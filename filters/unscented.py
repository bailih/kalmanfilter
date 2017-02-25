import numpy as np
import math
from lin_alg import extend_diagonal, extend_vector, check_pos_semi_def, matrix_square_root

ALPHA = 0.01
KI    = 0.0
BETA  = 2.0

DISCRETE_STRAT = 1

class Unscented:

    def __init__(self, order, Q, P, R, X, p, co):

        L = order + 1  # The dimension of the state variable
        self.L = L

        # We need 2L+1 sigma points to approximate the transformed gaussian distribution
        self.sigma_count = 1 + 2 * L

        self.Q = Q    # Process covariance
        self.P = P    # Estimation covariance
        self.R = R    # Measurement covariance
        self.X = X    # State information
        self.p = p    # Equivalent circuit paramaters
        self.co = co  # Polynomial coefficients

        self.c = (ALPHA * ALPHA * (L + KI))
        self.lamb = self.c - L

        # The zeroth weight for both mean and covariance are unique. All other weights are identical
        weights_mean = [self.lamb / self.c]
        weights_covariance = [(self.lamb / self.c) + 1 - (ALPHA * ALPHA) + BETA]
        weight = 1.0 / (2 * self.c)
        for i in range(1, 2 * L + 1):
            weights_mean.append(weight)
            weights_covariance.append(weight)

        self.wm = weights_mean
        self.wc = weights_covariance

        self._discrete_volts = {
            0 : self._volts_i_plus_a_t,
            1 : self._volts_inv_i_minus_a_t,
            2 : self._volts_bilinear
        }.get(DISCRETE_STRAT)


    def update(self, dt, voltage, current):

        current = 0 - current

        try:
            initial_sigmas = self._calculate_sigma_points(self.X, self.P)
        except TypeError:
            soc_index = self.L - 1
            if self.X[soc_index, 0] == None :
                # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                # in an OCV higher than the measured voltage; default to 1
                self.X[soc_index, 0] = 1
                for i in range(1, 10):
                    soc = i / 10
                    if voltage <= self._OCV(soc):
                        self.X[soc_index, 0] = soc
                        break
                initial_sigmas = self._calculate_sigma_points(self.X, self.P)
            else:
                raise TypeError("Initial state vector contains a None that's not at SOC position!")

        sigmas_star = self._priori_update(initial_sigmas, current, dt, self.p)
        # sigmas_star = []
        # for sigma in initial_sigmas:
        #     sigmas_star.append(self._priori_update(sigma, current, dt, self.p))

        priori_mean = self._calculate_mean(sigmas_star)

        priori_covariance = self._calculate_covariance(sigmas_star,
                                                       priori_mean,
                                                       sigmas_star,
                                                       priori_mean) \
                            + self.Q

        sigmas = self._calculate_sigma_points(priori_mean, priori_covariance)

        postiori_vt = []
        for sigma in sigmas:
            postiori_vt.append(self._calculate_terminal_voltage(sigma, current))

        postiori_mean = self._calculate_mean(postiori_vt)

        postiori_covariance = self._calculate_covariance(postiori_vt,
                                                         postiori_mean,
                                                         postiori_vt,
                                                         postiori_mean) \
                              + self.R

        combined_covariance = self._calculate_covariance(sigmas,
                                                         priori_mean,
                                                         postiori_vt,
                                                         postiori_mean)

        kalman_gain = combined_covariance * np.linalg.inv(postiori_covariance)

        self.X = priori_mean + kalman_gain * (voltage - postiori_mean)
        self.P = priori_covariance - kalman_gain * postiori_covariance * np.transpose(kalman_gain)


    def get_soc(self):
        return self.X[self.L - 1, 0]

    def get_estimated_voltage(self, current):
        return self._calculate_terminal_voltage(self.X, 0 - current)[0, 0]


    def _calculate_sigma_points(self, mean, covariance):
        '''
                mean, i = 0
        Xi = {  mean + ith column of sqrt((lambda + L) * covariance), 1 <= i <= L
                mean - (i - L)th column of sqrt((lambda + L) * covariance), L < i <= 2L

        '''
        sigma_points = [None] * self.sigma_count

        to_sqrt = self.c * covariance
        to_sqrt = check_pos_semi_def(to_sqrt, adjust=True)
        sqrt_matrix = np.linalg.cholesky(to_sqrt)

        sigma_points[0] = mean

        L = self.L
        for i in range(0, L):
            sigma_points[i + 1]     = mean + sqrt_matrix[:, [i]]
            sigma_points[i + 1 + L] = mean - sqrt_matrix[:, [i]]

        return sigma_points

    def _priori_update(self, points, current, dt, p):
        '''
        Vi_k+1  = (1 - dt / RiCi) * Vi_k + (dt / Ci) * current
        Soc_k+1 = Soc_k + (coulombic efficiency * (dt / capacity)) * current
        '''
        L = self.L
        sigmas = []
        for point in points:
            new_point = self._discrete_volts(self.L - 1, point, current, dt, p)
            new_point[L - 1, 0] = point[L - 1, 0] + (-1) * ((p['eta'] * dt) / p['cnom']) * current

            sigmas.append(new_point)
        return sigmas

    def _volts_i_plus_a_t(self, order, point, current, dt, p):
        new_point = np.matrix(np.zeros(point.shape))
        for i in range(order):
            cbr = p['cbranch'][i]
            new_point[i, 0] = (1 - dt / (p['rbranch'][i] * cbr)) * point[i, 0] \
                              + (dt / cbr) * current
        return new_point


    def _volts_inv_i_minus_a_t(self, order, point, current, dt, p):
        new_point = np.matrix(np.zeros(point.shape))
        for i in range(order) :
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            new_point[i, 0] = ((r * c) / (r * c + dt)) * point[i, 0] + (dt / c) * current

        return new_point

    def _volts_bilinear(self, order, point, current, dt, p):
        new_point = np.matrix(np.zeros(point.shape))
        for i in range(order):
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            new_point[i, 0] = (rc2 - dt) / (rc2 + dt) * point[i, 0]

        return new_point


    # Calculates the mean of a list of vectors, or a list of scalars.
    def _calculate_mean(self, vectors):
        '''mean = sum(Wi * vector_i) or sum(scalar_i)/2L for i = 0 to 2L'''
        if vectors[0].shape[0] > 1:
            mean = np.matrix(np.zeros((len(vectors[0]), 1)))
            is_matrix = True

        else:
            mean = 0
            is_matrix = False

        for i in range(len(vectors)):
            vector = vectors[i]

            if is_matrix:
                for j in range(len(vector)):
                    mean[j, 0] += self.wm[i] * vector[j, 0]

            else:
                mean += vector

        if not is_matrix:
            mean /= len(vectors)

        return mean


    def _calculate_covariance(self, points1, mean1, points2, mean2):
        '''cov = sum(Wi * (point1_i - mean1) * (point2_i - mean2)^T) for i = 0 to 2L'''
        covariance = np.matrix(np.zeros((len(points1[0]), len(points2[0]))))

        for i in range(len(points1)):
            weight, point1, point2 = self.wc[i], points1[i], points2[i]
            covariance += weight * ((point1 - mean1) * np.transpose(point2 - mean2))

        return covariance


    def _calculate_terminal_voltage(self, point, current):
        '''Vt = OCV(soc) - V1 - V2 - V3 - i * Ro'''
        volts, soc = (point[:-1], point[-1])
        ocv = self._OCV(soc)
        for volt in volts:
            ocv -= volt

        if current > 0:
            ocv -= current * self.p['rpos']
        else:
            ocv -= current * self.p['rneg']

        return ocv


    # Evaluated OCV_SOC curve
    def _OCV(self, soc):
        if soc > 1:
            b = self.co[-1]
            return (self._OCV(1) - b) * soc + b

        co_len = len(self.co)
        ocv = 0
        for i in range(co_len):
            ocv += self.co[i] * math.pow(soc, co_len - 1 - i)
        return ocv


class UnscentedSoh(Unscented):

    def __init__(self, order, Q, P, R, X, p, co):
        super(UnscentedSoh, self).__init__(order + 1, Q, P, R, X, p, co)

        self.Q = extend_diagonal(Q, 0.0001)
        self.P = extend_diagonal(P, 0.1)
        self.X = extend_vector(X, 0.5)

    def get_soc(self):
        return self.X[self.L - 2, 0]

    def get_soh(self):
        return self.X[self.L - 1, 0]

    def get_functional_soc(self):
        return self.X[self.L - 2, 0] / self.X[self.L - 1, 0]
        #return self.X[self.L - 1, 0]

    def _priori_update(self, points, current, dt, p):
        '''
        Vi_k+1  = (1 - dt / RiCi) * Vi_k + (dt / Ci) * current
        Soc_k+1 = Soc_k + (coulombic efficiency * (dt / capacity)) * current
        '''
        L = self.L
        sigmas = []
        for point in points:
            new_point = self._discrete_volts(self.L - 2, point, current, dt, p)
            new_point[L - 2, 0] = point[L - 2, 0] + (-1) * ((p['eta'] * dt) / p['cnom']) * current
            new_point[L - 1, 0] = point[L - 1, 0]

            sigmas.append(new_point)
        return sigmas

    def _calculate_terminal_voltage(self, point, current):
        '''Vt = OCV(soc) - V1 - V2 - V3 - i * Ro'''
        volts, soc, soh = (point[:-2], point[-2], point[-1])
        ocv = self._soh_OCV(soc, soh)
        for volt in volts:
            ocv -= volt

        if current > 0:
            ocv -= current * self.p['rpos']
        else:
            ocv -= current * self.p['rneg']

        return ocv

    # Evaluated OCV_SOC curve
    def _soh_OCV(self, soc, soh):
        co_len = len(self.co)
        ocv = 0
        for i in range(co_len):
            ocv += self.co[i] * math.pow(soc / soh, co_len - 1 - i)
        return ocv


class UnscentedLut(Unscented):

    def __init__(self, L, Q, P, R, X, p, co):
        super(UnscentedLut, self).__init__(L, Q, P, R, X, p, co)
        self.soc_list = []
        self.ocv_list = []
        with open('volt_bst.txt', 'r') as f:
            for line in f:
                soc, ocv = map(self._float_or_none, line.split(' '))
                self.soc_list.append(soc)
                self.ocv_list.append(ocv / 1000 if ocv is not None else None)

    def _OCV(self, soc):

        index = 1
        while self.soc_list[index] != soc:

            left_exists = (2 * index + 1 <= len(self.soc_list)
                           and self.soc_list[2 * index] is not None)
            right_exists = (2 * index + 2 <= len(self.soc_list)
                            and self.soc_list[2 * index + 1] is not None)

            if left_exists and self.soc_list[index] > soc:
                index = 2 * index
            elif right_exists and self.soc_list[index] < soc:
                index = 2 * index + 1

            else:
                break

        return self.ocv_list[index]


    def _float_or_none(self, string):
        try:
            value = float(string)
        except ValueError:
            value = None
        return value