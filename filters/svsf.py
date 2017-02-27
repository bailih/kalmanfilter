import numpy as np
from math import pow, e
from lin_alg import element_wise_multiplication

DISCRETE_STRAT = 1

class SVSF:
    """
        Smooth variable structure filter specifically for soc estimation.

        Written by Connor Rambo & Nelson Schultz for Cadex Electronics, 2016
    """

    def __init__(self, order, X, p, co, psi, gamma):
        """
            parameters
                order - (int) the order of the equivalent circuit model
                X   - (numpy.matrix) initial state vector
                p   - (kalman.Parameters) model parameters
                co  - (tuple) coefficients for the OCV(soc) curve
                psi - (float) the convergence rate of the svsf
                gamma - (float) the width of the smoothing boundary layer
            """
        L = order + 1
        self.L = L
        self.X = X
        self.p = p
        self.co = co
        self.psi = np.eye(L) * psi
        self.gamma = np.eye(L) * gamma
        self.e_k = 0
        self._discrete_f = {
            0 : self._volts_i_plus_a_t,
            1 : self._volts_inv_i_minus_a_t,
            2 : self._volts_bilinear
        }.get(DISCRETE_STRAT)



    def update(self, dt, voltage, current):

        # The model provided by the McMaster report uses backwards current.
        current *= -1

        try:
            priori_x = self._priori_update(self.X, current, dt, self.p)
        except TypeError:
            soc_index = self.L - 1
            if self.X[soc_index, 0] == None:
                # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                # in an OCV higher than the measured voltage; default to 1
                self.X[soc_index, 0] = 1
                for i in range(1, 100):
                    soc = i / 100
                    if voltage < self._OCV(soc):
                        self.X[soc_index, 0] = soc
                        break
                priori_x = self._priori_update(self.X, current, dt, self.p)
            else:
                raise TypeError("Initial state vector contains a None that's not at SOC position!")

        # We find the percent error between the voltage we measure and the voltage we would expect
        # based on our estimations
        priori_voltage = self._get_simulated_voltage(priori_x, current)
        e_k_1 = (voltage - priori_voltage)
        e_percent = e_k_1 / voltage

        # We then build a vector of state estimation errors by assuming the states differ by the
        # same percentage as the measurement
        error_vector = np.mat(np.zeros((self.L, 1)))
        for i in range(self.L):
            error_vector[i, 0] = priori_x[i, 0] * e_percent

        # eg. H = [-1, -1, -1, dOCV/soc]
        H = np.mat(np.zeros((1, self.L)))
        for i in range(self.L - 1):
            H[0, i] = -1
        H[0, self.L - 1] = self._dOCV(priori_x[self.L - 1, 0]) / priori_x[self.L - 1, 0]

        # The corrective gain is the pseudo inverse of the measurement matrix H, multiplied with
        # the weighted error term, all element wise multiplied with the vector returned by the sat
        # function, which provides the boundary layer.
        H_pseudo_inverted = np.linalg.pinv(H)

        pre_sat_gain = np.transpose(np.transpose(H_pseudo_inverted) *
                                    (abs(e_k_1) * np.eye(self.L) + self.gamma * abs(self.e_k)))
        sat = self._sat(error_vector)
        corrective_gain = element_wise_multiplication(pre_sat_gain, sat)

        self.X = priori_x + corrective_gain

        # We use PRESENT measurements to create the start points for FUTURE estimations
        sim_voltage = self._get_simulated_voltage(self.X, current)
        self.e_k = voltage - sim_voltage


    def get_soc(self):
        return self.X[self.L - 1, 0]

    def get_estimated_voltage(self, current):
        return self._get_simulated_voltage(self.X, 0 - current)

    # Takes in a point, the state vector, current, delta t and parameters to create the priori
    # estimate for the next state vector.
    def _priori_update(self, point, current, dt, p):

        f = self._discrete_f(self.L - 1, np.matrix(np.zeros((self.L, self.L))), dt, p)
        f[self.L - 1, self.L - 1] = 1  # Last element = 1

        # Control matrix
        b = np.matrix(np.zeros((self.L, 1)))
        for i in range(self.L - 1):
            b[i, 0] = dt / p['cbranch'][i]
        b[self.L - 1, 0] = (-1) * (p['eta'] * dt) / p['cnom']  # Last element

        new_point = f * point + b * current

        return new_point

    def _volts_i_plus_a_t(self, order, f, dt, p):
        for i in range(order):  # Reference all but last diagonal
            f[i, i] = 1 - dt / (p['rbranch'][i] * p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t(self, order, f, dt, p):
        for i in range(order):  # Reference all but last diagonal
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear(self, order, f, dt, p):
        for i in range(order) :  # Reference all but last diagonal
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f

    # Saturation function that provides the smoothing boundary layer
    def _sat(self, error):
        saturated = np.mat(np.zeros((self.L, 1)))

        for i in range(self.psi.shape[0]):
            saturated[i, 0] = error[i, 0] / self.psi[i, i]

            if saturated[i, 0] >= 1:
                saturated[i, 0] = 1

            elif saturated[i, 0] <= -1:
                saturated[i, 0] = -1

        return saturated

    # The less complicated SVSF only uses the sign of the error instead of the sat()
    # Left here for posterity
    def _signum(self, vector):
        sign_vector = np.mat(np.zeros((self.L, 1)))
        for i in range(sign_vector.shape[0]):
            if vector > 0:
                sign_vector[i, 0] = 1
            elif vector < 0:
                sign_vector[i, 0] = -1
        return sign_vector

    # Derivative of OCV(soc) curve at soc
    def _dOCV(self, soc):
        co_len = len(self.co)
        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * self.co[i] * pow(soc, co_len - 2 - i)
        return docv

    # Evaluated OCV_SOC curve
    def _OCV(self, soc):
        co_len = len(self.co)
        ocv = 0
        for i in range(co_len):
            ocv += self.co[i] * pow(soc, co_len - 1 - i)
        return ocv

    # This method assumes the current has already been properly reversed
    def _get_simulated_voltage(self, state_vector, current):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        soc = state_vector[self.L - 1]
        simulated_voltage = self._OCV(soc)
        for i in range(self.L - 1):
            simulated_voltage -= state_vector[i, 0]
        if current > 0:
            simulated_voltage -= current * self.p['rpos']
        else:
            simulated_voltage -= current * self.p['rneg']
        return simulated_voltage