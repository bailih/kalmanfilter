"""
    Header File for Kalman Filter SoC Estimator

    Contains Kalman Filter class and methods for 1st,
    2nd and 3rd order circuits.

"""

# Python Imports
# =============================================================================
from copy import deepcopy
from math import pow

import numpy as np

from classes.matrixHolder import matrixHolder, matrixHolder_Base
from lin_alg import extend_diagonal, extend_vector

# Constants
# =============================================================================
SOC_HISTORY_COUNT = 50
PD_HISTORY_COUNT = 150

# Chooses how we discretize our model.
# 0 ~ (I + AT), 1 ~ (I - AT)^-1, 2 ~ (I + AT/2)(I - AT/2)^-1
# Typically, for our purposes, strat 1 is best, because F won't go negative
# with large sampling times.
DISCRETE_STRAT = 1


class Kalman:
    """ Class to hold kalman state data and methods """

    def __init__(self, order, Q, P, R, X, p, c, inisoc, iniV):
        self.order = order  # Order of equivalent circuit
        vals = matrixHolder_Base(p['cnom'], inisoc, iniV)

        self.order = order  # Order of equivalent circuit
        self.dim = len(vals.getX())
        print("dim: {}".format(self.dim))

        self.case = vals.getcase()
        print("case: {}".format(self.case))

        baseH = np.matrix(np.zeros((1, order + 2)))
        baseK = np.matrix(np.zeros((order + 2, order + 2)))

        self.p = []
        self.c = []

        self.R = vals.getR()
        self.Q = vals.getQ()
        self.P = vals.getP()
        self.A = 1
        self.X = vals.getX()

        self.K = []
        self.H = []

        self.initial_x = []
        self.predicted_x = []
        self.new_x = []
        self.predicted_p = []
        self.new_p = []

        self.soc_error = []
        self.recent_socs = []
        self.recent_soc_diffs = []
        self.conver_step = []
        self.voltage_error = []
        self.simulated_voltage = []

        self.temp = []
        self.f = []
        self.b = []
        self.absolute_index = []

        self.ceiling_floor_hits = []
        for i in range(self.dim):
            self.absolute_index.append(self.case[i])
            self.ceiling_floor_hits.append(0)
            self.K.append(baseK)
            self.H.append(baseH)
            self.p.append(p)
            self.c.append(c)
            self.initial_x.append([])
            self.predicted_x.append([])
            self.new_x.append([])
            self.predicted_p.append([])
            self.new_p.append([])
            self.soc_error.append(0)
            self.recent_socs.append([])
            self.recent_soc_diffs.append([])
            self.conver_step.append(0)
            self.voltage_error.append(0)
            self.simulated_voltage.append(0)
            self.temp.append([])
            self.f.append([])
            self.b.append([])

        self._discrete_f = {
            0: self._volts_i_plus_a_t_new,
            1: self._volts_inv_i_minus_a_t_new,
            2: self._volts_bilinear_new
        }.get(DISCRETE_STRAT)

        # Kalman gain matrix
        self.iter = 0
        self.list_to_pop = []
        self.inisoc = inisoc

    def update(self, dt, voltage, current):
        order = self.order
        dim = self.dim
        self.list_to_pop = []

        if self.iter == 0:
            for i in range(dim):
                if self.X[i][self.order + 1, 0] == -1:
                    if abs(current) < 0.02:
                        self.X[i][self.order + 1, 0] = self._get_init_soc_new(voltage, i)
                    else:
                        self.X[i][self.order + 1, 0] = self._get_init_soc_new(voltage - 0.2, i)

        self.iter += 1

        # Reverse current direction
        current *= -1

        for i in range(dim):
            self.f[i], self.b[i] = self._create_f_and_b_new(dt, i)

        # Priori estimate (PREDICT)
        # Predict state:
        # X = F*X + B*u
        try:
            for i in range(dim):
                self.X[i] = np.add(self.f[i] * self.X[i], self.b[i] * current)
        except TypeError:
            for i in range(dim):
                if self.X[i][order, 0] == None:
                    # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                    # in an OCV higher than the measured voltage; default to 1
                    self.X[i][order, 0] = self.inisoc
                    # for j in range(1, 10):
                    #     soc = j / 10
                    #     if voltage < self._OCV(soc):
                    #         self.X[i] = soc * self.p[i]['cnom']
                    #         break
                    self.X[i] = np.add(self.f[i] * self.X[i], self.b[i] * current)
            else:
                raise TypeError("")

        # self.A = 0.98

        # Predict covariance:
        # P = F*P*Ft + Q
        for i in range(dim):
            self.predicted_p[i] = np.add(self.A * (self.f[i] * self.P[i] * self.f[i].T), self.Q[i])

        # Calculate H for posteriori estimate
        temp = np.matrix(np.zeros((1, order + 2)))
        for i in range(dim):
            for j in range(order):
                temp[0, j] = -1
            temp[0, order] = self._dOCVdIN(self.X[i][order, 0], self.X[i][order + 1, 0], self.p[i]['cnom'], i)
            temp[0, order + 1] = self._dOCVdINISOC(self.X[i][order, 0], self.X[i][order + 1, 0], self.p[i]['cnom'], i)

            self.H[i] = deepcopy(temp)

        # Calculate measurement error
        for i in range(dim):
            self.voltage_error[i] = (
            voltage - self._get_simulated_voltage_new(current, self.X[i][order, 0], self.X[i][order + 1, 0],
                                                      self.p[i]['cnom'], self.X[i], self.p[i], i))

            if abs(self.voltage_error[i]) < 0.1 and self.conver_step[i] == 0:
                self.conver_step[i] = self.iter
                print("converge step ({}): {}".format(self.absolute_index[i], self.conver_step[i]))

        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        for i in range(dim):
            self.temp[i] = self.predicted_p[i] * self.H[i].T
            self.K[i] = self.temp[i] * np.reciprocal(np.add((self.H[i] * self.temp[i]), self.R[i]))

        for i in range(dim):
            if self.X[i][self.order + 1, 0] < 0:
                self.X[i][self.order + 1, 0] = 0.0001

            if self.X[i][self.order + 1, 0] > self.p[i]['cnom']:
                self.X[i][self.order + 1, 0] = self.p[i]['cnom']

        # Update state:
        # X = X + error*K
        for i in range(dim):
            self.new_x[i] = np.add(self.X[i], (self.K[i] * self.voltage_error[i]))

        # Update covariance:
        # P = (I - K*H)*P*(I - K*H).T + K*R*K.T
        for i in range(dim):
            self.new_p[i] = np.subtract(np.eye(order + 2), (self.K[i] * self.H[i])) * self.predicted_p[i] * np.subtract(
                np.eye(order + 2), (self.K[i] * self.H[i])).T
            self.new_p[i] += self.K[i] * self.R[i] * self.K[i].T
            self.new_p[i] = (self.new_p[i] + self.new_p[i].T) / 2
            self.X[i] = self.new_x[i]
            self.P[i] = self.new_p[i]

        for i in range(dim):
            if self.X[i][self.order + 1, 0] < 0:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 1, 0] = 0.0001

            if self.X[i][self.order + 1, 0] > self.p[i]['cnom']:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 1, 0] = self.p[i]['cnom']

        # Recalculate simulated_voltage for printing purposes
        for i in range(dim):
            self.simulated_voltage[i] = self._get_simulated_voltage_new(current, self.X[i][order, 0],
                                                                        self.X[i][order + 1, 0], self.p[i]['cnom'],
                                                                        self.X[i], self.p[i], i)

            sim_soc = self.get_inputsoc(i)
            self.recent_socs[i].append(sim_soc)
            if len(self.recent_socs[i]) > SOC_HISTORY_COUNT:
                self.recent_socs[i].pop(0)
            mean = 0
            for j in range(len(self.recent_socs[i])):
                mean += self.recent_socs[i][j]
            mean = mean / (j + 1)

            self.recent_soc_diffs[i].append(abs(((mean - self.get_inputsoc(i)) / mean) * 100))
            if len(self.recent_soc_diffs[i]) > PD_HISTORY_COUNT:
                self.recent_soc_diffs[i].pop(0)

    def _create_f_and_b_new(self, dt, version):
        p = self.p[version]
        order = self.order
        # State transistion matrix
        f = self._discrete_f(np.matrix(np.zeros((order + 2, order + 2))), dt, version)
        f[order, order] = 1  # Last element = 1
        f[order + 1, order + 1] = 1

        # Control matrix
        b = np.zeros((order + 2, 1))
        for i in range(self.order):
            b[i, 0] = dt / p['cbranch'][i]
        b[order, 0] = (-1) * (p['eta'] * dt)  # Last element
        b[order + 1, 0] = 0

        return f, b

    def _volts_i_plus_a_t_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):  # Reference all but last diagonal
            f[i, i] = 1 - dt / (p['rbranch'][i] * p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):  # Reference all but last diagonal
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):  # Reference all but last diagonal
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f
        # Evaluated OCV_SOC curve

    def get_amp_seconds_new(self, version):
        return self.X[version][self.order, 0]

    def get_inputsoc(self, version):
        return self.X[version][self.order, 0]

    def get_inputsoc_2(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order, 0])
        return temp

    def get_P(self):

        temp = []
        for i in range(len(self.P)):
            temp.append(self.P[i])
        return temp

    def get_inisoc(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order + 1, 0])
        return temp

    def get_estimated_voltage_new(self, version):
        return self.simulated_voltage[version]

    def get_estimated_voltage_new_2(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.simulated_voltage[i])
        return temp

    def get_mean_percent_diff_new(self, version):
        mean = 0
        for i in range(len(self.recent_soc_diffs[version])):
            mean += self.recent_soc_diffs[version][i]
        return mean / (i + 1)

    # Evaluate derivative of OCV_SOC curve
    def _dOCVdIN(self, ini, input, soh, version):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        co_len = len(self.c[version])
        docv = 0
        for i in range(co_len - 1):
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) + (co_len - 2 - i))
            try:
                docv += (co_len - 1 - i) * self.c[version][i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (
                co_len - 1 - i))
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return docv

    # Evaluate derivative of OCV_SOC curve
    def _dOCVdINISOC(self, ini, input, soh, version):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        co_len = len(self.c[version])
        docv = 0
        for i in range(co_len - 1):
            # print ("i : {}, {}".format(i, self.c[i]))
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) * 2)
            try:
                docv += (co_len - 1 - i) * self.c[version][i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (
                co_len - 1 - i))
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return docv

    # Evaluated OCV_SOC curve
    def _OCV_C(self, ini, input, soh, version):
        co_len = len(self.c[version])
        ocv = 0
        for i in range(co_len):
            # ocv += self.c[i] * pow((ini/self.p['cnom']/soh/soh + input/self.p['cnom']/soh), co_len - 1 - i)
            try:
                ocv += self.c[version][i] * pow((ini / soh + input / soh), co_len - 1 - i)
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return ocv

    # Gets voltage error between predicted and simulated voltage
    def _get_simulated_voltage_new(self, current, ini, input, soh, X, p, version):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        simulated_voltage = self._OCV_C(ini, input, soh, version)
        for i in range(self.order):
            simulated_voltage -= X[i, 0]
        if current > 0:
            simulated_voltage -= current * p['rpos']
        else:
            simulated_voltage -= current * p['rneg']
        return simulated_voltage

    def get_simulated_voltage_new(self):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.simulated_voltage[i])
        return temp

    def get_simulated_soc_new(self):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order, 0])
        return temp

    def _get_init_soh(self, voltage, version):
        init_voltage = voltage
        smallest_dif = 0.1
        best_soh = self.p[0]['cnom']
        ini = self.get_inisoc(version)
        input = self.get_inputsoc(version)

        for i in range(100):
            soh = i * 90 + 0.01
            co_len = len(self.c[version])
            ocv = 0
            for j in range(co_len):
                ocv += self.c[version][j] * pow((ini / soh + input / soh), co_len - 1 - j)

            temp = (ocv - init_voltage)
            if temp < smallest_dif and temp > 0:
                best_soh = soh  # soc (%) * capacity
                break

        print('Init SoC_input guess: {}'.format(input))
        print('Init SoC_ini guess: {}'.format(ini))

        return best_soh  # for testing, return 0

    def _OCV_init_new(self, i, soc):
        co_len = len(self.c[i])
        ocv = 0
        for j in range(co_len):
            ocv += self.c[i][j] * pow(soc, co_len - 1 - j)
        return ocv

    def _get_init_soc_new(self, voltage, k):
        init_voltage = voltage
        temp_soc = 0.00001
        smallest_dif = 99999
        best_soc = 0.000000123
        for i in range(10000):
            temp = abs(self._OCV_init_new(k, temp_soc) - init_voltage)
            if temp < smallest_dif:
                smallest_dif = temp
                best_soc = temp_soc * self.p[k]['cnom']  # soc (%) * capacity
            temp_soc = (i + 1) / 10000

        print('Init SoC guess: {} ({})'.format(best_soc, best_soc / self.p[k]['cnom']))
        return best_soc


class KalmanSoC(Kalman):
    """ Class to hold kalman state data and methods
        This is the super filter
        state variable: V1, inputSoC, SoH, initSoC
    """

    def __init__(self, order, Q, P, R, X, p, c, inisoc, iniV):
        super(KalmanSoC, self).__init__(order, Q, P, R, X, p, c, inisoc, iniV)

        vals = matrixHolder(p['cnom'], inisoc, iniV)

        self.order = order  # Order of equivalent circuit
        self.dim = len(vals.getX())
        print("dim: {}".format(self.dim))

        self.case = vals.getcase()
        print("case: {}".format(self.case))
        # self.R = R  # Measurement covariance (R)
        # self.p = p  # Equivalent circuit paramaters (p)
        # self.c = c  # Polynomial coefficients

        baseH = np.matrix(np.zeros((1, order + 3)))
        baseK = np.matrix(np.zeros((order + 3, order + 3)))

        self.p = []
        self.c = []

        self.R = vals.getR()
        self.Q = vals.getQ()
        self.P = vals.getP()
        self.A = 1
        self.X = vals.getX()

        self.K = []
        self.H = []

        self.initial_x = []
        self.predicted_x = []
        self.new_x = []
        self.predicted_p = []
        self.new_p = []

        self.soc_error = []
        self.recent_socs = []
        self.recent_soc_diffs = []
        self.conver_step = []
        self.voltage_error = []
        self.simulated_voltage = []

        self.temp = []
        self.f = []
        self.b = []
        self.absolute_index = []

        self.ceiling_floor_hits = []
        for i in range(self.dim):
            self.absolute_index.append(self.case[i])
            self.ceiling_floor_hits.append(0)
            self.K.append(baseK)
            self.H.append(baseH)
            self.p.append(p)
            self.c.append(c)
            self.initial_x.append([])
            self.predicted_x.append([])
            self.new_x.append([])
            self.predicted_p.append([])
            self.new_p.append([])
            self.soc_error.append(0)
            self.recent_socs.append([])
            self.recent_soc_diffs.append([])
            self.conver_step.append(0)
            self.voltage_error.append(0)
            self.simulated_voltage.append(0)
            self.temp.append([])
            self.f.append([])
            self.b.append([])

        self._discrete_f = {
            0: self._volts_i_plus_a_t_new,
            1: self._volts_inv_i_minus_a_t_new,
            2: self._volts_bilinear_new
        }.get(DISCRETE_STRAT)

        self.iter = 0
        self.list_to_pop = []

    # Kalman update function
    def update(self, dt, voltage, current):
        order = self.order
        dim = self.dim
        self.list_to_pop = []

        if self.iter == 0:
            for i in range(dim):
                if self.X[i][self.order + 1, 0] == -1:
                    if abs(current) < 0.02:
                        self.X[i][self.order + 1, 0] = self._get_init_soc_new(voltage, i)
                    else:
                        self.X[i][self.order + 1, 0] = self._get_init_soc_new(voltage - 0.2, i)

        self.iter += 1

        # Reverse current direction
        current *= -1

        for i in range(dim):
            self.f[i], self.b[i] = self._create_f_and_b_new(dt, i)

        # Priori estimate (PREDICT)
        # Predict state:
        # X = F*X + B*u
        try:
            for i in range(dim):
                self.X[i] = np.add(self.f[i] * self.X[i], self.b[i] * current)
        except TypeError:
            for i in range(dim):
                if self.X[i][order, 0] == None:
                    # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                    # in an OCV higher than the measured voltage; default to 1
                    self.X[i][order, 0] = 1
                    for j in range(1, 10):
                        soc = j / 10
                        if voltage < self._OCV(soc):
                            self.X[i] = soc * self.p[i]['cnom']
                            break
                    self.X[i] = np.add(self.f[i] * self.X[i], self.b[i] * current)
            else:
                raise TypeError("")

        # esttimate soh and soc, we know soc_init + soc = soh
        # make it time varying
        # for i in range(dim):
        #     diff = self.X[i][order + 2, 0] - (self.X[i][order + 1, 0] + self.X[i][order, 0])
        #     weight = self.iter/10000
        #     weight2 = 0.3
        #
        #     self.X[i][order + 2, 0] -= diff * weight * weight2
        #     self.X[i][order + 1, 0] += diff * (1 - weight) * (1 - weight2)

        # self.A = 0.98

        # Predict covariance:
        # P = F*P*Ft + Q
        for i in range(dim):
            self.predicted_p[i] = np.add(self.A * (self.f[i] * self.P[i] * self.f[i].T), self.Q[i])

        # Calculate H for posteriori estimate
        temp = np.matrix(np.zeros((1, order + 3)))
        for i in range(dim):
            for j in range(order):
                temp[0, j] = -1
            temp[0, order] = self._dOCVdIN(self.X[i][order, 0], self.X[i][order + 1, 0], self.X[i][order + 2, 0], i)
            temp[0, order + 1] = self._dOCVdINISOC(self.X[i][order, 0], self.X[i][order + 1, 0],
                                                   self.X[i][order + 2, 0], i)
            temp[0, order + 2] = self._dOCVdSOH(self.X[i][order, 0], self.X[i][order + 1, 0], self.X[i][order + 2, 0],
                                                i)
            self.H[i] = deepcopy(temp)

        # Calculate measurement error
        for i in range(dim):
            self.voltage_error[i] = (
            voltage - self._get_simulated_voltage_new(current, self.X[i][order, 0], self.X[i][order + 1, 0],
                                                      self.X[i][order + 2, 0], self.X[i], self.p[i], i))

            if abs(self.voltage_error[i]) < 0.1 and self.conver_step[i] == 0:
                self.conver_step[i] = self.iter
                print("converge step ({}): {}".format(self.absolute_index[i], self.conver_step[i]))

        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        for i in range(dim):
            self.temp[i] = self.predicted_p[i] * self.H[i].T
            self.K[i] = self.temp[i] * np.reciprocal(np.add((self.H[i] * self.temp[i]), self.R[i]))

        for i in range(dim):
            if self.X[i][self.order + 1, 0] < 0:
                self.X[i][self.order + 1, 0] = 0.0001

            if self.X[i][self.order + 1, 0] > self.p[i]['cnom']:
                self.X[i][self.order + 1, 0] = self.p[i]['cnom']

            if self.X[i][self.order + 2, 0] > self.p[i]['cnom']:
                self.X[i][self.order + 2, 0] = self.p[i]['cnom']

            if self.X[i][self.order + 2, 0] < 0:
                self.X[i][self.order + 2, 0] = 0.0001
        # Update state:
        # X = X + error*K
        for i in range(dim):
            self.new_x[i] = np.add(self.X[i], (self.K[i] * self.voltage_error[i]))

        # Update covariance:
        # P = (I - K*H)*P*(I - K*H).T + K*R*K.T
        for i in range(dim):
            self.new_p[i] = np.subtract(np.eye(order + 3), (self.K[i] * self.H[i])) * self.predicted_p[i] * np.subtract(
                np.eye(order + 3), (self.K[i] * self.H[i])).T
            self.new_p[i] += self.K[i] * self.R[i] * self.K[i].T
            self.new_p[i] = (self.new_p[i] + self.new_p[i].T) / 2
            self.X[i] = self.new_x[i]
            self.P[i] = self.new_p[i]

        for i in range(dim):
            if self.X[i][self.order + 1, 0] < 0:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 1, 0] = 0.0001

            if self.X[i][self.order + 1, 0] > self.p[i]['cnom']:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 1, 0] = self.p[i]['cnom']

            if self.X[i][self.order + 2, 0] > self.p[i]['cnom']:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 2, 0] = self.p[i]['cnom']

            if self.X[i][self.order + 2, 0] < 0:
                self.ceiling_floor_hits[i] += 1
                self.X[i][self.order + 2, 0] = 0.0001

        # Recalculate simulated_voltage for printing purposes
        for i in range(dim):
            self.simulated_voltage[i] = self._get_simulated_voltage_new(current, self.X[i][order, 0],
                                                                        self.X[i][order + 1, 0],
                                                                        self.X[i][order + 2, 0], self.X[i], self.p[i], i)

            sim_soc = self.get_inputsoc(i)
            self.recent_socs[i].append(sim_soc)
            if len(self.recent_socs[i]) > SOC_HISTORY_COUNT:
                self.recent_socs[i].pop(0)
            mean = 0
            for j in range(len(self.recent_socs[i])):
                mean += self.recent_socs[i][j]
            mean = mean / (j + 1)

            self.recent_soc_diffs[i].append(abs(((mean - self.get_inputsoc(i)) / mean) * 100))
            if len(self.recent_soc_diffs[i]) > PD_HISTORY_COUNT:
                self.recent_soc_diffs[i].pop(0)


                # ''' popping from the list should occur here'''
                # for i in range(dim):
                #     if self.ceiling_floor_hits[i] >= 1000:
                #         self.list_to_pop.append(i)

                # for i in sorted(self.list_to_pop, reverse=True):
                #     print('Popping cases: '),
                #     for set_ in self.list_to_pop:
                #         print(self.absolute_index[set_])
                #     self.dim -= len(self.list_to_pop)
                #     self.absolute_index.pop(i)
                #     self.ceiling_floor_hits.pop(i)
                #     self.R.pop(i)
                #     self.Q.pop(i)
                #     self.P.pop(i)
                #     self.A.pop(i)
                #     self.X.pop(i)
                #     self.K.pop(i)
                #     self.H.pop(i)
                #     self.p.pop(i)
                #     self.c.pop(i)
                #     self.initial_x.pop(i)
                #     # self.predicted_x.pop(i)
                #     self.new_x.pop(i)
                #     self.predicted_p.pop(i)
                #     self.new_p.pop(i)
                #     self.soc_error.pop(i)
                #     self.recent_socs.pop(i)
                #     self.recent_soc_diffs.pop(i)
                #     self.conver_step.pop(i)
                #     self.voltage_error.pop(i)
                #     self.simulated_voltage.pop(i)
                #     self.temp.pop(i)
                #     self.f.pop(i)
                #     self.b.pop(i)

    def _create_f_and_b_new(self, dt, version):
        p = self.p[version]
        order = self.order
        # State transistion matrix
        f = self._discrete_f(np.matrix(np.zeros((order + 3, order + 3))), dt, version)
        f[order, order] = 1  # Last element = 1
        f[order + 1, order + 1] = 1
        f[order + 2, order + 2] = 1

        # Control matrix
        b = np.zeros((order + 3, 1))
        for i in range(self.order):
            b[i, 0] = dt / p['cbranch'][i]
        b[order, 0] = (-1) * (p['eta'] * dt)  # Last element
        b[order + 1, 0] = 0
        b[order + 2, 0] = 0
        return f, b

    def _volts_i_plus_a_t_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):     # Reference all but last diagonal
            f[i, i] = 1 - dt / (p['rbranch'][i] * p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):     # Reference all but last diagonal
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear_new(self, f, dt, version):
        p = self.p[version]
        for i in range(self.order):     # Reference all but last diagonal
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f

    def get_amp_seconds_new(self, version):
        return self.X[version][self.order, 0]

    def get_inputsoc(self, version):
        return self.X[version][self.order, 0]

    def get_inputsoc_2(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order, 0])
        return temp

    def get_P(self):

        temp = []
        for i in range(len(self.P)):
            temp.append(self.P[i])
        return temp

    def get_inisoc(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order + 1, 0])
        return temp

    def get_soh(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order + 2, 0])
        return temp

    def get_estimated_voltage_new(self, version):
        return self.simulated_voltage[version]

    def get_estimated_voltage_new_2(self, version):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.simulated_voltage[i])
        return temp

    def get_mean_percent_diff_new(self, version):
        mean = 0
        for i in range(len(self.recent_soc_diffs[version])):
            mean += self.recent_soc_diffs[version][i]
        return mean / (i + 1)


    # Evaluate derivative of OCV_SOC curve
    def _dOCVdIN(self, ini, input, soh, version):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        co_len = len(self.c[version])
        docv = 0
        for i in range(co_len - 1):
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) + (co_len - 2 - i))
            try:
                docv += (co_len - 1 - i) * self.c[version][i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (
                co_len - 1 - i))
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return docv

    # Evaluate derivative of OCV_SOC curve
    def _dOCVdINISOC(self, ini, input, soh, version):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        co_len = len(self.c[version])
        docv = 0
        for i in range(co_len - 1):
            # print ("i : {}, {}".format(i, self.c[i]))
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) * 2)
            try:
                docv += (co_len - 1 - i) * self.c[version][i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (
                co_len - 1 - i))
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return docv

        # Evaluate derivative of OCV_SOC curve

    def _dOCVdSOH(self, ini, input, soh, version):
        # x = (self.get_inputsoc() + self.get_inisoc()) / self.p['cnom']  # Get state of charge, convert to %
        co_len = len(self.c[version])
        docv = 0
        for i in range(co_len - 1):
            # docv -= (co_len - i) * self.c[i] * pow((input * soh) + ini, co_len - 2 - i) * ((input * soh) + ini * 2) / pow(self.p['cnom'], co_len - i) / pow(soh, (co_len - 1 - i)*2 + 1)
            try:
                docv -= (co_len - 1 - i) * self.c[version][i] * pow((ini + input), (co_len - 1 - i)) / pow(soh,
                                                                                                           (co_len - i))
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return docv

    # Evaluated OCV_SOC curve
    def _OCV_C(self, ini, input, soh, version):
        co_len = len(self.c[version])
        ocv = 0
        for i in range(co_len):
            # ocv += self.c[i] * pow((ini/self.p['cnom']/soh/soh + input/self.p['cnom']/soh), co_len - 1 - i)
            try:
                ocv += self.c[version][i] * pow((ini / soh + input / soh), co_len - 1 - i)
            except OverflowError or ZeroDivisionError:
                if version not in self.list_to_pop:
                    self.list_to_pop.append(version)
                break
        return ocv

    # Gets voltage error between predicted and simulated voltage
    def _get_simulated_voltage_new(self, current, ini, input, soh, X, p, version):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        simulated_voltage = self._OCV_C(ini, input, soh, version)
        for i in range(self.order):
            simulated_voltage -= X[i, 0]
        if current > 0:
            simulated_voltage -= current * p['rpos']
        else:
            simulated_voltage -= current * p['rneg']
        return simulated_voltage

    def get_simulated_voltage_new(self):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.simulated_voltage[i])
        return temp

    def get_simulated_soc_new(self):
        temp = []
        for i in range(len(self.X)):
            temp.append(self.X[i][self.order, 0])
        return temp

    def _get_init_soh(self, voltage, version):
        init_voltage = voltage
        smallest_dif = 0.1
        best_soh = self.p[0]['cnom']
        ini = self.get_inisoc(version)
        input = self.get_inputsoc(version)

        for i in range(100):
            soh = i * 90 + 0.01
            co_len = len(self.c[version])
            ocv = 0
            for j in range(co_len):
                ocv += self.c[version][j] * pow((ini / soh + input / soh), co_len - 1 - j)

            temp = (ocv - init_voltage)
            if temp < smallest_dif and temp > 0:
                best_soh = soh  # soc (%) * capacity
                break

        print('Init SoC_input guess: {}'.format(input))
        print('Init SoC_ini guess: {}'.format(ini))
        print('Init SoH guess: {}'.format(best_soh))
        return best_soh  # for testing, return 0

    def _OCV_init_new(self, i, soc):
        co_len = len(self.c[i])
        ocv = 0
        for j in range(co_len):
            ocv += self.c[i][j] * pow(soc, co_len - 1 - j)
        return ocv

    def _get_init_soc_new(self, voltage, k):
        init_voltage = voltage
        temp_soc = 0.00001
        smallest_dif = 99999
        best_soc = 0.000000123
        for i in range(10000):
            temp = abs(self._OCV_init_new(k, temp_soc) - init_voltage)
            if temp < smallest_dif:
                smallest_dif = temp
                best_soc = temp_soc * self.p[k]['cnom']  # soc (%) * capacity
            temp_soc = (i + 1) / 10000

        print('Init SoC guess: {} ({})'.format(best_soc, best_soc / self.p[k]['cnom']))
        return best_soc


class KalmanSoh(Kalman):
    '''
        Alternate implementation of the Kalman filter for SOC estimation which attempts to also
        estimate the soh of the battery.

        For batteries of low soh, the estimated soc will be much closer to the real soc.
        Without keeping track of soh, the filter will estimate closer to the 'practical soc':
        The soc with it's current maximum capacity as 100%, not its best maximum capacity.
    '''

    def __init__(self, order, Q, P, R, X, p, c):
        super(KalmanSoh, self).__init__(order, Q, P, R, X, p, c)

        self.Q = extend_diagonal(Q, 0.00001)
        self.P = extend_diagonal(P, 100)
        self.X = extend_vector(X, 1)
        self.A = 1

        self.Temp = np.zeros((order + 2, order + 2))

        self.iter = 0

    # Kalman update function
    def update(self, dt, voltage, current):

        if (self.iter == 0):
            self.X[self.order, 0] = self._get_init_soc(voltage)

        self.iter += 1

        order = self.order
        # Reverse current direction
        current *= -1
        f = self._soh_f(dt, order)

        # Priori estimate (PREDICT)
        # Predict state:
        initial_x = self.X
        predicted_x = f * initial_x

        # if predicted_x[order + 1, 0] > 1:
        #     predicted_x[order + 1, 0] = 0.99

        # x = x + bi
        for i in range(order):
            predicted_x[i, 0] += (dt * current) / self.p['cbranch'][i]
        predicted_x[order, 0] -= (dt * current)

        # predict covariance
        predicted_p = np.add((f * self.P * f.T), self.Q)

        # calculate h
        h = np.matrix(np.zeros((1, order + 2)))
        for i in range(order):
            h[0, i] = -1
        h[0, order] = self._dOCVdSOC(predicted_x, self.p['cnom'])
        h[0, order + 1] = self._dOCVdSOH(predicted_x, self.p['cnom'])

        # get voltage error
        self.voltage_error = (voltage - self._get_simulated_voltage(current))

        # will memory fade the filter
        # if abs(self.voltage_error) <= 0.01:
        #
        #     self.A = 1
        # else:
        if current > 0:
            self.A = 0.98
        else:
            self.A = 1.01

        # Posteriori estimate (UPDATE)

        # Update state:
        # X = X + error*K
        self.Temp = self.P * h.T
        temp = h.T * np.reciprocal(np.add((h * self.Temp), self.R))
        k = self.P * temp

        new_x = np.add(predicted_x, k * self.voltage_error)

        temp = predicted_p * (np.subtract(np.eye(order + 2), (k * h))).T
        updated_p = np.subtract(np.eye(self.order + 2), (k * h)) * temp + (k * self.R * k.T)

        new_p = (updated_p + updated_p.T) / 2

        self.P = new_p
        self.X = new_x

        self.simulated_voltage = self._get_simulated_voltage(current)

        sim_soc = self.get_simulated_soc()
        self.recent_socs.append(sim_soc)
        if len(self.recent_socs) > SOC_HISTORY_COUNT:
            self.recent_socs.pop(0)
        mean = 0
        for i in range(len(self.recent_socs)):
            mean += self.recent_socs[i]
        mean = mean / (i + 1)

        self.recent_soc_diffs.append(abs(((mean - self.get_simulated_soc()) / mean) * 100))
        if len(self.recent_soc_diffs) > PD_HISTORY_COUNT:
            self.recent_soc_diffs.pop(0)

    def _dOCVdSOC(self, X, capacity):
        order = self.order
        soc, soh = X[order, 0], X[order + 1, 0]
        c = self.c
        co_len = len(c)
        soc /= capacity

        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * c[i] * (pow(soc, co_len - 2 - i)) / pow(soh, co_len - 1 - i)
        return docv

    def _dOCVdSOH(self, X, capacity):
        order = self.order
        soc, soh = X[order, 0], X[order + 1, 0]
        soc /= capacity
        c = self.c
        co_len = len(c)

        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * c[i] * (pow(soc, co_len - 1 - i)) / pow(soh, co_len - i)
        return docv

    # def _dOCVdSOC(self):
    #     order = self.order
    #     soc, soh = self.X[order, 0], self.X[order + 1, 0]
    #     c = self.c
    #     co_len = len(c)
    #
    #     docv = 0
    #     for i in range(co_len - 1):
    #         docv += (co_len - 1 - i) * c[i] * (pow(soc, co_len - 2 - i) / pow(soh, co_len - 1 - i))
    #     return docv
    #
    # def _dOCVdSOH(self):
    #     order = self.order
    #     soc, soh = self.X[order, 0], self.X[order + 1, 0]
    #     c = self.c
    #     co_len = len(c)
    #
    #     docv = 0
    #     for i in range(co_len - 1):
    #         docv -= (co_len - 1 - i) * c[i] * (pow(soc, co_len - 1 - i) / pow(soh, co_len - i))
    #     return docv


    # Evaluated OCV_SOC curve ONLY ACCEPTS % AS INPUT
    def _OCV(self, soc=None):
        x = self.get_simulated_soc() / self.p['cnom'] if soc == None else soc  # Get state of charge
        soh = self.X[self.order + 1, 0]

        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            ocv += self.c[i] * pow(x / soh, co_len - 1 - i)
        return ocv

    def _soh_f(self, dt, order):
        # For this subclass, we don't actually create b
        # State transistion matrix
        f = self._discrete_f(np.matrix(np.zeros((order + 2, order + 2))), dt)
        f[order, order] = 1  # Last element = 1
        f[order + 1, order + 1] = 1

        return f

    def get_soh(self):
        return self.X[self.order + 1, 0]
