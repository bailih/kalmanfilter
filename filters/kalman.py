"""
    Header File for Kalman Filter SoC Estimator

    Contains Kalman Filter class and methods for 1st,
    2nd and 3rd order circuits.

"""

# Python Imports
# =============================================================================
import math
from math import pow

import matplotlib.pyplot as plt
import numpy as np

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

    def __init__(self, order, Q, P, R, X, p, c):
        self.order = order  # Order of equivalent circuit
        self.Q = Q  # Process covariance  (Q)
        self.P = P  # Estimation covariance (P)
        self.R = R  # Measurement covariance (R)
        self.X = X  # State information (X)
        self.p = p  # Equivalent circuit paramaters (p)
        self.c = c  # Polynomial coefficients

        self._discrete_f = {
            0: self._volts_i_plus_a_t,
            1: self._volts_inv_i_minus_a_t,
            2: self._volts_bilinear
        }.get(DISCRETE_STRAT)

        # Kalman gain matrix
        self.K = np.zeros((order + 1, 1))
        # Measurement matrix
        self.H = np.zeros((1, order + 1))
        for i in range(order):
            self.H[0, i] = -1  # All but last column = '-1'
        # Temporary matrix (for calculation)
        self.Temp = np.zeros((order + 1, order + 1))
        # Useful scalar values
        self.measured_voltage = 0
        self.measured_current = 0
        self.simulated_voltage = 0
        self.voltage_error = 0
        self.soc_error = 0
        self.timer = 0

        self.recent_socs = []
        self.recent_soc_diffs = []

        self.iter = 0

    # Kalman update function
    def update(self, dt):

        if self.iter == 1:
            self.X[self.order, 0] = self._get_init_soc()

        self.iter += 1
        # Reverse current direction
        self.measured_current = self.measured_current * (-1)

        self._create_f_and_b(dt)
        # Priori estimate (PREDICT)
        # Predict state:
        # X = F*X + B*u
        try:
            self.X = np.add((self.F * self.X), (self.B * self.measured_current))
        except TypeError:
            if self.X[self.order, 0] == None:
                # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                # in an OCV higher than the measured voltage; default to 1
                self.X[self.order, 0] = 1
                for i in range(1, 10):
                    soc = i / 10
                    if self.measured_voltage < self._OCV(soc):
                        self.X[self.order, 0] = soc
                        break
            else:
                raise TypeError("Initial state vector contains a None that's not at SOC position!")
        # Predict covariance:
        # P = F*P*Ft + Q

        self.P = np.add((1.0 * self.F * self.P * self.F.T), self.Q)
        # Calculate H for posteriori estimate
        self.H[0, self.order] = self._dOCV()  # Set final element of H
        # Calculate measurement error
        self.voltage_error = (self.measured_voltage
                              - self._get_simulated_voltage())

        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        self.Temp = self.P * self.H.T
        self.K = self.Temp * np.reciprocal(np.add((self.H * self.Temp), self.R))
        # Update state:
        # X = X + error*K
        self.X = np.add(self.X, (self.voltage_error * self.K))
        # Update covariance:
        # P = (I - K*H)*P
        posteriori = np.subtract(np.eye(self.order + 1), (self.K * self.H)) * self.P * \
                     (np.subtract(np.eye(self.order + 1), (self.K * self.H)).T) + (self.K * self.R * self.K.T)

        self.P = (posteriori + posteriori.T) / 2
        # Recalculate simulated_voltage for printing purposes
        self.simulated_voltage = self._get_simulated_voltage()
        self.timer += self.p['dt']

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

    def _create_f_and_b(self, dt):
        # State transistion matrix
        self.F = self._discrete_f(np.matrix(np.zeros((self.order + 1, self.order + 1))), dt)
        self.F[self.order, self.order] = 1  # Last element = 1

        # Control matrix
        self.B = np.zeros((self.order + 1, 1))
        for i in range(self.order):
            self.B[i, 0] = dt / self.p['cbranch'][i]
        self.B[self.order, 0] = (-1) * (self.p['eta'] * dt) / self.p['cnom']  # Last element

    def _volts_i_plus_a_t(self, f, dt):
        for i in range(self.order):  # Reference all but last diagonal
            f[i, i] = 1 - dt / (self.p['rbranch'][i] * self.p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t(self, f, dt):
        for i in range(self.order):  # Reference all but last diagonal
            r = self.p['rbranch'][i]
            c = self.p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear(self, f, dt):
        for i in range(self.order):  # Reference all but last diagonal
            rc2 = self.p['rbranch'][i] * self.p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f

    def read_inputs(self, current, voltage):
        self.measured_current = current
        self.measured_voltage = voltage

    def get_simulation_time(self):
        return self.timer

    def get_simulated_soc(self):
        # if self.X[self.order, 0] < 0:
        #     print("SOC is negative")
        return self.X[self.order, 0]

    def get_simulated_voltage(self):
        return self.simulated_voltage

    def get_measured_voltage(self):
        return self.measured_voltage

    def get_mean_percent_diff(self):
        mean = 0
        for i in range(len(self.recent_soc_diffs)):
            mean += self.recent_soc_diffs[i]
        return mean / (i + 1)

    # Evaluate derivative of OCV_SOC curve
    def _dOCV(self):
        x = self.get_simulated_soc()  # Get state of charge
        if x > 1:
            return self._OCV(1) - self.c[-1]
        co_len = len(self.c)
        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * self.c[i] * math.pow(x, co_len - 2 - i)
        return docv

    # Evaluated OCV_SOC curve
    def _OCV(self, soc=None):
        x = self.get_simulated_soc() if soc == None else soc  # Get state of charge
        overage = 0
        if x > 1:
            b = self.c[-1]
            return (self._OCV(1) - b) * x + b
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            ocv += self.c[i] * math.pow(x, co_len - 1 - i)
        return ocv + overage

    # Gets voltage error between predicted and simulated voltage
    def _get_simulated_voltage(self):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        simulated_voltage = self._OCV()
        # if simulated_voltage < 0:
        #     print("simulated voltage is negative: {}".format(simulated_voltage))

        for i in range(self.order):
            simulated_voltage -= self.X[i, 0]
        if self.measured_current > 0:
            simulated_voltage -= self.measured_current * self.p['rpos']
        else:
            simulated_voltage -= self.measured_current * self.p['rneg']

        return simulated_voltage

    def _OCV_init(self, soc=None):
        x = self.get_simulated_soc() if soc == None else soc  # Get state of charge
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            ocv += self.c[i] * math.pow(x, co_len - 1 - i)
        return ocv

    def _get_init_soc(self):
        init_voltage = self.get_measured_voltage()
        temp_soc = 0.00001
        smallest_dif = 99999
        best_soc = 0.000000123
        for i in range(10000):
            temp = abs(self._OCV_init(temp_soc) - init_voltage)
            if temp < smallest_dif:
                smallest_dif = temp
                best_soc = temp_soc
            temp_soc = (i + 1) / 10000

        print('Init SoC guess: {}'.format(best_soc))
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

        self.Q = extend_diagonal(Q, 100000)
        self.P = extend_diagonal(P, 0.001)
        self.X = extend_vector(X, 0.5)
        self.K = np.zeros((order + 2, ))
        self.A = 1
        self.H = np.zeros((1, order + 2))
        for i in range(order):
            self.H[0, i] = -1  # The elements of the vector that represent voltages are -1

        self.Temp = np.zeros((order + 2, order + 2))

        self.iter = 0

    # Kalman update function
    def update(self, dt):

        # self.plot_poly()

        if (self.iter == 0):
            self.X[self.order, 0] = self._get_init_soc()

        self.iter += 1

        order = self.order
        # Reverse current direction
        self.measured_current = self.measured_current * (-1)
        self._soh_f(dt, self.measured_current)
        # Priori estimate (PREDICT)
        # Predict state:

        self.X = self.F * self.X
        for i in range(order):
            self.X[i, 0] += (dt * self.measured_current) / self.p['cbranch'][i]
        self.X[order, 0] -= (dt * self.measured_current) / (self.p['cnom'])

        # We expect the last element of X (soh) to stay the same, so we don't bother changing it.

        # Predict covariance:
        temp = self.P * self.F.T
        # P = F*P*Ft + Q
        self.P = np.add((self.A * self.F * temp), self.Q)
        # Calculate H for posteriori estimate
        self.H[0, order] = self._dOCVdSOC()  # Set final element of H
        self.H[0, order + 1] = self._dOCVdSOH()
        # Calculate measurement error
        self.voltage_error = (self.measured_voltage
                              - self._get_simulated_voltage())

        # if abs(self.voltage_error) <= 0.01:
        #
        #     self.A = 1
        # else:
        if self.measured_current > 0:
            self.A = 0.98
        else:
            self.A = 1.01
        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        self.Temp = self.P * self.H.T
        temp = self.H.T * np.reciprocal(np.add((self.H * self.Temp), self.R))
        self.K = self.P * temp
        # Update state:
        # X = X + error*K
        self.X = np.add(self.X, (self.K * self.voltage_error))
        # Update covariance:
        # P = (I - K*H) * P * (I - K*H)t + K * R * Kt
        temp = self.P * \
               (np.subtract(np.eye(self.order + 2), (self.K * self.H))).T
        posteriori = np.subtract(np.eye(self.order + 2), (self.K * self.H)) * temp + (self.K * self.R * self.K.T)

        self.P = (posteriori + posteriori.T) / 2
        # P = (I - K*H)*P
        # self.P = np.subtract(np.eye(self.order + 2), (self.K * self.H)) * self.P
        # Recalculate simulated_voltage for printing purposes
        self.simulated_voltage = self._get_simulated_voltage()
        self.timer += self.p['dt']


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

    def get_soh(self):
        return self.X[self.order + 1, 0]

    def get_functional_soc(self):
        return self.X[self.order, 0] / self.X[self.order + 1, 0]

    def _soh_f(self, dt, current):
        # For this subclass, we don't actually create b
        order = self.order
        # State transistion matrix
        self.F = self._discrete_f(np.matrix(np.zeros((self.order + 2, self.order + 2))), dt)
        self.F[order, order] = 1  # Last element = 1
        self.F[order + 1, order + 1] = 1

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
    #
    # # Evaluated OCV_SOC curve
    # def _OCV(self, soc=None):
    #     x = self.get_simulated_soc() if soc == None else soc  # Get state of charge
    #     soh = self.X[self.order + 1, 0]
    #     co_len = len(self.c)
    #     ocv = 0
    #     for i in range(co_len):
    #         ocv += self.c[i] * math.pow(x / soh, co_len - 1 - i)
    #     return ocv

    def _dOCVdSOC(self):
        order = self.order
        soc, soh = self.X[order, 0], self.X[order + 1, 0]
        c = self.c
        co_len = len(c)

        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * c[i] * (pow(soc * soh, co_len - 2 - i))
        return docv

    def _dOCVdSOH(self):
        order = self.order
        soc, soh = self.X[order, 0], self.X[order + 1, 0]
        c = self.c
        co_len = len(c)

        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * c[i] * (pow(soc * soh, co_len - 2 - i))
        return docv

    # Evaluated OCV_SOC curve
    def _OCV(self, soc=None):
        x = self.get_simulated_soc() if soc == None else soc  # Get state of charge
        soh = self.X[self.order + 1, 0]
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            try:
                ocv += self.c[i] * math.pow(x * soh, co_len - 1 - i)
            except:
                print("x: {}; soh: {}".format(x, soh))
        return ocv

    def plot_poly(self):
        x = []
        y = []
        soc = 0.00001

        for i in range(6300):
            temp = self._OCV(soc)
            soc = i/6300
            y.append(temp)
            x.append(soc)

        plt.figure(2)
        plt.plot(x, y)
        plt.ylim(2.5, 4.5)
        # plt.show()
        x.pop(0)    #so i can breakpoint


