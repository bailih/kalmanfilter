"""
    Header File for Kalman Filter SoC Estimator

    Contains Kalman Filter class and methods for 1st,
    2nd and 3rd order circuits.

"""

# Python Imports
# =============================================================================
from math import pow
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

    def __init__(self, order, Q, P, R, X, p, c) :
        self.order = order  # Order of equivalent circuit
        self.Q = Q  # Process covariance  (Q)
        self.P = P  # Estimation covariance (P)
        self.R = R  # Measurement covariance (R)
        self.X = X  # State information (X)
        self.p = p  # Equivalent circuit paramaters (p)
        self.c = c  # Polynomial coefficients

        self._discrete_f = {
            0 : self._volts_i_plus_a_t,
            1 : self._volts_inv_i_minus_a_t,
            2 : self._volts_bilinear
        }.get(DISCRETE_STRAT)

        # Kalman gain matrix
        self.K = np.zeros((order + 1, 1))

        # Temporary matrix (for calculation)
        self.Temp = np.zeros((order + 1, order + 1))
        # Useful scalar values
        self.simulated_voltage = 0
        self.voltage_error = 0
        self.soc_error = 0

        self.recent_socs = []
        self.recent_soc_diffs = []

        self.iter = 0


    # Kalman update function
    def update(self, dt, voltage, current):

        if self.iter == 0:
            self.X[self.order, 0] = self._get_init_soc(voltage)

        self.iter += 1

        # Reverse current direction
        current *= -1

        f, b = self._create_f_and_b(dt)
        # Priori estimate (PREDICT)
        # Predict state:
        # X = F*X + B*u
        order = self.order
        initial_x = self.X
        try :
            predicted_x = np.add(f * initial_x, b * current)
        except TypeError:
            if initial_x[order, 0] == None:
                # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                # in an OCV higher than the measured voltage; default to 1
                initial_x[order, 0] = 1
                for i in range(1, 10):
                    soc = i / 10
                    if voltage < self._OCV(soc):
                        initial_x[order, 0] = soc
                        break
                predicted_x = np.add(f * initial_x, b * current)
            else:
                raise TypeError("")

        # Predict covariance:
        # P = F*P*Ft + Q
        predicted_p = np.add((f * self.P * f.T), self.Q)

        # Calculate H for posteriori estimate
        h = np.matrix(np.zeros((1, order + 1)))
        for i in range(order):
            h[0, i] = -1
        h[0, order] = self._dOCV(predicted_x[order, 0], self.p['cnom'])  # Set final element of H

        # Calculate measurement error
        self.voltage_error = (voltage - self._get_simulated_voltage(current))

        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        temp = predicted_p * h.T
        k = temp * np.reciprocal(np.add((h * temp), self.R))

        # Update state:
        # X = X + error*K
        new_x = np.add(predicted_x, (self.voltage_error * k))

        # Update covariance:
        # P = (I - K*H)*P
        new_p = np.subtract(np.eye(order + 1), (k * h)) * predicted_p

        self.X = new_x
        self.P = new_p

        # Recalculate simulated_voltage for printing purposes
        self.simulated_voltage = self._get_simulated_voltage(current)

        sim_soc = self.get_soc()
        self.recent_socs.append(sim_soc)
        if len(self.recent_socs) > SOC_HISTORY_COUNT :
            self.recent_socs.pop(0)
        mean = 0
        for i in range(len(self.recent_socs)) :
            mean += self.recent_socs[i]
        mean /= (i + 1)

        self.recent_soc_diffs.append(abs(((mean - self.get_soc()) / mean) * 100))
        if len(self.recent_soc_diffs) > PD_HISTORY_COUNT :
            self.recent_soc_diffs.pop(0)


    def _create_f_and_b(self, dt):
        p = self.p
        order = self.order
        # State transistion matrix
        f = self._discrete_f(np.matrix(np.zeros((order + 1, order + 1))), dt)
        f[order, order] = 1  # Last element = 1

        # Control matrix
        b = np.zeros((order + 1, 1))
        for i in range(self.order):
            b[i, 0] = dt / p['cbranch'][i]
        b[order, 0] = (-1) * (p['eta'] * dt)  # Last element

        return f, b

    def _volts_i_plus_a_t(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            f[i, i] = 1 - dt / (p['rbranch'][i] * p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f


    def get_amp_seconds(self):
        return self.X[self.order, 0]

    def get_soc(self):
        return self.X[self.order, 0]


    def get_estimated_voltage(self):
        return self.simulated_voltage


    def get_mean_percent_diff(self):
        mean = 0
        for i in range(len(self.recent_soc_diffs)):
            mean += self.recent_soc_diffs[i]
        return mean / (i + 1)


    # Evaluate derivative of OCV_SOC curve
    def _dOCV(self, charge, capacity):
        x = self.get_soc()/capacity  # Get state of charge, convert to %
        if x > 1:
            return self._OCV(1) - self.c[-1]
        co_len = len(self.c)
        docv = 0
        for i in range(co_len - 1):
            docv += (co_len - 1 - i) * self.c[i] * pow(x, co_len - 2 - i)# / pow(capacity, co_len - 1 - i)
        return docv


    # Evaluated OCV_SOC curve
    def _OCV(self, soc=None):
        x = self.get_soc()/self.p['cnom'] if soc == None else soc  # Get state of charge, convert to %
        if x > 1:
            b = self.c[-1]
            return (self._OCV(1) - b) * x + b
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            ocv += self.c[i] * pow(x, co_len - 1 - i)
        return ocv

    # Gets voltage error between predicted and simulated voltage
    def _get_simulated_voltage(self, current):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        simulated_voltage = self._OCV()
        for i in range(self.order):
            simulated_voltage -= self.X[i, 0]
        if current > 0:
            simulated_voltage -= current * self.p['rpos']
        else:
            simulated_voltage -= current * self.p['rneg']
        return simulated_voltage

    def get_simulated_voltage(self):
        return self.simulated_voltage

    def get_simulated_soc(self):
        return self.X[self.order, 0]

    def _OCV_init(self, soc=None):
        x = self.get_simulated_soc() if soc == None else soc  # Get state of charge
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            ocv += self.c[i] * pow(x, co_len - 1 - i)
        return ocv

    def _get_init_soc(self, voltage):
        init_voltage = voltage
        temp_soc = 0.00001
        smallest_dif = 99999
        best_soc = 0.000000123
        for i in range(10000):
            temp = abs(self._OCV_init(temp_soc) - init_voltage)
            if temp < smallest_dif:
                smallest_dif = temp
                best_soc = temp_soc * self.p['cnom']    #soc (%) * capacity
            temp_soc = (i + 1) / 10000

        print('Init SoC guess: {} ({})'.format(best_soc, best_soc / self.p['cnom']))
        return best_soc

class KalmanSoC(Kalman):
    """ Class to hold kalman state data and methods """

    def __init__(self, order, Q, P, R, X, p, c) :
        super(KalmanSoC, self).__init__(order, Q, P, R, X, p, c)
        self.order = order  # Order of equivalent circuit

        self.R = R  # Measurement covariance (R)
        self.p = p  # Equivalent circuit paramaters (p)
        self.c = c  # Polynomial coefficients

        self.Q = extend_diagonal(Q, 100)
        self.P = extend_diagonal(P, 999999999)
        self.X = extend_vector(X, 0)

        self.Q = extend_diagonal(self.Q, 50000)
        self.P = extend_diagonal(self.P, 9999)
        self.X = extend_vector(self.X, self.p['cnom'])

        self.K = np.zeros((order + 3, 1))
        self.A = 1
        self.H = np.zeros((1, order + 3))
        for i in range(order):
            self.H[0, i] = -1  # The elements of the vector that represent voltages are -1

        self.Temp = np.zeros((order + 3, order + 3))


        #

        self._discrete_f = {
            0 : self._volts_i_plus_a_t,
            1 : self._volts_inv_i_minus_a_t,
            2 : self._volts_bilinear
        }.get(DISCRETE_STRAT)

        # Kalman gain matrix
        # self.K = np.zeros((order + 1, 1))

        # Temporary matrix (for calculation)
        # self.Temp = np.zeros((order + 1, order + 1))

        # Useful scalar values
        self.simulated_voltage = 0
        self.voltage_error = 0
        self.soc_error = 0

        self.recent_socs = []
        self.recent_soc_diffs = []

        self.iter = 0

        self.conver_step = 0

    # Kalman update function
    def update(self, dt, voltage, current):
        order = self.order



        # if self.iter == 1:
            # if abs(current) < 0.02:
            #     self.X[self.order + 2, 0] = self._get_init_soh(voltage)
            # else:
            #     self.X[self.order + 2, 0] = self._get_init_soh(voltage - 0.2)
            # self.X[order + 2, 0] = 0
        self.iter += 1

        # Reverse current direction
        current *= -1

        f, b = self._create_f_and_b(dt)

        # Priori estimate (PREDICT)
        # Predict state:
        # X = F*X + B*u



        initial_x = self.X
        try:
            predicted_x = np.add(f * initial_x, b * current)

        except TypeError:
            if initial_x[order, 0] == None:
                # Set the initial guess of soc to the smallest multiple of 0.1 that would result
                # in an OCV higher than the measured voltage; default to 1
                initial_x[order, 0] = 1
                for i in range(1, 10):
                    soc = i / 10
                    if voltage < self._OCV(soc):
                        initial_x[order, 0] = soc * self.p['cnom']
                        break
                predicted_x = np.add(f * initial_x, b * current)
            else:
                raise TypeError("")

        # self.A = 0.98

        # Predict covariance:
        # P = F*P*Ft + Q
        predicted_p = np.add((self.A * f * self.P * f.T), self.Q)

        # Calculate H for posteriori estimate
        h = np.matrix(np.zeros((1, order + 3)))
        for i in range(order):
            h[0, i] = -1
        h[0, order] = self._dOCVdIN()  # Set final element of H
        h[0, order + 1] = self._dOCVdINISOC()
        h[0, order + 2] = self._dOCVdSOH()

        # Calculate measurement error
        self.voltage_error = (voltage - self._get_simulated_voltage(current))

        if abs(self.voltage_error) < 0.1 and self.conver_step == 0:

            # self.Q[order+1, order+1] = 0.00001
            # self.P[order+1, order+1] = 0.0001
            # self.Q[order + 2, order + 2] = 1000000
            # self.P[order + 2, order + 2] = 10000000
            self.conver_step = self.iter
            print("converge step: {}".format(self.conver_step))

        # Posteriori estimate (UPDATE)
        # Update kalman gain:
        # K = P*Ht*Inverse(H*P*Ht+R)
        temp = predicted_p * h.T
        k = temp * np.reciprocal(np.add((h * temp), self.R))


        if predicted_x[self.order + 1, 0] < 0:
            predicted_x[self.order + 1, 0] = 0.0001

        if predicted_x[self.order + 1, 0] > 9000:
            predicted_x[self.order + 1, 0] = 9000

        if predicted_x[self.order + 2, 0] > 9000:
            predicted_x[self.order + 2, 0] = 9000

        if predicted_x[self.order + 2, 0] < 0:
            predicted_x[self.order + 2, 0] = 0.0001
        # Update state:
        # X = X + error*K
        new_x = np.add(predicted_x, (k * self.voltage_error))

        # Update covariance:
        # P = (I - K*H)*P*(I - K*H).T + K*R*K.T
        new_p = np.subtract(np.eye(order + 3), (k * h)) * predicted_p * np.subtract(np.eye(order + 3), (k * h)).T
        new_p += k * self.R * k.T
        new_p = (new_p + new_p.T) / 2


        self.X = new_x
        self.P = new_p
        print(self.X)


        if self.X[self.order + 1, 0] < 0:
            self.X[self.order + 1, 0] = 0.0001

        if self.X[self.order + 1, 0] > 9000:
            self.X[self.order + 1, 0] = 9000

        if self.X[self.order + 2, 0] > 9000:
            self.X[self.order + 2, 0] = 9000

        if self.X[self.order + 2, 0] < 0:
            self.X[self.order + 2, 0] = 0.0001


        # Recalculate simulated_voltage for printing purposes
        self.simulated_voltage = self._get_simulated_voltage(current)

        sim_soc = self.get_inputsoc()
        self.recent_socs.append(sim_soc)
        if len(self.recent_socs) > SOC_HISTORY_COUNT :
            self.recent_socs.pop(0)
        mean = 0
        for i in range(len(self.recent_socs)) :
            mean += self.recent_socs[i]
        mean = mean / (i + 1)

        self.recent_soc_diffs.append(abs(((mean - self.get_inputsoc()) / mean) * 100))
        if len(self.recent_soc_diffs) > PD_HISTORY_COUNT :
            self.recent_soc_diffs.pop(0)


    def _create_f_and_b(self, dt):
        p = self.p
        order = self.order
        # State transistion matrix
        f = self._discrete_f(np.matrix(np.zeros((order + 3, order + 3))), dt)
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

    def _volts_i_plus_a_t(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            f[i, i] = 1 - dt / (p['rbranch'][i] * p['cbranch'][i])

        return f

    def _volts_inv_i_minus_a_t(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            r = p['rbranch'][i]
            c = p['cbranch'][i]
            f[i, i] = (r * c) / (r * c + dt)

        return f

    def _volts_bilinear(self, f, dt):
        p = self.p
        for i in range(self.order):     # Reference all but last diagonal
            rc2 = p['rbranch'][i] * p['cbranch'][i] * 2
            f[i, i] = (rc2 - dt) / (rc2 + dt)

        return f


    def get_amp_seconds(self):
        return self.X[self.order, 0]

    def get_inputsoc(self):
        return self.X[self.order, 0]

    def get_inisoc(self):
        return self.X[self.order+1, 0]

    def get_soh(self):
        return self.X[self.order + 2, 0]

    def get_estimated_voltage(self):
        return self.simulated_voltage


    def get_mean_percent_diff(self):
        mean = 0
        for i in range(len(self.recent_soc_diffs)):
            mean += self.recent_soc_diffs[i]
        return mean / (i + 1)


    # Evaluate derivative of OCV_SOC curve
    def _dOCVdIN(self):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        ini = self.get_inisoc()
        input = self.get_inputsoc()
        soh = self.get_soh()
        # if soh < 0:
        #     soh = 0.000000001
        # elif soh > 1:
        #     soh = 1
        # if x > 1:
        #     return self._OCV(1) - self.c[-1]
        co_len = len(self.c)
        docv = 0
        for i in range(co_len - 1):
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) + (co_len - 2 - i))
            docv += (co_len - 1 - i) * self.c[i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (co_len - 1 - i))
        return docv

    # Evaluate derivative of OCV_SOC curve
    def _dOCVdINISOC(self):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom']  # Get state of charge, convert to %
        ini = self.get_inisoc()
        input = self.get_inputsoc()
        soh = self.get_soh()
        # if soh < 0:
        #     soh = 0.000000001
        # elif soh > 1:
        #     soh = 1
        # if x > 1:
        #     return self._OCV(1) - self.c[-1]
        co_len = len(self.c)
        docv = 0
        for i in range(co_len - 1):
            # print ("i : {}, {}".format(i, self.c[i]))
            # docv += (co_len - 1 - i) * self.c[i] * pow((ini + input * soh), co_len - 2 - i) / pow(self.p['cnom'], co_len - 1 - i) / pow(soh, (co_len - 1 - i) * 2)
            docv += (co_len - 1 - i) * self.c[i] * pow((ini + input), (co_len - 2 - i)) / pow(soh, (co_len - 1 - i))
        return docv

        # Evaluate derivative of OCV_SOC curve

    def _dOCVdSOH(self):
        # x = (self.get_inputsoc() + self.get_inisoc()) / self.p['cnom']  # Get state of charge, convert to %
        ini = self.get_inisoc()
        input = self.get_inputsoc()
        soh = self.get_soh()
        # if soh < 0:
        #     soh = 0.00000001
        # elif soh > 1:
        #     soh = 1
        # if x > 1:
        #     return self._OCV(1) - self.c[-1]
        co_len = len(self.c)
        docv = 0
        for i in range(co_len - 1):
            # docv -= (co_len - i) * self.c[i] * pow((input * soh) + ini, co_len - 2 - i) * ((input * soh) + ini * 2) / pow(self.p['cnom'], co_len - i) / pow(soh, (co_len - 1 - i)*2 + 1)
            docv -= (co_len - 1 - i) * self.c[i] * pow((ini + input), (co_len - 1 - i)) / pow(soh, (co_len - i))
        return docv

    # Evaluated OCV_SOC curve
    def _OCV(self, soc=None):
        # x = (self.get_inputsoc() + self.get_inisoc())/self.p['cnom'] if soc == None else soc  # Get state of charge, convert to %
        ini = self.get_inisoc()
        input = self.get_inputsoc()
        soh = self.get_soh()
        # if soh < 0:
        #     soh = 0.00000001
        # elif soh > 1:
        #     soh = 1
        # while x > 1:
        #     self.X[self.order + 1, 0] -= 1000
        #     x = (self.get_inputsoc() + self.get_inisoc()) / self.p['cnom']
        #     b = self.c[-1]
        #     return (self._OCV(1) - b) * x + b
        co_len = len(self.c)
        ocv = 0
        for i in range(co_len):
            # ocv += self.c[i] * pow((ini/self.p['cnom']/soh/soh + input/self.p['cnom']/soh), co_len - 1 - i)
            ocv += self.c[i] * pow((ini / soh + input / soh), co_len - 1 - i)
        return ocv

    # Gets voltage error between predicted and simulated voltage
    def _get_simulated_voltage(self, current):
        # Simulated voltage = ocv(soc) - v1 - (v2 - v3)- i*R0
        simulated_voltage = self._OCV()
        for i in range(self.order):
            simulated_voltage -= self.X[i, 0]
        if current > 0:
            simulated_voltage -= current * self.p['rpos']
        else:
            simulated_voltage -= current * self.p['rneg']
        return simulated_voltage

    def get_simulated_voltage(self):
        return self.simulated_voltage

    def get_simulated_soc(self):
        return self.X[self.order, 0]

    def _get_init_soh(self, voltage):
        init_voltage = voltage
        smallest_dif = 0.1
        best_soh = 9000
        ini = self.get_inisoc()
        input = self.get_inputsoc()


        for i in range(100):
            soh = i * 90 + 0.01
            co_len = len(self.c)
            ocv = 0
            for j in range(co_len):
                ocv += self.c[j] * pow((ini / soh + input / soh), co_len - 1 - j)

            temp = (ocv - init_voltage)
            if temp < smallest_dif and temp > 0:
                best_soh = soh   #soc (%) * capacity
                break

        print('Init SoC_input guess: {}'.format(input))
        print('Init SoC_ini guess: {}'.format(ini))
        print('Init SoH guess: {}'.format(best_soh))
        return best_soh #for testing, return 0




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

        #predict covariance
        predicted_p = np.add((f * self.P * f.T), self.Q)

        #calculate h
        h = np.matrix(np.zeros((1, order + 2)))
        for i in range(order):
            h[0, i] = -1
        h[0, order] = self._dOCVdSOC(predicted_x, self.p['cnom'])
        h[0, order + 1] = self._dOCVdSOH(predicted_x, self.p['cnom'])

        #get voltage error
        self.voltage_error = (voltage - self._get_simulated_voltage(current))

        #will memory fade the filter
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
            docv += (co_len - 1 - i) * c[i] * (pow(soc, co_len - 1 - i))/ pow(soh, co_len - i)
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
        x = self.get_simulated_soc()/self.p['cnom'] if soc == None else soc  # Get state of charge
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
