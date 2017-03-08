'''
    IC curve handler

    Input: voltage   - can be either raw voltage or kalman filter processed voltage
           inputsoc  - from SOC tracker or kalman filter estimated soc

    Description:
            The IC handler takes the input data and processes with different
            algorithms. The purpose of this class is to generate the IC curve
            graph.

'''
# Python Imports
# =============================================================================
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys


class IC_handler:
    def __init__(self, voltage: list, inputsoc: list):
        self.voltage = voltage
        self.inputsoc = inputsoc
        self.QV = []

    def get_IC_volt(self):

        return self.voltage

    def get_IC_QV(self):
        self._dQdV2()
        lenQV = len(self.QV)

        fill = lenQV - len(self.voltage)

        if fill < 0:
            for _ in range(abs(fill)):
                self.QV.append(0)
        elif fill > 0:
            for _ in range(fill):
                self.QV.pop()

        return self.QV

    def get_IC_soc(self):
        return self.inputsoc

    def get_IC_peak(self):
        return max(self.QV)


    def _dQdV(self):
        Q = self.inputsoc
        V = self.voltage

        lenQ = len(Q)

        for i in range(lenQ - 1):
            if i == 0:
                self.QV.append(0)
            else:
                try:
                    dqdv = float((Q[i] - Q[i-1]) / V[i] - V[i-1])
                except ZeroDivisionError:
                    dqdv = 10000000000000000000000000000
                    print("uhhh u effed up cuz this isnt differentiable")
                self.QV.append(dqdv)

    def _get_d(self, c, order, x):
        result = 0

        for i in range(order - 1):
            result += (order - i) * x * c[i]

        return result

    def _dQdV2(self):
        Q = self.inputsoc
        V = self.voltage

        order = 2

        lenQ = len(Q)

        window_width = 500

        y = []
        x = []
        q = []
        endd = Q.index(Q[-1]) - window_width
        for q_ in Q:
            if q_ >= 4900:
                endd = Q.index(q_)
                break

        _V = V
        V = [_V[0]]

        for voltage in _V[1:]:
            if V[-1] < voltage:
                V.append(voltage)
            else:
                V.append(V[-1] + 0.00001)

        for j in range(endd):
            coeffs = np.polyfit(V[j:j+int(window_width)], Q[j:j+int(window_width)], order)
            p = np.poly1d(coeffs)
            y.append(p(V[j + int(window_width / 2)]))
            x.append(V[j + int(window_width / 2)])
            q.append(Q[j + int(window_width / 2)])
            # d = np.polyder(np.poly1d(coeffs))
            # dqdv = d(window_width/2)
            dp = np.polyder(p)
            dqdv = dp(V[j + int(window_width / 2)])
            # dqdv = 2 * coeffs[0] * V[int(window_width/2)] + coeffs[1]
            self.QV.append(1 / dqdv)

        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(x, y, '-', label='{} order poly fit V'.format(order))
        ax.plot(x, q, '-', label='Raw Data')
        ax2 = ax.twinx()
        ax2.plot(x, self.QV, '-r', label='temp')
        ax.legend(loc=0)
        ax.grid()
        ax.set_xlabel("Votlage")
        ax.set_ylabel("Charge A.s")
        ax2.set_ylabel("dQ / dV")
        ax2.set_ylim(min(self.QV), max(self.QV))
        ax.set_ylim(min(y), max(y))
        plt.show()

        return

    def _dQdV3(self):
        Q = self.inputsoc
        V = self.voltage

        order = 3

        coeffs = np.polyfit(V, Q, order)
        x = []
        y = []
        p = np.poly1d(coeffs)
        dp = np.polyder(p)
        for point in V:
            x.append(point)
            y.append(p(point))
            dqdv = dp(point)
            self.QV.append(1 / dqdv)

        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(x, y, '-', label='{} order poly fit V'.format(order))
        ax.plot(x, Q, '-', label='Raw Data')
        ax2 = ax.twinx()
        ax2.plot(x, self.QV, '-r', label='temp')
        ax.legend(loc=0)
        ax.grid()
        ax.set_xlabel("Votlage")
        ax.set_ylabel("Charge A.s")
        ax2.set_ylabel("dQ / dV")
        ax2.set_ylim(min(self.QV), max(self.QV))
        ax.set_ylim(min(y), max(y))
        plt.show()

        # plt.figure(3)
        # plt.plot(x,y)
        # plt.plot(x,self.QV)
        # plt.show()

        return

    def _dQdV4(self):
        Q = self.inputsoc
        V = self.voltage

        order = 2

        lenQ = len(Q)

        window_width = 1000

        y = []
        x = []
        q = []
        time = [j for j in range(len(V))]
        counter = 0
        for j in range(len(V) - window_width):
            if V[j] >= 4.2:
                break
            counter += 1
            coeffs = np.polyfit(Q[j:j + int(window_width)], time[j:j + int(window_width)], order)
            dq = self._get_d(coeffs, order, time[j + int(window_width / 2)])
            coeffs = np.polyfit(V[j:j + int(window_width)], time[j:j + int(window_width)], order)
            dv = self._get_d(coeffs, order, time[j + int(window_width / 2)])
            p = np.poly1d(coeffs)
            y.append(p(V[j + int(window_width / 2)]))
            x.append(V[j + int(window_width / 2)])
            q.append(Q[j + int(window_width / 2)])
            # d = np.polyder(np.poly1d(coeffs))
            # dqdv = d(window_width/2)
            dqdv = dq / dv
            # dqdv = 2 * coeffs[0] * V[int(window_width/2)] + coeffs[1]
            self.QV.append(dqdv)

        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(x, y, 'x', label='{} order poly fit V'.format(order))
        ax.plot(x, q, '+', label='Raw Data')
        ax2 = ax.twinx()
        ax2.plot(x, self.QV, '*', color='green', label='temp')
        ax.legend(loc=0)
        ax.grid()
        ax.set_xlabel("Votlage")
        ax.set_ylabel("Charge A.s")
        ax2.set_ylabel("dQ / dV")
        ax2.set_ylim(min(self.QV), max(self.QV))
        ax.set_ylim(min(y), max(y))

        fig = plt.figure(4)
        uh = fig.add_subplot(111)
        uh.plot(time[:len(q)], q)

        fig = plt.figure(5)
        uh = fig.add_subplot(111)
        uh.plot(time[:len(q)], x)

        plt.show()

        return

    def _dQdV5(self):
        Q = self.inputsoc
        V = self.voltage

        order = 13

        time = [j for j in range(len(V))]

        voltage_cutoff_index = 999

        for voltage in V:
            if voltage >= 4.2:
                voltage_cutoff_index = V.index(voltage)
                break

        coeffs = np.polyfit(time[:voltage_cutoff_index], V[:voltage_cutoff_index], order)

        dV = []
        v = []
        p = np.poly1d(coeffs)
        for i in range(voltage_cutoff_index):
            v.append(p(V[i]))
            dV.append(1 / self._get_d(coeffs, order, V[i]))

        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(time[:voltage_cutoff_index], V[:voltage_cutoff_index], label='actual voltage')
        ax.plot(time[:voltage_cutoff_index], v, label='modeled voltage')
        ax2 = ax.twinx()
        ax2.plot(time[:voltage_cutoff_index], dV, label='dV / dt ^ -1', color='green')
        ax.legend(loc=0)
        ax.grid()

        # plt.figure()
        #
        # plt.plot(time[:voltage_cutoff_index], V[:voltage_cutoff_index], marker='+', label='actual voltage')
        # plt.plot(time[:voltage_cutoff_index], v, marker='*', label='modeled voltage')
        # plt.plot(time[:voltage_cutoff_index], dV, label='dV / dt')
        #
        # plt.legend()

        plt.show()

        return
