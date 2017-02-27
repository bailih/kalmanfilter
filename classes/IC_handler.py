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
        self._dQdV3()
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

        window_width = 6000


        for j in range(len(V) - window_width):
            coeffs = np.polyfit(V[j:j+int(window_width)], Q[j:j+int(window_width)], order)
            # d = np.polyder(np.poly1d(coeffs))
            # dqdv = d(window_width/2)
            dqdv = self._get_d(coeffs, order, V[j + int(window_width/2)])
            # dqdv = 2 * coeffs[0] * V[int(window_width/2)] + coeffs[1]
            self.QV.append(dqdv)

        plt.plot()

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
            self.QV.append(dqdv)

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
