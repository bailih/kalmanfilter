'''
    IC curve handler

    Input: voltage   - can be either raw voltage or kalman filter processed voltage
           inputsoc  - from SOC tracker or kalman filter estimated soc

    Description:
            The IC handler takes the input data and processes with different
            algorithms. The purpose of this class is to generate the IC curve
            graph.

'''


class IC_handler:
    def __init__(self, voltage: list, inputsoc: list):
        self.voltage = voltage
        self.inputsoc = inputsoc
        self.QV = []

    def get_IC_volt(self):

        return self.voltage

    def get_IC_QV(self):
        self._dQdV()
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
                dqdv = float((Q[i] - Q[i - 1]) / V[i])
                self.QV.append(dqdv)
