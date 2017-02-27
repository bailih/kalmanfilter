'''
    This class is to hold a list of all the different matrix's for the Kalman Filter

    matrix must hold data for voltage, soc, init soc, and init soh
'''

import numpy as np


class matrixHolder:

    '''
    0: for 3.0-3.5 Vinit
    1: 3.5-3.8 Vinit, low SoH
    2: 3.5-3.8 Vinit, high SoH
    3: 3.8~3.9 Vinit 50 < SoH < 70
    4: 3.8~3.9 Vinit 90 < SoH
    5: 3.8~3.9 Vinit SoH < 50
    6: 3.8~3.9 Vinit 70 < SoH < 90
    '''

    def __init__(self, capacity, inisoc, iniV):

        case0_X = np.matrix([[0],   #case 0
                           [0],
                           [0],
                           [capacity]])

        case1_X = np.matrix([[0],   #case 1
                           [0],
                           [0],
                           [capacity*0.5]])

        case2_X = np.matrix([[0],   #case 2
                           [0],
                           [capacity/2],
                           [capacity]])

        case3_X = np.matrix([[0],   #case 3
                           [0],
                           [capacity * 0.3],
                           [capacity * 0.5]])

        case4_X = np.matrix([[0],   #case 4
                           [0],
                           [capacity * 0.8],
                           [capacity]])

        case5_X = np.matrix([[0],   #case 5
                           [0],
                           [0],
                           [capacity / 2]])

        case6_X = np.matrix([[0],   #case 6
                           [0],
                           [capacity * inisoc],
                           [capacity]])

        case0_Q = np.matrix([[0.00000001, 0, 0, 0],   #case 0
                               [0, 0.0001, 0, 0],
                               [0, 0, 0.0001, 0],
                               [0, 0, 0, 100000]])

        case1_Q = np.matrix([[0.00001, 0, 0, 0],   #case 1
                               [0, 0.00001, 0, 0],
                               [0, 0, 1000, 0],
                               [0, 0, 0, 100000]])

        case2_Q = np.matrix([[0.00000001, 0, 0, 0],   #case2
                               [0, 0.0001, 0, 0],
                               [0, 0, 1000, 0],
                               [0, 0, 0, 100]])

        case3_Q = np.matrix([[0.00000001, 0, 0, 0],   #case 3
                               [0, 0.0001, 0, 0],
                               [0, 0, 1000, 0],
                               [0, 0, 0, 500]])

        case4_Q = np.matrix([[0.00000001, 0, 0, 0],   #case 4
                               [0, 0.0001, 0, 0],
                               [0, 0, 1000, 0],
                               [0, 0, 0, 100]])

        case5_Q = np.matrix([[0.00000001, 0, 0, 0],   #case 5
                               [0, 0.0001, 0, 0],
                               [0, 0, 0.001, 0],
                               [0, 0, 0, 0.001]])

        case6_Q = np.matrix([[0.00000001, 0, 0, 0],   #case 5
                               [0, 0.0001, 0, 0],
                               [0, 0, 80, 0],
                               [0, 0, 0, 1000]])

        case0_P = np.matrix([[0.00001, 0, 0, 0], #case 0
                               [0, 0.0001, 0, 0],
                               [0, 0, 100, 0],
                               [0, 0, 0, 10000]])

        case1_P = np.matrix([[0.001, 0, 0, 0], #case 1
                               [0, 0.0001, 0, 0],
                               [0, 0, 100, 0],
                               [0, 0, 0, 0.001]])

        case2_P = np.matrix([[0.000001, 0, 0, 0], #case 2
                               [0, 0.0001, 0, 0],
                               [0, 0, 100, 0],
                               [0, 0, 0, 10000]])

        case3_P = np.matrix([[0.000001, 0, 0, 0], #case 3
                               [0, 0.0001, 0, 0],
                               [0, 0, 100, 0],
                               [0, 0, 0, 100]])

        case4_P = np.matrix([[0.000001, 0, 0, 0], #case 4
                           [0, 0.0001, 0, 0],
                           [0, 0, 1000, 0],
                           [0, 0, 0, 1000000]])

        case5_P = np.matrix([[0.000001, 0, 0, 0], #case 5
                           [0, 0.0001, 0, 0],
                           [0, 0, 99, 0],
                           [0, 0, 0, 999]])

        case6_P = np.matrix([[0.000001, 0, 0, 0], #case 5
                           [0, 0.0001, 0, 0],
                           [0, 0, 999, 0],
                           [0, 0, 0, 999]])
        self.X = [

        ]

        self.Q = [

        ]

        self.P = [

        ]

        self.R = [

        ]

        self.case = []
        '''
        0: 3.0-3.5 Vinit
        1: 3.5-3.8 Vinit, low SoH
        2: 3.5-3.8 Vinit, high SoH
        3: 3.8~3.9 Vinit 50 < SoH < 90
        4: 3.8~3.9 Vinit 90 < SoH
        5: 3.8~3.9 Vinit SoH < 50
        '''
        if True:
            self.X.append(case0_X)
            self.Q.append(case0_Q)
            self.P.append(case0_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(0)

            self.X.append(case1_X)
            self.Q.append(case1_Q)
            self.P.append(case1_P)
            self.R.append(np.matrix([10]))
            self.case.append(1)

            self.X.append(case2_X)
            self.Q.append(case2_Q)
            self.P.append(case2_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(2)

            self.X.append(case3_X)
            self.Q.append(case3_Q)
            self.P.append(case3_P)
            self.R.append(np.matrix([10]))
            self.case.append(3)

            self.X.append(case4_X)
            self.Q.append(case4_Q)
            self.P.append(case4_P)
            self.R.append(np.matrix([0.001]))
            self.case.append(4)

            self.X.append(case5_X)
            self.Q.append(case5_Q)
            self.P.append(case5_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(5)

            self.X.append(case6_X)
            self.Q.append(case6_Q)
            self.P.append(case6_P)
            self.R.append(np.matrix([10]))
            self.case.append(6)

            length = len(self.X)
            if any(len(lst) != length for lst in [self.Q, self.P, self.R]):
                raise RuntimeError(
                    "One or more of the EKF matrix parameters are missing! check matrixHolder.py -> __init__")
        elif iniV < 3.5:
            self.X.append(case0_X)
            self.Q.append(case0_Q)
            self.P.append(case0_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(0)
            length = len(self.X)
            if any(len(lst) != length for lst in [self.Q, self.P, self.R]):
                raise RuntimeError(
                    "One or more of the EKF matrix parameters are missing! check matrixHolder.py -> __init__")


        elif 3.8 >= iniV >= 3.5:

            self.X.append(case1_X)
            self.Q.append(case1_Q)
            self.P.append(case1_P)
            self.R.append(np.matrix([10]))
            self.case.append(1)

            self.X.append(case2_X)
            self.Q.append(case2_Q)
            self.P.append(case2_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(2)
            length = len(self.X)
            if any(len(lst) != length for lst in [self.Q, self.P, self.R]):
                raise RuntimeError(
                    "One or more of the EKF matrix parameters are missing! check matrixHolder.py -> __init__")

        elif iniV > 3.8:
            self.X.append(case3_X)
            self.Q.append(case3_Q)
            self.P.append(case3_P)
            self.R.append(np.matrix([10]))
            self.case.append(3)

            self.X.append(case4_X)
            self.Q.append(case4_Q)
            self.P.append(case4_P)
            self.R.append(np.matrix([0.001]))
            self.case.append(4)

            self.X.append(case5_X)
            self.Q.append(case5_Q)
            self.P.append(case5_P)
            self.R.append(np.matrix([0.0001]))
            self.case.append(5)

            self.X.append(case6_X)
            self.Q.append(case6_Q)
            self.P.append(case6_P)
            self.R.append(np.matrix([10]))
            self.case.append(6)

            length = len(self.X)
            if any(len(lst) != length for lst in [self.Q, self.P, self.R]):
                raise RuntimeError(
                    "One or more of the EKF matrix parameters are missing! check matrixHolder.py -> __init__")

        else:
            print('require more matrix setting')



    def getP(self):
        return self.P

    def getR(self):
        return self.R

    def getQ(self):
        return self.Q

    def getX(self):
        return self.X

    def getcase(self):
        return self.case

