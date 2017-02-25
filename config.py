"""
    Config File for Kalman Filter SoC Estimator
    
    Contains 1st, 2nd and 3rd order specific parameters, as well
    as independent parameters. Also contains a 'get' function that
    returns a list of the appropriate parameters.

"""

from copy import deepcopy

import numpy as np

ORDER = 1 # Order of equivalent circuit model

# First Order Parameters
# =============================================================================
o1 = dict(
    # Process covariance (Q)
    # Q = np.mat([[0.000001, 0],
    #             [0, 0.0000001]]),
    # # Estimation covariance (P)
    # P = np.mat([[0.000001, 0],
    #             [0, 0.01]]),

    # Q=np.mat([[0.005, 0],
    #           [0, 0.0000001]]),
    # # Estimation covariance (P)
    # P=np.mat([[0.000001, 0],
    #           [0, 0.01]]),
    #
    # # Measurement covariance (R)
    # R = 0.00008,
    Q=np.mat([[0.00001, 0],
              [0, 0.0001]]),
    # Estimation covariance (P)
    P=np.mat([[0.0001, 0],
              [0, 0.0001]]),

    # Measurement covariance (R)
    R=0.001,
    # State information (X)
    X=np.mat([[0],
              [0]]),

    rpos=0.0619247,
    rneg     = 0.09832708333311384,
    R_branch = [1.9143590268846616e-06],
    C_branch = [593.0601435936358],

    psi=5,
    gamma=0.5
)


# Second Order Parameters
# =============================================================================
o2 = dict(
    # Process covariance (Q)
    Q = np.mat([[0.000001, 0, 0], 
                [0, 0.000001, 0], 
                [0, 0, 0.000001]]),
    # Estimation covariance (P)
    P = np.mat([[0.1, 0, 0],
                [0, 0.1, 0],
                [0, 0, 100]]),
    # Measurement covariance (R)
    R = 5,
    # State information (X)
    X = np.mat([[0], 
                [0],
                [0.50]]),

    rpos = 0.0688898,
    rneg     = 0.0029966717107376955,
    R_branch = [0.01802230842925986, 0.08370486406120509],
    C_branch = [5145.494524174826, 719.8220950763777],

    psi   = 0.01,
    gamma = 0.5
)


# Third Order Parameters
# =============================================================================
o3 = dict(
    # Process covariance (Q)
    Q = np.mat([[0.00001, 0, 0, 0],
                [0, 0.00001, 0, 0],
                [0, 0, 0.00001, 0],
                [0, 0, 0, 0.00000001]]),
    # Estimation covariance (P)
    P = np.mat([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 100]]),
    # Measurement covariance (R)
    R = 8,
    # State information (X)
    X = np.mat([[0],
                [0],
                [0],
                [0.50]]),

    rpos = 0.0619369,
    rneg     = 0.03606348881123421,
    R_branch = [0.00044108508894862095, 0.020342613271534582, 0.032071047769841064],
    C_branch = [188.39473796483145, 13015.801680069379, 466.862686312624],

    psi   = 0.01,
    gamma = 0.01
)


# Independant Parameters
# =============================================================================
# Equivalent circuit
dt = 1
cnom = 4800 * 3.6 * 0.72
eta = 1.0
ASOC = 0.016 # Actual SoC - Used as reference for v_sim (SHOULD BE 0.907)
# Tourmaline Coefficients(p1, p2, ... p10, p11)
c = (-33951.366832094,
    184587.18448293468,
    -426275.99729642476,
    545339.423925515,
    -422650.6366402695,
    203921.42413937073,
    -60691.47015957726,
    10736.447069785818,
    -1068.3369113988201,
    61.69682265017936,
    50.51950639830661,
    )


                 
# Functions
# =============================================================================
def get():
    '''
    Returns the configs required by both the ekf and ukf, and the svsf, based on the specified order

    :return: kalman_config: [order, process noise, initial covariance, measurement noise,
                             initial state, model parameters, polynomial coefficients]

             svsf_config:   [order, initial state, model parameters, polynomial coefficients,
                             svsf convergence rate, smoothing boundary layer width]
    '''
    kalman_config = [ORDER]
    svsf_config   = [ORDER]

    o = {
        1 : o1,
        2 : o2,
        3 : o3
    }.get(ORDER, None)

    kalman_config.append(o['Q'])

    # We need deepcopies of anything the filters will modify, to avoid filters affecting each other
    kalman_config.append(deepcopy(o['P']))
    kalman_config.append(o['R'])

    for config in (kalman_config, svsf_config):
        config.append(deepcopy(o['X']))
        config.append({
            'dt'      : dt,
            'rpos'    : o['rpos'],
            'rneg'    : o['rneg'],
            'cnom'    : cnom,
            'rbranch' : o['R_branch'],
            'cbranch' : o['C_branch'],
            'eta'     : eta
        })
        config.append(c)

    svsf_config.append(o['psi'])
    svsf_config.append(o['gamma'])

    return kalman_config, svsf_config