import random
from copy import deepcopy
from math import sqrt

import numpy as np

from . import SoCTracker, ProgressIndicator
from filters import Kalman

VELOCITY = 0.1

def float_or_none(string):
    try:
        value = float(string)
    except:
        value = None
    return value

class KalmanOptimizer:

    def __init__(self, config, groupSize, inputFile, socTracker=None, initSoc=0):
        self.config = config
        self.k_wrappers = [OptimizingKalmanWrapper(*deepcopy(config)) for _ in range(groupSize)]
        self.inputFile = inputFile

        # If we aren't passed a pre-ran socTracker, create one ourselves
        if socTracker is None:
            tracker_kwargs = {
                'cell_count' : 14,
                'cell_capacity' : config[5]['cnom'],
                'init_soc' : initSoc,
                'balance_current' : 2000
            }

            # Overwrite defaults with input file values, if they exist
            for key in tracker_kwargs :
                try :
                    tracker_kwargs[key] = float(inputFile[key])
                except KeyError :
                    print('Using default value for key {}'.format(key))

            socTracker = SoCTracker(**tracker_kwargs)
            for line in inputFile:
                voltage, current, dt = map(float, line[:3])
                balance_info = (float_or_none(line[3]), float_or_none(line[4])) if len(line) == 5 \
                    else (None, None)
                socTracker.update(voltage * 1000, current * 1000, dt, balance_info)
        self.socTracker = socTracker

    def optimize(self, p: dict) -> (float, dict):
        '''
            Run Kalman filters with random variations of parameters passed.
            :returns: tuple of lowest rmse and best parameters associated with it
        '''
        inputLines = self.inputFile
        avg_soc = self.socTracker.avg_soc
        progressIndicator = ProgressIndicator.get_instance()
        progressIndicator.init_bar(len(self.k_wrappers), len(self.k_wrappers))

        lowestRmse = None
        bestParams = None
        for wrapper in self.k_wrappers:
            wrapper.start(p, inputLines, avg_soc)
            rmse = sum(wrapper.get_rmses())

            try:
                if rmse < lowestRmse:
                    lowestRmse = rmse
                    bestParams = wrapper.p
            except TypeError:
                # Ensure the type error is caused by lowestRmse, so we don't catch other issues
                if lowestRmse is None:
                    lowestRmse = rmse
                    bestParams = wrapper.p

            progressIndicator.update_bar()

        return(lowestRmse, bestParams)


class OptimizingKalmanWrapper:

    def __init__(self, order, Q, P, R, X, p, c):
        self.order = order
        self.Q = Q
        self.P = P
        self.R = R
        self.X = X

        self.c = c

    def start(self, p: dict, lines: list, avg_socs: list):
        _p = deepcopy(p)
        _p['rpos'] = _p['rpos'] + VELOCITY * _p['rpos'] * random.uniform(-1, 1)
        _p['rneg'] = _p['rneg'] + VELOCITY * _p['rneg'] * random.uniform(-1, 1)

        for i in range(len(_p['rbranch'])):
            _p['rbranch'][i] = _p['rbranch'][i] + VELOCITY * _p['rbranch'][i] * random.uniform(-1, 1)
            _p['cbranch'][i] = _p['cbranch'][i] + VELOCITY * _p['cbranch'][i] * random.uniform(-1, 1)

        self.p = _p

        kalman = Kalman(self.order, self.Q, self.P, self.R, self.X, _p, self.c)

        n = 0
        sum_voltage_error = 0
        sum_soc_error = 0
        max_error = 0

        for i in range(len(lines)) :
            line, avg_soc = lines[i], avg_socs[i]
            kalman.read_inputs(float(line[1]), float(line[0]))
            kalman.update(float(line[2]))

            n += 1
            sum_voltage_error += np.absolute(kalman.get_simulated_voltage()
                                             - kalman.get_measured_voltage())

            soc_error = abs(kalman.get_simulated_soc() - avg_soc)
            sum_soc_error += soc_error

            if soc_error > max_error :
                max_error = soc_error

        self.voltRmse = sqrt(sum_voltage_error / n)
        self.socRmse  = sqrt(sum_soc_error / n)
        self.maxError = max_error

    def get_rmses(self):
        return (self.voltRmse, self.socRmse)

    def get_max_error(self):
        return self.maxError