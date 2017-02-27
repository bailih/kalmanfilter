'''
    This is script generates model parameters used by an extended Kalman filter for state of
    charge estimation, and optimizes them for variations from real voltage and state of charge
    values.

    It first generates a large number of random model parameters. It then creates a corresponding
    Kalman filter for each set of parameters, and compares the filters' estimations to real values
    to determine performance quality. It takes the set of lowest error filters, and further
    optimizes their parameters through the use of iterative random perturbation.

    NOTE: Extreme optimizer does not currently use InputFile; you must change config.py

    Written by Nelson Schultz for Cadex Electronics, 2016
'''

import csv
import random
import sys
from copy import deepcopy
from math import sqrt

from classes import SoCTracker, KalmanOptimizer, ProgressIndicator
from filters import Kalman

import config

FILE_TO_OPTIMIZE = 'c_by_4_input_001'
RNEG_RANGE = (0.01, 3)
R_RANGE    = (0.000001, 0.5)
C_RANGE    = (100, 20000)

INIT_POINT_COUNT = 5000 # How many random points to test at the start
TOP_CUTOFF       = 10   # How many of the 'best' random points should we use
GROUP_SIZE       = 5    # How many kalman filters for each optimizer per iteration
MAX_ITERATIONS   = 100  # How many iterations the optimizer should go through

def float_or_none(string):
    '''Returns string cast as float if possible, or None if not'''
    try:
        value = float(string)
    except ValueError:
        value = None
    return value

def print_to_csv(*args, file):
    values = '{},' * len(args) + '\n'
    file.write(values.format(*args))


def insert_by_first_val(insertee: list, object: 'iterable'):
    '''
        Inserts object into list such that the list is in ascending order of
        the objects' first values
    '''
    # Try/except avoids checking if the list is empty every time, and is therefore faster
    try:
        if object[0] <= insertee[0][0]:
            insertee.insert(0, object)
            return

        for i in range(1, len(insertee)):
            if object[0] <= insertee[i][0]:
                insertee.insert(i, object)
                return

        insertee.append(object)

    except IndexError:
        # Make sure the list being empty is what caused the error, to avoid catching
        # unforeseen exceptions
        if not insertee:
            insertee.append(object)
        else:
            raise Exception('Index out of bounds')


# Pull the lines of the input file into a list so we have easy access to them at all times
fileLines = []
with open('inputs/{}.txt'.format(FILE_TO_OPTIMIZE), 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for line in reader:
        fileLines.append(line)


# Initialize the set of Kalman filters with pure random model parameters
kalmans = []
base_config = config.get()[0]
for _ in range(INIT_POINT_COUNT):
    config_copy = deepcopy(base_config)

    # Don't bother with rpos, as it doesn't affect the charging performance
    config_copy[5]['rneg'] = random.uniform(*RNEG_RANGE)
    for i in range(config_copy[0]):
        config_copy[5]['rbranch'][i] = random.uniform(*R_RANGE)
        config_copy[5]['cbranch'][i] = random.uniform(*C_RANGE)

    kalmans.append(Kalman(*config_copy))

# Build the list of actual soc values
sys.stdout.write('\rSocTracker')
socTracker = SoCTracker(14, base_config[5]['cnom'], 0.06, 2000)
for line in fileLines:

    voltage, current, dt = map(float, line[:3])
    balance_cell, balance_dir = float_or_none(line[3]), float_or_none(line[4])

    # socTracker expects milliamps and millivolts
    socTracker.update(voltage * 1000, current * 1000, dt, (balance_cell, balance_dir))

# We run the file information through each Kalman filter; we find the voltage and SOC
# root-mean-squared error, add them together and add them to a list, order by lowest
# (ie. best) rmse values.
k_count = 0
performances = []
for kalman in kalmans:
    k_count += 1
    sys.stdout.write('\rRANDOM BATCH:\t{}{}/{}'.format((4 - len(str(k_count))) * ' ',
                                                       k_count,
                                                       INIT_POINT_COUNT))

    soc_err_sum  = 0
    volt_err_sum = 0
    for i in range(len(fileLines)):

        line = fileLines[i]
        voltage, current, dt = map(float, line[:3])

        kalman.read_inputs(current, voltage)
        kalman.update(dt)

        volt_err_sum += abs(kalman.get_simulated_voltage() - kalman.get_measured_voltage())
        soc_err_sum  += abs(kalman.get_simulated_soc() - socTracker.avg_soc[i])

    total_rmse = sqrt(volt_err_sum / (i + 1)) + sqrt(soc_err_sum / (i + 1))
    insert_by_first_val(performances, (total_rmse, kalman.p))

# backup results to file, just in case
with open('best_rmses.csv', 'w') as f:
    print_to_csv('rmse','rneg', 'r1', 'r2', 'r3', 'c1', 'c2', 'c3', file=f)
    for performance in performances:
        rmse, params = performance
        print_to_csv(rmse, params['rneg'], *(params['rbranch'] + params['cbranch']), file=f)


progressIndicator = ProgressIndicator.get_instance()

# Take the TOP_CUTOFF best values, and further optimize each of them, through random
# perturbations.
optimizers = [KalmanOptimizer(base_config, GROUP_SIZE, fileLines, socTracker)
              for _ in range(TOP_CUTOFF)]

best_performances = performances[:TOP_CUTOFF]
for _ in range(MAX_ITERATIONS):
    for i in range(TOP_CUTOFF):
        progressIndicator.set_tag('{}{}/{} {}{}/{}'.format(' ' * (2 - len(str(i))),
                                                           i,
                                                           TOP_CUTOFF,
                                                           ' ' * (4 - len(str(_))),
                                                           _,
                                                           MAX_ITERATIONS))
        bestRmse, bestParams = best_performances[i]
        candRmse, candParams = optimizers[i].optimize(bestParams)
        if candRmse < bestRmse:
            best_performances[i] = (candRmse, candParams)

with open('best_of_the_best.txt', 'w') as f:
    for performance in best_performances:
        f.write('rmse: {} params: {}\n'.format(performance[0], performance[1]))
