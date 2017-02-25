import csv
import signal

from classes import ProgressIndicator, KalmanOptimizer, InputFile
import config

MAX_ITERATIONS = 2500
GROUP_SIZE = 7

def float_or_none(string):
    try:
        value = float(string)
    except:
        value = None
    return value

min_rmse        = 9999
max_soc_error   = 9999
best_parameters = None

def finish(signum=None, frame=None):
    raise Exception('\nLowest rmse: {}\tMax error: {}\n'.format(min_rmse, max_soc_error) +
          'Best Parameters:\n\tRpos\t:\t{}\n\tRneg\t:\t{}\n\tRbr\t:\t{}\n\tCbr\t:\t{}\n'
                    .format(best_parameters['rpos'],
                            best_parameters['rneg'],
                            best_parameters['rbranch'],
                            best_parameters['cbranch']))

signal.signal(signal.SIGINT, finish)
signal.signal(signal.SIGTERM, finish)



inputF = InputFile('inputs/note1_input_36.txt')
kalman_config = inputF.adjustConfig(*config.get())[0]



best_parameters = kalman_config[5]
optimizer = KalmanOptimizer(kalman_config, GROUP_SIZE, inputF, initSoc=0.06)


progressIndicator = ProgressIndicator.get_instance()
for i in range(0, MAX_ITERATIONS):
    # If the candidate returned by the optimizer is better than our best, replace our best
    # with the candidate
    progressIndicator.set_tag('{}{}/{}'.format(' ' * (5 - len(str(i))),
                                               i,
                                               MAX_ITERATIONS))
    candRmse, candParams = optimizer.optimize(best_parameters)
    if candRmse < min_rmse:
        min_rmse = candRmse
        best_parameters = candParams

finish()