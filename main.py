"""
    Kalman Filter SoC Estimator Program

    This program reads an input file of voltage 
    and current measurements. It performs a Kalman
    update on each line of the data file and prints
    the results to a graph. 

"""


# Python Imports
# =============================================================================
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

import config
from classes import SoCTracker, InputFile, IC_handler
from filters import FilterHandler


# Functions
# =============================================================================
def plot(max_pos, pos, t, d_real, d_other, title, xlabel, ylabel, l_sim, l_real, lb, ub):
    plt.subplot(max_pos, 1, pos)
    plt.plot(t, d_real, 'b-', *d_other)
    plt.axis([t[0], t[len(t) - 1], lb, ub])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    patch_dict = {
        '-e' : [mpatches.Patch(color='red'   , label='ekf')],
        '-u' : [mpatches.Patch(color='green' , label='ukf')],
        '-s' : [mpatches.Patch(color='yellow', label='svsf')],
        '-h' : [mpatches.Patch(color='#aa5555', label='ekf with soh'), mpatches.Patch(color='#ff9999', label='soc/soh')],
        '-uh': [mpatches.Patch(color='#005500', label='ukf with soh')],
        '-ck': [mpatches.Patch(color='#555500', label='ekf w/ Coulombs')],
        '-ckh': [mpatches.Patch(color='#005555', label='ekf w/ soh & Coulombs')],
        '-cki': [mpatches.Patch(color='#550055', label='ekf w/ Coulombs & initial SoC')]
    }
    blue_patch = mpatches.Patch(color='blue', label='real')
    patches = [blue_patch]
    for key in patch_dict:
        if key in sys.argv:
            patches += patch_dict.get(key)

    plt.legend(handles=patches, loc=4, prop={'size': 7})


def floatOrNone(string):
    try:
        value = float(string)
    except ValueError:
        value = None
    return value


# Main Program
# =============================================================================

# Initialization
# fileName = 'note1_input_43_2'

fileName = 'note1_input_94_200ms'
# fileName = 'note1_input_43'
# fileName = 'note1_input_54'
# fileName = 'note1_input_87'
# fileName = 'note1_input_57'
# fileName = 'note1_input_69_unknown_init_soc'
# fileName = 'note1_input_74'
# fileName = 'note1_input_94_unknown_init_soc'

inputF = InputFile('Good_Inputs/200ms/' + fileName + '.txt')

kalman_config, svsf_config = inputF.adjustConfig(*config.get()) # Read config file
voltage_actual = []
ini_soc = 0
iniV = 0
# If no options specified, default to ekf
if {'-e', '-u', '-s', '-h', '-uh', '-ck', '-ckh', '-cki'}.isdisjoint(sys.argv):
    sys.argv.append('-e')
    
sum_voltage_error = 0.0
sum_soc_error = 0.0
n = 0.0



# Init dict with default values
tracker_kwargs = {
    'cell_count'      : 1,
    'cell_capacity'   : kalman_config[5]['cnom'],
    'init_soc'        : 0.1,
    'balance_current' : 0
}


# Overwrite defaults with input file values, if they exist
for key in tracker_kwargs:
    try:
        tracker_kwargs[key] = float(inputF[key])
        print('key {}: {}'.format(key, tracker_kwargs[key]))
        if key == 'init_soc':
            print(ini_soc)
            ini_soc = tracker_kwargs[key]
    except KeyError:
        print('Using default value for key {}: {}'.format(key, tracker_kwargs[key]))

soc_tracker = SoCTracker(**tracker_kwargs)
time = []

counter = 0
tempx = []
tempy = []
if '-ic' in sys.argv:
    ekf_voltage = []
    ekf_soc = []

''' version is which set of parameters we want to look at the outputs for '''
version = 1
for line in inputF:
    iniV, current, dt = map(float, line[:3])
    if iniV != 0:
        print("iniV: {}".format(iniV))
        break

handler = FilterHandler(sys.argv, kalman_config, svsf_config, ini_soc, iniV)
for line in inputF:
    counter += 1
    tempx.append(counter)

    # if counter == 3000:
    #     break
    # Read current/voltage, then perform kalman update
    # line[1] is current and line[0] is voltage
    voltage, current, dt = map(float, line[:3])
    tempy.append(-current)

    balance_info = (floatOrNone(line[3]), floatOrNone(line[4])) if len(line) > 4 else (None, None)

    handler.update(dt, voltage, current, version)

    soc_tracker.update(voltage * 1000,  # Voltage
                       current * 1000,  # Current
                       dt,              # delta_t
                       balance_info)    # balance_info (cell, dir)

    time.append(time[-1] + dt if time else dt)
    voltage_actual.append(voltage)

    # for rmse calculation
    n += 1
    if '-cki' in sys.argv:
        sets_to_del = handler.filters['ekfi_coul']['filter'].list_to_pop
        if sets_to_del:
            sum_voltage_error = np.delete(sum_voltage_error,
                                          sets_to_del)  # if we elliminated some possible arrays, we need to drop them from these lists as well
            sum_soc_error = np.delete(sum_soc_error, sets_to_del)
        sum_voltage_error += np.absolute(handler.get_voltage_error(voltage))
        sum_soc_error += np.absolute(handler.get_soc_error(soc_tracker.avg_soc[-1]))
    elif '-ck' in sys.argv:
        sets_to_del = handler.filters['ekf_coul']['filter'].list_to_pop
        if sets_to_del:
            sum_voltage_error = np.delete(sum_voltage_error,
                                          sets_to_del)  # if we elliminated some possible arrays, we need to drop them from these lists as well
            sum_soc_error = np.delete(sum_soc_error, sets_to_del)
        sum_voltage_error += np.absolute(handler.get_voltage_error(voltage))
        sum_soc_error += np.absolute(handler.get_soc_error(soc_tracker.avg_soc[-1]))

    if '-ic' in sys.argv:
        ekf_voltage.append(handler.get_ekf_voltage()[0])
        ekf_soc.append(handler.get_ekf_soc()[0])

# fig = plt.figure(1)
# fig.canvas.set_window_title(fileName)
# # plt.subplot(2, 1, 1)
# x_to_plot = tempx[len(tempx) - 100: len(tempx)]
# y_to_plot = np.poly1d(np.polyfit(tempx[len(tempx) - 100: len(tempx)], tempy[len(tempy) - 100: len(tempy)], 1))(np.unique(tempx[0:100]))
# plt.plot(np.unique(x_to_plot), y_to_plot, tempx, tempy)
# print((y_to_plot[1] - y_to_plot[0]) / (x_to_plot[1] - x_to_plot[0]))
# # plt.subplot(2, 1, 2)
# # plt.plot(tempx, tempy)
# # plt.xlim(0, 15000)
# # plt.ylim(-2.5, 2.5)
# plt.show()
# exit()

# Determine how many subplots we need
max_pos = 2

# Show rms errors
if '-cki' in sys.argv:
    max_pos = 4
    num_of_params = handler.filters['ekfi_coul']['filter'].dim
    abs_case = handler.filters['ekfi_coul']['filter'].absolute_index
    for i in range(num_of_params):
        voltage_rmsd = np.sqrt([sum_voltage_error[i] / n])
        soc_rmsd = np.sqrt([sum_soc_error[i] / n])
        print("for case ", abs_case[i], "\tVOLTAGE\t\tSOC\n\t", voltage_rmsd, "\t", soc_rmsd)
elif '-ck' in sys.argv:
    max_pos = 4
    num_of_params = handler.filters['ekf_coul']['filter'].dim
    abs_case = handler.filters['ekf_coul']['filter'].absolute_index
    for i in range(num_of_params):
        voltage_rmsd = np.sqrt([sum_voltage_error[i] / n])
        soc_rmsd = np.sqrt([sum_soc_error[i] / n])
        print("for case ", abs_case[i], "\tVOLTAGE\t\tSOC\n\t", voltage_rmsd, "\t", soc_rmsd)
else:
    voltage_rmsd = np.sqrt([sum_voltage_error / n])
    soc_rmsd = np.sqrt([sum_soc_error / n])
    print("{\t\tRMSD\n\tVOLTAGE\t\tSOC\n\t", voltage_rmsd, "\t", soc_rmsd)

if '-q' in sys.argv:
    max_pos += 1
if '-h' in sys.argv or '-uh' in sys.argv or '-ckh' in sys.argv:  # or '-cki' in sys.argv:
    max_pos += 1
if '-ic' in sys.argv:
    max_pos += 1

# if '-cki' in sys.argv:
#     max_pos += 2

# plt_soc, plt_volts = handler.get_plotables(time)

# # Plot SoC comparison
# if '-ck' in sys.argv or '-ckh' in sys.argv or '-cki' in sys.argv:
#     counter = 0
#     for soc in soc_tracker.avg_soc:
#         soc_tracker.avg_soc[counter] = soc * tracker_kwargs['cell_capacity']
#         counter += 1
#     plot(max_pos, 1, time, soc_tracker.avg_soc, plt_soc, 'SoC Comparison', 'Time (s)', 'Capacity (As)',
#          'Simulated SoC', 'Actual SoC', 0, 1.1 * tracker_kwargs['cell_capacity'])
# else:
#     plot(max_pos, 1, time, soc_tracker.avg_soc, plt_soc, 'SoC Comparison',
#          'Time (s)', 'SoC',
#          'Simulated SoC', 'Actual SoC', 0, 1.1)
#
# # Plot V comparison
# plot(max_pos, 2, time, voltage_actual, plt_volts, 'Voltage Comparison', 'Time (s)',
#      'Voltage (V)', 'Simulated Voltage', 'Actual Voltage', 3, 5)

if '-h' in sys.argv or '-uh' in sys.argv or '-ckh' in sys.argv:
    plt.subplot(max_pos, 1, 3)
    plt.plot(*handler.get_soh_plots(time))
    plt.title('SoH Estimation')
    plt.xlabel('Time (s)')
    plt.ylabel('SoH')
    plt.xlim(0, len(voltage_actual))
    soh_patch = mpatches.Patch(color='green', label='SoH Estimation')
    plt.legend(handles=[soh_patch], loc=4, prop={'size': 7})

if '-q' in sys.argv:
    # Graph 'estimation quality' and quality thresholds
    jaj = handler.get_quality_plots(time)
    plt.subplot(max_pos, 1, max_pos)
    plt.plot((0, time[-1]), (15, 15), 'y-',     # okay threshold
             (0, time[-1]), (5, 5), 'g-',       # good threshold
             *handler.get_quality_plots(time))
    plt.title('Estimation Quality')
    plt.xlabel('Time (s)')
    plt.ylabel('Quality')

if '-cki' in sys.argv:

    for i in range(num_of_params):

        plt_soc, plt_volts, inputsoc = handler.get_plotables_2(time, abs_case[i] - abs_case[0],
                                                               abs_case[0])  # abs_case offset

        plt.figure(i + 2).canvas.set_window_title('Case {}'.format(abs_case[i]))
        # Plot SoC comparison
        counter = 0
        for soc in soc_tracker.avg_soc:
            soc_tracker.avg_soc[counter] = soc * tracker_kwargs['cell_capacity']
            counter += 1

        position = 1
        plot(max_pos, position, time, soc_tracker.avg_soc, plt_soc, fileName + ' SoC Comparison', 'Time (s)',
             'Capacity (As)',
             'Simulated SoC', 'Actual SoC', 0, 1.1 * tracker_kwargs['cell_capacity'])

        position += 1

        # Plot V comparison
        plot(max_pos, position, time, voltage_actual, plt_volts, 'Voltage Comparison', 'Time (s)',
             'Voltage (V)', 'Simulated Voltage', 'Actual Voltage', 3, 5)

        position += 1

        if '-ic' in sys.argv:

            ICcurve = IC_handler(ekf_voltage, ekf_soc)

            ICvoltage = ICcurve.get_IC_volt()

            ICdQdV = ICcurve.get_IC_QV()

            print("IC peak: {}".format(ICcurve.get_IC_peak()))

            for j in ICvoltage:
                if j >= 3.6:
                    l = int(ICvoltage.index(j))
                    break

            for j in ICvoltage:
                if j >= 4.2:
                    h = int(ICvoltage.index(j))
                    break

            # plt.subplot(max_pos, 1, position)
            plt.figure(99)
            plt.plot(ICvoltage[l:h], ICdQdV[l:h], 'r')
            plt.title('IC Curve')
            plt.xlabel('Voltage (V)')
            plt.ylabel('dQ / dV')
            plt.xlim(min(ICvoltage), max(ICvoltage))
            soc_init_patch = mpatches.Patch(color='red', label='dQ / dV')
            plt.legend(handles=[soc_init_patch], loc=4, prop={'size': 7})
            position += 1
            plt.figure(1)

        # plot soh and soc ini for the 1st case
        plt.subplot(max_pos, 1, position)
        plotinisoc, inisoc = handler.get_soc_init_plots(time, i, abs_case[i] - abs_case[0],
                                                        abs_case[0])  # abs_case offset
        plt.plot(*plotinisoc)
        plt.title('Initial SoC Estimation')
        plt.xlabel('Time (s)')
        plt.ylabel('SoC (A.s)')
        plt.xlim(0, len(time) * dt)
        soc_init_patch = mpatches.Patch(color='purple', label='Initial SoC Estimation')
        plt.legend(handles=[soc_init_patch], loc=4, prop={'size': 7})
        position += 1

        plt.subplot(max_pos, 1, position)
        plotsoh, soh = handler.get_soh_plots_2(time, i, abs_case[i] - abs_case[0], abs_case[0])  # abs_case offset
        plt.plot(*plotsoh)
        plt.title('SoH Estimation')
        plt.xlabel('Time (s)')
        plt.ylabel('SoH (A.s)')
        plt.xlim(0, len(time) * dt)
        soh_patch = mpatches.Patch(color='purple', label='SoH Estimation')
        plt.legend(handles=[soh_patch], loc=4, prop={'size': 7})

        plt.tight_layout()
        plt.subplots_adjust(hspace=1)

        if soh > inisoc + inputsoc or soh < inputsoc:
            print("SoH% (inisoc + inputsoc): {}".format((inisoc + inputsoc) / kalman_config[5]['cnom']))
        else:
            print("SoH% (soh): {}".format(soh / kalman_config[5]['cnom']))

        print(
            "compare ini%: {}".format(inisoc / (inisoc + inputsoc) * ((inisoc + inputsoc) / kalman_config[5]['cnom'])))

if '-ck' in sys.argv:

    for i in range(num_of_params):
        plt_soc, plt_volts, inputsoc = handler.get_plotables_2(time, abs_case[i] - abs_case[0],
                                                               abs_case[0])  # abs_case offset

        plt.figure(i + 2).canvas.set_window_title('Case {}'.format(abs_case[i]))
        # Plot SoC comparison
        counter = 0
        for soc in soc_tracker.avg_soc:
            soc_tracker.avg_soc[counter] = soc * tracker_kwargs['cell_capacity']
            counter += 1
        plot(max_pos, 1, time, soc_tracker.avg_soc, plt_soc, fileName + ' SoC Comparison', 'Time (s)',
             'Capacity (As)',
             'Simulated SoC', 'Actual SoC', 0, 1.1 * tracker_kwargs['cell_capacity'])

        # Plot V comparison
        plot(max_pos, 2, time, voltage_actual, plt_volts, 'Voltage Comparison', 'Time (s)',
             'Voltage (V)', 'Simulated Voltage', 'Actual Voltage', 3, 5)

        # plot soh and soc ini for the 1st case
        plt.subplot(max_pos, 1, 3)
        plotinisoc, inisoc = handler.get_soc_init_plots(time, i, abs_case[i] - abs_case[0],
                                                        abs_case[0])  # abs_case offset
        plt.plot(*plotinisoc)
        plt.title('Initial SoC Estimation')
        plt.xlabel('Time (s)')
        plt.ylabel('SoC (A.s)')
        plt.xlim(0, len(time) * dt)
        soc_init_patch = mpatches.Patch(color='green', label='Initial SoC Estimation')
        plt.legend(handles=[soc_init_patch], loc=4, prop={'size': 7})
        plt.tight_layout()
        print("SoH% (inisoc + inputsoc): {}".format((inisoc + inputsoc) / kalman_config[5]['cnom']))

        print("compare ini%: {}".format(
            inisoc / (inisoc + inputsoc) * ((inisoc + inputsoc) / kalman_config[5]['cnom'])))

# # Show results, or save figure
# fig = plt.figure(1)
# fig.canvas.set_window_title(fileName)
# plt.tight_layout()
# plt.savefig('plots/' + fileName + '_kalman_plot.png', dpi=800, bbox_inches='tight')
plt.show()

print('Done')
