import datetime
import os
import sys
from ast import literal_eval
from copy import deepcopy
from math import acos, atan, pi
from shutil import copyfile

import matplotlib
from matplotlib import pyplot as plt

import config
from classes import InputFile, FrobeniusWrapper, IC_handler
from filters import FilterHandler


class Logger(object):
    def __init__(self, dir):
        self.terminal = sys.stdout
        self.log = open(dir + "/logFile.log", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


now = datetime.datetime.now()
curr = now.strftime("%Y%m%d/%H_%M")
curr = 'plots2/' + curr

if not os.path.exists(curr):
    os.makedirs(curr)
copyfile('classes/matrixHolder.py', curr + '/matrixHolder.py')

logger = Logger(curr)

sys.stdout = logger

sys.argv.append('-cki')

def norm(x, y=0, z=0):
    return (x**2 + y**2 + z**2)**0.5

def spherical_coords(x=0, y=0, z=0):
    '''
        Returns spherical coordinates from cartesian coordinates.

        Coordinates are of the form (r, theta, phi) where:
        r is the radius of the sphere,
        theta is the azimuthal angle,
        phi is the inclination angle
    '''

    r = norm(x, y, z)
    phi = acos(z / r)

    if x:
        ref_theta = atan(y / x)

        theta = ref_theta
        if x < 0:
            theta += pi

        if theta < 0:
            theta += 2 * pi

    else:
        theta = 3 * pi / 2 if y < 0 else pi / 2

    return (r, theta, phi)


# clusters = ['input_clusters/high_health_note1.txt',
#             ]
clusters = ['input_clusters/note1_3.9_to_4.2.txt',
            ]

files_to_compare, known_sohs, colour_ranges= [], [], []
for cluster in clusters:
    with open(cluster, 'r') as f:
        for line in f:
            key, value = line.replace(' ', '').split('=')
            if key == 'input_files':
                files_to_compare.append(literal_eval(value))
            if key == 'states_of_health':
                known_sohs += literal_eval(value)
            if key == 'colour_range':
                colour_ranges.append(literal_eval(value))



input_files = [InputFile(filename) for filesets in files_to_compare for filename in filesets]

adjusted_configs = [input_files[i].adjustConfig(*deepcopy(config.get()))
                    for i in range(len(input_files))]
# Init dict with default values
tracker_kwargs = {
    'cell_count': 1,
    'cell_capacity': 9000,
    'init_soc': 0.1,
    'balance_current': 0
}

ini_soc = []
for i in range(len(input_files)):

    for key in tracker_kwargs:
        try:
            tracker_kwargs[key] = input_files[i][key]
            if key == 'init_soc':
                ini_soc.append(float(tracker_kwargs[key]))
                print("init_soc: {}".format(ini_soc))
        except KeyError:
            print('Using default value for key {}: {}'.format(key, tracker_kwargs[key]))

iniV = []
for i in range(len(input_files)):
    for line in (input_files[i]):
        tempV, current, dt = map(float, line[:3])
        if tempV != 0:
            iniV.append(tempV)
            print("iniV: {}".format(iniV))
            break

handlers = [FilterHandler(['-cki'], *adjusted_configs[i], ini_soc[i], iniV[i]) for i in range(len(input_files))]

colours = []
# Build a list of colours
for i in range(len(files_to_compare)):
    fileset = files_to_compare[i]
    c_start, c_end = colour_ranges[i]
    diff_divisor = len(fileset) - 1 if len(fileset) > 1 else 1

    r_incr = (((c_end & 0xFF0000) >> 16) - ((c_start & 0xFF0000) >> 16)) // diff_divisor
    g_incr = (((c_end & 0x00FF00) >> 8)  - ((c_start & 0x00FF00) >> 8) ) // diff_divisor
    b_incr = (( c_end & 0x0000FF)        - ( c_start & 0x0000FF)       ) // diff_divisor


    for j in range(len(fileset)):
        # Create hex value from RGB values.
        hex_value = hex(c_start + ((r_incr * j) << 16) + ((g_incr * j) << 8) + (b_incr * j))
        # For any colour with a red value < 16, we need to pad to 6 hex characters
        hex_string = '#{:0>6}'.format(str(hex_value).replace('0x', ''))
        colours.append(hex_string)

totalsoh_fig = plt.figure()
totalsoh_fig.suptitle('Input+Initial SOH')
totalsoh_plot = totalsoh_fig.add_subplot(111)

if '-h' in sys.argv:
    soh_fig = plt.figure()
    soh_fig.suptitle('Estimated SOH')
    soh_plot = soh_fig.add_subplot(111)

if '-fb' in sys.argv:
    frob_fig = plt.figure()
    frob_fig.suptitle('Frobenius Norm')
    frob_plot = frob_fig.add_subplot(111)

    frob_slope_fig = plt.figure()
    frob_slope_fig.suptitle('Frobenius Min slope')
    frob_slope_plot = frob_slope_fig.add_subplot(111)

    att_mag_fig = plt.figure()
    att_mag_fig.suptitle('Attribute magnitudes')
    att_mag_plot = att_mag_fig.add_subplot(111)

    att_angle_fig = plt.figure()
    att_angle_fig.suptitle('Attribute angles')
    att_angle_plot = att_angle_fig.add_subplot(111)

if '-ic' in sys.argv:
    ic_fig = plt.figure()
    ic_fig.suptitle('IC curve')
    ic_plot = ic_fig.add_subplot(111)

markers = ['.', 'o', '*', '+', 'x', 'X', 'D', 'p', 'h', 's']
markers_used = []
_file_handles = []

index = 0
for input_file, handler in zip(input_files, handlers):
    print('\nExpected SoH: {}% -> {} A.s'.format(known_sohs[index],
                                                 known_sohs[index] * handler.filters['ekfi_coul']['filter'].p[0]['cnom'] / 100))

    cases = handler.filters['ekfi_coul']['filter'].absolute_index
    time = []
    num_of_params = 0
    sohs = []
    icv = []
    icq = []
    ekf_voltage = []
    ekf_soc = []
    totalsohs = []

    xs, ys, zs = [[0] for _ in range(3)]

    if '-cki' in sys.argv:
        abs_case = handler.filters['ekfi_coul']['filter'].absolute_index
        num_of_params = handler.filters['ekfi_coul']['filter'].dim

    for line in input_file:
        voltage, current, dt = map(float, line[:3])
        handler.update(dt, voltage, current, 1)
        time.append(time[-1] + dt if time else dt)

    for i in range(num_of_params):
        totalsoh = []

        soh, lastsoh = handler.get_soh_plots_2(time, i, abs_case[i] - abs_case[0], abs_case[0])

        ini_soc, lastinisoc = handler.get_soc_init_plots(time, i, abs_case[i] - abs_case[0], abs_case[0])

        inputsoc, volt, lastinputsoc = handler.get_plotables_2(time, abs_case[i] - abs_case[0], abs_case[0])

        for j in range(len(ini_soc[1])):
            totalsoh.append(ini_soc[1][j] + inputsoc[1][j])
        totalsohs.append(totalsoh)

        totalsoh_plot.plot(time, totalsoh, colours[index])
        totalsoh_plot.annotate('{}'.format(abs_case[i]), xy=(ini_soc[0][-1], totalsohs[i][-1]),
                               xytext=(ini_soc[0][-1] + 200, totalsohs[i][-1]),
                               arrowprops=dict(facecolor='black', headwidth=2))

        if '-h' in sys.argv:
            sohs.append(soh)

            soh_plot.plot(time, soh[1], colours[index])
            soh_plot.annotate('{}'.format(abs_case[i]), xy=(sohs[i][0][-1], sohs[i][1][-1]),
                              xytext=(sohs[i][0][-1] + 200, sohs[i][1][-1]),
                              arrowprops=dict(facecolor='black', headwidth=2))

        if '-fb' in sys.argv:
            normal_soh = FrobeniusWrapper(time, soh[1]).normalized
            zs.append(normal_soh.min_slope[0])

        if '-ic' in sys.argv:
            ICvoltage = []
            ICdQdV = []

            ICcurve = IC_handler(volt[1], inputsoc[1])

            ICvoltage = ICcurve.get_IC_volt()

            ICdQdV = ICcurve.get_IC_QV()

            ICpeak = ICcurve.get_IC_peak()

            ic_plot.plot(ICvoltage, ICdQdV, colours[index])

            if not cases[i] in markers_used:
                markers_used.append(cases[i])

    if '-fb' in sys.argv:

        frobs = handler.get_frob_plots_2(time, abs_case)

        for i in range(0, len(frobs[1][0])):
            frob = FrobeniusWrapper(time, frobs[1][i])

            normal_frob = frob.normalized

            gaus_frob = frob.normalized.gaussianized

            frob_plot.plot(gaus_frob.time, gaus_frob.data, colours[index])

            ys.append(gaus_frob.min_slope[0])

            mini_slope, steepest_point = frob.min_slope

            if not cases[i] in markers_used:
                markers_used.append(cases[i])

            frob_slope_plot.plot(known_sohs[index], mini_slope, colours[index], marker=markers[cases[i]])
            frob_plot.plot(frob.time[steepest_point], frob[steepest_point], colours[index], marker=markers[cases[i]])

            if len(normal_frob.trough_value) > 0:
                xs.append(normal_frob.trough_value[0])
            else:
                xs.append(0)

        # att_plot.plot(xs=xs, ys=ys, zs=zs, c=colours[index])
        for i in range(1, len(frobs[1][0]) + 1):

            _min = 1000
            _max = -1
            for j in known_sohs:
                _min = j if j < _min else _min
                _max = j if j > _max else _max

            x_len = _max - _min

            r, theta, phi = spherical_coords(xs[i], ys[i], zs[i])

            att_mag_plot.plot(known_sohs[index], r, colours[index], marker='.')
            att_angle_plot.plot(known_sohs[index], phi, colours[index], marker='^')
            att_angle_plot.plot(known_sohs[index], theta, colours[index], marker='v')
            att_angle_plot.plot(known_sohs[index], phi + theta, colours[index], marker='D')

            att_mag_plot.annotate('{}'.format(abs_case[i - 1]), xy=(known_sohs[index] + x_len * 0.005, r),
                                  xytext=(known_sohs[index] + x_len * 0.05, r),
                                  arrowprops=dict(facecolor='black', headwidth=2, width=1))

            att_angle_plot.annotate('{}'.format(abs_case[i - 1]), xy=(known_sohs[index] + x_len * 0.005, phi),
                                    xytext=(known_sohs[index] + x_len * 0.05, phi),
                                    arrowprops=dict(facecolor='black', headwidth=2, width=1))
            att_angle_plot.annotate('{}'.format(abs_case[i - 1]), xy=(known_sohs[index] + x_len * 0.005, theta),
                                    xytext=(known_sohs[index] + x_len * 0.05, theta),
                                    arrowprops=dict(facecolor='black', headwidth=2, width=1))
            att_angle_plot.annotate('{}'.format(abs_case[i - 1]), xy=(known_sohs[index] + x_len * 0.005, phi + theta),
                                    xytext=(known_sohs[index] + x_len * 0.05, phi + theta),
                                    arrowprops=dict(facecolor='black', headwidth=2, width=1))

    _file_handles.append(matplotlib.lines.Line2D([], [], color=colours[index], label='File {}%'.format(known_sohs[index])))

    index += 1

_file_labels = [h.get_label() for h in _file_handles]

_handles = []

for i in markers_used:
    mark = markers[i]
    _handles.append(matplotlib.lines.Line2D([], [], marker=mark, color='black', label='Case {}'.format(i)))

_labels = [h.get_label() for h in _handles]
fnumber = 1


plt.xlabel('time', figure=plt.figure(fnumber))
plt.ylabel('input + initial SOC', figure=plt.figure(fnumber))
lgd = totalsoh_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

fnumber += 1

if '-h' in sys.argv:
    plt.xlabel('time', figure=plt.figure(fnumber))
    plt.ylabel('Estimated SOH', figure=plt.figure(fnumber))

    lgd2 = soh_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    # soh_fig.savefig(curr + '/SoH', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
    fnumber += 1


if '-ic' in sys.argv:
    plt.xlabel('voltage', figure=plt.figure(fnumber))
    plt.ylabel('dQ / dV', figure=plt.figure(fnumber))

    lgd3 = ic_plot.legend(handles=_file_handles, labels=_file_labels, loc=4)

    ic_fig.gca().add_artist(lgd3)

    fnumber += 1

if '-fb' in sys.argv:
    _special_angle_handles = []
    _special_angle_handles.append(matplotlib.lines.Line2D([], [], marker='^', color='black', label='Phi'))
    _special_angle_handles.append(matplotlib.lines.Line2D([], [], marker='v', color='black', label='Theta'))
    _special_angle_handles.append(matplotlib.lines.Line2D([], [], marker='D', color='black', label='Phi + Theta'))
    _special_angle_labels = [h.get_label() for h in _special_angle_handles]

    plt.xlabel('Time', figure=plt.figure(fnumber))
    plt.ylabel('Minimum Norm', figure=plt.figure(fnumber))
    fnumber += 1
    plt.xlabel('Known SoH', figure=plt.figure(fnumber))
    plt.ylabel('Minimum Slope', figure=plt.figure(fnumber))
    fnumber += 1
    plt.xlabel('Known SoH', figure=plt.figure(fnumber))
    plt.ylabel('Magnitude', figure=plt.figure(fnumber))
    fnumber += 1
    plt.xlabel('Known SoH', figure=plt.figure(fnumber))
    plt.ylabel('Angle', figure=plt.figure(fnumber))
    fnumber += 1
    plt.grid(figure=plt.figure(3))

    lgd1 = frob_slope_plot.legend(handles=_handles, labels=_labels, bbox_to_anchor=(1.05, 0.25), loc=2,
                                  borderaxespad=0.)
    lgd2 = frob_slope_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1.05, 1), loc=2,
                                  borderaxespad=0.)
    frob_slope_fig.gca().add_artist(lgd1)
    frob_slope_fig.savefig(curr + '/Slopes', bbox_extra_artists=(lgd1, lgd2), bbox_inches='tight', dpi=600)

    lgd1 = frob_plot.legend(handles=_handles, labels=_labels, bbox_to_anchor=(1.05, 0.25), loc=2, borderaxespad=0.)
    lgd2 = frob_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1.05, 1), loc=2,
                            borderaxespad=0.)
    frob_fig.gca().add_artist(lgd1)
    frob_fig.savefig(curr + '/Frobenius', bbox_extra_artists=(lgd1, lgd2), bbox_inches='tight', dpi=600)

    lgd = att_mag_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1.05, 1), loc=2,
                              borderaxespad=0.)
    att_mag_fig.savefig(curr + '/Magnitude', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600)

    lgd1 = att_angle_plot.legend(handles=_special_angle_handles, labels=_special_angle_labels,
                                 bbox_to_anchor=(1.05, 0.25),
                                 loc=2, borderaxespad=0.)
    lgd2 = att_angle_plot.legend(handles=_file_handles, labels=_file_labels, bbox_to_anchor=(1.05, 1), loc=2,
                                 borderaxespad=0.)
    att_angle_fig.gca().add_artist(lgd1)
    att_angle_fig.savefig(curr + '/Angle', bbox_extra_artists=(lgd1, lgd2), bbox_inches='tight', dpi=600)



logger.log.close()

plt.show()
