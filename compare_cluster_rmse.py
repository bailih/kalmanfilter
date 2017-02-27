from filters import FilterHandler
from classes import InputFile, SoCTracker
from math import sqrt
from matplotlib import pyplot as plt
from matplotlib import patches
import config
from ast import literal_eval

cluster = 'input_clusters/c8000_note1_87.txt'

files_to_compare, prime_sohs = [], []
with open(cluster, 'r') as clust_f:
    for line in clust_f:
        key, value = line.replace(' ', '').split('=')
        if key == 'input_files':
            files_to_compare += (literal_eval(value))
        if key == 'states_of_health':
            prime_sohs += literal_eval(value)

tracker_kwargs = {}

drum = []
for file, soh in zip(files_to_compare, prime_sohs):
    inF = InputFile(file)
    kalman_config = inF.adjustConfig(*config.get())[0]

    for key in ('cell_count', 'cell_capacity', 'init_soc', 'balance_current'):
        try:
            tracker_kwargs[key] = float(inF[key])
        except KeyError:
            tracker_kwargs[key] = None

    socTracker = SoCTracker(**tracker_kwargs)
    handler = FilterHandler(['-h'], kalman_config, None)

    time = []
    sum_soc_error = 0
    n = 0
    for line in inF:
        voltage, current, dt = map(float, line[:3])

        handler.update(dt, voltage, current)
        socTracker.update(voltage * 1000, current * 1000, dt, (None, None))

        time.append(time[-1] + dt if time else dt)

        soc_error = abs(handler.filters['ekfh']['filter'].get_simulated_voltage() - voltage)
        sum_soc_error += soc_error

        n += 1

    soc_rmse = sqrt(sum_soc_error / n)
    print(soc_rmse)


    fig = plt.figure()
    plot = plt.axes()

    title = 'Soc vs Time - {} SoH'.format(soh)
    plt.title(title)
    plt.xlabel('Time (s)')
    plt.ylabel('SoC')

    soc_plots =  handler.get_plotables(time)[0]
    plt.plot(time, socTracker.avg_soc, 'b-', *soc_plots)

    plt.axis([0, time[-1], 0, 1])

    legend_patches = [
        patches.Patch(color='Blue', label='Coulomb counting'),
        patches.Patch(color=handler.filters['ekfh']['main_colour'], label='Estimated SoC'),
        patches.Patch(color=handler.filters['ekfh']['alt_colour'], label='Estimated SoC / Estimated SoH')
    ]
    plot.legend(handles=legend_patches, loc=0, prop={'size':12})

plt.show()
