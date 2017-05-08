from filters.kalman import Kalman, KalmanSoh
from filters.coulomb_kalman import Kalman as KalmanC
from filters.coulomb_kalman import KalmanSoh as KalmanCSoh
from filters.coulomb_kalman import KalmanSoC as KalmanCSoC
from filters.svsf import SVSF
from filters.unscented import Unscented, UnscentedSoh
from numpy import linalg
from operator import sub
norm = linalg.norm


class FilterHandler:
    '''
    This handler allows us to plot any combination of extended kalman filter, unscented
    kalman filter and smooth variable structure filter. It uses flags from the command line to
    choose which filters to use, and keeps track of all necessary information.
    '''

    def __init__(self, argv: list, kalman_config: list, svsf_config: list, inisoc: float, iniV: float, classes):

        filters = {}

        if '-e' in argv:
            filters['ekf'] = {
                'filter'      : Kalman(*kalman_config),
                'soc'         : [],
                'volts'       : [],
                'quality'     : [],
                'frob'        : [],
                'main_colour' : 'r-'
            }

        if '-h' in argv:
            filters['ekfh'] = {
                'filter'      : KalmanSoh(*kalman_config),
                'soc'         : [],
                'volts'       : [],
                'soh'         : [],
                'frob'        : [],
                'main_colour' : '#aa5555',
                'alt_colour'  : '#ff9999'
            }

        if '-u' in argv:
            filters['ukf'] = {
                'filter'      : Unscented(*kalman_config),
                'soc'         : [],
                'volts'       : [],
                'frob'        : [],
                'main_colour' : 'g-'
            }

        if '-uh' in argv:
            filters['ukfh'] = {
                'filter'      : UnscentedSoh(*kalman_config),
                'soc'         : [],
                'volts'       : [],
                'soh'         : [],
                'frob'        : [],
                'main_colour' : '#005500',
                'alt_colour'  : '#99ff99'
            }

        if '-s' in argv:
            filters['svsf'] = {
                'filter'      : SVSF(*svsf_config),
                'soc'         : [],
                'volts'       : [],
                'main_colour' : 'y-'
            }

        if '-ck' in argv:
            filters['ekf_coul'] = {
                'filter': KalmanC(*kalman_config, inisoc, iniV),
                'soc'         : [],
                'volts'       : [],
                'soc_init': [],
                'frob': [],
                'main_colour' : 'm-'
            }

        if '-ckh' in argv:
            filters['ekfh_coul'] = {
                'filter': KalmanCSoh(*kalman_config),
                'soc': [],
                'soh': [],
                'volts': [],
                'main_colour': '#005500',
                'alt_colour': '#99ff99'
            }

        if '-cki' in argv:
            filters['ekfi_coul'] = {
                'filter': KalmanCSoC(*kalman_config, inisoc, iniV, classes),
                'soc': [],
                'soh': [],
                'soc_init': [],
                'volts': [],
                'frob': [],
                'main_colour': '#550055',
                'alt_colour': '#99ff99'
            }

        self.filters = filters
        self.for_coul_total_popped = []

    def update(self, dt: float, voltage: float, current: float, version):
        ''' Updates all instantiated filters and collects pertinent information '''

        filters = self.filters
        if 'ekf' in filters:
            dict = filters['ekf']
            filter = dict['filter']
            filter.read_inputs(current, voltage)
            filter.update(dt)
            dict['soc'].append(filter.get_simulated_soc())
            dict['volts'].append(filter.get_simulated_voltage())

            dict['frob'].append(norm(filter.P))

            dict['quality'].append(filter.get_mean_percent_diff())

        if 'ekfh' in filters:
            dict = filters['ekfh']
            filter = dict['filter']
            filter.read_inputs(current, voltage)
            filter.update(dt)
            dict['soc'].append(filter.get_simulated_soc())
            dict['volts'].append(filter.get_simulated_voltage())

            dict['frob'].append(norm(filter.P))

            dict['soh'].append(filter.get_soh())

        if 'ukf' in filters:
            dict = filters['ukf']
            filter = dict['filter']
            filter.update(dt, voltage, current)
            dict['soc'].append(filter.get_soc())
            dict['volts'].append(filter.get_estimated_voltage(current))

            dict['frob'].append(norm(filter.P))


        if 'ukfh' in filters:
            dict = filters['ukfh']
            filter = dict['filter']
            filter.update(dt, voltage, current)
            dict['soc'].append(filter.get_soc())
            dict['volts'].append(filter.get_estimated_voltage(current))

            dict['frob'].append(norm(filter.P))

            dict['soh'].append(filter.get_soh())

        if 'svsf' in filters:
            dict = filters['svsf']
            filter = dict['filter']
            filter.update(dt, voltage, current)
            dict['soc'].append(filter.get_soc())
            dict['volts'].append(filter.get_estimated_voltage(current))


        if 'ekf_coul' in filters:
            dict = filters['ekf_coul']
            filter = dict['filter']
            filter.update(dt, voltage, current)

            for i in self.filters['ekf_coul']['filter'].list_to_pop:
                self.for_coul_total_popped.append(i)

            temp = filter.get_inputsoc_2(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['soc'].append(temp)

            frobv = []
            temp = filter.get_P()
            for i in range(len(temp)):
                frobv.append(norm(temp[i]))
            dict['frob'].append(frobv)

            temp = filter.get_estimated_voltage_new_2(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['volts'].append(temp)

            temp = filter.get_inisoc(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['soc_init'].append(temp)

        if 'ekfh_coul' in filters:
            dict = filters['ekfh_coul']
            filter = dict['filter']
            filter.update(dt, voltage, current)
            dict['soc'].append(filter.get_soc())
            dict['volts'].append(filter.get_estimated_voltage())

            dict['soh'].append(filter.get_soh())

        if 'ekfi_coul' in filters:
            dict = filters['ekfi_coul']
            filter = dict['filter']
            filter.update(dt, voltage, current)

            for i in self.filters['ekfi_coul']['filter'].list_to_pop:
                self.for_coul_total_popped.append(i)

            temp = filter.get_inputsoc_2(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['soc'].append(temp)

            frobv = []
            temp = filter.get_P()
            for i in range(len(temp)):
                frobv.append(norm(temp[i]))
            dict['frob'].append(frobv)

            temp = filter.get_soh(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['soh'].append(temp)

            temp = filter.get_estimated_voltage_new_2(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['volts'].append(temp)

            temp = filter.get_inisoc(0)
            for i in self.for_coul_total_popped:
                temp.insert(i, -1)
            dict['soc_init'].append(temp)


    def get_plotables(self, t: list) -> '2 tuples of n*(t, y, line-type), n = number of filters':
        socs, volts = (), ()

        for key in self.filters:
            filter = self.filters[key]
            main_colour = filter['main_colour']

            socs += (t, filter['soc'], main_colour)
            if 'soh' in filter:
                socf = []
                for soc, soh in zip(filter['soc'], filter['soh']):
                    socf.append(soc / soh)
                socs += (t, socf, filter['alt_colour'])
            elif 'soc_init' in filter:
                soc_init_f = []
                for soc, soc_init in zip(filter['soc'], filter['soc_init']):
                    soc_init_f.append(soc_init)
                socs += (t, soc_init_f, filter['alt_colour'])

            volts += (t, filter['volts'], main_colour)

        return socs, volts

    def get_plotables_2(self, t: list, abs_version,
                        version_offset) -> '2 tuples of n*(t, y, line-type), n = number of filters':
        socs, volts = (), ()

        for key in self.filters:
            filter = self.filters[key]
            main_colour = filter['main_colour']

            socs += (t, [filter['soc'][i][abs_version] for i in range(len(t))], main_colour)
            volts += (t, [filter['volts'][i][abs_version] for i in range(len(t))], main_colour)
        print("\n\nSoC Input Case {}: {}".format(abs_version + version_offset, socs[-2][-1]))
        return socs, volts, socs[1][-1]

    def get_quality_plots(self, t: list) -> '1 tuple of n*(t, y, line-type), n = number of filters':
        # Currently, only ekf has quality implemented
        qualities = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'quality' in filter_dict:
                qualities += (t, filter_dict['quality'], filter_dict['main_colour'])

        return qualities

    def get_soh_plots(self, t: list):
        soh = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'soh' in filter_dict:
                soh += (t, filter_dict['soh'], filter_dict['main_colour'])
        print("SOH Case 1: {}".format(soh[-2][-1]))
        return soh

    def get_soc_init_plots(self, t: list, version, abs_version, version_offset):
        soc_init = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'soc_init' in filter_dict:
                temp = []
                for i in range(len(filter_dict['soc_init'])):
                    temp.append(filter_dict['soc_init'][i][abs_version])
                soc_init += (t, temp, filter_dict['main_colour'])
        print("SoC Init Case {}: {}".format(abs_version + version_offset, soc_init[-2][-1]))
        return soc_init, soc_init[-2][-1]

    def get_soh_plots_2(self, t: list, version, abs_version, version_offset):
        soh = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'soh' in filter_dict:
                temp = []
                for i in range(len(filter_dict['soh'])):
                    temp.append(filter_dict['soh'][i][abs_version])
                soh += (t, temp, filter_dict['main_colour'])
        print("SoH Case {}: {} -> {}".format(abs_version + version_offset, soh[-2][-1], soh[-2][-1] / 9000))
        lastsoh = soh[-2][-1]
        return soh, lastsoh

    #
    # def get_soc_init_plots_2(self, t: list):
    #     soc_init = ()
    #     filters = self.filters
    #
    #     for key in filters:
    #         filter_dict = filters[key]
    #         if 'soc_init_2' in filter_dict:
    #             soc_init += (t, filter_dict['soc_init_2'], filter_dict['main_colour'])
    #      print("SoC Init Case 2: {}".format(soc_init[-2][-1]))
    #     return soc_init

    def get_frob_plots(self, t: list):
        frobs = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'frob' in filter_dict:
                frobs += (t, filter_dict['frob'], filter_dict['main_colour'])

        return frobs

    def get_frob_plots_2(self, t: list, abs_version):
        frobs = ()
        filters = self.filters

        for key in filters:
            filter_dict = filters[key]
            if 'frob' in filter_dict:
                temp = []
                for i in range(len(filter_dict['frob'])):
                    temp.append(filter_dict['frob'][i])
                frobs += (t, temp, filter_dict['main_colour'])

        return frobs

    def get_voltage_error(self, voltage):
        # Used for RMSE calculation.
        if 'ekf' in self.filters:
            return abs(self.filters['ekf']['filter'].get_simulated_voltage() - voltage)
        elif 'ekfh' in self.filters:
            return abs(self.filters['ekfh']['filter'].get_simulated_voltage() - voltage)
        elif 'ekf_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekf_coul']['filter'].get_simulated_voltage_new()
            for i in range(self.filters['ekf_coul']['filter'].dim):
                temp.append(abs(temp2[i] - voltage))
            return temp

        elif 'ekfh_coul' in self.filters:
            return abs(self.filters['ekfh_coul']['filter'].get_simulated_voltage() - voltage)
        elif 'ekfi_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekfi_coul']['filter'].get_simulated_voltage_new()
            for i in range(self.filters['ekfi_coul']['filter'].dim):
                temp.append(abs(temp2[i] - voltage))
            return temp
        else:
            return 0

    def get_soc_error(self, soc):
        # Used for RMSE calculation.
        if 'ekf' in self.filters:
            return abs(self.filters['ekf']['filter'].get_simulated_soc() - soc)
        elif 'ekfh' in self.filters:
            return abs(self.filters['ekfh']['filter'].get_simulated_soc() - soc)
        elif 'ekf_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekf_coul']['filter'].get_simulated_soc_new()
            for i in range(self.filters['ekf_coul']['filter'].dim):
                temp.append(abs(temp2[i] - soc * self.filters['ekf_coul']['filter'].p[0]['cnom']))
            return temp
        elif 'ekfh_coul' in self.filters:
            return abs(
                self.filters['ekfh_coul']['filter'].get_simulated_soc() - soc * self.filters['ekfh_coul']['filter'].p[
                    'cnom'])
        elif 'ekfi_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekfi_coul']['filter'].get_simulated_soc_new()
            for i in range(self.filters['ekfi_coul']['filter'].dim):
                temp.append(abs(temp2[i] - soc * self.filters['ekfi_coul']['filter'].p[0]['cnom']))
            return temp
        else:
            return 0

    def get_ekf_voltage(self):

        if 'ekf' in self.filters:
            return self.filters['ekf']['filter'].get_simulated_voltage()
        elif 'ekfh' in self.filters:
            return self.filters['ekfh']['filter'].get_simulated_voltage()
        elif 'ekf_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekf_coul']['filter'].get_simulated_voltage_new()
            for i in range(self.filters['ekf_coul']['filter'].dim):
                temp.append(temp2[i])
            return temp

        elif 'ekfh_coul' in self.filters:
            return self.filters['ekfh_coul']['filter'].get_simulated_voltage()
        elif 'ekfi_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekfi_coul']['filter'].get_simulated_voltage_new()
            for i in range(self.filters['ekfi_coul']['filter'].dim):
                temp.append(temp2[i])
                return temp
        else:
            return 0

    def get_ekf_soc(self):

        if 'ekf' in self.filters:
            return self.filters['ekf']['filter'].get_simulated_soc()
        elif 'ekfh' in self.filters:
            return self.filters['ekfh']['filter'].get_simulated_soc()
        elif 'ekf_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekf_coul']['filter'].get_simulated_soc_new()
            for i in range(self.filters['ekf_coul']['filter'].dim):
                temp.append(temp2[i])
            return temp
        elif 'ekfh_coul' in self.filters:
            return self.filters['ekfh_coul']['filter'].get_simulated_soc()
        elif 'ekfi_coul' in self.filters:
            temp = []
            temp2 = self.filters['ekfi_coul']['filter'].get_simulated_soc_new()

            for i in range(self.filters['ekfi_coul']['filter'].dim):
                temp.append(temp2[i])
            return temp
        else:
            return 0
