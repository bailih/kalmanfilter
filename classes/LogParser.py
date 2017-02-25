import csv

'''
    Class to parse the log files spit out by the tourmaline battery charger.

    The strategy is to build 4 lists we then plug into the state of charge tracker:
        1. Voltage          (volts)
        2. Current          (amps)
        3. Change in time   (delta_t)
        4. The current      (balance_info)
           balancing cell
           and direction

    The two public methods are:

        - parse:    takes in the filename of the log to parse
                             the charging current at full power
                             the approximate sampling time used
                             the filename for the kalman filter input file (OPTIONAL)

        - get_info: returns info from the 4 lists at the specified index
'''

class LogParser:
    def __init__(self, millivolt_tolerance):
        self.MILLIVOLT_TOLERANCE = millivolt_tolerance

        self.volts        = []
        self.currents     = []
        self.delta_t      = []
        self.balance_info = []
        self.len          = 0

        self.lines = []

        self.line_volts_prev = None
        self.line_time_prev  = None

        self.bad_point_count = 0

    def parse(self, filename, full_charging_current, sample_time, file_to_write=None):

        self._read_csv(filename)

        prev_volt    = None
        prev_current = None
        prev_time    = None

        bad_volts_count    = 0
        bad_current_count  = 0
        tolerance_increase = 0

        starting = True
        stopping = False
        for line in self.lines:

            time, volt, current, balance_str = (line[0],                         # eg. 09-17 10:45:55 19.00231
                                                float(line[1]) * 14 * 1000 // 1, # line[1] is average cell voltage
                                                float(line[17]),                 # charging current in mA
                                                line[19].split(' '))             # [ balance_cell, balance_dir }

            # No actual support here for testings on a month switch. Add it later if it's needed.
            month_day, hms, t_clk   = time.split(' ')
            month, day              = month_day.split('-')
            hours, minutes, seconds = hms.split(':')
            time = 86400 * int(day) + 3600 * int(hours) + 60 * int(minutes) + int(seconds)



            volts_diff = volt - prev_volt if prev_volt else 0
            if volts_diff < self.MILLIVOLT_TOLERANCE + tolerance_increase:

                self.volts.append(volt)

                # Assume the voltage increased linearly through bad points.
                for i in range(bad_volts_count):
                    index = len(self.volts) - bad_volts_count - 1 + i
                    self.volts[index] = prev_volt + (i + 1) * ((volt - prev_volt) / (prev_volt + 1))

                prev_volt = volt
                tolerance_increase = 0
                bad_volts_count = 0

            else:
                # For bad points, we append None so we can come back and fix it later
                self.volts.append(None)
                tolerance_increase += 50
                bad_volts_count += 1



            # Once we get above the normal charging current, we know we are done the starting step
            if starting and current and current >= full_charging_current:
                starting = False

            # If we dip below the charging current, start checking the remaining currents. If we don't see
            # any currents above the charging current, that means we are stopping the charge.
            if current < full_charging_current and not starting and not stopping:
                stopping = True
                for i in range(self.lines.index(line) + 1, len(self.lines) - 1):
                    if float(self.lines[i][17]) >= full_charging_current:
                        stopping = False
                        break

            if current >= full_charging_current or starting or stopping:

                self.currents.append(current)

                # The current should be relatively stable, so replace bad points with the average of the currents
                # on either side
                for i in range(bad_current_count):
                    index = len(self.currents) - bad_current_count - 1 + i
                    self.currents[index] = (current + prev_current) / 2

                prev_current = current
                bad_current_count = 0


            else:
                self.currents.append(None)
                bad_current_count += 1



            # We calculate the delta t for all but the first point, where we assume it's the exact sample time
            self.delta_t.append(time - prev_time if prev_time != None else sample_time)
            prev_time = time



            if 'No' in balance_str:
                balance_cell, balance_dir = None, None
            else:
                balance_cell = balance_str[1]
                balance_dir = {
                    'charge': 0,
                    'discharge': 1
                }.get(balance_str[0])
            self.balance_info.append((balance_cell, balance_dir))



        # If we get have some straggling bad points at the end, just remove them.
        if bad_volts_count or bad_current_count:
            if bad_volts_count > bad_current_count:
                pops = range(bad_volts_count)
            else:
                pops = range(bad_current_count)

            for _ in pops:
                self.volts.pop()
                self.currents.pop()



        # Ensure the lengths of all lists are the same
        length_tuple = (len(self.volts), len(self.currents), len(self.delta_t), len(self.balance_info))
        if length_tuple.count(length_tuple[0]) != len(length_tuple):
            raise Exception('Lengths of log lists are not all equal!' +
                            '\tVolts : {}\tCurrents : {}\tDelta_t : {}\tBalance_info : {}'.format(*length_tuple))
        self.len = length_tuple[0]



        #Debug only, chumpi
        for i in range(self.len):
            if None in (self.volts[i], self.currents[i], self.delta_t[i], self.balance_info[i]):
                raise Exception('Missed a None! Index = {}'.format(i))



        # Allow us to build a file for use with the kalman filter.
        if file_to_write:
            with open(file_to_write, 'w') as kalman_file:
                time_sum = 0
                for i in range(len(self.volts)):
                    kalman_file.write('{}\t{}\t{}\t{}\t{}\n'.format(self.volts[i] / 1000,
                                                                    self.currents[i] / 1000,
                                                                    self.delta_t[i],
                                                                    self.balance_info[i][0],
                                                                    self.balance_info[i][1]))


    # Adds every line of the specified csv to parser's list of lines
    def _read_csv(self, filename):
        with open(filename, 'r') as f:
            rc = csv.reader(f, delimiter=',')
            next(rc, None)  # Skip header

            try:
                for line in rc:
                    self.lines.append(line)
            except:
                pass


    def get_info(self, index):
        return (self.volts[index], self.currents[index], self.delta_t[index], self.balance_info[index])