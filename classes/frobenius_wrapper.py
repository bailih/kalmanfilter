from math import pow, e, pi, sqrt
from numpy import convolve

class FrobeniusWrapper:

    def __init__(self, time, data):
        self.time = tuple(time)
        self.data = tuple(data)

    @property
    def trough_indices(self):
        ''' Returns index of the first local minimum. '''

        try:
            return self._trough_indices

        # If _trough_index doesn't exist, we need to create it.
        except AttributeError:
            self._find_troughs_and_peaks()
            return self._trough_indices

    @property
    def peak_indices(self):
        ''' Returns the index of the first local maximum after the first local minimum. '''
        try:
            return self._peak_indices

        except AttributeError:
            self._find_troughs_and_peaks()
            return self._peak_indices

    def _find_troughs_and_peaks(self):
        data = self.data
        time = self.time

        peaks   = []
        troughs = []

        slope = 0
        for i in range(len(data) - 1) :
            next_slope = (data[i + 1] - data[i]) / (time[i + 1] - time[i])
            if slope * next_slope < 0:
                if slope > 0:
                    peaks.append(i)
                else:
                    troughs.append(i)
            slope = next_slope

        self._peak_indices   = tuple(peaks)
        self._trough_indices = tuple(troughs)

        return

    @property
    def min_slope(self):
        '''
            Searches through data to find minimum slope, after the first trough, if there is one.
            Returns the value of the minimum slope and the index where it occurred.
        '''

        time = self.time
        data = self.data

        n = 0

        min = 2**64
        index = None
        search_start = self.peak_indices[-1] if self.peak_indices else 0# self.trough_indices[-1]# self.peak_indices[n]     if len(self.peak_indices) > 1 else 0
        search_end   = len(data) - 1# self.peak_indices[n + 1] if len(self.peak_indices) > 2  else len(data) - 1
        for i in range(search_start, search_end):
            slope = (data[i + 1] - data[i]) / (time[i + 1] - time[i])
            if slope < min:
                min = slope
                index = i
        if index == None:
            index = 0
            min = 0
        return min, index

    @property
    def peak_time(self):
        return [self.time[index] for index in self.peak_indices]

    @property
    def peak_value(self):
        return [self.data[index] for index in self.peak_indices]

    @property
    def trough_time(self):
        return [self.time[index] for index in self.trough_indices]


    @property
    def trough_value(self):
        return [self.data[index] for index in self.trough_indices]

    @property
    def normalized(self):
        try:
            return self._normalized

        except:
            data = self.data
            mean_sum = 0
            max = 0 - (2**64)
            min = 2**64

            for datum in data:
                # print(datum)
                mean_sum += datum
                if datum > max:
                    max = datum
                if datum < min:
                    min = datum
            mean = mean_sum / len(data)

            deviation_sum = 0
            for datum in data:
                deviation_sum += (datum - mean)**2
            deviation = (deviation_sum / (len(data) - 1))**0.5

            standard_min = (min - mean) / deviation
            standard_max = (max - mean) / deviation

            norm_data = []
            for datum in data:
                norm_data.append(((datum - mean) / deviation - standard_min)
                                 / (standard_max - standard_min))

            self._normalized = FrobeniusWrapper(self.time, norm_data)
            return self._normalized

    @property
    def gaussianized(self):
        try:
            return self._gaussianized

        except AttributeError:
            data = self.data

            mean = 0
            std_deviation = 40
            data_len = len(data)

            gaussian = [1 / (std_deviation * sqrt(2 * pi))
                        * e**(((x - mean)**2) / (-2 * std_deviation**2))
                        for x in range(-data_len//2, data_len//2)]
            new_data = convolve(data, gaussian)
            self._gaussianized = FrobeniusWrapper([i for i in range(len(new_data))], new_data)
            return self._gaussianized


    def __getitem__(self, item):
        if isinstance(item, int):
            return self.data[item]

        elif isinstance(item, slice):
            return FrobeniusWrapper(self.time[item], self.data[item])

def _test():
    from matplotlib import pyplot as plt
    time = [i for i in range(1000)]
    frob = FrobeniusWrapper(time, [x * x for x in range(-500, 500)])
    grib_nob = frob.gaussianized
    plt.plot(time, frob.data, 'k', [i for i in range(len(grib_nob))], grib_nob, 'r')
    plt.show()


