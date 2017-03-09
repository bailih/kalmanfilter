import numpy as np
import matplotlib.pyplot as plt


class HealthEstimater():
    def __init__(self, pulse_data, sample_rate=0.01):
        self.sample_rate = sample_rate
        self.pulse_data = pulse_data
        self.v_init = self.pulse_data[0]

        self.current_threshold = 0.01  # miliamples
        self.targets = [[0, 40, -0.0000999],
                        [41, 60, -0.000163],
                        [61, 100, -0.0002291]]
        self.window_size = 2500

        self.pulse_start_point = self._get_pulse_start_point()
        self.pulse_end_point = self._get_pulse_end_point()

        self.pulse_start_time = self.sample_rate * self.pulse_start_point
        self.pulse_end_time = self.sample_rate * self.pulse_end_point

        self.wait_time = 30
        self.evaluation_time = self.pulse_end_time + self.wait_time
        self.model_points = []
        self.derivative_points = []
        self.derivative_point = None

        self.soh_guess = None

    def _get_pulse_end_point(self, threshold=10.0):
        point = -1
        currents = []
        offset = self.pulse_start_point + 500
        counter = offset
        for line in self.pulse_data:
            currents.append(line[2])
        for current in currents[offset:]:
            if current <= threshold:
                point = counter
                break
            counter += 1

        print('Current did not go below threshold of {} mA!'.format(threshold)) if point == -1 else None
        return point

    def _get_pulse_start_point(self, threshold=10.0):
        point = -1
        currents = []
        for line in self.pulse_data:
            currents.append(line[2])
        for current in currents:
            if current >= threshold:
                point = currents.index(current)
                break

        print('Current did not exceed threshold of {} mA!'.format(threshold)) if point == -1 else None
        return point

    def _get_model_points(self, order=1):
        window_size = self.window_size
        model_points = []
        voltage = []
        time = []
        for line in self.pulse_data[self.pulse_end_point:]:
            time.append(line[0])
            voltage.append(line[1])

        for i in range(len(time) - window_size):
            coeffs = np.polyfit(time[i:i + window_size],
                                voltage[i:i + window_size],
                                order)
            p = np.poly1d(coeffs)
            model_points.append(p(int(i + window_size / 2)))
        return model_points

    def _get_model_derivative_points(self, order=1):
        window_size = self.window_size
        der_points = []
        voltage = []
        time = []
        for line in self.pulse_data[self.pulse_end_point:]:
            time.append(line[0])
            voltage.append(line[1])

        w = 100

        smoothed_voltage = voltage[:int(w / 2)]
        for i in range(len(voltage) - w):
            avg = 0
            for j in range(w):
                avg += voltage[i + j]
            avg /= w
            smoothed_voltage.append(avg)

        smoothed_voltage += voltage[-int(w / 2):-1]
        smoothed_voltage = smoothed_voltage[:len(time)]

        for i in range(len(time) - window_size):
            coeffs = np.polyfit(time[i:i + window_size],
                                smoothed_voltage[i:i + window_size],
                                order)
            p = np.poly1d(coeffs)
            dp = np.polyder(p)
            der_points.append(dp(int(i + window_size / 2)))
        return der_points

    def _plot_derivative(self):
        time = []
        for line in self.pulse_data[self.pulse_end_point:]:
            time.append(line[0])
        plt.plot(time[int(self.window_size / 2): int(-self.window_size / 2)], self.derivative_points)
        plt.show()

    def _estimate_health(self):
        diff = 1000000
        low, high = 0, 0

        for vals in self.targets:
            if abs(self.derivative_point - vals[-1]) < diff:
                diff = abs(self.derivative_point - vals[-1])
                low = vals[0]
                high = vals[1]

        return [low, high]

    def get_health(self):
        return self.soh_guess

    def calculate_model(self):
        self.model_points = self._get_model_points()

    def calculate_model_derivatives(self):
        self.derivative_points = self._get_model_derivative_points()

    def calculate_evaluation_point(self):
        self.derivative_point = self.derivative_points[int(self.wait_time / self.sample_rate)]

    def calculate_health(self):
        self.soh_guess = self._estimate_health()


def _test():
    filename = 'C:/Users/jacob.law/Desktop/test1/all/13_3/6.txt'
    with open(filename, 'r') as content_file:
        content = content_file.read()

    lines = content.split('\n')
    lines.pop(0)  # first is header
    lines.pop(-1)
    raw_data = [line.split('\t') for line in lines]

    print("Testing HealthEstimator...")

    float_raw_data = []
    for line in raw_data:
        float_raw_data.append([float(line[0]), float(line[1]), float(line[2])])

    he = HealthEstimater(float_raw_data)
    he.calculate_model_derivatives()
    he.calculate_evaluation_point()
    he.calculate_health()
    health = he.get_health()
    print('The derivative at {}s is {}'.format(he.evaluation_time, he.derivative_point))
    print('From pulse recovery derivative, I think SoH is between {}% and {}%'.format(*health))
    he._plot_derivative()


_test()
