import numpy as np
import matplotlib.pyplot as plt


class HealthEstimater():
    def __init__(self, pulse_data, sample_rate=0.01):
        self.sample_rate = sample_rate
        self.pulse_data = pulse_data
        self.v_init = self.pulse_data[0]

        self.current_threshold = 0.01  # miliamples
        self.targets = [[0, 40, ],
                        [41, 50, ],
                        [51, 100, ]]

        self.pulse_start_point = self._get_pulse_start_point()
        self.pulse_end_point = self._get_pulse_end_point()

        self.pulse_start_time = self.pulse_data[0][self.pulse_start_point]
        self.pulse_end_time = self.pulse_data[self.pulse_end_point]

        self.evaluation_time = self.pulse_end_time + 30
        self.model_points = self._get_model_points()
        self.derivative_points = self._get_model_derivative_points()
        self.derivative_point = self.derivative_points[self.pulse_data.index(self.evaluation_time)]

        self.soh_guess = self._estimate_health()

    def _get_pulse_end_point(self, threshold=0.01):
        point = -1
        for voltage in self.pulse_data[1][self.pulse_start_point:]:
            if voltage <= threshold:
                point = self.pulse_data[1].index(voltage)

        print('Current did not go below threshold of {} mA!'.format(threshold)) if point == -1 else None
        return point

    def _get_pulse_start_point(self, threshold=0.01):
        point = -1
        for voltage in self.pulse_data[1]:
            if voltage >= threshold:
                point = self.pulse_data[1].index(voltage)

        print('Current did not exceed threshold of {} mA!'.format(threshold)) if point == -1 else None
        return point

    def _get_model_points(self, window_size=2500, order=1):
        model_points = []
        for i in range(len(self.pulse_data) - window_size):
            coeffs = np.polyfit(self.pulse_data[0][self.pulse_start_point:],
                                self.pulse_data[1][self.pulse_start_point:], order)
            p = np.poly1d(coeffs)
            model_points.append(p(int(i + window_size / 2)))
        return model_points

    def _get_model_derivative_points(self, window_size=2500, order=1):
        der_points = []
        for i in range(len(self.model_points) - window_size):
            coeffs = np.polyfit(self.model_points[self.pulse_start_point:], self.model_points[self.pulse_start_point:])
            p = np.poly1d(coeffs)
            dp = np.polyder(p, )
            der_points.append(dp(int(i + window_size / 2)))
        return der_points

    def _plot_derivative(self):
        plt.plot([i for i in self.pulse_data[0][self.pulse_start_point:]], self.derivative_points)
        plt.show()

    def _estimate_health(self):
        diff = 1000000
        low, high = 0, 0

        for vals in self.targets:
            if abs(self.derivative_point - vals[-1]) < diff:
                diff = abs(self.derivative_point - vals[-1])
                low = vals[0]
                high = vals[1]

        return low, high

    def get_health(self):
        return self.soh_guess

    def _test(self):
        print("Testing HealthEstimator...")
