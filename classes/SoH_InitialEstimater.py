import numpy as np
import matplotlib.pyplot as plt


class HealthEstimater():
    def __init__(self, pulse_data, sample_rate):
        self.soh_guess = 0
        self.sample_rate = sample_rate
        self.pulse_data = pulse_data
        self.v_init = self.pulse_data[0]

        self.pulse_start_time = 0
        self.pulse_end_time = 0

        self.sample_point = self.pulse_end_time + 30
        self.derivative = []
        self.derivative_point = 0

    def get_model_points(self, window_size=2500, order=1):
        model_points = []
        for i in range(len(self.pulse_data) - window_size)
            coeffs = np.polyfit()
            p = np.poly1d()
            model_points.append(p(i))
