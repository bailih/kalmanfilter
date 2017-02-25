from filters.unscented import Unscented
import numpy as np
import unittest

# Constants for testing, so we don't break the tests if we adjust the real values.
ALPHA = 0.01
KI    = 0.0
BETA  = 2.0

class UnscentedTest(unittest.TestCase):

    def gen_filter(self, order=1):
        q = np.mat(np.zeros((order + 1, order + 1)))
        p = np.mat(np.zeros((order + 1, order + 1)))
        r = 5
        x = np.mat(np.zeros((order + 1, 1)))
        for i in range(order + 1):
            q[i, i] = 10**-7 + 10**-7 * i
            p[i, i] = 1 + 1 * i
            x[i, 0] = 0.01 + 0.01 * i
        x[order, 0] = 0.5

        param = {
            'dt'      : 1,
            'rpos'    : 0.5,    # Discharging resistance
            'rneg'    : 0.6,    # Charging resistance
            'cnom'    : 10000,
            'rbranch' : [0.1 + 0.1 * i for i in range(order)],
            'cbranch' : [1000 + 1000 * i for i in range(order)],
            'eta'     : 1
        }

        c = (110, -100, 90, -80, 70, -60, 50, -40, 30, -20, 10)

        unsc = Unscented(order, q, p, r, x, param, c)
        unsc.c = (ALPHA * ALPHA * (order + 1 + KI))
        unsc.lamb = unsc.c - (order + 1)
        return unsc


    def test_estimated_voltage(self):
        unsc = self.gen_filter()
        # set soc to 0. OCV(soc) is tested separately
        unsc.X[1, 0] = 0
        for current, expected in zip((1, 0, -1), (10.59, 9.99, 9.49)):
            self.assertEqual(unsc.get_estimated_voltage(current), expected)


    def test_ocv(self):
        unsc = self.gen_filter()
        for soc, expected in zip((0, 0.5, 1, 1.5), (10, 4.482421875, 60, 85)):
            self.assertEqual(unsc._OCV(soc), expected)


    def test_priori_update(self):
        # Tests using a specific update function to avoid breaking
        unsc = self.gen_filter()
        unsc._discrete_volts = unsc._volts_inv_i_minus_a_t
        sigmas = [np.mat([[0], [0.5]]), np.mat([[0.01], [0.1]])]
        current = 1
        dt = 20
        p = unsc.p
        updated_sigmas = unsc._priori_update(sigmas, current, dt, p)
        expected_sigmas = [np.mat([[0.02], [0.498]]), np.mat([[0.02833333], [0.098]])]
        for i, (updated, expected) in enumerate(zip(updated_sigmas, expected_sigmas)):
            for j in range(2):
                self.assertAlmostEqual(updated[j, 0], expected[j, 0], places=7,
                                       msg='Point: {}\tRow: {}'.format(i, j))


    def test_sigma_points(self):
        unsc = self.gen_filter()
        mean = unsc.X
        cov  = unsc.P
        got_list = unsc._calculate_sigma_points(mean, cov)
        expected_list = [np.mat([[0.01], [0.5]]),
                         np.mat([[(1 + 2**0.5) * 10**-2], [0.5]]),
                         np.mat([[0.01], [0.52]]),
                         np.mat([[(1 - 2**0.5) * 10**-2], [0.5]]),
                         np.mat([[0.01], [0.48]]),
                        ]
        for i, (got, expected) in enumerate(zip(got_list, expected_list)):
            for j in range(2):
                self.assertAlmostEqual(got[j, 0], expected[j, 0], places=12,
                                       msg='Point: {}\tRow: {}'.format(i, j))

    def test_sigma_negative_cov(self):
        # Since the actual sigma calculations are tested elsewhere, to prevent this
        # test from breaking, it only makes sure the function actually completes.
        unsc = self.gen_filter()
        mean = unsc.X
        cov = np.mat([[-1, 0],
                      [0, 2]])
        got_list = unsc._calculate_sigma_points(mean, cov)

    def test_soc(self):
        for unsc, order in ((self.gen_filter(order), order) for order in (1, 3, 25)):
            self.assertEqual(unsc.get_soc(), 0.5,
                             'get_soc() failed for order {}'.format(order))
