import unittest
from classes.frobenius_wrapper import FrobeniusWrapper

class FrobeniusWrapperTest(unittest.TestCase):

    def generate_analyzers(self):
        pos_line = FrobeniusWrapper(
            [10, 30, 50, 70, 90, 110],
            [20, 40, 60, 80, 100, 120]
        )
        neg_line = FrobeniusWrapper(
            [10, 30, 50, 70, 90, 110],
            [120, 100, 80, 60, 40, 20]
        )
        triangle = FrobeniusWrapper(
            [10, 30, 50, 70, 90, 110],
            [600, 500, 400, 450, 550, 650]
        )
        return pos_line, neg_line, triangle

    def test_getitem(self):
        long_frob = FrobeniusWrapper(
            [-99, -99, -99, -99, 10, 30, 50, 70, 90, 110, -99, -99, -99, -99],
            [-87, -87, -87, -87, 20, 40, 60, 80, 100, 120, -87, -87, -87, -87]
        )
        self.assertEqual(long_frob[4], 20)
        self.assertEqual(long_frob[4:10].normalized.data, (0, 0.2, 0.4, 0.6, 0.8, 1))


    def test_inflection_point(self):
        concave_switch = FrobeniusWrapper(
            [10, 30, 50, 70, 90, 110],
            [20, 60, 80, 80, 100, 140]
        )
        self.assertEqual(concave_switch.trough_index, 0)
        self.assertEqual(concave_switch.peak_index, 5)


    def test_min_slope(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.min_slope, (1, 0))
        self.assertEqual(neg_line.min_slope, (-1, 0))
        self.assertEqual(triangle.min_slope, (2.5, 2))



    def test_normalized_data(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.normalized.data, (0, 0.2, 0.4, 0.6, 0.8, 1))
        self.assertEqual(neg_line.normalized.data, (1, 0.8, 0.6, 0.4, 0.2, 0))

        # There is a rounding error for 0.2, so we use almost equal for each member individually
        for calculated, expected in zip(triangle.normalized.data, (0.8, 0.4, 0, 0.2, 0.6, 1)):
            self.assertAlmostEqual(calculated, expected, places=16)



    def test_peak_index(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.peak_index, 5)
        self.assertEqual(neg_line.peak_index, 0)
        self.assertEqual(triangle.peak_index, 5)

    def test_peak_time(self) :
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.peak_time, 110)
        self.assertEqual(neg_line.peak_time, 10)
        self.assertEqual(triangle.peak_time, 110)

    def test_peak_value(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.peak_value, 120)
        self.assertEqual(neg_line.peak_value, 120)
        self.assertEqual(triangle.peak_value, 650)



    def test_trough_index(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.trough_index, 0)
        self.assertEqual(neg_line.trough_index, 5)
        self.assertEqual(triangle.trough_index, 2)

    def test_trough_time(self) :
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.trough_time, 10)
        self.assertEqual(neg_line.trough_time, 110)
        self.assertEqual(triangle.trough_time, 50)

    def test_trough_value(self):
        pos_line, neg_line, triangle = self.generate_analyzers()
        self.assertEqual(pos_line.trough_value, 20)
        self.assertEqual(neg_line.trough_value, 20)
        self.assertEqual(triangle.trough_value, 400)
