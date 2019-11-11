# -*- coding: utf-8 -*-

import unittest
import freqerica.operator.util


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_absolute_truth_and_meaning(self):
        f = freqerica.operator.ham.construct_exact_ham
        assert True


if __name__ == '__main__':
    unittest.main()
