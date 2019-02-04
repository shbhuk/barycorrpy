import unittest
from .barycorrpy import get_BC_vel , exposure_meter_BC_vel


def fun(x):
    return x + 1

class MyTest(unittest.TestCase):
    def test(self):
        self.assertEqual(fun(3), 4)



if __name__ == '__main__':
    unittest.main()
