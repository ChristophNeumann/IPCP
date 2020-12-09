import unittest
from testinstances.ex_paper import *
import sys
import importlib.util
from model_manipulation import *
from model_information import *
from ipcp import *
from algorithm_analysis import *
from pyomo.environ import *
sys.path.append('testinstances')  # For loading models as modules
import logging

class MyTestCase(unittest.TestCase):

    def test_cutting_plane_method_for_paper_example(self):
        testinstance1 = importlib.import_module('ex_paper').model
        result = IPCP(testinstance1)
        print(result['x'])
        self.assertTrue(all(result['x']==np.array([3,3])))

    def test_reversed_ips_cuts(self):
        testinstance1 = importlib.import_module('ex_paper').model
        result = IPCP(testinstance1)
        add_reversed_IPS_cuts(testinstance1, result['model'], result['obj_pp'], result['x_pp'], False)
        model_vars = get_model_vars(testinstance1)
        coeff = get_coeff(testinstance1.reversed_cuts[1],model_vars)
        ub = testinstance1.reversed_cuts[1].upper()-get_const(testinstance1.reversed_cuts[1])
        self.assertEqual(coeff[0],0)
        self.assertEqual(2.375,ub/coeff[1])

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()