from algorithm_analysis import *
import logging
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
from datetime import datetime

sys.setrecursionlimit(10000) # For larger Problems
logging.basicConfig(level=logging.WARNING)
testbed = 'granular_instances'
test_problems = read_test_instances(testbed)

problem_data = get_data_matrix(test_problems)
result_matrix_FCPA = run_IPCP(test_problems)
result_matrix_FCPA.index = problem_data.index

results_overview = pd.concat([problem_data, result_matrix_FCPA], axis=1)
print(results_overview)
#save_obj(results_overview, testbed + '_results_IPCP')
