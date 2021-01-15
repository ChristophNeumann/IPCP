from algorithm_analysis import *
import logging
from datetime import datetime
sys.setrecursionlimit(10000) # For larger Problems
logging.basicConfig(level=logging.INFO)
testbed = 'granular_instances'
#write_output_to_file('IPCP_debug')
test_problems = read_test_instances(testbed)
problem_data = get_data_matrix(test_problems)
result_matrix_FCPA = run_IPCP(test_problems)
result_matrix_FCPA.index = problem_data.index

results_overview = pd.concat([problem_data, result_matrix_FCPA], axis=1)
cols = results_overview.columns
results_overview[cols[2]] = results_overview[cols[2]].apply(pd.to_numeric, errors='coerce')
print(results_overview)
save_obj(results_overview, testbed + '_results_IPCP')
