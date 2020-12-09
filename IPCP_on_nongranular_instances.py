from algorithm_analysis import *
pd.set_option('display.max_rows', None)

save_results_to_file = False
sys.setrecursionlimit(10000) # For larger Problems
sys.path.append('minlplib_instances')  # For loading models as modules

'''Results for consistent nongranular instances'''
testbed = 'nongranular'
test_problems = read_test_instances(testbed)
res_linear = evaluate_IPCP_empty_IPS(test_problems)
print(res_linear)
if save_results_to_file:
    save_obj(res_linear, testbed + '_results')

'''Results for inconsistent instances'''
test_problems = ['ball_mk3_10',
                 'ball_mk3_20','ball_mk3_30','ball_mk4_05',
                 'ball_mk4_10', 'ball_mk4_15']

res_linear = evaluate_IPCP_empty_IPS(test_problems)
print(res_linear)
if save_results_to_file:
    save_obj(res_linear, 'inconsistent_results')