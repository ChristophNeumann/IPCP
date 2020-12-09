import importlib.util
import sys
import csv
import importlib.util
import logging
from ipcp import *
from matplotlib import pyplot as plt
import os
import pickle
import pandas as pd
from matplotlib import rc
sys.path.append('minlplib_instances')  # For loading models as modules
import params
#sys.path.append('problem_information')

def save_obj(obj, name ):
    with open('solver_log/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('solver_log/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def read_optimal_values(test_problems):
    with open('problem_information/primal_bounds.csv',
              mode='r') as infile:
        reader = csv.reader(infile)
        mydict = {rows[0]: rows[1] for rows in reader}
        opt_vals = [mydict[item] for item in test_problems]
    return opt_vals

def read_only_epi_data(test_problems):
    with open('problem_information/only_epi.csv',
              mode='r') as infile:
        reader = csv.reader(infile)
        mydict = {rows[0]: rows[1] for rows in reader}
        only_epi = [mydict[item] == 'True' for item in test_problems]
        return only_epi

def get_data_matrix(test_problems, check_epi = True, opt_vals = True):
    print('############################')
    print('Extracting problem data table. Performing Epigraph reformulations on some instances. '
          'This may take some time')
    print('############################')
    data_matrix = []
    names = []
    for idx, name in enumerate(test_problems):
        current_model = load_pyomo_model(name)
        if check_epi & (not objective_is_linear(current_model)):
            name = name + "$^*$"
            current_model, _, _ = preprocess_model(current_model,printout=False)
        names.append(name)
        datalist = get_model_data_for_print(current_model)
        data_matrix.append(datalist)
    pd_data_matrix = pd.DataFrame(data_matrix, columns = ['vars','constrs'])
    pd_data_matrix.index = names
    if opt_vals:
        best_known = pd.DataFrame(read_optimal_values(test_problems), columns = ['$v^\star$'])
        best_known.index = names
        pd_data_matrix = pd.concat([pd_data_matrix,best_known],axis=1)
    return pd_data_matrix

def evaluate_IPCP_empty_IPS(test_problems, mode ='INNER', epsilon = params.epsilon, max_iter = 1000):

    result_list = []
    result_data = get_data_matrix(test_problems, opt_vals = False)
    for idx, name in enumerate(test_problems):
        print('############################')
        print('Testing problem ', name)
        current_model = load_pyomo_model(name)
        # Get model Data for result file later
        current_model, red_successful, is_epi = preprocess_model(current_model)

        if red_successful:
            result = IPCP(current_model, mode, epsilon, max_iter)
            print('Objective value: ', result['obj'])
            print('Iterations: ', result['iterations'])
            print('Runtime: ', result['runtime'])
            result_list.append([result['iterations'],result['runtime']])
            print('############################ \n')

    results_FCPA = pd.DataFrame(np.array(result_list), columns=['iter', 'time'])
    results_FCPA.index = result_data.index
    result_matrix = pd.concat([result_data,results_FCPA],axis=1)
    return result_matrix

def run_IPCP(test_problems):

    result_matrix = []
    only_epi = read_only_epi_data(test_problems)

    for idx, name in enumerate(test_problems):

        print('############################')
        print('Testing problem ', name)
        original_model = load_pyomo_model(name)
        current_model = original_model.clone()
        current_model, red_successful, is_epi = preprocess_model(current_model)

        logging.info('Only Epi: ' + str(only_epi[idx]))

        if red_successful:
            result = IPCP(current_model, params.mode, params.epsilon, params.max_iter, only_epi[idx])
            print_IPCP_results(result)
            del result['model']
            result_matrix.append([result['obj_pp'], result['obj'], result['iterations'],
                                 result['runtime_pp'], result['runtime']])
            del original_model
            del current_model
            del result

    result_matrix = pd.DataFrame(np.array(result_matrix), columns=['obj_pp', 'obj', 'iter', 'time_pp', 'time'])
    result_matrix.index = test_problems
    result_matrix['iter'] = result_matrix['iter'].astype('int')
    return result_matrix

def get_model_data_for_print(m):
    p = len(list(m.component_objects(Constraint)))
    p_nl = number_nonlinear_constrs(m)
    is_int = bool_vec_is_int(m)
    m = np.sum(is_int)
    return [str(len(is_int)) + '(' + str(m)  + ')',str(p) + '(' + str(p_nl) + ')']

def load_pyomo_model(problem_name):
    testinstance = importlib.import_module(problem_name)
    return testinstance.m

def plot_iterations(model_name):
    rc('text', usetex=True)
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)
    sys.path.append('minlplib_instances')
    m = load_pyomo_model(problem_name=model_name)
    result = IPCP(m)
    plt.plot(result['g_vals'],color = 'gray', linewidth=3.0)
    plt.xlabel('\huge $k$')
    plt.ylabel('\huge $g_{\ell_{k}}(\check x^k, \check y^k)$')
    plt.savefig(os.path.join('Z:\\hg2412\\Research\\03_FeasibleCuttingPlaneAlgorithm\\Latex', 'plot_' + model_name +'.pdf'), bbox_inches='tight')


def to_str(f,mode = 'float'):
    if mode == 'float' and f != "$*$":
        result = str("%0.2f" %float(f))
    else:
        result = str(f)
    return result

def read_test_instances(filename):
    '''Read test instance Data '''
    dir = os.path.dirname(__file__)
    problem_data = os.path.join(dir, 'problem_information', filename + '.txt')
    f = open(problem_data, 'r')
    test_problems = f.read().splitlines()
    f.close()
    return test_problems

def print_IPCP_results(result_IPCP):
    print('Iterations: ', result_IPCP['iterations'])
    print('Objective value: ', result_IPCP['obj'])
    print('Objective after PP: ', result_IPCP['obj_pp'])
    print('Run time: ', result_IPCP['runtime'])
    print('Run time after PP: ', result_IPCP['runtime_pp'])
    print('############################','\n')



def write_output_to_file(name):
    f = open('solver_log/' + name + '.txt', 'w')
    sys.stdout = f


