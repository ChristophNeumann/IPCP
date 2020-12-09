from pyomo.repn import generate_standard_repn
from pyomo.core.expr.visitor import identify_variables
from math_functions import *
from pyomo.core.base.symbolic import differentiate
from pyomo.opt import SolverStatus, TerminationCondition
from pyomo.environ import *
import logging
import time

int_type = ['Integers', 'PositiveIntegers', 'NonPositiveIntegers',
            'NegativeIntegers', 'NonNegativeIntegers', 'Boolean', 'Binary']


def bool_vec_is_int(m):
    """Returns an numpy-array isInt with isInt[i] = True iff variable i is some Integral variable"""
    variables = get_model_vars(m)
    isInt = np.array([str(v.domain) in int_type for v in variables])
    return isInt


def enlargement_param(coeff, is_int):
    if any(coeff[np.logical_not(is_int)]) | any(coeff[is_int] - np.floor(coeff[is_int])):
        g = 0
    else:
        g = gcd_vec(coeff[is_int])
    return g


def get_coeff(constr, model_vars):
    """Returns the vector of coefficients coeff of a linear pyomo constraint constr.
    Note that coeff[i] = 0 holds, if the variable does not appear in the constraint"""

    # assert(constr.body.polynomial_degree()==1)
    names_model_vars = [v.name for v in model_vars]  # all variablenames from the model
    coeff = np.zeros(len(model_vars))
    repn = generate_standard_repn(constr.body)
    coefs = repn.linear_coefs
    vars = repn.linear_vars
    for i, var in enumerate(vars or []):
        coeff[names_model_vars.index(var.name)] = coefs[i]
    return coeff

def get_const(constr):
    repn = generate_standard_repn(constr.body)
    return repn.constant

def is_leq_constr(constr):
    """Returns true, if and only if a pyomo constraint constr is of the form body <= upper"""
    if constr.lower is None:
        return True
    else:
        return False


def constr_value(constr):
    """Rearranges a constraint at some point x to g_i(x) <= 0 and computes g_i(x"""
    if is_leq_constr(constr):
        g_val = constr.body() - constr.upper()
    else:
        g_val = constr.lower() - constr.body()
    return g_val


def get_vars_from_constr(constr):
    """Returns all variables appearing in a constraint constr"""
    return list(identify_variables(constr.body))


def contains_eq_constrs_on_int_vars(m):
    """Checks if a model m contains some equality constraint ax+by=c with y Integer"""
    result = False
    for c in m.component_objects(Constraint, active=True):
        if c.equality:
            assert (c.body.polynomial_degree() in [0, 1])
            my_vars = get_vars_from_constr(c)
            for v in my_vars:
                if str(v.domain) in int_type:
                    result = True
                    break
    return result


def model_status(result_solver):
    if (result_solver.solver.status == SolverStatus.ok) and \
            (result_solver.solver.termination_condition == TerminationCondition.optimal):
        result = 'optimal'
    elif result_solver.solver.termination_condition == TerminationCondition.infeasible:
        result = 'infeasible'
    elif result_solver.solver.termination_condition == TerminationCondition.unbounded:
        result = 'unbounded'
    else:
        result = result_solver.solver.termination_message

    return result

def g_max_lin(m):

    linear_constrs = []
    for constr in m.component_objects(Constraint):
        if isinstance(constr,ConstraintList):
            for i in range(len(constr)):
                linear_constrs.append(constr[i+1])
        elif constr.body.polynomial_degree() in [0, 1]:
            linear_constrs.append(constr)
    # Note that equality constraints are always fulfilled, as they are assumed
    # to be linear and to only contain continuous variables and are not changed when rounding.
    if linear_constrs:
        if len(linear_constrs) >= 2:
            g_vals = np.zeros(len(linear_constrs))
            for idx, constr in enumerate(linear_constrs):
                g_vals[idx] = constr_value(constr)
            return max(g_vals)
        else:
            return constr_value(linear_constrs[0])
    else:
        return 0


def get_max_violated_constr(m):
    """Returns the maximum violated constraint of all constraint functions of a pyomo model m for a given point x,
    implicitly rearranging all constraints to g_i(x) <= 0 and computing max(g_i(x))"""

    nonlinear_constrs = []
    for constr in m.component_objects(Constraint):
        if not (constr.body.polynomial_degree() in [0, 1]):
            nonlinear_constrs.append(constr)
    # Note that equality constraints are always fulfilled, as they are assumed
    # to be linear and to only contain continuous variables and are not changed when rounding.
    if len(nonlinear_constrs) >= 2:
        g_vals = np.zeros(len(nonlinear_constrs))
        for idx, constr in enumerate(nonlinear_constrs):
            g_vals[idx] = constr_value(constr)
        idx_max = np.argmax(g_vals)
        return nonlinear_constrs[idx_max]
    else:
        return nonlinear_constrs[0]


def gradient(constr, variables):
    t = time.time()
    grad_num = np.array([value(partial) for partial in differentiate(constr.body, wrt_list=variables)])
    if (not (is_leq_constr(constr))):
        grad_num = -grad_num
    elapsed = time.time()-t
    logging.debug('gradient time is: ' + str(elapsed) + ' seconds')
    return grad_num


def var_value(model_vars):
    return np.array([value(v) for v in model_vars])


def get_model_vars(m):
    """Returns a list of all variables in a model. Distinguishes between
    indexed and non-indexed variables"""
    variable_list = []
    for v in m.component_objects(Var, active=True):
        if type(v) == pyomo.core.base.var.IndexedVar:
            for i in list(v.keys()):
                variable_list.append(v[i])
        else:
            variable_list.append(v)
    return variable_list

def get_int_vars(m):
    """Returns a list of all integer variables from the model"""
    model_vars = get_model_vars(m)
    int_vars =[v for v in model_vars if (str(v.domain) in int_type)]
    return int_vars


def ipcp_successful(g_nu, eps, solver_message, only_epi):
    result = False
    if ((constr_value(g_nu) <= eps) or only_epi) and (model_status(solver_message) == 'optimal'):
        print('IPCP successful.')
        result = True
    elif (model_status(solver_message) == 'infeasible'):
        print('Empty inner parallel set')
    elif (model_status(solver_message) == 'optimal' and constr_value(g_nu) > eps):
        print('Too many iterations')
    else:
        print('An unexpected case has occurred')
    return result


def print_cutting_planes(model):
    for i in range(1, len(model.kelley_ips_cuts._data) + 1):
        print(model.kelley_ips_cuts._data[i].expr)

def check_sense_on_minimize(model):
    return model.obj.sense == 1