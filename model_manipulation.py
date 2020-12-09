from model_information import *
from pyomo.core.expr.visitor import clone_expression
import params


ABS_BOUND = 1E12  # Lower bound for objective value/ variables in case of unboundedness
delta = 1 - 10 ** -4  # Enlargement parameter


def add_objective_bound(m):
    m.obj_con = Constraint(expr=m.obj.expr >= -ABS_BOUND)


def build_enlarged_IPS(m, mode='INNER'):
    eips = m.clone()
    deactivate_nonlinear_constrs(eips)
    constrs = [c for c in eips.component_objects(Constraint, active=True)]
    model_vars = get_model_vars(eips)
    is_int = bool_vec_is_int(m)
    for constr in constrs:
        coeff = get_coeff(constr, model_vars)
        beta = np.linalg.norm(coeff[is_int], ord=1)
        g = enlargement_param(coeff, is_int)
        if not constr.equality:
            if is_leq_constr(constr):
                # modify the constraint to obtain the inner parallel set of the linear constraints
                if mode != 'OUTER':
                    constr.set_value(constr.body <= floor_g(constr.upper(), g) - 1 / 2 * beta + delta * g)
                else:
                    # Version where the IPS of linear constrs is build due to cutting planes
                    constr.set_value(constr.body <= floor_g(constr.upper(), g) + delta * g)
            else:
                if mode != 'OUTER':
                    constr.set_value(-constr.body <= floor_g(-constr.lower(), g) - 1 / 2 * beta + delta * g)
                else:
                    # Version where the IPS is build due to cutting planes
                    constr.set_value(-constr.body <= floor_g(-constr.lower(), g) + delta * g)
    cont_relax_int_vars(model_vars, enlarge_box_constrs=True)
    if mode == 'CONVEX':
        activate_nonlinear_constrs(eips)
    return eips


def box_constrs_to_expr(m, in_vars):
    for var in in_vars:
        if var.bounds[0] is not None:
            m.add_component(var.name + "_lb", Constraint(expr=var >= var.bounds[0]))
        if var.bounds[1] is not None:
            m.add_component(var.name + "_ub", Constraint(expr=var <= var.bounds[1]))


def enlarge_box_constraints(var):
    if (var.bounds is not None) and (var.bounds[0] is not None):
        var.setlb(math.ceil(var.bounds[0]) - delta / 2)
    if (var.bounds is not None) and (var.bounds[1] is not None):
        var.setub(math.floor(var.bounds[1]) + delta / 2)


def cont_relax_int_vars(model_vars, enlarge_box_constrs=True):
    for var in model_vars:
        if str(var.domain) in int_type:
            if enlarge_box_constrs:
                enlarge_box_constraints(var)
            var.domain = Reals


def deactivate_constraints(constrs, m):
    for c in constrs:
        m.del_component(c)


def substitute_variables(m, substitution_map, old_int_vars):
    m.obj = clone_expression(m.obj.expr, substitute=substitution_map)
    for c in m.component_data_objects(Constraint, active=True):
        c.set_value((c.lower, clone_expression(c.body, substitute=substitution_map), c.upper))
    for v in old_int_vars:
        m.del_component(v)


def number_nonlinear_constrs(m):
    result = 0
    for constr in m.component_objects(Constraint, active=True):
        if not (constr.body.polynomial_degree() in [0, 1]):
            result = result + 1
    return result


def deactivate_nonlinear_constrs(m):
    """Deactivates all nonlinear constraints of a pyomo model m"""
    for constr in m.component_objects(Constraint, active=True):
        if not (constr.body.polynomial_degree() in [0, 1]):
            constr.deactivate()


def activate_nonlinear_constrs(m):
    """Activates all constraints of a pyomo model m"""
    for constr in m.component_objects(Constraint):
        if not (constr.body.polynomial_degree() in [0, 1]):
            constr.activate()


def add_box_constraints(m):
    for var in get_model_vars(m):
        if (var.bounds is None) or ((var.bounds[0] == None) and (var.bounds[1] == None)):
            logging.info('Adding lower and upper bounds for variable ' + str(var.name))
            var.setlb(-ABS_BOUND)
            var.setub(ABS_BOUND)
        elif var.bounds[0] is None:
            var.setlb(-ABS_BOUND)
            logging.info('Adding lower bound for variable ' + str(var.name))
        elif var.bounds[1] is None:
            var.setub(ABS_BOUND)
            logging.info('Adding upper bound for variable ' + str(var.name))


def set_var_vals(var_list, value_list, is_int):
    for idx, var in enumerate(var_list):
        if (is_int[idx]):
            var.set_value(int(value_list[idx]))
        else:
            var.set_value(value_list[idx])


def add_ips_cut(cut_list, g_nu, x_nu, vars_model, grad_num, is_int):
    beta = np.linalg.norm(grad_num[is_int], ord=1)
    cut_list.add(sum(
        [(x_i - x_nu[idx]) * grad_num[idx] for idx, x_i in enumerate(vars_model)]) + beta / 2 + constr_value(
        g_nu) <= 0)


def invert_sense(model):
    model.obj.sense = model.obj.sense * (-1)
    model.obj.expr = model.obj.expr * (-1)


def objective_is_linear(model):
    if (model.obj.expr.polynomial_degree() in [0, 1]):
        result = True
    else:
        result = False
    return result


def epgraph_reformulation(model, printout=True):
    model_lb_epi = model.clone()
    model_vars = get_model_vars(model_lb_epi)
    cont_relax_int_vars(model_vars)
    opt = SolverFactory(params.NLP_solver)
    try:
        # get lower bound for alpha
        solver_message = opt.solve(model_lb_epi)
        lb = value(model_lb_epi.obj)
        if printout:
            print('runtime for lower bound epi: ' + str(solver_message.solver.time))
            print('lower bound is: ' + str(lb))
    except:
        print('cannot compute lower bound with ipopt')
        lb = -ABS_BOUND
    epi_model = model.clone()
    epi_model.alpha_epi = Var(within=Reals, bounds=(lb, ABS_BOUND))
    epi_model.con_epi = Constraint(expr=epi_model.obj.expr <= epi_model.alpha_epi)
    epi_model.del_component('obj')
    epi_model.obj = Objective(expr=epi_model.alpha_epi, sense=minimize)
    return epi_model


def preprocess_model(model, printout=True):
    red_successful = True
    is_epi = False
    if (contains_eq_constrs_on_int_vars(model)):
#        red_successful = reduce_model(model)
        logging.error('Instance contains equality constraints on integer variables.')
    if (not check_sense_on_minimize(model)):
        invert_sense(model)
    if (not objective_is_linear(model)):
        model = epgraph_reformulation(model,printout)
        is_epi = True
        if printout:
            print('Performing IPCP on epigraph reformulation')
    return model, red_successful, is_epi


def fix_int_variables(model, y_k):
    y = get_int_vars(model)
    for i in range(len(y)):
        y[i].fix(y_k[i])
