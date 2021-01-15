from model_manipulation import *
import params


def IPCP(m, mode=params.mode, eps=params.epsilon, maxiter=params.max_iter, only_epi=False):
    """Implements the Inner Parallel Cutting Plane Method.
    Args:
    m -- Optimization model for which we want a feasible point
    v_convex -- possible modes 'INNER', 'CONVEX' and 'OUTER'. The article only discusses the "INNER" case. The
                other two modes are only experimental.
                "OUTER" corresponds to the case where we leave convex constraints in the model and solve a
                convex problem in each iteration. 'CONVEX' corresponds to the case where we do compute the inner
                parallel set of linear constraints also with cutting planes on the fly.
    eps -- Termination tolerance
    maxiter -- Maximum number of iterations
    only_epi -- True if the only source of nonlinearity stems from the reformulation of a nonlinear objective function.

    Returns:
        A dictionary with the following entries:
        'successful' -- True/False, indicating whether the method was able to generate a feasible point
        'model' -- the model which contains the last enlarged inner parallel set, including all cutting planes
        'message' -- the last solver message
        'iterations -- Number of iterations - corresponds to the number of sub LPs
        'runtime'-- cumulative runtimes of all sub LPs
        'runtime_pp' -- Cumulative runtimes of all sub LPs plus that of the postprocessing step
        'g_vals' -- the values of the maximum violated constraints of the iterates
        'obj' -- objective value of the generated feasible point before postprocessing
        'obj_pp' -- objective value of the generated feasible point after postprocessing
        'x' -- last iterate (feasible point) before postprocessing
        'x_pp' -- last iterate (feasible point) after postprocessing
        'obj_vals' -- objective values of the iterates.
        'primal_feas' -- maximum constraint violation of the generated feasible point before postprocessing
        'primal_feas_pp' -- maximum constraint violation of the generated feasible point after postprocessing
    """

    vars_master = get_model_vars(m)
    if only_epi:
        index_alpha_epi = get_alpha_epi_position(m)
    eips = build_enlarged_IPS(m, mode)
    vars_eips = get_model_vars(eips)
    unbounded_problem = False
    opt = SolverFactory(params.LP_solver)
    solver_message = opt.solve(eips)
    runtime = float(solver_message.Solver[0]['Time'])
    i = 1
    g_vals = []
    obj_vals_only_epi = []
    unsucc_epi_iter = 0

    if model_status(solver_message) == 'unbounded':
        logging.info('Add Box Constraints due to unboundedness')
        add_box_constraints(eips)
        unbounded_problem = True
        solver_message = opt.solve(eips)
        runtime += solver_message.solver.time

    if model_status(solver_message) != 'optimal':
        logging.info('First Problem is not solvable')
        successful = False

    else:
        logging.info('First Problem is solvable')
        x = var_value(vars_eips)
        is_int = bool_vec_is_int(m)
        x_nu = rounding(x, is_int)
        x_nu_epi = np.copy(x_nu)  # for evaluating if we're making progress in the objective value in epi-formulations
        set_var_vals(vars_master, x_nu, is_int)
        g_nu = get_max_violated_nonlinear_constr(m)
        if only_epi:
            obj_vals_only_epi.append(value(m.obj) + constr_value(g_nu))
            unsucc_epi_iter = 0
            x_nu_epi = np.copy(x_nu)
            x_nu_epi[index_alpha_epi] = obj_vals_only_epi[-1]
        eips.kelley_ips_cuts = ConstraintList()
        g_vals.append(constr_value(g_nu))

        while (constr_value(g_nu) > eps) and (i < maxiter) and \
                (model_status(solver_message) == 'optimal') and (unsucc_epi_iter < 20):

            i = i + 1
            logging.debug(str(constr_value(g_nu)))
            grad_num = gradient(g_nu, vars_master)
            add_ips_cut(eips.kelley_ips_cuts, g_nu, x_nu, vars_eips, grad_num, is_int)
            solver_message = opt.solve(eips, tee=params.print_solver_output)
            x = var_value(vars_eips)
            x_nu = rounding(x, is_int)
            set_var_vals(vars_master, x_nu, is_int)
            g_nu = get_max_violated_nonlinear_constr(m)
            runtime += float(solver_message.Solver[0]['Time'])
            g_vals.append(constr_value(g_nu))

            # Check if we make progress with feasible points
            if only_epi:
                logging.debug('Only epi: Appending objective value %f' % (value(m.obj) + constr_value(g_nu)))
                obj_vals_only_epi.append(value(m.obj) + constr_value(g_nu))
                unsucc_epi_iter = unsucc_epi_iter + 1
                if obj_vals_only_epi[-1] <= min(obj_vals_only_epi):
                    unsucc_epi_iter = 0
                    x_nu_epi = np.copy(x_nu)
                    x_nu_epi[index_alpha_epi] = obj_vals_only_epi[-1]
                logging.debug('Unsuccessful iter is %i' % unsucc_epi_iter)

        logging.debug('last constraint value is %f' % constr_value(g_nu))
        successful = ipcp_successful(g_nu, eps, solver_message, only_epi)

    if unbounded_problem:
        g_vals.pop(0)

    if successful:

        if only_epi:
            x_nu = x_nu_epi

        set_var_vals(vars_eips, x_nu, is_int)  # For query of optimal point and objective
        set_var_vals(vars_master, x_nu, is_int)
        obj_value = value(m.obj)
        max_constr_violation = get_max_violated_constr_value(m)

        x_pp, v_pp, runtime_PP = solve_CP_with_fixed_y(m, x_nu, only_epi)
        set_var_vals(vars_master,x_pp,is_int)
        max_constr_violation_pp = get_max_violated_constr_value(m)

        if is_epsilon_feasible_nonlinear(m):
            set_var_vals(vars_eips, x_pp, is_int)
        else:
            print('Using original feasible point because IPOPT-point is more than epsilon infeasible')
            set_var_vals(vars_master, x_nu, is_int)

        logging.info('Objective value before Post Processing: %f' % obj_value)
        logging.info('Objective value after Post Processing: %f' % v_pp)

        if (get_max_constr_value_lin(m) > 10 * eps):
            logging.warning('WARNING: Linear constraint violated, value is: ' + str(get_max_constr_value_lin(m)))

    else:
        x_pp = []
        max_constr_violation = runtime_PP = v_pp = np.inf
        obj_value = float('inf')
        x_nu = []

    return {'successful': successful, 'model': eips, 'iterations': i, 'message': solver_message, 'runtime': runtime,
            'runtime_pp': runtime_PP + runtime, 'g_vals': g_vals, 'obj': obj_value, 'obj_pp': v_pp,
            'x': x_nu, 'x_pp': x_pp, 'obj_vals': obj_vals_only_epi,
            'primal_feas' : max_constr_violation, 'primal_feas_pp': max_constr_violation_pp}


def add_reversed_IPS_cuts(original_model, result_IPCP, is_epi):

    eips = result_IPCP['model']
    v_ipcp = result_IPCP['obj_pp']
    x_ipcp = result_IPCP['x_pp']
    if (result_IPCP['primal_feas_pp']>params.epsilon) & (result_IPCP['primal_feas']<result_IPCP['primal_feas_pp']):
        print("Using point before post processing for reversed inner parallel cuts")
        v_ipcp = result_IPCP['obj']
        x_ipcp = result_IPCP['x']

    original_vars = get_model_vars(original_model)
    is_int = bool_vec_is_int(original_model)
    set_var_vals(original_vars, x_ipcp, is_int)

    g_max = get_max_violated_nonlinear_constr(original_model)
    grad_num = gradient(g_max, original_vars)
    if any(grad_num):
        original_model.optimal_point_constr = Constraint(expr=sum([(x_i - x_ipcp[idx]) * grad_num[idx] for idx, x_i in
                                                                   enumerate(original_vars)]) + constr_value(
            g_max) <= 0)
    added_cuts = eips.kelley_ips_cuts
    original_model.reversed_cuts = ConstraintList()

    eips_vars = get_model_vars(eips)
    substitution_map = {}
    var_size = len(is_int)
    n_of_cuts = len(added_cuts)
    first_cut = 0
    if var_size < n_of_cuts:
        first_cut = len(added_cuts) - var_size

    for i in range(first_cut, n_of_cuts):
        coeff = get_coeff(added_cuts[i + 1], eips_vars)
        beta = np.linalg.norm(coeff[is_int], ord=1)
        current_cut = added_cuts[i + 1]
        logging.debug(current_cut.body)
        logging.debug(current_cut.upper)
        logging.debug(str(beta))
        for j, y_i in enumerate(eips_vars):
            assert (y_i.name == original_vars[j].name)
            substitution_map[id(y_i)] = original_vars[j]
        original_model.reversed_cuts.add(
            (current_cut.lower, clone_expression(current_cut.body, substitution_map), current_cut.upper + beta / 2))
        assert not (current_cut.lower)

    if constr_value(g_max)<=params.solver_feas_tol: #Set an objective bound if solver considers the point to be feasible
        if is_epi:
            original_model.alpha_epi.setub(v_ipcp + 10*params.epsilon)
            logging.info("set ub of epi var to" + str(v_ipcp + 10*params.epsilon))
        else:
            original_model.objective_constraint = Constraint(expr=original_model.obj.expr <= v_ipcp + 10*params.epsilon)


def solve_CP_with_fixed_y(original_model, x_nu, only_epi=False):

    is_int = bool_vec_is_int(original_model)
    y_nu = x_nu[is_int]
    post_processing_model = original_model.clone()
    vars_pp = get_model_vars(post_processing_model)
    set_var_vals(vars_pp, x_nu, is_int)
    if only_epi:
        con_epi = get_max_violated_nonlinear_constr(post_processing_model)
    fix_int_variables(post_processing_model, y_nu)
    cont_relax_int_vars(vars_pp, enlarge_box_constrs=False)
    opt = SolverFactory(params.NLP_solver)
    opt.options['warm_start_init_point'] = "yes"
    opt.options['constr_viol_tol'] = 1E-8
    solver_message = opt.solve(post_processing_model)
    if model_status(solver_message) == 'optimal':
        x_pp = var_value(vars_pp)
        if only_epi:
            updated_alpha_epi_value = value(post_processing_model.obj) + constr_value(con_epi)
            x_pp[get_alpha_epi_position(post_processing_model)] = updated_alpha_epi_value
        set_var_vals(vars_pp,x_pp,is_int)
        v_pp = value(post_processing_model.obj)
        if not is_epsilon_feasible(post_processing_model):
            primal_feas = get_max_violated_constr_value(post_processing_model)
            print('Ipopt exceeds IPCPs feasibility tolerance by %.7f.'%primal_feas)
    else:
        x_pp = x_nu
        v_pp = value(post_processing_model.obj)
        print('Ipopt did not converge with an optimal point')

    runtime_PP = float(solver_message.Solver[0]['Time'])
    # pp_constr = get_max_violated_constr(post_processing_model)
    # eips.ips_cuts.add(constr_value())
    return x_pp, v_pp, runtime_PP

def compute_RIPCs_lower_bound(preprocessed_model, result_IPCP, is_epi):
    lb_model = preprocessed_model.clone()
    deactivate_nonlinear_constrs(lb_model)
    opt = SolverFactory(params.milp_solver)
    solver_message = opt.solve(lb_model)
    if model_status(solver_message) == 'optimal':
        lower_bound_initial = value(lb_model.obj)
    else:
        lower_bound_initial = -np.inf

    add_reversed_IPS_cuts(lb_model, result_IPCP, is_epi)
    opt = SolverFactory(params.milp_solver)
    solver_message = opt.solve(lb_model)
    runtime_lb = float(solver_message.Solver[0]['Time'])
    if model_status(solver_message) == 'optimal':
        lower_bound = value(lb_model.obj)
    else:
        lower_bound = -np.inf
    print("Difference in bounds is",lower_bound-lower_bound_initial)
    return lower_bound, lower_bound_initial, runtime_lb
