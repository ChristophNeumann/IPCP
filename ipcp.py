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
        'runtime'-- cumulative runtimes of all sub LPs
        'runtime_pp' -- Cumulative runtimes of all sub LPs plus that of the postprocessing step
        'g_vals' -- the values of the maximum violated constraints of the iterates
        'obj' -- objective value of the generated feasible point before postprocessing
        'obj_pp' -- objective value of the generated feasible point after postprocessing
        'x' -- last iterate (feasible point) before postprocessing
        'x_pp' -- last iterate (feasible point) after postprocessing
        'obj_vals' -- objective values, only relevant in the "only-epi" mode.}
    """

    vars_master = get_model_vars(m)
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
        x_nu_epi = x_nu  # for evaluating if we're making progress in the objective value in epi-formulations
        set_var_vals(vars_master, x_nu, is_int)
        g_nu = get_max_violated_constr(m)
        eips.kelley_ips_cuts = ConstraintList()
        g_vals.append(constr_value(g_nu))

        while (constr_value(g_nu) > eps) and (i < maxiter) and \
                (model_status(solver_message) == 'optimal') and (unsucc_epi_iter < 20):

            i = i + 1
            logging.debug(str(constr_value(g_nu)))
            grad_num = gradient(g_nu, vars_master)
            add_ips_cut(eips.kelley_ips_cuts, g_nu, x_nu, vars_eips, grad_num, is_int)
            solver_message = opt.solve(eips, tee=params.print_solver_output)
            logging.info(str(x))
            x = var_value(vars_eips)
            x_nu = rounding(x, is_int)
            set_var_vals(vars_master, x_nu, is_int)
            g_nu = get_max_violated_constr(m)
            runtime += float(solver_message.Solver[0]['Time'])
            g_vals.append(constr_value(g_nu))

            # Check if we make progress with feasible points
            if only_epi:
                logging.debug('Only epi: Appending objective value %f' % (value(m.obj) + constr_value(g_nu)))
                obj_vals_only_epi.append(value(m.obj) + constr_value(g_nu))
                unsucc_epi_iter = unsucc_epi_iter + 1
                if obj_vals_only_epi[-1] <= min(obj_vals_only_epi):
                    unsucc_epi_iter = 0
                    x_nu_epi = x_nu
                logging.debug('Unsuccessful iter is %i' % unsucc_epi_iter)

        logging.debug('last constraint value is %f' % constr_value(g_nu))
        successful = ipcp_successful(g_nu, eps, solver_message, only_epi)

    if unbounded_problem:
        g_vals.pop(0)

    if successful:
        obj_value = value(eips.obj)
        if only_epi:
            obj_value = min(obj_vals_only_epi)
            x_nu = x_nu_epi
        set_var_vals(vars_eips, x_nu, is_int)  # For query of optimal point and objective
        x_PP, v_PP, runtime_PP = solve_CP_with_fixed_y(m, x_nu)
        set_var_vals(vars_eips, x_PP, is_int)
        if v_PP != value(eips.obj):
            logging.warning('Objective value of post processing model does not match that of the feasible point')

        logging.info('Objective value before Post Processing: %f' % obj_value)
        logging.info('Objective value after Post Processing: %f' % v_PP)

        if (g_max_lin(m) > 10 * eps):
            print('WARNING: Linear constraint violated, value is: ' + str(g_max_lin(m)))
    else:
        x_PP = []
        runtime_PP = v_PP = np.inf
        obj_value = float('inf')
        x_nu = []

    return {'successful': successful, 'model': eips, 'iterations': i, 'message': solver_message, 'runtime': runtime,
            'runtime_pp': runtime_PP + runtime, 'g_vals': g_vals, 'obj': obj_value, 'obj_pp': v_PP,
            'x': x_nu, 'x_pp': x_PP, 'obj_vals': obj_vals_only_epi}


def add_reversed_IPS_cuts(original_model, eips, v_ipcp, x_ipcp, is_epi):
    original_vars = get_model_vars(original_model)
    is_int = bool_vec_is_int(original_model)

    set_var_vals(original_vars, x_ipcp, is_int)
    g_max = get_max_violated_constr(original_model)
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

    print("Value of maximum constraint violation of linear constraints:",g_max_lin(original_model))
    print("Value of maximum constraint violation of the nonlinear constraints:", constr_value(g_max))

    if is_epi:
        original_model.alpha_epi.setub(v_ipcp)
        logging.info("set ub of epi var to" + str(v_ipcp))
    else:
        original_model.objective_constraint = Constraint(expr=original_model.obj.expr <= v_ipcp)


def solve_CP_with_fixed_y(original_model, x_nu):
    is_int = bool_vec_is_int(original_model)
    y_nu = x_nu[is_int]
    post_processing_model = original_model.clone()
    vars = get_model_vars(post_processing_model)
    set_var_vals(vars, x_nu, is_int)
    fix_int_variables(post_processing_model, y_nu)
    cont_relax_int_vars(vars, enlarge_box_constrs=False)
    opt = SolverFactory(params.NLP_solver)
    opt.options['warm_start_init_point'] = "yes"
    opt.options['acceptable_tol'] = 0.1
    solver_message = opt.solve(post_processing_model)
    x_PP = var_value(vars)
    v_PP = value(post_processing_model.obj)
    runtime_PP = float(solver_message.Solver[0]['Time'])
    # pp_constr = get_max_violated_constr(post_processing_model)
    # eips.ips_cuts.add(constr_value())
    return x_PP, v_PP, runtime_PP
