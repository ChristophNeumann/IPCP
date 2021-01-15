'''IPCP'''
mode = 'INNER' # or 'OUTER', 'CONVEX'
epsilon = 2*10 ** -6
max_iter = 2500
print_solver_output = False
LP_solver = 'cbc' #CBC will use clp for solving the LPs
NLP_solver = 'ipopt'
milp_solver='cbc'
ABS_BOUND = 1E12  # Lower bound for objective value/ variables in case of unboundedness
delta = 1 - 10 ** -4  # Enlargement parameter

'''Bonmin'''
time_limit = 1800
algorithm = 'b-ifp' #b-hyb or b-oa
preprocess = True
save_to_file = True
solver_feas_tol = 10**(-6)
