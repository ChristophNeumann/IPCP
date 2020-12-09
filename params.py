'''IPCP'''
mode = 'INNER' # or 'OUTER', 'CONVEX'
epsilon = 2*10 ** -6
max_iter = 1500
print_solver_output = False
LP_solver = 'cbc' #CBC will use clp for solving the LPs
NLP_solver = 'ipopt'
