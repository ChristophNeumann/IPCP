from ipcp import *
from algorithm_analysis import *
model = ConcreteModel(name = "Example 1")

model.x1 = Var(bounds=(-1e3,1e3), within=Reals)
model.x2 = Var(bounds=(-1e3,1e3), within=Integers)

model.obj = Objective(expr=-model.x1-model.x2)

model.g1 = Constraint(expr=model.x1**4+model.x2**4-8 <=0)


result = IPCP(model)
print(result['obj'])
print(result['x'])
opt = SolverFactory('scipampl')
solver_message = opt.solve(model, tee=True)

