from ipcp import *
from algorithm_analysis import *
model = ConcreteModel(name = "Example 2")

model.x1 = Var(bounds=(-1e3,1e3), within=Reals)
model.x2 = Var(bounds=(-1e3,1e3), within=Integers)
model.obj = Objective(expr=-model.x1-model.x2)
model.g1 = Constraint(expr=model.x1**4+model.x2**4-8 <=0)
result = IPCP(model)
print(result['obj_pp'])
print(result['x_pp'])

