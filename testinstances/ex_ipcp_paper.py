import sys
sys.path.insert(0,'..')
from ipcp import *
from algorithm_analysis import *
model = ConcreteModel(name = "Example IPCP-Paper")
model.x1 = Var(bounds=(0,3), within=Reals)
model.x2 = Var(bounds=(None,None), within=Integers)
model.obj = Objective(expr=-model.x1+model.x2)
model.g1 = Constraint(expr= (model.x1-3)**2+(model.x2-3)**2<= 1/4)
model.l1 = Constraint(expr= -2/3*model.x1+ model.x2 <= 2)
model.l2 = Constraint(expr= 1/3*(model.x1) + model.x2>=2)
result = IPCP(model)
print_IPCP_results(result)




