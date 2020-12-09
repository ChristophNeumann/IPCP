import sys
sys.path.insert(0,'..')
from ipcp import *
model = ConcreteModel(name = "Example extended supporting hyperplane")
model.x1 = Var(bounds=(1,20), within=Reals)
model.x2 = Var(bounds=(1,20), within=Integers)
model.obj = Objective(expr=-model.x1-model.x2)
model.g1 = Constraint(expr=0.15 * ((model.x1 - 8) ** 2) + 0.1 * ((model.x2 - 6) ** 2) + 0.025 * exp(model.x1) * (
                       (model.x2) ** (-2)) - 5 <= 0)
model.g2 = Constraint(
    expr=(model.x1) ** (-1) + (model.x2) ** (-1) - (model.x1) ** (0.5) * (model.x2) ** (0.5) + 4 <= 0)
model.l1 = Constraint(expr=2 * (model.x1) - 3 * (model.x2) -2<=0)
result = IPCP(model)
print("Objective value is ",result['obj_pp'])
print("The last iterate after the postprocessing step is: ", result['x_pp'])
print("It took", result['iterations'], "LPs to compute this point.")
