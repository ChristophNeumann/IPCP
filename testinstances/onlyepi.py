from ipcp import *
model = ConcreteModel(name = "Example only epi")
model.y1 = Var(bounds=(0,3), within=Integers)
model.y2 = Var(bounds=(0,3), within=Integers)
model.obj = Objective(expr=model.y1**2 + model.y2**2)
model.g1 = Constraint(expr=model.y1+model.y2 >=1)
pp_model = preprocess_model(model)[0]
result = IPCP(pp_model, only_epi=True)
print("Objective value is ",result['obj_pp'])
print("The last iterate after the postprocessing step is: ", result['x_pp'])
print("It took", result['iterations'], "LPs to compute this point.")