# Inner Parallel Cutting Plane Algorithm
> This project implements the Inner Parallel Cutting Plane Method (IPCP) for convex MINLPs which is outlined in the article [Generating feasible points for mixed-integer convex optimization problems by inner parallel cuts](http://www.optimization-online.org/DB_HTML/2018/11/6947.html).

## Table of contents
* [General info](#general-info)
* [Setup](#setup)
* [Code Examples](#code-examples)
* [Reproducing the experiments](#reproducing-the-experiments)
* [Contact](#contact)

## General info
This project is a prototype of the IPCP that can be applied to Pyomo models. It is mainly intended for researchers to experiment and compare computational results on their MICP instances and to enable an easy reproducibility of the computational results of the above article. To enable reproducibility, we provide the minlplib instances as pyomo models that were used for the numerical comparison.

As a gerenal note, under a Gams license, you may convert any (other) Gams-model to a Pyomo model and then apply the IPCP. This process is described [here](https://www.gams.com/latest/docs/S_CONVERT.html). 

## Setup
The method has been tested under Python 3.7 using Pyomo 5.7 as the main framework. For solving the sub-LPs and the NLP outlined in the postprocessing step, it needs an LP-solver and an NLP-solver that can be accessed via Pyomo. 

The current version uses Cbc to access Clp as an LP-solver, which is available in [this github repository](https://github.com/coin-or/Cbc). As NLP-solver we currently use IPOPT, which is available in [this github repository](https://github.com/coin-or/Ipopt) (we recomment using coinbrew for the installation).

As an alternative, [Bonmin](https://projects.coin-or.org/Bonmin/wiki/GettingStarted), which can also be obtained by using the [coinbrew script](https://coin-or.github.io/coinbrew/), contains both solvers so that it suffices to install Bonmin (instead of Cbc and IPOPT). 

Note that when using a different LP-solver (e.g. Cplex, Gurobi), you probably need to change the lines where the solver time is queried (as the pyomo interface is different for different LP-solvers).

Further dependencies are ```numpy``` and ```pandas``` which can both be obtained via pip.

You can run the test 
```
python test_modules.py
```
to see if you installed all required packages and if the method is running correctly.

## Code Examples
The following example runs the IPCP on a testinstance from 

Jan Kronqvist, Andreas Lundell, and Tapio Westerlund. 
"The extended supporting hyperplane algorithm for convex mixed-integer nonlinear programming." 
Journal of Global Optimization 64.2 (2016): 249-272.

```
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
```

which returns:

```
IPCP successful.
Objective value is  -20.90361501458559
The last iterate after the postprocessing step is:  [ 8.90361501 12.        ]
It took 15 LPs to compute this point.
```

This example can also be run from the folder ```/testinstances``` using the command ```python ex_ex_sup_hyp_1.py```. In this folder, three other examples are also available.

Changing algorithm settings can be done in the params.py file.

## Reproducing the experiments

The computational experiments from  [Generating feasible points for mixed-integer convex optimization problems by inner parallel cuts](http://www.optimization-online.org/DB_HTML/2018/11/6947.html) may be reproduced by running the scripts 

```
IPCP_on_granular_instances
IPCP_on_nongranular_instances
```
These scripts load the models from ```/minlplib_instances```, run the IPCP on these models and print and summarize the results. 

Note that mainly due to interfacing times between Pyomo and the LP/NLP solver, applying the IPCP to problems may take significantly longer than the reported run time. The latter corresponds to the time spent in the LP/NLP solver and thus excludes interfacing times and is hence closer to the time the method would actually take when integrated into a solver.

In the folder ```/solver_log``` we further report the log files of B-Hyb, B-OA and SCIP that we used as a comparison. This folder also contains a xlsx-file with aggregated results that we used to compare solving times and lower bounds of all solvers and models with and without reversed inner parallel cuts.

## Contact
Christoph.Neumann@kit.edu
