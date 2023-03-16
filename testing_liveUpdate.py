import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from waveguideSolver_funcs import *

#lattice constant
a = 2.737678125000000e-07
#hole diameter prefactor 1
d1 = 0.666
#hole diameter prefactor 2
d2 = 0.773

# Define a function to update the plot
def update_plot(xk,optimizationResult):
    result = getattr(optimizationResult,'x')
    print(result) 
    # x_vals.append(xk)
    # fitness_vals.append(objective(xk))
    # plt.clf()
    # plt.plot(fitness_vals)
    # plt.title('Fitness over Time')
    # plt.xlabel('Iteration')
    # plt.ylabel('Fitness')
    # plt.draw()
    # plt.pause(1)

# Initialize the starting point
sim_run = [0]
p0 = [a,d1,d2]
# Initialize the boundary values 
bnd = ((2.0e-07,3.0e-07),(0.1,0.8),(0.67,1.69))

# Run the optimization and update the plot at each iteration
fitness_vals = [unitCellOptimization_SiC_elliptical(p0)]
res = minimize(unitCellOptimization_SiC_elliptical, p0, method='Nelder-Mead', callback=update_plot)