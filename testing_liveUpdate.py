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


def callback_func(xk,**kwargs):
    print(xk)
    print(kwargs)
    # print(f"Iteration {kwargs['jac']}: x={xk}, f(x)={kwargs['fun']}")

# Define objective function
def objective(x):
    return (x[0] - 1)**2 + (x[1] - 2.5)**2

# Define initial guess
x0 = [0, 0]

# Call optimization function with callback
res = minimize(objective, x0, method='BFGS', callback=callback_func)


########## code to run on the cluster 
# # Define a function to update the plot
# def update_plot(xk):
#     print("Debugging the live update plot ===================")
#     print(f"Iteration {len(res.x_iters)}: x={xk}, f(x)={res.fun}")
#     # x_vals.append(xk)
#     # fitness_vals.append(objective(xk))
#     # plt.clf()
#     # plt.plot(fitness_vals)
#     # plt.title('Fitness over Time')
#     # plt.xlabel('Iteration')
#     # plt.ylabel('Fitness')
#     # plt.draw()
#     # plt.pause(1)
    
    
# # Initialize the starting point
# sim_run = [0]
# p0 = [a,d1,d2]
# # Initialize the boundary values 
# bnd = ((2.0e-07,3.0e-07),(0.1,0.8),(0.67,1.69))

# # Run the optimization and update the plot at each iteration
# fitness_vals = [unitCellOptimization_SiC_elliptical(p0)]
# res = minimize(unitCellOptimization_SiC_elliptical, p0, method='Nelder-Mead', callback=update_plot)