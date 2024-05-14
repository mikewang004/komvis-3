import numpy as np 
import matplotlib.pyplot as plt 


# Equation implemented from Chaigne and Askenfeld 1993 
# Define constants 

eps = 1.0 # string stiffness parameter
b_1 = 1.0 # damping coefficient 
b_3 = 1.0 # damping coefficient 
T = 1.0 # string tension 
mu = 1.0 # linear mass density of string 
c = np.sqrt(T/mu)

N = 1.0e3 # no. of string segments 
dx = 0.001 
dt = 0.001

D = 1 + b_1 * dt + 2 * b_3 / dt 
r = c * dt / dx

a_1 = (2 - 2*r**2 + b_3/dt - 6 * eps * N**2 * r**2)/D
a_2 = (-1 + b_1 * dt + 2*b3/dt)/D
a_3 = (r**2 * (1 + 4 * eps * N**2))/D
a_4 = (b_3/dt - eps * N**2 * r**2)/D
a_5 = (-b_3 / dt) / D



#Implement FDM here 


def 
