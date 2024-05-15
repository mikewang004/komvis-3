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

N = 100 # no. of string segments 
dx = 0.01 
dt = 0.01
max_t = 1
fund_freq = 440 # Hz, assume standard tuning frequency A4 
samp_freq = 32000 # Hz

D = 1 + b_1 * dt + 2 * b_3 / dt 
r = c * dt / dx

a_1 = (2 - 2*r**2 + b_3/dt - 6 * eps * N**2 * r**2)/D
a_2 = (-1 + b_1 * dt + 2*b_3/dt)/D
a_3 = (r**2 * (1 + 4 * eps * N**2))/D
a_4 = (b_3/dt - eps * N**2 * r**2)/D
a_5 = (-b_3 / dt) / D


def construct_solution_matrix(N, max_t):
    # Implement FDM here 
    end_begin_conditions = 3 # index
    A = np.zeros([N, int(max_t/dt)])
    A[0, 0] = 0; A[N-1, max_t-1] = 0 # Define boundary conditions 
    A[:end_begin_conditions, :end_begin_conditions] = 0 # Implement easiest solution that first 3 timesteps are silent 
    for n in range(end_begin_conditions, int(max_t/dt)-1): # time index
        for i in range(end_begin_conditions-1, N-2): # place index 
            A[i, n + 1] = a_1 * A[i, n] + a_2 * A[i, n-1] + \
            a_3 * (A[i+1, n] + A[i -1, n]) + a_4 * (A[i + 2, n] + A[i - 2, n]) + \
            a_5 * (A[i+1, n-1] + A[i -1, n-1,] + A[i, n - 2])
    return A

A = construct_solution_matrix(N, max_t)
print(A)





