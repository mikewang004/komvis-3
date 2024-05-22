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
v_h0 = 1.0 # m / s**2
p = 2.5 # ideally between 2 and 3 
K_h = 1.0 #generalised hammer stiffness 
k_0 = 5 #Hammer location index 
M_s = 1.0 #string mass 
M_h = 1.0 #hammer length 

N = 100 # no. of string segments 
dx = 0.01 
dt = 0.01
deta = 0.01 
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

def sum_both_neighbours(A):
    """In 1D array, returns sum of both left and right neighbour."""
    return np.roll(A, -1) + np.roll(A, + 1)

def hammer_force(eta, y_x0, K_h=K_h, p=p):
    """Note eta should be put in as a scalar"""
    return K_h * np.abs(eta - y_x0)**p

# def runge_kutta_hammer_force(h,m_H, eta, A_t, n, i):
#     """dt follows implicitly from time index. Force is function, F_h current value for Force.""" 
#     h = dt
#     eta_n = eta[n]
#     deta = eta[n] - eta[n-1]
#     k11 =  h * deta
#     k12 =  hammer_force(eta_n, A_t[k_0])/m_H

#     #k21 = h * deta + 0.5 * k12
#     #k22 = h * hammer_force()


def solve_eta_hammer_force(m_H, eta, A_t, n):
    n_1 = 2 * eta[n] - eta[n-1] + dt**2 / m_H * hammer_force(eta[n], A_t[k_0])
    print(n_1)
    return n_1


def construct_solution_matrix(N, max_t):
    """Assuming discretisation in N string segments and max_t/dt time segments"""
    # Implement FDM here 
    end_begin_conditions = 3 # index
    A = np.zeros([N, int(max_t/dt)]); eta = np.zeros(int(max_t/dt)); F_h = np.zeros(int(max_t/dt))
    A[0, 0] = 0; A[N-1, max_t-1] = 0 # Define boundary conditions 
    A[:end_begin_conditions, :end_begin_conditions] = 0 # Implement easiest solution that first 3 timesteps are silent 
    eta[1] = v_h0 * dt; F_h[1] = K_h * np.abs(eta[1] - A[k_0, 1])**p
    # Follwing three lines based on Salo's review on FDM, page 8 
    A[:, 2] = sum_both_neighbours(A[:, 1]) - A[:, 0] + dt**2 * K_h * F_h[1] / M_s # eq (39)
    eta[2] = 2 * eta[1] - eta[0] - dt**2 * F_h[1] / M_h # eq (40)
    F_h[2] = K_h * np.abs(eta[2] - A[k_0, 2])**p # eq (41)
    for n in range(end_begin_conditions, int(max_t/dt)-1): # time index
        eta[n] = solve_eta_hammer_force(M_h, eta, A[:, n-1], n-1)
        F_h[n] = hammer_force(eta[n], A[k_0, n])
        for i in range(1, N-2): # place index 
            A[i, n + 1] = a_1 * A[i, n] + a_2 * A[i, n-1] + \
            a_3 * (A[i+1, n] + A[i -1, n]) + a_4 * (A[i + 2, n] + A[i - 2, n]) + \
            a_5 * (A[i+1, n-1] + A[i -1, n-1,] + A[i, n - 2]) + \
            (dt**2 * N * K_h) / (D * M_s) * F_h[n] 





    return A

A = construct_solution_matrix(N, max_t)
#print(A)





