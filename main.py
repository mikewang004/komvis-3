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
v_h0 = 10 # m / s**2
p = 2 # ideally between 2 and 3 
K_h = 1.0 #generalised hammer stiffness 
k_0 = 2 #Hammer location index 
M_s = 1.0 #string mass 
M_h = 1000 #hammer length 

N = 10 # no. of string segments 
dx = 0.1 
dt = 0.01
max_t = 1.0
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
    #print("k_0 string extension as follows: " + str(y_x0))
    return K_h * np.abs(eta - y_x0)**p



def solve_eta_hammer_force(m_H, eta, A_t, n):
    #print("hammer force as follows:" + str(hammer_force(eta[n], A_t[k_0])))
    n_1 = 2 * eta[n] - eta[n-1] + (dt**2 / m_H) * hammer_force(eta[n], A_t[k_0])
    #print(n_1)
    return n_1


def construct_solution_matrix(N, max_t):
    """Assuming discretisation in N string segments and max_t/dt time segments"""
    # Implement FDM here 
    hammer_force_not_needed = False
    end_begin_conditions = 2 # index
    A = np.zeros([N, int(max_t/dt)]); eta = np.zeros(int(max_t/dt)); F_h = np.zeros(int(max_t/dt))
    A[0, 0] = 0; A[N-1, int(max_t/dt)-1] = 0 # Define boundary conditions 
    eta[1] = v_h0 * dt; F_h[1] = K_h * np.abs(eta[1] - A[k_0, 1])**p
    # Follwing three lines based on Salo's review on FDM, page 8 
    A[:, 2] = sum_both_neighbours(A[:, 1]) - A[:, 0] + dt**2 * K_h * F_h[1] / M_s # eq (39)
    eta[2] = 2 * eta[1] - eta[0] - dt**2 * F_h[1] / M_h # eq (40)
    F_h[2] = K_h * np.abs(eta[2] - A[k_0, 2])**p # eq (41)
    for n in range(end_begin_conditions, int(max_t/dt)-1): # time index
        for i in range(1, N-2): # place index 
            if 0 < i < 3:
                hammer_window = 1
            else:
                hammer_window = 0
            A[i, n + 1] = a_1 * A[i, n] + a_2 * A[i, n-1] + \
            a_3 * (A[i+1, n] + A[i -1, n]) + a_4 * (A[i + 2, n] + A[i - 2, n]) + \
            a_5 * (A[i+1, n-1] + A[i -1, n-1,] + A[i, n - 2]) + \
            ((dt**2 * N * K_h) / (M_s)) * F_h[n] * hammer_window
            #((dt**2 * N * K_h) / (D * M_s)) * F_h[n] * hammer_window

            
        if hammer_force_not_needed != True:
            eta[n+1] = solve_eta_hammer_force(M_h, eta, A[:, n], n)
            F_h[n+1] = hammer_force(eta[n], A[k_0, n])
        if eta[n+1] < A[k_0, n+1]:
            F_h[n+1] = 0 
            hammer_force_not_needed = True
        print(n)
        print(eta[n+1])
        print(A[k_0, n+1])



    return A

A = construct_solution_matrix(N, max_t)
#print(A)





