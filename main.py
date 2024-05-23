import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# Equation implemented from Chaigne and Askenfeld 1993 
# Define constants 

N = 50 # no. of string segments 
dx = 0.01 
dt = 0.0001
max_t = 1.0

# Values follow for C4 as listed in Chainge and Askenfelt 1993 
fund_freq = 262 # Hz, assume standard tuning frequency A4 
L = 0.62 # m
M_s = 3.93 * 3 #string mass 
T = 670 # string tension 
E = 2.0e11 # Youngs modulus 
eps = 3.82e-5 # string stiffness parameter
b_1 = 0.5 # damping coefficient 
b_3 = 6.25e-9 # damping coefficient 

M_h = 2.97 * 3 #hammer mass 
mu = M_h/L # Transversal string desntiy, Placeholder value 
#mu = 10e6
c = np.sqrt(T/mu)
v_h0 = 0.5 # m / s**2; note 4.0 1.5 0.5 for resp. forte mezzo forte piano 
p = 2.5 # ideally between 2 and 3 
K_h = 4.5e9 #generalised hammer stiffness 
k_0 = int(0.12*N) #Hammer location index 




samp_freq = 32000 # Hz

D = 1 + b_1 * dt + 2 * b_3 / dt 
r = c * dt / dx

a_1 = (2 - 2*r**2 + b_3/dt - 6 * eps * N**2 * r**2)/D
a_2 = (-1 + b_1 * dt + 2*b_3/dt)/D
a_3 = (r**2 * (1 + 4 * eps * N**2))/D
a_4 = (b_3/dt - eps * N**2 * r**2)/D
a_5 = (-b_3 / dt) / D
print(a_4)

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
            if k_0 - 3 < i < k_0 + 3:
                hammer_window = 1
            else:
                hammer_window = 0
            A[i, n + 1] = a_1 * A[i, n] + a_2 * A[i, n-1] + \
            a_3 * (A[i+1, n] + A[i -1, n]) + a_4 * (A[i + 2, n] + A[i - 2, n]) + \
            a_5 * (A[i+1, n-1] + A[i -1, n-1,] + A[i, n - 2]) + \
            ((dt**2 * N * K_h) / (M_s)) * F_h[n] * hammer_window
            #((dt**2 * N * K_h) / (D * M_s)) * F_h[n] * hammer_window
        A[N-2, n+1] = -1 * A[1, n+1]
        #print(A[1, n], A[N-2, n])
        if hammer_force_not_needed != True:
            eta[n+1] = solve_eta_hammer_force(M_h, eta, A[:, n], n)
            F_h[n+1] = hammer_force(eta[n], A[k_0, n])
        if eta[n+1] < A[k_0, n+1]:
            F_h[n+1] = 0 
            hammer_force_not_needed = True
        #print(n)
        #print(eta[n+1])
        #print(A[:, n+1])



    return A

A = construct_solution_matrix(N, max_t)


#plt.scatter(np.arange(N), A[:, -1])
#plt.plot(A[:, -1])
#plt.show()
print(A)


# Plot wave 

# Create the figure and axis

fig= plt.figure()
ax = plt.axes(xlim=(0, N), ylim = (-10e12, 10e12))
line, = ax.plot([], [], lw=3)
x = np.arange(0, N)

def plot_A(frame):
    print(frame)
    frame = frame * 50 
    line.set_data(x, A[:, int(frame)])
    return A[:, int(frame)]

ani = animation.FuncAnimation(fig, plot_A, frames = 100)
plt.show()
ani.save("myvideo.mp4")



