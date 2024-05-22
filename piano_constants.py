import numpy as np 

N = 50 # no. of string segments 
dx = 0.1 
dt = 0.01
max_t = 2

# Values follow for C4 as listed in Chainge and Askenfelt 1993 
fund_freq = 262 # Hz, assume standard tuning frequency A4 
L = 0.62
M_s = 3.93 * 3 #string mass 
T = 670 # string tension 
E = 2.0e11 # Youngs modulus 
eps = 3.82e-5 # string stiffness parameter
b_1 = 0.5 # damping coefficient 
b_3 = 6.25e-9 # damping coefficient 

mu = E # linear mass density of string ## Not sure on this one 
c = np.sqrt(T/mu)
v_h0 = 10 # m / s**2
p = 2.5 # ideally between 2 and 3 
K_h = 4.5e9 #generalised hammer stiffness 
k_0 = int(0.12*N) #Hammer location index 

M_h = 2.97 * 3 #hammer mass 



samp_freq = 32000 # Hz



# Placeholder values listed below 


eps = 1.0 # string stiffness parameter
b_1 = 1.0 # damping coefficient 
b_3 = 1.0 # damping coefficient 
T = 1.0 # string tension 
mu = 1.0 # linear mass density of string 
c = np.sqrt(T/mu)
v_h0 = 10 # m / s**2
p = 2 # ideally between 2 and 3 
K_h = 1.0 #generalised hammer stiffness 
k_0 = 10 #Hammer location index 
M_s = 1.0 #string mass 
M_h = 1000 #hammer length 

N = 50 # no. of string segments 
dx = 0.1 
dt = 0.01
max_t = 2
fund_freq = 440 # Hz, assume standard tuning frequency A4 
samp_freq = 32000 # Hz