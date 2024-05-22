# Test file for the hammer force described by 
# M d2h/dt2 = -Fh(t)
# Based on the FDM review by Salo 


import numpy as np 
import matplotlib.pyplot as plt 

dt = 0.01 
max_t = 100

eta = np.zeros(int(max_t/dt))
hammer_force = np.zeros(int(max_t/dt))

def hammer_force(K_h, eta):
    pass