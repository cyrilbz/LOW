# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:17:31 2024

@author: cbozonnet
"""
        
# This is a program for a Lockhart/Ortega model
# with elongation
# with wall thickening related to sugar availability
# with sugar uptake
# all parameters are in function.py

import numpy as np
from scipy.integrate import odeint
import pickle
from functions import data2save, parameters
import time
import os

# Start the timer.
start_time = time.time()

# Specify the output filename
filename ='complete-Cs_ext200-omegax0-osmoreg2.pkl'
overwrite = True # True if you want to overwrite the existing file 

if os.path.exists(filename) and overwrite==False: # Check if the file already exists
    raise Exception("This output file already exists.") 
   
print("Simulation running: " + str(filename))
        
def dydt(y,t,p): # all physical equations for time integration
    # extract data from y vector
    ns = y[0]# extract sugar content
    sig_l = y[1] # extract longitudinal stress
    L = y[2] # extract cell length
    W = y[3] # extract wall thickness
    
    # Compute useful quantities
    Vh = np.pi*(p.R0-W)**2*L # water volume in the cell (m3)
    Cs = ns/Vh # Sugar concentration (mol/m3)
    PI =  Cs*p.Rg*p.T # Cell osmotic potential (Pa)
    dMmax = p.omega*Vh*np.exp(p.Eaw/p.kb*(1/p.T0-1/p.T)) # Maximal speed of mass increment 
    dMdt = dMmax*Cs/(Cs+p.Km) # Mass growth (kg/s)
    myphi_w = p.phi_w # max(p.phi_w*(1-t/(p.t_end)),0) #custom extensibility
    fs = p.eta_s*(p.Cs_ext - Cs) # sugar flux
    #fs = max(fs,0) # take only influx
    
    # Compute time changes
    dsig_ldt = p.E*p.phi_h/(1+2*W/p.R0)*(p.Psi_src - 2*W*sig_l/p.R0 -p.P_ext + PI) \
                + p.E/(Vh*p.rho_w*(1+2*W/p.R0))*dMdt \
                - p.E*myphi_w*max(sig_l-p.sig_Y,0) # longitudinal wall stress
    dLdt = L*myphi_w*max(sig_l-p.sig_Y,0) + L/p.E*dsig_ldt # cell length    
    dWdt = 1/(p.rho_w*2*np.pi*p.R0*L)*dMdt - W/L*dLdt # cell thickness
    #dnsdt = - 1/p.MMs*dMdt + fs # Cell sugar content (mol/s)
    dnsdt = ns*(1/L*dLdt+2*np.pi*L/Vh*(W*dWdt-p.R0*dWdt)) # Imposes dPI/dt=0

    # return the vector of time changes
    return np.array([dnsdt, dsig_ldt, dLdt, dWdt]) 
    
######### Main program ##########
p = parameters() # get all parameters       
y0 = [p.ns0,p.sig_l0,p.L0,p.W0] # initial values
t = np.arange(p.t0,p.t_end+p.dt,p.dt) # time vector
sol = odeint(dydt, y0, t, args=(p,)) # resolution

######## Store solution #########
data = data2save(p,t,sol) # create data structure
with open(filename, "wb") as file: # open file
    pickle.dump(data, file) # save data

# Calculate the elapsed time.
elapsed_time = time.time() - start_time

# Display the elapsed time to the user.
print("The elapsed time is {} seconds.".format(elapsed_time))
