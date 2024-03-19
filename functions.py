# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:28:55 2024

@author: cbozonnet
"""

import numpy as np

class parameters():
    def __init__(self):
        # Mechanical parameters
        self.E = 30e6 # Young Modulus [Dumais, new Phytol., 2021]
        self.nu = 0.8e10 # wall viscosity in Pa.s [Dumais, new Phytol., 2021]
        self.rho_w = 1500 # cell wall density in kg/m3 [Friend et al, Nature com., 2022]
        self.sig_Y = 1e6 # Plastic threshold in Pa [user choice]
        
        # Hydraulic parameters
        self.kh = 1e-16 # hydraulic conduct in m/Pa/s [Dumais, new Phytol., 2021]
        self.Psi_src = 0 # (=P-PI)_{ext} : external (xylem) potential in Pa
        self.P_ext = 0 # external pressure      
        self.Pi0 = 0.5e6 # initial osmotic potential in Pa [Uggla et al, Plant phy, 2001]
        
        # Wall synthesis parameters [Friend et al, Nature com., 2022]
        self.omega = 0*2.2e-4 # normalised rate of mass growth (kg/m3/s) at T0
        self.Eaw = 1.43 # activation energy for wall building (eV)
        self.kb = 8.617e-5 # Bolztmann's constant (eV/K)
        self.Km = 14.9 # Michaelis constant for wall synthesis (mol/m3)
        self.MMs = 0.342 # molar mass for sucrose (kg/mol)
               
        # Temperature
        self.T0 = 273.15 # Reference temperature in K
        self.T = self.T0 + 15 # Actual temperature in K
        
        # Geometry
        self.R0 = 10e-6 # initial radius in m
        self.L0 = 10e-6 # initial cell length in m 
        self.W0 = 0.5e-6 # initial wall thickness in m
        
        # Simulation parameters
        self.t0 = 0
        self.t_end = 1000*3600 # final time (s)
        self.dt = 0.5*3600 # time step (s)
      
        # Physical constant
        self.Rg = 8.314 # Ideal gas constant
        
        # Sugar transport
        self.eta_s = 2.6e-16 # sugar diffusion constant (m3/s) [Friend et al, Nature com., 2022]
        self.Cs_ext = 200 # sugar source conentration (mol/m3) [Uggla et al, Plant phy, 2001]

        # derived quantities
        self.phi_w = 1/self.nu # wall extensibility (1/(Pa.s))
        self.phi_h = 2*self.kh/self.R0 # =A*k_h/V (2k_w/R for a cylinder) 
        self.alpha = self.phi_h/(self.phi_h+self.phi_w)
        self.Vh0 = np.pi*(self.R0-self.W0)**2*self.L0 # initial water volume
        self.ns0 = self.Vh0*self.Pi0/(self.Rg*self.T) # intial sugar content
        self.P0 = 0 # self.sig_Y*2*self.W0/self.R0 # initial turgor pressure in Pa
        self.sig_l0 = self.P0*self.R0/(2*self.W0) # intial wall stress
        
class data2save: # this creates a structure to save all datas
    def __init__(self, p, t, sol):
        self.p = p
        self.t = t
        self.sol = sol