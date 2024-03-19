# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:21:22 2024

@author: cbozonnet
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import pickle
from functions import data2save, parameters

# List filenames and legends
# list_files = ['ortega.pkl','no_osmoreg.pkl','wall_thinning.pkl','with_wall_building.pkl','complete.pkl','complete-Cs_ext200.pkl','complete-Cs_ext200-psi_ext-0.5.pkl','complete-osmoreg.pkl','complete-dpidt_eq_zero.pkl']
# my_legend = ['With osmoregulation','No osmoregulation','with wall thinning', 'with wall building','complete case','complete case - Cs_ext=200 mol/m3', 'same - low Psi','complete with dndt=0','complete with dPIdt=0']
# list_files = ['ortega.pkl','no_osmoreg.pkl','wall_thinning.pkl','with_wall_building.pkl','complete.pkl']
# my_legend = ['With osmoregulation','No osmoregulation','with wall thinning', 'with wall building','complete case']
list_files = ['complete-Cs_ext80-low_sigY.pkl','complete-Cs_ext200-high_sigY.pkl','complete-Cs_ext80-low_sigY-high_nu.pkl']
my_legend = ['Cs_ext=80-low_sigY','Cs=200-high_sigY','Cs=80-low_sigY-high_nu']
list_files = ['complete-Cs_ext200.pkl','complete-Cs_ext200-omegax2.pkl','complete-Cs_ext200-omegax0-osmoreg2.pkl']
my_legend = ['Cs_ext=200','Cs_ext=200 omegax2','Cs_ext=200 omegax0']
my_color =['b','orange','green']
# list_files = ['complete-Cs_ext200.pkl']
# my_legend = ['Cs_ext=200']
# list_files = ['complete-Cs_ext200-long.pkl','complete-Cs_ext200-n_constant.pkl','complete-Cs_ext200-PI_constant.pkl']
# my_legend = ['Cs_ext=200','n=cste','PI=cste']
# list_files = ['complete-Cs_ext200-dt0_5h.pkl','complete-Cs_ext200-dt2h.pkl','complete-Cs_ext200-dt4h.pkl','complete-Cs_ext200-dt10h.pkl','complete-Cs_ext200-dt50h.pkl']
# my_legend = ['dt=0.5h','dt=2h','dt=4h','dt=10h','dt=50h']
# list_files = ['complete-Cs_ext200-dt4h.pkl']
# my_legend = ['dt=4h']
size = len(list_files) # get the number of files

for i in range(size): # loop to open the files one by one and plot things
    with open(list_files[i], "rb") as file:
        data = pickle.load(file)
    sol = data.sol
    p = data.p
    t = data.t
    
    ####### Post process #########
    # extract data from solution
    ns=sol[:,0]
    sig_L=sol[:,1]
    L=sol[:,2]
    W=sol[:,3]
    
    # data treatment
    P = 2*W*sig_L/p.R0 # pressure from equilibrium condition
    Vh = np.pi*(p.R0-W)**2*L # water volume in the cell (m3)
    Cs = ns/Vh # Sugar concentration (mol/m3)
    PI = Cs*p.Rg*p.T # Cell osmotic potential (Pa)
    th=t/3600 # time in hours
    
    # Compute Lockhart's solution
    # NB: to match with it 
    # you should take W -> 0 in the parameters
    # In addition to osmoregulation & no wall synthesis
    Py0 = p.sig_Y*2*p.W0/p.R0 # yield pressure
    # corrected wall synthesis rate
    phi_w_star = p.phi_w*p.R0/(2*p.W0)*(1+2*p.W0/p.R0)
    # Compute growth rate
    Gth_cor = phi_w_star*p.phi_h/(phi_w_star+p.phi_h)*(p.Psi_src+p.Pi0-Py0)
    # exponential growth
    L_th = p.L0*np.exp(Gth_cor*th*3600)
    
    ######### Plots ##########
    plt.figure(1)
    plt.plot(th,L*1000000,label=my_legend[i],color=my_color[i])
    plt.plot(th,L_th*1000000,'--k')
    plt.plot([150, 150], [0, 84],':b',linewidth=1.5)
    # plt.plot([90, 90], [0, 42],color='orange',linewidth=1.5,linestyle=':')
    plt.xlim((0,300))
    plt.ylim((0,150))
    #plt.legend(loc='best')
    #plt.yscale('log')
    plt.xlabel(r"\textbf{t [h]}", fontsize=16)
    plt.ylabel(r"\textbf{L [$\mu$m]}", fontsize=16)
    plt.title(r"$\textbf{Cell length}$", fontsize=16)
    # Set grid and minor ticks
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    # Use LaTeX for tick labels (optional)
    plt.tick_params(labelsize=12, which='both', top=True, bottom=True, left=True, right=True)
    plt.tight_layout()
    #plt.arrow(x=40, y=87.5, dx=19, dy=-3.5, width=0.1, color='black',head_width=5,head_length=12)
    # plt.text(5,135, r"\bf{No\\  synthesis}", color="green", fontsize=20,rotation=0)
    plt.text(100,120, r"\bf{Lockhart}", color="k", fontsize=20,rotation=0)
    plt.text(200,85, r"\bf{Base case}", color="b", fontsize=20,rotation=0)
    # plt.text(170,25, r"\bf{Wall synthesis \\ 2x faster}", color="orange", fontsize=20,rotation=0, va="center")
    # Save the current figure as a PNG in the specified directory
    # Set a higher resolution for the JPEG (adjust as needed)
    dpi = 300  # Dots per inch
    # Save the figure as a high-resolution JPEG
    plt.savefig("growth_base_cases.jpeg", dpi=dpi)
    # Optional: Clear the figure to avoid saving multiple plots on top of each other
    #plt.clf()  # Clears the figure
    
    # # Specify the desired directory (replace with your actual path)
    # target_dir = "/path/to/your/directory"  # Adjust the path accordingly

    # # Create the directory if it doesn't exist (recommended)
    # os.makedirs(target_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # # Combine the directory and filename for the full path
    # filename = "length_base_case.png"  # Replace with your desired filename
    # full_path = os.path.join(target_dir, filename)



    
    plt.figure(2)
    plt.plot(th,sig_L/1e6,label=my_legend[i],color=my_color[i])
    plt.plot([150, 150], [0, 1],':b',linewidth=1.5)
    # plt.plot([90, 90], [0, 1],color='orange',linewidth=1.5,linestyle=':')
    plt.plot([0, 300], [1, 1],':k',linewidth=1.5)
    plt.xlim((0,300))
    # plt.ylim((0,1.1))
    plt.xlabel(r"\textbf{t [h]}", fontsize=16)
    plt.ylabel(r"\textbf{$\sigma$ [MPa]}", fontsize=16)
    plt.title(r"$\textbf{Wall stress}$", fontsize=16)
    # Set grid and minor ticks
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    # Use LaTeX for tick labels (optional)
    plt.tick_params(labelsize=12, which='both', top=True, bottom=True, left=True, right=True)
    plt.tight_layout()
    plt.text(250,0.92, r"\bf{$\sigma_Y$}", color="k", fontsize=20,rotation=0)
    # plt.savefig("sigma_base_cases.jpeg", dpi=dpi)
    
    plt.figure(3)
    plt.plot(th,P/1e6,color=my_color[i])
    plt.xlabel('t [h]')
    plt.ylabel('P [MPa]')
    
    plt.figure(4)
    Wp=W*1000000
    plt.plot(th,Wp,label=my_legend[i],color=my_color[i])
    plt.plot([150, 150], [0, 3.15],':b',linewidth=1.5)
    plt.plot([0, 150], [3.15, 3.15],':b',linewidth=1.5)
    plt.plot([90, 90], [0, 4],color='orange',linewidth=1.5,linestyle=':')
    plt.plot([0, 90], [4, 4],color='orange',linewidth=1.5,linestyle=':')
    plt.xlim((0,300))
    plt.ylim((0,8))
    plt.xlabel(r"\textbf{t [h]}", fontsize=16)
    plt.ylabel(r"\textbf{W [$\mu$m]}", fontsize=16)
    plt.title(r"$\textbf{Wall thickness}$", fontsize=16)
    # Set grid and minor ticks
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    # Use LaTeX for tick labels (optional)
    plt.tick_params(labelsize=12, which='both', top=True, bottom=True, left=True, right=True)
    plt.tight_layout()
    #plt.text(250,0.92, r"\bf{$\sigma_Y$}", color="k", fontsize=20,rotation=0)
    plt.savefig("W_all_cases.jpeg", dpi=dpi)

    # # plt.figure(5)
    # # plt.plot(th,ns/p.ns0,label=my_legend[i])
    # # plt.legend(loc='best')
    # # plt.xlabel('t [h]')
    # # plt.ylabel('ns/ns0 [-]')
    
    plt.figure(6)
    plt.plot(th,PI/1e6,color=my_color[i])
    # plt.legend(loc='best')
    plt.xlabel('t [h]')
    plt.ylabel('PI [MPa]')
    
    # plt.figure(7)
    # plt.plot(th,Cs,label=my_legend[i])
    # plt.legend(loc='best')
    # plt.xlabel('t [h]')
    # plt.ylabel('Cs [mol/m3]')
    
    # plt.figure(8)
    # plt.plot(th[1:],np.diff(W),label=my_legend[i])
    # plt.legend(loc='best')
    # plt.xlabel('t [h]')
    # plt.ylabel('DW [m]')
    