#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: evaloughridge
"""

import numpy as np
import matplotlib.pyplot as plt

import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# 1500nm data with L = 181, f = 0.78 au
data_IR15 = np.loadtxt("Harm_vel_0001_z_IR1.5_only")
data_xuv_only = np.loadtxt("Harm_vel_0001_z_xuv_0.1_only") 

data_xuv_IR15_d0 = np.loadtxt("Harm_vel_0001_z_xuv_IR1.5") #f = 0.78 au, d = 0 fs, I = 0.1%
data_xuv_IR15_dm2 = np.loadtxt("Harm_vel_0001_z_dm2") #f = 0.78 au, d = -2 fs
data_xuv_IR15_a181_dm13 = np.loadtxt("Harm_vel_0001_z_dm13") #f = 0.78 au, d = -13 fs

# Detune, d = 0 fs
data_xuv_detune_16 = np.loadtxt("Harm_vel_0001_z_16.06eV") 
data_xuv_detune_40 = np.loadtxt("Harm_vel_0001_z_40eV")
data_xuv_detune_5 = np.loadtxt("Harm_vel_0001_z_4.89eV")

# Intensities
data_xuv_IR15_x1 = np.loadtxt("Harm_vel_0001_z_xuv_IRx1.0") #f = 0.78 au, I = 1%
data_xuv_IR15_x01 = np.loadtxt("Harm_vel_0001_z_xuv_IRx0.01") #f = 0.78 au, I = 0.01%
data_xuv_IR15_x001 = np.loadtxt("Harm_vel_0001_z_xuv_IRx0.001") #f = 0.78 au, I = 0.001%
data_xuv_IR15_x0001 = np.loadtxt("Harm_vel_0001_z_xuv_IRx0.0001") #f = 0.78 au, I = 0.0001%

#-------------------------------------------------------------------------------------------

Int = 1e14 #Intensity
c = 3e8 #speed of light 
p_e = 27.212 #eV au conversion
I_p = 24.59 #Ionisation potential 
l_l = 1500e-9 #1500nm laser
U_p = 9.33e16*Int*((l_l/1e9)**2) #pondermotive energy
cutoff_2 = I_p + 3.17*U_p #empirical cutoff formula
f_au = ((c/l_l)/(4.13567e16))*2*np.pi #frequency in atomic units 

datasets = [
    (data_IR15,  'cornflowerblue',  'IR only, 1.5 $\mu$m', 0, 'dashed'),
    (data_xuv_only,  'orange',  'XUV only, 21.23 eV', 0, (0, (3, 3, 1, 3, 1, 3))),
    (data_xuv_IR15_d0,  'navy',  'XUV, d = 0 fs', 0, 'dashdot'),
  
    (data_xuv_IR15_dm2, 'aqua', 'XUV', -2.01, '-'),
    (data_xuv_IR15_a181_dm13, 'magenta', 'XUV', -13.16, (0, (1, 1))),
    
    #(data_xuv_detune_5, 'lawngreen', 'Detuned 5 eV XUV', 0, '-'),
    #(data_xuv_detune_16, 'red', 'Detuned 16 eV XUV', 0, '-'),
    #(data_xuv_detune_40, 'orange', 'Detuned 40 eV XUV', 0, '-'),
    
    #(data_xuv_IR15_x1,  'red',  'XUV, Int. = 1%', 0, '-'),
    #(data_xuv_IR15_d0,  'navy',  'XUV, Int. = 0.1%', 0, '-'),
    #(data_xuv_IR15_x01,  'lawngreen',  'XUV, Int. = 0.01%', 0, '-'),
    #(data_xuv_IR15_x001,  'orange',  'XUV, Int. = 0.001%', 0, '-'),
    #(data_xuv_IR15_x0001,  'violet',  'XUV, Int. = 0.0001%', 0, '-'),
]
plt.figure(figsize=(8, 5))

for data, color, intensity_label, delay, linestyle in datasets:
    energy = data[:, 0]
    spectrum = data[:, 1]

    window_size = 7
    window = np.ones(window_size) / window_size
    spectrum_smooth = np.convolve(spectrum, window, mode='same')

    if delay:
        label = f"{intensity_label}, d = {delay} fs"
    else:
        label = intensity_label

    #plt.semilogy(energy, spectrum, ':', color=color,linestyle=linestyle, label=f"{label} Raw HHG")
    plt.semilogy(energy, spectrum_smooth, color=color,linestyle=linestyle, linewidth=1.5, label=f"{label}")


plt.xlabel("Photon Energy (eV)")
plt.ylabel("HHG Intensity (arb. units)")
#plt.title("Helium HHG Spectrum (IR=1500nm)")
plt.xlim(0, 100)
plt.ylim(1e-12,1e4)

plt.vlines(92, 0, 10e3, 'green', '--', label='Calculated Cutoff')
plt.axvline(79, color='grey', linestyle=':', label='Reduced Cutoff')
plt.legend(loc ='upper right')
plt.savefig("helium_hhg_IR15_xuv_delay_dash.png")
plt.show()


# ============================================================
# For Efield Plot
# ============================================================
"""
timefile = "EField.he_IR1.5_xuv_idm2"    
vel_col = 1                            

Efield_data = np.loadtxt(timefile, skiprows=1)
time = Efield_data[:,0]
E_col = 3
E_t = Efield_data[:,E_col]
dt = time[1] - time[0]
N = len(time)

xuv_25 = 0.195*2.5
IR_25 = 5.067*2.5
delay = -2

plt.figure(figsize=(8,4))
plt.plot((time/41.34) - IR_25, E_t,  color = 'purple')
plt.xlabel("Time (fs)")
plt.ylabel("E-field (au)")
x1 = delay - (xuv_25)
x2 = x1 + (xuv_25 * 2)
plt.axvspan(x1, x2, color='firebrick', alpha=0.2)
plt.title(f"Combined Laser Electric Field: 2.5 cycle ramp on/off pulse, with d = {delay} fs")
plt.grid()
#plt.savefig('Laser_EField_IR15_dm2.png')
plt.show()
"""



