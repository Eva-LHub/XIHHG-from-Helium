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

# 1500nm data with L = 181
data_IR15 = np.loadtxt("Harm_vel_0001_z_IR1.5_only")
data_xuv_IR15_d0 = np.loadtxt("Harm_vel_0001_z_xuv_IR1.5") # d = 0 fs

data_xuv_IR15_d2 = np.loadtxt("Harm_vel_0001_z_d2") # d = 2 fs 
data_xuv_IR15_d4 = np.loadtxt("Harm_vel_0001_z_d4") # d = 4 fs
data_xuv_IR15_d6 = np.loadtxt("Harm_vel_0001_z_d6") # d = 6 fs

data_xuv_IR15_dm05 = np.loadtxt("Harm_vel_0001_z_dm0.5") # d = -0.5 fs
data_xuv_IR15_dm1 = np.loadtxt("Harm_vel_0001_z_dm1") # d = -1 fs
data_xuv_IR15_dm15 = np.loadtxt("Harm_vel_0001_z_dm1.5") # d = -1.5 fs
data_xuv_IR15_dm2 = np.loadtxt("Harm_vel_0001_z_dm2") # d = -2 fs
data_xuv_IR15_dm4 = np.loadtxt("Harm_vel_0001_z_dm4") # d = -4 fs
data_xuv_IR15_dm6 = np.loadtxt("Harm_vel_0001_z_dm6") # d = -6 fs
data_xuv_IR15_dm8 = np.loadtxt("Harm_vel_0001_z_dm8") # d = -8 fs
data_xuv_IR15_dm10 = np.loadtxt("Harm_vel_0001_z_dm10") # d = -10 fs
data_xuv_IR15_dm13 = np.loadtxt("Harm_vel_0001_z_dm13") # d = -13 fs
data_xuv_IR15_dm15 = np.loadtxt("Harm_vel_0001_z_dm15") # d = -15 fs
data_xuv_IR15_dm17 = np.loadtxt("Harm_vel_0001_z_dm17") # d = -17 fs

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
    (data_xuv_IR15_d0,   0),  
    (data_xuv_IR15_dm05, -0.5),
    
    #(data_xuv_IR15_dm1, -1),
    #(data_xuv_IR15_dm15, -1.5),
    
    (data_xuv_IR15_dm2, -2.01),
    (data_xuv_IR15_d2,  2),
    
    
    (data_xuv_IR15_dm4, -4),
    (data_xuv_IR15_d4,  4),
    
    (data_xuv_IR15_dm6,  -6),
    (data_xuv_IR15_d6,  6),
    
    (data_xuv_IR15_dm8,  -8),
    (data_xuv_IR15_dm10, -10),
    (data_xuv_IR15_dm13,  -13),
    (data_xuv_IR15_dm15,  -15),


]

#------------------------------------------------------------------------------------------

from scipy.interpolate import interp1d

data_with_delays = [(data, delay) for data, delay in datasets[1:]]
data_with_delays.sort(key=lambda x: x[1])

sorted_datasets = [d[0] for d in data_with_delays]
sorted_delays = np.array([d[1] for d in data_with_delays])

E_common = np.linspace(0, 100, 800)  
spectra_interp = [] 

for data in sorted_datasets:
    f_interp = interp1d(data[:, 0], data[:, 1], bounds_error=False, fill_value=0)
    spectra_interp.append(f_interp(E_common))


spectra_2D = np.array(spectra_interp)
spectra_2D = np.log10(np.clip(spectra_2D, 1e-12, None))  

plt.figure(figsize=(8,5))
extent = [E_common[0], E_common[-1], sorted_delays[0], sorted_delays[-1]]

plt.imshow(spectra_2D, aspect='auto', extent=extent,origin='lower', cmap='jet')
plt.colorbar(label='HHG Intensity (arb. units)')
plt.xlabel('Photon Energy (eV)')
plt.ylabel('XUV–IR Delay (fs)')
#plt.title('Helium HHG Spectrum vs XUV–IR Delay')

plt.axhline(0, color='white', linestyle='-')
#plt.axhline(-2.01, color='white', linestyle=':')

ax = plt.gca()

# 32 eV
x1 = 32 + I_p
ax.axvline(x1, linestyle=':', color='black')
ax.text(x1, 1.02, 'P1 + I$_p$',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 54.5 eV (slightly higher)
x4 = 54.5 + I_p
ax.axvline(x4, linestyle=':', color='black')
ax.text(x4-2, 1.02, 'P4 + I$_p$',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 58 eV (slightly lower)
x2 = 58 + I_p
ax.axvline(x2, linestyle=':', color='black')
ax.text(x2+1, 1.08, 'P2 + I$_p$',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 68.5 eV
x3 = 68.5 + I_p
ax.axvline(x3, linestyle=':', color='black')
ax.text(x3, 1.02, 'P3 + I$_p$',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')


plt.plot(x3, -2.01, marker='o', color = 'black') 
plt.plot(x4, 0, marker='o', color = 'black') 


plt.tight_layout()
plt.savefig("helium_hhg_smooth_xuv_delayheat_line.png")
plt.show()


