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
data_xuv_IR15_d0 = np.loadtxt("Harm_vel_0001_z_xuv_IR1.5") #f = 0.78 au, d = 0 fs, I = 0.1%
data_xuv_IR15_dm2 = np.loadtxt("Harm_vel_0001_z_dm2") #f = 0.78 au, d = -2 fs

# Pulse Train
data_IR_pulse = np.loadtxt("Harm_vel_0001_z_train")
data_xuv_pulse_d0 = np.loadtxt("Harm_vel_0001_z_train_d0")
data_xuv_pulse_dm2 = np.loadtxt("Harm_vel_0001_z_train_dm2")

#-------------------------------------------------------------------------------------------

Int = 1e14 #Intensity
c = 3e8 #speed of light 
p_e = 27.212 #eV au conversion
I_p = 24.59 #Ionisation potential 
l_l = 1500e-9 #1500nm laser
U_p = 9.33e16*Int*((l_l/1e9)**2) #pondermotive energy
cutoff_2 = I_p + 3.17*U_p #empirical cutoff formula
f_au = ((c/l_l)/(4.13567e16))*2*np.pi #frequency in atomic units 
harm = f_au*p_e #Harmonics


#-------------------------------------------------------------------------------------------

datasets = [
    
    (data_IR_pulse,  'cornflowerblue',  'Pulse train: IR only, 1.5 $\mu$m', 0),
    (data_xuv_pulse_d0,  'navy',  'Pulse train: XUV, d = 0 fs', 0),
    
    (data_xuv_pulse_dm2, 'aqua', 'Pulse train: XUV', -2.01),

    (data_IR15,  'lawngreen',  'IR only, 1.5 $\mu$m', 0),
    (data_xuv_IR15_d0,  'red',  'XUV, d = 0 fs', 0),
  
    (data_xuv_IR15_dm2, 'orange', 'XUV', -2.01),
]

plt.figure(figsize=(8, 5))

for data, color, intensity_label, delay in datasets:
    energy = data[:, 0]
    spectrum = data[:, 1]

    window_size = 7
    window = np.ones(window_size) / window_size
    spectrum_smooth = np.convolve(spectrum, window, mode='same')

    if delay:
        label = f"{intensity_label}, d = {delay} fs"
    else:
        label = intensity_label

    #plt.semilogy(energy, spectrum, ':', color=color, label=f"{label} Raw HHG")
    plt.semilogy(energy, spectrum_smooth, '-', color=color, linewidth=1.5, label=f"{label}")


plt.xlabel("Photon Energy (eV)")
plt.ylabel("HHG Intensity (arb. units)")
plt.xlim(0, 120)
#plt.title("Helium HHG Spectrum (IR Pulse Train)")


"""
# Harmonic lines
plt.vlines(harm, 0, 1e3, 'olive', ':', label='0.83 eV Harmonics')
for n in range(1, 120, 2):
    plt.vlines(harm * n, 0, 1e3, 'olive', ':')
"""

# Key energies

plt.vlines(cutoff_2, 0, 10e3, 'green', '--', label='Calculated Cutoff')
plt.legend(loc='upper center', bbox_to_anchor=(0.78, 1.2))

#plt.legend()
plt.savefig("helium_hhg_xuv_train.png")
plt.show()

"""
#----------------------------------------------------------------------------------
timefile = "EField.he_IR_train"    
vel_col = 1                            

Efield_data = np.loadtxt(timefile, skiprows=1)
time = Efield_data[:,0]
E_col = 3
E_t = Efield_data[:,E_col]
dt = time[1] - time[0]
N = len(time)

# ============================================================
# For Efield Plot
# ============================================================
xuv_25 = 0.195*2.5
IR_25 = 5.067*2.5
delay = 0

# Plot E-field (time domain)
plt.figure(figsize=(8,4))
plt.plot((time/41.34), E_t,  color = 'purple')
plt.xlabel("Time (fs)")
plt.ylabel("E-field (au)")
plt.title("Laser Electric Field: 2.5 cycle ramp on/off IR pulse, 15 cycles")
plt.grid()
plt.savefig('Laser_EField_IR15_train.png')
plt.show()

timefile_2 = "EField.he_xuv_train_d0"    
vel_col = 1                            

Efield_data = np.loadtxt(timefile_2, skiprows=1)
time = Efield_data[:,0]
E_col = 3
E_t = Efield_data[:,E_col]
dt = time[1] - time[0]
N = len(time)

# ============================================================
# For Efield Plot
# ============================================================
xuv_25 = 0.6425
IR_25 = 6.675
delay = 0

plt.figure(figsize=(8,4))
plt.plot((time/41.34), E_t,  color = 'purple')
plt.xlabel("Time (fs)")
plt.ylabel("E-field (au)")
#plt.axvline(0.4875, color='firebrick', linestyle=':')
#plt.axvline(0.4875+7.2, color='firebrick', linestyle=':')
#plt.xlim(-xuv_25, xuv_25)
plt.title(f"Laser Electric Field: 2.5 cycle ramp on/off IR pulse, 15 cycles with d={delay}")
plt.grid()
plt.savefig('Laser_EField_IR15_xuv_train_d0.png')
plt.show()

timefile_3 = "EField.he_xuv_train_dm2"    
vel_col = 1                            

Efield_data = np.loadtxt(timefile_3, skiprows=1)
time = Efield_data[:,0]
E_col = 3
E_t = Efield_data[:,E_col]
dt = time[1] - time[0]
N = len(time)

# ============================================================
# For Efield Plot
# ============================================================
xuv_25 = 0.6425
IR_25 = 6.675
delay = -2

plt.figure(figsize=(8,4))
plt.plot((time/41.34), E_t,  color = 'purple')
plt.xlabel("Time (fs)")
plt.ylabel("E-field (au)")
#plt.axvline(0.4875, color='firebrick', linestyle=':')
#plt.axvline(0.4875+7.2, color='firebrick', linestyle=':')
#plt.xlim(-xuv_25, xuv_25)
plt.title(f"Laser Electric Field: 2.5 cycle ramp on/off IR pulse, 15 cycles with d={delay}")
plt.grid()
plt.savefig('Laser_EField_IR15_xuv_train_dm2.png')
plt.show()
"""