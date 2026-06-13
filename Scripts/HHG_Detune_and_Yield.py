#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 23:29:45 2026

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
data_xuv_detune_2 = np.loadtxt("Harm_vel_0001_z_2eV")
data_xuv_detune_4 = np.loadtxt("Harm_vel_0001_z_4eV")
data_xuv_detune_6 = np.loadtxt("Harm_vel_0001_z_6eV")
data_xuv_detune_8 = np.loadtxt("Harm_vel_0001_z_8eV")
data_xuv_detune_10 = np.loadtxt("Harm_vel_0001_z_10eV")
data_xuv_detune_12 = np.loadtxt("Harm_vel_0001_z_12eV")
data_xuv_detune_14 = np.loadtxt("Harm_vel_0001_z_14eV")
data_xuv_detune_16 = np.loadtxt("Harm_vel_0001_z_16.06eV") 
data_xuv_detune_18 = np.loadtxt("Harm_vel_0001_z_18eV")
data_xuv_detune_20 = np.loadtxt("Harm_vel_0001_z_20eV")
data_xuv_detune_22 = np.loadtxt("Harm_vel_0001_z_22eV")
data_xuv_detune_24 = np.loadtxt("Harm_vel_0001_z_24eV")

data_xuv_detune_5 = np.loadtxt("Harm_vel_0001_z_4.89eV")
data_xuv_detune_40 = np.loadtxt("Harm_vel_0001_z_40eV")

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
    (data_IR15,  'cornflowerblue',  'IR only, 1.5 $\mu$m', 0, '-'),
    #(data_xuv_only,  'orange',  'XUV only, 21.23 eV', 0, (0, (3, 3, 1, 3, 1, 3))),
    (data_xuv_IR15_d0,  'navy',  'XUV, d = 0 fs', 0, 'dashdot'),
  
    #(data_xuv_IR15_dm2, 'aqua', 'XUV', -2.01, '-'),
    #(data_xuv_IR15_a181_dm13, 'magenta', 'XUV', -13.16, (0, (1, 1))),
    
    (data_xuv_detune_5, 'lawngreen', 'Detuned 5 eV XUV', 0, (0, (3, 3, 1, 3, 1, 3))),
    (data_xuv_detune_16, 'red', 'Detuned 16 eV XUV', 0, (0, (1, 1))),
    (data_xuv_detune_40, 'orange', 'Detuned 40 eV XUV', 0, '--'),
    
    #(data_xuv_IR15_x1,  'red',  'XUV, Int. = 1%', 0, '--'),
    #(data_xuv_IR15_d0,  'navy',  'XUV, Int. = 0.1%', 0, (0, (1, 1))),
    #(data_xuv_IR15_x01,  'lawngreen',  'XUV, Int. = 0.01%', 0, 'dashdot'),
    #(data_xuv_IR15_x001,  'orange',  'XUV, Int. = 0.001%', 0, (0, (3, 3, 1, 3, 1, 3))),
    #(data_xuv_IR15_x0001,  'violet',  'XUV, Int. = 0.0001%', 0, (0, (3, 4, 1, 4, 1, 4, 1, 4))),
]
plt.figure(figsize=(8, 5))

# -------------------------------------------------------------------------

def integrated_plateau_yield(data, Ip, cutoff):

    energy = data[:,0]
    spectrum = data[:,1]

    mask = (energy >= Ip) & (energy <= cutoff)

    dE = energy[1] - energy[0]

    return np.sum(spectrum[mask]) * dE

# -------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(8,8),
    gridspec_kw={'height_ratios':[3,1]},
    sharex=False
)

# -------------------------------------------------------------------------
# TOP PANEL: spectra

for data, color, intensity_label, delay, linestyle in datasets:

    energy = data[:,0]
    spectrum = data[:,1]

    window_size = 7
    window = np.ones(window_size)/window_size
    spectrum_smooth = np.convolve(
        spectrum,
        window,
        mode='same'
    )

    if delay:
        label = f"{intensity_label}, d = {delay} fs"
    else:
        label = intensity_label

    ax1.semilogy(
        energy,
        spectrum_smooth,
        color=color,
        linestyle=linestyle,
        linewidth=1.5,
        label=label
    )

ax1.set_ylabel("HHG Intensity (arb. units)")
ax1.set_xlim(0,100)
ax1.set_ylim(1e-12,1e4)

ax1.legend(loc='upper right')

# -------------------------------------------------------------------------
# BOTTOM PANEL: integrated yield vs XUV energy

detuning_energies = np.array([
    2,
    4,
    4.89,
    6,
    8,
    10,
    12,
    14,
    16.06,
    18,
    20,
    22,
    24
])

detuning_datasets = [
    data_xuv_detune_2,
    data_xuv_detune_4,
    data_xuv_detune_5,
    data_xuv_detune_6,
    data_xuv_detune_8,
    data_xuv_detune_10,
    data_xuv_detune_12,
    data_xuv_detune_14,
    data_xuv_detune_16,
    data_xuv_detune_18,
    data_xuv_detune_20,
    data_xuv_detune_22,
    data_xuv_detune_24
]

integrated_yields = []

for data in detuning_datasets:

    Y = integrated_plateau_yield(
        data,
        I_p,
        cutoff_2
    )

    integrated_yields.append(Y)

integrated_yields = np.array(integrated_yields)

ax2.plot(detuning_energies, integrated_yields,'o-',color='purple', linewidth=1.5, markersize=5)

ax2.set_xlabel("XUV Photon Energy (eV)")
ax2.set_ylabel("Plateau Yield")
ax2.set_yscale('log')

ax2.grid(alpha=0.3)

# -------------------------------------------------------------------------

plt.tight_layout()

plt.savefig("helium_hhg_IR15_xuv_detune_yield.png")

plt.show()