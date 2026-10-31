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
    #(data_xuv_IR15_d0,  'navy',  'XUV, d = 0 fs', 0, 'dashdot'),
  
    #(data_xuv_IR15_dm2, 'aqua', 'XUV', -2.01, '-'),
    #(data_xuv_IR15_a181_dm13, 'magenta', 'XUV', -13.16, (0, (1, 1))),
    
    (data_xuv_IR15_x1,    'red',      r'XUV, $I = 10^{12}\,\mathrm{{W\,cm^{{-2}}}}$', 0, '--'),
    (data_xuv_IR15_d0,   'navy',      r'XUV, $I = 10^{11}\,\mathrm{W\,cm^{-2}}$', 0, (0, (1, 1))),
    (data_xuv_IR15_x01,  'lawngreen', r'XUV, $I = 10^{10}\,\mathrm{W\,cm^{-2}}$', 0, 'dashdot'),
    (data_xuv_IR15_x001, 'orange',    r'XUV, $I = 10^{9}\,\mathrm{W\,cm^{-2}}$', 0, (0, (3, 3, 1, 3, 1, 3))),
    (data_xuv_IR15_x0001,'violet',    r'XUV, $I = 10^{8}\,\mathrm{W\,cm^{-2}}$', 0, (0, (3, 4, 1, 4, 1, 4, 1, 4))),
]
plt.figure(figsize=(8, 5))


# -------------------------------------------------------------------------
def harmonic_yield(data, centre_energy, width=1.0):

    energy = data[:,0]
    spectrum = data[:,1]

    mask = (
        (energy >= centre_energy - width/2)
        &
        (energy <= centre_energy + width/2)
    )

    dE = energy[1] - energy[0]

    return np.sum(spectrum[mask]) * dE

def integrated_plateau_yield(data, Ip, cutoff):

    energy = data[:,0]
    spectrum = data[:,1]

    mask = (energy >= Ip) & (energy <= cutoff)

    dE = energy[1] - energy[0]

    return np.sum(spectrum[mask]) * dE

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(8,8),
    gridspec_kw={'height_ratios':[3,1]},
    sharex=False
)

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

ax1.axvline(50, color='grey', linestyle=':', label='50 eV')
ax1.set_ylabel("HHG Intensity (arb. units)")
ax1.set_xlim(0,100)
ax1.set_ylim(1e-12,1e4)

ax1.legend(loc='upper right')

# BOTTOM PANEL: integrated yield vs XUV energy

intensity_range = np.array([
    1, 0.1, 0.01, 0.001, 0.0001
])

intensity_data = [
    data_xuv_IR15_x1,
    data_xuv_IR15_d0,
    data_xuv_IR15_x01,
    data_xuv_IR15_x001,
    data_xuv_IR15_x0001
]

harmonic_energy = 50.0      
integration_width = 1.0     

integrated_yields = []

for data in intensity_data:

    Y = harmonic_yield(
        data,
        centre_energy=harmonic_energy,
        width=integration_width
    )

    integrated_yields.append(Y)

integrated_yields = np.array(integrated_yields)

ax2.plot(
    intensity_range,
    integrated_yields,
    'o-',
    color='purple',
    linewidth=1.5,
    markersize=5
)

ax2.set_xlabel("XUV Intensity (%)")
ax2.set_ylabel(f"Yield at {harmonic_energy:.0f} eV")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.grid(alpha=0.3)

ax2.set_xlabel("Intensity (%)")
ax2.set_ylabel("Plateau Yield")
ax2.set_yscale('log')

ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("helium_hhg_IR15_xuv_int_yield.png")

plt.show()