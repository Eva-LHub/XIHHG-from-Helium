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
data_xuv_IR15_d0 = np.loadtxt("Harm_vel_0001_z_xuv_IR1.5")  # f = 0.78 au, d = 0 fs, I = 0.1%
data_xuv_IR15_dm2 = np.loadtxt("Harm_vel_0001_z_dm2")       # f = 0.78 au, d = -2 fs

# Long Pulse
data_IR_pulse = np.loadtxt("Harm_vel_0001_z_train")
data_xuv_pulse_d0 = np.loadtxt("Harm_vel_0001_z_train_d0")
data_xuv_pulse_dm2 = np.loadtxt("Harm_vel_0001_z_train_dm2")

# -------------------------------------------------------------------------------------------

Int = 1e14  # Intensity
c = 3e8     # speed of light
p_e = 27.212  # eV au conversion
I_p = 24.59   # Ionisation potential
l_l = 1500e-9  # 1500 nm laser

U_p = 9.33e16 * Int * ((l_l / 1e9) ** 2)  # ponderomotive energy
cutoff_2 = I_p + 3.17 * U_p               # empirical cutoff formula

f_au = ((c / l_l) / (4.13567e16)) * 2 * np.pi  # frequency in atomic units
harm = f_au * p_e                              # Harmonics

# -------------------------------------------------------------------------------------------

datasets = [

    (data_IR_pulse, 'cornflowerblue',
     'Long pulse: IR only, 1.5 $\\mu$m', 0, 'dashdot'),

    (data_xuv_pulse_d0, 'navy',
     'Long pulse: XUV, d = 0 fs', 0, '-'),

    (data_xuv_pulse_dm2, 'aqua',
     'Long pulse: XUV', -2.01, (0, (1, 1))),

    (data_IR15, 'lawngreen',
     'Short pulse: IR only, 1.5 $\\mu$m', 0, '--'),

    (data_xuv_IR15_d0, 'orange',
     'Short pulse: XUV, d = 0 fs', 0, '-'),

    (data_xuv_IR15_dm2, 'red',
     'Short pulse: XUV', -2.01, (0, (1, 1))),
]

fig, ax = plt.subplots(figsize=(8, 5))

line_handles = {}

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

    line, = ax.semilogy(
        energy,
        spectrum_smooth,
        linestyle=linestyle,
        color=color,
        linewidth=1.5,
        label=label
    )

    line_handles[label] = line

ax.set_xlabel("Photon Energy (eV)")
ax.set_ylabel("HHG Intensity (arb. units)")
ax.set_xlim(0, 120)

"""
# Harmonic lines
ax.vlines(harm, 0, 1e3, 'olive', ':', label='0.83 eV Harmonics')
for n in range(1, 120, 2):
    ax.vlines(harm * n, 0, 1e3, 'olive', ':')
"""
cutoff_line = ax.vlines(
    cutoff_2,
    0,
    1e4,
    colors='green',
    linestyles='--',
    label='Calculated Cutoff'
)


xuv_legend = ax.legend(
    handles=[
        line_handles['Long pulse: XUV, d = 0 fs'],
        line_handles['Long pulse: XUV, d = -2.01 fs'],
        line_handles['Short pulse: XUV, d = 0 fs'],
        line_handles['Short pulse: XUV, d = -2.01 fs']
    ],
    loc='upper right'
)

ax.add_artist(xuv_legend)


ax.legend(
    handles=[
        line_handles['Long pulse: IR only, 1.5 $\\mu$m'],
        line_handles['Short pulse: IR only, 1.5 $\\mu$m'],
        cutoff_line
    ],
    loc='lower left'
)

plt.tight_layout()
plt.savefig("helium_hhg_xuv_long_pulse.png")
plt.show()
