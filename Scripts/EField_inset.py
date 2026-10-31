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

#-------------------------------------------------------------------------------------------

Int = 1e14 #Intensity
c = 3e8 #speed of light 
p_e = 27.212 #eV au conversion
I_p = 24.59 #Ionisation potential 
l_l = 1500e-9 #1500nm laser
U_p = 9.33e16*Int*((l_l/1e9)**2) #pondermotive energy
cutoff_2 = I_p + 3.17*U_p #empirical cutoff formula
f_au = ((c/l_l)/(4.13567e16))*2*np.pi #frequency in atomic units 

timefile = "EField_XUVIR.test"  # Replace with your new EField filename
Efield_data = np.loadtxt(timefile, skiprows=1)

timefile_dm2 = "EField_XUVIR_dm2.test"  # Replace with your new EField filename
Efield_data_dm2 = np.loadtxt(timefile_dm2, skiprows=1)

time = Efield_data[:, 0]

E_XUV = Efield_data[:, 2]  # Third column: XUV (y-direction)
E_IR = Efield_data[:, 3]   # Fourth column: IR (z-direction)
E_XUV_dm2 = Efield_data_dm2[:, 2]  # Third column: XUV (y-direction)

idx_max = np.argmax(np.abs(E_IR))
t_mid = time[idx_max]

time_fs = (time - t_mid) / 41.34


plt.figure(figsize=(11, 5))

plt.plot(time_fs, E_IR, '--', linewidth='2', color="cornflowerblue", label="IR")
plt.plot(time_fs, E_XUV*20,'-.', color="navy", label="XUV, d = 0 fs")
plt.plot(time_fs, E_XUV_dm2*20, color="cyan", label="XUV, d = -2.01 fs ")


plt.axvline(0, color="black", linestyle=":", alpha=0.5)

plt.xlabel("Time (fs)", fontsize=28)
#plt.ylabel("Electric field (a.u.)", fontsize=16)
#plt.title("Combined IR and XUV Electric Fields")
#plt.grid(alpha=0.3)
plt.xticks(fontsize=26)
plt.yticks([])
plt.legend(loc='upper right',fontsize=24)
plt.tight_layout()

plt.savefig("Laser_EField_IR_XUV_dm2.png", dpi=300)

print(f"IR midpoint (a.u.): {t_mid}")

plt.show()


