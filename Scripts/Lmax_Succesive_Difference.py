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

harm = 0.057*27.212
I_p= 24.59
cutoff = 43.5
c = 3e8
l= 1800e6
f = c/l


file_info = {
    121: "Harm_vel_0001_z_a121",
    131: "Harm_vel_0001_z_a131",
    141: "Harm_vel_0001_z_a141",
    151: "Harm_vel_0001_z_a151",
    161: "Harm_vel_0001_z_a161",
    171: "Harm_vel_0001_z_a171",    
    181: "Harm_vel_0001_z_a181",
    191: "Harm_vel_0001_z_a191",
    201: "Harm_vel_0001_z_a201",
    }


data_vel = {L: np.loadtxt(fname) for L, fname in file_info.items()}

ref_L = 201
E_min, E_max = 0.0, cutoff      # eV window for the averages
rel_thresh = 1e-6               # exclude points below this fraction of window max
highlight_L = 181               # the value adopted in the paper

Ls_all = sorted(data_vel)
energy = data_vel[ref_L][:, 0]
    
regions = {
    "0-24.6 eV":     (0.0, I_p),
    "24.6-43.5 eV":  (I_p, cutoff),
    "All (0-43.5 eV)": (0.0, 43.5),
}

def regional_diff(I_lo, I_hi):
    valid = np.isfinite(I_hi) & (I_hi > 0)
    out = {}
    for name, (lo, hi) in regions.items():
        sel = valid & (energy >= lo) & (energy < hi)
        pct = 100 * np.abs(I_lo[sel] - I_hi[sel]) / I_hi[sel]
        out[name] = {"mean": pct.mean(), "median": np.median(pct), "max": pct.max()}
    return out

succ_reg = {}
for L_lo, L_hi in zip(Ls_all[:-1], Ls_all[1:]):
    succ_reg[L_hi] = regional_diff(data_vel[L_lo][:, 1], data_vel[L_hi][:, 1])

print(f"{'Pair':>10} " + " ".join(f"{n:>15}" for n in regions))
for L in sorted(succ_reg):
    print(f"{f'{L-10}->{L}':>10} " + " ".join(f"{succ_reg[L][n]['mean']:15.3f}" for n in regions))

fig, ax = plt.subplots(figsize=(7, 4))
for (name, col, mk) in [("0-24.6 eV", "cornflowerblue", "o"),
                        ("24.6-43.5 eV", "magenta", "s"),
                        ("All (0-43.5 eV)", "k", "d")]:
    Ls = sorted(succ_reg)
    ax.semilogy(Ls, [succ_reg[L][name]["mean"] for L in Ls],color=col, marker=mk, label=name)
    
ax.axvline(highlight_L, color='green', linestyle='--', alpha=0.6)
ax.set_xlabel(r"Maximum angular momentum $L$")
ax.set_ylabel("Mean relative difference (%)")
ax.legend()
plt.tight_layout()
plt.savefig("successive_difference.png")
plt.show()
    
