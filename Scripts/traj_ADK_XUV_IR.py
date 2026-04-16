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


start_index          = 158
index_stride         = 1

def do_integ1(yprime,y,dt):
    y.append(y[-1] - np.trapezoid(yprime[-2:],dx=dt))
    return(y)

def do_integ2(yprime,y,dt):
    y.append(y[-1] + np.trapezoid(yprime[-2:],dx=dt))
    return(y)

def read_field(fieldfile='EField'):
    field=[]
    times=[]
    with open(fieldfile, 'r') as f:
        for line in f.readlines()[1:]:
            field.append(float((line.split()[-1])))
            times.append(float((line.split()[0])))
#    field=[fd/100. for fd in field]
    return times,field

def stage(field,vel,pos,dt):
    vel=do_integ1(field,vel,dt)
    pos=do_integ2(vel,pos,dt)
    return vel,pos

def cross_zero(pos):
    if pos[-1]*pos[-2] < 0.:
#    if np.abs(pos[-1]) < 1 :
        return True
    else:
        return False

def get_energy(vel,Ip):
    return 0.5*vel*vel

def get_field():
    A = 10  # Amplitude
    phi = 0  # (3/5)*np.pi # Phase
    omega = 2  # Frequency
    t = np.linspace(0, 5, 1000)
    E = A*np.sin(omega*t+phi)  # *np.exp(-1*(t-2)**2)
    return t, [x for x in E]

def get_return(index,field,birth_vel,Ip):
    if birth_vel != 0:
        birth_array=[birth_vel, -birth_vel]
    else:
        birth_array=[0]

    for init_vel in birth_array:
        vel= [0 for ii in range(index)]
        pos= [0 for ii in range(index)]
        ff=field[index-2 : index]
        vel[-1]=init_vel

        vel,pos = stage(ff,vel,pos,dt)

        for ii in range(index,len(field)-1):
            ff.append(field[ii+1])
            vel,pos = stage(ff,vel,pos,dt)
            if cross_zero(pos):
                cross_energy= get_energy(vel[-1],Ip)
                if cross_energy>196:
                    print ("excursion length ="+str(max(pos)))
                return cross_energy,ii

    return None, None


def get_local_max(time,array):
    max_list=[]
    time_list=[]
    for i1,i2,i3,tt in zip(array[0:-2],array[1:-1],array[2:],time[1:-1]):
        if i2> i1 and i2 > i3:
            max_list.append(i2)
            time_list.append(tt)
    return(time_list,max_list)



def traj_and_energy(field, index):
    vel= [0 for ii in range(index)]
    pos= [0 for ii in range(index)]
    ff = list(field[index-2:index])


    vel,pos = stage(ff,vel,pos,dt)

    for ii in range(index,len(field)-1):
        ff.append(field[ii+1])
        vel,pos = stage(ff,vel,pos,dt)
        if cross_zero(pos):
            cross_energy= get_energy(vel[-1],0)
            return pos, cross_energy
    return None, None


time, field = read_field(fieldfile='EField.he_IR1.5_xuv_d0')
time_ir, field_ir = read_field(fieldfile='EField.he_IR1.5')

field_total = np.asarray(field)
field_ir = np.asarray(field_ir)
field_xuv = field_total - field_ir

label = 'd = -2 fs'
dt=time[1]-time[0]

traj_list = []
KElist = []

E_h = 5.1422067e11  # atomic unit of electric field (V/m)
euler = np.e
pi = np.pi

# Helium parameters 
Ip_au = 24.587 / 27.211386
Z = 1.0

# Ionisation potentials (eV) 
Ip_ground_eV  = 24.59
Ip_excited_eV = 24.59 - 21.23  
Ip_ground_au  = Ip_ground_eV / 27.211386
Ip_excited_au = Ip_excited_eV / 27.211386

au_to_s = 2.418884e-17


def ADK_rate_au(E, Ip):
    """
    ADK rate in atomic units for a single field value and Ip.
    """
    E = abs(E)

    if E < 1e-12:
        return 0.0

    kappa = np.sqrt(2 * Ip)

    rate = (4 * kappa**3 / E) * np.exp(-2 * kappa**3 / (3 * E))

    return rate


def ADK_probability(E_au, Ip_au, dt_au):
    """
    Converts ADK rate into timestep ionisation probability.
    """
    rate_au = ADK_rate_au(E_au, Ip_au)
    rate_SI = rate_au / au_to_s
    dt_SI   = dt_au * au_to_s

    return 1 - np.exp(-rate_SI * dt_SI)

prob_list = []

xuv_threshold = 0.01 * np.max(np.abs(field_xuv))
xuv_weight = np.abs(field_xuv) / np.max(np.abs(field_xuv))

for ii in range(start_index, len(field)-2):

    #  Determine effective Ip 
    if np.abs(field_xuv[ii]) > xuv_threshold:
        Ip_eff = Ip_excited_au
    else:
        Ip_eff = Ip_ground_au

    P = ADK_probability(field_ir[ii], Ip_eff, dt)
    P *= xuv_weight[ii]

    xx, en = traj_and_energy(field_ir, ii)

    if xx:
        traj_list.append(xx)
        KElist.append(en)
        prob_list.append(P)
    else:
        traj_list.append(np.zeros(10))
        KElist.append(0)
        prob_list.append(0)


print(ii, P)
print(max(KElist))

#  unit conversion for plotting 
IR_25_fs = 2.5 * 5.067
xuv_25 = 0.195*2.5
time_fs = (np.asarray(time) / 41.43) - IR_25_fs
delay = 0

KElist_eV = np.asarray(KElist) * 27.211386

W_array = []

for i in range(len(field_ir)):
    if np.abs(field_xuv[i]) > xuv_threshold:
        Ip_eff = Ip_excited_au
    else:
        Ip_eff = Ip_ground_au

    W_array.append(ADK_rate_au(field_ir[i], Ip_eff))

W_array = np.array(W_array)

fig, ax1 = plt.subplots()

#  Return Energy 
index = start_index
for en in KElist_eV:
    ax1.plot(time_fs[index], en, '.', 
             color='deepskyblue', markersize=3)
    index += 1

ax1.set_xlabel('Birth Time (fs)')
ax1.set_ylabel('Return Energy (eV)', color='deepskyblue')
ax1.tick_params(axis='y', labelcolor='deepskyblue')

#  Electric Field 
ax2 = ax1.twinx()
ax2.plot(time_fs, field, color='purple', alpha=0.6)
ax2.set_ylabel('E-Field (a.u.)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

x1 = delay - (xuv_25)
x2 = x1 + (xuv_25 * 2)

plt.axvspan(x1, x2, color='limegreen', alpha=0.2)

#  ADK rate (normalised) 
ax3 = ax1.twinx()
ax3.spines.right.set_position(("axes", 1.15))

W_norm = W_array / np.max(W_array)
ax3.plot(time_fs, W_norm, '--', color='crimson', linewidth=1)
ax3.set_ylabel('Normalised ADK Rate', color='crimson')
ax3.tick_params(axis='y', labelcolor='crimson')
#ax3.axhline(max(W_norm), color='red', linewidth=0.5)

ax = plt.gca()

# 32 eV
x1 = -7.7
#ax.axvline(x1, linestyle=':', color='black')
ax.text(x1, 1.02, 'Peak 1',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 54.5 eV (slightly higher)
x2 = -5
#ax.axvline(x2, linestyle=':', color='black')
ax.text(x2, 1.08, 'Peak 2',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 58 eV (slightly lower)
x3 = -2.2
#ax.axvline(x3, linestyle=':', color='black')
ax.text(x3, 1.02, 'Peak 3',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

# 68.5 eV
x4 = 0.2
#ax.axvline(x4, linestyle=':', color='black')
ax.text(x4, 1.08, 'Peak 4',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom')

#plt.title('1500 nm, 10$^{14}$ W/cm$^2$ — Classical Return Energies and ADK Window')
plt.savefig("helium_hhg_ADK_XUV_peak_d0.png")
plt.show()

print("Peak field (a.u.):", max(np.abs(field)))
print("Peak ADK rate (a.u.):", max(W_array))
print("Peak ADK rate (s^-1):", max(W_array)/au_to_s)


