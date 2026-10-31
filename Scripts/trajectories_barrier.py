"""Classical return energies vs birth time for the five-cycle 1.5 um IR pulse used in the
RMT calculations, with the XUV pulse windows for d = 0 and d = -2.01 fs and the field
strength above which the 1s2p electron is released over the barrier.

Field: E(t) = E0 sin^2(pi t'/Tp) cos(w t), with t' measured from the start of the pulse and
t from the centre of the envelope; this reproduces the return energies in the original
trajectory figure (Peak 3: 68 eV, Peak 4: 54 eV)."""
import numpy as np
import matplotlib.pyplot as plt

AU_FS = 2.4188843e-2
HARTREE = 27.211386
w = 0.030                      # IR frequency (a.u.)
E0 = np.sqrt(1e14 / 3.50945e16)  # peak field for 1e14 W/cm^2 (a.u.)
Tcyc = 2 * np.pi / w
Tp = 5 * Tcyc                  # five-cycle sin^2 pulse
Ip_gs = 24.587 / HARTREE
Ip_2p = 3.36 / HARTREE
F_bsi_2p = Ip_2p**2 / 4        # barrier-suppression field, Z = 1
F_bsi_gs = Ip_gs**2 / 4
w_xuv = 0.78
T_xuv = 5 * 2 * np.pi / w_xuv  # 0.98 fs XUV pulse

dt = 0.05
tp = np.arange(0, Tp + dt, dt)               # time from pulse start
t = tp - Tp / 2                               # time relative to envelope centre
E = E0 * np.sin(np.pi * tp / Tp)**2 * np.cos(w * t)
A = -np.cumsum(E) * dt                        # vector potential, E = -dA/dt

# classical trajectories: born at rest at x=0, v(t) = A(t) - A(t0)
births, energies = [], []
for i0 in range(0, len(tp) - 1, 4):
    v = A[i0:] - A[i0]
    x = np.cumsum(v) * dt
    cross = np.where(np.sign(x[1:]) != np.sign(x[:-1]))[0]
    cross = cross[cross > 2]
    if cross.size:
        j = cross[0]
        births.append(t[i0]); energies.append(0.5 * v[j]**2 * HARTREE)
births = np.array(births) * AU_FS
energies = np.array(energies)

if __name__ == '__main__':
    print(f"E0 = {E0:.4f} a.u.; F_BSI(1s2p) = {F_bsi_2p:.4f} a.u. "
          f"({F_bsi_2p**2*3.50945e16:.2e} W/cm2); F_BSI(1s^2) = {F_bsi_gs:.3f} a.u. "
          f"({F_bsi_gs**2*3.50945e16:.2e} W/cm2)")
    print(f"Ponderomotive energy Up = {E0**2/(4*w**2)*HARTREE:.2f} eV; half cycle = {Tcyc/2*AU_FS:.3f} fs")
    # report lobes
    lobes = []
    for c in np.arange(-12, 12, 2.533 / 1):
        m = (births > c - 1.27) & (births < c + 1.27)
        if m.any() and energies[m].max() > 1:
            k = np.argmax(np.where(m, energies, -1))
            lobes.append((births[k], energies[k]))
    for b, e in lobes:
        print(f"lobe: max return {e:5.1f} eV at birth {b:6.2f} fs -> cutoff {e + Ip_gs*HARTREE:5.1f} eV")
    frac = np.mean(np.abs(E[np.abs(t * AU_FS) < 2.53 * 1.5]) > F_bsi_2p)
    print(f"fraction of time |E|>F_BSI(1s2p) within +/-3.8 fs of centre: {frac:.3f}")

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax2 = ax.twinx()
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    for d, col, lab in [(0.0, '#2a9d4b', r'XUV, $d=0$ fs'), (-2.01, '#e08a1e', r'XUV, $d=-2.01$ fs')]:
        ax2.axvspan(d - T_xuv / 2 * AU_FS, d + T_xuv / 2 * AU_FS, color=col, alpha=0.25, lw=0, label=lab)
    ax2.plot(t * AU_FS, E, color='#7b3f9e', lw=1.4, alpha=0.8, label='IR field')
    ax2.axhspan(-F_bsi_2p, F_bsi_2p, color='0.5', alpha=0.25, lw=0,
                label=r'$|E|<F_\mathrm{BSI}$ for $1s2p$')
    ax2.set_ylabel('IR electric field (a.u.)', color='#7b3f9e')
    ax2.tick_params(axis='y', colors='#7b3f9e')
    ax2.set_ylim(-0.065, 0.065)
    ax.plot(births, energies, '.', ms=2.5, color='#1f78c8', label='Return energy')
    ax.set_xlabel('Birth time (fs)')
    ax.set_ylabel('Return energy (eV)', color='#1f78c8')
    ax.tick_params(axis='y', colors='#1f78c8')
    ax.set_ylim(0, 80)
    ax.set_xlim(-12.7, 12.7)
    top4 = sorted(lobes, key=lambda x: -x[1])[:4]
    ycoord = top4[0][1] + 5
    for n, (b, e) in enumerate(sorted(top4), start=1):
        ax.annotate(f'Peak {n}', (b, ycoord), textcoords='offset points', xytext=(0, 5),
                    ha='center', fontsize=8)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=7.5, frameon=False, markerscale=4)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'Images/helium_traj_barrier.{ext}', dpi=200)
