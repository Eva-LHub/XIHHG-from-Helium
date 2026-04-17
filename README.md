# XIHHG from Helium - Research Data & Code Scripts

This repository contains raw experimental data and analysis code for:

**"Strongly Enhanced and Controllable Plateau in Extreme Ultraviolet Initiated High Harmonic Generation from Helium"**

Eva Loughridge, QUB MSci Thesis April 2026

## Data Directory
Raw data generated using RMT (see https://gitlab.com/Uk-amor/RMT/rmt) on the Kelvin2 high-performance computing facility. 

RMT provides the time-dependent wavefunction, from which can be computed the time-dependent expectation value of the dipole, d(t). The harmonic spectra are then calculated as the square modulus of the Fourier transformed dipole, d(ω). In practice, RMT produces both the length and velocity forms of the dipole expectation value, and the two can be compared to verify the accuracy of the calculation. All spectra shown in this paper use the velocity form, but the results are the same from either form. 

- `data/IR Only/` - EField and Harmonic Spectrum of the IR Laser Only
- `data/XUV Only/` - EField and Harmonic Spectrum of the XUV Laser Only
- `data/XUV-IR/` - EField and Harmonic Spectra of the combined XUV-IR Laser field 
   - `data/XUV-IR/Delays/` - EField and Harmonic Spectra for various time delays (time elapsed between the peak of the IR and the peak of the XUV laser pulse).
   - `data/XUV-IR/Detune/` - Harmonic Spectra for various time XUV energies that are detuned from the 1s2p state in Helium
   - `data/XUV-IR/Intensities/` - Harmonic Spectra for various intensities of the XUV laser
   - `data/XUV-IR/Pulse-Train/` - EField and Harmonic Spectra for an XUV-IR pulse train

## Analysis & Visualisation Scripts

Four visualisation/analysis scripts for generating thesis figures:

- `Scripts/HHG_XUV_IR.py` - Standard Script for visualising HHG Spectra, including various delays, detuned XUV energies and XUV intensities
- `Scripts/traj_ADK_XUV_IR.py` - Return Energies Trajectory Script with modified ADK Ionisation Potential
- `Scripts/HHG_XUV_IR_Delays.py`- Colour-map Script for the effect of the time delay on the HHG spectrum
- `Scripts/HHG_Pulse_Train.py`- Script for visualising the HHG Spectrum of a Pulse Train
- Momentum distribution figures were generated using the `RMT_plot_mom` function from the RMT code, executed on the Kelvin2 system.


See the thesis for figure references and detailed analysis.
