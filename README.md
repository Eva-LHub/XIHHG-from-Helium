# XIHHG from Helium - Research Data & Code Scripts

This repository contains raw experimental data and analysis code for:

**"Strongly Enhanced and Controllable Plateau in Extreme Ultraviolet Initiated High Harmonic Generation from Helium"**

Eva Loughridge, QUB MSci Thesis April 2026

## Data Directory
Raw data generated using RMT (see https://gitlab.com/Uk-amor/RMT/rmt) on the Kelvin2 high-performance computing facility

- `data/IR Only/` - EField and Harmonic Velocity Dipole of the IR Laser Only
- `data/XUV Only/` - EField and Harmonic Velocity Dipole of the XUV Laser Only
- `data/XUV-IR/` - EField and Harmonic Velocity Dipole of the combined XUV-IR Laser field 
   - `data/XUV-IR/Delays/` - EField and Harmonic Velocity Dipole for various time delays (time elapsed between the peak of the IR and the peak of the XUV laser pulse).
   - `data/XUV-IR/Detune/` - Harmonic Velocity Dipole for various time XUV energies that are detuned from the 1s2p state in Helium
   - `data/XUV-IR/Intensities/` - Harmonic Velocity Dipole for various intensities of the XUV laser
   - `data/XUV-IR/Pulse-Train/` - EField and Harmonic Velocity Dipole for an XUV-IR pulse train

## Analysis & Visualisation Scripts

Four visualisation/analysis scripts for generating thesis figures:

- `Scripts/HHG_XUV_IR.py` - Standard Script for producing HHG Spectrum, including various delays, detuned XUV energies and XUV intensities
- `Scripts/traj_ADK_XUV_IR.py` - Return Energies Trajectory Script with modified ADK Ionisation Potential
- `Scripts/HHG_XUV_IR_Delays.py`- Colour-map of the effect of the time delay on the HHG spectrum
- `Scripts/HHG_Pulse_Train.py`- HHG Spectrum for a Pulse Train
- Momentum distribution figures were generated using the `RMT_plot_mom` function from the RMT code, executed on the Kelvin2 system.


See the thesis for figure references and detailed analysis
