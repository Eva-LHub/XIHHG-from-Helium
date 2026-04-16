# Data Directory

## Structure
Raw data generated using RMT (see https://gitlab.com/Uk-amor/RMT/rmt)

- `IR Only/` - EField and Harmonic Velocity Dipole of the IR Laser Only
- `XUV Only/` - EField and Harmonic Velocity Dipole of the XUV Laser Only
- `XUV-IR/` - EField and Harmonic Velocity Dipole of the combined XUV-IR Laser field 
   - `Delays/` - EField and Harmonic Velocity Dipole for various time delays (time elapsed between the peak of the IR and the peak of the XUV laser pulse).
   - `Detune/` - Harmonic Velocity Dipole for various time XUV energies that are detuned from the 1s2p state in Helium
   - `Intensities/` - Harmonic Velocity Dipole for various intensities of the XUV laser
   - `Pulse-Train/` - EField and Harmonic Velocity Dipole for an XUV-IR pulse train
