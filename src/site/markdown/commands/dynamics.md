### Molecular Dynamics

Run NVE, NVT, NPT or stochastic dynamics.

---
```
Usage: ffxc Dynamics [-hox] [-b=Bussi] [-d=1.0] [-F=XYZ] [-i=Verlet] [-k=1.0] [-n=1000000] [-p=0] [-r=0.25] [-t=298.15] [-w=10.0] files...
Run dynamics on a system.
    files...                    XYZ or PDB input files.
-b, --thermostat=Bussi          Thermostat: [Adiabatic / Berendsen / Bussi].
-d, --dt=1.0                    Time discretization step in femtoseconds.
-F, --fileFormat=XYZ            Choose file type to write [PDB/XYZ].
-h, --help                      Print this help message.
-i, --integrator=Verlet         Integrator: [Beeman / Respa / Stochastic / Verlet].
-k, --checkpoint=1.0            Interval to write out restart files (.dyn, .his, etc).
-n, --numberOfSteps=1000000     Number of molecular dynamics steps.
-o, --optimize                  Optimize and save low-energy snapshots.
-p, --npt=0                     Specify use of a MC Barostat at the given pressure; the default 0 disables NPT (atm).
-r, --report=0.25               Interval to report thermodynamics (psec).
-t, --temperature=298.15        Temperature (Kelvin).
-w, --write=10.0                Interval to write out coordinates (psec).
-x, --repEx                     Execute temperature replica exchange.
```
---
