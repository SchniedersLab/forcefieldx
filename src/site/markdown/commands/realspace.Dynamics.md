### Real Space Molecular Dynamics

Molecular dynamics on a Real Space target.

---
```
Usage: ffxc realspace.Dynamics [-ho] [--wA=1.0] [-b=Bussi] [-d=1.0] [-i=Verlet] [-k=1.0] [-n=1000000] [-r=0.25] [-t=298.15] [-w=10.0] [-X=<data> <data>]... files...
Molecular dynamics on a Real Space target.
    files...                    PDB and Real Space input files.
    --wA, --dataWeight=1.0      The weight of the real space data (wA).
-b, --thermostat=Bussi          Thermostat: [Adiabatic / Berendsen / Bussi].
-d, --dt=1.0                    Time discretization step in femtoseconds.
-h, --help                      Print this help message.
-i, --integrator=Verlet         Integrator: [Beeman / Respa / Stochastic / Verlet].
-k, --checkpoint=1.0            Interval to write out restart files (.dyn, .his, etc).
-n, --numberOfSteps=1000000     Number of molecular dynamics steps.
-o, --optimize                  Optimize and save low-energy snapshots.
-r, --report=0.25               Interval to report thermodynamics (psec).
-t, --temperature=298.15        Temperature (Kelvin).
-w, --write=10.0                Interval to write out coordinates (psec).
-X, --data=<data> <data>        Specify input data filename, and its weight (wA) (e.g. -X filename 1.0).
```
---