### Simulated Annealing

Run molecular dynamics at a series of temperatures to optimize a structure.

---
```
Usage: ffxc Anneal [-ho] [--tl=10.0] [--tu=1000.0] [-b=Bussi] [-d=1.0] [-i=Verlet] [-k=1.0] [-n=1000000] [-r=0.25] [-t=298.15] [-w=10.0] [-W=10] files...
Run simulated annealing on a system.
     files...                           XYZ or PDB input files.
     --tl, --temperatureLow=10.0        Low temperature limit (Kelvin).
     --tu, --temperatureUpper=1000.0    High temperature limit (Kelvin).
 -b, --thermostat=Bussi                 Thermostat: [Adiabatic / Berendsen / Bussi].
 -d, --dt=1.0                           Time discretization step in femtoseconds.
 -h, --help                             Print this help message.
 -i, --integrator=Verlet                Integrator: [Beeman / Respa / Stochastic / Verlet].
 -k, --checkpoint=1.0                   Interval to write out restart files (.dyn, .his, etc).
 -n, --numberOfSteps=1000000            Number of molecular dynamics steps.
 -o, --optimize                         Optimize and save low-energy snapshots.
 -r, --report=0.25                      Interval to report thermodynamics (psec).
 -t, --temperature=298.15               Temperature (Kelvin).
 -w, --write=10.0                       Interval to write out coordinates (psec).
 -W, --windows=10                       Number of annealing windows.
 ```
 ---
 