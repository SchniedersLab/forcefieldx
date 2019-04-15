### L-BFGS Minimization

Run L-BFGS minimization on a system.

---
```
Usage: ffxc Minimize [-h] [--af=-1] [--as=1] [--ef1=-1] [--ef2=-1] [--es1=1] [--es2=1] [--f1=-1] [--f2=-1] [--la1=-1] [--la2=-1] [--np=1] [--s1=0] [--s2=0] [--sf=1.0] [--uaA=-1] [--uaB=-1] [-e=1.0] [-I=Unlimited] [-l=-1] files...
Run L-BFGS minimization on a system.
    files...                        Atomic coordinate files in PDB or XYZ format.
    --af, --activeFinal=-1          Final active atom (single-topology only).
    --as, --activeStart=1           Starting active atom (single-topology only).
    --ef1, --noElecFinal1=-1        Final no-electrostatics atom for 1st topology.
    --ef2, --noElecFinal2=-1        Final no-electrostatics atom for 2nd topology
    --es1, --noElecStart1=1         Starting no-electrostatics atom for 1st topology.
    --es2, --noElecStart2=1         Starting no-electrostatics atom for 2nd topology
    --f1, --final1=-1               Final ligand atom for the 1st topology.
    --f2, --final2=-1               Final ligand atom for the 2nd topology
    --la1, --ligAtoms1=-1           Period-separated ranges of 1st topology ligand atoms (e.g. 40-50.72-83).
    --la2, --ligAtoms2=-1           Period-separated ranges of 2nd toplogy ligand atoms (e.g. 40-50.72-83)
    --np, --nParallel=1             Number of topologies to evaluate in parallel
    --s1, --start1=0                Starting ligand atom for 1st topology.
    --s2, --start2=0                Starting ligand atom for 2nd topology
    --sf, --switchingFunction=1.0   Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)
    --uaA, --unsharedA=-1           Unshared atoms in the A dual topology (period-separated hyphenated ranges)
    --uaB, --unsharedB=-1           Unshared atoms in the B dual topology (period-separated hyphenated ranges)
-e, --eps=1.0                       Convergence criteria.
-h, --help                          Print this help message.
-I, --iterations=Unlimited          Number of minimization steps.
-l, --lambda=-1                     Initial lambda value.
 ```
 ---
 