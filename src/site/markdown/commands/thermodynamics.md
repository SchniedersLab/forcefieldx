### Calculate Thermodynamic Differences

Use Tempered Orthogonal Space Random Walk algorithm to estimate a free energy difference.

---
```
Usage: ffxc Thermodynamics [-hoy] [--mc] [--rn] [--ts] [--af=-1] [--as=1] [--bM=0.05] [--dw=OFF] [--ef1=-1] [--ef2=-1] [--es1=1] [--es2=1] [--f1=-1] [--f2=-1] [--la1=-1] [--la2=-1] [--lf=1.0E-18] [--lm=1.0E-18] [--lw=1.0] [--mcL=0.1] [--mcMD=100] [--np=1] [--rsym=-1.0] [--ruc=-1.0] [--s1=0] [--s2=0] [--sf=1.0] [--tp=8.0] [--uaA=-1] [--uaB=-1] [-b=Bussi] [-C=10] [-d=1.0] [-F=XYZ] [-i=Verlet] [-k=1.0] [-l=-1] [-n=1000000] [-p=0] [-Q=1000] [-r=0.25] [-t=298.15] [-w=10.0]
files...
Use the Transition-Tempered Orthogonal Space Random Walk algorithm to estimate a free energy.
    files...                        The atomic coordinate file in PDB or XYZ format.
    --af, --activeFinal=-1          Final active atom (single-topology only).
    --as, --activeStart=1           Starting active atom (single-topology only).
    --bM, --biasMag=0.05            OST Gaussian bias magnitude (kcal/mol).
    --dw, --distributeWalkers=OFF   AUTO: Pick up per-walker configurations as [filename. pdb]_[num], or specify a residue to distribute on.
    --ef1, --noElecFinal1=-1        Final no-electrostatics atom for 1st topology.
    --ef2, --noElecFinal2=-1        Final no-electrostatics atom for 2nd topology
    --es1, --noElecStart1=1         Starting no-electrostatics atom for 1st topology.
    --es2, --noElecStart2=1         Starting no-electrostatics atom for 2nd topology
    --f1, --final1=-1               Final ligand atom for the 1st topology.
    --f2, --final2=-1               Final ligand atom for the 2nd topology
    --la1, --ligAtoms1=-1           Period-separated ranges of 1st topology ligand atoms (e.g. 40-50.72-83).
    --la2, --ligAtoms2=-1           Period-separated ranges of 2nd toplogy ligand atoms (e.g. 40-50.72-83)
    --lf, --lambdaFriction=1.0E-18  Friction on the lambda particle.
    --lm, --lambdaMass=1.0E-18      Mass of the lambda particle.
    --lw, --lambdaWritOut=1.0       Only write out snapshots if lambda is greater than the value specified.
    --mc, monteCarlo                Monte Carlo OST
    --mcL, --mcLambdaStd=0.1        Standard deviation for lambda move.
    --mcMD, --mcTraj=100            Number of dynamics steps to take for each MD trajectory for Monte Carlo OST
    --np, --nParallel=1             Number of topologies to evaluate in parallel
    --rn, --resetNumSteps           Ignore prior steps logged in .lam or similar files
    --rsym, --randomSymOp=-1.0      Apply a random Cartesian symmetry operator with a random translation in the range -X .. X; less than 0 disables.
    --ruc, --randomUnitCell=-1.0    Apply random unit cell axes to achieve the specified density (g/cc).
    --s1, --start1=0                Starting ligand atom for 1st topology.
    --s2, --start2=0                Starting ligand atom for 2nd topology
    --sf, --switchingFunction=1.0   Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)
    --tp, --temperingParam=8.0      Dama et al tempering rate parameter in multiples of kBT
    --ts, --twoStep                 Sample MC-OST using separate lambda and MD moves.
    --uaA, --unsharedA=-1           Unshared atoms in the A dual topology (period-separated hyphenated ranges)
    --uaB, --unsharedB=-1           Unshared atoms in the B dual topology (period-separated hyphenated ranges)
-b, --thermostat=Bussi              Thermostat: [Adiabatic / Berendsen / Bussi].
-C, --count=10                      Time steps between MD-OST counts.
-d, --dt=1.0                        Time discretization step in femtoseconds.
-F, --fileFormat=XYZ                Choose file type to write [PDB/XYZ].
-h, --help                          Print this help message.
-i, --integrator=Verlet             Integrator: [Beeman / Respa / Stochastic / Verlet].
-k, --checkpoint=1.0                Interval to write out restart files (.dyn, .his, etc).
-l, --lambda=-1                     Initial lambda value.
-n, --numberOfSteps=1000000         Number of molecular dynamics steps.
-o, --optimize                      Optimize and save low-energy snapshots.
-p, --npt=0                         Specify use of a MC Barostat at the given pressure; the default 0 disables NPT (atm).
-Q, --equilibrate=1000              Number of equilibration steps before evaluation of thermodynamics.
-r, --report=0.25                   Interval to report thermodynamics (psec).
-t, --temperature=298.15            Temperature (Kelvin).
-w, --write=10.0                    Interval to write out coordinates (psec).
-y, --synchronous                   Walker communication is synchronous
 ```
 ---
 