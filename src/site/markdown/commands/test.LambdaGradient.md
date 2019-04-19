### Test Lambda Derivatives

Test the potential energy gradient using finite-differences.

---
```
Usage: ffxc test.LambdaGradient [-hv] [--ls] [--sdX] [--sk2] [--af=-1] [--as=1] [--ef1=-1] [--ef2=-1] [--es1=1] [--es2=1] [--f1=-1] [--f2=-1] [--la=-1] [--la1=-1] [--la2=-1] [--lm=0.01] [--np=1] [--s1=0] [--s2=0] [--sf=1.0] [--tol=1.0e-3] [--uaA=-1] [--uaB=-1] [-a=1] [-d=1.0e-5] [-l=-1] files...
Test potential energy derivatives with respect to Lambda.
    files...                        The atomic coordinate file in PDB or XYZ format.
    --af, --activeFinal=-1          Final active atom (single-topology only).
    --as, --activeStart=1           Starting active atom (single-topology only).
    --ef1, --noElecFinal1=-1        Final no-electrostatics atom for 1st topology.
    --ef2, --noElecFinal2=-1        Final no-electrostatics atom for 2nd topology
    --es1, --noElecStart1=1         Starting no-electrostatics atom for 1st topology.
    --es2, --noElecStart2=1         Starting no-electrostatics atom for 2nd topology
    --f1, --final1=-1               Final ligand atom for the 1st topology.
    --f2, --final2=-1               Final ligand atom for the 2nd topology
    --la, --lastAtomID=-1           The last atom to test (default is to test all Atoms, unless a first atom is specified).
    --la1, --ligAtoms1=-1           Period-separated ranges of 1st topology ligand atoms (e.g. 40-50.72-83).
    --la2, --ligAtoms2=-1           Period-separated ranges of 2nd toplogy ligand atoms (e.g. 40-50.72-83)
    --lm, --lambdaMoveSize=0.01     Size of the lambda moves during the test.
    --ls, --lambdaScan              Scan lambda values.
    --np, --nParallel=1             Number of topologies to evaluate in parallel
    --s1, --start1=0                Starting ligand atom for 1st topology.
    --s2, --start2=0                Starting ligand atom for 2nd topology
    --sdX, --skipdX                 Skip calculating per-atom dUdX values and only test lambda gradients.
    --sf, --switchingFunction=1.0   Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)
    --sk2, --skip2                  Skip 2nd derivatives.
    --tol, --tolerance=1.0e-3       The analytic vs. finite-difference gradient error tolerance.
    --uaA, --unsharedA=-1           Unshared atoms in the A dual topology (period-separated hyphenated ranges)
    --uaB, --unsharedB=-1           Unshared atoms in the B dual topology (period-separated hyphenated ranges)
-a, --atomID=1                      The first atom to test (default is Atom 1)
-d, --dx=1.0e-5                     The finite-difference step size.
-h, --help                          Print this help message.
-l, --lambda=-1                     Initial lambda value.
-v, --verbose                       Print out the energy for each step.
```
---
 