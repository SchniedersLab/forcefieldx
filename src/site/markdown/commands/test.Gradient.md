### Test the Gradient

Test the potential energy gradient using finite-differences.

---
```
Usage: ffxc test.Gradient [-hv] [--la=-1] [--tol=1.0e-3] [-a=1] [-d=1.0e-5] files...
Test the potential energy gradient.
  files...                      The atomic coordinate file in PDB or XYZ format.
    --la, --lastAtomID=-1       The last atom to test (default is to test all Atoms, unless a first atom is specified).
    --tol, --tolerance=1.0e-3   The analytic vs. finite-difference gradient error tolerance.
-a, --atomID=1                  The first atom to test (default is Atom 1)
-d, --dx=1.0e-5                 The finite-difference step size.
-h, --help                      Print this help message.
-v, --verbose                   Print out the energy for each step.
 ```
 ---
 