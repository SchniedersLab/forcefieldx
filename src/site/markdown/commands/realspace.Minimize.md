### Real Space Refinement

Minimize a Real Space target.

---
```
Usage: ffxc realspace.Minimize [-h] [--wA=1.0] [-e=1.0] [-I=Unlimited] [-X=<data> <data>]... files...
Minimization on a Real Space target.
    files...                PDB and Real Space input files.
    --wA, --dataWeight=1.0  The weight of the real space data (wA).
-e, --eps=1.0               Convergence criteria.
-h, --help                  Print this help message.
-I, --iterations=Unlimited  Number of minimization steps.
-X, --data=<data> <data>    Specify input data filename, and its weight (wA) (e.g. -X filename 1.0).
```
---
 