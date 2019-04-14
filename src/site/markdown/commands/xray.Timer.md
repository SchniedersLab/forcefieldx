### X-ray Timer

Time calculation of the X-ray target.

---
```
Usage: ffxc xray.Timer [-AghSUv] [--rmo] [--sf] [--aRadBuffer=0.75] [--FSigFCutoff=-1.0] [--nBins=10] [--nResidueBFactor=0] [--nt=0] [--sigmaATol=0.05] [--sol=POLYNOMIAL] [--wA=1.0] [--xrayScaleTol=1.0e-4] [-B=1.0] [-G=0.6] [-m=coordinates] [-n=5] [-R=-1] [-X=<data> <data> <data>]... files...
Time calculation of the X-ray target.
    files...                        PDB and Diffraction input files.
    --aRadBuffer=0.75               Set the distance beyond the atomic radius to evaluate scattering (A).
    --FSigFCutoff=-1.0              F / SigF cutoff (-1.0 is no cutoff).
    --nBins=10                      The number of refection bins.
    --nResidueBFactor=0             Number of residues per B-factor. 0 uses atomic B-factors (default).
    --nt, --threads=0               Number of SMP threads (0 specifies use of all CPU cores).
    --rmo, --refineMolOcc           Refine on molecules.
    --sf, --splineFit               Use a resolution dependent spline scale.
    --sigmaATol=0.05                Sigma A optimization tolerance.
    --sol, --solvent=POLYNOMIAL     Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]
    --wA, --dataWeight=1.0          The weight of the real space data (wA).
    --xrayScaleTol=1.0e-4           X-ray scale optimization tolerance.
-A, --allGaussians                  Use all defined Gaussiansfor atomic scattering density (the default is to use the top 3).
-B, --bSimWeight=1.0                B-Factor similarity weight.
-g, --noGradient                    Ignore computation of the atomic coordinates noGradient.
-G, --sampling=0.6                  The number of grid spaces per Angstrom for the scattering FFT grid.
-h, --help                          Print this help message.
-m, --mode=coordinates              Refinement mode: coordinates, bfactors and/or occupancies.
-n, --iterations=5                  Number of iterations.
-R, --rFreeFlag=-1                  R-Free Flag value (-1 attempts to auto-determine from the data).
-S, --solventGridSearch             Perform a grid search for optimal bulk solvent parameters.
-U, --addAnisoU                     Add Anisotropic B-Factors to refinement.
-v, --verbose                       Print the energy for each iteration.
-X, --data=<data> <data> <data>     Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -X filename 1.0 false).
```
---
