### X-ray Refinement

Refine an X-ray, Neutron or Joint X-ray/Neutron target.

---
```
Usage: ffxc xray.Minimize [-AhStU] [--rmo] [--sf] [--aRadBuffer=0.75] [--FSigFCutoff=-1.0] [--nBins=10] [--nResidueBFactor=0] [--sigmaATol=0.05] [--sol=POLYNOMIAL] [--suffix=_refine] [--wA=1.0] [--xrayScaleTol=1.0e-4] [-B=1.0] [-e=1.0] [-G=0.6] [-I=Unlimited] [-m=coordinates] [-R=-1] [-E=-1.0 -1.0 -1.0]... [-X=<data> <data> <data>]... files...
Refine an X-ray/Neutron target.
  files...                          PDB and Diffraction input files.
  --aRadBuffer=0.75                 Set the distance beyond the atomic radius to evaluate scattering (A).
  --FSigFCutoff=-1.0                F / SigF cutoff (-1.0 is no cutoff).
  --nBins=10                        The number of refection bins.
  --nResidueBFactor=0               Number of residues per B-factor. 0 uses atomic B-factors (default).
  --rmo, --refineMolOcc             Refine on molecules.
  --sf, --splineFit                 Use a resolution dependent spline scale.
  --sigmaATol=0.05                  Sigma A optimization tolerance.
  --sol, --solvent=POLYNOMIAL       Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]
  --suffix=_refine                  Suffix to apply to files written out by minimization.
  --wA, --dataWeight=1.0            The weight of the real space data (wA).
  --xrayScaleTol=1.0e-4             X-ray scale optimization tolerance.
-A, --allGaussians                  Use all defined Gaussiansfor atomic scattering density (the default is to use the top 3).
-B, --bSimWeight=1.0                B-Factor similarity weight.
-e, --eps=1.0                       Convergence criteria.
-E, --eps3=-1.0 -1.0 -1.0           RMS gradient convergence criteria for three stage refinement (default of -1.0 automatically determines eps for each stage).
-G, --sampling=0.6                  The number of grid spaces per Angstrom for the scattering FFT grid.
-h, --help                          Print this help message.
-I, --iterations=Unlimited          Number of minimization steps.
-m, --mode=coordinates              Refinement mode: coordinates, bfactors and/or occupancies.
-R, --rFreeFlag=-1                  R-Free Flag value (-1 attempts to auto-determine from the data).
-S, --solventGridSearch             Perform a grid search for optimal bulk solvent parameters.
-t, --threeStage                    Perform refinement in 3 stages: coordinates, b-factors, and then occupancies (overrides mode setting if true)
-U, --addAnisoU                     Add Anisotropic B-Factors to refinement.
-X, --data=<data> <data> <data>     Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -X filename 1.0 false).
```
---
