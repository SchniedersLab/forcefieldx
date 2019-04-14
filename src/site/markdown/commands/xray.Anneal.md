### X-ray Target Simulated Annealing

Simulated annealing on a X-ray refinement target.

---
```
Usage: ffxc xray.Anneal [-AhoSU] [--rmo] [--sf] [--aRadBuffer=0.75] [--FSigFCutoff=-1.0] [--nBins=10] [--nResidueBFactor=0] [--sigmaATol=0.05] [--sol=POLYNOMIAL] [--tl=10.0] [--tu=1000.0] [--wA=1.0] [--xrayScaleTol=1.0e-4] [-b=Bussi] [-B=1.0] [-d=1.0] [-G=0.6] [-i=Verlet] [-k=1.0] [-m=coordinates] [-n=1000000] [-r=0.25] [-R=-1] [-t=298.15] [-w=10.0] [-W=10] [-X=<data> <data> <data>]... files...
Simulated annealing on an X-ray target.
    files...                        PDB and Diffraction input files.
    --aRadBuffer=0.75               Set the distance beyond the atomic radius to evaluate scattering (A).
    --FSigFCutoff=-1.0              F / SigF cutoff (-1.0 is no cutoff).
    --nBins=10                      The number of refection bins.
    --nResidueBFactor=0             Number of residues per B-factor. 0 uses atomic B-factors (default).
    --rmo, --refineMolOcc           Refine on molecules.
    --sf, --splineFit               Use a resolution dependent spline scale.
    --sigmaATol=0.05                Sigma A optimization tolerance.
    --sol, --solvent=POLYNOMIAL     Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]
    --tl, --temperatureLow=10.0     Low temperature limit (Kelvin).
    --tu, --temperatureUpper=1000.0 High temperature limit (Kelvin).
    --wA, --dataWeight=1.0          The weight of the real space data (wA).
    --xrayScaleTol=1.0e-4           X-ray scale optimization tolerance.
-A, --allGaussians                  Use all defined Gaussiansfor atomic scattering density (the default is to use the top 3).
-b, --thermostat=Bussi              Thermostat: [Adiabatic / Berendsen / Bussi].
-B, --bSimWeight=1.0                B-Factor similarity weight.
-d, --dt=1.0                        Time discretization step in femtoseconds.
-G, --sampling=0.6                  The number of grid spaces per Angstrom for the scattering FFT grid.
-h, --help                          Print this help message.
-i, --integrator=Verlet             Integrator: [Beeman / Respa / Stochastic / Verlet].
-k, --checkpoint=1.0                Interval to write out restart files (.dyn, .his, etc).
-m, --mode=coordinates              Refinement mode: coordinates, bfactors and/or occupancies.
-n, --numberOfSteps=1000000         Number of molecular dynamics steps.
-o, --optimize                      Optimize and save low-energy snapshots.
-r, --report=0.25                   Interval to report thermodynamics (psec).
-R, --rFreeFlag=-1                  R-Free Flag value (-1 attempts to auto-determine from the data).
-S, --solventGridSearch             Perform a grid search for optimal bulk solvent parameters.
-t, --temperature=298.15            Temperature (Kelvin).
-U, --addAnisoU                     Add Anisotropic B-Factors to refinement.
-w, --write=10.0                    Interval to write out coordinates (psec).
-W, --windows=10                    Number of annealing windows.
-X, --data=<data> <data> <data>     Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -X filename 1.0 false).
```
---
