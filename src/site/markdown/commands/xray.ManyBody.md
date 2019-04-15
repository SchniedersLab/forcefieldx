### X-ray Target Rotamer Optimization

Discrete optimization using a many-body expansion of an X-ray target function and elimination expressions.

---
```
Usage: ffxc xray.ManyBody [-AEhOSTUz] [--dee] [--out] [--rmo] [--sf] [--aRadBuffer=0.75] [--bB=0.0] [--bC=1] [--bL=20.0] [--ch=-1] [--eR=none] [--fi=-1] [--fR=-1,-1] [--FSigFCutoff=-1.0] [--increment=3] [--Ln=Richardson] [--lR=none] [--mC=-1] [--nB=3,3,3] [--nBins=10] [--nResidueBFactor=0] [--pr=1] [--radius=2.0] [--sigmaATol=0.05] [--sol=POLYNOMIAL] [--tC=3.0] [--thC=3.0] [--wA=1.0] [--window=7] [--xrayScaleTol=1.0e-4] [-a=0] [-B=1.0] [-G=0.6] [-L=2] [-m=coordinates] [-R=-1] [-s=-1] [-x=-1] [-X=<data> <data> <data>]... files...
Discrete optimization using a many-body expansion and elimination expressions.
    files...                            PDB and Real Space input files.
    --aRadBuffer=0.75                   Set the distance beyond the atomic radius to evaluate scattering (A).
    --bB, --boxBorderSize=0.0           Extent of overlap between optimization boxes in Angstroms.
    --bC, --boxInclusionCriterion=1     Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer
    --bL, --approxBoxLength=20.0        Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).
    --ch, --chain=-1                    Single character chain ID of the residues to optimize.
    --dee, --deadEnd                    Use dead-end elimination criteria instead of Goldstein criteria.
    --eR, --energyRestart=none          Load energy restart file from a previous run (requires that all parameters are the same).
    --fi, --final=-1                    Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.
    --fR, --forceResidues=-1,-1         Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.
    --FSigFCutoff=-1.0                  F / SigF cutoff (-1.0 is no cutoff).
    --increment=3                       Sliding window increment.
    --Ln, --libraryNucleic=Richardson   Nucleic acid library to select: [Richardson]
    --lR, --listResidues=none           Choose a list of individual residues to optimize (eg. A11,A24,B40).
    --mC, --monteCarlo=-1               Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.
    --nB, --numBoxes=3,3,3              Specify number of boxes along X, Y, and Z (default: 3,3,3)
    --nBins=10                          The number of refection bins.
    --nResidueBFactor=0                 Number of residues per B-factor. 0 uses atomic B-factors (default).
    --out, --output                     Save eliminated singles and eliminated pairs to a text file.
    --pr, --prune=1                     Prune no clashes (0), only single clashes (1), or all clashes (2)
    --radius=2.0                        The sliding window cutoff radius (Angstroms).
    --rmo, --refineMolOcc               Refine on molecules.
    --sf, --splineFit                   Use a resolution dependent spline scale.
    --sigmaATol=0.05                    Sigma A optimization tolerance.
    --sol, --solvent=POLYNOMIAL         Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]
    --tC, --twoBodyCutoff=3.0           Cutoff distance for two body interactions.
    --thC, --threeBodyCutoff=3.0        Cutoff distance for three-body interactions.
    --wA, --dataWeight=1.0              The weight of the real space data (wA).
    --window=7                          Size of the sliding window with respect to adjacent residues.
    --xrayScaleTol=1.0e-4               X-ray scale optimization tolerance.
-a, --algorithm=0                       Algorithm: default automatic settings (0), independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)
-A, --allGaussians                      Use all defined Gaussiansfor atomic scattering density (the default is to use the top 3).
-B, --bSimWeight=1.0                    B-Factor similarity weight.
-E, --decompose                         Print energy decomposition for the input structure (no optimization!).
-G, --sampling=0.6                      The number of grid spaces per Angstrom for the scattering FFT grid.
-h, --help                              Print this help message.
-L, --library=2                         Ponder and Richards (1) or Richardson (2) rotamer library.
-m, --mode=coordinates                  Refinement mode: coordinates, bfactors and/or occupancies.
-O, --noOriginal                        Do not include starting coordinates as their own rotamer.
-R, --rFreeFlag=-1                      R-Free Flag value (-1 attempts to auto-determine from the data).
-s, --start=-1                          Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.
-S, --solventGridSearch                 Perform a grid search for optimal bulk solvent parameters.
-T, --threeBody                         Include 3-Body interactions in the elimination criteria.
-U, --addAnisoU                         Add Anisotropic B-Factors to refinement.
-x, --all=-1                            Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.
-X, --data=<data> <data> <data>         Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -X filename 1.0 false).
-z, --revert                            Revert unfavorable changes.
```
---
