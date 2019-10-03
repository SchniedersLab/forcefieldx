package ffx.algorithms.optimize.manybody;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import edu.rit.pj.WorkerTeam;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.utils.EnergyException;

public class EnergyExpansion {

    private static final Logger logger = Logger.getLogger(EnergyExpansion.class.getName());

    /**
     * The potential energy of the system with all side-chains to be optimized
     * turned off.
     */
    private double backboneEnergy;
    /**
     * Self-energy of each residue for each rotamer. [residue][rotamer]
     */
    private double[][] selfEnergy;
    /**
     * Two-body energies for each pair of residues and pair of rotamers.
     * [residue1][rotamer1][residue2][rotamer2]
     */
    private double[][][][] twoBodyEnergy;
    /**
     * Trimer-energies for each trimer of rotamers.
     * [residue1][rotamer1][residue2][rotamer2][residue3][rotamer3]
     */
    private double[][][][][][] threeBodyEnergy;
    /**
     * Map of self-energy values to compute.
     */
    private final HashMap<Integer, Integer[]> selfEnergyMap = new HashMap<>();
    /**
     * Map of 2-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> twoBodyEnergyMap = new HashMap<>();
    /**
     * Map of 3-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> threeBodyEnergyMap = new HashMap<>();
    /**
     * Map of 4-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> fourBodyEnergyMap = new HashMap<>();

    private RotamerOptimization rO;
    private DistanceMatrix dM;
    private EliminatedRotamers eR;
    /**
     * MolecularAssembly to perform rotamer optimization on.
     */
    private MolecularAssembly molecularAssembly;
    /**
     * RotamerLibrary instance.
     */
    private RotamerLibrary library;
    /**
     * AlgorithmListener who should receive updates as the optimization runs.
     */
    private AlgorithmListener algorithmListener;
    /**
     * A list of all residues being optimized. Note that Box and Window
     * optimizations operate on subsets of this list.
     */
    private ArrayList<Residue> allResiduesList;
    /**
     * Interaction partners of a Residue that come after it.
     */
    private int[][] resNeighbors;
    /**
     * Flag to control use of 3-body terms.
     */
    private boolean threeBodyTerm;
    /**
     * Flag to indicate computation of a many-body expansion for original
     * coordinates.
     */
    private boolean decomposeOriginal;
    /**
     * Flag to indicate use of box optimization.
     */
    private boolean usingBoxOptimization;
    /**
     * Flag to indicate verbose logging.
     */
    private boolean verbose;
    /**
     * Flag to prune clashes.
     */
    private boolean pruneClashes;
    /**
     * Flag to prune pair clashes.
     */
    private boolean prunePairClashes;
    /**
     * Flag to indicate if this is the master process.
     */
    private final boolean master;

    /**
     * Maximum number of 4-body energy values to compute.
     */
    private int max4BodyCount;

    /**
     * Default value for the ommRecalculateThreshold in kcal/mol.
     */
    private static final double DEFAULT_OMM_RECALCULATE_THRESHOLD = -200;
    /**
     * If OpenMM is in use, recalculate any suspiciously low self/pair/triple
     * energies using the reference Java implementation.
     */
    private final double ommRecalculateThreshold;
    /**
     * Default value for the singularityThreshold in kcal/mol.
     */
    private static final double DEFAULT_SINGULARITY_THRESHOLD = -1000;
    /**
     * Reject any self, pair, or triple energy beneath singularityThreshold.
     * <p>
     * This is meant to check for running into singularities in the energy
     * function, where an energy is "too good to be true" and actually
     * represents a poor conformation.
     */
    private final double singularityThreshold;
    /**
     * Indicates if the Potential is an OpenMMForceFieldEnergy.
     */
    private final boolean potentialIsOpenMM;

    public EnergyExpansion(RotamerOptimization rO, DistanceMatrix dM, EliminatedRotamers eR,
                           MolecularAssembly molecularAssembly, Potential potential,
                           RotamerLibrary library, AlgorithmListener algorithmListener, ArrayList<Residue> allResiduesList,
                           int[][] resNeighbors, boolean threeBodyTerm, boolean decomposeOriginal,
                           boolean usingBoxOptimization, boolean verbose,
                           boolean pruneClashes, boolean prunePairClashes, boolean master) {
        this.rO = rO;
        this.dM = dM;
        this.eR = eR;
        this.molecularAssembly = molecularAssembly;
        this.library = library;
        this.algorithmListener = algorithmListener;
        this.allResiduesList = allResiduesList;
        this.resNeighbors = resNeighbors;
        this.threeBodyTerm = threeBodyTerm;
        this.decomposeOriginal = decomposeOriginal;
        this.usingBoxOptimization = usingBoxOptimization;
        this.verbose = verbose;
        this.pruneClashes = pruneClashes;
        this.prunePairClashes = prunePairClashes;
        this.master = master;

        CompositeConfiguration properties = molecularAssembly.getProperties();
        max4BodyCount = properties.getInt("ro-max4BodyCount", Integer.MAX_VALUE);
        if (max4BodyCount != Integer.MAX_VALUE) {
            logger.info(format(" Max 4Body Count: %d", max4BodyCount));
        }
        singularityThreshold = properties.getDouble("ro-singularityThreshold", DEFAULT_SINGULARITY_THRESHOLD);
        potentialIsOpenMM = potential instanceof ForceFieldEnergyOpenMM;
        if (potentialIsOpenMM) {
            ommRecalculateThreshold = properties.getDouble("ro-ommRecalculateThreshold", DEFAULT_OMM_RECALCULATE_THRESHOLD);
        } else {
            ommRecalculateThreshold = -1E200;
        }
    }

    public void setBackboneEnergy(double backboneEnergy) {
        this.backboneEnergy = backboneEnergy;
    }

    /**
     * Computes a self energy, defined as energy with all sidechains but one
     * turned off, minus the backbone energy.
     *
     * @param residues Residues under optimization.
     * @param i        A residue index.
     * @param ri       A rotamer index for residue i.
     * @return Eself(ri)=E1(ri)-Eenv/bb.
     */
    public double computeSelfEnergy(Residue[] residues, int i, int ri) {
        rO.turnOffAllResidues(residues);
        rO.turnOnResidue(residues[i], ri);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            energy = rO.currentEnergy(residues) - backboneEnergy;
            if (potentialIsOpenMM && energy < ommRecalculateThreshold) {
                logger.warning(format(" Experimental: re-computing self energy %s-%d using Force Field X", residues[i], ri));
                energy = rO.currentFFXPE() - backboneEnergy;
            }
            if (energy < singularityThreshold) {
                String message = format(" Rejecting self energy for %s-%d is %10.5g << %10f, "
                        + "likely in error.", residues[i], ri, energy, singularityThreshold);
                logger.warning(message);
                throw new EnergyException(message);
            }
        } finally {
            rO.turnOffResidue(residues[i]);
        }

        return energy;
    }

    /**
     * Computes a pair energy, defined as energy with all sidechains but two
     * turned off, minus the sum of backbone and component self energies.
     *
     * @param residues Residues under optimization.
     * @param i        A residue index.
     * @param ri       A rotamer index for residue i.
     * @param j        A residue index j!=i.
     * @param rj       A rotamer index for residue j.
     * @return Epair(ri, rj)=E2(ri,rj)-Eself(ri)-Eself(rj)-Eenv/bb.
     */
    public double compute2BodyEnergy(Residue[] residues, int i, int ri, int j, int rj) {
        rO.turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            double subtract = -backboneEnergy - getSelf(i, ri) - getSelf(j, rj);
            energy = rO.currentEnergy(residues) + subtract;
            if (potentialIsOpenMM && energy < ommRecalculateThreshold) {
                logger.warning(format(" Experimental: re-computing pair energy %s-%d %s-%d using Force Field X", residues[i], ri, residues[j], rj));
                energy = rO.currentFFXPE() + subtract;
            }
            if (energy < singularityThreshold) {
                String message = format(" Rejecting pair energy for %s-%d %s-%d is %10.5g << %10f, "
                        + "likely in error.", residues[i], ri, residues[j], rj, energy, singularityThreshold);
                logger.warning(message);
                throw new EnergyException(message);
            }
        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
        }
        return energy;
    }

    /**
     * Computes a 3-body energy, defined as the energy with all sidechains but
     * three turned off, minus the sum of backbone and component self/2-Body
     * energies.
     *
     * @param residues Residues under optimization.
     * @param i        A residue index.
     * @param ri       A rotamer index for residue i.
     * @param j        A residue index j!=i.
     * @param rj       A rotamer index for residue j.
     * @param k        A residue index k!=j k!=i.
     * @param rk       A rotamer index for residue k.
     * @return Etri(ri,
     *rj)=E3(ri,rj,rk)-Epair(ri,rj)-Epair(ri,rk)-Epair(rj,rk)-Eself(ri)-Eself(rj)-Eself(rk)-Eenv/bb.
     */
    public double compute3BodyEnergy(Residue[] residues, int i, int ri, int j, int rj, int k, int rk) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        turnOnResidue(residues[k], rk);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            double subtract = -backboneEnergy - getSelf(i, ri) - getSelf(j, rj) - getSelf(k, rk)
                    - get2Body(i, ri, j, rj) - get2Body(i, ri, k, rk) - get2Body(j, rj, k, rk);
            energy = rO.currentEnergy(residues) + subtract;
            if (potentialIsOpenMM && energy < ommRecalculateThreshold) {
                logger.warning(format(" Experimental: re-computing triple energy %s-%d %s-%d %s-%d using Force Field X", residues[i], ri, residues[j], rj, residues[k], rk));
                energy = rO.currentFFXPE() + subtract;
            }
            if (energy < singularityThreshold) {
                String message = format(" Rejecting triple energy for %s-%d %s-%d %s-%d is %10.5g << %10f, "
                        + "likely in error.", residues[i], ri, residues[j], rj, residues[k], rk, energy, singularityThreshold);
                logger.warning(message);
                throw new EnergyException(message);
            }
        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
            turnOffResidue(residues[k]);
        }
        return energy;
    }

    /**
     * Computes a 4-body energy, defined as the energy with all sidechains but
     * four turned off, minus the sum of backbone and component
     * self/2-Body/3-body energies.
     *
     * @param residues Residues under optimization.
     * @param i        A residue index.
     * @param ri       A rotamer index for residue i.
     * @param j        A residue index j!=i.
     * @param rj       A rotamer index for residue j.
     * @param k        A residue index k!=j k!=i.
     * @param rk       A rotamer index for residue k.
     * @param l        A residue index l!=i l!=j l!=k.
     * @param rl       A rotamer index for residue l.
     * @return The 4-body energy.
     */
    public double compute4BodyEnergy(Residue[] residues, int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        turnOnResidue(residues[k], rk);
        turnOnResidue(residues[l], rl);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            double subtract = -backboneEnergy
                    - getSelf(i, ri) - getSelf(j, rj) - getSelf(k, rk) - getSelf(l, rl)
                    - get2Body(i, ri, j, rj) - get2Body(i, ri, k, rk) - get2Body(i, ri, l, rl)
                    - get2Body(j, rj, k, rk) - get2Body(j, rj, l, rl) - get2Body(k, rk, l, rl)
                    - get3Body(residues, i, ri, j, rj, k, rk) - get3Body(residues, i, ri, j, rj, l, rl) - get3Body(residues, i, ri, k, rk, l, rl) - get3Body(residues, j, rj, k, rk, l, rl);
            energy = rO.currentEnergy(residues) + subtract;

            if (potentialIsOpenMM && energy < ommRecalculateThreshold) {
                logger.warning(format(" Experimental: re-computing quad energy %s-%d %s-%d %s-%d %s-%d using Force Field X", residues[i], ri, residues[j], rj, residues[k], rk, residues[l], rl));
                energy = rO.currentFFXPE() + subtract;
            }
            if (energy < singularityThreshold) {
                String message = format(" Rejecting quad energy for %s-%d %s-%d %s-%d %s-%d is %10.5g << %10f, "
                        + "likely in error.", residues[i], ri, residues[j], rj, residues[k], rk, residues[l], rl, energy, singularityThreshold);
                logger.warning(message);
                throw new EnergyException(message);
            }

        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
            turnOffResidue(residues[k]);
            turnOffResidue(residues[l]);
        }
        return energy;
    }

    public void setSelf(int i, int ri, double e) {
        setSelf(i, ri, e, false);
    }

    /**
     * Stores a self energy in the self energy matrix.
     *
     * @param i     A residue index.
     * @param ri    A rotamer for residue i.
     * @param e     Computed energy to store.
     * @param quiet Silence warnings about exceptions.
     */
    public void setSelf(int i, int ri, double e, boolean quiet) {
        try {
            selfEnergy[i][ri] = e;
        } catch (NullPointerException | ArrayIndexOutOfBoundsException ex) {
            if (!quiet) {
                logger.warning(format(" NPE or array index error for (%3d,%2d)", i, ri));
            }
            throw ex;
        }
    }

    public void set2Body(int i, int ri, int j, int rj, double e) {
        set2Body(i, ri, j, rj, e, false);
    }

    /**
     * Stores a pair energy in the pairs energy matrix.
     *
     * @param i     A residue index.
     * @param ri    A rotamer for residue i.
     * @param j     A residue index j != i.
     * @param rj    A rotamer for residue j.
     * @param e     Computed energy to store.
     * @param quiet Silence warnings about exceptions.
     */
    public void set2Body(int i, int ri, int j, int rj, double e, boolean quiet) {
        // Ensure i < j.
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        try {
            // Find where j is in i's neighbor list (and thus the 2-body energy matrix).
            int[] nI = resNeighbors[i];
            int indJ = -1;
            for (int l = 0; l < nI.length; l++) {
                if (nI[l] == j) {
                    indJ = l;
                    break;
                }
            }
            if (indJ == -1) {
                throw new IllegalArgumentException(format(" Residue %d not found in neighbors of %d; assumed past cutoff.", j, i));
            } else {
                twoBodyEnergy[i][ri][indJ][rj] = e;
            }
        } catch (NullPointerException npe) {
            if (!quiet) {
                logger.info(format(" NPE for 2-body energy (%3d,%2d) (%3d,%2d).", i, ri, j, rj));
            }
            throw npe;
        }

    }

    /**
     * <p>
     * set3Body.</p>
     *
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     * @param i        a int.
     * @param ri       a int.
     * @param j        a int.
     * @param rj       a int.
     * @param k        a int.
     * @param rk       a int.
     * @param e        a double.
     */
    public void set3Body(Residue[] residues, int i, int ri, int j, int rj, int k, int rk, double e) {
        set3Body(residues, i, ri, j, rj, k, rk, e, false);
    }

    /**
     * Stores a triple energy in the triples energy matrix.
     *
     * @param i        A residue index.
     * @param ri       A rotamer for residue i.
     * @param j        A residue index j != i.
     * @param rj       A rotamer for residue j.
     * @param k        A residue index k != j, k != i.
     * @param rk       A rotamer for residue k.
     * @param e        Computed energy to store.
     * @param quiet    Silence warnings about exceptions.
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     * @throws java.lang.IllegalStateException If threeBodyTerm is false.
     */
    public void set3Body(Residue[] residues, int i, int ri, int j, int rj, int k, int rk, double e, boolean quiet)
            throws IllegalStateException {
        if (!threeBodyTerm) {
            throw new IllegalStateException(" Attempting to set a 3-body energy when threeBodyTerm is false!");
        }
        // Ensure i < j and j < k.
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        if (k < i) {
            int ii = i;
            int iri = ri;
            i = k;
            ri = rk;
            k = ii;
            rk = iri;
        }
        if (k < j) {
            int jj = j;
            int jrj = rj;
            j = k;
            rj = rk;
            k = jj;
            rk = jrj;
        }

        // Find where j is in i's neighbor list, and where k is in j's neighbor list.
        int[] nI = resNeighbors[i];
        int indJ = -1;
        for (int l = 0; l < nI.length; l++) {
            if (nI[l] == j) {
                indJ = l;
                break;
            }
        }

        int[] nJ = resNeighbors[j];
        int indK = -1;
        for (int l = 0; l < nJ.length; l++) {
            if (nJ[l] == k) {
                indK = l;
                break;
            }
        }

        // i,j,k: Indices in the current Residue array.
        // indJ, indK: Index of j in i's neighbor list, index of k in j's neighbor list.
        // indexI, indexJ, indexK: Indices in allResiduesList.
        int indexI = allResiduesList.indexOf(residues[i]);
        int indexJ = allResiduesList.indexOf(residues[j]);
        int indexK = allResiduesList.indexOf(residues[k]);
        if (dM.checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
            throw new IllegalArgumentException(format(" Residue %d not found in neighbors of %d; assumed past cutoff.", j, i));
        } else {
            try {
                threeBodyEnergy[i][ri][indJ][rj][indK][rk] = e;
            } catch (NullPointerException | ArrayIndexOutOfBoundsException ex) {
                if (!quiet) {
                    String message = format(" Could not access threeBodyEnergy array for 3-body energy (%3d,%2d) (%3d,%2d) (%3d,%2d)", i, ri, j, rj, k, rk);
                    logger.warning(message);
                }
                throw ex;
            }
        }
    }

    public double getBackboneEnergy() {
        return backboneEnergy;
    }

    /**
     * Return a previously computed self-energy.
     *
     * @param i  Residue i.
     * @param ri Rotamer ri of residue i.
     * @return The self-energy.
     */
    public double getSelf(int i, int ri) {
        try {
            return selfEnergy[i][ri];
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for self energy (%3d,%2d).", i, ri));
            throw npe;
        }
    }

    /**
     * Return a previously computed 2-body energy.
     *
     * @param i  Residue i.
     * @param ri Rotamer ri of residue i.
     * @param j  Residue j.
     * @param rj Rotamer rj of residue j.
     * @return The 2-Body energy.
     */
    public double get2Body(int i, int ri, int j, int rj) {
        // Ensure i < j.
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        try {
            // Find where j is in i's neighbor list (and thus the 2-body energy matrix).
            int[] nI = resNeighbors[i];
            int indJ = -1;
            for (int l = 0; l < nI.length; l++) {
                if (nI[l] == j) {
                    indJ = l;
                    break;
                }
            }
            if (indJ == -1) {
                logger.fine(format(" Residue %d not found in neighbors of %d; assumed past cutoff.", j, i));
                return 0;
            } else {
                return twoBodyEnergy[i][ri][indJ][rj];
            }
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for 2-body energy (%3d,%2d) (%3d,%2d).", i, ri, j, rj));
            throw npe;
        }
    }

    /**
     * Return a previously computed 3-body energy.
     *
     * @param i        Residue i.
     * @param ri       Rotamer ri of residue i.
     * @param j        Residue j.
     * @param rj       Rotamer rj of residue j.
     * @param k        Residue k.
     * @param rk       Rotamer rk of residue k.
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     * @return The 3-Body energy.
     */
    public double get3Body(Residue[] residues, int i, int ri, int j, int rj, int k, int rk) {
        if (!threeBodyTerm) {
            return 0.0;
        }
        // Ensure i < j and j < k.
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        if (k < i) {
            int ii = i;
            int iri = ri;
            i = k;
            ri = rk;
            k = ii;
            rk = iri;
        }
        if (k < j) {
            int jj = j;
            int jrj = rj;
            j = k;
            rj = rk;
            k = jj;
            rk = jrj;
        }

        // Find where j is in i's neighbor list, and where k is in j's neighbor list.
        int[] nI = resNeighbors[i];
        int indJ = -1;
        for (int l = 0; l < nI.length; l++) {
            if (nI[l] == j) {
                indJ = l;
                break;
            }
        }

        int[] nJ = resNeighbors[j];
        int indK = -1;
        for (int l = 0; l < nJ.length; l++) {
            if (nJ[l] == k) {
                indK = l;
                break;
            }
        }

        // i,j,k: Indices in the current Residue array.
        // indJ, indK: Index of j in i's neighbor list, index of k in j's neighbor list.
        // indexI, indexJ, indexK: Indices in allResiduesList.
        int indexI = allResiduesList.indexOf(residues[i]);
        int indexJ = allResiduesList.indexOf(residues[j]);
        int indexK = allResiduesList.indexOf(residues[k]);
        if (dM.checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
            return 0;
        } else {
            try {
                return threeBodyEnergy[i][ri][indJ][rj][indK][rk];
            } catch (NullPointerException | ArrayIndexOutOfBoundsException ex) {
                String message = format(" Could not find an energy for 3-body energy (%3d,%2d) (%3d,%2d) (%3d,%2d)", i, ri, j, rj, k, rk);
                logger.warning(message);
                throw ex;
            }
        }
    }

    /**
     * Return the lowest self-energy for residue i.
     *
     * @param residues
     * @param i
     * @return
     */
    public double lowestSelfEnergy(Residue[] residues, int i) {
        if (residues == null) {
            return 0.0;
        }
        int n = residues.length;
        if (i < 0 || i >= n) {
            return 0.0;
        }
        Rotamer[] rotamers = residues[i].getRotamers(library);
        int nr = rotamers.length;
        double energy = Double.MAX_VALUE;
        for (int ni = 0; ni < nr; ni++) {
            try {
                double e = getSelf(i, ni);
                if (e < energy) {
                    energy = e;
                }
            } catch (Exception e) {
                continue;
            }

        }
        return energy;
    }

    /**
     * Return the lowest pair-energy for residue (i,ri) with residue j.
     *
     * @param residues
     * @param i
     * @param ri
     * @param j
     * @return
     */
    public double lowestPairEnergy(Residue[] residues, int i, int ri, int j) {
        if (residues == null) {
            return 0.0;
        }
        int n = residues.length;
        if (i < 0 || i >= n) {
            return 0.0;
        }
        if (j < 0 || j >= n) {
            return 0.0;
        }

        Rotamer[] rotamers = residues[j].getRotamers(library);
        int nr = rotamers.length;
        double energy = Double.MAX_VALUE;
        for (int jr = 0; jr < nr; jr++) {
            try {
                double e = get2Body(i, ri, j, jr);
                if (e < energy) {
                    energy = e;
                }
            } catch (Exception e) {
                continue;
            }
        }
        return energy;
    }

    /**
     * Find the min/max of the 2-body energy.
     *
     * @param residues The residue array.
     * @param minMax   The bound on the 3-body energy (minMax[0] = min, minMax[1]
     *                 = max.
     * @param i        Residue i
     * @param ri       Rotamer ri of Residue i
     * @param j        Residue j
     * @param rj       Rotamer rj of Residue j
     * @return true if this term is valid.
     */
    private boolean minMax2BodySum(Residue[] residues, double[] minMax, int i, int ri, int j, int rj) {
        int nres = residues.length;
        double minSum = 0.0;
        double maxSum = 0.0;
        for (int k = 0; k < nres; k++) {
            if (k == i || k == j) {
                continue;
            }
            Residue residuek = residues[k];
            Rotamer[] romatersk = residuek.getRotamers(library);
            int lenrk = romatersk.length;
            double currentMin = Double.MAX_VALUE;
            double currentMax = Double.MIN_VALUE;
            for (int rk = 0; rk < lenrk; rk++) {
                if (eR.check(k, rk)) {
                    // k,rk is part of no valid phase space, so ignore it.
                    continue;
                }
                if (eR.check(i, ri, k, rk) || eR.check(j, rj, k, rk)) {
                    // Not implemented: check(i, ri, j, rj, k, rk).
                    // k,rk conflicts with i,ri or j,rj, so the max is now Double.NaN. No effect on minimum.
                    currentMax = Double.NaN;
                } else {
                    double current = get3Body(residues, i, ri, j, rj, k, rk);
                    if (Double.isFinite(current) && current < currentMin) {
                        currentMin = current;
                    } // Else, no new minimum found.
                    if (Double.isFinite(current) && Double.isFinite(currentMax)) {
                        if (current > currentMax) {
                            currentMax = current;
                        } // Else, we have failed to find a new finite maximum.
                    } else {
                        // The maximum is NaN.
                        currentMax = Double.NaN;
                    }
                }
            }
            if (currentMin == Double.MAX_VALUE || !Double.isFinite(minSum)) {
                // We have failed to find a viable configuration for i,ri,j,rj, as it conflicts with all rk for this k.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            } else {
                // Add finite current min to minSum.
                minSum += currentMin;
            }
            if (Double.isFinite(maxSum) && Double.isFinite(currentMax)) {
                maxSum += currentMax;
            } else {
                maxSum = Double.NaN;
            }
        }
        minMax[0] = minSum;
        minMax[1] = maxSum;
        return true;
    }

    /**
     * Computes the maximum and minimum energy i,ri might have with j, and
     * optionally (if three-body energies in use) third residues k.
     * <p>
     * The return value should be redundant with minMax[0] being NaN.
     *
     * @param residues Array of residues under consideration.
     * @param minMax   Index 0 to be filled by minimum energy, index 1 filled by
     *                 maximum energy.
     * @param i        Some residue i under consideration.
     * @param ri       A rotamer for residue i.
     * @param j        Some arbitrary residue i!=j.
     * @return If a valid configuration between i,ri and j could be found.
     */
    public boolean minMaxPairEnergy(Residue[] residues, double[] minMax, int i, int ri, int j) {
        Residue residuej = residues[j];
        Rotamer[] rotamersj = residuej.getRotamers(library);
        int lenrj = rotamersj.length;
        boolean valid = false;
        minMax[0] = Double.MAX_VALUE;
        minMax[1] = Double.MIN_VALUE;

        // Loop over the 2nd residues' rotamers.
        for (int rj = 0; rj < lenrj; rj++) {
            // Check for an eliminated single or eliminated pair.
            if (eR.check(i, ri) || eR.check(j, rj) || eR.check(i, ri, j, rj)) {
                continue;
            }

            double currMax = get2Body(i, ri, j, rj);
            double currMin = currMax; // Will remain identical if truncating at 2-body.

            if (threeBodyTerm) {
                double[] minMaxTriple = new double[2];
                // Loop over residue k to find the min/max 3-Body energy.
                boolean validPair = minMax2BodySum(residues, minMaxTriple, i, ri, j, rj);
                if (!validPair) {
                    // Eliminate Rotamer Pair
                    Residue residuei = residues[i];
                    rO.logIfMaster(format(" Inconsistent Pair: %8s %2d, %8s %2d.",
                            residuei.toFormattedString(false, true), ri,
                            residuej.toFormattedString(false, true), rj), Level.INFO);
                    continue;
                }

                if (Double.isFinite(currMin) && Double.isFinite(minMaxTriple[0])) {
                    currMin += minMaxTriple[0];
                } else {
                    currMin = Double.NaN;
                }

                if (Double.isFinite(currMax) && Double.isFinite(minMaxTriple[1])) {
                    currMax += minMaxTriple[1];
                } else {
                    currMax = Double.NaN;
                }
            }

            valid = true;
            if (Double.isFinite(currMin) && currMin < minMax[0]) {
                minMax[0] = currMin;
            } // Else, we do not have a new minimum.

            if (Double.isFinite(currMax) && Double.isFinite(minMax[1])) {
                if (currMax > minMax[1]) {
                    // We have a new, finite maximum.
                    minMax[1] = currMax;
                }  // Else, if currMax is finite and less than minMax[1], we do not have a new maximum.
            } else {
                // We have a non-finite maximum.
                minMax[1] = Double.NaN;
            }
        }

        // minMax[0] being set to NaN should be redundant with valid being false.
        // It would indicate i,ri clashes with something in every possible configuration.
        minMax[0] = (minMax[0] == Double.MAX_VALUE) ? Double.NaN : minMax[0];
        // minMax[1] always gets set, unless somehow everything turns up as Double.MIN_VALUE.
        return valid;
    }

    /**
     * Calculates the minimum and maximum summations over additional residues
     * for some pair ri-rj.
     *
     * @param residues Residues under consideration.
     * @param minMax   Result array: 0 is min summation, 1 max summation.
     * @param i        Residue i.
     * @param ri       Rotamer for residue i.
     * @param j        Residue j!=i.
     * @param rj       Rotamer for residue j.
     * @return False if ri-rj always clashes with other residues.
     * @throws IllegalArgumentException If ri, rj, or ri-rj eliminated.
     */
    public boolean minMaxE2(Residue[] residues, double[] minMax, int i, int ri, int j, int rj)
            throws IllegalArgumentException {
        Residue resi = residues[i];
        Residue resj = residues[j];
        if (eR.check(i, ri) || eR.check(j, rj) || eR.check(i, ri, j, rj)) {
            throw new IllegalArgumentException(format(" Called for minMaxE2 on an eliminated pair %s-%d %s-%d", resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj));
        }

        // Minimum summation over third residues k.
        minMax[0] = 0;
        // Maximum summation over third residues k.
        minMax[1] = 0;

        int nRes = residues.length;
        for (int k = 0; k < nRes; k++) {
            if (k == i || k == j) {
                continue;
            }
            Residue resk = residues[k];
            Rotamer[] rotsk = resk.getRotamers(library);
            int lenrk = rotsk.length;
            double[] minMaxK = new double[2];
            minMaxK[0] = Double.MAX_VALUE;
            minMaxK[1] = Double.MIN_VALUE;

            for (int rk = 0; rk < lenrk; rk++) {
                if (eR.check(k, rk)) {
                    // Not a valid part of phase space.
                    continue;
                }
                if (eR.check(i, ri, k, rk) || eR.check(j, rj, k, rk)) {
                    // Not implemented: check(i, ri, j, rj, k, rk).

                    // i,ri or j,rj clashes with this rotamer, max will be NaN.
                    // Minimum for this rk will be a clash, which is never a minimum.
                    minMaxK[1] = Double.NaN;
                } else {

                    // Min and max summations over 4th residues l, plus the ri-rk and rj-rk interactions.
                    // If no 3-body term, just the ri-rk and rj-rk interactions.
                    double currentMin = get2Body(i, ri, k, rk) + get2Body(j, rj, k, rk);
                    double currentMax = currentMin;
                    if (threeBodyTerm) {
                        // If the 3-Body eliminated, would fill max to Double.NaN.
                        currentMin += get3Body(residues, i, ri, j, rj, k, rk);
                        currentMax = currentMin;

                        // Obtain min and max summations over l.
                        double[] minMaxTriple = new double[2];
                        if (minMaxE3(residues, minMaxTriple, i, ri, j, rj, k, rk)) {
                            // A non-finite triples minimum should have the code taking the else branch.
                            assert (Double.isFinite(minMaxTriple[0]) && minMaxTriple[0] != Double.MAX_VALUE);

                            // Add the min and max summations over all 4th residues l.
                            currentMin += minMaxTriple[0];

                            if (Double.isFinite(currentMax) && Double.isFinite(minMaxTriple[1])) {
                                currentMax += minMaxTriple[1];
                            } else {
                                currentMax = Double.NaN;
                            }
                        } else {
                            // i, ri, j, rj, k, rk creates an inevitable clash with some residue l.
                            currentMin = Double.NaN;
                            currentMax = Double.NaN;
                        }
                    }

                    assert (threeBodyTerm || currentMax == currentMin);

                    // Now check if rk displaces previously searched rk for min/max over this k.
                    if (Double.isFinite(currentMin) && currentMin < minMaxK[0]) {
                        // rk has a more favorable minimum than previously searched rk.
                        minMaxK[0] = currentMin;
                    } // Else, no new minimum found.

                    if (Double.isFinite(currentMax) && Double.isFinite(minMaxK[1])) {
                        // rk has a less favorable maximum than previously searched rk.
                        minMaxK[1] = Math.max(currentMax, minMaxK[1]);
                    } else {
                        // Our maximum is a NaN.
                        minMaxK[1] = Double.NaN;
                    }
                }
            }

            if (Double.isFinite(minMaxK[0])) {
                // Add the minimum contribution from this k to the summation.
                minMax[0] += minMaxK[0];
            } else {
                // Else, ri-rj conflicts with all rk for this k, and can be swiftly eliminated.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            }
            if (Double.isFinite(minMaxK[1]) && Double.isFinite(minMax[1])) {
                // Add the max contribution from this k to the summation.
                minMax[1] += minMaxK[1];
            } else {
                // Otherwise, the max for ri-rj is a clash.
                minMax[1] = Double.NaN;
            }
        }

        return Double.isFinite(minMax[0]);
    }

    /**
     * Calculates the minimum and maximum summations over additional residues
     * for some rotamer triples ri-rj-rk.
     *
     * @param residues Residues under consideration.
     * @param minMax   Result array: 0 is min summation, 1 max summation.
     * @param i        Residue i.
     * @param ri       Rotamer for residue i.
     * @param j        Residue j!=i.
     * @param rj       Rotamer for residue j.
     * @param k        Residue k!=j and k!=i.
     * @param rk       Rotamer for residue k.
     * @return False if ri-rj-rk always clashes with other residues.
     * @throws IllegalArgumentException if there are pre-existing eliminations
     *                                  in ri-rj-rk.
     */
    private boolean minMaxE3(Residue[] residues, double[] minMax, int i, int ri, int j, int rj, int k, int rk)
            throws IllegalArgumentException {
        Residue resi = residues[i];
        Residue resj = residues[j];
        Residue resk = residues[k];
        if (eR.check(i, ri) || eR.check(j, rj) || eR.check(k, rk) ||
                eR.check(i, ri, j, rj) || eR.check(i, ri, k, rk) || eR.check(j, rj, k, rk)) {
            // Not implemented: check(i, ri, j, rj, k, rk).
            throw new IllegalArgumentException(format(" Called for minMaxE2 on an eliminated triple %s-%d %s-%d %s-%d", resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk));
        }

        // These two are a summation of mins/maxes over all fourth residues l.
        minMax[0] = 0;
        minMax[1] = 0;
        int nRes = residues.length;
        for (int l = 0; l < nRes; l++) {
            if (l == i || l == j || l == k) {
                continue;
            }
            Residue resl = residues[l];
            Rotamer[] rotsl = resl.getRotamers(library);
            int lenrl = rotsl.length;

            // Find min/max rl for residue l.
            double currentMax = Double.MIN_VALUE;
            double currentMin = Double.MAX_VALUE;

            for (int rl = 0; rl < lenrl; rl++) {
                if (eR.check(l, rl) || eR.check(k, rk, l, rl)) {
                    // Not valid phase space for anything.
                    continue;
                }

                double current;
                if (eR.check(i, ri, l, rl) || eR.check(j, rj, l, rl)) {
                    // Not implemented: checking ri-rj-rl, ri-rk-rl, rj-rk-rl, or ri-rj-rk-rl.
                    current = Double.NaN;
                } else {
                    // ri-rj-rl is accounted for at a different part of the summation as ri-rj-rk.
                    current = get3Body(residues, i, ri, k, rk, l, rl) + get3Body(residues, j, rj, k, rk, l, rl);
                }

                // TODO: Add quads to the DEE summation.
                // Would have to replace "current" with array "currentQuads".
                //double[] minMaxQuads;
                // minMaxE4(args)
                if (Double.isFinite(current) && current < currentMin) {
                    // rl forms a more favorable 3-body than any prior rl for this residue l.
                    currentMin = current;
                }

                if (Double.isFinite(current) && Double.isFinite(currentMax)) {
                    if (current > currentMax) {
                        currentMax = current;
                    } // Else, no new finite max found.
                } else {
                    currentMax = Double.NaN;
                }
            }

            if (Double.isFinite(currentMin)) {
                minMax[0] += currentMin;
            } else {
                // Else, ri-rj-rk inevitably conflicts with l.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            }

            if (Double.isFinite(currentMax) && Double.isFinite(minMax[1])) {
                minMax[1] += currentMax;
            } else {
                minMax[1] = Double.NaN;
            }
            // Finished with this residue l.
        }
        return Double.isFinite(minMax[0]);
    }

    public HashMap<Integer, Integer[]> getSelfEnergyMap() {
        return selfEnergyMap;
    }

    public HashMap<Integer, Integer[]> getTwoBodyEnergyMap() {
        return twoBodyEnergyMap;
    }

    public HashMap<Integer, Integer[]> getThreeBodyEnergyMap() {
        return threeBodyEnergyMap;
    }

    public HashMap<Integer, Integer[]> getFourBodyEnergyMap() {
        return fourBodyEnergyMap;
    }

    public HashMap<String, Integer> allocateSelfJobMap(Residue[] residues, int nResidues, boolean reverseMap) {
        selfEnergyMap.clear();
        // allocate selfEnergy array and create self jobs
        HashMap<String, Integer> reverseJobMapSingles = new HashMap<>();
        int singleJobIndex = 0;
        selfEnergy = new double[nResidues][];
        for (int i = 0; i < nResidues; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            selfEnergy[i] = new double[roti.length];
            for (int ri = 0; ri < roti.length; ri++) {
                if (!eR.check(i, ri)) {
                    Integer selfJob[] = {i, ri};
                    if (decomposeOriginal && ri != 0) {
                        continue;
                    }
                    selfEnergyMap.put(singleJobIndex, selfJob);
                    if (reverseMap) {
                        String revKey = format("%d %d", i, ri);
                        reverseJobMapSingles.put(revKey, singleJobIndex);
                    }
                    singleJobIndex++;
                }
            }
        }
        return reverseJobMapSingles;
    }

    public HashMap<String, Integer> allocate2BodyJobMap(Residue[] residues, int nResidues, boolean reverseMap) {
        twoBodyEnergyMap.clear();
        // allocated twoBodyEnergy array and create pair jobs
        HashMap<String, Integer> reverseJobMapPairs = new HashMap<>();
        int pairJobIndex = 0;
        twoBodyEnergy = new double[nResidues][][][];
        for (int i = 0; i < nResidues; i++) {
            Residue resi = residues[i];
            int indexI = allResiduesList.indexOf(resi);
            Rotamer roti[] = resi.getRotamers(library);
            int[] nI = resNeighbors[i];
            int lenNI = nI.length;
            twoBodyEnergy[i] = new double[roti.length][lenNI][];

            for (int ri = 0; ri < roti.length; ri++) {
                if (eR.check(i, ri)) {
                    continue;
                }
                //for (int j = i + 1; j < nResidues; j++) {
                for (int indJ = 0; indJ < lenNI; indJ++) {
                    int j = nI[indJ];
                    if (rO.checkNeighboringPair(i, j)) {
                        Residue resj = residues[j];
                        int indexJ = allResiduesList.indexOf(resj);
                        Rotamer rotj[] = resj.getRotamers(library);
                        twoBodyEnergy[i][ri][indJ] = new double[rotj.length];
                        for (int rj = 0; rj < rotj.length; rj++) {
                            if (eR.checkToJ(i, ri, j, rj)) {
                                continue;
                            }

                            // Skip creating a job if the pair is outside pair cut-off.
                            if (dM.checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                                continue;
                            }

                            Integer[] pairJob = {i, ri, j, rj};
                            if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                continue;
                            }
                            twoBodyEnergyMap.put(pairJobIndex, pairJob);
                            if (reverseMap) {
                                String revKey = format("%d %d %d %d", i, ri, j, rj);
                                reverseJobMapPairs.put(revKey, pairJobIndex);
                            }
                            pairJobIndex++;
                        }
                    }
                }
            }
        }
        return reverseJobMapPairs;
    }

    public HashMap<String, Integer> allocate3BodyJobMap(Residue[] residues, int nResidues, boolean reverseMap) {
        HashMap<String, Integer> reverseJobMapTrimers = new HashMap<>();
        threeBodyEnergyMap.clear();
        // fill in 3-Body energies from the restart file.
        threeBodyEnergy = new double[nResidues][][][][][];
        int trimerJobIndex = 0;
        for (int i = 0; i < nResidues; i++) {
            Residue resi = residues[i];
            int indexI = allResiduesList.indexOf(resi);
            Rotamer[] roti = resi.getRotamers(library);
            int lenri = roti.length;
            int[] nI = resNeighbors[i];
            int lenNI = nI.length;
            threeBodyEnergy[i] = new double[lenri][lenNI][][][];

            for (int ri = 0; ri < lenri; ri++) {
                if (eR.check(i, ri)) {
                    continue;
                }
                for (int indJ = 0; indJ < lenNI; indJ++) {
                    //for (int j = i + 1; j < nResidues; j++) {
                    int j = nI[indJ];
                    Residue resj = residues[j];
                    int indexJ = allResiduesList.indexOf(resj);
                    Rotamer[] rotj = resj.getRotamers(library);
                    int lenrj = rotj.length;
                    int[] nJ = resNeighbors[j];
                    int lenNJ = nJ.length;
                    threeBodyEnergy[i][ri][indJ] = new double[lenrj][lenNJ][];

                    for (int rj = 0; rj < lenrj; rj++) {
                        if (eR.checkToJ(i, ri, j, rj)) {
                            continue;
                        }
                        //for (int k = j + 1; k < nResidues; k++) {
                        for (int indK = 0; indK < lenNJ; indK++) {
                            int k = nJ[indK];
                            Residue resk = residues[k];
                            int indexK = allResiduesList.indexOf(resk);
                            Rotamer[] rotk = resk.getRotamers(library);
                            int lenrk = rotk.length;
                            threeBodyEnergy[i][ri][indJ][rj][indK] = new double[lenrk];

                            for (int rk = 0; rk < lenrk; rk++) {
                                if (eR.checkToK(i, ri, j, rj, k, rk)) {
                                    continue;
                                }
                                if (dM.checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
                                    continue;
                                }
                                Integer[] trimerJob = {i, ri, j, rj, k, rk};
                                if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                    continue;
                                }
                                threeBodyEnergyMap.put(trimerJobIndex, trimerJob);
                                if (reverseMap) {
                                    String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                                    reverseJobMapTrimers.put(revKey, trimerJobIndex);
                                }
                                trimerJobIndex++;
                            }
                        }
                    }
                }
            }
        }
        return reverseJobMapTrimers;
    }

    public void allocate4BodyJobMap(Residue[] residues, int nResidues) {
        logger.info(" Creating 4-Body jobs...");
        fourBodyEnergyMap.clear();
        boolean maxedOut = false;
        // create 4-Body jobs (no memory allocation)
        int quadJobIndex = 0;
        for (int i = 0; i < nResidues; i++) {
            Residue resi = residues[i];
            Rotamer[] roti = resi.getRotamers(library);
            for (int ri = 0; ri < roti.length; ri++) {
                if (eR.check(i, ri)) {
                    continue;
                }
                for (int j = i + 1; j < nResidues; j++) {
                    Residue resj = residues[j];
                    Rotamer[] rotj = resj.getRotamers(library);
                    for (int rj = 0; rj < rotj.length; rj++) {
                                /*if (check(j, rj) || check(i, ri, j, rj)) {
                                 continue;
                                 }*/
                        if (eR.checkToJ(i, ri, j, rj)) {
                            continue;
                        }
                        for (int k = j + 1; k < nResidues; k++) {
                            Residue resk = residues[k];
                            Rotamer[] rotk = resk.getRotamers(library);
                            for (int rk = 0; rk < rotk.length; rk++) {
                                        /*if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk) || check(i, ri, j, rj, k, rk)) {
                                         continue;
                                         }*/
                                if (eR.checkToK(i, ri, j, rj, k, rk)) {
                                    continue;
                                }
                                for (int l = k + 1; l < nResidues; l++) {
                                    Residue resl = residues[l];
                                    Rotamer[] rotl = resl.getRotamers(library);
                                    for (int rl = 0; rl < rotl.length; rl++) {
                                        if (eR.checkToL(i, ri, j, rj, k, rk, l, rl)) {
                                            continue;
                                        }
                                        Integer[] quadJob = {i, ri, j, rj, k, rk, l, rl};
                                        if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0 || rl != 0)) {
                                            continue;
                                        }
                                        fourBodyEnergyMap.put(quadJobIndex++, quadJob);
                                        if (fourBodyEnergyMap.size() > max4BodyCount) {
                                            maxedOut = true;
                                            break;
                                        }
                                    }
                                    if (maxedOut) {
                                        break;
                                    }
                                }
                                if (maxedOut) {
                                    break;
                                }
                            }
                            if (maxedOut) {
                                break;
                            }
                        }
                        if (maxedOut) {
                            break;
                        }
                    }
                    if (maxedOut) {
                        break;
                    }
                }
                if (maxedOut) {
                    break;
                }
            }
            if (maxedOut) {
                break;
            }
        }
    }

    public int loadEnergyRestart(File restartFile, Residue[] residues) {
        return loadEnergyRestart(restartFile, residues, -1, null);
    }

    public int loadEnergyRestart(File restartFile, Residue[] residues, int boxIteration, int[] cellIndices) {
        try {
            int nResidues = residues.length;
            Path path = Paths.get(restartFile.getCanonicalPath());
            List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
            List<String> linesThisBox = new ArrayList<>();

            try {
                backboneEnergy = rO.computeBackboneEnergy(residues);
            } catch (ArithmeticException ex) {
                logger.severe(format(" Exception %s in calculating backbone energy; FFX shutting down.", ex.toString()));
            }
            rO.logIfMaster(format(" Backbone energy:  %s\n", rO.formatEnergy(backboneEnergy)));

            if (usingBoxOptimization && boxIteration >= 0) {
                boolean foundBox = false;
                for (int i = 0; i < lines.size(); i++) {
                    String line = lines.get(i);
                    if (line.startsWith("Box")) {
                        String tok[] = line.replaceAll("Box", "").replaceAll(":", ",").replaceAll(" ", "").split(",");
                        int readIteration = Integer.parseInt(tok[0]);
                        int readCellIndexX = Integer.parseInt(tok[1]);
                        int readCellIndexY = Integer.parseInt(tok[2]);
                        int readCellIndexZ = Integer.parseInt(tok[3]);
                        if (readIteration == boxIteration
                                && readCellIndexX == cellIndices[0]
                                && readCellIndexY == cellIndices[1]
                                && readCellIndexZ == cellIndices[2]) {
                            foundBox = true;
                            for (int j = i + 1; j < lines.size(); j++) {
                                String l = lines.get(j);
                                if (l.startsWith("Box")) {
                                    break;
                                }
                                linesThisBox.add(l);
                            }
                            break;
                        }
                    }
                }
                if (!foundBox) {
                    rO.logIfMaster(format(" Didn't find restart energies for Box %d: %d,%d,%d",
                            boxIteration, cellIndices[0], cellIndices[1], cellIndices[2]));
                    return 0;
                } else if (linesThisBox.size() == 0) {
                    return 0;
                } else {
                    lines = linesThisBox;
                }
            }

            List<String> singleLines = new ArrayList<>();
            List<String> pairLines = new ArrayList<>();
            List<String> tripleLines = new ArrayList<>();
            for (String line : lines) {
                String tok[] = line.split("\\s");
                if (tok[0].startsWith("Self")) {
                    singleLines.add(line);
                } else if (tok[0].startsWith("Pair")) {
                    pairLines.add(line);
                } else if (tok[0].startsWith("Triple")) {
                    tripleLines.add(line);
                }
            }
            int loaded = 0;
            if (tripleLines.size() > 0) {
                loaded = 3;
            } else if (pairLines.size() > 0) {
                loaded = 2;
            } else if (singleLines.size() > 0) {
                loaded = 1;
            } else {
                logger.warning(format(" Empty or unreadable energy restart file: %s.", restartFile.getCanonicalPath()));
            }
            if (loaded >= 1) {
                boolean reverseMap = true;
                HashMap<String, Integer> reverseJobMapSingles = allocateSelfJobMap(residues, nResidues, reverseMap);
                // fill in self-energies from file while removing the corresponding jobs from selfEnergyMap
                for (String line : singleLines) {
                    try {
                        String tok[] = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        double energy = Double.parseDouble(tok[3]);
                        try {
                            setSelf(i, ri, energy);
                            if (verbose) {
                                rO.logIfMaster(format(" From restart file: Self energy %3d (%8s,%2d): %s",
                                        i, residues[i].toFormattedString(false, true), ri,
                                        rO.formatEnergy(energy)));
                            }
                        } catch (Exception e) {
                            if (verbose) {
                                rO.logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d", i, ri);
                        Integer ret[] = selfEnergyMap.remove(reverseJobMapSingles.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format(" Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                rO.logIfMaster(" Loaded self energies from restart file.");

                //Pre-Prune if self-energy is Double.NaN.
                eR.prePruneSelves(residues);

                // prune singles
                if (pruneClashes) {
                    eR.pruneSingleClashes(residues);
                }
            }

            // Remap to sequential integer keys.
            condenseEnergyMap(selfEnergyMap);

            if (loaded >= 2) {
                if (selfEnergyMap.size() > 0) {
                    rO.logIfMaster(" Double-check that parameters match original run due to missing self-energies.");
                }
                boolean reverseMap = true;
                HashMap<String, Integer> reverseJobMapPairs = allocate2BodyJobMap(residues, nResidues, reverseMap);
                // fill in pair-energies from file while removing the corresponding jobs from twoBodyEnergyMap
                for (String line : pairLines) {
                    try {
                        String tok[] = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        int j;
                        if (tok[3].contains("-")) {
                            j = nameToNumber(tok[3], residues);
                        } else {
                            j = Integer.parseInt(tok[3]);
                        }
                        int rj = Integer.parseInt(tok[4]);
                        double energy = Double.parseDouble(tok[5]);
                        try {
                            // When a restart file is generated using a large cutoff, but a new simulation is being done
                            // with a smaller cutoff, the two-body distance needs to be checked. If the two-body
                            // distance is larger than the cutoff, then the two residues are not considered 'neighbors'
                            // so that pair should not be added to the pairs map.
                            if (rO.checkNeighboringPair(i, j)) {
                                //If inside the cutoff, set energy to previously computed value.
                                //Gather distances and indices for printing.
                                Residue residueI = residues[i];
                                Residue residueJ = residues[j];
                                int indexI = allResiduesList.indexOf(residueI);
                                int indexJ = allResiduesList.indexOf(residueJ);
                                if (!dM.checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                                    set2Body(i, ri, j, rj, energy);

                                    double resDist = dM.getResidueDistance(indexI, ri, indexJ, rj);
                                    String resDistString = format("large");
                                    if (resDist < Double.MAX_VALUE) {
                                        resDistString = format("%5.3f", resDist);
                                    }

                                    double dist = dM.checkDistMatrix(indexI, ri, indexJ, rj);
                                    String distString = format("     large");
                                    if (dist < Double.MAX_VALUE) {
                                        distString = format("%10.3f", dist);
                                    }

                                    logger.fine(format(" Pair %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                                            residueI.toFormattedString(false, true), ri,
                                            residueJ.toFormattedString(false, true), rj,
                                            rO.formatEnergy(get2Body(i, ri, j, rj)), distString, resDistString));
                                }
                            } else {
                                logger.fine(format("Ignoring a pair-energy from outside the cutoff: 2-energy [(%8s,%2d),(%8s,%2d)]: %12.4f", residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj, energy));
                            }

                            if (verbose) {
                                rO.logIfMaster(format(" From restart file: Pair energy [(%8s,%2d),(%8s,%2d)]: %12.4f",
                                        residues[i].toFormattedString(false, true), ri,
                                        residues[j].toFormattedString(false, true), rj, energy));
                            }
                        } catch (Exception e) {
                            if (verbose) {
                                rO.logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d", i, ri, j, rj);
                        Integer[] ret = twoBodyEnergyMap.remove(reverseJobMapPairs.get(revKey));
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                rO.logIfMaster(" Loaded 2-body energies from restart file.");

                // Pre-Prune if pair-energy is Double.NaN.
                eR.prePrunePairs(residues);

                // prune pairs
                if (prunePairClashes) {
                    eR.prunePairClashes(residues);
                }
            }

            // Remap to sequential integer keys.
            condenseEnergyMap(twoBodyEnergyMap);

            if (loaded >= 3) {
                if (twoBodyEnergyMap.size() > 0) {
                    if (master) {
                        logger.warning("Double-check that parameters match original run!  Found trimers in restart file, but pairs job queue is non-empty.");
                    }
                }
                boolean reverseMap = true;
                HashMap<String, Integer> reverseJobMapTrimers = allocate3BodyJobMap(residues, nResidues, reverseMap);

                // fill in 3-Body energies from file while removing the corresponding jobs from threeBodyEnergyMap
                for (String line : tripleLines) {
                    try {
                        String[] tok = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        int j;
                        if (tok[3].contains("-")) {
                            j = nameToNumber(tok[3], residues);
                        } else {
                            j = Integer.parseInt(tok[3]);
                        }
                        int rj = Integer.parseInt(tok[4]);
                        int k;
                        if (tok[5].contains("-")) {
                            k = nameToNumber(tok[5], residues);
                        } else {
                            k = Integer.parseInt(tok[5]);
                        }
                        int rk = Integer.parseInt(tok[6]);
                        double energy = Double.parseDouble(tok[7]);

                        try {
                            //threeBodyEnergy[i][ri][j][rj][k][rk] = energy;
                            //IntegerKeyset ijk = new IntegerKeyset(i, ri, j, rj, k, rk);
                            //threeBodyEnergies.put(ijk, energy);

                            // When a restart file is generated using a large cutoff, but a new simulation is being done
                            // with a smaller cutoff, the three-body distance needs to be checked. If the three-body
                            // distance is larger than the cutoff, then the three residues are not considered 'neighbors'
                            // so that triple should not be added to the pairs map.
                            if (rO.checkNeighboringTriple(i, j, k)) {
                                //If within the cutoff, the energy should be set to the previously calculated energy.
                                Residue residueI = residues[i];
                                Residue residueJ = residues[j];
                                Residue residueK = residues[k];
                                int indexI = allResiduesList.indexOf(residueI);
                                int indexJ = allResiduesList.indexOf(residueJ);
                                int indexK = allResiduesList.indexOf(residueK);
                                if (!dM.checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
                                    set3Body(residues, i, ri, j, rj, k, rk, energy);

                                    double rawDist = dM.getRawNBodyDistance(indexI, ri, indexJ, rj, indexK, rk);
                                    double resDist = dM.get3BodyResidueDistance(indexI, ri, indexJ, rj, indexK, rk);

                                    String resDistString = "     large";
                                    if (resDist < Double.MAX_VALUE) {
                                        resDistString = format("%5.3f", resDist);
                                    }

                                    String distString = "     large";
                                    if (rawDist < Double.MAX_VALUE) {
                                        distString = format("%10.3f", rawDist);
                                    }

                                    logger.fine(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                                            residueI.toFormattedString(false, true), ri,
                                            residueJ.toFormattedString(false, true), rj,
                                            residueK.toFormattedString(false, true), rk,
                                            rO.formatEnergy(get3Body(residues, i, ri, j, rj, k, rk)), distString, resDistString));
                                }
                            } else {
                                logger.fine(format("Ignoring a triple-energy from outside the cutoff: 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s",
                                        residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj, residues[k].toFormattedString(false, true), rk,
                                        rO.formatEnergy(get3Body(residues, i, ri, j, rj, k, rk))));
                            }
                        } catch (ArrayIndexOutOfBoundsException ex) {
                            if (verbose) {
                                rO.logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        } catch (NullPointerException npe) {
                            if (verbose) {
                                rO.logIfMaster(format(" NPE in loading 3-body energies: pruning "
                                                + "likely changed! 3-body %s-%d %s-%d %s-%d", residues[i].toFormattedString(false, true),
                                        ri, residues[j], rj, residues[k], rk));
                            }
                        }
                        if (verbose) {
                            rO.logIfMaster(format(" From restart file: Trimer energy %3d %-2d, %3d %-2d, %3d %-2d: %s",
                                    i, ri, j, rj, k, rk, rO.formatEnergy(energy)));
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                        Integer ret[] = threeBodyEnergyMap.remove(reverseJobMapTrimers.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                rO.logIfMaster(" Loaded trimer energies from restart file.");
            }

            // Remap to sequential integer keys.
            condenseEnergyMap(threeBodyEnergyMap);

            return loaded;
        } catch (IOException ex) {
            logger.log(Level.WARNING, "Exception while loading energy restart file.", ex);
        }

        return 0;
    }

    public void turnOnResidue(Residue residue, int ri) {
        turnOnAtoms(residue);
        Rotamer rotamers[] = residue.getRotamers(library);
        RotamerLibrary.applyRotamer(residue, rotamers[ri]);
    }

    public void turnOffResidue(Residue residue) {
        turnOffAtoms(residue);
        applyDefaultRotamer(residue);
    }

    /**
     * Applies the "default" rotamer: currently the 0'th rotamer.
     *
     * @param residue Residue to apply a default rotamer for.
     */
    private void applyDefaultRotamer(Residue residue) {
        RotamerLibrary.applyRotamer(residue, residue.getRotamers(library)[0]);
    }

    public void turnOffAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            turnOffResidue(residues[i]);
        }
    }

    public void turnOnAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
    }

    /**
     * Set the "use" flag to true for all variable atoms in a residue.
     *
     * @param residue The Residue to turn on.
     */
    private static void turnOnAtoms(Residue residue) {
        switch (residue.getResidueType()) {
            case NA:
            case AA:
                List<Atom> atomList = residue.getVariableAtoms();
                for (Atom atom : atomList) {
                    atom.setUse(true);
                }
                break;
            default:
                atomList = residue.getAtomList();
                for (Atom atom : atomList) {
                    atom.setUse(true);
                }
                break;
        }
    }

    /**
     * Set the "use" flag to true for all variable atoms in a residue.
     *
     * @param residue The residue to turn off.
     */
    private static void turnOffAtoms(Residue residue) {
        switch (residue.getResidueType()) {
            case NA:
            case AA:
                List<Atom> atomList = residue.getVariableAtoms();
                for (Atom atom : atomList) {
                    atom.setUse(false);
                }
                break;
            default:
                atomList = residue.getAtomList();
                for (Atom atom : atomList) {
                    atom.setUse(false);
                }
                break;
        }
    }

    private int nameToNumber(String residueString, Residue[] residues) throws NumberFormatException {
        int ret = -1;
        for (int x = 0; x < residues.length; x++) {
            if (residueString.equals(residues[x].toString())) {
                ret = x;
                break;
            }
        }
        if (ret == -1) {
            throw new NumberFormatException();
        }
        return ret;
    }

    private void condenseEnergyMap(HashMap<Integer, Integer[]> energyMap) {
        Set<Integer> keys = energyMap.keySet();
        HashMap<Integer, Integer[]> tempMap = new HashMap<>();
        int count = 0;
        for (int key : keys) {
            tempMap.put(count, energyMap.get(key));
            count++;
        }
        energyMap.clear();
        energyMap.putAll(tempMap);
    }
}
