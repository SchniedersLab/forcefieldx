/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.PolarizeType;

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link NeighborList} for any
 * {@link Crystal} space group. The real space contribution is contained within
 * this class and the reciprocal space contribution is delegated to the
 * {@link ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders
 */
public abstract class ParticleMeshEwald implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    protected Atom atoms[];
    /**
     * The number of atoms in the system.
     */
    protected int nAtoms;
    /**
     * Polarization groups.
     */
    protected int ip11[][];
    protected int ip12[][];
    protected int ip13[][];

    /**
     * Total electrostatic energy == permanentMultipole + polarizationEnergy +
     * generalizedKirkwood.
     */
    protected double totalElectrostaticEnergy;
    /**
     * Total multipole energy == permanentMultipole + polarizationEnergy (no
     * GK).
     */
    protected double totalMultipoleEnergy;
    /**
     * Permanent multipole energy in kcal/mol == permanentRealSpace +
     * permanentSelf + permanentReciprocal.
     */
    protected double permanentMultipoleEnergy;
    protected double permanentRealSpaceEnergy;
    protected double permanentSelfEnergy;
    protected double permanentReciprocalEnergy;
    /**
     * Polarization energy in kcal/mol == inducedRealSpace + inducedSelf +
     * inducedReciprocal.
     */
    protected double polarizationEnergy;
    protected double inducedRealSpaceEnergy;
    protected double inducedSelfEnergy;
    protected double inducedReciprocalEnergy;
    /**
     * Solvation energy due to implicit GK solvent.
     */
    protected double generalizedKirkwoodEnergy;
    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public Polarization polarization;

    protected SCFAlgorithm scfAlgorithm = SCFAlgorithm.CG;

    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    public double coordinates[][][];
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    public int neighborLists[][][];

    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    public double globalMultipole[][][];

    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double inducedDipole[][][];
    public double inducedDipoleCR[][][];


    /**
     * Log the induced dipole magnitudes and directions. Use the cgo_arrow.py
     * script (available from the wiki) to draw these easily in PyMol.
     */
    public boolean printInducedDipoles = Boolean.valueOf(System.getProperty("pme.printInducedDipoles", "false"));
    /**
     * Log the seven components of total electrostatic energy at each
     * evaluation: (Permanent) PermanentRealSpace, PermanentSelf, PermanentRecip
     * (Induced) InducedRealSpace, InducedSelf, InducedRecip, and
     * GeneralizedKirkwood. Self, Recip terms apply only to periodic systems; GK
     * applies only when requested and aperiodic.
     */
    public boolean printDecomposition;
    /**
     * Disables windowed lambda ranges by setting permLambdaStart =
     * polLambdaStart = 0.0 and permLambdaEnd = polLambdaEnd = 1.0.
     */
    protected final boolean noWindowing = Boolean.valueOf(System.getProperty("pme.noWindowing", "false"));

    /**
     * Default cutoff values for PME and aperiodic systems.
     */
    public static final double PERIODIC_DEFAULT_EWALD_CUTOFF = 7.0;
    public static final double APERIODIC_DEFAULT_EWALD_CUTOFF = 1000.0;

    public enum Polarization {
        MUTUAL, DIRECT, NONE
    }

    public void setPolarization(Polarization set) {
        this.polarization = set;
    }

    public enum ELEC_FORM {
        PAM, FIXED_CHARGE
    }

    public enum LambdaMode {
        OFF, CONDENSED, CONDENSED_NO_LIGAND, VAPOR
    }

    /**
     * Describes available SCF algorithms, and whether they are supported by the
     * FFX and/or CUDA implementations.
     */
    public enum SCFAlgorithm {
        // I actually don't know if OpenMM does SOR or CG, but they both just become "Mutual".
        SOR(true, true), CG(true, true), EPT(false, true);

        private final List<Platform> supportedPlatforms;

        SCFAlgorithm(boolean ffx, boolean openMM, Platform... otherPlatforms) {
            List<Platform> platforms = new ArrayList<>();
            if (ffx) {
                platforms.add(Platform.FFX);
            }
            if (openMM) {
                platforms.add(Platform.OMM);
                platforms.add(Platform.OMM_CUDA);
                platforms.add(Platform.OMM_REF);
                // Inapplicable because the OpenCL and optimized C implementations don't have AMOEBA yet.
                //platforms.add(ForceFieldEnergy.Platform.OMM_OPENCL);
                //platforms.add(ForceFieldEnergy.Platform.OMM_OPTCPU);
            }
            platforms.addAll(Arrays.asList(otherPlatforms));
            supportedPlatforms = Collections.unmodifiableList(platforms);
            // Short of reflection, a final unmodifiable list is, well, unmodifiable.
        }

        /**
         * Returns the list of supported Platforms.
         *
         * @return The supported platform List. Unmodifiable.
         */
        public List<Platform> getSupportedPlatforms() {
            return supportedPlatforms;
        }

        /**
         * Checks if this platform is supported
         *
         * @param platform To check
         * @return Supported
         */
        public boolean isSupported(Platform platform) {
            return supportedPlatforms.contains(platform);
        }
    }

    public enum SCFPredictor {
        NONE, LS, POLY, ASPC
    }

    public enum Mask {
        Permanent, PolarGroup, PolarEnergy, Polar, All;
    }

    public abstract GeneralizedKirkwood getGK();

    public abstract double getEwaldCutoff();

    protected abstract double[][][] getGradient();

    protected abstract double[][][] getTorque();

    protected abstract double[][][] getLambdaGradient();

    protected abstract double[][][] getLambdaTorque();

    public abstract void setAtoms(Atom atoms[], int molecule[]);

    public abstract void setFixedCharges(Atom atoms[]);

    public abstract double energy(boolean gradient, boolean print);

    public abstract int getInteractions();

    public abstract int getGKInteractions();

    @Override
    public abstract void setLambda(double lambda);

    @Override
    public abstract double getdEdL();

    @Override
    public abstract void getdEdXdL(double[] gradients);

    @Override
    public abstract double getd2EdL2();

    public abstract void destroy() throws Exception;

    public abstract void setCrystal(Crystal crystal);

    public abstract double getCavitationEnergy(boolean throwError);

    public abstract double getDispersionEnergy(boolean throwError);

    public abstract double[][][] getCoordinates();

    public abstract double getPolarEps();

    public abstract int[][] getPolarization11();

    public abstract int[][] getPolarization12();

    public abstract int[][] getPolarization13();

    public abstract Polarization getPolarizationType();

    public abstract int[][] getAxisAtoms();

    public abstract double getScale14();

    public abstract double getEwaldCoefficient();

    public abstract ReciprocalSpace getReciprocalSpace();

    public abstract void setLambdaMultipoleScale(double scale);

    public abstract ELEC_FORM getElecForm();

    public abstract String getName();

    public double getTotalElectrostaticEnergy() {
        return permanentMultipoleEnergy + polarizationEnergy + generalizedKirkwoodEnergy;
    }

    public double getTotalMultipoleEnergy() {
        return permanentMultipoleEnergy + polarizationEnergy;
    }

    public double getPermanentEnergy() {
        return permanentMultipoleEnergy;
    }

    public double getPermRealEnergy() {
        return permanentRealSpaceEnergy;
    }

    public double getPermSelfEnergy() {
        return permanentSelfEnergy;
    }

    public double getPermRecipEnergy() {
        return permanentReciprocalEnergy;
    }

    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    public double getIndRealEnergy() {
        return inducedRealSpaceEnergy;
    }

    public double getIndSelfEnergy() {
        return inducedSelfEnergy;
    }

    public double getIndRecipEnergy() {
        return inducedReciprocalEnergy;
    }

    public double getGKEnergy() {
        return generalizedKirkwoodEnergy;
    }

    /**
     * Returns the SCF algorithm in use.
     *
     * @return The SCF algorithm used.
     */
    public SCFAlgorithm getScfAlgorithm() {
        return scfAlgorithm;
    }

    protected void assignPolarizationGroups() {
        /**
         * Find directly connected group members for each atom.
         */
        List<Integer> group = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.getIndex() - 1 != i) {
                logger.severe(" Atom indexing is not consistent in PME.");
            }
        }
        for (Atom ai : atoms) {
            group.clear();
            Integer index = ai.getIndex() - 1;
            group.add(index);
            PolarizeType polarizeType = ai.getPolarizeType();
            if (polarizeType != null) {
                if (polarizeType.polarizationGroup != null) {
                    growGroup(group, ai);
                    Collections.sort(group);
                    ip11[index] = new int[group.size()];
                    int j = 0;
                    for (int k : group) {
                        ip11[index][j++] = k;
                    }
                } else {
                    ip11[index] = new int[group.size()];
                    int j = 0;
                    for (int k : group) {
                        ip11[index][j++] = k;
                    }
                }
            } else {
                String message = "The polarize keyword was not found for atom "
                        + (index + 1) + " with type " + ai.getType();
                logger.severe(message);
            }
        }
        /**
         * Find 1-2 group relationships.
         */
        int mask[] = new int[nAtoms];
        List<Integer> list = new ArrayList<>();
        List<Integer> keep = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {
            mask[i] = -1;
        }
        for (int i = 0; i < nAtoms; i++) {
            list.clear();
            for (int j : ip11[i]) {
                list.add(j);
                mask[j] = i;
            }
            keep.clear();
            for (int j : list) {
                Atom aj = atoms[j];
                ArrayList<Bond> bonds = aj.getBonds();
                for (Bond b : bonds) {
                    Atom ak = b.get1_2(aj);
                    int k = ak.getIndex() - 1;
                    if (mask[k] != i) {
                        keep.add(k);
                    }
                }
            }
            list.clear();
            for (int j : keep) {
                for (int k : ip11[j]) {
                    list.add(k);
                }
            }
            Collections.sort(list);
            ip12[i] = new int[list.size()];
            int j = 0;
            for (int k : list) {
                ip12[i][j++] = k;
            }
        }
        /**
         * Find 1-3 group relationships.
         */
        for (int i = 0; i < nAtoms; i++) {
            mask[i] = -1;
        }
        for (int i = 0; i < nAtoms; i++) {
            for (int j : ip11[i]) {
                mask[j] = i;
            }
            for (int j : ip12[i]) {
                mask[j] = i;
            }
            list.clear();
            for (int j : ip12[i]) {
                for (int k : ip12[j]) {
                    if (mask[k] != i) {
                        if (!list.contains(k)) {
                            list.add(k);
                        }
                    }
                }
            }
            ip13[i] = new int[list.size()];
            Collections.sort(list);
            int j = 0;
            for (int k : list) {
                ip13[i][j++] = k;
            }
        }
    }

    /**
     * A recursive method that checks all atoms bonded to the seed atom for
     * inclusion in the polarization group. The method is called on each newly
     * found group member.
     *
     * @param group XYZ indeces of current group members.
     * @param seed The bonds of the seed atom are queried for inclusion in the
     * group.
     */
    private void growGroup(List<Integer> group, Atom seed) {
        List<Bond> bonds = seed.getBonds();
        for (Bond bi : bonds) {
            Atom aj = bi.get1_2(seed);
            int tj = aj.getType();
            boolean added = false;
            PolarizeType polarizeType = seed.getPolarizeType();
            for (int type : polarizeType.polarizationGroup) {
                if (type == tj) {
                    Integer index = aj.getIndex() - 1;
                    if (!group.contains(index)) {
                        group.add(index);
                        added = true;
                        break;
                    }
                }
            }
            if (added) {
                growGroup(group, aj);
            }
        }
    }

}
