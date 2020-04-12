//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Collections.sort;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField.ELEC_FORM;
import ffx.potential.parameters.PolarizeType;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.utilities.Constants.DEFAULT_ELECTRIC;
import static ffx.utilities.Constants.ELEC_ANG_TO_DEBYE;

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link ffx.potential.nonbonded.NeighborList} for any
 * {@link ffx.crystal.Crystal} space group. The real space contribution is contained within
 * this class and the reciprocal space contribution is delegated to the
 * {@link ffx.potential.nonbonded.ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders
 */
public abstract class ParticleMeshEwald implements LambdaInterface {

    /**
     * Default cutoff values for PME and aperiodic systems.
     */
    public static final double PERIODIC_DEFAULT_EWALD_CUTOFF = 7.0;
    /**
     * Constant <code>APERIODIC_DEFAULT_EWALD_CUTOFF=1000.0</code>
     */
    public static final double APERIODIC_DEFAULT_EWALD_CUTOFF = 1000.0;
    private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());
    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public Polarization polarization;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    public double[][][] coordinates;
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    public int[][][] neighborLists;
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    public double[][][] globalMultipole;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double[][][] inducedDipole;
    public double[][][] inducedDipoleCR;
    /**
     * Direct induced dipoles.
     */
    public double[][] directDipole;
    public double[][] directDipoleCR;
    /**
     * Vacuum induced dipoles
     */
    public double[][][] vacuumInducedDipole;
    public double[][][] vacuumInducedDipoleCR;
    /**
     * Vacuum induced dipoles
     */
    public double[][] vacuumDirectDipole;
    public double[][] vacuumDirectDipoleCR;
    /**
     * Log the seven components of total electrostatic energy at each
     * evaluation: (Permanent) PermanentRealSpace, PermanentSelf, PermanentRecip
     * (Induced) InducedRealSpace, InducedSelf, InducedRecip, and
     * GeneralizedKirkwood. Self, Recip terms apply only to periodic systems; GK
     * applies only when requested and aperiodic.
     */
    public boolean printDecomposition = false;
    /**
     * Disables windowed lambda ranges by setting permLambdaStart =
     * polLambdaStart = 0.0 and permLambdaEnd = polLambdaEnd = 1.0.
     */
    public boolean noWindowing = false;
    public double electric = DEFAULT_ELECTRIC;
    /**
     * An ordered array of atoms in the system.
     */
    protected Atom[] atoms;
    /**
     * The number of atoms in the system.
     */
    protected int nAtoms;
    /**
     * Polarization groups.
     */
    protected int[][] ip11;
    protected int[][] ip12;
    protected int[][] ip13;
    /**
     * Total multipole energy = permanentMultipoleEnergy + polarizationEnergy.
     * <br>
     * This does not include GK.
     */
    protected double totalMultipoleEnergy;
    /**
     * Permanent multipole energy = permanentRealSpaceEnergy + permanentSelfEnergy + permanentReciprocalEnergy.
     */
    protected double permanentMultipoleEnergy;
    protected double permanentRealSpaceEnergy;
    protected double permanentSelfEnergy;
    protected double permanentReciprocalEnergy;
    /**
     * Polarization energy = inducedRealSpaceEnergy + inducedSelfEnergy + inducedReciprocalEnergy.
     */
    protected double polarizationEnergy;
    protected double inducedRealSpaceEnergy;
    protected double inducedSelfEnergy;
    protected double inducedReciprocalEnergy;
    /**
     * Solvation energy due to implicit GK solvent.
     */
    protected double generalizedKirkwoodEnergy;
    protected SCFAlgorithm scfAlgorithm = SCFAlgorithm.CG;
    /**
     * Log the induced dipole magnitudes and directions. Use the cgo_arrow.py
     * script (available from the wiki) to draw these easily in PyMol.
     */
    protected boolean printInducedDipoles = false;

    /**
     * Compute multipole moments for an array of atoms.
     *
     * @param activeAtoms Atom array to consider.
     * @param forceEnergy Force calculation of the electrostatic energy (rotate multipoles, perform SCF).
     */
    public void computeMoments(Atom[] activeAtoms, boolean forceEnergy) {
        // Zero out total charge, dipole and quadrupole components.
        var netchg = 0.0;
        var netdpl = 0.0;
        var xdpl = 0.0;
        var ydpl = 0.0;
        var zdpl = 0.0;
        var xxqdp = 0.0;
        var xyqdp = 0.0;
        var xzqdp = 0.0;
        var yxqdp = 0.0;
        var yyqdp = 0.0;
        var yzqdp = 0.0;
        var zxqdp = 0.0;
        var zyqdp = 0.0;
        var zzqdp = 0.0;

        // Find the center of mass of the set of active atoms.
        double xmid = 0.0;
        double ymid = 0.0;
        double zmid = 0.0;
        double totalMass = 0;
        for (Atom atom : activeAtoms) {
            var m = atom.getMass();
            totalMass += m;
            xmid = xmid + atom.getX() * m;
            ymid = ymid + atom.getY() * m;
            zmid = zmid + atom.getZ() * m;
        }
        if (totalMass > 0) {
            xmid /= totalMass;
            ymid /= totalMass;
            zmid /= totalMass;
        }
        int n = activeAtoms.length;
        double[] xcm = new double[n];
        double[] ycm = new double[n];
        double[] zcm = new double[n];
        int index = 0;
        for (Atom atom : activeAtoms) {
            xcm[index] = atom.getX() - xmid;
            ycm[index] = atom.getY() - ymid;
            zcm[index] = atom.getZ() - zmid;
            index++;
        }

        if (forceEnergy) {
            energy(false, false);
        }

        // Account for charge, dipoles and induced dipoles.
        for (Atom atom : activeAtoms) {
            int i = atom.getIndex() - 1;
            double[] globalMultipolei = globalMultipole[0][i];
            double[] inducedDipolei = inducedDipole[0][i];

            var ci = globalMultipolei[t000];
            var dix = globalMultipolei[t100];
            var diy = globalMultipolei[t010];
            var diz = globalMultipolei[t001];
            var uix = inducedDipolei[0];
            var uiy = inducedDipolei[1];
            var uiz = inducedDipolei[2];

            netchg += ci;
            xdpl += xcm[i] * ci + dix + uix;
            ydpl += ycm[i] * ci + diy + uiy;
            zdpl += zcm[i] * ci + diz + uiz;
            xxqdp += xcm[i] * xcm[i] * ci + 2.0 * xcm[i] * (dix + uix);
            xyqdp += xcm[i] * ycm[i] * ci + xcm[i] * (diy + uiy) + ycm[i] * (dix + uix);
            xzqdp += xcm[i] * zcm[i] * ci + xcm[i] * (diz + uiz) + zcm[i] * (dix + uix);
            yxqdp += ycm[i] * xcm[i] * ci + ycm[i] * (dix + uix) + xcm[i] * (diy + uiy);
            yyqdp += ycm[i] * ycm[i] * ci + 2.0 * ycm[i] * (diy + uiy);
            yzqdp += ycm[i] * zcm[i] * ci + ycm[i] * (diz + uiz) + zcm[i] * (diy + uiy);
            zxqdp += zcm[i] * xcm[i] * ci + zcm[i] * (dix + uix) + xcm[i] * (diz + uiz);
            zyqdp += zcm[i] * ycm[i] * ci + zcm[i] * (diy + uiy) + ycm[i] * (diz + uiz);
            zzqdp += zcm[i] * zcm[i] * ci + 2.0 * zcm[i] * (diz + uiz);
        }

        // Convert the quadrupole from traced to traceless form.
        var qave = (xxqdp + yyqdp + zzqdp) / 3.0;
        xxqdp = 1.5 * (xxqdp - qave);
        xyqdp = 1.5 * xyqdp;
        xzqdp = 1.5 * xzqdp;
        yxqdp = 1.5 * yxqdp;
        yyqdp = 1.5 * (yyqdp - qave);
        yzqdp = 1.5 * yzqdp;
        zxqdp = 1.5 * zxqdp;
        zyqdp = 1.5 * zyqdp;
        zzqdp = 1.5 * (zzqdp - qave);

        // Add the traceless atomic quadrupoles to total quadrupole.
        for (Atom atom : activeAtoms) {
            int i = atom.getIndex() - 1;
            double[] globalMultipolei = globalMultipole[0][i];
            var qixx = globalMultipolei[t200];
            var qiyy = globalMultipolei[t020];
            var qizz = globalMultipolei[t002];
            var qixy = globalMultipolei[t110];
            var qixz = globalMultipolei[t101];
            var qiyz = globalMultipolei[t011];
            xxqdp += qixx;
            xyqdp += qixy;
            xzqdp += qixz;
            yxqdp += qixy;
            yyqdp += qiyy;
            yzqdp += qiyz;
            zxqdp += qixz;
            zyqdp += qiyz;
            zzqdp += qizz;
        }

        // Convert dipole to Debye and quadrupole to Buckingham.
        xdpl = xdpl * ELEC_ANG_TO_DEBYE;
        ydpl = ydpl * ELEC_ANG_TO_DEBYE;
        zdpl = zdpl * ELEC_ANG_TO_DEBYE;
        xxqdp = xxqdp * ELEC_ANG_TO_DEBYE;
        xyqdp = xyqdp * ELEC_ANG_TO_DEBYE;
        xzqdp = xzqdp * ELEC_ANG_TO_DEBYE;
        yxqdp = yxqdp * ELEC_ANG_TO_DEBYE;
        yyqdp = yyqdp * ELEC_ANG_TO_DEBYE;
        yzqdp = yzqdp * ELEC_ANG_TO_DEBYE;
        zxqdp = zxqdp * ELEC_ANG_TO_DEBYE;
        zyqdp = zyqdp * ELEC_ANG_TO_DEBYE;
        zzqdp = zzqdp * ELEC_ANG_TO_DEBYE;

        // Get dipole magnitude and diagonalize quadrupole tensor.
        netdpl = sqrt(xdpl * xdpl + ydpl * ydpl + zdpl * zdpl);
        double[][] a = new double[3][3];
        a[0][0] = xxqdp;
        a[0][1] = xyqdp;
        a[0][2] = xzqdp;
        a[1][0] = yxqdp;
        a[1][1] = yyqdp;
        a[1][2] = yzqdp;
        a[2][0] = zxqdp;
        a[2][1] = zyqdp;
        a[2][2] = zzqdp;
        EigenDecomposition e = new EigenDecomposition(new Array2DRowRealMatrix(a));
        // Eigenvalues are returned in descending order, but logged below in ascending order.
        var netqdp = e.getRealEigenvalues();

        logger.info("\n Electric Moments\n");
        logger.info(format("  Total Electric Charge:    %13.5f Electrons\n", netchg));
        logger.info(format("  Dipole Moment Magnitude:  %13.5f Debye\n", netdpl));
        logger.info(format("  Dipole X,Y,Z-Components:  %13.5f %13.5f %13.5f\n", xdpl, ydpl, zdpl));
        logger.info(format("  Quadrupole Moment Tensor: %13.5f %13.5f %13.5f", xxqdp, xyqdp, xzqdp));
        logger.info(format("       (Buckinghams)        %13.5f %13.5f %13.5f", yxqdp, yyqdp, yzqdp));
        logger.info(format("                            %13.5f %13.5f %13.5f\n", zxqdp, zyqdp, zzqdp));
        logger.info(format("  Principal Axes Quadrupole %13.5f %13.5f %13.5f\n", netqdp[2], netqdp[1], netqdp[0]));
    }

    /**
     * <p>destroy.</p>
     *
     * @throws java.lang.Exception if any.
     */
    public abstract void destroy() throws Exception;

    /**
     * <p>energy.</p>
     *
     * @param gradient a boolean.
     * @param print    a boolean.
     * @return a double.
     */
    public abstract double energy(boolean gradient, boolean print);

    /**
     * <p>getAxisAtoms.</p>
     *
     * @return an array of {@link int} objects.
     */
    public abstract int[][] getAxisAtoms();

    /**
     * <p>getCavitationEnergy.</p>
     *
     * @return a double.
     */
    public abstract double getCavitationEnergy();

    /**
     * <p>Getter for the field <code>coordinates</code>.</p>
     *
     * @return an array of {@link double} objects.
     */
    public abstract double[][][] getCoordinates();

    /**
     * <p>getDispersionEnergy.</p>
     *
     * @return a double.
     */
    public abstract double getDispersionEnergy();

    /**
     * <p>getElecForm.</p>
     *
     * @return a {@link ffx.potential.parameters.ForceField.ELEC_FORM} object.
     */
    public abstract ELEC_FORM getElecForm();

    /**
     * <p>getEwaldCoefficient.</p>
     *
     * @return a double.
     */
    public abstract double getEwaldCoefficient();

    /**
     * <p>getEwaldCutoff.</p>
     *
     * @return a double.
     */
    public abstract double getEwaldCutoff();

    /**
     * <p>getGK.</p>
     *
     * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood} object.
     */
    public abstract GeneralizedKirkwood getGK();

    /**
     * <p>getDispersionEnergy.</p>
     *
     * @return a double.
     */
    public abstract double getGKEnergy();

    /**
     * <p>getGKInteractions.</p>
     *
     * @return a int.
     */
    public abstract int getGKInteractions();

    /**
     * <p>getIndRealEnergy.</p>
     *
     * @return a double.
     */
    public double getIndRealEnergy() {
        return inducedRealSpaceEnergy;
    }

    /**
     * <p>getIndRecipEnergy.</p>
     *
     * @return a double.
     */
    public double getIndRecipEnergy() {
        return inducedReciprocalEnergy;
    }

    /**
     * <p>getIndSelfEnergy.</p>
     *
     * @return a double.
     */
    public double getIndSelfEnergy() {
        return inducedSelfEnergy;
    }

    /**
     * <p>getInteractions.</p>
     *
     * @return a int.
     */
    public abstract int getInteractions();

    /**
     * <p>getName.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public abstract String getName();

    /**
     * <p>getPermRealEnergy.</p>
     *
     * @return a double.
     */
    public double getPermRealEnergy() {
        return permanentRealSpaceEnergy;
    }

    /**
     * <p>getPermRecipEnergy.</p>
     *
     * @return a double.
     */
    public double getPermRecipEnergy() {
        return permanentReciprocalEnergy;
    }

    /**
     * <p>getPermSelfEnergy.</p>
     *
     * @return a double.
     */
    public double getPermSelfEnergy() {
        return permanentSelfEnergy;
    }

    /**
     * <p>getPermanentEnergy.</p>
     *
     * @return a double.
     */
    public double getPermanentEnergy() {
        return permanentMultipoleEnergy;
    }

    /**
     * <p>getPolarEps.</p>
     *
     * @return a double.
     */
    public abstract double getPolarEps();

    /**
     * <p>getPolarization11.</p>
     *
     * @return an array of {@link int} objects.
     */
    public abstract int[][] getPolarization11();

    /**
     * <p>getPolarization12.</p>
     *
     * @return an array of {@link int} objects.
     */
    public abstract int[][] getPolarization12();

    /**
     * <p>getPolarization13.</p>
     *
     * @return an array of {@link int} objects.
     */
    public abstract int[][] getPolarization13();

    /**
     * <p>Getter for the field <code>polarizationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    /**
     * <p>getPolarizationType.</p>
     *
     * @return a {@link ffx.potential.nonbonded.ParticleMeshEwald.Polarization} object.
     */
    public abstract Polarization getPolarizationType();

    /**
     * <p>getReciprocalSpace.</p>
     *
     * @return a {@link ffx.potential.nonbonded.ReciprocalSpace} object.
     */
    public abstract ReciprocalSpace getReciprocalSpace();

    /**
     * <p>getScale14.</p>
     *
     * @return a double.
     */
    public abstract double getScale14();

    /**
     * Returns the SCF algorithm in use.
     *
     * @return The SCF algorithm used.
     */
    public SCFAlgorithm getScfAlgorithm() {
        return scfAlgorithm;
    }

    /**
     * <p>getGKEnergy.</p>
     *
     * @return a double.
     */
    public double getSolvationEnergy() {
        return generalizedKirkwoodEnergy;
    }

    /**
     * <p>Getter for the field <code>totalElectrostaticEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getTotalElectrostaticEnergy() {
        return permanentMultipoleEnergy + polarizationEnergy + generalizedKirkwoodEnergy;
    }

    /**
     * <p>Getter for the field <code>totalMultipoleEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getTotalMultipoleEnergy() {
        return permanentMultipoleEnergy + polarizationEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract double getd2EdL2();

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract double getdEdL();

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract void getdEdXdL(double[] gradients);

    /**
     * <p>Setter for the field <code>atoms</code>.</p>
     *
     * @param atoms    an array of {@link ffx.potential.bonded.Atom} objects.
     * @param molecule an array of {@link int} objects.
     */
    public abstract void setAtoms(Atom[] atoms, int[] molecule);

    /**
     * <p>setCrystal.</p>
     *
     * @param crystal a {@link ffx.crystal.Crystal} object.
     */
    public abstract void setCrystal(Crystal crystal);

    /**
     * <p>setFixedCharges.</p>
     *
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public abstract void setFixedCharges(Atom[] atoms);

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract void setLambda(double lambda);

    /**
     * <p>setLambdaMultipoleScale.</p>
     *
     * @param scale a double.
     */
    public abstract void setLambdaMultipoleScale(double scale);

    /**
     * <p>Setter for the field <code>polarization</code>.</p>
     *
     * @param set a {@link ffx.potential.nonbonded.ParticleMeshEwald.Polarization} object.
     */
    public void setPolarization(Polarization set) {
        this.polarization = set;
    }

    public enum Polarization {
        MUTUAL, DIRECT, NONE
    }

    public enum LambdaMode {
        OFF, CONDENSED, CONDENSED_NO_LIGAND, VAPOR
    }

    /**
     * Describes available SCF algorithms, and whether they are supported by the
     * FFX and/or CUDA implementations.
     */
    public enum SCFAlgorithm {
        SOR(true, true), CG(true, true), EPT(true, true);

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

    /**
     * <p>assignPolarizationGroups.</p>
     */
    protected void assignPolarizationGroups() {
        // Find directly connected group members for each atom.
        List<Integer> group = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.getIndex() - 1 != i) {
                logger.severe(" Atom indexing is not consistent in PME.");
            }
        }
        for (Atom ai : atoms) {
            group.clear();
            int index = ai.getIndex() - 1;
            group.add(index);
            PolarizeType polarizeType = ai.getPolarizeType();
            if (polarizeType != null) {
                if (polarizeType.polarizationGroup != null) {
                    growGroup(group, ai);
                    sort(group);
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
        // Find 1-2 group relationships.
        int[] mask = new int[nAtoms];
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
                List<Bond> bonds = aj.getBonds();
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
            sort(list);
            ip12[i] = new int[list.size()];
            int j = 0;
            for (int k : list) {
                ip12[i][j++] = k;
            }
        }

        // Find 1-3 group relationships.
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
            sort(list);
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
     * @param seed  The bonds of the seed atom are queried for inclusion in the
     *              group.
     */
    private void growGroup(List<Integer> group, Atom seed) {
        List<Bond> bonds = seed.getBonds();
        for (Bond bond : bonds) {
            Atom atom = bond.get1_2(seed);
            int tj = atom.getType();
            boolean added = false;
            PolarizeType polarizeType = seed.getPolarizeType();
            for (int type : polarizeType.polarizationGroup) {
                if (type == tj) {
                    Integer index = atom.getIndex() - 1;
                    if (!group.contains(index)) {
                        group.add(index);
                        added = true;
                        break;
                    }
                }
            }
            if (added) {
                growGroup(group, atom);
            }
        }
    }

}
