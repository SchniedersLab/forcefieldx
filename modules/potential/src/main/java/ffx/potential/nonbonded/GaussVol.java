//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.rsq;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;

/**
 * A class that implements the Gaussian description of an object (molecule) made of a overlapping spheres.
 * <p>
 * Ported from C++ code by Emilio Gallicchio.
 * GaussVol is part of the AGBNP/OpenMM implicit solvent model.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GaussVol {

    /* -------------------------------------------------------------------------- *
     *                                 GaussVol                                   *
     * -------------------------------------------------------------------------- *
     * This file is part of the AGBNP/OpenMM implicit solvent model software      *
     * implementation funded by the National Science Foundation under grant:      *
     * NSF SI2 1440665  "SI2-SSE: High-Performance Software for Large-Scale       *
     * Modeling of Binding Equilibria"                                            *
     *                                                                            *
     * copyright (c) 2016 Emilio Gallicchio                                       *
     * Authors: Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>                 *
     * Contributors:                                                              *
     *                                                                            *
     *  AGBNP/OpenMM is free software: you can redistribute it and/or modify      *
     *  it under the terms of the GNU Lesser General Public License version 3     *
     *  as published by the Free Software Foundation.                             *
     *                                                                            *
     *  AGBNP/OpenMM is distributed in the hope that it will be useful,           *
     *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
     *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
     *  GNU General Public License for more details.                              *
     *                                                                            *
     *  You should have received a copy of the GNU General Public License         *
     *  along with AGBNP/OpenMM.  If not, see <http://www.gnu.org/licenses/>      *
     *                                                                            *
     * -------------------------------------------------------------------------- */

    private static final Logger logger = Logger.getLogger(GaussVol.class.getName());

    /**
     * Conversion factors from spheres to Gaussians.
     */
    private static double KFC = 2.2269859253;
    private static double PFC = 2.5;
    /**
     * Minimum volume.
     */
    private static double MIN_GVOL = Double.MIN_VALUE;
    /**
     * Maximum overlap level.
     */
    private static int MAX_ORDER = 8;
    /**
     * Finite-Difference step size to compute surface area.
     */
    private static final double offset = 0.00005;

    // Volume cutoffs in switching function in nanometers.
    // private static double ANG3 = 0.001;
    // private static double VOLMINA = 0.01 * ANG3;
    // private static double VOLMINB = 0.1 * ANG3;
    /**
     * Volume cutoffs for switching function in Angstroms.
     */
    private static double ANG3 = 1.0;
    private static double VOLMINA = 0.01 * ANG3;
    private static double VOLMINB = 0.1 * ANG3;

    /**
     * Number of atoms.
     */
    private int nAtoms;
    /**
     * Atomic radii for all atoms. To approximate solvent excluded volume, or solvent accessible surface area,
     * a probe radius can be added to each radius.
     */
    private double[] radii;
    /**
     * All atomic radii with a small offset added, which are used to compute surface area using a
     * finite-difference approach.
     */
    private double[] radiiOffset;
    /**
     * Atomic "self" volumes for all atoms, which are computed from atomic radii.
     */
    private double[] volumes;
    /**
     * Atomic "self" volumes for all atoms, which are computed from the radii plus offset array.
     */
    private double[] volumeOffset;
    /**
     * Surface tension -- in our implementation these values are kept at 1.0.
     */
    private double[] gammas;
    /**
     * Flag to denote if an atom is a hydrogen, which results in it not contributing to the volume or S.A.
     */
    private boolean[] ishydrogen;
    /**
     * The Gaussian Overlap Tree.
     */
    private GaussianOverlapTree tree;
    /**
     * Surface area (Ang^2).
     */
    private double surfaceArea;
    /**
     * Surface area energy (kcal/mol).
     */
    private double surfaceAreaEnergy;
    /**
     * Volume (Ang^3).
     */
    private double volume;
    /**
     * Volume energy (kcal/mol).
     */
    private double volumeEnergy;
    /**
     * Cavitation energy, which is a function of volume and/or surface area.
     */
    private double cavitationEnergy;
    /**
     * The volume gradient, which does not include the solvent pressure constant (i.e. this is in Ang^2).
     */
    private double[] volumeGradient;
    /**
     * The surface area gradient, which does not include the surface tension constant (i.e. this is in Ang).
     */
    private double[] surfaceAreaGradient;
    /**
     * Solvent pressure in kcal/mol/Ang^3.
     */
    private static final double solventPressure = 0.0327;
    /**
     * Surface tension in kcal/mol/Ang^2.
     */
    private static final double surfaceTension = 0.08;
    /**
     * Radius where volume dependence crosses over to surface area dependence (approximately at 1 nm).
     */
    private double crossOver = 3.0 * surfaceTension / solventPressure;
    /**
     * Begin turning off the Volume term.
     */
    private double volumeOff = crossOver - 3.5;
    /**
     * Volume term is zero at the cut-off.
     */
    private double volumeCut = crossOver + 3.5;
    /**
     * Begin turning off the SA term.
     */
    private double surfaceAreaOff = crossOver + 3.9;
    /**
     * SA term is zero at the cut-off.
     */
    private double surfaceAreaCut = crossOver - 3.5;
    /**
     * Volume multiplicative switch.
     */
    private MultiplicativeSwitch volumeSwitch = new MultiplicativeSwitch(volumeCut, volumeOff);
    /**
     * Surface area multiplicative switch.
     */
    private MultiplicativeSwitch surfaceAreaSwitch = new MultiplicativeSwitch(surfaceAreaCut, surfaceAreaOff);

    /**
     * Creates/Initializes a GaussVol instance.
     *
     * @param nAtoms     The number of atoms.
     * @param ishydrogen If each atom is a hydrogen or not.
     */
    public GaussVol(int nAtoms, boolean[] ishydrogen) {
        tree = new GaussianOverlapTree(nAtoms);
        this.nAtoms = nAtoms;
        this.ishydrogen = ishydrogen;

        radii = new double[nAtoms];
        radiiOffset = new double[nAtoms];
        fill(radii, 1.0);
        fill(radiiOffset, 1.0 + offset);

        volumes = new double[nAtoms];
        volumeOffset = new double[nAtoms];
        fill(volumes, 0.0);
        fill(volumeOffset, 0.0);

        gammas = new double[nAtoms];
        fill(gammas, 0.0);
    }

    /**
     * Creates/Initializes a GaussVol instance.
     *
     * @param nAtoms     The number of atoms.
     * @param radii      Atomic radii.
     * @param volumes    Atomic volumes.
     * @param gammas     Atomic surface tensions.
     * @param ishydrogen True if the atom is a hydrogen.
     */
    public GaussVol(int nAtoms, double[] radii, double[] volumes,
                    double[] gammas, boolean[] ishydrogen) {
        tree = new GaussianOverlapTree(nAtoms);
        this.nAtoms = nAtoms;
        this.radii = radii;
        this.volumes = volumes;
        this.gammas = gammas;
        this.ishydrogen = ishydrogen;

        radiiOffset = new double[nAtoms];
        volumeOffset = new double[nAtoms];
        double fourThirdsPI = 4.0 / 3.0 * PI;
        for (int i = 0; i < nAtoms; i++) {
            radiiOffset[i] = radii[i] + offset;
            volumeOffset[i] = fourThirdsPI * pow(radiiOffset[i], 3);
        }
    }

    /**
     * Return the cavitation energy.
     *
     * @return
     */
    public double getEnergy() {
        return cavitationEnergy;
    }

    public double getVolumeEnergy() {
        return volumeEnergy;
    }

    public double getSurfaceAreaEnergy() {
        return surfaceAreaEnergy;
    }

    public double getVolume() {
        return volume;
    }

    public double getSurfaceArea() {
        return surfaceArea;
    }

    /**
     * Set the isHydrogen flag.
     *
     * @param isHydrogen
     * @return
     * @throws Exception
     */
    int setIsHydrogen(boolean[] isHydrogen) throws Exception {
        if (nAtoms == isHydrogen.length) {
            this.ishydrogen = isHydrogen;
            return nAtoms;
        } else {
            throw new Exception(" setIsHydrogen: number of atoms does not match");
        }
    }

    /**
     * Set radii.
     *
     * @param radii Atomic radii (Angstroms).
     * @return
     * @throws Exception
     */
    int setRadii(double[] radii) throws Exception {
        if (nAtoms == radii.length) {
            this.radii = radii;
            radiiOffset = new double[nAtoms];
            volumeOffset = new double[nAtoms];
            double fourThirdsPI = 4.0 / 3.0 * PI;
            for (int i = 0; i < nAtoms; i++) {
                radiiOffset[i] = radii[i] + offset;
                volumeOffset[i] = fourThirdsPI * pow(radiiOffset[i], 3);
            }
            return nAtoms;
        } else {
            throw new Exception(" setRadii: number of atoms does not match");
        }
    }

    /**
     * Set volumes.
     *
     * @param volumes Atomic volumes (Angstroms^3).
     * @return
     * @throws Exception
     */
    int setVolumes(double[] volumes) throws Exception {
        if (nAtoms == volumes.length) {
            this.volumes = volumes;
            return nAtoms;
        } else {
            throw new Exception(" setVolumes: number of atoms does not match");
        }
    }

    /**
     * Set gamma values.
     *
     * @param gammas Gamma values (kcal/mol/A^2).
     * @return
     * @throws Exception
     */
    int setGammas(double[] gammas) throws Exception {
        if (nAtoms == gammas.length) {
            this.gammas = gammas;
            return nAtoms;
        } else {
            throw new Exception(" setGammas: number of atoms does not match");
        }
    }

    /**
     * Compute molecular volume and surface area.
     *
     * @param positions Atomic positions to use.
     */
    public double computeVolumeAndSA(double[][] positions) {
        return energyAndGradient(positions, new double[3][positions.length]);
    }

    /**
     * Compute molecular volume and surface area.
     *
     * @param positions Atomic positions to use.
     */
    public double energyAndGradient(double[][] positions, double[][] gradient) {

        // Output
        double[][] grad = new double[nAtoms][3];
        double[] totalVolume = new double[1];
        double[] energy = new double[1];
        double[] gradV = new double[nAtoms];
        double[] freeVolume = new double[nAtoms];
        double[] selfVolume = new double[nAtoms];

        if (volumeGradient == null || volumeGradient.length != nAtoms * 3) {
            volumeGradient = new double[nAtoms * 3];
            surfaceAreaGradient = new double[nAtoms * 3];
        }

        computeTree(positions);
        computeVolume(positions, totalVolume, energy, grad, gradV, freeVolume, selfVolume);

        double selfVolumeSum = 0;
        for (int i = 0; i < nAtoms; i++) {
            selfVolumeSum += selfVolume[i];

            // Volume based gradient
            int index = i * 3;
            volumeGradient[index++] = grad[i][0];
            volumeGradient[index++] = grad[i][1];
            volumeGradient[index] = grad[i][2];

            // Surface Area based gradient
            index = i * 3;
            surfaceAreaGradient[index++] = grad[i][0];
            surfaceAreaGradient[index++] = grad[i][1];
            surfaceAreaGradient[index] = grad[i][2];
        }

        // Save the total molecular volume.
        volume = totalVolume[0];

        // Calculate a purely volume based cavitation energy.
        volumeEnergy = volume * solventPressure;

        // Save a reference to non-offset radii and volumes.
        double[] radiiBak = this.radii;
        double[] volumeBak = this.volumes;

        // Load the offset radii and volumes.
        this.radii = radiiOffset;
        this.volumes = volumeOffset;

        // Run Volume calculation on radii that are slightly offset in order to do finite difference to get back surface area
        rescanTreeVolumes(positions);
        computeVolume(positions, totalVolume, energy, grad, gradV, freeVolume, selfVolume);

        // Set the radii and volumes to their non-offset values.
        this.radii = radiiBak;
        this.volumes = volumeBak;

        double selfVolumeOffsetSum = 0;
        for (int i = 0; i < nAtoms; i++) {
            selfVolumeOffsetSum += selfVolume[i];

            // Surface Area based gradient
            int index = i * 3;
            surfaceAreaGradient[index] = (grad[i][0] - surfaceAreaGradient[index++]) / offset;
            surfaceAreaGradient[index] = (grad[i][1] - surfaceAreaGradient[index++]) / offset;
            surfaceAreaGradient[index] = (grad[i][2] - surfaceAreaGradient[index]) / offset;
        }

        // Calculate the surface area.
        surfaceArea = (selfVolumeOffsetSum - selfVolumeSum) / offset;

        // Calculate a purely surface area based cavitation energy.
        surfaceAreaEnergy = surfaceArea * surfaceTension;

        // Use SA to find an effective cavity radius.
        double reff = 0.5 * sqrt(surfaceArea / PI);
        double reff2 = reff * reff;
        double reff3 = reff2 * reff;
        double reff4 = reff3 * reff;
        double reff5 = reff4 * reff;
        double dreff = reff / (2.0 * surfaceArea);

        // Find the cavitation energy using a combination of volume and surface area dependence.
        if (reff < volumeOff) {
            // Volume only.
            // logger.info(format(" Volume Only: %8.6f < %8.6f", reff, volumeOff));
            // Find cavity energy from only the molecular volume.
            cavitationEnergy = volumeEnergy;
            addVolumeGradient(1.0, 0.0, gradient);
        } else if (reff <= volumeCut) {
            // Include a tapered molecular volume.
            double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
            double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dreff;
            // logger.info(format(" Tapered Volume: %8.6f <= %8.6f <= %8.6f -> %8.6f", volumeOff, reff, volumeCut, taper));
            cavitationEnergy = taper * volumeEnergy;
            addVolumeGradient(taper, dtaper, gradient);
        }

        if (reff > surfaceAreaOff) {
            // Find cavity energy from only SA.
            cavitationEnergy = surfaceAreaEnergy;
            // logger.info(format(" SA Only: %8.6f > %8.6f", surfaceAreaCut, reff));
            addSurfaceAreaGradient(1.0, 0.0, gradient);
        } else if (reff >= surfaceAreaCut) {
            // Include a tapered surface area term.
            double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
            double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dreff;
            // logger.info(format(" Tapered SA: %8.6f <= %8.6f <= %8.6f -> %8.6f", surfaceAreaCut, reff, surfaceAreaOff, taperSA));
            cavitationEnergy += taperSA * surfaceAreaEnergy;
            addSurfaceAreaGradient(taperSA, dtaperSA, gradient);
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("\n Volume:              %8.3f (Ang^3)", totalVolume[0]));
            logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
            logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
            logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
            logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
        }

        return cavitationEnergy;
    }

    /**
     * Collect the volume base cavitation energy and gradient contributions.
     *
     * @param taperVolume
     * @param dTaperVolumedR
     * @param gradient
     */
    private void addVolumeGradient(double taperVolume, double dTaperVolumedR, double[][] gradient) {
        taperVolume *= solventPressure;
        dTaperVolumedR *= volumeEnergy;
        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            gradient[0][i] += taperVolume * volumeGradient[index] + dTaperVolumedR * surfaceAreaGradient[index++];
            gradient[1][i] += taperVolume * volumeGradient[index] + dTaperVolumedR * surfaceAreaGradient[index++];
            gradient[2][i] += taperVolume * volumeGradient[index] + dTaperVolumedR * surfaceAreaGradient[index];
        }
    }

    /**
     * Collect the surface area based cavitation energy and gradient contributions.
     *
     * @param taperSA
     * @param dTaperSAdR
     * @param gradient
     */
    private void addSurfaceAreaGradient(double taperSA, double dTaperSAdR, double[][] gradient) {
        double factor = surfaceTension * taperSA + surfaceAreaEnergy * dTaperSAdR;
        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            gradient[0][i] += factor * surfaceAreaGradient[index++];
            gradient[1][i] += factor * surfaceAreaGradient[index++];
            gradient[2][i] += factor * surfaceAreaGradient[index];
        }
    }

    /**
     * Constructs the tree.
     *
     * @param positions Current atomic positions.
     */
    private void computeTree(double[][] positions) {
        tree.computeOverlapTreeR(positions, radii, volumes, gammas, ishydrogen);
    }

    /**
     * Returns GaussVol volume energy function and forces.
     * Also returns gradients with respect to atomic volumes and atomic free-volumes and self-volumes.
     *
     * @param positions
     * @param totalVolume
     * @param totalEnergy
     * @param grad
     * @param gradV
     * @param free_volume
     * @param self_volume
     */
    private void computeVolume(double[][] positions,
                               double[] totalVolume, double[] totalEnergy,
                               double[][] grad, double[] gradV,
                               double[] free_volume, double[] self_volume) {
        tree.computeVolume2R(positions, totalVolume, totalEnergy, grad, gradV, free_volume, self_volume);
        for (int i = 0; i < nAtoms; ++i) {
            if (volumes[i] > 0) {
                gradV[i] = gradV[i] / volumes[i];
            }
        }

    }

    /**
     * Rescan the tree after resetting gammas, radii and volumes
     *
     * @param positions
     */
    private void rescanTreeVolumes(double[][] positions) {
        tree.rescanTreeV(positions, radii, volumes, gammas, ishydrogen);
    }

    /**
     * Rescan the tree resetting gammas only with current values.
     */
    private void rescanTreeGammas() {
        tree.rescanTreeG(gammas);
    }

    /**
     * Returns number of overlaps for each atom.
     *
     * @param nov Number of overlaps.
     */
    private void getStats(int[] nov) {
        if (nov.length != nAtoms) return;

        for (int i = 0; i < nAtoms; i++) nov[i] = 0;
        for (int atom = 0; atom < nAtoms; atom++) {
            int slot = atom + 1;
            nov[atom] = tree.nchildrenUnderSlotR(slot);
        }

    }

    /**
     * Print the tree.
     */
    private void printTree() {
        tree.printTree();
    }

    /**
     * 3D Gaussian, V,c,a representation.
     */
    private class GaussianVca {
        // Gaussian volume
        public double v;
        // Gaussian exponent
        public double a;
        // Center
        public double[] c = new double[3];
    }

    /**
     * Overlap between two Gaussians represented by a (V,c,a) triplet
     * V: volume of Gaussian
     * c: position of Gaussian
     * a: exponential coefficient
     * g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
     * this version is based on V=V(V1,V2,r1,r2,alpha)
     * alpha = (a1 + a2)/(a1 a2)
     * dVdr is (1/r)*(dV12/dr)
     * dVdV is dV12/dV1
     * dVdalpha is dV12/dalpha
     * d2Vdalphadr is (1/r)*d^2V12/dalpha dr
     * d2VdVdr is (1/r) d^2V12/dV1 dr
     */
    private class GaussianOverlap implements Comparable<GaussianOverlap> {
        /**
         * level (0=root, 1=atoms, 2=2-body, 3=3-body, etc.)
         */
        public int level;
        /**
         * Gaussian representing overlap
         */
        public GaussianVca g = new GaussianVca();
        /**
         * Volume of overlap (also stores Psi1..i in GPU version)
         */
        public double volume;
        /**
         * Derivative wrt volume of first atom (also stores F1..i in GPU version)
         */
        double dvv1;
        /**
         * Derivative wrt position of first atom (also stores P1..i in GPU version)
         */
        double[] dv1 = new double[3];
        ;
        /**
         * Sum gammai for this overlap
         */
        double gamma1i;
        /**
         * Self volume accumulator (also stores Psi'1..i in GPU version)
         */
        double selfVolume;
        /**
         * Switching function derivatives
         */
        double sfp;
        /**
         * The atomic index of the last atom of the overlap list (i, j, k, ..., atom)
         */
        public int atom;
        /**
         * = (Parent, atom)
         * index in tree list of parent overlap
         */
        int parentIndex;
        /**
         * Start index in tree array of children
         */
        int childrenStartindex;
        /**
         * Number of children.
         */
        int childrenCount;

        /**
         * Print overlaps.
         */
        public void printOverlap() {
            logger.info(format(" Gaussian Overlap %d: Atom: %d, Parent: %d, ChildrenStartIndex: %d, ChildrenCount: %d," +
                            "Volume: %6.3f, Gamma: %6.3f, Gauss.a: %6.3f, Gauss.v: %6.3f, Gauss.center (%6.3f,%6.3f,%6.3f)," +
                            "dedx: %6.3f, dedy: %6.3f, dedz: %6.3f, sfp: %6.3f",
                    level, atom, parentIndex, childrenStartindex, childrenCount, volume, gamma1i, g.a, g.v,
                    g.c[0], g.c[1], g.c[2], dv1[0], dv1[1], dv1[2], sfp));
        }

        /**
         * Order by volume, larger first.
         *
         * @param o GaussianOverlap to compare to.
         * @return Compare by volume.
         */
        @Override
        public int compareTo(GaussianOverlap o) {
            return Double.compare(volume, o.volume);
        }
    }

    /**
     * Gaussian Overlap Tree.
     */
    private class GaussianOverlapTree {

        /**
         * Number of atoms.
         */
        int nAtoms;
        /**
         * The root is at index 0, atoms are at 1..natoms+1
         */
        List<GaussianOverlap> overlaps;

        GaussianOverlapTree(int nAtoms) {
            this.nAtoms = nAtoms;
        }


        /**
         * @param pos        Atomic positions.
         * @param radii      Atomic radii.
         * @param volumes    Atomic volumes.
         * @param gammas     Atomic surface tensions.
         * @param ishydrogen True if the atom is a hydrogen.
         * @return
         */
        int initOverlapTree(double[][] pos, double[] radii, double[] volumes,
                            double[] gammas, boolean[] ishydrogen) {


            // Reset tree
            overlaps = new ArrayList<>();

            // Slot 0 contains the master tree information, children = all of the atoms.
            GaussianOverlap overlap = new GaussianOverlap();
            overlap.level = 0;
            overlap.volume = 0;
            //overlap.dv1 = new double[3];
            overlap.dvv1 = 0.;
            overlap.selfVolume = 0;
            overlap.sfp = 1.;
            overlap.gamma1i = 0.;
            overlap.parentIndex = -1;
            overlap.atom = -1;
            overlap.childrenStartindex = 1;
            overlap.childrenCount = nAtoms;

            overlaps.add(overlap);

            //logger.info(String.format("Initializing Tree Using Following Inputs"));
            // list of atoms start at slot #1
            for (int iat = 0; iat < nAtoms; iat++) {
                overlap = new GaussianOverlap();
                double a = KFC / (radii[iat] * radii[iat]);
                double vol = ishydrogen[iat] ? 0. : volumes[iat];
                overlap.level = 1;
                overlap.g.v = vol;
                overlap.g.a = a;
                overlap.g.c = pos[iat];
                overlap.volume = vol;
                //logger.info(String.format("Atom: %s, Gaussian Volume: %8.6f, Center: (%8.6f,%8.6f,%8.6f)", iat, overlap.g.v,
                //overlap.g.c[0],overlap.g.c[1],overlap.g.c[2]));

                //overlap.dv1 = new double[3];
                overlap.dvv1 = 1.; //dVi/dVi
                overlap.selfVolume = 0.;
                overlap.sfp = 1.;
                overlap.gamma1i = gammas[iat];// gamma[iat]/SA_DR;
                overlap.parentIndex = 0;
                overlap.atom = iat;
                overlap.childrenStartindex = -1;
                overlap.childrenCount = -1;
                overlaps.add(overlap);
            }

            return 1;
        }

        /**
         * @param parent_index      Parent index.
         * @param children_overlaps Children overlaps.
         * @return
         */
        int addChildren(int parent_index, ArrayList<GaussianOverlap> children_overlaps) {
            int i, slot;

            /* adds children starting at the last slot */
            int start_index = overlaps.size();

            int noverlaps = children_overlaps.size();

            /* retrieves address of root overlap */
            GaussianOverlap root = overlaps.get(parent_index);

            /* registers list of children */
            root.childrenStartindex = start_index;
            root.childrenCount = noverlaps;

            // Sort neighbors by overlap volume.
            Collections.sort(children_overlaps);

            int root_level = root.level;

            // Now copies the children overlaps from temp buffer.
            for (GaussianOverlap child : children_overlaps) {
                child.level = root_level + 1;

                // Connect overlap to parent
                child.parentIndex = parent_index;

                // Reset its children indexes
                child.childrenStartindex = -1;
                child.childrenCount = -1;

                // Add to tree.
                // note that the 'root' pointer may be invalidated by the push back below
                overlaps.add(child);
            }

            return start_index;
        }

        /**
         * scans the siblings of overlap identified by "root_index" to create children overlaps,
         * returns them into the "children_overlaps" buffer: (root) + (atom) -> (root, atom)
         *
         * @param root_index        Root index.
         * @param children_overlaps Children overlaps.
         * @return
         */
        int computeChildren(int root_index, ArrayList<GaussianOverlap> children_overlaps) {
            int parent_index;
            int sibling_start, sibling_count;
            int j;

            // Reset output buffer
            children_overlaps.clear();

            // Retrieves overlap.
            GaussianOverlap root = overlaps.get(root_index);

            // Retrieves parent overlap.
            parent_index = root.parentIndex;

            //master root? can't do computeChildren() on master root
            if (parent_index < 0) return 1;

            if (root.level >= MAX_ORDER) return 1;

            GaussianOverlap parent = overlaps.get(parent_index);

            /* Retrieves start index and count of siblings */
            sibling_start = parent.childrenStartindex;
            sibling_count = parent.childrenCount;

            // Parent is not initialized?
            if (sibling_start < 0 || sibling_count < 0) return -1;

            // This overlap somehow is not the child of registered parent.
            if (root_index < sibling_start && root_index > sibling_start + sibling_count - 1) return -1;

            // Now loops over "younger" siblings (i<j loop) to compute new overlaps.
            for (int slotj = root_index + 1; slotj < sibling_start + sibling_count; slotj++) {

                GaussianVca g12 = new GaussianVca();

                GaussianOverlap sibling = overlaps.get(slotj);
                double gvol;
                double[] dVdr = new double[1];
                double[] dVdV = new double[1];
                double[] sfp = new double[1];

                // Atomic gaussian of last atom of sibling.

                int atom2 = sibling.atom;
                GaussianVca g1 = root.g;
                // Atoms are stored in the tree at indexes 1...N
                GaussianVca g2 = overlaps.get(atom2 + 1).g;
                gvol = overlapGaussianAlpha(g1, g2, g12, dVdr, dVdV, sfp);

                /* create child if overlap volume is not zero */
                if (gvol > MIN_GVOL) {
                    GaussianOverlap ov = new GaussianOverlap();
                    ov.g = g12;
                    ov.volume = gvol;
                    ov.selfVolume = 0;
                    ov.atom = atom2;

                    // dv1 is the gradient of V(123..)n with respect to the position of 1
                    // ov.dv1 = ( g2.c - g1.c ) * (-dVdr);
                    diff(g2.c, g1.c, ov.dv1);
                    scalar(ov.dv1, -dVdr[0], ov.dv1);

                    //dvv1 is the derivative of V(123...)n with respect to V(123...)
                    ov.dvv1 = dVdV[0];
                    ov.sfp = sfp[0];
                    ov.gamma1i = root.gamma1i + overlaps.get(atom2 + 1).gamma1i;
                    children_overlaps.add(ov);
                }
            }
            return 1;
        }

        /**
         * Grow the tree with more children starting at the given root slot (recursive).
         *
         * @param root The root index.
         * @return
         */
        int computeAndAddChildrenR(int root) {
            ArrayList<GaussianOverlap> children_overlaps = new ArrayList<>();
            computeChildren(root, children_overlaps);
            int noverlaps = children_overlaps.size();
            if (noverlaps > 0) {
                int start_slot = addChildren(root, children_overlaps);
                for (int ichild = start_slot; ichild < start_slot + noverlaps; ichild++) {
                    computeAndAddChildrenR(ichild);
                }
            }
            return 1;
        }

        /**
         * @param pos        Atomic positions.
         * @param radii      Atomic radii.
         * @param volumes    Atomic volumes.
         * @param gammas     Atomic surface tensions.
         * @param ishydrogen True if the atom is a hydrogen.
         * @return
         */
        int computeOverlapTreeR(double[][] pos, double[] radii,
                                double[] volumes, double[] gammas, boolean[] ishydrogen) {
            initOverlapTree(pos, radii, volumes, gammas, ishydrogen);
            for (int slot = 1; slot <= nAtoms; slot++) {
                computeAndAddChildrenR(slot);
            }
            return 1;
        }

        /**
         * Compute volumes, energy of the overlap at slot and calls itself recursively to get
         * the volumes of the children.
         *
         * @param slot
         * @param psi1i       Subtree accumulator for free volume.
         * @param f1i         Subtree accumulator for free volume.
         * @param p1i         Subtree accumulator for free volume.
         * @param psip1i      Subtree accumulator for self volume.
         * @param fp1i        Subtree accumulators for self volume.
         * @param pp1i        Subtree accumulators for self volume.
         * @param energy1i    Subtree accumulator for volume-based energy.
         * @param fenergy1i   Subtree accumulator for volume-based energy.
         * @param penergy1i   subtree accumulator for volume-based energy.
         * @param dr          Gradient of volume-based energy wrt to atomic positions.
         * @param dv          Gradient of volume-based energy wrt to atomic volumes.
         * @param free_volume Atomic free volumes.
         * @param self_volume Atomic self volumes.
         * @return
         */
        int computeVolumeUnderSlot2R(
                int slot,
                double[] psi1i, double[] f1i, double[] p1i,
                double[] psip1i, double[] fp1i, double[] pp1i,
                double[] energy1i, double[] fenergy1i, double[] penergy1i,
                double[][] dr, double[] dv, double[] free_volume, double[] self_volume) {

            GaussianOverlap ov = overlaps.get(slot);
            double cf = ov.level % 2 == 0 ? -1.0 : 1.0;
            double volcoeff = ov.level > 0 ? cf : 0;
            double volcoeffp = ov.level > 0 ? volcoeff / (double) ov.level : 0;

            int atom = ov.atom;
            double ai = overlaps.get(atom + 1).g.a;
            double a1i = ov.g.a;
            double a1 = a1i - ai;

            psi1i[0] = volcoeff * ov.volume; //for free volumes
            //logger.info(format("Atom: %d, psi1i: %8.6f",ov.atom, psi1i[0]));
            f1i[0] = volcoeff * ov.sfp;

            psip1i[0] = volcoeffp * ov.volume; //for self volumes
            fp1i[0] = volcoeffp * ov.sfp;

            energy1i[0] = volcoeffp * ov.gamma1i * ov.volume; //EV energy
            fenergy1i[0] = volcoeffp * ov.sfp * ov.gamma1i;

            // These arrays must be allocated prior to entering the method to return their values.
            // p1i[0] = new double[3];
            // pp1i[0] = new double[3];
            // penergy1i[0] = new double[3];

            if (ov.childrenStartindex >= 0) {
                for (int sloti = ov.childrenStartindex; sloti < ov.childrenStartindex + ov.childrenCount; sloti++) {
                    double[] psi1it = new double[1];
                    double[] f1it = new double[1];
                    double[] p1it = new double[3];
                    double[] psip1it = new double[1];
                    double[] fp1it = new double[1];
                    double[] pp1it = new double[3];
                    double[] energy1it = new double[1];
                    double[] fenergy1it = new double[1];
                    double[] penergy1it = new double[3];
                    computeVolumeUnderSlot2R(sloti,
                            psi1it, f1it, p1it,
                            psip1it, fp1it, pp1it,
                            energy1it, fenergy1it, penergy1it,
                            dr, dv, free_volume, self_volume);
                    psi1i[0] += psi1it[0];
                    f1i[0] += f1it[0];
                    sum(p1i, p1it, p1i);

                    psip1i[0] += psip1it[0];
                    fp1i[0] += fp1it[0];
                    sum(pp1i, pp1it, pp1i);

                    energy1i[0] += energy1it[0];
                    fenergy1i[0] += fenergy1it[0];
                    sum(penergy1i, penergy1it, penergy1i);
                }
            }

            if (ov.level > 0) {
                // Contributions to free and self volume of last atom
                free_volume[atom] += psi1i[0];
                self_volume[atom] += psip1i[0];

                // Contributions to energy gradients
                double c2 = ai / a1i;

                // dr[atom] += (-ov.dv1) * fenergy1i + penergy1i * c2;
                double[] work1 = new double[3];
                double[] work2 = new double[3];
                double[] work3 = new double[3];
                scalar(penergy1i, c2, work1);
                scalar(ov.dv1, -fenergy1i[0], work2);
                sum(work1, work2, work3);
                sum(dr[atom], work3, dr[atom]);

                // ov.g.v is the unswitched volume
                // dv[atom] += ov.g.v * fenergy1i; //will be divided by Vatom later
                dv[atom] += ov.g.v * fenergy1i[0];

                // Update subtree P1..i's for parent
                c2 = a1 / a1i;

                // p1i = (ov.dv1) * f1i + p1i * c2;
                scalar(ov.dv1, f1i[0], work1);
                scalar(pp1i, c2, work2);
                sum(work1, work2, p1i);

                // pp1i = (ov.dv1) * fp1i + pp1i * c2;
                scalar(ov.dv1, fp1i[0], work1);
                scalar(pp1i, c2, work2);
                sum(work1, work2, pp1i);

                // penergy1i = (ov.dv1) * fenergy1i + penergy1i * c2;
                scalar(ov.dv1, fenergy1i[0], work1);
                scalar(penergy1i, c2, work2);
                sum(work1, work2, penergy1i);

                // Update subtree F1..i's for parent
                f1i[0] = ov.dvv1 * f1i[0];
                fp1i[0] = ov.dvv1 * fp1i[0];
                fenergy1i[0] = ov.dvv1 * fenergy1i[0];
            }
            return 1;


        }

        /**
         * Recursively traverses tree and computes volumes, etc.
         *
         * @param pos
         * @param volume
         * @param energy
         * @param dr
         * @param dv
         * @param free_volume
         * @param self_volume
         * @return
         */
        int computeVolume2R(double[][] pos, double[] volume, double[] energy,
                            double[][] dr, double[] dv,
                            double[] free_volume, double[] self_volume) {

            double[] psi1i = new double[1];
            double[] f1i = new double[1];
            double[] p1i = new double[3];
            double[] psip1i = new double[1];
            double[] fp1i = new double[1];
            double[] pp1i = new double[3];
            double[] energy1i = new double[1];
            double[] fenergy1i = new double[1];
            double[] penergy1i = new double[3];

            // Reset volumes, gradient
            for (int i = 0; i < dr.length; ++i) {
                dr[i][0] = 0.0;
                dr[i][1] = 0.0;
                dr[i][2] = 0.0;
            }
            fill(dv, 0.0);
            fill(free_volume, 0.0);
            fill(self_volume, 0.0);

            computeVolumeUnderSlot2R(0,
                    psi1i, f1i, p1i,
                    psip1i, fp1i, pp1i,
                    energy1i, fenergy1i, penergy1i,
                    dr, dv, free_volume, self_volume);

            volume[0] = psi1i[0];
            energy[0] = energy1i[0];

            return 1;

        }

        /**
         * Rescan the sub-tree to recompute the volumes, does not modify the tree.
         *
         * @param slot
         * @return
         */
        int rescanR(int slot) {
            int parent_index;

            // This overlap.
            GaussianOverlap ov = overlaps.get(slot);

            // Recompute its own overlap by merging parent and last atom.
            parent_index = ov.parentIndex;
            if (parent_index > 0) {
                GaussianVca g12 = new GaussianVca();
                double[] dVdr = new double[1];
                double[] dVdV = new double[1];
                double[] sfp = new double[1];

                int atom = ov.atom;
                GaussianOverlap parent = overlaps.get(parent_index);
                GaussianVca g1 = parent.g;

                // Atoms are stored in the tree at indexes 1...N
                GaussianVca g2 = overlaps.get(atom + 1).g;
                double gvol = overlapGaussianAlpha(g1, g2, g12, dVdr, dVdV, sfp);
                ov.g = g12;
                ov.volume = gvol;

                // dv1 is the gradient of V(123..)n with respect to the position of 1
                // ov.dv1 = ( g2.c - g1.c ) * (-dVdr);
                diff(g2.c, g1.c, ov.dv1);
                scalar(ov.dv1, -dVdr[0], ov.dv1);

                // dvv1 is the derivative of V(123...)n with respect to V(123...)
                ov.dvv1 = dVdV[0];
                ov.sfp = sfp[0];
                ov.gamma1i = parent.gamma1i + overlaps.get(atom + 1).gamma1i;
            }

            /* calls itself recursively on the children */
            for (int slot_child = ov.childrenStartindex; slot_child < ov.childrenStartindex + ov.childrenCount; slot_child++) {
                rescanR(slot_child);
            }

            return 1;

        }

        /**
         * Rescan the tree to recompute the volumes, does not modify the tree.
         *
         * @param pos        Atomic positions.
         * @param radii      Atomic radii.
         * @param volumes    Atomic volumes.
         * @param gammas     Atomic surface tensions.
         * @param ishydrogen True if the atom is a hydrogen.
         * @return
         */
        int rescanTreeV(double[][] pos, double[] radii,
                        double[] volumes, double[] gammas, boolean[] ishydrogen) {

            int slot = 0;
            GaussianOverlap ov = overlaps.get(slot);

            ov.level = 0;
            ov.volume = 0;
            ov.dv1 = new double[3];
            ov.dvv1 = 0.;
            ov.selfVolume = 0;
            ov.sfp = 1.;
            ov.gamma1i = 0.;

            slot = 1;
            for (int iat = 0; iat < nAtoms; iat++, slot++) {
                double a = KFC / (radii[iat] * radii[iat]);
                double vol = ishydrogen[iat] ? 0. : volumes[iat];
                ov = overlaps.get(slot);
                ov.level = 1;
                ov.g.v = vol;
                ov.g.a = a;
                ov.g.c = pos[iat];
                ov.volume = vol;
                ov.dv1 = new double[3];
                //dVi/dVi
                ov.dvv1 = 1.;
                ov.selfVolume = 0.;
                ov.sfp = 1.;
                // gamma[iat]/SA_DR
                ov.gamma1i = gammas[iat];
            }

            rescanR(0);

            return 1;
        }

        /**
         * Rescan the sub-tree to recompute the gammas, does not modify the volumes nor the tree.
         *
         * @param slot
         * @return
         */
        int rescanGammaR(int slot) {

            // This overlap.
            GaussianOverlap ov = overlaps.get(slot);

            // Recompute its own overlap by merging parent and last atom.
            int parent_index = ov.parentIndex;
            if (parent_index > 0) {
                int atom = ov.atom;
                GaussianOverlap parent = overlaps.get(parent_index);
                ov.gamma1i = parent.gamma1i + overlaps.get(atom + 1).gamma1i;
            }

            // Calls itself recursively on the children.
            for (int slot_child = ov.childrenStartindex; slot_child < ov.childrenStartindex + ov.childrenCount; slot_child++) {
                rescanGammaR(slot_child);
            }

            return 1;

        }

        /**
         * Rescan the tree to recompute the gammas only, does not modify volumes and the tree.
         *
         * @param gammas
         * @return
         */
        int rescanTreeG(double[] gammas) {

            int slot = 0;
            GaussianOverlap ov = overlaps.get(slot);
            ov.gamma1i = 0.;

            slot = 1;
            for (int iat = 0; iat < nAtoms; iat++, slot++) {
                ov = overlaps.get(slot);
                ov.gamma1i = gammas[iat];
            }

            rescanGammaR(0);
            return 1;
        }

        /**
         * Print the contents of the tree.
         */
        void printTree() {
            //logger.info("slot level LastAtom parent ChStart ChCount SelfV V gamma a x y z dedx dedy dedz sfp");
            for (int i = 1; i <= nAtoms; i++) {
                printTreeR(i);
            }

        }

        /**
         * Print the contents of the tree (recursive).
         *
         * @param slot
         */
        void printTreeR(int slot) {
            GaussianOverlap ov = overlaps.get(slot);
            logger.info(String.format("tg:      %d ", slot));
            ov.printOverlap();
            for (int i = ov.childrenStartindex; i < ov.childrenStartindex + ov.childrenCount; i++) {
                printTreeR(i);
            }
        }

        /**
         * Counts number of overlaps under the one given.
         *
         * @param slot
         * @return
         */
        int nchildrenUnderSlotR(int slot) {
            int n = 0;
            if (overlaps.get(slot).childrenCount > 0) {
                n += overlaps.get(slot).childrenCount;
                //now calls itself on the children
                for (int i = 0; i < overlaps.get(slot).childrenCount; i++) {
                    n += nchildrenUnderSlotR(overlaps.get(slot).childrenStartindex + i);
                }
            }
            return n;
        }
    }

    /**
     * Overlap volume switching function + 1st derivative
     *
     * @param gvol
     * @param volmina
     * @param volminb
     * @param sp
     * @return
     */
    private double polSwitchFunction(double gvol, double volmina, double volminb, double[] sp) {

        double swf = 0.0f;
        double swfp = 1.0f;
        double swd, swu, swu2, swu3, s;
        if (gvol > volminb) {
            swf = 1.0f;
            swfp = 0.0f;
        } else if (gvol < volmina) {
            swf = 0.0f;
            swfp = 0.0f;
        }
        swd = 1.f / (volminb - volmina);
        swu = (gvol - volmina) * swd;
        swu2 = swu * swu;
        swu3 = swu * swu2;
        s = swf + swfp * swu3 * (10.f - 15.f * swu + 6.f * swu2);
        sp[0] = swfp * swd * 30.f * swu2 * (1.f - 2.f * swu + swu2);

        //turn off switching function
        //*sp = 0.0;
        //s = 1.0;

        return s;
    }

    /**
     * Overlap between two Gaussians represented by a (V,c,a) triplet
     * <p>
     * V: volume of Gaussian
     * c: position of Gaussian
     * a: exponential coefficient
     * <p>
     * g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
     * <p>
     * this version is based on V=V(V1,V2,r1,r2,alpha)
     * alpha = (a1 + a2)/(a1 a2)
     * <p>
     * dVdr is (1/r)*(dV12/dr)
     * dVdV is dV12/dV1
     * dVdalpha is dV12/dalpha
     * d2Vdalphadr is (1/r)*d^2V12/dalpha dr
     * d2VdVdr is (1/r) d^2V12/dV1 dr
     *
     * @param g1   Gaussian 1.
     * @param g2   Gaussian 2.
     * @param g12  Overlap Gaussian.
     * @param dVdr is (1/r)*(dV12/dr)
     * @param dVdV is dV12/dV1
     * @param sfp
     * @return
     */
    private double overlapGaussianAlpha(GaussianVca g1, GaussianVca g2, GaussianVca g12,
                                        double[] dVdr, double[] dVdV, double[] sfp) {

        double d2, deltai, gvol, a12;
        double s, df, dgvol, dgvolv, ef;
        double[] c1 = g1.c;
        double[] c2 = g2.c;
        double[] dist = new double[3];
        double[] sp = new double[1];

        diff(c2, c1, dist);
        d2 = rsq(dist);

        a12 = g1.a + g2.a;
        deltai = 1. / a12;

        // 1/alpha
        df = (g1.a) * (g2.a) * deltai;

        ef = exp(-df * d2);
        gvol = ((g1.v * g2.v) / pow(PI / df, 1.5)) * ef;

        // (1/r)*(dV/dr) w/o switching function
        dgvol = -2.f * df * gvol;

        // (1/r)*(dV/dr) w/o switching function
        dgvolv = g1.v > 0 ? gvol / g1.v : 0.0;

        // Parameters for overlap gaussian. Note that c1 and c2 are Vec3's and the "*" operator wants
        // the vector first and scalar second vector2 = vector1 * scalar
        // g12.c = ((c1 * g1.a) + (c2 * g2.a)) * deltai;

        double[] c1a = new double[3];
        double[] c2a = new double[3];
        scalar(c1, g1.a * deltai, c1a);
        scalar(c2, g2.a * deltai, c2a);
        sum(c1a, c2a, g12.c);

        g12.a = a12;
        g12.v = gvol;

        // Switching function
        s = polSwitchFunction(gvol, VOLMINA, VOLMINB, sp);
        sfp[0] = sp[0] * gvol + s;
        dVdr[0] = dgvol;
        dVdV[0] = dgvolv;

        return s * gvol;
    }
}
