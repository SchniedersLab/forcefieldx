/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.numerics.MultipoleTensor;
import ffx.potential.bonded.Atom;

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link NeighborList} for any
 * {@link Crystal} space group. The real space contribution is contained within
 * this class and the reciprocal space contribution is delegated to the
 * {@link ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders
 */
public abstract class ParticleMeshEwald {

    protected final Logger logger = Logger.getLogger(this.getClass().getName());
    /**
     * Unit cell and spacegroup information.
     */
    protected Crystal crystal;
    /**
     * Number of symmetry operators.
     */
    protected int nSymm;
    /**
     * An ordered array of atoms in the system.
     */
    protected Atom atoms[];
    /**
     * The number of atoms in the system.
     */
    protected int nAtoms;

    /**
     * Axis defining atoms.
     */
    protected int axisAtom[][];

    /**
     * Polarization covalent lists
     */
    protected int ip11[][];
    protected int ip12[][];
    protected int ip13[][];

    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public enum Polarization {

        MUTUAL, DIRECT, NONE
    }

    protected Polarization polarization;

    public enum ELEC_FORM {

        PAM, FIXED_CHARGE
    }

    protected enum LambdaMode {
        OFF, CONDENSED, CONDENSED_NO_LIGAND, VAPOR
    }
    protected LambdaMode lambdaMode = LambdaMode.OFF;

    /**
     * Optionally predict induced dipoles prior to the SCF calculation.
     */
    protected ScfPredictor scfPredictor = null;
    /**
     * Dimensions of [nsymm][3][nAtoms].
     */
    protected double coordinates[][][];
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    protected int neighborLists[][][];

    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    protected double globalMultipole[][][];
    protected double globalMultipoleDot[][][];

    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    protected double inducedDipole[][][];
    protected double inducedDipoleCR[][][];

    /**
     * Number of unique tensors for given order.
     */
    protected static final int tensorCount = MultipoleTensor.tensorCount(3);
    protected static final double oneThird = 1.0 / 3.0;

    /**
     * PME real space cut-off.
     */
    protected double off;

    /**
     * Ewald coefficient.
     */
    protected double aewald;

    /**
     * SCF convergence criteria.
     */
    protected double poleps;

    /**
     * Reciprocal Space
     */
    protected ReciprocalSpace reciprocalSpace;

    public class EnergyForceTorque {

        public double energy;
        public double[] permFi = new double[3];
        public double[] permTi = new double[3];
        public double[] permFk = new double[3];
        public double[] permTk = new double[3];
        public double dPermdZ;
    }

    public class LambdaFactors {

        public double sc1 = 0.0;
        public double dsc1dL = 0.0;
        public double d2sc1dL2 = 0.0;
        public double sc2 = 1.0;
        public double dsc2dL = 0.0;
        public double d2sc2dL2 = 0.0;

        public LambdaFactors(double sc1, double dsc1dL, double d2sc1dL2,
                double sc2, double dsc2dL, double d2sc2dL2) {
            this.sc1 = sc1;
            this.dsc1dL = dsc1dL;
            this.d2sc1dL2 = d2sc1dL2;
            this.sc2 = sc2;
            this.dsc2dL = dsc2dL;
            this.d2sc2dL2 = d2sc2dL2;
        }
    }

    public int[][] getAxisAtoms() {
        return axisAtom;
    }

    public double getEwaldCutoff() {
        return off;
    }

    public double getEwaldCoefficient() {
        return aewald;
    }

    public ReciprocalSpace getReciprocalSpace() {
        return reciprocalSpace;
    }

    public Polarization getPolarizationType() {
        return polarization;
    }

    public double getPolarEps() {
        return poleps;
    }

    public int[][] getPolarization11() {
        return ip11;
    }

    public int[][] getPolarization12() {
        return ip12;
    }

    public int[][] getPolarization13() {
        return ip13;
    }

    /**
     * <p>
     * Getter for the field <code>gradient</code>.</p>
     *
     * @return an array of double.
     */
    protected abstract double[][][] getGradient();

    /**
     * <p>
     * Getter for the field <code>torque</code>.</p>
     *
     * @return an array of double.
     */
    protected abstract double[][][] getTorque();

    protected abstract double[][][] getLambdaGradient();

    protected abstract double[][][] getLambdaTorque();

    public abstract void setAtoms(Atom atoms[], int molecule[]);

    public abstract void setFixedCharges(Atom atoms[]);

    public abstract double energy(boolean gradient, boolean print);

    public abstract double getPermanentEnergy();

    public abstract double getPermanentRealSpaceEnergy();

    public abstract double getPermanentReciprocalEnergy();

    public abstract String getDecomposition();

    public abstract double getPolarizationEnergy();

    public abstract int getInteractions();

    public abstract double getGKEnergy();

    public abstract GeneralizedKirkwood getGK();

    public abstract int getGKInteractions();

    public abstract void setLambda(double lambda);

    public abstract double getdEdL();

    public abstract void getdEdXdL(double[] gradients);

    public abstract double getd2EdL2();

    public abstract double[] getdEdEsv();

    public abstract void destroy() throws Exception;

    public abstract void setCrystal(Crystal crystal);

    public abstract double getCavitationEnergy(boolean throwError);

    public abstract double getDispersionEnergy(boolean throwError);

    private void log(int i, int k, double r, double eij) {
        log("ELEC", i, k, r, eij);
    }

    /**
     * Log the real space electrostatics interaction.
     *
     * @param i Atom i.
     * @param k Atom j.
     * @param r The distance rij.
     * @param eij The interaction energy.
     */
    private void log(String type, int i, int k, double r, double eij) {
        logger.info(String.format("%s %6d-%s %6d-%s %10.4f  %10.4f",
                type, atoms[i].getIndex(), atoms[i].getAtomType().name, atoms[k].getIndex(), atoms[k].getAtomType().name, r, eij));
    }

}
