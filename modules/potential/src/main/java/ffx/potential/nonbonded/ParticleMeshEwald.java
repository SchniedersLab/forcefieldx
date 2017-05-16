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
import ffx.potential.bonded.Atom;

import static ffx.potential.extended.ExtUtils.formatArray;

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
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public Polarization polarization;

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

    public enum SCFAlgorithm {
        SOR, CG
    }

    public enum SCFPredictor {
        NONE, LS, POLY, ASPC
    }

    public abstract double getTotalMultipoleEnergy();

    public abstract double getPermanentEnergy();

    public abstract double getPermRealEnergy();

    public abstract double getPermSelfEnergy();

    public abstract double getPermRecipEnergy();

    public abstract double getPolarizationEnergy();

    public abstract double getIndRealEnergy();

    public abstract double getIndSelfEnergy();

    public abstract double getIndRecipEnergy();

    public abstract double getGKEnergy();

    public abstract GeneralizedKirkwood getGK();

    public static class LambdaFactors {

        public final double sc1;
        public final double dsc1dL;
        public final double d2sc1dL2;
        public final double sc2;
        public final double dsc2dL;
        public final double d2sc2dL2;
        public static final LambdaFactors Defaults
                = new LambdaFactors(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

        public LambdaFactors(double sc1, double dsc1dL, double d2sc1dL2,
                double sc2, double dsc2dL, double d2sc2dL2) {
            this.sc1 = sc1;
            this.dsc1dL = dsc1dL;
            this.d2sc1dL2 = d2sc1dL2;
            this.sc2 = sc2;
            this.dsc2dL = dsc2dL;
            this.d2sc2dL2 = d2sc2dL2;
        }

        @Override
        public String toString() {
            return formatArray(new double[]{sc1, dsc1dL, d2sc1dL2, sc2, dsc2dL, d2sc2dL2});
        }
    }

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

    public abstract void setLambda(double lambda);

    public abstract double getdEdL();

    public abstract void getdEdXdL(double[] gradients);

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
}
