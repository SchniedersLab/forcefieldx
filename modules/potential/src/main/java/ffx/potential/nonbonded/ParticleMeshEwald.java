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

import java.util.List;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedVariable;

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

    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public enum Polarization {

        MUTUAL, DIRECT, NONE
    }

    public enum ELEC_FORM {

        PAM, FIXED_CHARGE
    }

    /**
     * Polarization mode.
     */
    protected Polarization polarization;

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

    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    protected double inducedDipole[][][];
    protected double inducedDipoleCR[][][];

    public abstract double getEwaldCutoff();

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

    public abstract double getPolarizationEnergy();

    public abstract int getInteractions();

    public abstract double getGKEnergy();

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
}
