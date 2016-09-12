/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.nonbonded;

import java.util.List;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedVariable;

/**
 *
 * @author mjschnie
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
    
    public abstract double[] getdEdLdh();
    
    public abstract double[] getd2EdLdh2();
    
    public abstract double[][] getdEdXdLdh();

    public abstract void setESVList(List<ExtendedVariable> esvList);

    public abstract void destroy() throws Exception;

    public abstract void setCrystal(Crystal crystal);

    public abstract double getCavitationEnergy(boolean throwError);

    public abstract double getDispersionEnergy(boolean throwError);
}
