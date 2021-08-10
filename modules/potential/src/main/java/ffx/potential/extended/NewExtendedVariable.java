// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.extended;

import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.extended.NewExtendedSystem.NewExtendedSystemConfig;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.*;
import java.util.logging.Logger;

import static ffx.utilities.Constants.kB;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * A generalized extended system variable. Treatment of ESVs: a. Bonded terms interpolate linearly
 * between end states. ESVs based on MultiResidue (e.g. TitrationESV) place a multiplier in the term
 * objects themselves. b. PME and vdW scaling and derivatives are handled inside these classes'
 * inner loops. Softcoring follows the same form used by OST lambda.
 *
 * @author Andrew Thiel
 * @since 1.0
 */
public abstract class NewExtendedVariable {

    /**
     * Constant <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(NewExtendedVariable.class.getName());

    private final boolean isConstantESV;
    /**
     * Index of this ESV in the list of ESVs.
     */
    public final int esvIndex;
    /**
     * Atoms affected by this extended variable
     */
    protected final List<Atom> atomsExtended;
    /**
     * Maps multipole end points of this ESV's lambda path.
     */
    protected final AminoAcidUtils.AminoAcid3 aminoAcid3;
    /**
     * TODO: Replace ExtendedSystemConfig with an instance of normal
     * FFX Properties.
     */
    protected final NewExtendedSystemConfig config;
    /**
     * ESVs travel on {0,1}.
     */
    protected double lambda = 0.0;
    /**
     * Propagates lambda particle via:
     * lambda = sin(theta)^2
     */
    protected double theta = 0.0;
    /**
     * Velocity of the theta particle (radian/psec).
     */
    protected double thetaVelocity = 0.0;
    /**
     * Acceleration of the theta particle (radian/psec^2).
     */
    protected double thetaAccel = 0.0;
    /**
     * Mass of the theta particle (amu).
     */
    protected final double thetaMass;
    /**
     * Magnitude of the discretization bias (kcal/mol).
     */
    protected final double discrBiasMag;
    /**
     * Discretization bias.
     */
    protected double discrBias;
    /**
     * Discretization bias chain rule derivative.
     */
    protected double dDiscrBiasdL;
    /**
     * Sigmoidal switching function.
     *
     * Maps lambda -> S(lambda) which has a flatter derivatiuve near zero/unity.
     * Properties: {S(0)=0, S(1)=1, dS(0)=0, dS(1)=0, S(1-L)=1-S(L)}.
     */
    protected final MultiplicativeSwitch switchingFunction;
    /**
     * Switched lambda value.
     */
    protected double lSwitch;
    /**
     * Switched lambda chain rule derivative.
     */
    protected double dlSwitch;
    /**
     * Utils for setting ESV multipoles and derivs
     */
    ConstantPhUtils constantPhUtils;

    /**
     * Prefer ExtendedSystem::populate to manual ESV creation.
     *
     * @param residue:      from TitrationUtils.titrationFactory()
     * @param initialLambda: (optional) starting position of the extended particle
     * @param esvSystem      a {@link ExtendedSystem} object.
     */
    public NewExtendedVariable(NewExtendedSystem esvSystem, Residue residue, double initialLambda) {
        this.esvIndex = esvSystem.requestIndexing();
        this.config = esvSystem.config;
        this.discrBiasMag = config.discrBias;
        this.switchingFunction = new MultiplicativeSwitch(1.0, 0.0);
        this.aminoAcid3 = residue.getAminoAcid3();
        if(initialLambda < 0.0){this.isConstantESV = true;}
        else{this.isConstantESV = false;}

        MolecularAssembly molecularAssembly = esvSystem.getMolecularAssembly();
        ForceField forceField = molecularAssembly.getForceField();
        constantPhUtils = new ConstantPhUtils(forceField);
        CompositeConfiguration properties = molecularAssembly.getProperties();
        this.thetaMass = properties.getDouble("esv.mass",ExtendedSystem.THETA_MASS);
        setInitialLambda(initialLambda);

        atomsExtended = new ArrayList<>();
        atomsExtended.addAll(residue.getAtomList());

        /* Background atoms don't get automatically typed by PME since they're
         * disconnected from the molecular assembly; must be done manually. */
        if (config.electrostatics) {
            updateMultipoleTypes();
        }
    }

    /**
     * Getter for the field <code>esvIndex</code>.
     *
     * @return a int.
     */
    public final int getEsvIndex() {
        return esvIndex;
    }

    /**
     * The unswitched lambda value, ie input to S(L).
     *
     * @return a double.
     */
    public final double getLambda() {
        return lambda;
    }

    boolean isConstantESV(){
        return isConstantESV;
    }
    /**
     * Should only be called by ExtendedSystem since an updateListeners() call is required afterwards
     * to notify VdW and PME.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        setLambda(lambda, true);
    }

    /**
     * getName.
     *
     * @return a {@link String} object.
     */
    public String getName() {
        return String.format("Esv%d", esvIndex);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return String.format("%s (%4.2f)", this.getName(), getLambda());
    }

    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     *
     * @param temperature a double.
     * @param print       a boolean.
     * @return a double.
     */
    protected abstract double getTotalBias(double temperature, boolean print);

    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     *
     * @param temperature a double.
     * @param print       a boolean.
     * @return a double.
     */
    protected abstract double getTotalBiasDeriv(double temperature, boolean print);

    /**
     * First pass setting of lambda to initialize theta and theta velocity
     * @param lambda a double
     */
    private void setInitialLambda(double lambda) {
        setTheta(Math.asin(Math.sqrt(lambda)));
        Random random = new Random();
        setThetaVelocity(random.nextGaussian() * sqrt(kB * 298.15 / thetaMass));
        this.lambda = lambda;
        this.lSwitch = (config.allowLambdaSwitch) ? switchingFunction.taper(lambda) : lambda;
        this.dlSwitch = (config.allowLambdaSwitch) ? switchingFunction.dtaper(lambda) : 1.0;
        discrBias = - 4.0 * discrBiasMag * (lambda - 0.5) * (lambda - 0.5);
        dDiscrBiasdL = -8.0 * discrBiasMag * (lambda - 0.5);
    }

    /**
     * Used when a manual setting of lambda is needed; not during Langevin dynamics
     * @param lambda a double
     * @param updateComponents a boolean
     */
    private void setLambda(double lambda, boolean updateComponents) {
        setTheta(Math.asin(Math.sqrt(lambda)));
        updateLambda(lambda, updateComponents);
    }

    /**
     * Used to update lambda during dynamics, does not update theta value.
     * Bias update from Eq. 3 of "All-Atom Continuous Constant pH Molecular Dynamics..." J. Shen 2016.
     * @param lambda a double
     * @param updateComponents a boolean
     */
    protected void updateLambda(double lambda, boolean updateComponents) {
        this.lambda = lambda;
        this.lSwitch = (config.allowLambdaSwitch) ? switchingFunction.taper(lambda) : lambda;
        this.dlSwitch = (config.allowLambdaSwitch) ? switchingFunction.dtaper(lambda) : 1.0;
        discrBias = - 4.0 * discrBiasMag * (lambda - 0.5) * (lambda - 0.5);
        dDiscrBiasdL = -8.0 * discrBiasMag * (lambda - 0.5);
        if (updateComponents) {
            updateMultipoleTypes();
        }
    }

    protected void setTheta(double theta) {
        this.theta = theta;
    }

    protected double getTheta() {
        return theta;
    }

    protected void setThetaVelocity(double thetaVelocity) {
        this.thetaVelocity = thetaVelocity;
    }

    protected double getThetaVelocity() {
        return thetaVelocity;
    }

    protected void setThetaAccel(double thetaAccel) {
        this.thetaAccel = thetaAccel;
    }

    protected double getThetaAccel() {
        return thetaAccel;
    }

    protected double getThetaMass(){return thetaMass; }

    /**
     * Invoked by ExtendedSystem after lambda changes.
     */
    protected abstract void updateMultipoleTypes();

    /**
     * getLambdaSwitch.
     *
     * @return a double.
     */
    protected final double getLambdaSwitch() {
        return (config.allowLambdaSwitch) ? lSwitch : lambda; // S(L)
    }

    /**
     * getSwitchDeriv.
     *
     * @return a double.
     */
    protected final double getSwitchDeriv() {
        return (config.allowLambdaSwitch) ? dlSwitch : 1.0; // dS(L)dL
    }

    /**
     * From Shen and Huang 2016; drives ESVs to zero/unity. bias = 4B*(L-0.5)^2
     *
     * @return a double.
     */
    protected double getDiscrBias() {
        return discrBias;
    }

    /**
     * dBiasdL = -8B*(L-0.5)
     *
     * @return a double.
     */
    protected double getDiscrBiasDeriv() {
        return dDiscrBiasdL;
    }

    /**
     * viewUnsharedAtoms.
     *
     * @return a {@link List} object.
     */
    protected List<Atom> viewExtendedAtoms() {
        return Collections.unmodifiableList(atomsExtended);
    }
}
