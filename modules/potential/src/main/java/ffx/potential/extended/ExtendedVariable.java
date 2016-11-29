package ffx.potential.extended;

import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.Atom;

import static ffx.potential.extended.ExtUtils.prop;

/**
 * A generalized extended system variable.
 * Treatment of ESVs: 
 *  a. Bonded terms interpolate linearly between end states.
 *      ESVs based on MultiResidue (e.g. TitrationESV) place a multiplier in the term objects themselves.
 *  b. PME and vdW scaling and derivatives are handled inside these classes' inner loops.
 *      Softcoring follows the same form used by OSRW lambda.
 * @author slucore
 */
public abstract class ExtendedVariable {
    
    // System handles
    private static final Logger logger = Logger.getLogger(ExtendedVariable.class.getName());
    private static int esvIndexer = 0;
    public final int index;
    protected final List<Atom> atoms = new ArrayList<>();
    protected final List<Atom> backbone = new ArrayList<>();
    boolean ready = false;
    
    // Static Properties
    protected static final boolean esvPropagation = prop("esv-propagation", true);
    protected static final boolean scaleBondedTerms = prop("esv-scaleBonded", true);
    protected static final OptionalDouble biasOverride = prop("esv-bias", OptionalDouble.empty());
    
    // Instance Properties
    private final double thetaMass = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
    private final double thetaFriction = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps
    
    // Lamedh variables
    private double lambda;                        // ESVs travel on {0,1}
    private double theta;                           // Propagates lamedh particle via "lamedh=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final Random stochasticRandom = ThreadLocalRandom.current();
    
    private final double betat;
    private double discrBias, dDiscrBiasdL;
    
    public ExtendedVariable(double biasMag, double initialLamedh) {
        index = esvIndexer++;
        betat = biasOverride.isPresent() ? biasOverride.getAsDouble() : biasMag;
        lambda = initialLamedh;
        theta = Math.asin(Math.sqrt(lambda));
    }
    
    public ExtendedVariable(double biasMag) {
        this(biasMag, 1.0);
    }    
    public ExtendedVariable() {
        this(0.0, 1.0);
    }
        
    /**
     * Propagate lambda using Langevin dynamics.
     * Check that temperature goes to the value used below (when set as a constant) even when sim is decoupled.
     */
    public void propagate(double dEdLdh, double dt, double setTemperature) {
        if (!esvPropagation) {
            return;
        }
        double rt2 = 2.0 * ExtConstants.BOLTZMANN * setTemperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / ExtConstants.randomConvert;
        double dEdL = -dEdLdh * sin(2.0 * theta);
        halfThetaVelocity = (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
                + ExtConstants.randomConvert2 * 2.0 * dt * (dEdL + randomForce))
                / (2.0 * thetaMass + thetaFriction * dt);
        theta = theta + dt * halfThetaVelocity;

        if (theta > PI) {
            theta -= 2.0 * PI;
        } else if (theta <= -PI) {
            theta += 2.0 * PI;
        }

        double sinTheta = sin(theta);
        setLambda(sinTheta * sinTheta);
    }
    
    public abstract void updateBondedLambdas();
    
    public void setLambda(double lambda) {
        this.lambda = lambda;
        theta = Math.asin(Math.sqrt(lambda));
        discrBias = betat - 4*betat*(lambda-0.5)*(lambda-0.5);
        dDiscrBiasdL = -8*betat*(lambda-0.5);
        updateBondedLambdas();
    }
    
    public final double getLambda() {
        return lambda;
    }
    
    public final int getIndex() {
        return index;
    }
    
    @Override
    public String toString() {
        return String.format("ESV%d", index);
    }
    
    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double totalBiasEnergy(double temperature, StringBuilder sb);
    public abstract double totalBiasDeriv(double temperature, StringBuilder sb);
    
    /**
     * From Shen&Huang 2016; drives ESVs to zero/unity.
     * bias = 4B*(L-0.5)^2
     */
    public double discretizationBiasEnergy() {
        return discrBias;
    }
    
    /**
     * dBiasdL = -8B*(L-0.5)
     */
    public double discretizationBiasDeriv() {
        return dDiscrBiasdL;
    }
    
    /**
     * Implementations should fill the protected atoms[] array, set the 
     * esvLambda value of any affected bonded terms, and set the ready flag.
     */
    public abstract void readyup();
    
    public boolean isReady() {
        return ready;
    }
    
}
