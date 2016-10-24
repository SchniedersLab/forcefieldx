package ffx.potential.extended;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.ROLS;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * A generalized extended system variable.
 * Notes:
 *  1. Bonded terms interpolate linearly between the two lamedh end states.
 *  2. VdW and PME use softcoring (nonlinear) as with lambda.
 *  3. Treatment of ESVs: 
 *      a. Bonded terms are handled at the level of ForceFieldEnergy (or, generically, at Potential).
 *      b. PME and vdW scaling and derivatives are handled inside these classes' inner loops.
 *  4. Interaction of lamedh with lambda:
 *      a. dU/dLambda/dLamedh is taken to be nil.
 *      b. Lambda (OSRW) statistics are to be collected only at zero or unity lamedh.
 *      c. Lambda-scaling on atoms and their potential terms stacks multiplicatively with lamedh scaling.
 *              e.g. An ASH-ALA transition is coupled to lambda, while lamedh couples ASH-ASP.
 *                   ASH-HD2 vdW interactions in this case are scaled by lambda*lamedh.
 *              TODO: Decide on the form of the softcore vdW denominator scaling.
 *      d. Bonded term potential interpolation along lambda is handled at the DualTopology level
 *              and thus naturally encapsulates lamedh.
 *      e. PME and vdW interaction is handled explicitly in these classes.
 *              Since there are no lambda derivs at intermediate lamedh and since these
 *              derivatives take the same form, the vdW code for dUdLambda is reused: 
 *              albeit with an extra dimension (for multiple lamedhs).
 * @author slucore
 */
public abstract class ExtendedVariable {
    
    // System handles
    private static int esvIndexer = 0;
    public final int index;
    protected List<Atom> atoms = new ArrayList<>();
    
    // Lamedh variables
    protected double lambda;                        // ESVs travel on {0,1}
    private double theta;                           // Propagates lamedh particle via "lamedh=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final double thetaMass = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
    private final double thetaFriction = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps
    private final Random stochasticRandom = ThreadLocalRandom.current();
    
    private final double betat;
    
    public ExtendedVariable(double biasMag, double initialLamedh) {
        index = esvIndexer++;
        betat = System.getProperty("esv-bias") == null ? biasMag : 
                Double.parseDouble(System.getProperty("esv-bias"));
        lambda = initialLamedh;
        theta = Math.asin(Math.sqrt(lambda));
    }
    
    public ExtendedVariable(double biasMag) {
        this(biasMag, 1.0);
    }    
    public ExtendedVariable() {
        this(0.0, 1.0);
    }
    
    public static int prop(String key, int def) {
        return (System.getProperty(key) != null) ? Integer.parseInt(System.getProperty(key)) : def;
    }
    public static double prop(String key, double def) {
        return (System.getProperty(key) != null) ? Double.parseDouble(System.getProperty(key)) : def;
    }
    
    public List<Atom> getAtoms() {
        List<Atom> ret = new ArrayList<>();
        ret.addAll(atoms);
        return ret;
    }
    
    public boolean containsAtom(Atom atom) {
        return atoms.contains(atom);
    }
    
    /**
     * Propagate lamedh using Langevin dynamics.
     * Check that temperature goes to the value used below (when set as a constant) even when sim is decoupled.
     */
    public void propagateLamedh(double dEdLamedh, double currentTemperature, double dt) {
        double rt2 = 2.0 * ThermoConstants.R * currentTemperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / ThermoConstants.randomConvert;
//        double dEdLamedh = getdEdLdh();   // TODO ASSESS architecture that assigns getdEdLdh() to Potential
        double dEdL = -dEdLamedh * sin(2.0 * theta);
        halfThetaVelocity = (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
                + ThermoConstants.randomConvert2 * 2.0 * dt * (dEdL + randomForce))
                / (2.0 * thetaMass + thetaFriction * dt);
        theta = theta + dt * halfThetaVelocity;

        if (theta > PI) {
            theta -= 2.0 * PI;
        } else if (theta <= -PI) {
            theta += 2.0 * PI;
        }

        double sinTheta = sin(theta);
        lambda = sinTheta * sinTheta;
    }
    
    public final void setLambda(double lambda) {
        this.lambda = lambda;
        theta = Math.asin(Math.sqrt(lambda));
    }
    public final double getLambda() {
        return lambda;
    }
    public final int getIndex() {
        return index;
    }
    
    /**
     * From Shen&Huang 2016; drives ESVs to zero/unity.
     * bias = 4B*(L-0.5)^2
     */
    public double getBiasEnergy() {
        return (4*betat - (lambda-0.5)*(lambda-0.5));
    }
    
    /**
     * dBiasdL = -8B*(L-0.5)
     */
    public double getdBiasdLdh() {
        return (-8*betat*(lambda-0.5));
    }
    
    /**
     * d2BiasdL2 = -8B
     */
    public double getd2BiasdLdh2() {
        return -8*betat;
    }
    
    /**
     * Declared abstract as a reminder to both fill local array and update Atom fields.
     */
    protected abstract void finalize();
    /**
     * Called by ForceFieldEnergy to apply ESVs to bonded terms.
     */
    public abstract OptionalDouble getROLSScaling(ROLS rols);
    
}
