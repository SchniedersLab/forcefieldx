package ffx.potential.extended;

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
 *  1. Potential terms interpolate linearly between the two lamedh end states.
 *  2. Treatment of ESVs: 
 *      a. Bonded terms are handled at the level of ForceFieldEnergy (or, generically, at Potential).
 *      b. PME and vdW scaling and derivatives are handled inside these classes' inner loops.
 *  3. Interaction of lamedh with lambda:
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
    private final MolecularAssembly mola;
    private static int esvIndexer = 0;
    public final int index;
    protected List<Atom> atoms;
    
    // Thermo variables
    private double temperature = 298.15;            // Kelvin
    private double dt = 1.0;                        // femtoseconds
    
    // Lamedh variables
    protected double lamedh;                        // ESVs travel on {0,1}
    private double theta;                           // Propagates lamedh particle via "lamedh=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final double thetaMass = 1.0e-18;       // from OSRW, reasonably 100 a.m.u.
    private final double thetaFriction = 1.0e-19;   // from OSRW, reasonably 60/ps
    private final Random stochasticRandom = ThreadLocalRandom.current();
    
    public ExtendedVariable(MolecularAssembly mola, double temperature, double dt) {
        this.mola = mola;
        this.temperature = temperature;
        this.dt = dt;
        this.index = esvIndexer++;
        setLamedh(1.0);
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
     */
    public void propagateLamedh() {
        double rt2 = 2.0 * ThermoConstants.R * temperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / ThermoConstants.randomConvert;
        double dEdLamedh = getdEdLamedh();
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
        lamedh = sinTheta * sinTheta;
    }
    
    public final void setLamedh(double lamedh) {
        this.lamedh = lamedh;
        theta = Math.asin(Math.sqrt(lamedh));
    }
    public final double getLamedh() {
        return lamedh;
    }
    public final int getIndex() {
        return index;
    }
    
    /**
     * Declared abstract as a reminder to fill local array and update Atom fields.
     */
    protected abstract void setAtoms();
    /**
     * Necessary for particle propagation.
     */
    public abstract double getdEdLamedh();
    /**
     * Called by ForceFieldEnergy to apply ESVs.
     */
    public abstract OptionalDouble getROLSScaling(ROLS rols);
    
}
