package ffx.potential.extended;

import java.util.Arrays;
import java.util.List;

import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * 
 * @author slucore
 */
public class ExtConstants {
    
    public static final List<String> titratableHydrogenNames = 
            Arrays.asList("HH", "HG", "HE2", "HD1", "HE2", "HD2", "HZ3");
    public static final List<String> backboneNames = Arrays.asList("N","CA","C","O","HA","H");
    
    /**
     * Boltzmann's constant is kcal/mol/Kelvin.
     */
    public static final double Boltzmann = 0.0019872041;
    public static final double beta = 1 / Boltzmann;
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K.
     */
    public static final double kB = 0.83144725;
    /**
     * Conversion from kcal/mole to g*Ang**2/ps**2.
     */
    public static final double convert = 4.1840e2;
    /**
     * Gas constant (in Kcal/mole/Kelvin).
     */
    public static final double R = 1.9872066e-3;
    /**
     * Random force conversion to kcal/mol/A; formerly randomForce.
     */
    public static final double forceToKcal = sqrt(4.184) / 10e9;
    /**
     * Random force conversion to (kcal/mol/A)^2; formerly randomForce2.
     */
    public static final double forceToKcalSquared = forceToKcal * forceToKcal;
    
    public static final double roomTemperature = 298.15;
    
    public static final double log10 = Math.log(10);
    
}
