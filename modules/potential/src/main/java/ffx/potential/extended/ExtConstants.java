package ffx.potential.extended;

import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.potential.extended.ExtUtils.prop;

/**
 * 
 * @author slucore
 */
public class ExtConstants {
    
    public static final boolean DEBUG = prop("esv-debug", false);
    
    public static final boolean VERBOSE = prop("esv-verbose", false);
    
    /**
     * Boltzmann's constant is kcal/mol/Kelvin.
     */
    public static final double BOLTZMANN = 0.0019872041;
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K
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
     * One over kbT.
     */
    public static final double beta = 1 / BOLTZMANN;
    /**
     * Random force conversion to kcal/mol/A;
     */
    public static final double randomConvert = sqrt(4.184) / 10e9;
    public static final double randomConvert2 = randomConvert * randomConvert;
    
    public static final double roomTemperature = 298.15;
    
    public static final double log10 = Math.log(10);
    
}
