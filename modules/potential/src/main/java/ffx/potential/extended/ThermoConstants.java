package ffx.potential.extended;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * 
 * @author slucore
 */
public class ThermoConstants {

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
    public static final double beta = 1 / R;
    /**
     * Random force conversion to kcal/mol/A;
     */
    public static final double randomConvert = sqrt(4.184) / 10e9;
    public static final double randomConvert2 = randomConvert * randomConvert;
    private static final Random random = ThreadLocalRandom.current();
    
    private static final double roomTemperature = 298.15;
    
    public static final double log10 = Math.log(10);
    
    /**
     * Return room-temperature velocities from a Maxwell-Boltzmann distribution of momenta.
     */
    public static double[] singleRoomtempMaxwell(double mass) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = random.nextGaussian() * sqrt(kB * roomTemperature / mass);
        }
        return vv;
    }
    
    /**
     * Return velocities from a Maxwell-Boltzmann distribution of momenta.
     * The variance of each independent momentum component is kT * mass.
     */
    public double[] singleMaxwell(double mass, double temperature) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = random.nextGaussian() * sqrt(kB * temperature / mass);
        }
        return vv;
    }
    
}
