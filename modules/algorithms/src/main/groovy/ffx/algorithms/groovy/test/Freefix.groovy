package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import org.apache.commons.math3.special.Erf
import static java.lang.String.format
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

class Freefix extends AlgorithmsScript {


    @Option(names = ['--ri'], paramLabel = '0.00',
            description = 'Constant Inner Radius value.')
    double innerRadius = 0.00

    @Option(names = ['--fi'], paramLabel = '0.0',
            description = 'Constant Inner Force value.')
    double innerForce = 0.0

    @Option(names = ['--ro'], paramLabel = '0.00',
            description = 'Constant Outer Radius value.')
    double outerRadius = 0.00

    @Option(names = ['--fo'], paramLabel = '15.0',
            description = 'Constant Outer Force value.')
    double outerForce = 15.0

    @Option(names = ['--temp'], paramLabel = '298.0',
            description = 'Constant Temperature value.')
    double temperature = 298.0

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate files in PDB or XYZ format.')
    private String filename = null

    private static final double AVOGADRO = 6.02214076e23;
    private static final double STD_CONVERSION = 1.0e27 / AVOGADRO;

    /**
     * Freefix constructor.
     */
    Freefix() {
        super()
    }

    /**
     * Freefix constructor.
     * @param binding The Groovy Binding to use.
     */
    Freefix(Binding binding) {
        super(binding)
    }

    @Override
    Freefix run() {

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        if (innerForce == 0.0) {
            innerForce = 1.0;
            innerRadius = 0.0;
        }
        if (outerForce == 0.0) {
            outerForce = 1.0;
        }

        System.out.printf("%-35s %.4f Ang%n", "Inner Flat-Bottom Radius:", innerRadius);
        System.out.printf("%-35s %.4f Kcal/mole/Ang^2%n", "Inner Force Constant:", innerForce);
        System.out.printf("%-35s %.4f Ang%n", "Outer Flat-Bottom Radius:", outerRadius);
        System.out.printf("%-35s %.4f Kcal/mole/Ang^2%n", "Outer Force Constant:", outerForce);
        System.out.printf("%-35s %.4f Kelvin%n%n", "System Temperature Value:", temperature);

        double kB = 0.0019872041;
        double kt = kB * temperature;

        double[] volumeResults = computeVolumeIntegrals(innerRadius, innerForce, outerRadius, outerForce, temperature);
        double vol = volumeResults[0];
        double dvol = volumeResults[1];

        logger.info(
                format("%-35s %.4f Ang^3%n", "Analytical Volume Integral:", vol));
        logger.info(
                format("%-35s %.4f Ang^3/K%n", "Analytical dVol/dT Value:", dvol));

        double dg = -kt * Math.log(vol / STD_CONVERSION);
        double ds = -dg / temperature + kt * dvol / vol;
        double dh = dg + temperature * ds;

        logger.info(
                format("%-35s %.4f Kcal/mole%n", "Restraint Free Energy:", dg));
        logger.info(
                format("%-35s %.4f Kcal/mole/K%n", "Restraint Entropy Value:", ds));
        logger.info(
                format("%-35s %.4f Kcal/mole%n", "Restraint Enthalpy Value:", dh));
        logger.info(
                format("%-35s %.4f Kcal/mole%n", "Restraint -T deltaS Value:", -temperature * ds));

        return this
    }

    static double[] computeVolumeIntegrals(double innerRadius, double innerForce, double outerRadius, double outerForce, double temperature) {
        double kB = 0.0019872041;
        double kt = kB * temperature;
        double volume1 = 0.0, volume2, volume3 = 0.0;

        if (innerRadius != 0.0) {
            double term1_v1 = 2.0 * Math.PI * innerRadius * (-2.0 + Math.exp(-Math.pow(innerRadius, 2) * innerForce / kt)) * kt / innerForce;
            double term2_v1 = Math.sqrt(kt * Math.pow(Math.PI / innerForce, 3)) *
                    (2.0 * innerForce * Math.pow(innerRadius, 2) + kt) * Erf.erf(innerRadius * Math.sqrt(innerForce / kt));
            volume1 = term1_v1 + term2_v1;
        }

        volume2 = (4.0 / 3.0) * Math.PI * (Math.pow(outerRadius, 3) - Math.pow(innerRadius, 3));

        if (outerRadius == 0.0) {
            double term1_v3 = Math.sqrt(kt * (Math.PI / outerForce)**3) * kt
            volume3 = term1_v3
        } else {
            double term1_v3 = Math.sqrt(kt * Math.pow(Math.PI / outerForce, 3)) *
                    (2.0 * outerForce * Math.pow(outerRadius, 2) + kt + 4.0 * outerRadius * Math.sqrt(kt * outerForce / Math.PI));
            volume3 = term1_v3;
        }


        double volume = volume1 + volume2 + volume3;

        double dv1 = (innerRadius != 0.0) ?
                computeDv1(innerRadius, innerForce, kt, temperature) : 0.0;
        double dv3 = (outerRadius != 0.0) ?
                computeDv3(outerRadius, outerForce, kt, temperature) :
                1.5 * Math.pow(kt * Math.PI / outerForce, 1.5) / temperature;

        double dvolume = dv1 + dv3;
        return new double[]{volume, dvolume};
    }

    private static double computeDv1(double innerRadius, double innerForce, double kt, double temperature) {
        double term1_dv1 = 2.0 * Math.PI * Math.pow(innerRadius, 3) * Math.exp(-Math.pow(innerRadius, 2) * innerForce / kt) / temperature;
        double term2_dv1 = 2.0 * Math.PI * innerRadius * (-2.0 + Math.exp(-Math.pow(innerRadius, 2) * innerForce / kt)) * kt / (innerForce * temperature);
        double term3_dv1 = 0.5 * Math.sqrt(Math.pow(Math.PI / innerForce, 3)) * Math.sqrt(kt) *
                (2.0 * innerForce * Math.pow(innerRadius, 2) + kt) * Erf.erf(innerRadius * Math.sqrt(innerForce / kt)) / temperature;
        double term4_dv1 = -Math.PI * innerRadius * Math.exp(-Math.pow(innerRadius, 2) * innerForce / kt) *
                (2.0 * innerForce * Math.pow(innerRadius, 2) + kt) / (innerForce * temperature);
        double term5_dv1 = Math.sqrt(Math.pow(kt * Math.PI / innerForce, 3)) * Erf.erf(innerRadius * Math.sqrt(innerForce / kt)) / temperature;
        return term1_dv1 + term2_dv1 + term3_dv1 + term4_dv1 + term5_dv1;
    }

    private static double computeDv3(double outerRadius, double outerForce, double kt, double temperature) {
        double term1_dv3 = Math.sqrt(kt * Math.pow(Math.PI / outerForce, 3)) * outerForce * Math.pow(outerRadius, 2) / temperature;
        double term2_dv3 = 4.0 * kt * (Math.PI / outerForce) * outerRadius / temperature;
        double term3_dv3 = 1.5 * Math.sqrt(Math.pow(kt * Math.PI / outerForce, 3)) / temperature;
        return term1_dv3 + term2_dv3 + term3_dv3;
    }
}