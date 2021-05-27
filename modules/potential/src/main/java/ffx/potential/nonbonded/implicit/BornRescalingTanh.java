package ffx.potential.nonbonded.implicit;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.PI;

public class BornRescalingTanh {
    private static final Logger logger = Logger.getLogger(BornRescalingTanh.class.getName());

    // Implements tanh function for Born radii rescaling
    // A couple of different forms described in Mongan 2007 and Aguliar/Onufriev 2010

    // Ri^-1 = (rhoi^-3 - rhoi^-3*tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3)^(1/3)
    // Psi = (3/4*pi)*int(|r|^-6 dV)
    // beta parameters will be fit using an optimizer

    public static final double MaxBornRadius = 50.0;
    /**
     * 1/50^3 where 50 Angstroms is the maximum Born radius
     */
    private static final double recipMaxBornRadius3 = 1.0 / Math.pow(MaxBornRadius, 3.0);
    private static final double PI4_3 = 4.0 / 3.0 * PI;
    private static final double oneThird = 1.0 / 3.0;

    // Set defaults to Aguliar/Onufriev original numbers: beta0 = 1.0 beta1 = 18.4377 beta2 = 343.7171
    private static double beta0 = 0.4720;
    private static double beta1 = 1.0252;
    private static double beta2 = 0.2488;

    // Pass in Psi and rhoi for rescaling
    public static double rescale(double Ii, double rhoi) {

        // tanh function:
        // Ri_inverse = ((1/rhoi^3) - ((1/rhoi^3)-(1/50^3))*tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3))^1/3
        // Set up tanh function components

        double rhoi3 = rhoi * rhoi * rhoi;
        // Here, Psi = Ii, which is passed in
        double rhoi3Psi = rhoi3 * Math.abs(Ii);
        double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
        double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
        // If the output of the tanh function is 1.0, then the Born radius will be MaxBornRadius
        double tanh_constant = PI4_3*((1.0 / rhoi3) - recipMaxBornRadius3);

        double input = beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3;
        double output = Math.copySign(Math.tanh(input),Ii);
        //logger.info("Ii: "+Ii+" input: "+input+" output: "+output+" returned: "+tanh_constant*output);
        return tanh_constant * output;
    }

    public static double getTanhContribution(double Ii, double rhoi){

        double rhoi3 = rhoi * rhoi * rhoi;
        // Here, Psi = Ii, which is passed in
        double rhoi3Psi = rhoi3 * Math.abs(Ii);
        double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
        double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;

        double rhoi6Psi = rhoi3 * rhoi3 * Math.abs(Ii);
        double rhoi9Psi2 = rhoi6Psi2 * rhoi3;
        // If the output of the tanh function is 1.0, then the Born radius will be MaxBornRadius
        double tanh_constant = PI4_3 * ((1.0 / rhoi3) - recipMaxBornRadius3);

        double input = beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3;
        double chainRuleTerm = beta0 * rhoi3 - 2 * beta1 * rhoi6Psi + 3 * beta2 * rhoi9Psi2;

        return tanh_constant * chainRuleTerm * (1 - Math.pow(Math.tanh(input),2));
    }

    // Getters and setters for optimizing constants
    public static double getBeta0() {
        return beta0;
    }

    public static double getBeta1() {
        return beta1;
    }

    public static double getBeta2() {
        return beta2;
    }

    public static void setBeta0(double beta0_input) {
        beta0 = beta0_input;
    }

    public static void setBeta1(double beta1_input) {
        beta1 = beta1_input;
    }

    public static void setBeta2(double beta2_input) {
        beta2 = beta2_input;
    }

}
