package ffx.potential.nonbonded.implicit;

public class BornRescalingTanh {

    // Implements tanh function for Born radii rescaling
    // A couple of different forms described in Mongan 2007 and Aguliar/Onufriev 2010

    // Ri^-1 = (rhoi^-3 - rhoi^-3*tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3)^(1/3)
    // Psi = (3/4*pi)*int(|r|^-6 dV)
    // beta parameters will be fit using an optimizer

    public static final double MaxBornRadius = 50.0;
    /**
     * 1/50^3 where 50 Angstroms is the maximum Born radius
     */
    private static final double recipMaxBornRadius3 = 1.0 / Math.pow(MaxBornRadius,3.0);

    //public static class Tanh {
        // Set defaults to Aguliar/Onufriev original numbers: beta1 = 18.4377 beta2 = 343.7171
        private static double beta0 = 1.0;
        private static double beta1 = 0.8;//1.0508;
        private static double beta2 = 2.75;//4.6766;
        // Getters and setters for optimizing constants
        public double getBeta0() {
            return beta0;
        }

        public double getBeta1() {
            return beta1;
        }

        public double getBeta2() {
            return beta2;
        }

        public static void setBeta0(double beta0) {
            beta0 = beta0;
        }

        public static void setBeta1(double beta1) {
            beta1 = beta1;
        }

        public static void setBeta2(double beta2) {
            beta2 = beta2;
        }

        // Pass in Psi and rhoi for rescaling
        public static double rescale(double Ii, double rhoi) {

            // tanh function:
            // Ri_inverse = ((1/rhoi^3) - ((1/rhoi^3)-(1/50^3))*tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3))^1/3
            // Set up tanh function components
            double rhoi3 = rhoi * rhoi * rhoi;
            // Here, Psi = Ii, which is passed in
            double rhoi3Psi = rhoi3 * Ii;
            double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
            double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
            double tanh_constant = (1.0 / rhoi3) - recipMaxBornRadius3;

            // In Aguilar/Onufriev 2010, beta0 is set based on rhoi and A (electrostatic size of the molecule)
            // If A = 1.0, then ci = 1 - (1/rhoi^-3) and beta0 is set this way
            // If we decide to use this radii-dependent beta0, we won't optimize that parameter
            return tanh_constant * Math.tanh(beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3);

        }
   // }
}
