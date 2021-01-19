package ffx.potential.nonbonded.implicit;

public class BornRescalingTanh {

    // Implements tanh function for Born radii rescaling
    // A couple of different forms described in Mongan 2007 and Aguliar/Onufriev 2010

    // Ri^-1 = (rhoi^-3 - rhoi^-3*tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3)^(1/3)
    // Psi = (3/4*pi)*int(|r|^-6 dV)
    // beta parameters will be fit using an optimizer

    public static class Tanh {
        // Set defaults to Aguliar/Onufriev original numbers
        private static double beta0 = 0.0;
        private static double beta1 = 18.4377;
        private static double beta2 = 313.7171;

        // Getters and setters for optimizing constants
        public double getBeta0() { return beta0; }

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
        public static double rescale(double Ri, double rhoi) {

            boolean localBeta = true;
            // Set up tanh function components
            double Psi = 1.0 / rhoi - 1.0 / Ri;
            double rhoi3 = rhoi * rhoi * rhoi;
            double rhoi3Psi = rhoi3 * Psi;
            double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
            double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
            double Imax = 1.0 / rhoi3;
            //double ci = 1.0 - rhoi3 / (Ri*Ri*Ri);
            double b0 = 1.0/(1.0 - rhoi3);
            double tanh = Ri;

            // In Aguilar/Onufriev 2010, beta0 is set based on rhoi and A (electrostatic size of the molecule)
            // If A = 1.0, then ci = 1 - (1/rhoi^-3) and beta0 is set this way
            // If we decide to use this radii-dependent beta0, we won't optimize that parameter
            if(localBeta) {
                tanh = Math.tanh(b0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3);
            } else {
                tanh = Math.tanh(beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3);
            }
            double rescaled3 = Imax - Imax*tanh;
            double onethird = 1.0 / 3.0;
            //System.out.println("Rescaled radius: " + Math.pow(rescaled3, onethird));

            // Return the tanh rescaled Born radius
            // Original tanh function is for Ri^-1
            return 1.0/Math.pow(rescaled3, onethird);
        }
    }
}
