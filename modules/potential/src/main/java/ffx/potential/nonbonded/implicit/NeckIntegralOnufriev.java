package ffx.potential.nonbonded.implicit;

import java.util.ArrayList;
import java.util.Arrays;

public class NeckIntegralOnufriev {
    // Inputs: radius of atom i (rho i) and radius of atom j (rho j)
    // Outputs: Aij and Bij (interpolated/extrapolated where necessary)

    private static final ArrayList<Double> rhoiRows = new ArrayList<>(Arrays.asList(1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8));
    private static final ArrayList<Double> rhojColumns = new ArrayList<>(Arrays.asList(1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8));
    private static final double[][] AijAguilarOnufriev = {
            {0.0000947808,0.0001137020,0.0001001180,0.0001195900,0.0001050630,0.0000922010,0.0001099390,0.0000961817,0.0001143420,0.0000999772,0.0001184650,0.0001035490,0.0001227270},
            {0.0000789821,0.0000948444,0.0000836245,0.0000733297,0.0000878314,0.0000769799,0.0000919335,0.0000804789,0.0000958057,0.0000839558,0.0000997966,0.0000870743,0.0001031750},
            {0.0000665366,0.0000799715,0.0000703720,0.0000842561,0.0000740921,0.0000649757,0.0000775663,0.0000681323,0.0000811720,0.0000708781,0.0000842888,0.0000735226,0.0000872208},
            {0.0000560640,0.0000673256,0.0000593409,0.0000523516,0.0000625504,0.0000549136,0.0000656357,0.0000782465,0.0000685443,0.0000599304,0.0000712967,0.0000622049,0.0000738641},
            {0.0000476864,0.0000573346,0.0000505739,0.0000606430,0.0000534346,0.0000638940,0.0000560721,0.0000668434,0.0000585494,0.0000511473,0.0000608988,0.0000531489,0.0000632062},
            {0.0000406501,0.0000489724,0.0000432430,0.0000519647,0.0000456150,0.0000545524,0.0000478982,0.0000571349,0.0000500254,0.0000437652,0.0000521263,0.0000455330,0.0000540793},
            {0.0000349541,0.0000421574,0.0000370920,0.0000445857,0.0000391597,0.0000469342,0.0000411581,0.0000491397,0.0000430438,0.0000512731,0.0000448652,0.0000392138,0.0000466779},
            {0.0000411296,0.0000362900,0.0000319527,0.0000384631,0.0000337429,0.0000404799,0.0000355183,0.0000424310,0.0000371662,0.0000443414,0.0000388286,0.0000461746,0.0000403668},
            {0.0000355236,0.0000313722,0.0000276089,0.0000332680,0.0000292288,0.0000350426,0.0000307676,0.0000367966,0.0000323063,0.0000385546,0.0000337092,0.0000400815,0.0000350018},
            {0.0000309168,0.0000272940,0.0000240680,0.0000289772,0.0000254577,0.0000305399,0.0000268793,0.0000321616,0.0000281553,0.0000336062,0.0000294027,0.0000349824,0.0000305618},
            {0.0000269536,0.0000237946,0.0000287199,0.0000252553,0.0000222642,0.0000266807,0.0000234688,0.0000280674,0.0000245787,0.0000293847,0.0000256658,0.0000305777,0.0000267274},
            {0.0000236280,0.0000208424,0.0000251832,0.0000221566,0.0000195124,0.0000234327,0.0000205702,0.0000246477,0.0000215557,0.0000257936,0.0000225548,0.0000268717,0.0000234796},
            {0.0000207603,0.0000182939,0.0000161523,0.0000195050,0.0000234519,0.0000206073,0.0000180540,0.0000216798,0.0000189606,0.0000226817,0.0000198338,0.0000236529,0.0000206561}
    };
    private static final double[][] BijAguilarOnufriev = {
            {0.00,0.15,0.10,0.25,0.20,0.15,0.30,0.25,0.40,0.35,0.50,0.45,0.60},
            {0.05,0.20,0.15,0.10,0.25,0.20,0.35,0.30,0.45,0.40,0.55,0.50,0.65},
            {0.10,0.25,0.20,0.35,0.30,0.25,0.40,0.35,0.50,0.45,0.60,0.55,0.70},
            {0.15,0.30,0.25,0.20,0.35,0.30,0.45,0.60,0.55,0.50,0.65,0.60,0.75},
            {0.20,0.35,0.30,0.45,0.40,0.55,0.50,0.65,0.60,0.55,0.70,0.65,0.80},
            {0.25,0.40,0.35,0.50,0.45,0.60,0.55,0.70,0.65,0.60,0.75,0.70,0.85},
            {0.30,0.45,0.40,0.55,0.50,0.65,0.60,0.75,0.70,0.85,0.80,0.75,0.90},
            {0.55,0.50,0.45,0.60,0.55,0.70,0.65,0.80,0.75,0.90,0.85,1.00,0.95},
            {0.60,0.55,0.50,0.65,0.60,0.75,0.70,0.85,0.80,0.95,0.90,1.05,1.00},
            {0.65,0.60,0.55,0.70,0.65,0.80,0.75,0.90,0.85,1.00,0.95,1.10,1.05},
            {0.70,0.65,0.80,0.75,0.70,0.85,0.80,0.95,0.90,1.05,1.00,1.15,1.10},
            {0.75,0.70,0.85,0.80,0.75,0.90,0.85,1.00,0.95,1.10,1.05,1.20,1.15},
            {0.80,0.75,0.70,0.85,1.00,0.95,0.90,1.05,1.00,1.15,1.10,1.25,1.20}
    };

    /**
     * NeckIntegralOnufrievConstants Static Class
     *
     */
    public static class NeckIntegralOnufrievConstants {
        private double Aij;
        private double Bij;

        public double getAij() {
            return this.Aij;
        }

        public double getBij() {
            return this.Bij;
        }

        public static double[] run(double rhoi, double rhoj) {

            double Aij = 0.0;
            double Bij = 0.0;

            // If both rho i and rho j are in Aguilar/Onufriev data, get Aij and Bij values directly from tables
            if (rhoiRows.contains(rhoi) && rhojColumns.contains(rhoj)) {
                int row = rhoiRows.indexOf(rhoi);
                int col = rhojColumns.indexOf(rhoj);
                Aij = AijAguilarOnufriev[row][col];
                Bij = BijAguilarOnufriev[row][col];
            } else {
                // Otherwise, interpolate/extrapolate as needed
                boolean calculatei = false;
                boolean counti = true;
                boolean calculatej = false;
                boolean countj = true;
                // Find which two values of rho i and rho j the inputs fall between
                int lowi = 0;
                int lowj = 0;
                int highi = lowi + 1;
                int highj = lowj + 1;

                if (!rhoiRows.contains(rhoi)) {
                    calculatei = true;
                    // If input rho i is smaller than all values in table, extrapolate down using first two table values
                    // These values are the defaults (set above)
                    if (rhoiRows.get(0) > rhoi) {
                        counti = false;
                    }
                    // If input rho i is larger than all values in table, extrapolate up using last two table values
                    if (rhoiRows.get(rhoiRows.size() - 1) < rhoi) {
                        lowi = rhoiRows.size() - 2;
                        highi = rhoiRows.size() - 1;
                        counti = false;
                    }
                    while (counti) {
                        // Find the two table values that the input rho i falls between
                        if (rhoiRows.get(lowi) < rhoi && rhoi < rhoiRows.get(lowi + 1)) {
                            highi = lowi + 1;
                            counti = false;
                        } else {
                            lowi++;
                        }
                        if (lowi >= rhoiRows.size()) {
                            counti = false;
                        }
                    }
                }

                if (!rhojColumns.contains(rhoj)) {
                    calculatej = true;
                    // If input rho j is smaller than all values in table, extrapolate down using first two table values
                    // These values are the defaults (set above)
                    if (rhojColumns.get(0) > rhoj) {
                        countj = false;
                    }
                    // If input j is larger than all values in table, extrapolate up using last two table values
                    if (rhojColumns.get(rhojColumns.size() - 1) < rhoj) {
                        lowj = rhojColumns.size() - 2;
                        highj = rhojColumns.size() - 1;
                        countj = false;
                    }
                    while (countj) {
                        // Find the two table values that the input rho j falls between
                        if (rhojColumns.get(lowj) < rhoj && rhoj < rhojColumns.get(lowj + 1)) {
                            highj = lowj + 1;
                            countj = false;
                        } else {
                            lowj++;
                        }
                        if (lowj >= rhojColumns.size()) {
                            countj = false;
                        }
                    }
                }

                // If the values of rho i and rho j aren't in table, interpolate/extrapolate
                double startInterp_i = rhoiRows.get(lowi);
                double endInterp_i = rhoiRows.get(highi);
                double startInterp_j = rhojColumns.get(lowj);
                double endInterp_j = rhojColumns.get(highj);

                if (calculatei && calculatej) {
                    // Rho i and rho j aren't table values
                    Aij = interpolateAij(startInterp_i, endInterp_i, startInterp_j, endInterp_j, rhoi, rhoj);
                    Bij = interpolateBij(startInterp_i, endInterp_i, startInterp_j, endInterp_j, rhoi, rhoj);
                }
                // Rho i is in table, but rho j isn't
                if (!calculatei && calculatej) {
                    Aij = interpolateAij(rhoi, rhoi, startInterp_j, endInterp_j, rhoi, rhoj);
                    Bij = interpolateBij(rhoi, rhoi, startInterp_j, endInterp_j, rhoi, rhoj);
                }
                // Rho i isn't in table, but rho j is
                if (calculatei && !calculatej) {
                    Aij = interpolateAij(startInterp_i, endInterp_i, rhoj, rhoj, rhoi, rhoj);
                    Bij = interpolateBij(startInterp_i, endInterp_i, rhoj, rhoj, rhoi, rhoj);
                }
            }

            // Never let Aij or Bij be negative
            if(Aij < 0.00){ Aij = 0.0; }
            if(Bij < 0.00){ Bij = 0.0; }

            return new double[]{Aij, Bij};
        }
    }

    private static double interpolateAij(double startInterp_i, double endInterp_i, double startInterp_j, double endInterp_j, double rhoi, double rhoj){

        double Aij;
        if(startInterp_i == endInterp_i){
            // 1D interpolation: Only interpolate between rho j values
            // System.out.println("Start j: "+startInterp_j+" End j: "+endInterp_j+" Start i: "+startInterp_i+" End i: "+endInterp_i);
            Aij = interpolate1D(startInterp_j, endInterp_j, rhoj,
                    AijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    AijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);

        } else if(startInterp_j == endInterp_j){
            // 1D interpolation: Only interpolate between rho i values
            // System.out.println("Start j: "+startInterp_j+" End j: "+endInterp_j+" Start i: "+startInterp_i+" End i: "+endInterp_i);
            Aij = interpolate1D(startInterp_i, endInterp_i, rhoi,
                    AijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    AijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);
        } else{
            // 2D interpolation: Interpolate both values
            // System.out.println("Start j: "+startInterp_j+" End j: "+endInterp_j+" Start i: "+startInterp_i+" End i: "+endInterp_i);
            Aij = interpolate2D(startInterp_i, endInterp_i, startInterp_j, endInterp_j,rhoi, rhoj,
                    AijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    AijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    AijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(endInterp_j)],
                    AijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);
        }

        // Test capping Aij at table values


        return Aij;
    }

    private static double interpolateBij(double startInterp_i, double endInterp_i, double startInterp_j, double endInterp_j, double rhoi, double rhoj){

        double Bij;

        if(startInterp_i == endInterp_i){
            // 1D interpolation: Only interpolate between rho j values
            Bij = interpolate1D(startInterp_j, endInterp_j, rhoj,
                    BijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    BijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);
        } else if(startInterp_j == endInterp_j){
            // 1D interpolation: Only interpolate between rho i values
            Bij = interpolate1D(startInterp_i, endInterp_i, rhoi,
                    BijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    BijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);
        } else{
            // 2D interpolation: Interpolate both values
            Bij = interpolate2D(startInterp_i, endInterp_i, startInterp_j, endInterp_j, rhoi, rhoj,
                    BijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    BijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(startInterp_j)],
                    BijAguilarOnufriev[rhoiRows.indexOf(startInterp_i)][rhojColumns.indexOf(endInterp_j)],
                    BijAguilarOnufriev[rhoiRows.indexOf(endInterp_i)][rhojColumns.indexOf(endInterp_j)]);
        }
        return Bij;
    }

    private static double interpolate1D(double y1, double y2, double y, double fxy1, double fxy2){
        double frac1 = (y2 - y)/(y2 - y1);
        double frac2 = (y - y1)/(y2 - y1);
        double product1 = frac1 * fxy1;
        double product2 = frac2 * fxy2;

        return product1 + product2;
    }

    private static double interpolate2D(double x1, double x2, double y1, double y2, double x, double y,
                                 double fx1y1, double fx2y1, double fx1y2, double fx2y2){
        double fxy = 0.0;
        double fxy1 = (x2 - x)/(x2 -x1) * fx1y1 + (x - x1)/(x2 - x1) * fx2y1;
        double fxy2 = (x2 -x)/(x2 - x1) * fx1y2 + (x - x1)/(x2 - x1) * fx2y2;

        fxy = (y2 - y)/(y2 - y1) * fxy1 + (y - y1)/(y2 - y1) * fxy2;

        return fxy;
    }

}
