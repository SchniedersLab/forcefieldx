/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

/**
 *
 * @author Michael J. Schnieders, Wei Yang and Pengyu Ren
 */
public class OSRW {

    /**
     * State variable lambda ranges from 0.0 .. 1.0.
     */
    private double lambda;
    /**
     * Flag to indicate that the total gradient dE/dL should be calculated.
     */
    private boolean lambdaGradient = false;
    /**
     * Flag to indicate that the number of energy evaluations is currently 
     * being counted.
     */
    private boolean doCounting = true;
    /**
     * Number of times the OSRW biasing potential has been evaluated with the
     * "doCounting" flag true.
     */
    private int energyCount;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005).
     * The final Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     * 
     * With this scheme, the maximum of biasing Gaussians
     * is at the edges. 
     */
    private int lambdaBins = 101;
    /**
     * It is useful to have an odd number of bins, so that there is
     * a bin from FL=-dFL/2 to dFL/2 so that as FL approaches zero its
     * contribution to thermodynamic integration goes to zero. Otherwise 
     * a contribution of zero from a L bin can only result from equal
     * sampling of the ranges -dFL to 0 and 0 to dFL.
     */
    private int FLambdaBins = 401;
    /**
     * The recursion kernel stores the number of visits to 
     * each [lambda][Flambda] bin.
     */
    private int recursionKernel[][];
    /**
     * When evaluating the biasing potential, contributions from Gaussians
     * centered more the "biasCutoff" away will be neglected.
     */
    private int biasCutoff = 5;
    /**
     * Width of the lambda bin.
     */
    private double dL = 0.01;
    /**
     * Half the width of a lambda bin.
     */
    private double dL_2 = dL / 2.0;
    /**
     * The width of the F_lambda bin.
     */
    private double dFL = 2.0;
    /**
     * Half the width of the F_lambda bin.
     */
    private double dFL_2 = dFL / 2.0;
    /**
     * The minimum value of the first lambda bin.  
     */
    private double minLambda = -0.005;
    /**
     * The minimum value of the first F_lambda bin.
     */
    private double minFLambda = -(dFL * FLambdaBins) / 2.0;
    /**
     * The maximum value of the last F_lambda bin.
     */
    private double maxFLambda = minFLambda + FLambdaBins * dFL;
    /**
     * Total partial derivative of the potential being sampled w.r.t. lambda.
     */
    private double dEdLambda = 0.0;
    /**
     * 2nd partial derivative of the potential being sampled w.r.t lambda.
     */
    private double d2EdLambda2 = 0.0;
    private double dUdXdL[] = null;
    private double gaussianMag = 0.005;
    private double FLambda[];
}
