package ffx.numerics.multipole;

import static java.lang.Math.abs;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 *  The AmoebaPlusDampTensorGlobal class computes derivatives of overlap via recursion to
 *  order &lt;= 6 for Cartesian multipoles defined in AMOEBA+. These are the multipole overlap
 *  interactions. The zeroth order overlap term is 1-A*exp(-alpha1*r)-B*exp(-alpha2*r).
 *
 *  @author Matthew J. Speranza
 *  @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 *      Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 *      computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 *      Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @see <a href="https://doi.org/10.1039/C6CP06017J" target="_blank"> Rackers, Wang, Liu,
 *      Piquemal, Ren and Ponder, An optimized charge penetration model for use with the AMOEBA
 *      force field. Phys. Chem. Chem. Phys., 2017,19, 276-291. </a>
 *  @since 1.0
 */
public class AmoebaPlusOverlapTensorGlobal extends CoulombTensorGlobal {

    private double aI;
    private double aK;
    private double A;
    private double B;


    /**
     * Constructor for CoulombTensorGlobal.
     *
     * @param order The order of the tensor.
     */
    public AmoebaPlusOverlapTensorGlobal(int order, double alphaI, double alphaK) {
        super(order);
        this.operator = Operator.AMOEBA_PLUS_OVERLAP_FIELD;
        this.aI = alphaI;
        this.aK = alphaK;
        this.A = aK*aK / ((aK*aK) - (aI * aI));
        this.B = aI*aI / ((aI*aI) - (aK * aK));
        // Multipole-Multipole overlap only defined up to order 6 (OST - Quadrupole-Quadrupole)
        assert (order <= 6);
    }

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param T000 Location to store the source terms.
     */
    @Override
    protected void source(double[] T000) {
        // Compute the normal Coulomb auxiliary term.
        super.source(T000);

        // Add the Thole damping terms: edamp = exp(-thole*u^3).
        overlapSource(aI, aK, A, B, R, T000);
    }

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param R The separation distance.
     * @param T000 Location to store the source terms.
     */
    protected static void overlapSource(double aI, double aK, double A, double B, double R, double[] T000) {
        // Add the damping terms: edamp = 1-A*exp(-alphaI*r)-B*exp(-alphaK*r).
        // 1/r^11 terms from tinker code & are not included in supp info of paper
        if (abs(aI-aK) < 1e-7) { // ai = ak
            double alpha = aI;
            double alphaR = alpha * R;
            double alphaR2 = alphaR * alphaR;
            double alphaR4 = alphaR2 * alphaR2;
            double expAR = exp(-alphaR);

            T000[0] *= 1 - (1 + alphaR/2.0) * expAR; // 1/r
            T000[1] *= 1 - (1 + alphaR + alphaR2 / 2.0) * expAR; // 1/r^3
            T000[2] *= 1 - (1.0 + alphaR + alphaR2 / 2.0 + alphaR2 * alphaR / 6.0) * expAR; // 1/r^5
            T000[3] *= 1 - (1.0 + alphaR + alphaR2 / 2.0 + alphaR2 * alphaR / 6.0 +
                    alphaR2 * alphaR2 / 30.0) * expAR; // 1/r^7
            T000[4] *= 1 - (1.0 + alphaR + alphaR2 / 2.0 + alphaR2 * alphaR / 6.0 +
                    4.0 * alphaR4 / 105.0 + alphaR4 * alphaR / 210.0) * expAR; // 1/r^9
            T000[5] *= 1 - (1.0 + alphaR + alphaR2 * .5 + alphaR2*alphaR/6.0 +
                    5.0*alphaR4/126.0 + 2.0*alphaR4*alphaR/315.0 +
                    alphaR4*alphaR2/1890.0)*expAR; // 1/r^11

        } else { // ai != ak
            assert(A != Double.POSITIVE_INFINITY && B != Double.POSITIVE_INFINITY);
            double aIR = aI * R;
            double aIR2 = aIR * aIR;
            double aKR = aK * R;
            double aKR2 = aKR * aKR;
            double expAI = A * exp(-aIR);
            double expAK = B * exp(-aKR);

            T000[0] *= 1 - expAI - expAK; // 1/r
            T000[1] *= 1 - (1 + aIR) * expAI - (1 + aKR) * expAK; // 1/r^3
            T000[2] *= 1 - (1.0 + aIR + aIR2 / 3.0) * expAI - (1.0 + aKR + aKR2 / 3.0) * expAK; // 1/r^5
            T000[3] *= 1 - (1.0 + aIR + 2.0 * aIR2 / 5.0 + aIR * aIR2 / 15.0) * expAI - // 1/r^7
                    (1.0 + aKR + 2.0 * aKR2 / 5.0 + aKR * aKR2 / 15.0) * expAK;
            T000[4] *= 1 - (1.0 + aIR + 3.0 * aIR2 / 7.0 + 2.0 * aIR * aIR2 / 21.0 + aIR2 * aIR2 / 105.0) * expAI - // 1/r^9
                    (1.0 + aKR + 3.0 * aKR2 / 7.0 + 2.0 * aKR * aKR2 / 21.0 + aKR2 * aKR2 / 105.0) * expAK;
            T000[5] *= 1.0 - (1.0 + aIR + 4.0*aIR2/9.0 + aIR2*aI/9.0 +
                             aIR2*aIR2/63.0 + aIR2*aIR2*aIR/945.0)*expAI
                    - (1.0 + aKR + 4.0*aKR2/9.0 + aKR2*aKR/9.0 + aKR2*aKR2/63.0 + aKR2*aKR2*aKR/945.0)*expAK; // TODO: Get 6th order terms
        }
    }
}
