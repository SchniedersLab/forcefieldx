package ffx.numerics.multipole;

import static java.lang.Math.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 *  The AmoebaPlusDampTensorGlobal class computes derivatives of damping via recursion to
 *  order &lt;= 2 for Cartesian multipoles defined in AMOEBA+. These are the nuclear charge
 *  interactions with permanent multipoles. The zeroth order damp term is 1-exp(-alpha*r).
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
public class AmoebaPlusDampTensorGlobal extends CoulombTensorGlobal {

    private double alpha;
    private double alpha2;
    private static final double oneThird = 1.0 / 3.0;

    /**
     * Constructor for CoulombTensorGlobal.
     *
     * @param order The order of the tensor.
     */
    public AmoebaPlusDampTensorGlobal(int order, double alpha1, double alpha2) {
        super(order);
        this.alpha = alpha1;
        this.alpha2 = alpha2;
        this.operator = abs(alpha - alpha2) < 1e-5 ?
                Operator.AMOEBA_PLUS_SYM_DAMP_FIELD : Operator.AMOEBA_PLUS_DAMP_FIELD;
        assert (order <= 2); // Nuclear charge is a point charge
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

        // Add the damping term: edamp = 1-exp(-alpha*r).
        dampSource(alpha, R, T000);

        // Swap alphas for next tensor gen in required
        if (this.operator != Operator.AMOEBA_PLUS_SYM_DAMP_FIELD) {
            double temp = alpha;
            this.alpha = alpha2;
            this.alpha2 = temp;
        }
    }

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param alpha adjustable damping parameter
     * @param R The separation distance.
     * @param T000 Location to store the source terms.
     */
    protected static void dampSource(double alpha, double R, double[] T000) {
        // Add the damping terms: edamp = 1-exp(-alpha*r).
        double alphaR = alpha * R;
        double alphaR2 = alphaR * alphaR;
        double expAR = exp(-alphaR);

        T000[0] *= 1 - expAR;
        T000[1] *= 1 - (1 + alphaR) * expAR;
        T000[2] *= 1 - (1.0 + alphaR + oneThird*alphaR2) * expAR;
        T000[3] *= 1 - (1.0 + alphaR + .4*alphaR2 + alphaR2*alphaR/15) * expAR;
    }
}
