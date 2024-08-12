package ffx.numerics.multipole;

import java.util.Arrays;

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

    private double beta = 0.0;
    private double alpha;
    private double alpha2;
    private static final double oneThird = 1.0 / 3.0;
    private double[] ewaldSource;

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
        assert (order <= 3); // Nuclear charge is a point charge
    }

    public AmoebaPlusDampTensorGlobal(int order, double alpha1, double alpha2, double ewaldA) {
        this(order, alpha1, alpha2);
        this.beta = ewaldA;
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
        double[] copy = Arrays.copyOf(T000, o1);
        // Add the damping term: edamp = 1-exp(-alpha*r).
        dampSource(alpha, R, T000);
        if(beta > 1e-3){
            this.ewaldSource = new double[this.order+1];
            EwaldTensorGlobal.fillEwaldSource(this.order, beta,
                    EwaldTensorGlobal.initEwaldSource(this.order, beta, new double[o1]),
                    this.R, ewaldSource);
            // T000 = Ewald - Coulomb + Core
            for(int i = 0; i < ewaldSource.length; i++){ T000[i] += ewaldSource[i] - copy[i]; }
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


    /**
     * Terms 1, 2, 3 in Eq. 5 of AMOEBA+ paper. Uses a swap of alpha -> alpha2
     * that takes place in the first call to the source method to generate a new
     * tensor with the second alpha.
     * @param mI
     * @param mK
     * @return
     */
    public double coreInteraction(PolarizableMultipole mI, PolarizableMultipole mK) {
        // Coulomb energy of core charges (No damping)
        double energy = beta >= 1e-3 ? mK.Z*mI.Z * ewaldSource[0] : mK.Z*mI.Z/R;

        // Cores contract with multipole moments
        multipoleIPotentialAtK(mI, 0);
        energy += mK.Z * E000;
        if (this.operator != Operator.AMOEBA_PLUS_SYM_DAMP_FIELD){
            // Generate tensor with alpha2 for term 2 -> Zi * T(damp)ij * Mj
            double temp = this.alpha;
            this.alpha = this.alpha2;
            this.alpha2 = temp;
            this.generateTensor();
            temp = this.alpha;
            this.alpha = this.alpha2;
            this.alpha2 = temp;
        }
        multipoleKPotentialAtI(mK, 0);
        energy += mI.Z * E000;

        return energy;
    }

    public double coreInteractionAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
                                             double[] Gi, double[] Gk){
        // Coulomb energy of core charges (No damping -> cant use tensor terms)
        double energy = this.beta >= 1e-3 ?  mK.Z*mI.Z * this.ewaldSource[0] : mK.Z*mI.Z/R;
        double tensorElement = this.beta >= 1e-3 ? this.ewaldSource[1]: -pow(R,-3);
        Gk[0] = mK.Z * mI.Z * x * tensorElement;
        Gk[1] = mK.Z * mI.Z * y * tensorElement;
        Gk[2] = mK.Z * mI.Z * z * tensorElement;

        // Cores contract with multipole moments
        multipoleIPotentialAtK(mI, 1);
        energy += mK.Z * E000;
        Gk[0] += mK.Z * E100;
        Gk[1] += mK.Z * E010;
        Gk[2] += mK.Z * E001;
        if (this.operator != Operator.AMOEBA_PLUS_SYM_DAMP_FIELD){
            // Generate tensor with alpha2 for term 2 -> Zi * T(damp)ij * Mj
            // R = |Rj - Ri| = Rji --> T(damp)ij != T(damp)ji after first order?
            double temp = this.alpha;
            this.alpha = this.alpha2;
            this.alpha2 = temp;
            this.generateTensor();
            temp = this.alpha;
            this.alpha = this.alpha2;
            this.alpha2 = temp;
        }
        multipoleKPotentialAtI(mK, 1);
        energy += mI.Z * E000;
        // mI.Z * EXXX = derivative Gi[X], so flip sign
        Gk[0] -= mI.Z * E100;
        Gk[1] -= mI.Z * E010;
        Gk[2] -= mI.Z * E001;

        Gi[0] = -Gk[0];
        Gi[1] = -Gk[1];
        Gi[2] = -Gk[2];

        return energy;
    }

}
