package ffx.numerics.multipole;


import edu.rit.util.Arrays;

import static org.apache.commons.math3.util.FastMath.pow;

/**
 * This class allows for the source terms of tensors to be combined, and therefore
 * multiple interaction tensors can be computed simultaneously (as a sum).
 * <p></p>
 * Non-Ex. Amoeba = M-pole*(Ewald)*M-pole + self
 * Ex. Amoeba PolarizationE = Dipole*(Coulomb - Thole)*M-pol
 * Ex. Amoeba SCF Eval = Dipole*(permPolField + (Coulomb - Thole)*Dipole) (at convergence goes to A PolE)
 * <p></p>
 * Ex. Amoeba+ = M-pol*(Ewald - CPenOverlap)*M-pol + Core*(Ewald-CPenDamp)*(Mpole) + Core*(Ewald)*Core + self
 * Ex. Amoeba+ PolE = Dipole*(Ewald - TholeDirect)*M-pol
 * Ex. Amoeba+ SCF Eval = Dipole(permPolField + (Coulomb - Thole)*Dipole) (at convergence this goes to A+ PolE)
 */
public class CombinedTensorGlobal extends CoulombTensorGlobal {

    private double[] termsToAdd;
    private double multiplier = 1.0;
    private double[] termsToAddSep;

    /**
     * Constructor for CoulombTensorGlobal.
     *
     * @param order The order of the tensor.
     */
    public CombinedTensorGlobal(int order) {
        super(order);
        this.termsToAdd = new double[order+1];
        this.termsToAddSep = new double[order+1];
    }

    /**
     * Accumulates onto existing terms. Assumes we add all the elements
     * onto the respective index on existing array. Implicitly adds coulomb
     * interaction.
     * @param terms
     */
    public void addTerms(double[] terms) {
        assert(terms.length <= this.termsToAdd.length);
        for(int i = 0; i < terms.length; i++) {
            this.termsToAdd[i] += terms[i];
        }
    }

    public void resetSource(){
        this.termsToAdd = new double[this.termsToAdd.length];
        this.termsToAddSep = new double[this.termsToAdd.length];
    }

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param T000 Location to store the source terms.
     */
    @Override
    protected void source(double[] T000) {
        // Compute the normal Coulomb auxiliary term.
        super.source(T000); // Important that what is passed into here does not already have coulomb terms in it

        for(int i = 0; i < termsToAdd.length; i++) {
            T000[i] *= 1*multiplier + termsToAdd[i]; // multiplier may make coulomb term negative
            T000[i] += termsToAddSep[i];
        }
    }

    public void addCoulombMultiplier(double mult) {
        this.multiplier = mult;
    }

    public void addTermsSeparate(double[] terms) {
        assert(terms.length <= this.termsToAddSep.length);
        for(int i = 0; i < terms.length; i++) {
            this.termsToAddSep[i] += terms[i];
        }
    }
}
