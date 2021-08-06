package ffx.potential.extended;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.MultipoleType;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.List;
import java.util.logging.Level;

import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static ffx.potential.parameters.MultipoleType.zeroM;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.log;

public class NewTautomerESV extends NewExtendedVariable {
    /**
     * Simulation pH.
     */
    private final double constPh;

    private static final double LOG10 = log(10.0);

    private int pairedTitrationIndex = 0;

    public double[] pairedTitrationLambda = new double[]{1.0};

    public double[] tautomerLambda = new double[]{1.0};

    /**
     * Constructor for TautomerESV.
     *
     * @param esvSystem a {@link ExtendedSystem} object.
     * @param residue  a {@link MultiResidue} object.
     */
    public NewTautomerESV(NewExtendedSystem esvSystem, Residue residue) {
        super(esvSystem, residue, 1.0);
        pairedTitrationIndex = esvIndex - 1;
        MolecularAssembly mola = esvSystem.getMolecularAssembly();
        CompositeConfiguration properties = mola.getProperties();
        this.constPh = esvSystem.getConstantPh();
        this.pairTitrationESV(esvSystem, pairedTitrationIndex, this.getEsvIndex());
    }

    protected void pairTitrationESV(NewExtendedSystem esvSystem, int pairedTitrationIndex, int esvIndex){
        NewTitrationESV pairedTitrationESV = (NewTitrationESV) esvSystem.getEsv(pairedTitrationIndex);
        tautomerLambda = pairedTitrationESV.pairedTautomerLambda;
        pairedTitrationESV.pairTautomerESV(esvSystem, esvIndex);

    }

    @Override
    protected void updateLambda(double lambda, boolean updateComponents) {
        tautomerLambda[0] = lambda;
        super.updateLambda(lambda, updateComponents);
        logger.info(format("Paired titration ESV: %f", pairedTitrationLambda[0]));
    }

    /**
     * Invoked by ExtendedSystem after lambda changes and by PME after multipole rotation.
     * TODO: Currently a copy of TitrationESV updateMultipoleTypes--Must implement tautomer type
     */
    protected final void updateMultipoleTypes() {
        int atomIndex;
        MultipoleType atomMultipoleType;
        double[] esvMultipole;
        double[] esvTautomerDotMultipole;
        MultipoleType esvType;
        MultipoleType esvTautomerDotType;
        for (Atom atom : atomsExtended){
            atomIndex = atom.getArrayIndex();
            atomMultipoleType = atom.getMultipoleType();
            esvMultipole = constantPhUtils.getMultipole(aminoAcid3, atomIndex, pairedTitrationLambda[0],
                    tautomerLambda[0], atomMultipoleType.getMultipole());
            esvTautomerDotMultipole = constantPhUtils.getMultipoleTautomerDeriv(aminoAcid3, atomIndex,
                    pairedTitrationLambda[0], tautomerLambda[0], atomMultipoleType.getMultipole());
            esvType = new MultipoleType(esvMultipole,atomMultipoleType.frameAtomTypes,
                    atomMultipoleType.frameDefinition, false);
            esvTautomerDotType = new MultipoleType(esvTautomerDotMultipole, atomMultipoleType.frameAtomTypes,
                    atomMultipoleType.frameDefinition, false);
            //TODO: Detect hydrogen and scale alpha.
            atom.setEsv(this, esvType, esvTautomerDotType);
        }
    }

    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     *
     * @param temperature a double.
     * @param print       a boolean.
     * @return a double.
     */
    protected double getTotalBias(double temperature, boolean print){
        double eDiscr = getDiscrBias();
        return (eDiscr);
    }

    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     *
     * @param temperature a double.
     * @param print       a boolean.
     * @return a double.
     */
    protected double getTotalBiasDeriv(double temperature, boolean print){
        double dDiscr = getDiscrBiasDeriv();
        return (dDiscr);
    }
}
