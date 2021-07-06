package ffx.potential.extended;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.extended.TitrationUtils.Titration;
import ffx.potential.parameters.MultipoleType;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.List;
import java.util.logging.Level;

import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static ffx.potential.parameters.MultipoleType.zeroM;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.log;

public class TautomerESV extends ExtendedVariable {
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
     * @param esvSystem a {@link ffx.potential.extended.ExtendedSystem} object.
     * @param multiRes  a {@link ffx.potential.bonded.MultiResidue} object.
     */
    public TautomerESV(ExtendedSystem esvSystem, MultiResidue multiRes) {
        super(esvSystem, multiRes, 1.0);
        pairedTitrationIndex = esvIndex - 1;
        MolecularAssembly mola = esvSystem.getMolecularAssembly();
        CompositeConfiguration properties = mola.getProperties();
        this.constPh = esvSystem.getConstantPh();
    }

    protected void pairTitrationESV(ExtendedSystem esvSystem, int pairedTitrationIndex, int esvIndex){
        TitrationESV pairedTitrationESV = (TitrationESV) esvSystem.getEsv(pairedTitrationIndex);
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
        if (!config.electrostatics) {
            return;
        }
        /* If not softcoring unshared atoms, then scale them as well
         * (with an implied zero-multipole background atom).	  */
        List<Atom> iterate = ExtUtils.joinedListView(atomsShared, atomsUnshared);

        for (Atom fg : iterate) {
            //            MultipoleType Ptype, Utype;
            Atom bg = fg2bg.get(fg);
            final MultipoleType Ptype = fg.getMultipoleType();
            final MultipoleType Utype;
            if (bg == null) {
                if (atomsUnshared.contains(fg)) {
                    Utype = new MultipoleType(zeroM, Ptype.frameAtomTypes, Ptype.frameDefinition, false);
                } else {
                    logger.warning("Error @ESV.updateMultipoles: bg null && !unshared.");
                    Utype = null;
                }
            } else {
                Utype = bg.getMultipoleType();
            }
            MultipoleType types[] = new MultipoleType[]{Ptype, Utype};
            double mWeights[], mdotWeights[];
            if (config.allowLambdaSwitch && config.nonlinearMultipoles) {
                mWeights = new double[]{lSwitch, 1.0 - lSwitch};
                mdotWeights = new double[]{dlSwitch, -dlSwitch};
            } else {
                mWeights = new double[]{lambda, 1.0 - lambda};
                mdotWeights = new double[]{1.0, -1.0};
            }
            if (Ptype == null) {
                // Multipoles not yet defined.
                continue;
            }
            StringBuilder sb = new StringBuilder();
            if (Utype == null) {
                sb.append(format("Error @ESV.updateMultipoleTypes: bgType null."));
                sb.append(format("   fg: %s, '%s'", fg.toString(), fg.getName()));
                sb.append(format(" Background atoms available for match: "));
                for (Atom debug : atomsBackground) {
                    sb.append(format("   bg: %s, '%s'", debug.toString(), debug.getName()));
                }
                logger.warning(sb.toString());
                continue;
            }
            final double[] esvMultipole;
            final double[] esvMultipoleDot;
            final MultipoleType esvType;
            final MultipoleType esvDotType;
            try {
                esvType = MultipoleType.weightMultipoleTypes(types, mWeights, Ptype.frameAtomTypes);
                esvDotType = MultipoleType.weightMultipoleTypes(types, mdotWeights, Ptype.frameAtomTypes);
            } catch (IllegalArgumentException ex) {
                logger.warning(
                        format(
                                "Multipole scaling failed for fg,bg atoms: %s %s",
                                fg.describe(Atom.Descriptions.Trim), bg.describe(Atom.Descriptions.Trim)));
                throw ex;
            }
            if (esvType == null || esvDotType == null) {
                logger.severe("Error @ESV.updateMultipoleTypes: M or Mdot null.");
            }
            if (isTitratableHydrogen(fg)) {
                final double scaledAlpha = fg.getPolarizeType().polarizability * getLambda();
                fg.setEsv(this, esvType, esvDotType, scaledAlpha);
            } else {
                fg.setEsv(this, esvType, esvDotType);
            }

            sb.append(
                    format(
                            " Assigning ESV MultipoleTypes for atom %s\n",
                            fg.describe(Atom.Descriptions.Resnum_Name)));
            sb.append(format("  U: %.2f*%s\n", lambda, Ptype.toCompactBohrString()));
            sb.append(format("  P: %.2f*%s\n", 1.0 - lambda, Utype.toCompactBohrString()));
            sb.append(format("  M:      %s\n", fg.getEsvMultipole().toCompactBohrString()));
            sb.append(format("  Mdot:   %s\n", fg.getEsvMultipoleDot().toCompactBohrString()));
            if (logger.isLoggable(Level.FINE)) {
                logger.info(sb.toString());
            }
        }
    }

    /**
     * Scales bonded terms based on lambda. Currently off.
     */
    protected void updateBondedLambdas() {
        if (!config.bonded) {
            return;
        }
        bondedDeriv.set(0.0);
        double Sl = getLambdaSwitch();
        double dSldL = getSwitchDeriv();
        for (BondedTerm bt1 : bondedFg) {
            boolean attachExtended = true;
            Sl = getLambdaSwitch();
            dSldL = getSwitchDeriv();
            Atom[] atoms = bt1.getAtomArray();
            for (Atom unshared : atomsUnshared) {
                for (Atom atom : atoms) {
                    if (atom.equals(unshared)) {
                        // logger.info(" Unshared Atom Bonded Term: " + bt1);
                        //Sl = 1.0;
                        //dSldL = 1.0;
                        attachExtended = false;
                        break;
                    }
                }
            }
            if (attachExtended) {
                if (config.decomposeBonded) {
                    //bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv, fgBondedDerivDecomp);
                    fgBondedDerivDecomp.clear();
                } else {
                    //bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv);
                }
            }
        }

        for (BondedTerm bt0 : bondedBg) {
            Sl = getLambdaSwitch();
            dSldL = getSwitchDeriv();
            Atom[] atoms = bt0.getAtomArray();
            for (Atom unshared : atomsUnshared) {
                for (Atom atom : atoms) {
                    if (atom.equals(unshared)) {
                        logger.info(" Background Unshared Atom Bonded Term: " + bt0);
                        // Sl = 0.0;
                        // dSldL = -1.0;
                        break;
                    }
                }
            }
            if (config.decomposeBonded) {
                //bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv, bgBondedDerivDecomp);
                bgBondedDerivDecomp.clear();
            } else {
                //bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv);
            }
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
