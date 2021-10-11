// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.potential.extended;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
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

/**
 * An extended system variable that allows continuous fractional protonation of an amino acid. All
 * atomic charges and bonded terms scale linearly between prot and deprot states.
 *
 * <p>Possible expansions: (1) Add back the ability to interact with OST lambda and thereby combine
 * with protein design (QuadTop). (2) Allow triple-state systems such as histidine with 0.5 protons
 * per nitrogen or tautomeric ASP/GLU. (3) Assess bytecode-output implementation via eg. ASM.
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public final class TitrationESV extends ExtendedVariable {

    /**
     * Model PMF coefficient A in A*(lambda-B)^2
     */
    private final double referenceEnergy;
    /**
     *  Model PMF coefficient B in A*(lambda-B)^2
     */
    private final double lambdaIntercept;
    /**
     * Simulation pH.
     */
    private final double constPh;
    /**
     * Reference pKa value.
     */
    private final double pKaModel;

    private static final double LOG10 = log(10.0);

    /**
     *
     */
    public double[] pairedTautomerLambda = new double[]{1.0};

    public double[] titrationLambda = new double[]{1.0};

    /**
     * Constructor for TitrationESV.
     *
     * @param esvSystem a {@link ffx.potential.extended.ExtendedSystem} object.
     * @param multiRes  a {@link ffx.potential.bonded.MultiResidue} object.
     */
    public TitrationESV(ExtendedSystem esvSystem, MultiResidue multiRes) {
        super(esvSystem, multiRes, 1.0);

        MolecularAssembly mola = esvSystem.getMolecularAssembly();
        CompositeConfiguration properties = mola.getProperties();
        this.constPh = esvSystem.getConstantPh();
        Residue currentRes = multiRes.getActive();
        Titration titration = Titration.lookup(currentRes);
        AminoAcidUtils.AminoAcid3 currentAA3 = AminoAcidUtils.AminoAcid3.valueOf(currentRes.getName());
        if(currentAA3 == AminoAcidUtils.AminoAcid3.LYS || currentAA3 == AminoAcidUtils.AminoAcid3.LYD){
            this.referenceEnergy = properties.getDouble("PMF-LYS-ReferenceEnergy",titration.refEnergy);
            this.lambdaIntercept = properties.getDouble("PMF-LYS-LambdaIntercept",titration.lambdaIntercept);
            this.pKaModel = properties.getDouble("PMF-LYS-pkaModel",titration.pKa);
        }
        else{
            this.referenceEnergy = titration.refEnergy;
            this.lambdaIntercept = titration.lambdaIntercept;
            this.pKaModel = titration.pKa;
        }
    }

    /**
     * Invoked by ExtendedSystem after lambda changes and by PME after multipole rotation.
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
                sb.append("Error @ESV.updateMultipoleTypes: bgType null.");
                sb.append(format("   fg: %s, '%s'", fg.toString(), fg.getName()));
                sb.append(" Background atoms available for match: ");
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
     * {@inheritDoc}
     */
    @Override
    public String getName() {
        return format("Titr%d", esvIndex);
    }

    /**
     * Pairs with tautomer ESV by copying its 1-element titrationLambda array.
     * @param esvSystem
     * @param esvIndex
     */
    public void pairTautomerESV(ExtendedSystem esvSystem, int esvIndex){
        TautomerESV tautomerESV = (TautomerESV) esvSystem.getEsv(esvIndex);
        titrationLambda = tautomerESV.pairedTitrationLambda;
    }

    @Override
    protected void updateLambda(double lambda, boolean updateComponents) {
        titrationLambda[0] = lambda;
        super.updateLambda(lambda, updateComponents);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double getTotalBias(double temperature, boolean print) {
        double eDiscr = getDiscrBias();
        double ePh = getPhBias(temperature);
        return (eDiscr + ePh);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double getTotalBiasDeriv(double temperature, boolean print) {
        double dDiscr = getDiscrBiasDeriv();
        double dPh = getPhBiasDeriv(temperature);
        return (dDiscr + dPh);
    }

    /**
     * Eqs 5,6 from Wallace+Shen 2011 "Continuous constant pH M.D. in explicit..." U_pH(ldh) =
     * log(10)*kb*T*(pKa_model - pH)*ldh U_mod(ldh) = potential of mean force for protonation (or
     * -deprot) of model compound U_star = sum(ldh) { U_pH(ldh) + U_mod_prot(ldh) + U_barr(ldh) This
     * method returns U_pH + U_mod_prot.
     *
     * @param temperature a double.
     * @return a double.
     */
    protected double getPhBias(double temperature) {
        double lambda = getLambda();
        double uph = LOG10 * Constants.R * temperature * (pKaModel - constPh) * (1.0 - lambda);
        double umod = -referenceEnergy * (lambda - lambdaIntercept) * (lambda - lambdaIntercept); //PMF Equation: A*(lambda-B)^2
        // TODO Find PMFs for monomers/trimers/pentapeptides.
        return uph + umod;
    }

    /**
     * getPhBiasDeriv.
     *
     * @param temperature a double.
     * @return a double.
     */
    protected double getPhBiasDeriv(double temperature) {
        double lambda = getLambda();
        double duphdl = LOG10 * Constants.R * temperature * (pKaModel - constPh) * -1.0;
        double dumoddl = -2.0 * referenceEnergy * (lambda - lambdaIntercept);
        return duphdl + dumoddl;
    }
}
