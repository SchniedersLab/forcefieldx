/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.extended;

import java.util.logging.Logger;

import static java.lang.String.format;

import ffx.potential.bonded.MultiResidue;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.extended.ExtendedSystem.EsvConfiguration;
import ffx.potential.extended.TitrationUtils.Titr;

/**
 * An extended system variable that allows continuous fractional protonation of an amino acid.
 * All atomic charges and bonded terms scale linearly between prot and deprot states.
 *
 * TODO List:
 * TODO Assess need to finalize MultiResidues before adding and using them in TitrationESVs.
 *
 * Possible expansions:
 *  (1) Add back the ability to interact with OSRW lambda and thereby combine with protein design (QuadTop).
 *  (2) Allow triple-state systems such as histidine with 0.5 protons per nitrogen or tautomeric ASP/GLU.
 *  (3) Assess bytecode-output implementation via eg. ASM.
 *
 * @author slucore
 */
public final class TitrationESV extends ExtendedVariable {

    // System handles
    private static final Logger logger = Logger.getLogger(TitrationESV.class.getName());
    private final MultiResidue titrating;           // from TitrationUtils.titrationFactory()
    private final double referenceEnergy;           // deprotonation free energy of a model tripeptide

    private final double constPh;                   // Simulation pH.
    private final double pKaModel;                  // Reference pKa value.

    /**
     * @{inherit_doc}
     */
    public TitrationESV(EsvConfiguration esvConfig, MultiResidue multiRes, double constPh, double biasMag) {
        super(esvConfig, multiRes, biasMag, 1.0);
        this.constPh = constPh;
        this.titrating = multiRes;
        Titr titr = TitrationUtils.titrationLookup(titrating.getActive());
        this.referenceEnergy = titr.refEnergy;
        this.pKaModel = titr.pKa;
    }

    public TitrationESV(EsvConfiguration esvConfig, MultiResidue titrating, double constPh) {
        this(esvConfig, titrating, constPh, 1.0);
    }

    public MultiResidue getMultiRes() {
        return titrating;
    }

    @Override
    protected double getTotalBias(double temperature, boolean print) {
        double eDiscr = getDiscrBias();
        double ePh = getPhBias(temperature);
        if (print) {
            SB.logfn(" eDiscr %d: %g", index, eDiscr);
            SB.logfn(" epH    %d: %g", index, ePh);
        }
        return (eDiscr + ePh);
    }

    @Override
    protected double getTotalBiasDeriv(double temperature, boolean print) {
        double dDiscr = getDiscrBiasDeriv();
        double dPh = getPhBiasDeriv(temperature);
        if (print) {
            SB.logfn("  Discr  %d: %g", index, dDiscr);
            SB.logfn("  pH     %d: %g", index, dPh);
        }
        return (dDiscr + dPh);
    }

    /**
     * Eqs 5,6 from Wallace+Shen 2011 "Continuous constant pH M.D. in explicit..."
     * U_pH(ldh) = log(10)*kb*T*(pKa_model - pH)*ldh
     * U_mod(ldh) = potential of mean force for protonation (or -deprot) of model compound
     * U_star = sum(ldh) { U_pH(ldh) + U_mod_prot(ldh) + U_barr(ldh)
     * This method returns U_pH + U_mod_prot.
     */
    protected double getPhBias(double temperature) {
        double lambda = getLambdaSwitch();
        double uph = ExtConstants.log10*ExtConstants.Boltzmann*temperature*(pKaModel - constPh)*lambda;
//        buglog(" U(pH): 2.303kT*(pKa-pH)*L = %.4g * (%.2f - %.2f) * %.2f = %.4g",
//                ExtConstants.log10*ExtConstants.Boltzmann*temperature,
//                pKaModel, constPh, lambda, uph);
        double umod = referenceEnergy * lambda;     // TODO PRIO find PMFs for monomers/trimers/pentapeptides
//        return uph + umod;
        return umod;
    }

    protected double getPhBiasDeriv(double temperature) {
        double chain = getSwitchDeriv();
        double duphdl = ExtConstants.log10*ExtConstants.Boltzmann*temperature*(pKaModel - constPh);
//        buglog(" dU(pH)dL: 2.303kT*(pKa-pH) = %.4g * (%.2f - %.2f) = %.4g",
//                ExtConstants.log10*ExtConstants.Boltzmann*temperature,
//                pKaModel, constPh, duphdl);
        double dumoddl = referenceEnergy * chain;
//        return duphdl + dumoddl;
        return dumoddl;
    }

    @Override
    public String getName() {
        return format("ESV-Titr_%d", index);
    }

}
