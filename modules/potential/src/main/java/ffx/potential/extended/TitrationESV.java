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

import static java.lang.String.format;

import ffx.potential.bonded.MultiResidue;
import ffx.potential.extended.TitrationUtils.Titration;

/**
 * An extended system variable that allows continuous fractional protonation of an amino acid.
 * All atomic charges and bonded terms scale linearly between prot and deprot states.
 * 
 * Possible expansions:
 *  (1) Add back the ability to interact with OSRW lambda and thereby combine with protein design (QuadTop).
 *  (2) Allow triple-state systems such as histidine with 0.5 protons per nitrogen or tautomeric ASP/GLU.
 *  (3) Assess bytecode-output implementation via eg. ASM.
 *
 * @author slucore
 */
public final class TitrationESV extends ExtendedVariable {

    private final MultiResidue titrating;           // from TitrationUtils.titrationFactory()
    private final double referenceEnergy;           // deprotonation free energy of a model tripeptide
    private final double constPh;                   // Simulation pH.
    private final double pKaModel;                  // Reference pKa value.
    
    public TitrationESV(ExtendedSystem esvSystem, MultiResidue multiRes) {
        super(esvSystem, multiRes, 1.0);
        this.constPh = esvSystem.getConstantPh();
        this.titrating = multiRes;
        Titration titration = Titration.lookup(titrating.getActive());
        this.referenceEnergy = titration.refEnergy;
        this.pKaModel = titration.pKa;
    }

    @Override
    protected double getTotalBias(double temperature, boolean print) {
        double eDiscr = getDiscrBias();
        double ePh = getPhBias(temperature);
        return (eDiscr + ePh);
    }

    @Override
    protected double getTotalBiasDeriv(double temperature, boolean print) {
        double dDiscr = getDiscrBiasDeriv();
        double dPh = getPhBiasDeriv(temperature);
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
        double lambda = getLambda();
        double uph = ExtConstants.log10*ExtConstants.Boltzmann*temperature*(pKaModel - constPh)*lambda;
        double umod = referenceEnergy * lambda;     // TODO Find PMFs for monomers/trimers/pentapeptides.
        return uph + umod;
    }

    protected double getPhBiasDeriv(double temperature) {
        double duphdl = ExtConstants.log10*ExtConstants.Boltzmann*temperature*(pKaModel - constPh);
        double dumoddl = referenceEnergy;
        return duphdl + dumoddl;
    }

    @Override
    public String getName() {
        return format("Titr%d", esvIndex);
    }

}
