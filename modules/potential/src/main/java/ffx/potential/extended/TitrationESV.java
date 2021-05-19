// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.extended.TitrationUtils.Titration;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

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
        ResidueEnumerations.AminoAcid3 currentAA3 = ResidueEnumerations.AminoAcid3.valueOf(currentRes.getName());
        if(currentAA3 == ResidueEnumerations.AminoAcid3.LYS || currentAA3 == ResidueEnumerations.AminoAcid3.LYD){
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
     * {@inheritDoc}
     */
    @Override
    public String getName() {
        return format("Titr%d", esvIndex);
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
