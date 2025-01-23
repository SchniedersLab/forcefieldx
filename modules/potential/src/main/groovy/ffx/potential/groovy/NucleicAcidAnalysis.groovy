//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
//******************************************************************************
package ffx.potential.groovy

import ffx.numerics.math.DoubleMath
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Atom
import ffx.potential.bonded.NucleicAcidUtils
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.GetProteinFeatures
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.math3.util.FastMath.toDegrees

@Command(description = " Create a Feature Map for a given protein structure", name = "FeatureMap")
class NucleicAcidAnalysis extends PotentialScript {

//    @Option(names = ["--iP", "--includePolarity"], paramLabel = "false",
//            description = "Include polarity change in feature map.")
//    private boolean includePolarity = false

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format, variant list file, and a free energy file')
    private List<String> filenames = null

    private List<Residue> residues

    /**
     * ffx.potential.FeatureMap constructor.
     */
    NucleicAcidAnalysis() {
        this(new Binding())
    }

    /**
     * ffx.potential.FeatureMap constructor.
     * @param binding The Groovy Binding to use.
     */
    NucleicAcidAnalysis(Binding binding) {
        super(binding)
    }

    /**
     * ffx.potential.FeatureMap the script.
     */
    @Override
    NucleicAcidAnalysis run() {
        // Init the context and bind variables.
        if (!init()) {
            return null
        }

        // Load the MolecularAssembly.
        activeAssembly = getActiveAssembly(filenames[0])
        if (activeAssembly == null) {
            logger.info(helpString())
            return null
        }

        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)
        forceFieldEnergy.energy(x)

        residues = activeAssembly.getResidueList()

        for (Residue residue : residues) {
            NucleicAcidUtils.NucleicAcid3 nucleicAcid = residue.getNucleicAcid3()
            Atom o4p = residue.getAtomByName("O4'", true)
            Atom c1p = residue.getAtomByName("C1'", true)
            Atom baseN = null
            Atom baseC = null
            if (nucleicAcid == NucleicAcidUtils.NucleicAcid3.DAD || nucleicAcid == NucleicAcidUtils.NucleicAcid3.DGU) {
                baseN = residue.getAtomByName("N9", true)
                baseC = residue.getAtomByName("C4", true)
            } else if (nucleicAcid == NucleicAcidUtils.NucleicAcid3.DCY || nucleicAcid == NucleicAcidUtils.NucleicAcid3.DTY) {
                baseN = residue.getAtomByName("N1", true)
                baseC = residue.getAtomByName("C2", true)
            }

            if (o4p != null && c1p != null && baseN != null && baseC != null) {
                double angle = DoubleMath.dihedralAngle(o4p.getXYZ(null), c1p.getXYZ(null), baseN.getXYZ(null), baseC.getXYZ(null))
                double glyco = toDegrees(angle)
                logger.info("RES : " + residue.name + " = " + glyco)
            }
        }

        return this
    }
}



