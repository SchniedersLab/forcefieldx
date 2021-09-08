//******************************************************************************
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
//******************************************************************************

import ffx.potential.bonded.Polymer
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import ffx.potential.GetProteinFeatures
import ffx.potential.bonded.AminoAcidUtils
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import java.nio.file.Files
import java.nio.file.Paths
import ffx.potential.bonded.Torsion



@Command(description = " Create a Feature Map for a given protein structure", name = "ffxc ffx.potential.FeatureMap")
class FeatureMap extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private String filename = null

    private List<Residue> residues
    private int[] polarity
    private int[] acidity
    private int[] secondaryStructure
    private double[] surfaceArea


    /**
     * ffx.potential.FeatureMap constructor.
     */
    FeatureMap() {
        this(new Binding())
    }

    /**
     * ffx.potential.FeatureMap constructor.
     * @param binding The Groovy Binding to use.
     */
    FeatureMap(Binding binding) {
        super(binding)
    }

    /**
     * ffx.potential.FeatureMap the script.
     */
    @Override
    FeatureMap run() {
        // Init the context and bind variables.
        if (!init()) {
            return null
        }

        // Load the MolecularAssembly.
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return null
        }

        residues = activeAssembly.getResidueList()
        GetProteinFeatures getProteinFeatures = new GetProteinFeatures()


        String fileDir = FilenameUtils.getFullPath(filename).replace("filename", "")
        String baseName = FilenameUtils.getBaseName(filename)
        String featureFileName = fileDir + baseName + ".csv"
        try {
            FileWriter fos = new FileWriter(featureFileName)
            PrintWriter dos = new PrintWriter(fos)
            dos.println("Residue\tPosition\tPolarity\tAcidity\tSecondary Structure\tPhi\tPsi\tSurface Area")
            for (int i = 0; i < residues.size(); i++) {
                getProteinFeatures.saveFeatures(residues.get(i))
                String[] features = getProteinFeatures.getFeatures()
                for (int j=0; j<features.length; j++){
                    dos.print(features[j] + "\t")
                }
                dos.println()
            }
            dos.close()
            fos.close()
        } catch (IOException e) {
            logger.info("Could Not Write Tab Delimited File")
        }


    }


}



