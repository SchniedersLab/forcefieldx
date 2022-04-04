//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy

import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.GetProteinFeatures
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.io.FilenameUtils.getFullPath

@Command(description = " Create a Feature Map for a given protein structure", name = "ffxc FeatureMap")
class FeatureMap extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private List<String> filenames = null

    private List<Residue> residues

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
        System.setProperty("gkterm", "true")
        System.setProperty("cavmodel", "CAV")
        System.setProperty("surface-tension", "1.0")
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
        GetProteinFeatures getProteinFeatures = new GetProteinFeatures()

        // Use the current base directory, or update if necessary based on the given filename.
        String dirString = getBaseDirString(filenames[0])

        String baseName = getBaseName(filenames[0])
        String csvPath = getBaseDirString(filenames[0]).replace(baseName + '/', '')
        String[] geneSplit =  baseName.split('_')
        String csvFileName = geneSplit[0] + ".csv"

        logger.info(csvFileName)

        List<String[]> featureList = new ArrayList<>()
        //Store all features for each residue in an array list called Features
        for (int i = 0; i < residues.size(); i++) {
            double residueSurfaceArea =
                    forceFieldEnergy.getGK().getSurfaceAreaRegion().getResidueSurfaceArea(residues.get(i))
            featureList.add(getProteinFeatures.saveFeatures(residues.get(i), residueSurfaceArea))
        }

        BufferedReader txtReader = null;
        List<String> ddgunLines = new ArrayList<>()
        try {
            File txtfile = new File(dirString, "output.txt")
            txtReader = new BufferedReader(new FileReader(txtfile));
            String line = txtReader.readLine();
            while (line != null) {
                if(line.contains('.pdb')){
                    ddgunLines.add(line)
                }
                // read next line
                line = txtReader.readLine();
            }
            txtReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<String> npChanges = getProteinFeatures.ddgunToNPChange(ddgunLines)
        List<Double[]> ddGun = getProteinFeatures.getDDGunValues(ddgunLines)

        BufferedReader br = null;
        BufferedWriter bw = null;

        final String lineSep = System.getProperty("line.separator");
        logger.info(lineSep)

        try {
            File file = new File(csvPath, csvFileName)
            File updatedFile = new File(csvPath, "update_" + csvFileName)

            br = new BufferedReader(new InputStreamReader(new FileInputStream(file)))
            if (csvFileName.length() == 0){
                bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(updatedFile)))
            } else {
                FileWriter fw = new FileWriter(updatedFile,true)
                bw = new BufferedWriter(fw)
            }

            String line = null;
            int i = 0;
            for (line = br.readLine(); line != null; line = br.readLine(), i++) {
                if (i == 0 || i == 1) {

                    if (updatedFile.length() == 0 && i == 1) {
                        bw.write(line + ',\"Surface Area\",\"Normalized SA\",\"Confidence Score\",\"ddG\",\"|ddG|\"')
                    } else if (i == 0 && updatedFile.length() == 0){
                        bw.write(line+ '\n')
                    }
                }  else {
                    String[] splits = line.split('\",\"')
                    int length = splits.length
                    int position = splits[8].toInteger()
                    String proteinChange = splits[2]
                    String[] ddG = ""
                    if (npChanges.indexOf(proteinChange) != -1){
                        ddG = ddGun.get(npChanges.indexOf(proteinChange))
                    } else {
                        ddG = ["null","null"]
                    }
                    String isomer = proteinChange.split(':p.')[0]

                    if (length == 14 && isomer == geneSplit[1]+ '_' + geneSplit[2]) {
                        String[] feat = featureList.get(position - 1)
                        bw.newLine()
                        bw.write(line + '\",\"' + feat[0] + '\",\"' + feat[1] + '\",\"' + feat[2] + '\",\"' + String.valueOf(ddG[0]) +'\"' + String.valueOf(ddG[1])+ '\"')
                    }
                }

            }
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            if (br != null)
                br.close();
            if (bw != null)
                bw.close();
        }

        logger.info(format("\n Total SurfacAreaRegion Solvent Accessible Surface Area: %1.6f",
                forceFieldEnergy.getGK().getSurfaceAreaRegion().getEnergy()))
        logger.info(format("\n Total Calculated Solvent Accessible Surface Area: %1.6f",
                getProteinFeatures.getTotalSurfaceArea()))


    }


}



