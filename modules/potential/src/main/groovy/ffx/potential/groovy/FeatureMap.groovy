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

import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.GetProteinFeatures
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.io.FilenameUtils.getFullPath

@Command(description = " Create a Feature Map for a given protein structure", name = "FeatureMap")
class FeatureMap extends PotentialScript {

    @Option(names = ["-d", "--delimiter"], paramLabel = ",",
            description = "Delimiter of input variant list file")
    private String delimiter = ","

    @Option(names = ["--iP", "--includePolarity"], paramLabel = "false",
            description = "Include polarity change in feature map.")
    private boolean includePolarity = false

    @Option(names = ["--iA", "--includeAcidity"], paramLabel = "false",
            description = "Include acidity change in feature map.")
    private boolean includeAcidity = false

    @Option(names = ["--iAng", "--includeAngles"], paramLabel = "false",
            description = "Include ramachandran angles in feature map.")
    private boolean includeAngles = false

    @Option(names = ["--iS", "--includeStructure"], paramLabel = "false",
            description = "Include secondary structure annotations in feature map.")
    private boolean includeStructure = false

    @Option(names = ["--mI", "--multiple isomers"], paramLabel = "false",
            description = "Set this flag if the variant list contains variants from multiple isomers. Isomer should be in name of pdb file")
    private boolean multipleIsomers = false
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

        // Handles when variant files will have multiple isoforms and will need to isoform specific variants when
        // writing the file csv file
        String[] geneSplit
        String fileIsomer
        if(multipleIsomers){
            String baseName = getBaseName(filenames[0])
            if(baseName.contains("ENS")){
                geneSplit = baseName.replace(".pdb", "").split("_")
                fileIsomer = geneSplit
            } else {
                geneSplit = baseName.replace(".pdb", "").split("_")
                fileIsomer = geneSplit[1] + '_' + geneSplit[2]
            }
        }

        List<String[]> featureList = new ArrayList<>()
        //Store all features for each residue in an array list called Features
        for (int i = 0; i < residues.size(); i++) {
            double residueSurfaceArea =
                    forceFieldEnergy.getGK().getSurfaceAreaRegion().getResidueSurfaceArea(residues.get(i))
            featureList.add(getProteinFeatures.saveFeatures(residues.get(i), residueSurfaceArea, includeAngles, includeStructure))
        }

        BufferedReader txtReader = null

        // Parse the ddGun output file for extracting free energy values for the map
        List<String> ddgunLines = new ArrayList<>()
        try {
            File freeEnergyFile = new File(filenames[2])
            txtReader = new BufferedReader(new FileReader(freeEnergyFile));
            String line = txtReader.readLine();
            while (line != null) {
                if (line.contains('.pdb')) {
                    ddgunLines.add(line)
                }
                // read next line
                line = txtReader.readLine();
            }
            txtReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Get ddGun features
        List<String> npChanges = getProteinFeatures.ddgunToNPChange(ddgunLines)
        List<Double[]> ddGun = getProteinFeatures.getDDGunValues(ddgunLines)

        // Get Polarity and Acidity Features if flags are set
        List<String[]> polarityAndAcidityChange = new ArrayList<>()
        if(includePolarity || includeAcidity){
            polarityAndAcidityChange = getProteinFeatures.getPolarityAndAcidityChange(npChanges,
                    includePolarity, includeAcidity)
        }


        BufferedReader br = null;
        BufferedWriter bw = null;

        // Write a new csv file with all the original data and the features determined through this script
        try {
            File inputCSVFile = new File(filenames[1])
            String inputCSVPath = inputCSVFile.getAbsolutePath().replace( "/" +filenames[1], '')
            String newCSVFileName = "update_" + filenames[1]
            File updatedFile = new File(inputCSVPath, newCSVFileName)
            br = new BufferedReader(new InputStreamReader(new FileInputStream(inputCSVFile)))
            if (newCSVFileName.length() == 0) {
                bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(updatedFile)))
            } else {
                FileWriter fw = new FileWriter(updatedFile, true)
                bw = new BufferedWriter(fw)
            }

            String line = null
            int i = 0
            int npIndex
            int isoformIndex
            for (line = br.readLine(); line != null; line = br.readLine(), i++) {
                StringBuilder newCSVLine = new StringBuilder()
                if (i == 0 || i == 1) {
                    if (updatedFile.length() == 0 && i == 1) {
                        newCSVLine.append(line + delimiter +'Surface Area'+ delimiter + 'Normalized SA'+ delimiter +
                                'Confidence Score'+ delimiter + 'ddG' + delimiter + '|ddG|')
                        if(includeAcidity){
                            newCSVLine.append(delimiter + 'Acidity Change')
                        }
                        if(includePolarity){
                            newCSVLine.append(delimiter + 'Polarity Change')
                        }
                        if(includeAngles){
                            newCSVLine.append(delimiter + 'Phi' + delimiter + 'Psi' + delimiter + 'Omega')
                        }
                        if(includeStructure){
                            newCSVLine.append(delimiter + 'Secondary Structure Annotation')
                        }
                        bw.write(newCSVLine.toString())
                    } else if (i == 0 && updatedFile.length() == 0) {
                        bw.write(line + '\n')
                    }
                } else {
                    String[] splits = line.split(delimiter)
                    if(i == 1 || i==2){
                        for(int j=0; j < splits.length; j++){
                            if(splits[j].contains("p.")){
                                npIndex = j
                            }
                            if(multipleIsomers){
                                if(splits[j].contains("NP") || splits[j].contains("XP") || splits[j].contains("ENS")){
                                    isoformIndex = j
                                }
                            }
                        }
                    }
                    int length = splits.length
                    int position
                    String[] ddG = ""
                    String[] pA = ""
                    String npChange = splits[npIndex]

                    String[] feat = new String[featureList.get(0).length]
                    if (splits[npIndex].contains('del')) {
                        ddG = ["null", "null"]
                        Arrays.fill(feat,null)
                    } else {
                        String proteinChange = npChange.substring(npChange.indexOf('p'))
                        String splitstring = "p\\."
                        String sub = proteinChange.split(splitstring)[1]
                        position = sub.replace(sub.substring(0,3),'').replace(sub.substring(sub.length()-3, sub.length()), '').toInteger()
                        if (position <= residues.size()) {
                            if (npChanges.indexOf(proteinChange) != -1) {
                                ddG = ddGun.get(npChanges.indexOf(proteinChange))
                                if(includeAcidity || includePolarity){
                                    pA = polarityAndAcidityChange.get(npChanges.indexOf(proteinChange))
                                }
                            } else {
                                ddG = ["null", "null"]
                                pA = ["null", "null"]
                            }
                            feat = featureList.get(position - 1)
                        }
                    }
                    String isomer
                    if(multipleIsomers){
                        if(splits[isoformIndex].contains(":p.")){
                            isomer = splits[isoformIndex].split(":")[0]
                        } else {
                            isomer = splits[isoformIndex]
                        }

                    }

                    if (length == splits.length) {
                        if(multipleIsomers){
                            if(isomer != fileIsomer){
                                continue
                            }
                        }
                        newCSVLine.append(line + delimiter + feat[0] + delimiter + feat[1] + delimiter + feat[2] + delimiter
                                + String.valueOf(ddG[0]) + delimiter + String.valueOf(ddG[1]))
                        if(includeAcidity){
                            newCSVLine.append(delimiter + pA[0])
                        }
                        if(includePolarity){
                            newCSVLine.append(delimiter + pA[1])
                        }
                        if(includeAngles){
                            newCSVLine.append(delimiter + feat[3] + delimiter + feat[4] + delimiter + feat[5])
                            if(includeStructure){
                                newCSVLine.append(delimiter + feat[6])
                            }
                        } else if(includeStructure){
                            newCSVLine.append(delimiter + feat[3])
                        }
                        bw.newLine()
                        bw.write(newCSVLine.toString())
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



