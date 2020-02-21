package test

import ffx.potential.cli.PotentialScript

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

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The WriteSoluteLines script converts a TINKER *.prm file to Java properties.
 * <br>
 * Usage:
 * <br>
 * ffxc test.WriteSoluteLines &lt;filename&gt;
 */
@Command(description = "WriteSoluteLines writes out SOLUTE parameter lines using a map and input radii.", name = "ffxc WriteSoluteLines")
class WriteSoluteLines extends PotentialScript{

    /**
     * The final argument(s) should be two or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'SMARTs to atom type map and radii comparison file(s).')
    private List<String> filenames = null

    /**
     * -o or --option enables use of previously optmized radii to start script instead of vdw radii.
     */
    @CommandLine.Option(names = ['-o', '--option'], paramLabel = "0",
            description = "Specify whether a single input file will include both groups of radii (0), GK radii (1), or COSMO radii (2)")
    private int option = 0

    /**
     * -r or --actualRadii enables use of previously optmized radii to start script instead of vdw radii.
     */
    @CommandLine.Option(names = ['-r', '--actualRadii'], paramLabel = "false",
            description = "Specify if values entered are radii as opposed to diameters (usually diameters)")
    private boolean actualRadii = false

    ArrayList<String[]> mapList = new ArrayList<>()
    HashMap<String, String[]> gkRadiiList = new HashMap<>()
    HashMap<String, String[]> cosmoRadiiList = new HashMap<>()

    /**
     * Execute the WriteSoluteLines script.
     */
    @Override
    WriteSoluteLines run() {

        if (!init()) {
            return
        }

        // Read in the command line file.
        List<String> arguments = filenames
        String mapName = arguments.get(0)
        String gkRadiiFile
        String cosmoRadiiFile
        if(arguments.size() == 3){
            // Include both radii comp files, assume GK is first
            gkRadiiFile = arguments.get(1)
            cosmoRadiiFile = arguments.get(2)
        } else if(arguments.size() == 2){
            // Read in options flag and assign from there
            switch (option){
                case 1:
                    // GK radii given
                    gkRadiiFile = arguments.get(1)
                    break
                case 2:
                    // COSMO radii given
                    cosmoRadiiFile = arguments.get(1)
                    break
                default:
                    // Assume GK radii
                    gkRadiiFile = arguments.get(1)
            }
        }

        // Read map file
        // Map file should be header-less CSV file with three columns:
        // atomType,atomClass,SMARTS
        try {
            // input the (modified) file content to the StringBuffer "input"
            BufferedReader mapReader = new BufferedReader(new FileReader(mapName))
            String line

            while ((line = mapReader.readLine()) != null) {
                // Print lines to ArrayList
                // Index 0 of each array is atom type
                // Index 1 is atom class
                // Index 2 is SMARTS string description
                String[] splitLine = line.split("\\s+")
                mapList.add(splitLine)
            }
        } catch(IOException e){
            System.out.println("Could not read map file")
            e.printStackTrace()
        }

        // Read radii input file(s)
        // Currently set up to be radiiComp.csv files from optimizer output
        try{
            switch (option){
                case 0:
                    // Both radii groups
                    readGKRadiiFile(gkRadiiFile)
                    readCOSMORadiiFile(cosmoRadiiFile)
                    break
                case 1:
                    // GK radii only
                    readGKRadiiFile(gkRadiiFile)
                    break
                case 2:
                    // COSMO radii only
                    readCOSMORadiiFile(cosmoRadiiFile)
                    break
                default:
                    break
            }
        } catch(IOException e){
            System.out.println("Could not read radii file")
            e.printStackTrace()
        }

        // Write out SOLUTE lines for each atom type:
        // SOLUTE  atomTypeNum  SMARTS  GK-rad  COSMO-rad
        try {
            FileWriter soluteFile = new FileWriter("soluteLines.txt")

            for (int i = 0; i < mapList.size(); i++) {
                String[] mapEntry = mapList.get(i)
                String[] matchedGK = ["0.0000","0.0000"]
                String[] matchedCOSMO = ["0.0000","0.0000"]
                String gkRad = "0.0000"
                double gkRadNum = 0.0000
                String cosmoRad = "0.0000"
                double cosmoRadNum = 0.0000

                switch (option) {
                    case 0:
                        // Both radii groups
                        matchedGK = gkRadiiList.get(mapEntry[2])
                        matchedCOSMO = cosmoRadiiList.get(mapEntry[2])
                        if(matchedGK != null) {
                            gkRad = matchedGK[1]
                            gkRadNum = Double.parseDouble(gkRad)
                            gkRad = String.format("%.4f",gkRadNum)
                        }else {gkRad = "0.0000"}
                        if(matchedCOSMO != null) {
                            cosmoRad = matchedCOSMO[1]
                            cosmoRadNum = Double.parseDouble(cosmoRad)
                            cosmoRad = String.format("%.4f",cosmoRadNum)
                        }else {cosmoRad = "0.0000"}
                        //soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+mapEntry[2]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        break
                    case 1:
                        // GK only
                        matchedGK = gkRadiiList.get(mapEntry[2])
                        gkRad = matchedGK[1]
                        cosmoRad = matchedGK[0]
                        //soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+mapEntry[2]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        break
                    case 2:
                        // COSMO only
                        matchedCOSMO = cosmoRadiiList.get(mapEntry[2])
                        gkRad = matchedCOSMO[0]
                        cosmoRad = matchedCOSMO[1]
                        //soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+mapEntry[2]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        soluteFile.write("SOLUTE\t"+mapEntry[0]+"\t"+gkRad+"\t"+cosmoRad+"\n")
                        break
                    default:
                        break
                }
            }
            soluteFile.close()
        } catch(IOException e){
            System.out.println("Exception writing solute lines file")
            e.printStackTrace()
        }

        /*CompositeConfiguration properties = Keyword.loadProperties(null)
        properties.setProperty("parameters", mapName)
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties)

        ForceField forceField = forceFieldFilter.parse()

        int prms = arguments.size()
        for (int i = 1; i < prms; i++) {
            mapName = arguments.get(i)
            properties = Keyword.loadProperties(null)
            properties.setProperty("parameters", mapName)
            forceFieldFilter = new ForceFieldFilter(properties)
            ForceField forceField2 = forceFieldFilter.parse()
            forceField.append(forceField2)
        }

        if (forceField != null) {
            forceField.print()
        }*/

        return this
    }

    private void readGKRadiiFile(String gkRadiiFile){
        BufferedReader radiiReader = new BufferedReader(new FileReader(gkRadiiFile))
        String gkLine
        while((gkLine = radiiReader.readLine()) != null){
            String[] splitLine = gkLine.split("\\s+")
            String smarts = splitLine[0]
            String[] radii = new String[2]
            radii[0] = splitLine[2]
            radii[1] = splitLine[3]
            gkRadiiList.put(smarts,radii)
        }
    }

    private void readCOSMORadiiFile(String cosmoRadiiFile){
        BufferedReader radiiReader = new BufferedReader(new FileReader(cosmoRadiiFile))
        String cosmoLine
        while((cosmoLine = radiiReader.readLine()) != null){
            String[] splitLine = cosmoLine.split("\\s+")
            String smarts = splitLine[0]
            String[] radii = new String[2]
            radii[0] = splitLine[2]
            radii[1] = splitLine[3]
            cosmoRadiiList.put(smarts,radii)
        }
    }
}
