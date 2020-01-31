//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

package test

import ffx.potential.Utilities
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import org.apache.commons.math3.ml.clustering.CentroidCluster
import org.apache.commons.math3.ml.clustering.Clusterable
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer

import java.nio.file.Path
import java.util.logging.Level
import java.util.stream.Collectors;

/**
 * Wrapper command for the PACCOM algorithm
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc Cluster [options] &lt;filename&gt;
 */
@Command(description = " Compute RMSD for crystals using PACCOM.", name = "ffxc CallPAC")
class CallPAC extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'PDB, XYZ, ARC files to be compared.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    CallPAC run() {
        if (!init()) {
            return this
        }

        if (filenames == null || filenames.isEmpty()) {
            logger.info(helpString())
            return this
        }
        List<String> arguments = filenames
        // Check nArgs should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1
        nArgs = (nArgs < 1) ? 1 : nArgs


        // Segment of code for MultiDynamics and OST.
        List<File> allFiles = arguments.stream().
                map { fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath()) }.
                collect(Collectors.toList())
        boolean xyz = false;
        List<File> xyzFiles = new ArrayList<File>();
        boolean arc = false;
        List<File> arcFiles = new ArrayList<File>();
        boolean pdb = false;
        List<File> pdbFiles = new ArrayList<File>();
        for (File f in allFiles) {
            String ext = FilenameUtils.getExtension(f.getName());
            //TODO compatible with other file types?
            if (ext.toUpperCase().startsWith("XYZ")) {
                xyzFiles.add(f)
                xyz = true;
            }
            if (ext.toUpperCase().startsWith("ARC")){
                arcFiles.add(f)
                arc = true;
            }
            if (ext.toUpperCase().startsWith("PDB")){
                pdbFiles.add(f)
                pdb = true;
            }
        }
        if (xyzFiles.isEmpty() && arcFiles.isEmpty() && pdbFiles.isEmpty()) {
            logger.severe("No inputted structure files.")
        }

        String line;
        if (arc){ //if arc break down into xyz files
            File newXYZ = new File("temp")
            FileWriter fwNewXYZ = new FileWriter(newXYZ);
            BufferedWriter bwNewXYZ = new BufferedWriter(fwNewXYZ);

            for(File f in arcFiles){
                //Need to separate out arc files into individual xyz for PACCOM.

                int fileCount = 0;
                FileReader frArc = new FileReader(f);
                BufferedReader brArc = new BufferedReader(frArc);
                while((line = brArc.readLine()) != null){
                    if(line.contains("_opt.xyz_") || line.contains(".arc")){
                        bwNewXYZ.close();
                        if(fileCount>0){
                            xyzFiles.add(newXYZ)
                        }
                        newXYZ = new File(String.format(FilenameUtils.getBaseName(f.getName())+"_PAC.xyz_%d", ++fileCount))
                        fwNewXYZ = new FileWriter(newXYZ)
                        bwNewXYZ = new BufferedWriter(fwNewXYZ)
                        bwNewXYZ.write(line + "\n")
                    }else{
                        bwNewXYZ.write(line + "\n")
                    }
                }
            }
        }

        File data33 = new File("data33.in")

        File structs1 = new File("mol500_1.dat")
        File props1 = new File("mol500_1p.dat")
        File structs2 = new File("mol500_2.dat")
        File props2 = new File("mol500_2p.dat")

        File output = new File("output.dat")

        if (data33.exists()){
            data33.delete()
        }
        if (structs1.exists()) {
            structs1.delete()
        }
        if (props1.exists()) {
            props1.delete()
        }
        if (structs2.exists()) {
            structs2.delete()
        }
        if (props2.exists()) {
            props2.delete()
        }
        if (output.exists()){
            output.delete()
        }
        data33.createNewFile()
        structs1.createNewFile()
        props1.createNewFile()
        structs2.createNewFile()
        props2.createNewFile()
        output.createNewFile()

        try {
            FileWriter fwData33 = new FileWriter(data33, true)
            FileWriter fwStructs1 = new FileWriter(structs1, true)
            FileWriter fwProps1 = new FileWriter(props1, true)
            FileWriter fwStructs2 = new FileWriter(structs2, true)
            FileWriter fwProps2 = new FileWriter(props2, true)

            BufferedWriter bwData33 = new BufferedWriter(fwData33)
            BufferedWriter bwStructs1 = new BufferedWriter(fwStructs1)
            BufferedWriter bwProps1 = new BufferedWriter(fwProps1)
            BufferedWriter bwStructs2 = new BufferedWriter(fwStructs2)
            BufferedWriter bwProps2 = new BufferedWriter(fwProps2)

            //TODO Make inputs to data33.in flags.
            if(pdb) {
                bwData33.write("30" + "\n" + "8" + " " + "6" + " " + "8" + "\n" + "32" + "\n" + "32" + "\n" + "2" + "\n" + "6" + "\n" + "7" + "\n" + "PDB" + "\n" + "FFX")
            }else{
                bwData33.write("30" + "\n" + "8" + " " + "6" + " " + "8" + "\n" + "32" + "\n" + "32" + "\n" + "2" + "\n" + "6" + "\n" + "7" + "\n" + "FFX" + "\n" + "FFX")

            }

            for (File structFile in xyzFiles) {
                if(pdb){
                    for(File structPDB in pdbFiles){
                        bwStructs1.write(structPDB.getName() + "\n")
                    }
                }else {
                    bwStructs1.write(structFile.getName() + "\n")
                }
                bwStructs2.write(structFile.getName() + "\n")
                logger.info(structFile.getName())
                //TODO Better determination of properties/key file
                String baseName
                if(arc) {
                    baseName = structFile.getName().replaceAll("_PAC.xyz_.*", "")
                }
                if(xyz){
                    baseName = structFile.getName().replaceAll("_opt.xyz_.*", "")
                }
                baseName = baseName.replaceAll(".xyz.*", "")
                File tempProp1 = new File(baseName + ".properties")
                File tempProp2 = new File(baseName + ".key")
                if (tempProp1.exists()) {
                    if(!pdb) {
                        bwProps1.write(baseName + ".properties" + "\n")
                    }
                    bwProps2.write(baseName + ".properties" + "\n")
                } else if (tempProp2.exists()) {
                    if(!pdb) {
                        bwProps1.write(baseName + ".key" + "\n")
                    }
                    bwProps2.write(baseName + ".key" + "\n")
                } else {
                    logger.severe("Property/Key file not found...")
                }
            }
            bwData33.close()
            bwStructs1.close()
            bwProps1.close()
            bwStructs2.close()
            bwProps2.close()
        } catch (IOException e) {
            logger.severe(e.toString());
        }

        if(pdb){
            if (props1.exists()) {
                props1.delete()
            }
        }

        ProcessBuilder processBuilder = new ProcessBuilder()
//        Map<String, String> envMap = processBuilder.environment();
//        Set<String> keys = envMap.keySet();
//        for(String key:keys){
//            logger.info(String.format((key+" ==> "+envMap.get(key));
//        }
        // Run this on Windows, cmd, /c = terminate after this run
        processBuilder.command("/bin/bash", "-c", "/Users/anessler/Research/MTPC/PACCOM/src_bin_analysis_forFFX_33_Aaron_org/comp_rmsd_33_ffx_L");

        try {
            logger.info(String.format("Running PACCOM"));
            FileWriter fwOutput = new FileWriter(output, true)
            BufferedWriter bwOutput = new BufferedWriter(fwOutput)

            Process process = processBuilder.start();
            logger.info("PACCOM has Started.")
            BufferedReader reader =
                    new BufferedReader(new InputStreamReader(process.getInputStream()));
            ArrayList<String> distances = new ArrayList<>();
            logger.info("Reading input")
            String array = ""
            if(!pdb) {
                while ((line = reader.readLine()) != null) {
                    int rmsdEnd = line.indexOf("   "); // three spaces go behind rmsd
                    distances.add(line.substring(2, rmsdEnd)); //two blank spaces before rmsd...
                }
                for (int i = 0; i < xyzFiles.size(); i++) {
                    File cluster = new File(String.format("Cluster_%d", i+1))
                    FileWriter fwCluster = new FileWriter(cluster);
                    BufferedWriter bwCluster = new BufferedWriter(fwCluster);
                    String tempLine = ""
                    for (int j = 0; j < Math.pow(xyzFiles.size(), 2); j += xyzFiles.size()) {
                        tempLine += distances.get(i + j)
                        tempLine += "\t"

                        //Threshold?
                        if(distances.get(i + j).toDouble()<1.0){
                            bwCluster.write((i+1)+"\t"+((j/xyzFiles.size())+1)+"\n")
                        }
                        //?Threshold
                    }
                    array += tempLine
                    array += "\n"
                    logger.info(tempLine)
                    bwCluster.close();
                }

            }else{
                while ((line = reader.readLine()) != null){
                    logger.info(String.format(line))
                    array += line
                    array += "\n"
                }
            }
            bwOutput.write(array)
            int exitCode = process.waitFor();
//            logger.info(String.format("\nExited with error code : " + exitCode));

            bwOutput.close()
            if(exitCode == 0){
                logger.info(String.format("Results written to output.dat"))
            }

        } catch (IOException e) {
            logger.info(String.format("Be certain that PACCOM is setup using the following command: comp_rmsd_33_ffx_L"))
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        if (data33.exists()){
            data33.delete()
        }
        if (structs1.exists()) {
            structs1.delete()
        }
        if (props1.exists()) {
            props1.delete()
        }
        if (structs2.exists()) {
            structs2.delete()
        }
        if (props2.exists()) {
            props2.delete()
        }

        return this
    }
}
