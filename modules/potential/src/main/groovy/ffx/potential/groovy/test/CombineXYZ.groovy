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
 * Merge two xyz files.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc Cluster [options] &lt;filename&gt;
 */
@Command(description = " Merge two xyz files.", name = "ffxc CombineXYZ")
class CombineXYZ extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "2", paramLabel = "files",
            description = 'Two xyz files.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    CombineXYZ run() {
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
        nArgs = (nArgs < 2) ? 2 : nArgs


        // Segment of code for MultiDynamics and OST.
        List<File> allFiles = arguments.stream().
                map { fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath()) }.
                collect(Collectors.toList())

        File newXYZ = new File("out.xyz")

        if (newXYZ.exists()){
            newXYZ.delete()
        }

        newXYZ.createNewFile()

        try {
            FileWriter fwnewXYZ = new FileWriter(newXYZ, true)
            BufferedWriter bwnewXYZ = new BufferedWriter(fwnewXYZ)

            File xyz1 = allFiles.get(0)
            File xyz2 = allFiles.get(1)

            FileReader frXYZ1 = new FileReader(xyz1)
            FileReader frXYZ2 = new FileReader(xyz2)
            BufferedReader brXYZ1 = new BufferedReader(frXYZ1)
            BufferedReader brXYZ2 = new BufferedReader(frXYZ2)
            String line;
            String[] tokens;
            int atomCount = 0;
            int atomCount2 = 0;
            //find number of atoms in file one
            if((line = brXYZ1.readLine()) != null){
                tokens = line.split("\\s+")
                atomCount = tokens[1].toInteger(); //I though should be zero.. but there is a "" character in 0...
            }else{
                logger.warning("Number of atoms not found at head of first file. Be sure you are using XYZ format.");
            }
            //Find number atoms in file two.
            if((line = brXYZ2.readLine()) != null){
                tokens = line.split("\\s+")
                atomCount2 = tokens[1].toInteger(); //I though should be zero.. but there is a "" character in 0...
            }else{
                logger.warning("Number of atoms not found at head of second file. Be sure you are using XYZ format.");
            }
            if(atomCount+atomCount2<10){
                line="      "
            }else if(atomCount+atomCount2<100){
                line = "     "
            }else{
                line = "    "
            }
            line+=atomCount+atomCount2
            bwnewXYZ.write(line)

            while((line = brXYZ1.readLine()) != null){
                bwnewXYZ.write(line + "\n")
            }
            //Skip the Crystal information in second xyz
            brXYZ2.readLine()
            while((line = brXYZ2.readLine()) != null){
                String line2 = "";
                int num;

                tokens = line.split("\\s+")
                logger.info("printing tokens:")
                for(String token in tokens){
                    logger.info(token);
                }
                logger.info("DONE printing tokens:")
                for(int i = 0; i<tokens.length;i++){
                    switch(i){
                        case 1:
                            //Change the atom number to account for new file
                            num = tokens[i].toInteger();
                            num+=atomCount;
                            tokens[i]=num.toString()

                            if (tokens[i].toInteger() < 10) {
                                line2 += "      "
                            } else if (tokens[i].toInteger() < 100) {
                                line2 += "     "
                            } else {
                                line2 += "    "
                            }
                            break;
                        case 2:
                            line2 += "   "
                            break;
                        case 6:
                            //Change the atom type to account for new file.
                            num = tokens[i].toInteger();
                            num+=atomCount;
                            tokens[i]=num.toString()

                            line2+="   "
                            break;
                        default:
                            if(i==0){
                                line2+="";
                            } else if(i>6) {
                                num = tokens[i].toInteger();
                                num += atomCount;
                                tokens[i] = num.toString()

                                if (tokens[i].toInteger() < 10) {
                                    line2 += "       "
                                } else if (tokens[i].toInteger() < 100) {
                                    line2 += "      "
                                } else {
                                    line2 += "     "
                                }
                            }else {
                                if (tokens[i].toDouble() < 0.0) {
                                    line2 += "   "
                                } else {
                                    line2 += "    "
                                }
                            }
                    }
                    line2+=tokens[i]
                }
                bwnewXYZ.write(line2 + "\n")
            }

            bwnewXYZ.close()
        } catch (IOException e) {
            logger.severe(e.toString());
        }

        return this
    }
}