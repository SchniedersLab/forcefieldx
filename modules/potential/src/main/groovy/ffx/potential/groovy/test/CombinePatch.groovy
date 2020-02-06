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
 * Merge two patch files.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc Cluster [options] &lt;filename&gt;
 */
@Command(description = " Merge two patch files.", name = "ffxc CombinePatch")
class CombinePatch extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "2", paramLabel = "files",
            description = 'Two patch files.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    CombinePatch run() {
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

        File newPatch = new File("out.patch")

        if (newPatch.exists()){
            newPatch.delete()
        }

        newPatch.createNewFile()

        try {
            FileWriter fwNewPatch = new FileWriter(newPatch, true)
            BufferedWriter bwNewPatch = new BufferedWriter(fwNewPatch)

            File patch1 = allFiles.get(0)
            File patch2 = allFiles.get(1)

            FileReader frPatch1 = new FileReader(patch1)
            FileReader frPatch2 = new FileReader(patch2)
            BufferedReader brPatch1 = new BufferedReader(frPatch1)
            BufferedReader brPatch2 = new BufferedReader(frPatch2)
            String line;
            int atomCount = 0;
            while((line = brPatch1.readLine()) != null){
                if(line.contains("atom")){
                    bwNewPatch.write(line + "\n")
                    atomCount++;
                }else{
                    bwNewPatch.write(line + "\n")
                }
            }
            while((line = brPatch2.readLine()) != null) {
                String line2 = "";
                int num

                String[] tokens = line.split("\\s+")
                String tokensOut = ""
                for(String token in tokens){
                    tokensOut +=token;
                    tokensOut +=" ";
                }
                logger.info(tokensOut);
                switch(tokens[0]){
                    case "atom":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()
                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "vdw":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()
                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "bond":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()
                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "angle":
                    case "anglep":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        num = tokens[3].toInteger()
                        num+=atomCount;
                        tokens[3]=num.toString()
                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "strbnd":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        num = tokens[3].toInteger()
                        num+=atomCount;
                        tokens[3]=num.toString()

                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break
                    case "opbend":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "torsion":
                    case " torsion":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        num = tokens[3].toInteger()
                        num+=atomCount;
                        tokens[3]=num.toString()

                        num = tokens[4].toInteger()
                        num+=atomCount;
                        tokens[4]=num.toString()
                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "multipole":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        num = tokens[3].toInteger()
                        num+=atomCount;
                        tokens[3]=num.toString()

                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "polarize":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        for( int i =4; i<tokens.length; i++){
                            num = tokens[i].toInteger()
                            num+=atomCount;
                            tokens[i]=num.toString()
                        }

                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    case "pitors":
                        num = tokens[1].toInteger()
                        num+=atomCount;
                        tokens[1]=num.toString()

                        num = tokens[2].toInteger()
                        num+=atomCount;
                        tokens[2]=num.toString()

                        for(String part in tokens){
                            line2+=part;
                            line2+="    ";
                        }
                        bwNewPatch.write(line2 + "\n")
                        break;
                    default:
                        //Assume all other lines should be copied over identically.
                        bwNewPatch.write(line + "\n")
                }
            }

            bwNewPatch.close()
        } catch (IOException e) {
            logger.severe(e.toString());
        }

        return this
    }
}

