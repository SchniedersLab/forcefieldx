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

import ffx.potential.Utilities
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Cluster script clusters structures by RMSD.
 * <br>
 * Usage:
 * <br>
 * ffxc Cluster [options] &lt;filename&gt;
 */
@Command(description = " Cluster structures using an RMSD matrix.", name = "ffxc Cluster")
class Cluster extends PotentialScript {

    /**
     * -a or --algorithm Clustering algorithm to use.
     */
    @Option(names = ['-a', '--algorithm'], paramLabel = "kmeans",
            description = "Print out a file with density adjusted to match mean calculated density")
    private String algorithm = "kmeans";

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The RMSD matrix.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    Cluster run() {
        if (!init()) {
            return this
        }

        if (filenames == null || filenames.isEmpty()) {
            logger.info(helpString())
            return this
        }
        File file = new File(filenames.get(0));
        int nDim = 0;
        double[][] distMatrix;
        // TODO: Read in the RMSD matrix.
        try {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String data = br.readLine();
            // Check for blank lines at the top of the file
            while (data != null && data.trim().equals("")) {
                data = br.readLine();
            }
            if (data == null) {
                logger.severe("No Data in input file.");
            }
            String[] tokens = data.trim().split("\t");
            // Expect a n x n matrix of distance values.
            nDim = tokens.size();
            distMatrix = new double[nDim][nDim];
            for(int i=0; i<nDim; i++){
                for(int j=0; j<nDim; j++){
                    distMatrix[i][j] = tokens[j].toDouble();
                }
                data = br.readLine();
                if (data != null) {
                    logger.info("Next Line:");
                    tokens = data.trim().split("\t");
                }
            }
            br.close();
            fr.close();
        } catch (IOException e) {
            logger.severe(e.toString());
        }
        if(distMatrix == null){
            logger.severe("Input read attempt failed.");
        }
        String tempString="";
        for(int i =0; i<nDim; i++){
            for(int j = 0; j<nDim; j++){
                tempString+=String.format("%f \t", distMatrix[i][j]);
            }
            tempString+="\n";
        }
        logger.info(tempString);
        // TODO: Input the RMSD matrix to the clustering algorithm
        // Use the org.apache.commons.math3.ml.clustering package.

        // TODO: Output the clusters in a useful way.

        return this
    }
}
