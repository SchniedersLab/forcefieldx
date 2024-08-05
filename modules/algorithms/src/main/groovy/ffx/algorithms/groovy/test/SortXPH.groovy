package ffx.algorithms.groovy.test

import edu.rit.pj.ParallelTeam

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

import ffx.algorithms.cli.AlgorithmsScript
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The SortXPH command unwinds .ARC files from CpHMD runs.
 * <br>
 * Usage:
 * <br>
 * ffxc test.SortXPH [options] &lt;structures1&gt &lt;structures2&gt;
 */

@Command(description = " The SortXPH command unwinds .ARC files from CpHMD runs.", name = "test.SortXPH")
class SortXPH extends AlgorithmsScript {

    @Mixin
    private TopologyOptions topology

    @Option(names = ["--nw", "--nWindows"], paramLabel = "-1",
            description = "If set, auto-determine lambda values and subdirectories (overrides other flags).")
    private int nWindows = -1

    @Option(names = ["--bT", '--sortByTemp'], paramLabel = "false",
            description = "If set, sort archive files by temperature values")
    private boolean sortTemp = false

    @Option(names = ["--sT", '--startTemp'], paramLabel = "298.15",
            defaultValue = "298.15",
            description = "Sets the starting temperature for the exponential temperature ladder if sorting by temperature.")
    private double lowTemperature = 298.15

    @Option(names = ['--pH'], paramLabel = "7.4",
            description = "Sets the middle of the pH ladder")
    private double pH = 7.4

    @Option(names = ['--pHGaps'], paramLabel = "1",
            description = "Sets the size of the gaps in the pH ladder")
    private double pHGap = 1

    @Option(names = ["--ex", '--exponent'], paramLabel = "0.5",
            defaultValue = "0.5",
            description = "Sets the exponent for the exponential temperature ladder if sorting by temperature.")
    private double exponent = 0.05

    /**
     * The final argument(s) should be filenames for lambda windows in order..
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
    private String filename = null

    private double[] lambdaValues
    private double[] temperatureValues
    private double[] pHValues
    private XPHFilter opener
    private XPHFilter[] writers
    private CompositeConfiguration additionalProperties
    private List<String> windowFiles = new ArrayList<>()
    MolecularAssembly molecularAssembly
    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail
    ExtendedSystem extendedSystem


    /**
     * Sets an optional Configuration with additional properties.
     * @param additionalProps
     */
    void setProperties(CompositeConfiguration additionalProps) {
        this.additionalProperties = additionalProps
    }

    /**
     * SortArc Constructor.
     */
    SortXPH() {
        this(new Binding())
    }

    /**
     * SortArc Constructor.
     * @param binding The Groovy Binding to use.
     */
    SortXPH(Binding binding) {
        super(binding)
    }

    @Override
    SortXPH run() {
        logger.info(" Running")

        if (!init()) {
            return this
        }

        if (nWindows != -1) {
            for (int i = 0; i < nWindows; i++) {
                String fullPathToFile = FilenameUtils.getFullPath(filename)
                String directoryFullPath = fullPathToFile.replace(filename, "") + i
                windowFiles.add(directoryFullPath + File.separator + i)
            }

            lambdaValues = new double[nWindows]
            temperatureValues = new double[nWindows]
            pHValues = new double[nWindows]
            for (int i = 0; i < nWindows; i++) {
                if (sortTemp) {
                    temperatureValues[i] = lowTemperature * Math.exp(exponent * i);
                } else{
                    double range = nWindows * pHGap
                    double pHMin = pH - range/2
                    if(nWindows % 2 != 0){
                        pHMin += pHGap/2
                    }
                    pHValues[i] = pHMin + i * pHGap
                }
            }
        }

        String[] archiveFullPaths = new String[nWindows]
        File file = new File(filename)
        String directoryPath = file.getAbsoluteFile().getParent() + File.separator
        String[] archiveNewPath = new String[nWindows]
        File[] saveFile = new File[nWindows]
        File[] arcFiles = new File[nWindows]
        String archiveName = FilenameUtils.getBaseName(filename) + ".arc"

        writers = new XPHFilter[nWindows]
        for (int i = 0; i < nWindows; i++) {
            archiveFullPaths[i] = directoryPath + i + File.separator + archiveName
            File arcFile = new File(archiveFullPaths[i])
            arcFiles[i] = arcFile
            archiveNewPath[i] = directoryPath + i + File.separator + FilenameUtils.getBaseName(filename) + "_E" + i + ".arc"
            saveFile[i] = new File(archiveNewPath[i])

        }

        molecularAssembly = getActiveAssembly(filename)
        extendedSystem = new ExtendedSystem(molecularAssembly, pH, null)
        molecularAssembly.getPotentialEnergy().attachExtendedSystem(extendedSystem)


        int numParallel = topology.getNumParallel(threadsAvail, 1)
        threadsPer = (int) (threadsAvail / numParallel)

        for (int i = 0; i < nWindows; i++) {
            File arc = saveFile[i]
            writers[i] = new XPHFilter(arc, molecularAssembly, molecularAssembly.getForceField(),
                    additionalProperties, extendedSystem)
        }
        opener = new XPHFilter(algorithmFunctions.getFilter(), extendedSystem)

        double tolerance
        if (sortTemp){
            tolerance = 1.0e-2
        } else  {
            tolerance = 1.0e-1
        }


        for (int i = 0; i < nWindows; i++) {

            opener.setFile(arcFiles[i])
            molecularAssembly.setFile(arcFiles[i])
            logger.info("Set file to:" + arcFiles[i].toString())


            int snapshots = opener.countNumModels()
            logger.info(snapshots.toString())

            for (int n = 0; n < snapshots; n++) {
                boolean resetPosition = n == 0

                //TODO: Fix ReadNex to actually read in esv
                opener.readNext(resetPosition, false)

                String remarkLine = opener.getRemarkLines()


                double lambda = 0
                double temp = 0
                double pH = 0
                if (remarkLine.contains(" Lambda: ")) {
                    String[] tokens = remarkLine.split(" +")
                    for (int p = 0; p < tokens.length; p++) {
                        if (tokens[p].startsWith("Lambda")) {
                            lambda = Double.parseDouble(tokens[p + 1])
                        }
                        if (tokens[p].startsWith("Temp")) {
                            temp = Double.parseDouble(tokens[p + 1])
                        }
                        if (tokens[p].startsWith("pH")) {
                            pH = Double.parseDouble(tokens[p + 1])
                        }
                    }

                }

                double diff
                for (int k = 0; k < nWindows; k++) {
                    if (sortTemp) {
                        diff = Math.abs(temperatureValues[k] - temp)
                    } else  {
                        diff = Math.abs(pHValues[k] - pH)
                    }

                    if (diff < tolerance) {
                        double[] lambdas = opener.getExtendedSystem().getLambdaArray()
                        logger.info("************ Lambdas: " + lambdas)
                        writers[k].getExtendedSystem().setLambdaArray(lambdas)
                        logger.info(" Writing to XPH")
                        writers[k].writeFile(saveFile[k], true, remarkLine)
                        //set topology back to archive being read in
                        molecularAssembly.setFile(arcFiles[i])
                        break
                    }
                }
            }
        }
    }
}


