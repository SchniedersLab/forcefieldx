package ffx.algorithms.groovy
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
import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.CrystalPotential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFileFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.parsers.XPHFilter
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
import org.biojava.nbio.structure.chem.ResidueType
import picocli.CommandLine

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The SortArc script sort Monte Carlo archive files by lambda value. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc SortArc [options] &lt;structures1&gt &lt;structures2&gt;
 */

@Command(description = " Unwind .ARC files for nWindows", name = "ffxc SortArc")
class SortArc extends AlgorithmsScript {

    @Mixin
    private AlchemicalOptions alchemical

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

    @Option(names = ["--bpH", '--sortByPH'], paramLabel = "false",
            description = "If set, sort archive files by pH values")
    private boolean sortPh = false

    @Option(names = ['--pH'], paramLabel = "7.4",
            description = "Sets the middle of the pH ladder")
    private double pH = 7.4

    @Option(names = ['--pHGaps'], paramLabel = "1",
            description = "Sets the size of the gaps in the pH latter")
    private double pHGap = 1

    @Option(names = ["--ex", '--exponent'], paramLabel = "0.5",
            defaultValue = "0.5",
            description = "Sets the exponent for the exponential temperature ladder if sorting by temperature.")
    private double exponent = 0.05

    /**
     * The final argument(s) should be filenames for lambda windows in order..
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
    List<String> filenames = null

    private double[] lambdaValues
    private double[] temperatureValues
    private double[] pHValues
    private SystemFilter[] openers
    private SystemFilter[][] writers
    private String[] files
    private CompositeConfiguration additionalProperties
    private List<String> windowFiles = new ArrayList<>()
    MolecularAssembly[] molecularAssemblies
    MolecularAssembly ma
    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail


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
    SortArc() {
        this(new Binding())
    }

    /**
     * SortArc Constructor.
     * @param binding The Groovy Binding to use.
     */
    SortArc(Binding binding) {
        super(binding)
    }

    @Override
    SortArc run() {
        logger.info(" Running")

        if (!init()) {
            return this
        }

        int nMolAssemblies = filenames.size()
        files = new String[nMolAssemblies]
        for (int i = 0; i < nMolAssemblies; i++) {
            files[i] = filenames.get(i)
        }

        if (nWindows != -1) {
            for (int i = 0; i < nWindows; i++) {
                for (int j = 0; j < nMolAssemblies; j++) {
                    String fullPathToFile = FilenameUtils.getFullPath(files[j])
                    String directoryFullPath = fullPathToFile.replace(files[j], "") + i
                    windowFiles.add(directoryFullPath + File.separator + i)
                }
            }

            lambdaValues = new double[nWindows]
            temperatureValues = new double[nWindows]
            pHValues = new double[nWindows]
            for (int i = 0; i < nWindows; i++) {
                if (sortTemp) {
                    temperatureValues[i] = lowTemperature * Math.exp(exponent * i);
                } else if(sortPh){
                    double range = nWindows * pHGap
                    double pHMin = pH - range/2
                    if(nWindows % 2 != 0){
                        pHMin += pHGap/2
                    }
                    pHValues[i] = pHMin + i * pHGap

                } else {
                    lambdaValues[i] = alchemical.getInitialLambda(nWindows, i, false);
                }
            }
        }

        if (filenames == null) {
            return this
        }

        String[][] archiveFullPaths = new String[nWindows][nMolAssemblies]
        File file = new File(files[0])
        String directoryPath = file.getAbsoluteFile().getParent() + File.separator
        String[][] archiveNewPath = new String[nWindows][nMolAssemblies]
        File[][] saveFile = new File[nWindows][nMolAssemblies]
        File[][] arcFiles = new File[nWindows][nMolAssemblies]


        for (int j = 0; j < nMolAssemblies; j++) {
            String archiveName = FilenameUtils.getBaseName(files[j]) + ".arc"
            for (int i = 0; i < nWindows; i++) {
                archiveFullPaths[i][j] = directoryPath + i + File.separator + archiveName
                File arcFile = new File(archiveFullPaths[i][j])
                arcFiles[i][j] = arcFile
                archiveNewPath[i][j] = directoryPath + i + File.separator + FilenameUtils.getBaseName(files[j]) + "_E" + i + ".arc"
                saveFile[i][j] = new File(archiveNewPath[i][j])
            }
        }

        if(!sortPh) {
            openers = new XYZFilter[nMolAssemblies]
            writers = new XYZFilter[nWindows][nMolAssemblies]
        } else{
            openers = new XPHFilter[nMolAssemblies]
            writers = new XPHFilter[nWindows][nMolAssemblies]
        }
        int numParallel = topology.getNumParallel(threadsAvail, nMolAssemblies)
        threadsPer = (int) (threadsAvail / numParallel)


        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
    The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nMolAssemblies == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OST parameters than absolute free energies.
        if (nMolAssemblies >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        molecularAssemblies = new MolecularAssembly[nMolAssemblies]
        for (int j = 0; j < nMolAssemblies; j++) {

            if(filenames[j].contains(".pdb")){
                ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, archiveFullPaths[0][j], j)
            } else {
                ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[j], j)
            }

            molecularAssemblies[j] = ma
            ExtendedSystem extendedSystem = new ExtendedSystem(molecularAssemblies[j], null)
            //TODO: Figure out how to get esv to be properly read in
            logger.info(" Number of ESVs: " + extendedSystem.getExtendedResidueList())

            if(sortPh){
                openers[j] = new XPHFilter(algorithmFunctions.getFilter(), extendedSystem)
            } else{
                openers[j] = algorithmFunctions.getFilter()
            }

            for (int i = 0; i < nWindows; i++) {
                File arc = saveFile[i][j]
                if(sortPh) {
                    writers[i][j] = new XPHFilter(arc, molecularAssemblies[j], molecularAssemblies[j].getForceField(), additionalProperties, extendedSystem)
                } else{
                    writers[i][j] = new XYZFilter(arc, molecularAssemblies[j], molecularAssemblies[j].getForceField(), additionalProperties)
                }
            }
        }

        double tolerance
        if (sortTemp){
            tolerance = 1.0e-2
        } else if(sortPh) {
            tolerance = 1.0e-1
        } else {
            tolerance = 1.0e-4
        }

        for (int j = 0; j < nMolAssemblies; j++) {

            for (int i = 0; i < nWindows; i++) {

                logger.info(format(" Initializing %d topologies for each end", nMolAssemblies))
                openers[j].setFile(arcFiles[i][j])
                molecularAssemblies[j].setFile(arcFiles[i][j])
                logger.info("Set file to:" + arcFiles[i][j].toString())


                int snapshots = openers[j].countNumModels()
                logger.info(snapshots.toString())

                for (int n = 0; n < snapshots; n++) {
                    boolean resetPosition = n == 0

                    //TODO: Fix ReadNex to actually read in esv
                    openers[j].readNext(resetPosition, true)
                    if(sortPh) {
                        ExtendedSystem esv = (openers[j] as XPHFilter).getExtendedSystem()
                        logger.info(" ESV Residues: " + esv.getExtendedResidueList().toString())
                    }

                    String remarkLine = openers[j].getRemarkLines()


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
                            if (tokens[p].startsWith("pH")){
                                pH = Double.parseDouble(tokens[p + 1])
                            }
                        }

                    }

                    double diff
                    for (int k = 0; k < nWindows; k++) {
                        if (sortTemp) {
                            diff = Math.abs(temperatureValues[k] - temp)
                        }else if(sortPh){
                            diff = Math.abs(pHValues[k] - pH)
                        } else {
                            diff = Math.abs(lambdaValues[k] - lambda)
                        }

                        if (diff < tolerance) {
                            logger.info(" Writing to XYZ")
                            writers[k][j].writeFile(saveFile[k][j], true, remarkLine)
                            //set topology back to archive being read in
                            molecularAssemblies[j].setFile(arcFiles[i][j])
                            break
                        }
                    }
                }
            }
        }
    }
}


