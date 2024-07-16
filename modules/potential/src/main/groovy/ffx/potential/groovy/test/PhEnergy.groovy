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
package ffx.potential.groovy.test

import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.potential.utils.StructureMetrics.momentsOfInertia
import static ffx.potential.utils.StructureMetrics.radiusOfGyration
import static java.lang.Double.parseDouble
import static java.lang.Double.parseDouble
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*
import static org.apache.commons.math3.util.FastMath.pow

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc test.PhEnergy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy for a CpHMD system.",
    name = "test.PhEnergy")
class PhEnergy extends PotentialScript {

    @Mixin
    AtomSelectionOptions atomSelectionOptions

    /**
     * -m or --moments print out electrostatic moments.
     */
    @Option(names = ['-m', '--moments'], paramLabel = "false", defaultValue = "false",
            description = 'Print out electrostatic moments.')
    private boolean moments = false

    /**
     * --rg or --gyrate Print out the radius of gyration.
     */
    @Option(names = ['--rg', '--gyrate'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the radius of gyration.')
    private boolean gyrate = false

    /**
     * --in or --inertia Print out the moments of inertia.
     */
    @Option(names = ['--in', '--inertia'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the moments of inertia.')
    private boolean inertia = false

    /**
     * -g or --gradient to print out gradients.
     */
    @Option(names = ['-g', '--gradient'], paramLabel = "false", defaultValue = "false",
            description = 'Compute the atomic gradient as well as energy.')
    private boolean gradient = false

    /**
     * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
     */
    @Option(names = ['--fl', '--findLowest'], paramLabel = "0", defaultValue = "0",
            description = 'Return the n lowest energies from an ARC/PDB file.')
    private int fl = 0

    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @Option(names = ['-v', '--verbose'], paramLabel = "false", defaultValue = "false",
            description = "Print out all energy components for each snapshot.")
    private boolean verbose = false

    /**
     * --dc or --densityCutoff Collect structures above a specified density.
     */
    @Option(names = ['--dc', '--densityCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures above a specified density.")
    private double dCutoff = 0.0

    /**
     * --ec or --energyCutOff Collect structures below a specified energy range from the minimum energy.
     */
    @Option(names = ['--ec', '--energyCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
    private double eCutoff = 0.0

    /**
     * --pH or --constantPH Constant pH value for the test.
     */
    @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
            description = 'pH value for the energy evaluation. (Only applies when esvTerm is true)')
    double pH = 7.4

    @Option(names = ['--aFi', '--arcFile'], paramLabel = "traj",
            description = 'A file containing snapshots to evaluate on when using a PDB as a reference to build from. There is currently no default.')
    private String arcFileName = null

    @Option(names = ["--bar", "--mbar"], paramLabel = "false",
            description = "Run (restartable) energy evaluations for MBAR. Requires an ARC file to be passed in. Set the tautomer flag to true for tautomer parameterization.")
    boolean mbar = false

    @Option(names = ["--numLambda", "--nL", "--nw"], paramLabel = "-1",
            description = "Required for lambda energy evaluations. Ensure numLambda is consistent with the trajectory lambdas, i.e. gaps between traj can be filled easily. nL >> nTraj is recommended.")
    int numLambda = -1

    @Option(names = ["--outputDir", "--oD"], paramLabel = "",
            description = "Where to place MBAR files. Default is ../mbarFiles/energy_(window#).mbar. Will write out a file called energy_0.mbar.")
    String outputDirectory = ""

    @Option(names = ["--lambdaDerivative", "--lD"], paramLabel = "false",
            description = "Perform dU/dL evaluations and save to mbarFiles.")
    boolean derivatives = false

    @Option(names = ["--perturbTautomer"], paramLabel = "false",
            description = "Change tautomer instead of lambda state for MBAR energy evaluations.")
    boolean tautomer = false

    /**
     * --testEndStateEnergies
     */
    @Option(names = ['--testEndStateEnergies'], paramLabel = 'false',
            description = 'Test both ESV energy end states as if the polarization damping factor is initialized from the respective protonated or deprotonated state')
    boolean testEndstateEnergies = false

    @Option(names = ['--recomputeAverage'], paramLabel = 'false',
            description = 'Recompute average position and spit out said structure from trajectory')
    boolean recomputeAverage = false

    /**
     * The final argument is a PDB or XPH coordinate file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XPH format.')
    private String filename = null

    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null
    public MolecularAssembly mola = null

    /**
     * Energy constructor.
     */
    PhEnergy() {
        this(new Binding())
    }

    /**
     * Energy constructor.
     * @param binding The Groovy Binding to use.
     */
    PhEnergy(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    PhEnergy run() {
        // Init the context and bind variables.
        if (!init()) {
            return this
        }
        if (mbar){ // This is probably set to true since parameterization requires locked lambda states
            System.setProperty("lock.esv.states", "false")
        }

        // Load the MolecularAssembly.
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        // Set the filename.
        filename = activeAssembly.getFile().getAbsolutePath()

        logger.info("\n Running Energy on " + filename)
        forceFieldEnergy = activeAssembly.getPotentialEnergy()

        // Restart File
        File esv = new File(removeExtension(filename) + ".esv")
        if (!esv.exists()) {
            esv = null
        }

        ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, pH, esv)
        if(testEndstateEnergies && BigDecimal.valueOf(esvSystem.getExtendedLambdas()[0]) ==  0.0){
            for(Atom atom : esvSystem.getExtendedAtoms()){
                int atomIndex = atom.getArrayIndex()
                if (esvSystem.isTitratingHeavy(atomIndex)) {
                    double endstatePolar = esvSystem.titrationUtils.getPolarizability(atom, 0.0, 0.0, atom.getPolarizeType().polarizability)
                    double sixth = 1.0 / 6.0
                    atom.getPolarizeType().pdamp = pow(endstatePolar, sixth)
                }
            }
        } else if(testEndstateEnergies && BigDecimal.valueOf(esvSystem.getExtendedLambdas()[0]) ==  1.0){
            for(Atom atom : esvSystem.getExtendedAtoms()){
                int atomIndex = atom.getArrayIndex()
                if (esvSystem.isTitratingHeavy(atomIndex)) {
                    double endstatePolar = esvSystem.titrationUtils.getPolarizability(atom, 1.0, 0.0, atom.getPolarizeType().polarizability)
                    double sixth = 1.0 / 6.0
                    atom.getPolarizeType().pdamp = pow(endstatePolar, sixth)
                }
            }
        }
        esvSystem.setConstantPh(pH)
        int numESVs = esvSystem.extendedResidueList.size()
        forceFieldEnergy.attachExtendedSystem(esvSystem)
        logger.info(format(" Attached extended system with %d residues.", numESVs))

        // Apply atom selections
        atomSelectionOptions.setActiveAtoms(activeAssembly)

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)
        double[] averageCoordinates = Arrays.copyOf(x, x.length)
        if (gradient) {
            double[] g = new double[nVars]
            int nAts = (int) (nVars / 3)
            energy = forceFieldEnergy.energyAndGradient(x, g, true)
            logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
            for (int i = 0; i < nAts; i++) {
                int i3 = 3 * i
                logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
            }
        } else {
            energy = forceFieldEnergy.energy(x, true)
        }

        if (moments) {
            forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false)
        }

        if (gyrate) {
            double rg = radiusOfGyration(activeAssembly.getActiveAtomArray())
            logger.info(format(" Radius of gyration:           %10.5f A", rg))
        }

        if (inertia) {
            double[][] inertiaValue = momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true)
        }

        SystemFilter systemFilter
        if(arcFileName != null){
            File arcFile = new File(arcFileName)
            systemFilter = new XPHFilter(arcFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties(), esvSystem)
        } else{
            systemFilter = potentialFunctions.getFilter()
            if(systemFilter instanceof XYZFilter){
                systemFilter = new XPHFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties(), esvSystem)
                systemFilter.readFile()
                logger.info("Reading ESV lambdas from XPH file")
                forceFieldEnergy.getCoordinates(x)
                forceFieldEnergy.energy(x, true)
            } 
        }

        if (mbar){
            computeESVEnergiesAndWriteFile(systemFilter, esvSystem)
            return this
        }

        if (systemFilter instanceof XPHFilter || systemFilter instanceof PDBFilter) {
            int index = 1
            while (systemFilter.readNext()) {
                index++
                Crystal crystal = activeAssembly.getCrystal()
                forceFieldEnergy.setCrystal(crystal)
                forceFieldEnergy.getCoordinates(x)
                if(recomputeAverage){
                    for(int i = 0; i < x.length; i++){
                        averageCoordinates[i] += x[i]
                    }
                }
                if (verbose) {
                    logger.info(format(" Snapshot %4d", index))
                    if (!crystal.aperiodic()) {
                        logger.info(format("\n Density:                                %6.3f (g/cc)",
                                crystal.getDensity(activeAssembly.getMass())))
                    }
                    energy = forceFieldEnergy.energy(x, true)
                } else {
                    energy = forceFieldEnergy.energy(x, false)
                    logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))
                }
            }
            if(recomputeAverage){
                for(int i = 0; i < x.length; i++){
                    x[i] = averageCoordinates[i] / index
                }
                forceFieldEnergy.setCoordinates(x)
            }
        }
        if(recomputeAverage){
            String modelFilename = activeAssembly.getFile().getAbsolutePath()
            if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
                baseDir = new File(FilenameUtils.getFullPath(modelFilename))
            }
            String dirName = baseDir.toString() + File.separator
            String fileName = FilenameUtils.getName(modelFilename)
            String ext = FilenameUtils.getExtension(fileName)
            fileName = FilenameUtils.removeExtension(fileName)
            File saveFile
            SystemFilter writeFilter
            if (ext.toUpperCase().contains("XYZ")) {
                saveFile = new File(dirName + fileName + ".xyz")
                writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                        activeAssembly.getProperties())
                potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
            } else if (ext.toUpperCase().contains("ARC")) {
                saveFile = new File(dirName + fileName + ".arc")
                saveFile = potentialFunctions.versionFile(saveFile)
                writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                        activeAssembly.getProperties())
                potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
            } else {
                saveFile = new File(dirName + fileName + ".pdb")
                saveFile = potentialFunctions.versionFile(saveFile)
                writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                        activeAssembly.getProperties())
                int numModels = systemFilter.countNumModels()
                if (numModels > 1) {
                    writeFilter.setModelNumbering(0)
                }
                writeFilter.writeFile(saveFile, true, false, false)
            }
        }


        return this
    }

    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (forceFieldEnergy == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList((Potential) forceFieldEnergy)
        }
        return potentials
    }

    private class StateContainer implements Comparable<StateContainer> {

        private final AssemblyState state
        private final double e

        StateContainer(AssemblyState state, double e) {
            this.state = state
            this.e = e
        }

        AssemblyState getState() {
            return state
        }

        double getEnergy() {
            return e
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy())
        }
    }

    void computeESVEnergiesAndWriteFile(SystemFilter systemFilter, ExtendedSystem esvSystem) {
        // Find all run directories to determine lambda states & make necessary files
        File mbarFile;
        File mbarGradFile;
        double[] lambdas
        if (outputDirectory.isEmpty()) {
            File dir = new File(filename).getParentFile()
            File parentDir = dir.getParentFile()
            int thisRung = -1
            dir.getName().find(/(\d+)/) { match ->
                thisRung = match[0].toInteger()
            }
            assert thisRung != -1: "Could not determine the rung number from the directory name."
            mbarFile = new File(parentDir.getAbsolutePath() + File.separator + "mbarFiles" + File.separator + "energy_"
                    + thisRung + ".mbar")
            mbarGradFile = new File(parentDir.getAbsolutePath() + File.separator + "mbarFiles" + File.separator + "derivative_"
                    + thisRung + ".mbar")
            mbarFile.getParentFile().mkdir()
            File[] lsFiles = parentDir.listFiles()
            List<File> rungFiles = new ArrayList<>()
            for (File file : lsFiles) {
                if (file.isDirectory() && file.getName().matches(/\d+/)) {
                    rungFiles.add(file)
                }
            }
            if (numLambda == -1) {
                numLambda = rungFiles.size()
            }
            lambdas = new double[numLambda]
            for (int i = 0; i < numLambda; i++) {
                double dL = 1 / (numLambda - 1)
                lambdas[i] = i * dL
            }


            logger.info(" Computing energies for each lambda state for generation of mbar file.")
            logger.info(" MBAR File: " + mbarFile)
            logger.info(" Lambda States: " + lambdas)
        } else {
            mbarFile     = new File(outputDirectory + File.separator + "energy_0.mbar")
            mbarGradFile = new File(outputDirectory + File.separator + "derivative_0.mbar")
            lambdas = new double[numLambda]
            if (numLambda == -1){
                logger.severe("numLambda must be set when outputDirectory is set.")
            }
            for (int i = 0; i < numLambda; i++) {
                double dL = 1 / (numLambda - 1)
                lambdas[i] = i * dL
            }
        }

        int progress = 1
        if (mbarFile.exists()){ // Restartable
            // Count lines in mbarFile
            progress = mbarFile.readLines().size() - 1 // Subtract header
            for(int i = 0; i < progress; i++){
                systemFilter.readNext()
            }
            progress += 1 // Increment to next snapshot
            logger.info("\n Restarting MBAR file at snapshot " + progress)
        }

        if (systemFilter instanceof XPHFilter || systemFilter instanceof PDBFilter) {
            int index = progress
            double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
            try(FileWriter fw = new FileWriter(mbarFile, mbarFile.exists())
                BufferedWriter writer = new BufferedWriter(fw)
                FileWriter fwGrad = new FileWriter(mbarGradFile, mbarGradFile.exists())
                BufferedWriter writerGrad = new BufferedWriter(fwGrad)
            ){
                // Write header (Number of snaps, temp, and this.xyz)
                StringBuilder sb = new StringBuilder(systemFilter.countNumModels() + "\t" + "298.0" + "\t" + getBaseName(filename))
                StringBuilder sbGrad = new StringBuilder(systemFilter.countNumModels() + "\t" + "298.0" + "\t" + getBaseName(filename))
                logger.info(" MBAR file temp is hardcoded to 298.0 K. Please change if necessary.")
                sb.append("\n")
                sbGrad.append("\n")
                if (progress == 1) { // File didn't exist (more consistent than checking for existence)
                    writer.write(sb.toString())
                    writer.flush()
                    logger.info(" Header: " + sb.toString())
                    if(derivatives){
                       writerGrad.write(sbGrad.toString())
                       writerGrad.flush()
                       logger.info(" Header: " + sbGrad.toString())
                    }
                }
                while (systemFilter.readNext()) {
                    // MBAR lines (\t index\t lambda0 lambda1 ... lambdaN)
                    sb = new StringBuilder("\t" + index + "\t")
                    sbGrad = new StringBuilder("\t" + index + "\t")
                    index++
                    Crystal crystal = activeAssembly.getCrystal()
                    forceFieldEnergy.setCrystal(crystal)
                    for(double lambda : lambdas) {
                        if (tautomer){
                            setESVTautomer(lambda, esvSystem)
                        } else {
                            setESVLambda(lambda, esvSystem)
                        }
                        forceFieldEnergy.getCoordinates(x)
                        if (derivatives) {
                            energy = forceFieldEnergy.energyAndGradient(x, new double[x.length * 3])
                            double grad = esvSystem.getDerivatives()[0] // Only one residue
                            sbGrad.append(grad).append(" ")
                        } else {
                            energy = forceFieldEnergy.energy(x, false)
                        }
                        sb.append(energy).append(" ")
                    }
                    sb.append("\n")
                    writer.write(sb.toString())
                    writer.flush() // Flush after each snapshot, otherwise it doesn't do it consistently
                    logger.info(sb.toString())
                    if (derivatives) {
                        sbGrad.append("\n")
                        writerGrad.write(sbGrad.toString())
                        writerGrad.flush()
                        logger.info(sbGrad.toString())
                    }
                }
            } catch (IOException e) {
                logger.severe("Error writing to MBAR file.")
            }
        }
    }

    /**
     * Sets lambda values for the extended system. Note that it is expected that the
     * tautomer is set correctly from dynamics.
     * @param lambda
     * @param extendedSystem
     */
    static void setESVLambda(double lambda, ExtendedSystem extendedSystem) {
        List<Residue> residueList = extendedSystem.getExtendedResidueList()
        if (residueList.size() == 1 || (residueList.size() == 2 && extendedSystem.isTautomer(residueList.get(0)))){
            extendedSystem.setTitrationLambda(residueList.get(0), lambda, false);
        } else {
            if (residueList.size() == 0) {
                logger.severe("No residues found in the extended system.")
            } else {
                logger.severe("Only one lambda path is allowed for MBAR energy evaluations.")
            }
        }
    }

    /**
     * Sets tautomer values for the extended system. Note that it is expected that the
     * lambda is set correctly from dynamics.
     * @param tautomer
     * @param extendedSystem
     */
    static void setESVTautomer(double tautomer, ExtendedSystem extendedSystem) {
        List<Residue> residueList = extendedSystem.getExtendedResidueList()
        if (residueList.size() == 1 || (residueList.size() == 2 && extendedSystem.isTautomer(residueList.get(0)))) {
            extendedSystem.setTautomerLambda(residueList.get(0), tautomer, false);
        } else {
            if (residueList.size() == 0) {
                logger.severe("No residues found in the extended system.")
            } else {
                logger.severe("Only one lambda path is allowed for MBAR energy evaluations.")
            }
        }
    }
}

