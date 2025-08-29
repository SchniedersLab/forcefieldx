//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.ConformationScan
import ffx.algorithms.optimize.Minimize
import ffx.algorithms.optimize.TorsionSearch
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Molecule
import ffx.potential.utils.PotentialsUtils
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CoformerBindingSearch is a Groovy script that generates a set of molecular orientations in vacuum and
 * calculates the energy of each conformation.
 * <br>
 * Usage:
 * <br>
 * ffxc test.CoformerBindingSearch [options] &lt;filename&gt;
 */
@Command(description = " Calculates interaction energies of different molecular orientations and saves low energy orientations.",
        name = "test.CoformerBindingSearch")
class CoformerBindingSearch extends AlgorithmsScript {

    /**
     * --eps
     */
    @Option(names = ['--eps'], paramLabel = '.01',
            description = 'Gradient cutoff for minimization.')
    double eps = 0.01

    /**
     * --maxIter
     */
    @Option(names = ['--maxIter'], paramLabel = '10000',
            description = 'Max iterations for minimization.')
    int maxIter = 10000

    /**
     * --solventDielectric
     */
    @Option(names = ['--solventDielectric'], paramLabel = '78.4',
            description = 'Sets the gk solvent dielectric constant.')
    double gkSolventDielec = 78.4

    /**
     * --skipHomodimerNumber
     */
    @Option(names = ['--skipHomodimerNumber'], paramLabel = '-1',
            description = 'Skip conformation search on input dimer one or two.')
    int skipHomodimerNumber = -1

    /**
     * --torsionScan
     */
    @Option(names = ['--torsionScan', "--tscan"], paramLabel = 'false', defaultValue = 'false',
            description = 'During sampling, statically scan torsions after direct minimization to find the lowest energy conformation.')
    private boolean tscan = false

    /**
     * --noMinimize
     */
    @Option(names = ['--noMinimize', '--noMin'], paramLabel = 'false', defaultValue = 'false',
            description = 'Don\'t minimize or torsion scan after conformations are generated. Useful for testing.')
    private boolean noMinimize = false

    /**
     * --excludeH
     */
    @Option(names = ['--excludeH', "--eh"], paramLabel = 'false', defaultValue = 'false',
            description = 'Only include H bonded to electronegative atoms in conformations.')
    private boolean excludeH = false

    /**
     * --coformerOnly
     */
    @Option(names = ['--coformerOnly'], paramLabel = 'false', defaultValue = 'false',
            description = 'Only conformation search the coformer.')
    private boolean coformerOnly = false

    /**
     * --gk
     */
    @Option(names = ['--gk'], paramLabel = 'false', defaultValue = 'false',
            description = 'Use generalized kirkwood solvent.')
    private boolean gk = false

    /**
     * Filename.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ input file.")
    private List<String> filenames

    boolean minimize = null


    /**
     * Constructor.
     */
    CoformerBindingSearch() {
        super()
    }

    /**
     * Constructor.
     * @param binding The Groovy Binding to use.
     */
    CoformerBindingSearch(Binding binding) {
        super(binding)
    }

    /**
     * Constructor that sets the command line arguments.
     * @param args Command line arguments.
     */
    CoformerBindingSearch(String[] args) {
        super(args)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    CoformerBindingSearch run() {
        System.setProperty("direct-scf-fallback", "true")

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // Cat the key files together and set the -Dkey property to be the new file we created
        // Write default gk options to the key file
        // Only does a simple search for the patch file so it needs to be named accordingly with the .xyz
        setKeyAndPatchFiles(gk, gkSolventDielec, filenames)

        // Check the size of the filenames list
        if (!(filenames.size() == 1 || filenames.size() == 2)) {
            logger.severe("Must provide one or two filenames.")
            return this
        }

        minimize = !noMinimize
        boolean skipMoleculeOne = skipHomodimerNumber == 1
        boolean skipMoleculeTwo = skipHomodimerNumber == 2

        // Perform scan on monomer one w/ itself
        ConformationScan systemOneScan = null
        ConformationScan systemTwoScan = null
        if(!coformerOnly) {
            if (!skipMoleculeOne) {
                systemOneScan = runScan(filenames.get(0), filenames.get(0))
            } else{
                logger.info(" Skipping monomer one scan.")
            }

            if (!skipMoleculeTwo && filenames.size() == 2) {
                systemTwoScan = runScan(filenames.get(1), filenames.get(1))
            } else if (!skipMoleculeTwo && filenames.size() == 1) {
                logger.info(" Only one file provided, skipping second homodimer scan.")
            }
        } else {
            logger.info(" Skipping monomer one and two scan.")
        }

        if(filenames.size() == 2){
            ConformationScan bothSystems = runScan(filenames.get(0), filenames.get(1))
            if(systemOneScan != null && systemTwoScan != null){
                logger.info("\n System one (" + FilenameUtils.removeExtension(filenames.get(0)) + ") self scan energy information:")
                systemOneScan.logAllEnergyInformation()
                logger.info("\n System two (" + FilenameUtils.removeExtension(filenames.get(1)) + ") self scan energy information:")
                systemTwoScan.logAllEnergyInformation()
                ConformationScan.logBindingEnergyCalculation(systemOneScan, systemTwoScan, bothSystems)
            }
        }
        else {
            logger.info(" Only one file provided, skipping coformer scan.")
        }

        return this
    }

    ConformationScan runScan(String fileOne, String fileTwo){
        boolean coformer = fileOne != fileTwo
        PotentialsUtils potentialsUtils = new PotentialsUtils()
        MolecularAssembly[] molecularAssemblies = potentialsUtils.openAll(new String[]{fileOne, fileTwo})
        minimizeMolecularAssemblies(molecularAssemblies) // Only if noMinimize is false
        int secondSystemStartIndex = molecularAssemblies[0].getMoleculeArray().length
        MolecularAssembly combined = combineTwoMolecularAssemblies(molecularAssemblies[0], molecularAssemblies[1])
        // Split combined mola into two lists based on molecules in files
        ArrayList<Molecule> setOne = new ArrayList<>()
        ArrayList<Molecule> setTwo = new ArrayList<>()
        for(int i = 0; i < combined.getMoleculeArray().length; i++){
            if(i < secondSystemStartIndex){
                setOne.add(combined.getMoleculeArray()[i])
            } else{
                setTwo.add(combined.getMoleculeArray()[i])
            }
        }
        Molecule[] systemOne = setOne.toArray(new Molecule[0])
        Molecule[] systemTwo = setTwo.toArray(new Molecule[0])
        ConformationScan scan = new ConformationScan(
                combined,
                systemOne,
                systemTwo,
                eps,
                maxIter,
                tscan,
                excludeH,
                minimize
        )
        scan.scan()
        String fileName = coformer ? "coformerScan.arc" : FilenameUtils.removeExtension(fileOne) + ".arc"
        File file = new File(fileName)
        if(!scan.writeStructuresToXYZ(file)){
            logger.warning(" No structures saved from scan.")
            scan = null
        } else if (!coformer){
            logger.info("\n System (" + FilenameUtils.removeExtension(fileOne) + ") self scan energy information:")
            scan.logAllEnergyInformation()
        } else{
            logger.info("\n System (" + FilenameUtils.removeExtension(fileOne) +
                    ") and system (" + FilenameUtils.removeExtension(fileTwo) + ") scan energy information:")
            scan.logAllEnergyInformation()
        }
        return scan
    }

    static MolecularAssembly combineTwoMolecularAssemblies(MolecularAssembly mola1, MolecularAssembly mola2){
        MolecularAssembly mainMonomerAssembly = mola1
        MolecularAssembly feederAssembly = mola2
        Molecule[] assemblyTwoMolecules = feederAssembly.getMoleculeArray()
        int molNum = mainMonomerAssembly.getMoleculeArray().length
        for(Molecule m: assemblyTwoMolecules) {
            for (Atom a : m.getAtomList()) {
                a.setMoleculeNumber(molNum)
            }
            m.setName("AddedMol" + molNum)
            mainMonomerAssembly.addMSNode(m)
            molNum++
        }
        mainMonomerAssembly.update()
        mainMonomerAssembly.setPotential(null) // energyFactory doesn't do anything if it isn't null
        ForceFieldEnergy forceFieldEnergy = ForceFieldEnergy.energyFactory(mainMonomerAssembly)
        return mainMonomerAssembly
    }

    static void setKeyAndPatchFiles(boolean gk, double gkSolventDielec, List<String> filenames) {
        String key = "coformerScan.key"
        String patch = "coformerScan.patch"
        // Create the key file
        File keyFile = new File(key)
        File patchFile = new File(patch)
        logger.info(" Creating key file: " + key)
        keyFile.createNewFile()
        // concatenate the two files together with bufferedReader and bufferedWriter
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(keyFile))
            bw.write("patch " + FilenameUtils.removeExtension(key) + ".patch\n\n")
            if(gk){
                bw.write("gkterm true\n")
                bw.write("solvent-dielectric " + gkSolventDielec + "\n")
                bw.write("gk-radius solute\n")
                bw.write("cavmodel gauss-disp\n")
            }
            bw.close()
        } catch (IOException e) {
            e.printStackTrace()
        }
        logger.info(" Creating patch file: " + patch)
        patchFile.createNewFile()
        String patchOne = FilenameUtils.removeExtension(filenames.get(0)) + ".patch"
        String patchTwo = FilenameUtils.removeExtension(filenames.get(1)) + ".patch"
        String[] files = new String[]{patchOne, patchTwo}
        BufferedWriter bw = new BufferedWriter(new FileWriter(patchFile))
        for (String file : files) {
            BufferedReader br = new BufferedReader(new FileReader(file))
            String line = br.readLine()
            while (line != null) {
                bw.write(line + "\n")
                line = br.readLine()
            }
            br.close()
        }
        bw.close()
        System.setProperty("key", key)
    }

    void minimizeMolecularAssemblies(MolecularAssembly[] molecularAssemblies) {
        if (!noMinimize) {
            logger.info(" Minimizing molecular assemblies.")
            for (MolecularAssembly molecularAssembly : molecularAssemblies) {
                for(Molecule m: molecularAssembly.getMoleculeArray()) {
                    if (tscan) {
                        TorsionSearch ts = new TorsionSearch(molecularAssembly, m, 32, 1)
                        ts.staticAnalysis(0, 100)
                        if (!ts.getStates().isEmpty()) {
                            AssemblyState minState = ts.getStates().get(0)
                            minState.revertState()
                        }
                    }
                }
                Minimize minimizer = new Minimize(molecularAssembly, molecularAssembly.getPotentialEnergy(), algorithmListener)
                minimizer.minimize(eps, maxIter)
            }
            logger.info(" Done minimizing molecular assemblies.")
        }
    }
}

