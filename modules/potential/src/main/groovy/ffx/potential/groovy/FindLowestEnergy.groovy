package ffx.potential.groovy

import ffx.potential.ForceFieldEnergy
import ffx.potential.AssemblyState
import ffx.potential.bonded.Atom
import ffx.potential.parsers.XYZFileFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.utils.PotentialsUtils
import org.apache.commons.io.FilenameUtils
import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.Superpose
import org.checkerframework.checker.units.qual.K
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import com.google.common.collect.MinMaxPriorityQueue

import static java.lang.String.format

/**
 * The FindLowestEnergy script calculates energies for all assemblies in an arc file , finds the lowest energy assembly,
 * and saves that assembly to a pdb file.
 * <br>
 * Usage:
 * <br>
 * ffxc FindLowestEnergy [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "ffxc FindLowestEnergy")
class FindLowestEnergy extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    /**
     * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
     */
    @Option(names = ['-K', '--nlowest'], paramLabel = "1",
            description = 'allows you to get the n lowest energy structures in an arc file')
    int numSnaps = 1

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    private double energy = Double.MAX_VALUE;
    private AssemblyState assemblyState = null;

    private class StateContainer implements Comparable<StateContainer> {
        private final AssemblyState state;
        private final double e;

        StateContainer(AssemblyState state, double e) {
            this.state = state;
            this.e = e;

        }

        AssemblyState getState() {
            return state;
        }
        double getEnergy() {
            return e;
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy());
        }
    }

    /**
     * Execute the script.
     */
    @Override
    FindLowestEnergy run() {
        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames == null || filenames.size() != 1) {
            logger.warning(" Invalid arguments provided! Must have one non-null argument.")
            logger.info(helpString())
            return this
        } else {
            assemblies = potentialFunctions.open(filenames.get(0));
            activeAssembly = assemblies[0]
        }

        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)
        energy = forceFieldEnergy.energy(x, true)
        assemblyState = new AssemblyState(activeAssembly)
        SystemFilter systemFilter = potentialFunctions.getFilter()


        int maxnum= 1

        /**
         * Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
         */

        if (numSnaps < 1 ) {
            numSnaps = 1
            logger.info(String.format(" Warning!!! Cannot request 0 lowest enrgies! System will return 1 lowest energies and PDB file"))
        }

        MinMaxPriorityQueue<StateContainer> energyQueue;
        energyQueue = MinMaxPriorityQueue
                .maximumSize(numSnaps)
                .create()
        StateContainer firstContainer;
        firstContainer = new StateContainer(assemblyState, energy)
        energyQueue.add(firstContainer)

        if (systemFilter instanceof XYZFilter) {
            XYZFilter xyzFilter = (XYZFilter) systemFilter
            //calling the next assembly of the arc file
            while (xyzFilter.readNext()) {
                forceFieldEnergy.getCoordinates(x) // getting the coordinates for the next assembly
                energy= forceFieldEnergy.energy(x, true) //calculating energy for new assembly
                assemblyState = new AssemblyState(activeAssembly) //saving new assembly if the energy is less than the current energy
                energyQueue.add(new StateContainer(assemblyState, energy))
                maxnum = maxnum + 1

            }
        } else {
            logger.severe(String.format(" System %s does not appear to be a .arc or .xyz file!", filenames.get(0)));
        }

        /**File saveDir = baseDir
        *String modelFilename = assemblyState.mola.getFile().getAbsolutePath()
        *if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        *}
        *String dirName = saveDir.toString() + File.separator
        *String fileName = FilenameUtils.getName(modelFilename)
        *\\fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        *String arcFileName = fileName + ".arc"
        **/


        if (numSnaps > maxnum ) {
            logger.info(String.format(" Warning!!! System does not appear to contain enough entries! All %d energies will be reported", maxnum))
            numSnaps = maxnum

        }

        for (int i = 0; i < numSnaps-1; i++) {
            StateContainer savedState = energyQueue.removeLast()
            AssemblyState finalAssembly = savedState.getState()
            finalAssembly.revertState();
            double finalEnergy = savedState.getEnergy()
            logger.info(String.format("The potential energy found is %12.6g kcal/mol", finalEnergy))



            File saveDir = baseDir
            String modelFilename = assemblyState.mola.getFile().getAbsolutePath()
            if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                saveDir = new File(FilenameUtils.getFullPath(modelFilename))
            }
            String dirName = saveDir.toString() + File.separator
            String fileName = FilenameUtils.getName(modelFilename)
            //fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
            String arcFileName = fileName + ".arc"
            File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
            potentialFunctions.saveAsPDB(assemblyState.mola, saveFile)

        }


        StateContainer savedState = energyQueue.removeLast()
        AssemblyState lowestAssembly = savedState.getState()
        lowestEnergy = savedState.getEnergy()

        assemblyState.revertState()
        logger.info(String.format(" The lowest potential energy found is %12.6g kcal/mol", lowestEnergy))
        //prints our final energy (which will be the lowest energy

        File saveDir = baseDir
        String modelFilename = assemblyState.mola.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        //fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        String arcFileName = fileName + ".arc"
        File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
        potentialFunctions.saveAsPDB(assemblyState.mola, saveFile)

        return this
    }

    /**
     * Returns the lowest energy found.
     * @return Lowest potential energy.
     */
    public double getLowestEnergy() {
        return energy;
    }

    /**
     * Returns a copy of the lowest-energy state found.
     * @return Lowest-energy state found.
     */
    public AssemblyState getOptimumState() {
        return AssemblyState.copyState(assemblyState);
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
*/
