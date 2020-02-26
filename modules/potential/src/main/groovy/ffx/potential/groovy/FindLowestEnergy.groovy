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
package ffx.potential.groovy

import com.google.common.collect.MinMaxPriorityQueue

import org.apache.commons.io.FilenameUtils

import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The FindLowestEnergy script calculates energies for all assemblies in an arc file,
 * finds the lowest energy assembly, and saves that assembly to a pdb file.
 * <br>
 * Usage:
 * <br>
 * ffxc FindLowestEnergy [options] &lt;filename&gt;
 */
@Command(description = " Finds the lowest energy structures in an arc file.", name = "ffxc FindLowestEnergy")
class FindLowestEnergy extends PotentialScript {

    /**
     * -K or --nLowest Finds the K lowest energy structures in an arc file.
     */
    @Option(names = ['-K', '--nLowest'], paramLabel = "1", defaultValue = "1",
            description = 'Finds the K lowest energy structures in an arc file.')
    int numSnaps = 1

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    private double energy = Double.MAX_VALUE
    private AssemblyState assemblyState = null

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

    /**
     * Execute the script.
     */
    @Override
    FindLowestEnergy run() {
        if (!init()) {
            return null
        }

        if (filenames == null || filenames.size() != 1) {
            logger.info(helpString())
            return null
        } else {
            MolecularAssembly[] assemblies = [potentialFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
        }

        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)
        energy = forceFieldEnergy.energy(x, true)
        assemblyState = new AssemblyState(activeAssembly)
        SystemFilter systemFilter = potentialFunctions.getFilter()

        // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
        if (numSnaps < 1) {
            numSnaps = 1
        }

        MinMaxPriorityQueue<StateContainer> energyQueue
        energyQueue = MinMaxPriorityQueue.maximumSize(numSnaps).create()
        StateContainer firstContainer
        firstContainer = new StateContainer(assemblyState, energy)
        energyQueue.add(firstContainer)

        int maxNumber = 1
        if (systemFilter instanceof XYZFilter) {
            XYZFilter xyzFilter = (XYZFilter) systemFilter
            while (xyzFilter.readNext()) {
                forceFieldEnergy.getCoordinates(x)
                energy = forceFieldEnergy.energy(x, true)
                assemblyState = new AssemblyState(activeAssembly)
                // Save the new assembly if the energy is less than the current energy
                energyQueue.add(new StateContainer(assemblyState, energy))
                maxNumber = maxNumber + 1
            }
        } else {
            logger.severe(String.format(" System %s does not appear to be a .arc or .xyz file!", filenames.get(0)))
        }

        if (numSnaps > maxNumber) {
            logger.info(String.format(" The archive does not contain enough entries; all %d energies will be reported.", maxNumber))
            numSnaps = maxNumber
        }

        for (int i = 0; i < numSnaps - 1; i++) {
            StateContainer savedState = energyQueue.removeLast()
            AssemblyState finalAssembly = savedState.getState()
            MolecularAssembly molecularAssembly = finalAssembly.getMolecularAssembly()
            finalAssembly.revertState()
            double finalEnergy = savedState.getEnergy()
            logger.info(String.format("The potential energy found is %12.6g kcal/mol", finalEnergy))

            File saveDir = baseDir
            String modelFilename = molecularAssembly.getFile().getAbsolutePath()
            if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                saveDir = new File(FilenameUtils.getFullPath(modelFilename))
            }
            String dirName = saveDir.toString() + File.separator
            String fileName = FilenameUtils.getName(modelFilename)
            File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
            potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
        }

        StateContainer savedState = energyQueue.removeLast()
        double lowestEnergy = savedState.getEnergy()
        AssemblyState finalAssembly = savedState.getState()
        MolecularAssembly molecularAssembly = finalAssembly.getMolecularAssembly()
        finalAssembly.revertState()
        logger.info(String.format(" The lowest potential energy found is %12.6g kcal/mol", lowestEnergy))

        File saveDir = baseDir
        String modelFilename = molecularAssembly.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
        potentialFunctions.saveAsPDB(molecularAssembly, saveFile)

        return this
    }

    /**
     * Returns the lowest energy found.
     * @return Lowest potential energy.
     */
    double getLowestEnergy() {
        return energy
    }

    /**
     * Returns a copy of the lowest-energy state found.
     * @return Lowest-energy state found.
     */
    AssemblyState getOptimumState() {
        return AssemblyState.copyState(assemblyState)
    }
}
