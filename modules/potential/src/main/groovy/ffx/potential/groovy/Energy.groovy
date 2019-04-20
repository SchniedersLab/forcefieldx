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
package ffx.potential.groovy

import java.util.logging.Level
import static java.lang.String.format

import com.google.common.collect.MinMaxPriorityQueue

import org.apache.commons.io.FilenameUtils

import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc Energy")
class Energy extends PotentialScript {

    /**
     * -g or --gradient to print out gradients.
     */
    @Option(names = ['-g', '--gradient'], paramLabel = "false",
            description = 'Print out atomic gradients as well as energy.')
    private boolean gradient = false

    /**
     * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
     */
    @Option(names = ['--es1', '--noElecStart1'], paramLabel = "1",
            description = 'Starting no-electrostatics atom for 1st topology')
    private int es1 = 1
    /**
     * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
     */

    @Option(names = ['--fl', '--findLowest'], paramLabel = "0",
            description = 'Return the n lowest energy structures from an ARC or PDB file.')
    private int fl = 0

    /**
     * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
     */
    @Option(names = ['--ef1', '--noElecFinal1'], paramLabel = "-1",
            description = 'Final no-electrostatics atom for 1st topology')
    private int ef1 = -1

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private List<String> filenames = null


    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null


    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

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
    Energy run() {

        if (!init()) {
            return this
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info(" Running Energy on " + filename)

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        Atom[] atoms = activeAssembly.getAtomArray()

        // Apply the no electrostatics atom selection
        int noElecStart = es1
        noElecStart = (noElecStart < 1) ? 1 : noElecStart

        int noElecStop = ef1
        noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop

        if (noElecStart <= noElecStop) {
            logger.info(format(" Disabling electrostatics for atoms %d (%s) to %d (%s).",
                    noElecStart, atoms[noElecStart - 1], noElecStop, atoms[noElecStop - 1]))
        }
        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1]
            ai.setElectrostatics(false)
            ai.print(Level.FINE)
        }

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        if (gradient) {
            double[] g = new double[nVars]
            int nAts = nVars / 3
            energy = forceFieldEnergy.energyAndGradient(x, g, true)
            logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
            for (int i = 0; i < nAts; i++) {
                int i3 = 3 * i
                logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
            }
        } else {
            energy = forceFieldEnergy.energy(x, true)
        }

        SystemFilter systemFilter = potentialFunctions.getFilter()

        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {

            int numSnaps = fl
            double lowestEnergy = energy
            assemblyState = new AssemblyState(activeAssembly)
            int index = 1

            // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
            MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = null
            if (fl > 0) {
                lowestEnergyQueue = MinMaxPriorityQueue
                        .maximumSize(numSnaps)
                        .create()
                lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy))
            }

            while (systemFilter.readNext()) {
                index++
                forceFieldEnergy.getCoordinates(x)
                energy = forceFieldEnergy.energy(x, false)
                logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))

                if (fl > 0) {
                    lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), energy))
                }
            }

            if (fl > 0) {
                if (numSnaps > index) {
                    logger.warning(format(
                            " Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported",
                            numSnaps, filename, index, index))
                    numSnaps = index
                }

                File saveDir = baseDir
                String modelFilename = activeAssembly.getFile().getAbsolutePath()
                if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                    saveDir = new File(FilenameUtils.getFullPath(modelFilename))
                }
                String dirName = saveDir.toString() + File.separator
                String fileName = FilenameUtils.getName(modelFilename)

                for (int i = 0; i < numSnaps - 1; i++) {
                    StateContainer savedState = lowestEnergyQueue.removeLast()
                    AssemblyState finalAssembly = savedState.getState()
                    finalAssembly.revertState()
                    double finalEnergy = savedState.getEnergy()
                    logger.info(format(" The potential energy found is %16.8f (kcal/mol)", finalEnergy))
                    File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                    MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                    potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
                }

                StateContainer savedState = lowestEnergyQueue.removeLast()
                AssemblyState lowestAssembly = savedState.getState()
                lowestEnergy = savedState.getEnergy()

                assemblyState.revertState()
                logger.info(format(" The lowest potential energy found is %16.8f (kcal/mol)", lowestEnergy))

                // Prints our final energy (which will be the lowest energy
                File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
            }
        }

        return this

    }

    @Override
    List<Potential> getPotentials() {
        return forceFieldEnergy == null ? Collections.emptyList() : Collections.singletonList(forceFieldEnergy)
    }
}

