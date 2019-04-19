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

import com.google.common.collect.MinMaxPriorityQueue
import ffx.potential.AssemblyState
import org.apache.commons.io.FilenameUtils

import java.util.logging.Level
import static java.lang.String.format

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.utils.Superpose

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
     * * --fl or --findlowest returns the n lowest energy structures.
     */

    @Option(names = ['--fl', '--findlowest'], paramLabel = "0",
            description = 'Find this number of lowest-energy structures in a .arc file (does not function with .pdb)')
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
    Energy run() {

        if (!init()) {
            return
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
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
        if (systemFilter instanceof XYZFilter) {
            XYZFilter xyzFilter = (XYZFilter) systemFilter

            double[] x2 = new double[nVars]
            double[] mass = new double[nVars / 3]

            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                mass[i] = atoms[i].getMass()
            }

            //Get heavy atom masses.
            int nHeavyVars = forceFieldEnergy.getNumberOfHeavyAtomVariables()
            double[] massHeavy = new double[nHeavyVars / 3]
            for (int i = 0; i < nHeavyVars / 3; i++) {
                if (!atoms[i].isHydrogen()) {
                    massHeavy[i] = atoms[i].getMass()
                }
            }

            //Array containing heavy atom indices.
            int[] heavyAtomPositions = new int[nHeavyVars / 3];
            int j = 0;
            for (int i = 0; i < nVars / 3; i++) {
                if (!atoms[i].isHydrogen()) {
                    heavyAtomPositions[j] = i
                    j++
                }
            }

            int numSnaps = fl
            /*double lowestEnergy = Double.MAX_VALUE
            lowestEnergy = forceFieldEnergy.energy(x, false)*/
            double lowestEnergy = energy;
            assemblyState = new AssemblyState(activeAssembly)

            int maxnum = 1

            /**
             * Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
             */

            MinMaxPriorityQueue<StateContainer> lowestEnergyQueue;
            if (fl > 0) {
                lowestEnergyQueue = MinMaxPriorityQueue
                        .maximumSize(numSnaps)
                        .create()
                StateContainer firstContainer;
                firstContainer = new StateContainer(assemblyState, lowestEnergy)
                lowestEnergyQueue.add(firstContainer)
            }

            while (xyzFilter.readNext()) {
                //Arrays for holding coordinates of heavy atoms after rotation and translation.
                double[] xHeavy = new double[nHeavyVars]
                double[] x2Heavy = new double[nHeavyVars]

                forceFieldEnergy.getCoordinates(x2)
                energy = forceFieldEnergy.energy(x2, true)

                //Original RMSD.
                for (int i = 0; i < nHeavyVars / 3; i++) {
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i * 3] = x[positionOfHeavyAtom]
                    xHeavy[i * 3 + 1] = x[positionOfHeavyAtom + 1]
                    xHeavy[i * 3 + 2] = x[positionOfHeavyAtom + 2]
                    x2Heavy[i * 3] = x2[positionOfHeavyAtom]
                    x2Heavy[i * 3 + 1] = x2[positionOfHeavyAtom + 1]
                    x2Heavy[i * 3 + 2] = x2[positionOfHeavyAtom + 2]
                }
                double origRMSDHeavy = Superpose.rmsd(xHeavy, x2Heavy, massHeavy);

                //Translated RMSD.
                Superpose.translate(x, mass, x2, mass)
                for (int i = 0; i < nHeavyVars / 3; i++) {
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i * 3] = x[positionOfHeavyAtom]
                    xHeavy[i * 3 + 1] = x[positionOfHeavyAtom + 1]
                    xHeavy[i * 3 + 2] = x[positionOfHeavyAtom + 2]
                    x2Heavy[i * 3] = x2[positionOfHeavyAtom]
                    x2Heavy[i * 3 + 1] = x2[positionOfHeavyAtom + 1]
                    x2Heavy[i * 3 + 2] = x2[positionOfHeavyAtom + 2]
                }
                double transRMSDHeavy = Superpose.rmsd(xHeavy, x2Heavy, massHeavy)

                //Rotated RMSD.
                Superpose.rotate(x, x2, mass)
                for (int i = 0; i < nHeavyVars / 3; i++) {
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i * 3] = x[positionOfHeavyAtom]
                    xHeavy[i * 3 + 1] = x[positionOfHeavyAtom + 1]
                    xHeavy[i * 3 + 2] = x[positionOfHeavyAtom + 2]
                    x2Heavy[i * 3] = x2[positionOfHeavyAtom]
                    x2Heavy[i * 3 + 1] = x2[positionOfHeavyAtom + 1]
                    x2Heavy[i * 3 + 2] = x2[positionOfHeavyAtom + 2]
                }
                double rotRMSDHeavy = Superpose.rmsd(xHeavy, x2Heavy, massHeavy)
                logger.info(format(
                        "\n Coordinate RMSD Based On Heavy Atoms (Angstroms)\n Original:\t\t%7.3f\n After Translation:\t%7.3f\n After Rotation:\t%7.3f\n",
                        origRMSDHeavy, transRMSDHeavy, rotRMSDHeavy))

                if (fl > 0) {

                    forceFieldEnergy.getCoordinates(x) // getting the coordinates for the next assembly
                    lowestEnergy = forceFieldEnergy.energy(x, false) //calculating energy for new assembly
                    assemblyState = new AssemblyState(activeAssembly)
                    lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy))
                    ++maxnum;
                }

            }

            if (fl > 0) {

                if (numSnaps > maxnum) {
                    logger.warning(String.format(" Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported", numSnaps, filename, maxnum, maxnum))
                    numSnaps = maxnum
                }

                for (int i = 0; i < numSnaps - 1; i++) {
                    StateContainer savedState = lowestEnergyQueue.removeLast()
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
                    String arcFileName = fileName + ".pdb"
                    File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                    potentialFunctions.saveAsPDB(assemblyState.mola, saveFile)

                }


                StateContainer savedState = lowestEnergyQueue.removeLast()
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
                File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                potentialFunctions.saveAsPDB(assemblyState.mola, saveFile)

            }
        }
        return this
        
    }

    @Override
    List<Potential> getPotentials() {
        return forceFieldEnergy == null ? Collections.emptyList() : Collections.singletonList(forceFieldEnergy);
    }
}

