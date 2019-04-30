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

import java.util.stream.IntStream
import static java.lang.String.format

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
 * The Superpose script superposes all molecules in an arc file to the first molecule in the arc file and reports the RMSD.
 * <br>
 * Usage:
 * <br>
 * ffxc Superpose [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "ffxc SaveAsPDB")
class Superpose extends PotentialScript {

    /**
     * --aS or --atomSelection The atom selection [HEAVY (0) / ALL (1) / CALPHA (2)] for the RMSD calculation (CALPHA chooses N1 or N9 for nucleic acids).
     */
    @Option(names = ['--aS', '--atomSelection'], paramLabel = "0",
            description = 'The atom selection [HEAVY (0) / ALL (1) / CALPHA (2)] for the RMSD calculation (CALPHA chooses N1 or N9 for nucleic acids).')
    private String atomSelection = "0"

    /**
     * -A or --allvsAll Frames to be compared within the arc file. Select [true] for all versus all comparison; select [false] for one versus all comparison.
     */
    @Option(names = ['-A', '--allvsAll'], paramLabel = "false",
            description = 'Compare all snapshots versus all others, instead of the first snapshot versus all others.')
    private boolean frameComparison = false

    /**
     * -s or --start Atom number where RMSD calculation of structure will begin.
     */
    @Option(names = ['-s', '--start'], paramLabel = "1",
            description = 'Starting atom to include in the RMSD calculation.')
    private int start = 1

    /**
     * -f or --final Atom number where RMSD calculation of structure will end.
     */
    @Option(names = ['-f', '--final'], paramLabel = "nAtoms",
            description = 'Final atom to include in the RMSD calculation.')
    private int finish = Integer.MAX_VALUE

    /**
     * -w or --write Write out superposed snapshots.
     */
    @Option(names = ['-w', '--write'], paramLabel = "false",
            description = 'Write out superposed snapshots.')
    private boolean writeSnapshots = false

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

    public ForceFieldEnergy forceFieldEnergy = null
    private File outFile
    private XYZFilter outputFilter

    /**
     * Execute the script.
     */
    @Override
    Superpose run() {
        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        Atom[] atoms = activeAssembly.getAtomArray()

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        if (writeSnapshots) {
            String outFileName = activeAssembly.getFile().toString().replaceFirst(~/\.(?:xyz|pdb|arc).*$/, "")
            outFileName = outFileName + "_superposed.arc"
            outFileName = potentialFunctions.versionFile(outFileName)

            outFile = new File(outFileName)
            outputFilter = new XYZFilter(outFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
        }

        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (systemFilter instanceof PDBFilter || systemFilter instanceof XYZFilter) {
            double[] x2 = new double[nVars]
            double[] mass = new double[nVars / 3]

            int nAtoms = atoms.length
            for (int i = 0; i < nAtoms; i++) {
                mass[i] = atoms[i].getMass()
            }

            if (finish > nAtoms - 1) {
                finish = nAtoms - 1
            }
            if (start < 0 || start > finish) {
                start = 0
            }

            // Note that atoms are indexed from 0 to nAtoms - 1.
            logger.info(format(" Atoms from %d to %d will be considered.", start + 1, finish + 1))

            // Begin streaming the possible atom indices, filtering out inactive atoms.
            IntStream atomIndexStream = IntStream.range(start, finish + 1).filter({ int i -> return atoms[i].isActive() })

            // String describing the selection type.
            String selectionType = "All Atoms"

            // Switch on what type of atoms to select, filtering as appropriate. Support the old integer indices.
            switch (atomSelection.toUpperCase()) {
                case "HEAVY":
                case "0":
                    // Filter only for heavy (non-hydrogen) atoms.
                    atomIndexStream = atomIndexStream.filter({ int i -> atoms[i].isHeavy() })
                    selectionType = "Heavy Atoms"
                    break
                case "ALL":
                case "1":
                    // Unmodified stream; we have just checked for active atoms.
                    selectionType = "All Atoms"
                    break
                case "ALPHA":
                case "2":
                    // Filter only for reference atoms: carbons named CA (protein) or nitrogens named N1 or N9 (nucleic acids).
                    atomIndexStream = atomIndexStream.filter({ int i ->
                        Atom ati = atoms[i]
                        String atName = ati.getName().toUpperCase()
                        boolean proteinReference = atName.equals("CA") && ati.getAtomType().atomicNumber == 6
                        boolean naReference = (atName.equals("N1") || atName.equals("N9")) && ati.getAtomType().atomicNumber == 7
                        return proteinReference || naReference
                    })
                    selectionType = "C-Alpha Atoms (or N1/N9 for nucleic acids)"
                    break
                default:
                    logger.severe(format(" Could not parse %s as an atom selection! Must be ALL, HEAVY, or ALPHA", atomSelection))
                    break
            }

            logger.info(" Superpose selection criteria: " + selectionType)

            // Indices of atoms used in alignment and RMSD calculations.
            int[] usedIndices = atomIndexStream.toArray()
            int nUsed = usedIndices.length
            int nUsedVars = nUsed * 3
            double[] massUsed = Arrays.stream(usedIndices).
                    mapToDouble({ int i -> atoms[i].getAtomType().atomicWeight }).toArray()
            double[] xUsed = new double[nUsedVars]
            double[] x2Used = new double[nUsedVars]

            if (writeSnapshots) {
                AssemblyState origState = new AssemblyState(activeAssembly)
                forceFieldEnergy.getCoordinates(x2)
                copyCoordinates(nUsed, usedIndices, x2, x2Used)
                double[] translate = ffx.potential.utils.Superpose.calculateTranslation(x2Used, massUsed)
                ffx.potential.utils.Superpose.applyTranslation(x2, translate)
                forceFieldEnergy.setCoordinates(x2)
                outputFilter.writeFile(outFile, true)
                origState.revertState()
            }

            // Check which molecular assemblies to do RMSD comparisons among.
            if (!frameComparison) {
                // The first snapshot is being used for all comparisons here; therefore, snapshot = 1.
                rmsd(systemFilter, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, 1)
            } else {
                rmsd(systemFilter, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, 1)
                SystemFilter systemFilter1 = null
                if (systemFilter instanceof PDBFilter) {
                    systemFilter1 = new PDBFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                    systemFilter1.readFile()
                } else if (systemFilter instanceof XYZFilter) {
                    systemFilter1 = new XYZFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                }
                while (systemFilter1.readNext(false, false)) {
                    int snapshot1 = systemFilter1.getSnapshot()
                    forceFieldEnergy.getCoordinates(x)
                    SystemFilter systemFilter2 = null
                    if (systemFilter instanceof PDBFilter) {
                        systemFilter2 = new PDBFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                        systemFilter2.readFile()
                    } else if (systemFilter instanceof XYZFilter) {
                        systemFilter2 = new XYZFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                    }
                    rmsd(systemFilter2, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, snapshot1)
                }
            }
        }
        return this
    }

    /**
     * Copy coordinates from the entire system to the used subset.
     *
     * @param nUsed Number of atoms used.
     * @param usedIndices Mapping from the xUsed array to its source in x.
     * @param x All atomic coordinates.
     * @param xUsed The used subset of coordinates.
     */
    private static void copyCoordinates(int nUsed, int[] usedIndices, double[] x, double[] xUsed) {
        for (int i = 0; i < nUsed; i++) {
            int index3 = 3 * usedIndices[i]
            int i3 = 3 * i
            for (int j = 0; j < 3; j++) {
                xUsed[i3 + j] = x[index3 + j]
            }
        }
    }

    void rmsd(SystemFilter systemFilter, int nUsed, int[] usedIndices, double[] x, double[] x2, double[] xUsed, double[] x2Used, double[] massUsed, int snapshot1) {
        double[] xBak = Arrays.copyOf(x, x.length)
        while (systemFilter.readNext(false, false)) {
            int snapshot2 = systemFilter.getSnapshot()
            // Only calculate RMSD for snapshots if they aren't the same snapshot.
            // Also avoid double calculating snapshots in the matrix by only calculating the upper triangle.
            if (snapshot1 != snapshot2 && snapshot1 < snapshot2) {
                AssemblyState origStateB = new AssemblyState(activeAssembly)
                forceFieldEnergy.getCoordinates(x2)
                copyCoordinates(nUsed, usedIndices, x, xUsed)
                copyCoordinates(nUsed, usedIndices, x2, x2Used)

                double origRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)

                // Calculate the translation on only the used subset, but apply it to the entire structure.
                double[] tA = ffx.potential.utils.Superpose.calculateTranslation(xUsed, massUsed)
                ffx.potential.utils.Superpose.applyTranslation(x, tA)
                double[] tB = ffx.potential.utils.Superpose.calculateTranslation(x2Used, massUsed)
                ffx.potential.utils.Superpose.applyTranslation(x2, tB)
                // Copy the applied translation to xUsed and x2Used.
                copyCoordinates(nUsed, usedIndices, x, xUsed)
                copyCoordinates(nUsed, usedIndices, x2, x2Used)
                double translatedRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)

                // Calculate the rotation on only the used subset, but apply it to the entire structure.
                double[][] rotation = ffx.potential.utils.Superpose.calculateRotation(xUsed, x2Used, massUsed)
                ffx.potential.utils.Superpose.applyRotation(x2, rotation)
                // Copy the applied rotation to x2Used.
                copyCoordinates(nUsed, usedIndices, x2, x2Used)
                double rotatedRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)

                logger.info(format(
                        " Coordinate RMSD for %d and %d: Original %7.3f, After Translation %7.3f, After Rotation %7.3f",
                        snapshot1, snapshot2, origRMSD, translatedRMSD, rotatedRMSD))

                if (writeSnapshots) {
                    forceFieldEnergy.setCoordinates(x2)
                    outputFilter.writeFile(outFile, true)
                    origStateB.revertState()
                }
                System.arraycopy(xBak, 0, x, 0, x.length)
            }
        }
    }
}
