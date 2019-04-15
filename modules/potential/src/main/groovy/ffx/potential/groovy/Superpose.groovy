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
     * --atoms defines which atoms to calculate RMSD on.
     */
    @Option(names = ['--rA', '--rmsdAtoms'], paramLabel = "HEAVY",
            description = 'Atoms to be included in RMSD calculation. Select [ALL/HEAVY/ALPHA] to choose which atoms are used for RMSD (nucleic acids will use N1 or N9 in place of alpha carbons).')
    private String atomSelection = "HEAVY"

    /**
     * --atoms defines which atoms to calculate RMSD on.
     */
    @Option(names = ['--fC', '--frameComparison'], paramLabel = "true",
            description = 'Frames to be compared within the arc file. Select [true] for all versus all comparison; select [false] for one versus all comparison.')
    private boolean frameComparison = true

    /**
     * -s or --start defines which atoms in the structure will be used in RMSD calculation.
     */
    @Option(names = ['-s', '--start'], paramLabel = "-1",
            description = 'Atom number where RMSD calculation of structure will begin.')
    private int start = -1

    /**
     * --fi or --finish defines which atoms in the structure will be used in RMSD calculation.
     */
    @Option(names = ['-f', '--fi'], paramLabel = "Max_Value",
            description = 'Atom number where RMSD calculation of structure will end.')
    private int finish = Integer.MAX_VALUE

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

        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (systemFilter instanceof PDBFilter || systemFilter instanceof XYZFilter) {
            double[] x2 = new double[nVars]
            double[] mass = new double[nVars / 3]

            int nAtoms = atoms.length
            for (int i = 0; i < nAtoms; i++) {
                mass[i] = atoms[i].getMass()
            }

            // Begin streaming the possible atom indices, filtering out inactive atoms.
            // TODO: Decide if we only want active atoms.
            IntStream atomIndexStream = IntStream.range(0, atoms.length).
                    filter({ int i -> return atoms[i].isActive() }).
                    filter({ int i -> return atoms[i].xyzIndex >= start && atoms[i].xyzIndex <= finish })
            if (start > -1 || finish < Integer.MAX_VALUE) {
                logger.info(String.format("Calculating RMSD on residues %d to %d.", start, finish))
            }

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
                    selectionType = "Reference Atoms (i.e. alpha carbons and N1/N9 for nucleic acids)";
                    break

                case "ALL":
                case "1":
                    selectionType = "All Atoms"
                    // Unmodified stream; we have just checked for active atoms.
                    break

                default:
                    logger.severe(String.format(" Could not parse %s as an atom selection! Must be ALL, HEAVY, or ALPHA", atomSelection))
                    break
            }

            // Indices of atoms used in alignment and RMSD calculations.
            int[] usedIndices = atomIndexStream.toArray()
            System.out.println("used indices length: " + usedIndices.length)
            int nUsed = usedIndices.length
            int nUsedVars = nUsed * 3
            double[] massUsed = Arrays.stream(usedIndices).
                    mapToDouble({ int i -> atoms[i].getAtomType().atomicWeight }).
                    toArray()
            double[] xUsed = new double[nUsedVars]
            double[] x2Used = new double[nUsedVars]

            // Check which molecular assemblies to do RMSD comparisons among.
            if (frameComparison == false) {
                // The first snapshot is being used for all comparisons here; therefore, snapshot = 1.
                rmsd(systemFilter, nUsed, usedIndices, selectionType, x, x2, xUsed, x2Used, massUsed, 1)
            } else if (frameComparison == true) {
                rmsd(systemFilter, nUsed, usedIndices, selectionType, x, x2, xUsed, x2Used, massUsed, 1)
                SystemFilter systemFilter1
                if (systemFilter instanceof PDBFilter) {
                    systemFilter1 = new PDBFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                    systemFilter1.readFile()
                } else if (systemFilter instanceof XYZFilter) {
                    systemFilter1 = new XYZFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                }
                while (systemFilter1.readNext(false, false)) {
                    int snapshot1 = systemFilter1.getSnapshot()
                    forceFieldEnergy.getCoordinates(x)
                    SystemFilter systemFilter2
                    if (systemFilter instanceof PDBFilter) {
                        systemFilter2 = new PDBFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                        systemFilter2.readFile()
                    } else if (systemFilter instanceof XYZFilter) {
                        systemFilter2 = new XYZFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
                    }
                    rmsd(systemFilter2, nUsed, usedIndices, selectionType, x, x2, xUsed, x2Used, massUsed, snapshot1)
                }
            } else{
                logger.severe(String.format(" Could not parse %s as a frame comparison! Must be true or false", frameComparison))
            }
        }
        return this
    }

    void rmsd(SystemFilter systemFilter, int nUsed, int[] usedIndices, String selectionType, double[] x, double[] x2, double[] xUsed, double[] x2Used, double[] massUsed, int snapshot1) {
        while (systemFilter.readNext(false, false)) {
            int snapshot2 = systemFilter.getSnapshot()
            // Only calculate RMSD for snapshots if they aren't the same snapshot.
            // Also avoid double calculating snapshots in the matrix by only calculating the upper triangle.
            if (snapshot1 != snapshot2 && snapshot1 < snapshot2) {
                forceFieldEnergy.getCoordinates(x2)
                for (int i = 0; i < nUsed; i++) {
                    int index3 = 3 * usedIndices[i]
                    int i3 = 3 * i
                    for (int j = 0; j < 3; j++) {
                        xUsed[i3 + j] = x[index3 + j]
                        x2Used[i3 + j] = x2[index3 + j]
                    }
                }

                double origRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)
                ffx.potential.utils.Superpose.translate(xUsed, massUsed, x2Used, massUsed)
                double translatedRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)
                ffx.potential.utils.Superpose.rotate(xUsed, x2Used, massUsed)
                double rotatedRMSD = ffx.potential.utils.Superpose.rmsd(xUsed, x2Used, massUsed)

                logger.info(format(
                        "\n Coordinate RMSD Based On %s (Angstroms) on Model %d and Model %d\n Original:\t\t%7.3f\n After Translation:\t%7.3f\n After Rotation:\t%7.3f\n",
                        selectionType, snapshot1, snapshot2, origRMSD, translatedRMSD, rotatedRMSD))

            }
        }
    }
}
