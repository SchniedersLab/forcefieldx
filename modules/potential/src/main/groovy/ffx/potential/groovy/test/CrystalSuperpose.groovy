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
package ffx.potential.groovy.test

import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.ProgressiveAlignmentOfCrystals
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * Quantifies the similarity of input crystals based on progressive alignments.
 * This script is based off of the PACCOM code created by Okimasa Okada.
 *
 * @author Okimasa OKADA
 * @author Aaron J. Nessler and Michael J. Schnieders
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler and Micheal Schnieders 2020
 * revised by Aaron Nessler and Michael Schnieders 2021
 * <br>
 * Usage:
 * <br>
 * ffxc test.CrystalSuperpose &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Compare crystal polymorphs based on intermolecular distances of aligned crystals.", name = "ffxc test.CrystalSuperpose")
class CrystalSuperpose extends PotentialScript {

    /**
     * --nm or --numMolecules Number of asymmetric units to include from each crystal in RMSD comparison.
     */
    @Option(names = ['--na', '--numAU'], paramLabel = '20',
            description = 'Determines number of asymmetric units to include in final comparison.')
    int nAU = 20

    /**
     * --im or --inflatedAU Number of molecules in the inflated sphere.
     */
    @Option(names = ['--ia', '--inflatedAU'], paramLabel = '100',
            description = 'Specifies the minimum number of asymmetric units in the inflated crystal.')
    int inflatedAU = 200

    /**
     * --ns or --numSearch Number of molecules to loop through in first system.
     */
    @Option(names = ['--ns', '--numSearch'], paramLabel = '-1',
            description = 'Determines number of asymmetric units to search in first system (mirror check).')
    int numSearch = 5

    /**
     * --ns2 or --numSearch2 Number of molecules to loop through in second system.
     */
    @Option(names = ['--ns2', '--numSearch2'], paramLabel = '-1',
            description = 'Determines number of asymmetric units to search in second system (mirror check).')
    int numSearch2 = 5

    /**
     * --nms or --noMirrorSearch Internal loop over (--ns) molecules to check for mirrors.
     */
    @Option(names = ['--nms', '--noMirrorSearch'], paramLabel = "false", defaultValue = "false",
            description = 'Loop over structures looking for mirrors.')
    private static boolean noMirrorSearch = false

    /**
     * -w or --write Write out a distance matrix for each comparison.
     */
    @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
            description = 'Write a distance comparison matrix to a text file.')
    private static boolean write = false

    /**
     * --sp or --savePDB Save out a PDB.
     */
    @Option(names = ['--sp', '--savePDB'], paramLabel = "false", defaultValue = "false",
            description = 'Save out a PDB file for the comparison.')
    private static boolean save = false

    /**
     * -f or --force Perform comparison even if a single asymmetric unit from each crystal differs
     */
    @Option(names = ['-f', '--force'], paramLabel = "false", defaultValue = "false",
            description = 'Force comparison regardless of single asymmetric unit RMSD.')
    private static boolean force = false

    /**
     * --nh or --noHydrogens Perform comparison without hydrogen atoms.
     */
    @Option(names = ['--nh', '--noHydrogens'], paramLabel = "false", defaultValue = "false",
            description = 'Crystal RMSD calculated without hydrogen atoms.')
    private static boolean noHydrogens = false

    //TODO finish implementing symmetric flag.
    /**
     * --sym or --symmetric Attempt prioritization in a symmetric manner (should produce symmetric distance matrix).
     */
    @Option(names = ['--sym', '--symmetric'], paramLabel = "false", defaultValue = "false",
            description = 'Attempt to enforce symmetric output matrix.')
    private static boolean symmetric = false

    /**
     * Select individual atoms to calculate RMSD rather than using all atoms.
     */
    @Option(names = ['--da', '--desiredAtoms'], arity = "1..*", paramLabel = "atom integers",
            description = 'Specify atoms to calculate RMSD. Otherwise, all atom.')
    private int[] desiredAtoms = null

    /**
     * CrystalSuperpose Test requires a public variable containing observables to test.
     */
    public double[][] distMatrix

    /**
     * The final argument(s) should be two or more filenames (same file twice if comparing same structures).
     */
    @Parameters(arity = "1..2", paramLabel = "files",
            description = 'Atomic coordinate files to compare in XYZ format.')
    List<String> filenames = null

    /**
     * CrystalSuperpose Constructor.
     */
    CrystalSuperpose() {
        this(new Binding())
    }

    /**
     * CrystalSuperpose Constructor.
     * @param binding Groovy Binding to use.
     */
    CrystalSuperpose(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    @Override
    CrystalSuperpose run() {
        // Turn off non-bonded interactions.
        System.setProperty("vdwterm", "false")

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // Ensure file exists.
        if (filenames == null) {
            logger.info(helpString())
            return this
        }

        // Number of files to read in.
        int numFiles = filenames.size()

        // System Filter containing structures stored in file 0.
        SystemFilter systemFilter
        // System Filter containing structures stored in file 1.
        SystemFilter systemFilter2

        // If only one file compare structures within to self. Otherwise compare structures between files 0 and 1.
        if(numFiles==1){
            potentialFunctions.openAll(filenames.get(0))
            systemFilter = potentialFunctions.getFilter()
            //Need to reinitialize a new filter object
            potentialFunctions.openAll(filenames.get(0))
            systemFilter2 = potentialFunctions.getFilter()
        }else{
            potentialFunctions.openAll(filenames.get(0))
            systemFilter = potentialFunctions.getFilter()
            potentialFunctions.openAll(filenames.get(1))
            systemFilter2 = potentialFunctions.getFilter()
        }

        // Example atoms to determine single molecule characteristics (e.g. number of atoms, hydrogens, etc.)
        final List<Atom> exampleAtoms = systemFilter.getActiveMolecularSystem().getAtomList()
        // Number of atoms in asymmetric unit.
        int nAtoms = exampleAtoms.size()

        // If comparison atoms are not specified, use all atoms.
        Integer[] comparisonAtoms
        if (desiredAtoms != null) {
            comparisonAtoms = new int[desiredAtoms.length]
            for (int i = 0; i < desiredAtoms.size(); i++) {
                comparisonAtoms[i] = desiredAtoms[i] - 1
                if (comparisonAtoms[i] < 0 || comparisonAtoms[i] >= nAtoms) {
                    logger.severe(" Selected atoms are outside of molecular size.")
                    return this
                }
            }
        } else {
            //TODO implement no hydrgoens && centerAtoms
            if (noHydrogens) {
                int n = 0
                for (int i = 0; i < nAtoms; i++) {
                    if (!exampleAtoms.get(i).isHydrogen()) {
                        n++
                    }
                }
                comparisonAtoms = new int[n]
                n = 0
                for (int i = 0; i < nAtoms; i++) {
                    if (!exampleAtoms.get(i).isHydrogen()) {
                        comparisonAtoms[n++] = i
                    }
                }
            } else {
                comparisonAtoms = new int[nAtoms]
                for (int i = 0; i < nAtoms; i++) {
                    comparisonAtoms[i] = i
                }
            }
        }
        // Number of atoms being included for comparison.
        int compareAtomsSize = comparisonAtoms.size()
        logger.info(format(" Number of atoms being compared: %3d of %3d\n",
                compareAtomsSize, nAtoms))
        // Current method to save PDB only works when using all atoms...
        if(compareAtomsSize != nAtoms && save){
            save = false
            logger.warning(" Saving a PDB is not compatible with subsets of atoms (--nh and --da).")
        }

        // If search is desired ensure inner loop will be used. Else use single comparison.
        if(noMirrorSearch){
            numSearch = 1
            numSearch2 = 1
        }

        // Create object to perform comparison (I believe this was a necessary step to achieve parallelization...)
        ProgressiveAlignmentOfCrystals progressiveAlignmentOfCrystals = new ProgressiveAlignmentOfCrystals(systemFilter, systemFilter2)
        // Compare structures in SystemFilter and SystemFilter2.
        distMatrix = progressiveAlignmentOfCrystals.comparisons(nAtoms, Arrays.asList(comparisonAtoms), nAU, inflatedAU,
                numSearch, numSearch2, write, force, symmetric, save)

        return this
    }
}