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

import static java.lang.String.format

import org.apache.commons.io.FilenameUtils
import static org.apache.commons.math3.util.FastMath.pow

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.nonbonded.GaussVol

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The SaveAsP1 script expands a specified file to P1
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsP1 [options] &lt;filename&gt;
 */
@Command(description = " Calculate the surface area and volume of molecular system.", name = "ffxc Volume")
class Volume extends PotentialScript {

    private static final double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0)

    /**
     * -y or --includeHydrogen leaves in hydrogen when calculating the overlap tree.
     */
    @CommandLine.Option(names = ['-y', '--includeHydrogen'], paramLabel = "false",
            description = "Include Hydrogen in calculation of overlaps and volumes")
    private boolean includeHydrogen = false

    /**
     * -p or --probe Add a probe radius offset to all atomic radii.
     */
    @CommandLine.Option(names = ['-p', '--probe'], paramLabel = "0.0",
            description = "Add a probe radius offset to all atomic radii")
    private double probe = 0.0

    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @CommandLine.Option(names = ['-v', '--verbose'], paramLabel = "false",
            description = "Print out all components of volume of molecule and offset")
    private boolean verbose = false

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

    /**
     * JUnit Testing Variables
     */
    public double totalVolume = 0.0;
    public double totalSurfaceArea = 0.0;

    /**
     * Execute the script.
     */

    @Override
    Volume run() {
        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Calculating volume and surface area for " + modelFilename)

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        // Input
        boolean[] isHydrogen = new boolean[nAtoms]
        double[] radii = new double[nAtoms]
        double[] volume = new double[nAtoms]
        double[] gamma = new double[nAtoms]
        double[][] positions = new double[nAtoms][3]

        Arrays.fill(gamma, 1.0)
        double fourThirdsPI = 4.0 / 3.0 * Math.PI
        int index = 0
        for (Atom atom : atoms) {
            isHydrogen[index] = atom.isHydrogen()
            if (includeHydrogen) {
                isHydrogen[index] = false
            }
            radii[index] = atom.getVDWType().radius / 2.0// * rminToSigma
            radii[index] += probe
            volume[index] = fourThirdsPI * pow(radii[index], 3)
            positions[index][0] = atom.getX()
            positions[index][1] = atom.getY()
            positions[index][2] = atom.getZ()
            index++
        }

        // Run Volume calculation to get vdw volume of molecule
        GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen)
        gaussVol.computeVolumeAndSA(positions)
        logger.info(format("\n Maximum depth of overlaps in tree: %d", gaussVol.getMaximumDepth()))
        if(verbose){
            gaussVol.printTree()
        }
        logger.info(format("\n Total number of overlaps in tree: %d", gaussVol.getTotalNumberOfOverlaps()))

        // Calculate effective radius by assuming the GaussVol volume is the volume of a sphere
        double threeOverFourPi = 3.0/(4.0*Math.PI)
        double radical = gaussVol.getVolume()*threeOverFourPi
        double effectiveRadius = pow(radical, 1/3)

        logger.info(format("\n Volume:              %8.3f (Ang^3)", gaussVol.getVolume()))
        logger.info(format(" Volume Energy:       %8.3f (kcal/mol)", gaussVol.getVolumeEnergy()))
        logger.info(format(" Surface Area:        %8.3f (Ang^2)", gaussVol.getSurfaceArea()))
        logger.info(format(" Surface Area Energy: %8.3f (kcal/mol)", gaussVol.getSurfaceAreaEnergy()))
        logger.info(format(" Volume + SA Energy:  %8.3f (kcal/mol)", gaussVol.getEnergy()))
        logger.info(format(" Effective Radius:    %8.3f (Ang)", effectiveRadius))

        // Set JUnit testing variables based on output volume and surface area
        totalVolume = gaussVol.getVolume()
        totalSurfaceArea = gaussVol.getSurfaceArea()

        return this
    }
}
