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

import ffx.potential.bonded.Atom
import ffx.potential.nonbonded.GaussVol
import org.apache.commons.io.FilenameUtils

import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static org.apache.commons.math3.util.FastMath.pow

/**
 * The SaveAsP1 script expands a specified file to P1
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsP1 [options] &lt;filename&gt;
 */
@Command(description = " Calculate the surface area and volume of molecular system.", name = "ffxc Volume")
class Volume extends PotentialScript {

    private static final double angToNm =0.1
    private static final double nmToAng =10.0
    private static final double nmVolumeToAngVolume=1000
    private static final double kjToKcal=1.0/4.184
    private static final double offset=0.00005

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
        double[] radiiOffset = new double[nAtoms]
        double[] volume = new double[nAtoms]
        double[] volumeOffset = new double[nAtoms]
        double[] gamma = new double[nAtoms]
        double[][] positions = new double[nAtoms][3]

        // Output
        double[][] force = new double[nAtoms][3]
        double[] totalVolume = new double[1]
        double[] energy = new double[1]
        double[] gradV = new double[nAtoms]
        double selfVolumeSum = 0
        double[] freeVolume = new double[nAtoms]
        double[] selfVolume = new double[nAtoms]

        double[][] forceOffset = new double[nAtoms][3]
        double[] totalVolumeOffset = new double[1]
        double[] energyOffset = new double[1]
        double[] gradVOffset = new double[nAtoms]
        double selfVolumeOffsetSum = 0
        double[] freeVolumeOffset = new double[nAtoms]
        double[] selfVolumeOffset = new double[nAtoms]

        double surfaceArea = 0

        int index = 0
        for (Atom atom : atoms) {
            //isHydrogen[index] = atom.isHydrogen()
            isHydrogen[index] = false
            radii[index] = atom.getVDWType().radius / 2.0 - 0.2//* angToNm / 2.0
            radiiOffset[index] = radii[index] - offset
            volume[index] = 4.0 / 3.0 * Math.PI * pow(radii[index], 3)
            volumeOffset[index] = 4.0 / 3.0 * Math.PI * pow(radiiOffset[index], 3)
            gamma[index] = 1.0
            positions[index][0] = atom.getX() //* angToNm
            positions[index][1] = atom.getY() //* angToNm
            positions[index][2] = atom.getZ() //* angToNm
            index++
        }

        GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen)

        gaussVol.computeTree(positions)
        gaussVol.computeVolume(positions, totalVolume, energy, force, gradV, freeVolume, selfVolume)
        //gaussVol.printTree()


        totalVolume[0]=totalVolume[0]//*nmVolumeToAngVolume
        //TODO: Check Energy Conversion
        energy[0]=energy[0]//*kjToKcal*nmVolumeToAngVolume

        logger.info(String.format(" Total Volume of Molecule: %8.6f, Energy: %8.6f", totalVolume[0], energy[0]))

        int i = 0
        for (Atom atom : atoms) {
            //if (!atom.isHydrogen()) {
                logger.info(String.format(" Atom %s, Force: (%8.6f,%8.6f,%8.6f), GradV: %8.6f, FreeV: %8.6f, SelfV: %8.6f ",
                        atom.toString(), force[i][0], force[i][1], force[i][2], gradV[i], freeVolume[i], selfVolume[i]))
                selfVolumeSum += selfVolume[i] / offset
            //}
            i++
        }

        GaussVol gaussVolOffest = new GaussVol(nAtoms, radiiOffset, volumeOffset, gamma, isHydrogen)

        gaussVolOffest.computeTree(positions)
        gaussVolOffest.computeVolume(positions, totalVolumeOffset, energyOffset, forceOffset, gradVOffset, freeVolumeOffset, selfVolumeOffset)
        //gaussVolOffest.printTree()


        totalVolumeOffset[0]=totalVolumeOffset[0]//*nmVolumeToAngVolume
        //TODO: Check Energy Conversion
        energyOffset[0]=energyOffset[0]//*kjToKcal*nmVolumeToAngVolume

        logger.info(String.format(" Total Volume of Offset: %8.6f, Energy: %8.6f", totalVolumeOffset[0], energyOffset[0]))

        int j = 0
        for (Atom atom : atoms) {
            //if (!atom.isHydrogen()) {
                logger.info(String.format(" Atom %s, Force: (%8.6f,%8.6f,%8.6f), GradV: %8.6f, FreeV: %8.6f, SelfV: %8.6f ",
                        atom.toString(), forceOffset[j][0], forceOffset[j][1], forceOffset[j][2], gradVOffset[j], freeVolumeOffset[j], selfVolumeOffset[j]))
                selfVolumeOffsetSum += selfVolumeOffset[j] / offset
            //}
            j++
        }
        surfaceArea=selfVolumeSum-selfVolumeOffsetSum

        logger.info(String.format("Solvent Accessible Surface Area: %8.6f", surfaceArea))


        return this
    }
}