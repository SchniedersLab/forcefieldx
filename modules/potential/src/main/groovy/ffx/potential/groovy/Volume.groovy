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

import static java.lang.String.format

import org.apache.commons.io.FilenameUtils
import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.pow

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.nonbonded.GaussVol
import ffx.potential.nonbonded.GeneralizedKirkwood

import edu.rit.pj.ParallelTeam
import ffx.potential.nonbonded.implicit.VolumeRegion

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
@Command(description = " Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm.",
        name = "ffxc Volume")
class Volume extends PotentialScript {

    private static final double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0)

    /**
     * -c or --connolly Use the Connolly algorithm to compute volume and surface area (instead of GaussVol).
     */
    @CommandLine.Option(names = ['-c', '--connolly'], paramLabel = "false",
            description = "Use the Connolly algorithm to compute solvent excluded volume and solvent accessible surface area.")
    private boolean connolly = false

//   NOT SUPPORTED YET.
//    /**
//     * -m or --molecular If using the Connolly algorithm, compute molecular volume and surface area (instead of SEV/SASA).
//     */
//    @CommandLine.Option(names = ['-m', '--molecular'], paramLabel = "false",
//            description = "If using the Connolly algorithm, compute molecular volume and surface area (instead of SEV/SASA).")
//    private boolean molecular = false

    /**
     * -p or --probe If using the Connolly algorithm, set the exclude radius (SASA) or probe radius (molecular surface).
     */
    @CommandLine.Option(names = ['-p', '--probe'], paramLabel = "1.4",
            description = "If using the Connolly algorithm, set the exclude radius (SASA) or probe radius (molecular surface).")
    private double probe = 1.4

    /**
     * -y or --includeHydrogen Include Hydrogen in calculation volume and surface area.
     */
    @CommandLine.Option(names = ['-y', '--includeHydrogen'], paramLabel = "false",
            description = "Include Hydrogen in calculation volume and surface area.")
    private boolean includeHydrogen = false

    /**
     * -s or --sigma Use sigma radii instead of Rmin.
     */
    @CommandLine.Option(names = ['-s', '--sigma'], paramLabel = "false",
            description = "Use sigma radii instead of Rmin.")
    private boolean sigma = false

    /**
     * -o or --offset For GaussVol, add an offset to all atomic radii.
     */
    @CommandLine.Option(names = ['-o', '--offset'], paramLabel = "0.0",
            description = "Add an offset to all atomic radii for GaussVol volume and surface area.")
    private double offset = 0.0

    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @CommandLine.Option(names = ['-v', '--verbose'], paramLabel = "false",
            description = "Print out all components of volume of molecule and offset.")
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
    public double totalVolume = 0.0
    public double totalSurfaceArea = 0.0

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

        CompositeConfiguration properties = activeAssembly.getProperties()
        double solventPressure = properties.getDouble("solvent-pressure", GeneralizedKirkwood.DEFAULT_SOLVENT_PRESSURE)
        double surfaceTension = properties.getDouble("surface-tension", GeneralizedKirkwood.DEFAULT_CAVDISP_SURFACE_TENSION)
        double crossOver = properties.getDouble("cross-over", GeneralizedKirkwood.DEFAULT_CROSSOVER)

        if (!connolly) {
            // Input
            double[] radii = new double[nAtoms]
            boolean[] isHydrogen = new boolean[nAtoms]
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
                radii[index] = atom.getVDWType().radius / 2.0
                if (sigma) {
                    radii[index] *= rminToSigma;
                }
                radii[index] += offset
                volume[index] = fourThirdsPI * pow(radii[index], 3)
                positions[index][0] = atom.getX()
                positions[index][1] = atom.getY()
                positions[index][2] = atom.getZ()
                index++
            }

            // Run Volume calculation to get vdw volume of molecule
            GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen)
            gaussVol.setSolventPressure(solventPressure)
            gaussVol.setSurfaceTension(surfaceTension)
            gaussVol.setCrossOver(crossOver)

            gaussVol.computeVolumeAndSA(positions)

            if (verbose) {
                logger.info(format("\n Maximum depth of overlaps in tree: %d", gaussVol.getMaximumDepth()))
                logger.info(format(" Total number of overlaps in tree: %d", gaussVol.getTotalNumberOfOverlaps()))

                //gaussVol.printTree()

                index = 0
                for (Atom atom : atoms) {
                    logger.info(" Radius for atom " + atom.name + ": " + radii[index] + "\n")
                    index++
                }
            }

            logger.info("\n GaussVol Surface Area and Volume\n")
            if (sigma) {
                logger.info(format("  Radii:                  Sigma"))
            } else {
                logger.info(format("  Radii:                   Rmin"))
            }
            logger.info(format("  Radii offset:        %8.4f (Ang)", offset))
            logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
            logger.info(format("\n  Volume:              %8.4f (Ang^3)", gaussVol.getVolume()))
            logger.info(format("  Solvent Pressure:    %8.4f (kcal/mol/Ang^3)", gaussVol.getSolventPressure()))
            logger.info(format("  Volume Energy:       %8.4f (kcal/mol)", gaussVol.getVolumeEnergy()))
            logger.info(format("\n  Surface Area:        %8.4f (Ang^2)", gaussVol.getSurfaceArea()))
            logger.info(format("  Surface Tension:     %8.4f (kcal/mol/Ang^2)", gaussVol.getSurfaceTension()))
            logger.info(format("  Surface Area Energy: %8.4f (kcal/mol)", gaussVol.getSurfaceAreaEnergy()))
            logger.info(format("\n  Effective Radius:    %8.4f (Ang)", gaussVol.getEffectiveRadius()))
            logger.info(format("  Cross-over Radius:   %8.4f (Ang)", gaussVol.getCrossOver()))
            logger.info(format("  Volume + SA Energy:  %8.4f (kcal/mol)", gaussVol.getEnergy()))

            // Set JUnit testing variables based on output volume and surface area
            totalVolume = gaussVol.getVolume()
            totalSurfaceArea = gaussVol.getSurfaceArea()
        } else {
            // Connolly algorithm.

            // For molecular, use the chosen probe and set exclude to 0.0.
            double exclude = 0.0

            // For SEV/SASA, set exclude to the chosen probe, and zero the probe.
            // if (!molecular) {
                exclude = probe
                probe = 0.0
            // }

            // Connolly Volume and Surface Area.
            double[] radii = new double[nAtoms]
            double[] x = new double[nAtoms]
            double[] y = new double[nAtoms]
            double[] z = new double[nAtoms]

            int index = 0
            for (Atom atom : atoms) {
                radii[index] = atom.getVDWType().radius / 2.0
                if (sigma) {
                    radii[index] *= rminToSigma;
                }
                boolean hydrogen = atom.isHydrogen()
                if (!includeHydrogen && hydrogen) {
                    radii[index] = 0.0
                }
                x[index] = atom.getX()
                y[index] = atom.getY()
                z[index] = atom.getZ()
                index++
            }

            // Note that the VolumeRegion code is currently limited to a single thread.
            ParallelTeam parallelTeam = new ParallelTeam(1)
            VolumeRegion volumeRegion = new VolumeRegion(atoms, x, y, z, radii, parallelTeam.getThreadCount())
            volumeRegion.setSolventPressure(solventPressure)
            volumeRegion.setSurfaceTension(surfaceTension)
            volumeRegion.setCrossOver(crossOver)
            // For solvent excluded volume.
            volumeRegion.setExclude(exclude)
            // For molecular surface.
            volumeRegion.setProbe(probe)
            try {
                parallelTeam.execute(volumeRegion)
            } catch (Exception e) {
                logger.warning(" Exception executing Volume region");
            }

            // if (!molecular) {
                logger.info("\n Connolly Solvent Accessible Surface Area and Solvent Excluded Volume\n")
                logger.info(format("  Exclude radius:      %8.4f (Ang)", exclude))
            // } else {
            //    logger.info("\n Connolly Molecular Surface Area and Volume")
            //    logger.info(format("  Probe radius:       %8.4f (Ang)", probe))
            // }

            if (sigma) {
                logger.info(format("  Radii:                  Sigma"))
            } else {
                logger.info(format("  Radii:                   Rmin"))
            }
            logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
            logger.info(format("  Volume:              %8.4f (Ang^3)", volumeRegion.getVolume()))
            logger.info(format("  Solvent Pressure:    %8.4f (kcal/mol/Ang^3)", volumeRegion.getSolventPressure()))
            logger.info(format("  Volume Energy:       %8.4f (kcal/mol)", volumeRegion.getVolumeEnergy()))
            logger.info(format("\n  Surface Area:        %8.4f (Ang^2)", volumeRegion.getSurfaceArea()))
            logger.info(format("  Surface Tension:     %8.4f (kcal/mol/Ang^2)", volumeRegion.getSurfaceTension()))
            logger.info(format("  Surface Area Energy: %8.4f (kcal/mol)", volumeRegion.getSurfaceAreaEnergy()))
            logger.info(format("\n  Effective Radius:    %8.4f (Ang)", volumeRegion.getEffectiveRadius()))
            logger.info(format("  Cross-over Radius:   %8.4f (Ang)", volumeRegion.getCrossOver()))
            logger.info(format("  Volume + SA Energy:  %8.4f (kcal/mol)", volumeRegion.getEnergy()))

            // Set JUnit testing variables based on output volume and surface area
            totalVolume = volumeRegion.getVolume()
            totalSurfaceArea = volumeRegion.getSurfaceArea()
        }

        return this
    }
}
