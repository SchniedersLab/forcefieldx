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

import edu.rit.pj.ParallelTeam
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.nonbonded.implicit.ConnollyRegion
import ffx.potential.nonbonded.implicit.GaussVol
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.pow

/**
 * The SaveAsP1 script expands a specified file to P1
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsP1 [options] &lt;filename&gt;
 */
@Command(description = " Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm.",
        name = "ffxc VolumeArc")
class VolumeArc extends PotentialScript {

    private static final double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0)

    /**
     * -c or --connolly Use the Connolly algorithm to compute volume and surface area (instead of GaussVol).
     */
    @CommandLine.Option(names = ['-c', '--connolly'], paramLabel = "false",
            description = "Use the Connolly algorithm to compute solvent excluded volume and solvent accessible surface area.")
    private boolean connolly = false

    /**
     * -m or --molecular For Connolly, compute molecular volume and surface area (instead of SEV/SASA).
     */
    @CommandLine.Option(names = ['-m', '--molecular'], paramLabel = "false",
            description = "For Connolly, compute molecular volume and surface area (instead of SEV/SASA).")
    private boolean molecular = false

    /**
     * --vdW or --vanDerWaals For Connolly, compute van der Waals volume and surface area (instead of SEV/SASA).
     */
    @CommandLine.Option(names = ['--vdW', '--vanDerWaals'], paramLabel = "false",
            description = "For Connolly, compute van der Waals volume and surface area (instead of SEV/SASA)")
    private boolean vdW = false

    /**
     * -p or --probe For Connolly, set the exclude radius (SASA) or probe (molecular surface). Ignored for vdW.
     */
    @CommandLine.Option(names = ['-p', '--probe'], paramLabel = "1.4",
            description = "For Connolly, set the exclude radius (SASA) or probe radius (molecular surface). Ignored for vdW.")
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

    /**
     * JUnit Testing Variables
     */
    public double totalVolume = 0.0
    public double totalSurfaceArea = 0.0

    /**
     * Volume Constructor.
     */
    VolumeArc() {
        this(new Binding())
    }

    /**
     * Volume Constructor.
     * @param binding Groovy Binding to use.
     */
    VolumeArc(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    @Override
    VolumeArc run() {
        if (!init()) {
            return null
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = potentialFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return null
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Calculating volume and surface area for " + modelFilename)

        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (filenames.size() == 1) {
            String ext = FilenameUtils.getExtension(filenames.get(0))
            if (ext.toUpperCase().contains("ARC")) {
                if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
                    while (systemFilter.readNext()) {


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
                                    radii[index] *= rminToSigma
                                }
                                radii[index] += offset
                                volume[index] = fourThirdsPI * pow(radii[index], 3)
                                positions[index][0] = atom.getX()
                                positions[index][1] = atom.getY()
                                positions[index][2] = atom.getZ()
                                index++
                            }

                            // Run Volume calculation to get vdw volume of molecule
                            ParallelTeam parallelTeam = new ParallelTeam()
                            GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen, parallelTeam)
                            gaussVol.computeVolumeAndSA(positions)

                            if (verbose) {
                                logger.info(format("\n Maximum depth of overlaps in tree: %d", gaussVol.getMaximumDepth()))
                                logger.info(
                                        format(" Total number of overlaps in tree: %d", gaussVol.getTotalNumberOfOverlaps()))

                                //gaussVol.printTree()

                                index = 0
                                for (Atom atom : atoms) {
                                    logger.info(" Radius for atom " + atom.name + ": " + radii[index] + "\n")
                                    index++
                                }
                            }

                            logger.info("\n GaussVol Surface Area and Volume for Arc File\n")
                            if (sigma) {
                                logger.info(format("  Radii:                  Sigma"))
                            } else {
                                logger.info(format("  Radii:                   Rmin"))
                            }
                            logger.info(format("  Radii offset:        %8.4f (Ang)", offset))
                            logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
                            logger.info(format("  Volume:              %8.4f (Ang^3)", gaussVol.getVolume()))
                            logger.info(format("  Surface Area:        %8.4f (Ang^2)", gaussVol.getSurfaceArea()))
                            // Set JUnit testing variables based on output volume and surface area
                            totalVolume = gaussVol.getVolume()
                            totalSurfaceArea = gaussVol.getSurfaceArea()
                        } else {
                            // For Connolly molecular volume & surface area, use the chosen probe and set exclude to 0.0.
                            double exclude = 0.0

                            if (vdW) {
                                // For Connolly vdW, both exclude & probe are zero.
                                exclude = 0.0
                                probe = 0.0
                            } else if (!molecular) {
                                // For Connolly SEV/SASA, set exclude to the chosen probe, and zero the probe.
                                exclude = probe
                                probe = 0.0
                            }

                            double[] radii = new double[nAtoms]
                            int index = 0
                            for (Atom atom : atoms) {
                                radii[index] = atom.getVDWType().radius / 2.0
                                if (sigma) {
                                    radii[index] *= rminToSigma
                                }
                                boolean hydrogen = atom.isHydrogen()
                                if (!includeHydrogen && hydrogen) {
                                    radii[index] = 0.0
                                }
                                index++
                            }

                            // Note that the VolumeRegion code is currently limited to a single thread.
                            ParallelTeam parallelTeam = new ParallelTeam(1)
                            ConnollyRegion connollyRegion = new ConnollyRegion(atoms, radii,
                                    parallelTeam.getThreadCount())
                            // For solvent excluded volume.
                            connollyRegion.setExclude(exclude)
                            // For molecular surface.
                            connollyRegion.setProbe(probe)
                            connollyRegion.runVolume()

                            if (vdW || (probe == 0.0 && exclude == new Double(0.0))) {
                                logger.info("\n Connolly van der Waals Surface Area and Volume for Arc File\n")
                            } else if (!molecular) {
                                logger.info("\n Connolly Solvent Accessible Surface Area and Solvent Excluded Volume\n")
                                logger.info(format("  Exclude radius:      %8.4f (Ang)", exclude))
                            } else {
                                logger.info("\n Connolly Molecular Surface Area and Volume")
                                logger.info(format("  Probe radius:       %8.4f (Ang)", probe))
                            }
                            if (sigma) {
                                logger.info(format("  Radii:                  Sigma"))
                            } else {
                                logger.info(format("  Radii:                   Rmin"))
                            }
                            logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
                            logger.info(format("  Volume:              %8.4f (Ang^3)", connollyRegion.getVolume()))
                            logger.info(format("  Surface Area:        %8.4f (Ang^2)", connollyRegion.getSurfaceArea()))

                            // Set JUnit testing variables based on output volume and surface area
                            totalVolume = connollyRegion.getVolume()
                            totalSurfaceArea = connollyRegion.getSurfaceArea()
                        }
                    }
                }
            }
        }
        // END HERE
        return this
    }
}
