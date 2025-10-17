//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential.commands;

import edu.rit.pj.ParallelTeam;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.nonbonded.implicit.ConnollyRegion;
import ffx.potential.nonbonded.implicit.GaussVol;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm for each
 * snapshot in an ARC file.
 *
 * Usage:
 *   ffxc VolumeArc [options] &lt;filename&gt;
 */
@Command(name = "VolumeArc", description = " Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm.")
public class VolumeArc extends PotentialCommand {

  private static final double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0);

  @Option(names = {"-c", "--connolly"}, paramLabel = "false",
      description = "Use the Connolly algorithm to compute solvent excluded volume and solvent accessible surface area.")
  private boolean connolly = false;

  @Option(names = {"-m", "--molecular"}, paramLabel = "false",
      description = "For Connolly, compute molecular volume and surface area (instead of SEV/SASA).")
  private boolean molecular = false;

  @Option(names = {"--vdW", "--vanDerWaals"}, paramLabel = "false",
      description = "For Connolly, compute van der Waals volume and surface area (instead of SEV/SASA)")
  private boolean vdW = false;

  @Option(names = {"-p", "--probe"}, paramLabel = "1.4",
      description = "For Connolly, set the exclude radius (SASA) or probe radius (molecular surface). Ignored for vdW.")
  private double probe = 1.4;

  @Option(names = {"-y", "--includeHydrogen"}, paramLabel = "false",
      description = "Include Hydrogen in calculation volume and surface area.")
  private boolean includeHydrogen = false;

  @Option(names = {"-s", "--sigma"}, paramLabel = "false",
      description = "Use sigma radii instead of Rmin.")
  private boolean sigma = false;

  @Option(names = {"-o", "--offset"}, paramLabel = "0.0",
      description = "Add an offset to all atomic radii for GaussVol volume and surface area.")
  private double offset = 0.0;

  @Option(names = {"-v", "--verbose"}, paramLabel = "false",
      description = "Print out all components of volume of molecule and offset.")
  private boolean verbose = false;

  /** The final argument(s) should be one or more filenames. */
  @Parameters(arity = "1", paramLabel = "files",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private java.util.List<String> filenames = null;

  /** JUnit Testing Variables. */
  public double totalVolume = 0.0;
  public double totalSurfaceArea = 0.0;

  public VolumeArc() { super(); }
  public VolumeArc(FFXBinding binding) { super(binding); }
  public VolumeArc(String[] args) { super(args); }

  @Override
  public VolumeArc run() {
    if (!init()) {
      return this;
    }

    if (filenames != null && !filenames.isEmpty()) {
      MolecularAssembly[] assemblies = potentialFunctions.openAll(filenames.get(0));
      activeAssembly = assemblies[0];
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    String modelFilename = activeAssembly.getFile().getAbsolutePath();
    logger.info("\n Calculating volume and surface area for " + modelFilename);

    Atom[] atoms = activeAssembly.getAtomArray();
    int nAtoms = atoms.length;

    SystemFilter systemFilter = potentialFunctions.getFilter();
    if (filenames != null && filenames.size() == 1) {
      String ext = FilenameUtils.getExtension(filenames.get(0));
      if (ext.toUpperCase().contains("ARC")) {
        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
          while (systemFilter.readNext()) {
            if (!connolly) {
              // Inputs for GaussVol.
              double[][] positions = new double[nAtoms][3];
              int index = 0;
              for (Atom atom : atoms) {
                positions[index][0] = atom.getX();
                positions[index][1] = atom.getY();
                positions[index][2] = atom.getZ();
                index++;
              }

              ForceField forceField = activeAssembly.getForceField();
              if (includeHydrogen) {
                forceField.addProperty("GAUSSVOL_HYDROGEN", "true");
              }
              if (sigma) {
                forceField.addProperty("GAUSSVOL_USE_SIGMA", "true");
              }
              if (offset > 0.0) {
                forceField.addProperty("GAUSSVOL_RADII_OFFSET", Double.toString(offset));
              }
              if (!forceField.hasProperty("GAUSSVOL_RADII_SCALE")) {
                forceField.addProperty("GAUSSVOL_RADII_SCALE", Double.toString(1.0));
              }

              ParallelTeam parallelTeam = new ParallelTeam();
              GaussVol gaussVol = new GaussVol(atoms, forceField, parallelTeam);
              gaussVol.computeVolumeAndSA(positions);

              if (verbose) {
                logger.info(format("\n Maximum depth of overlaps in tree: %d", gaussVol.getMaximumDepth()));
                logger.info(format(" Total number of overlaps in tree: %d", gaussVol.getTotalNumberOfOverlaps()));
                double[] radii = gaussVol.getRadii();
                index = 0;
                for (Atom atom : atoms) {
                  logger.info(" " + String.format("%5d %4s: %6.4f", index, atom.getName(), radii[index]));
                  index++;
                }
              }

              logger.info("\n GaussVol Surface Area and Volume for Arc File\n");
              if (sigma) {
                logger.info(format("  Radii:                  Sigma"));
              } else {
                logger.info(format("  Radii:                   Rmin"));
              }
              logger.info(format("  Radii offset:        %8.4f (Ang)", offset));
              logger.info(format("  Include hydrogen:    %8b", includeHydrogen));
              logger.info(format("  Volume:              %8.4f (Ang^3)", gaussVol.getVolume()));
              logger.info(format("  Surface Area:        %8.4f (Ang^2)", gaussVol.getSurfaceArea()));
              totalVolume = gaussVol.getVolume();
              totalSurfaceArea = gaussVol.getSurfaceArea();
            } else {
              double exclude = 0.0;
              if (vdW) {
                exclude = 0.0;
                probe = 0.0;
              } else if (!molecular) {
                exclude = probe;
                probe = 0.0;
              }

              double[] radii = new double[nAtoms];
              int index = 0;
              for (Atom atom : atoms) {
                double r = atom.getVDWType().radius / 2.0;
                if (sigma) {
                  r *= rminToSigma;
                }
                if (!includeHydrogen && atom.isHydrogen()) {
                  r = 0.0;
                }
                radii[index++] = r;
              }

              // VolumeRegion currently limited to 1 thread.
              ParallelTeam parallelTeam = new ParallelTeam(1);
              ConnollyRegion connollyRegion = new ConnollyRegion(atoms, radii, parallelTeam.getThreadCount());
              connollyRegion.setExclude(exclude);
              connollyRegion.setProbe(probe);
              connollyRegion.runVolume();

              if (vdW || (probe == 0.0 && exclude == 0.0)) {
                logger.info("\n Connolly van der Waals Surface Area and Volume for Arc File\n");
              } else if (!molecular) {
                logger.info("\n Connolly Solvent Accessible Surface Area and Solvent Excluded Volume\n");
                logger.info(format("  Exclude radius:      %8.4f (Ang)", exclude));
              } else {
                logger.info("\n Connolly Molecular Surface Area and Volume");
                logger.info(format("  Probe radius:       %8.4f (Ang)", probe));
              }
              if (sigma) {
                logger.info(format("  Radii:                  Sigma"));
              } else {
                logger.info(format("  Radii:                   Rmin"));
              }
              logger.info(format("  Include hydrogen:    %8b", includeHydrogen));
              logger.info(format("  Volume:              %8.4f (Ang^3)", connollyRegion.getVolume()));
              logger.info(format("  Surface Area:        %8.4f (Ang^2)", connollyRegion.getSurfaceArea()));
              totalVolume = connollyRegion.getVolume();
              totalSurfaceArea = connollyRegion.getSurfaceArea();
            }
          }
        }
      }
    }

    return this;
  }
}
