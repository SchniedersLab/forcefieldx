//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import ffx.potential.GraalPyANI;
import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.logging.Level;

/**
 * The ANI command evaluates the ANI2x energy of a system.
 * This must be executed using GraalVM and a virtual environment with Torch installed.
 * <br>
 * Usage:
 * <br>
 * ffxc ANI &lt;filename&gt;
 */
@Command(description = " Compute the ANI-2x energy.", name = "ANI-2x")
public class ANI extends PotentialCommand {

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  @Option(names = {"-p", "--printPython"},
      description = "Print the Python code used to evaluate ANI-2x.")
  private boolean printPython;

  public ANI() {
    super();
  }

  public ANI(FFXBinding binding) {
    super(binding);
  }

  public ANI(String[] args) {
    super(args);
  }

  @Override
  public ANI run() {
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    logger.info("\n Running ANI-2x Energy on " + filename);

    if (printPython) {
      logger.info("\n Python code used to evaluate ANI-2x:\n");
      logger.info(GraalPyANI.getPythonSource());
      logger.info("");
    }

    Path graalpy;
    try {
      graalpy = GraalPyANI.getGraalPyPath();
    } catch (IllegalStateException e) {
      logger.severe(" " + e.getMessage());
      return this;
    }
    logger.info(" graalpy (-Dgraalpy=path.to.graalpy):             " + graalpy);
    if (!Files.isDirectory(graalpy)) {
      logger.severe(" GraalPy resource directory was not found: " + graalpy);
      return this;
    }

    Path torchScript = GraalPyANI.getTorchScriptPath();
    logger.info(" torchscript (-Dtorchscript=path.to.torchscript): " + torchScript);
    if (!Files.isRegularFile(torchScript)) {
      logger.severe(" TorchScript model was not found: " + torchScript);
      return this;
    }

    // Collect atomic number and coordinates for each atom.
    Atom[] atoms = activeAssembly.getAtomArray();
    int nAtoms = atoms.length;
    int[] species = new int[nAtoms];
    double[] coordinates = new double[nAtoms * 3];
    for (int i = 0; i < nAtoms; i++) {
      Atom a = atoms[i];
      species[i] = a.getAtomicNumber();
      int index = 3 * i;
      coordinates[index] = a.getX();
      coordinates[index + 1] = a.getY();
      coordinates[index + 2] = a.getZ();
    }
    double[] grad = new double[nAtoms * 3];
    try {
      GraalPyANI ani = new GraalPyANI(species);
      double energy = ani.energyAndGradient(coordinates, grad);
      logger.info(" ANI-2x Energy (Hartree): " + energy);
    } catch (org.graalvm.polyglot.PolyglotException e) {
      logger.log(Level.SEVERE, " ANI-2x evaluation failed.", e);
    }

    return this;
  }

}
