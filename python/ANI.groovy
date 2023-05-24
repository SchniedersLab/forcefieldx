//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import org.graalvm.polyglot.*
import org.graalvm.polyglot.proxy.*

import java.io.OutputStream
import java.io.BufferedOutputStream
import java.io.ByteArrayOutputStream
import java.nio.file.Path
import java.nio.file.Paths

import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The ANI script evaluates the ANI energy of a system.
 * This must be executed using GraalVM and a virtual environment with Torch installed.
 * <br>
 * Usage:
 * <br>
 * ffxc ANI &lt;filename&gt;
 */
@Command(description = " Compute the ANI-2x energy.", name = "ANI-2x")
class ANI extends PotentialScript {

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  private String filename = null

  /**
   * Energy constructor.
   */
  ANI() {
    this(new Binding())
  }

  /**
   * Energy constructor.
   * @param binding The Groovy Binding to use.
   */
  ANI(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  ANI run() {
    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Running Energy on " + filename)

    String JAVA_HOME = System.getProperty("java.home")
    Path graalpy = Paths.get(JAVA_HOME, "ffx_venv", "bin", "graalpy") 
    graalpyString = System.getProperty("graalpy", graalpy.toString());
    graalpy = Paths.get(graalpyString)
    logger.info(" graalpy (-Dgraalpy=path.to.graalpy):             " + graalpy)
    String torchScript = "ANI2x.pt"
    torchScript = System.getProperty("torchscript", torchScript);
    logger.info(" torchscript (-Dtorchscript=path.to.torchscript): " + torchScript);
    Atom[] atoms = activeAssembly.getAtomArray();
    int nAtoms = atoms.length;
    StringBuilder coords = new StringBuilder();
    StringBuilder species = new StringBuilder();
    for (int i=0; i < nAtoms; i++) {
       Atom atom = atoms[i]
       coords.append(format("[%17.15f, %17.15f, %17.15f]", atom.getX(), atom.getY(), atom.getZ())) 
       species.append(format("%d", atom.getAtomicNumber()))
       if (i < nAtoms - 1) {
          coords.append(",\n")
          species.append(",");
       }
    }
    String coordTensor = format("coordinates = torch.tensor([[\n%s]], dtype=torch.double)\n", coords.toString());
    String speciesTensor = format("species = torch.tensor([[%s]], dtype=torch.int64)\n", species.toString());
    double energy = 0.0;
    double[] grad = new double[nAtoms * 3]; 
    try (Context context = Context.newBuilder("python").allowAllAccess(true).
                            option("python.Executable", graalpy.toString()).build()) {
       String torch = "import site\n"
       torch += "import torch\n"
       torch += "ani = torch.jit.load('" + torchScript + "')\n"
       torch += coordTensor
       torch += speciesTensor
       torch += "gradient = ani(species, coordinates)\n"
       logger.fine(torch)

       Value result = context.eval("python", torch);
       Value bindings = context.getBindings("python");
       Value ret = bindings.getMember("gradient")
       for (int i=0; i < nAtoms * 3; i++) {
          grad[i] = ret.getArrayElement(i).asDouble();
       }
       energy = ret.getArrayElement(nAtoms * 3).asDouble();
    } catch (Exception e) {
       logger.info(" Exception:\n" + e.toString()) 
    }

    logger.info(" Energy: " + energy)
    logger.info(" Gradient:\n" + Arrays.toString(grad))

    return this
  }
}

