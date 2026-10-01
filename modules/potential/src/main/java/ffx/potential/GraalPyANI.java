// ******************************************************************************
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
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
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
// executable, regardless of the license terms of your choice, provided that
// you also meet, for each linked independent module, the terms and conditions
// of the license of that module. An independent module is a module which is
// not derived from or based on this library. If you modify this library, you
// may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// ******************************************************************************
package ffx.potential;

import org.graalvm.polyglot.Context;
import org.graalvm.polyglot.HostAccess;
import org.graalvm.polyglot.PolyglotAccess;
import org.graalvm.polyglot.PolyglotException;
import org.graalvm.polyglot.Value;
import org.graalvm.python.embedding.GraalPyResources;

import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Persistent GraalPy evaluator for the ANI-2x TorchScript model.
 */
public final class GraalPyANI {

  private static final Logger logger = Logger.getLogger(GraalPyANI.class.getName());

  private static final String TORCH_EVALUATOR = """
      import polyglot
      import torch

      species = [int(atomic_number) for atomic_number in polyglot.import_value('species')]
      torchScript = str(polyglot.import_value('torchScript'))
      speciesTensor = torch.tensor([species], dtype=torch.int64)
      ani = torch.jit.load(torchScript)

      def evaluate(coordinates):
          coordinates = [float(component) for component in coordinates]
          coordinatesTensor = torch.tensor([coordinates], dtype=torch.double).reshape(1, -1, 3)
          return ani(speciesTensor, coordinatesTensor)
      """.stripIndent();

  private final int coordinateCount;
  private final Context context;
  private final Value evaluate;

  /**
   * Initialize an ANI-2x evaluator and load its TorchScript model.
   *
   * @param species Atomic numbers in the model's atom order.
   */
  public GraalPyANI(int[] species) {
    this(species, getTorchScriptPath(), getGraalPyPath());
  }

  private GraalPyANI(int[] species, Path torchScript, Path graalpy) {
    if (!Files.isRegularFile(torchScript)) {
      throw new IllegalArgumentException("ANI-2x TorchScript model was not found: " + torchScript);
    }
    if (!Files.isDirectory(graalpy)) {
      throw new IllegalArgumentException("GraalPy resource directory was not found: " + graalpy);
    }

    coordinateCount = 3 * species.length;
    context = createGraalPyContext(graalpy);
    try {
      Value polyglotBindings = context.getPolyglotBindings();
      polyglotBindings.putMember("species", species.clone());
      polyglotBindings.putMember("torchScript", torchScript.toString());
      context.eval("python", TORCH_EVALUATOR);
      evaluate = context.getBindings("python").getMember("evaluate");
    } catch (PolyglotException e) {
      logger.log(Level.SEVERE, "Unable to initialize the ANI-2x evaluator.", e);
      throw new IllegalStateException("Unable to initialize the ANI-2x evaluator.", e);
    }
  }

  /**
   * Return the Python source used by the evaluator.
   *
   * @return ANI-2x GraalPy evaluator source.
   */
  public static String getPythonSource() {
    return TORCH_EVALUATOR;
  }

  /**
   * Return the GraalPy resource directory configured for ANI-2x.
   *
   * @return GraalPy resource directory.
   */
  public static Path getGraalPyPath() {
    String ffxHome = System.getProperty("basedir");
    if (ffxHome == null || ffxHome.isBlank()) {
      throw new IllegalStateException("FFX base directory is not configured.");
    }
    Path defaultGraalPy = Paths.get(ffxHome, "python-resources");
    return Paths.get(System.getProperty("graalpy", defaultGraalPy.toString())).toAbsolutePath();
  }

  /**
   * Return the ANI-2x TorchScript model configured for evaluation.
   *
   * @return ANI-2x TorchScript model.
   */
  public static Path getTorchScriptPath() {
    return Paths.get(System.getProperty("torchscript", "ANI2x.pt")).toAbsolutePath();
  }

  /**
   * Evaluate the ANI-2x energy.
   *
   * @param coordinates Cartesian coordinates in Angstroms.
   * @return Energy in Hartree.
   */
  public synchronized double energy(double[] coordinates) {
    return evaluate(coordinates, null);
  }

  /**
   * Evaluate the ANI-2x energy and gradient.
   *
   * @param coordinates Cartesian coordinates in Angstroms.
   * @param gradient Gradient destination.
   * @return Energy in Hartree.
   */
  public synchronized double energyAndGradient(double[] coordinates, double[] gradient) {
    if (gradient.length < coordinateCount) {
      throw new IllegalArgumentException(
          "ANI-2x requires a gradient array of at least " + coordinateCount + " elements.");
    }
    return evaluate(coordinates, gradient);
  }

  private double evaluate(double[] coordinates, double[] gradient) {
    if (coordinates.length != coordinateCount) {
      throw new IllegalArgumentException(
          "ANI-2x requires " + coordinateCount + " coordinates, but received " + coordinates.length + '.');
    }

    try {
      Value result = evaluate.execute(coordinates);
      int expectedElements = coordinateCount + 1;
      if (!result.hasArrayElements() || result.getArraySize() != expectedElements) {
        throw new IllegalStateException(
            "ANI-2x returned " + (result.hasArrayElements() ? result.getArraySize() : "a non-array")
                + "; expected " + expectedElements + " energy and gradient elements.");
      }

      if (gradient != null) {
        for (int i = 0; i < coordinateCount; i++) {
          gradient[i] = result.getArrayElement(i).asDouble();
        }
      }
      return result.getArrayElement(coordinateCount).asDouble();
    } catch (PolyglotException e) {
      logger.log(Level.SEVERE, "ANI-2x evaluation failed.", e);
      throw e;
    }
  }

  /**
   * Create the persistent GraalPy context used for ANI-2x evaluations.
   *
   * <p>The context intentionally remains open because PyTorch may retain native work after an
   * evaluation. Keeping it alive also avoids repeatedly loading the ANI-2x TorchScript model.</p>
   */
  private static Context createGraalPyContext(Path graalpy) {
    return Context.newBuilder().allowNativeAccess(true)
        .allowHostAccess(HostAccess.ALL)
        .allowPolyglotAccess(PolyglotAccess.ALL)
        .allowExperimentalOptions(true)
        .apply(GraalPyResources.forExternalDirectory(graalpy))
        .option("python.BackgroundGCTask", "false")
        .option("python.NoAsyncActions", "true")
        .option("python.WarnExperimentalFeatures", "false")
        .logHandler(OutputStream.nullOutputStream()).build();
  }
}
