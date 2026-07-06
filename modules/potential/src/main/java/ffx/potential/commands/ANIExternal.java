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

import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.parseDouble;

/**
 * The ANIExternal command evaluates the ANI2x energy of a system using an external GraalPy
 * subprocess.
 */
@Command(description = " Compute the ANI-2x energy via an external GraalPy subprocess.",
    name = "ANIExternal")
public class ANIExternal extends PotentialCommand {

  private static final String PYTHON_SCRIPT = """
      import os
      import sys
      import torch

      def flatten(value):
          if isinstance(value, (list, tuple)):
              flattened = []
              for item in value:
                  flattened.extend(flatten(item))
              return flattened
          return [value]

      input_path = sys.argv[1]
      output_path = sys.argv[2]

      with open(input_path, encoding="utf-8") as handle:
          torch_script = handle.readline().rstrip("\\n")
          n_atoms = int(handle.readline().strip())
          species = [int(value) for value in handle.readline().split()]
          coords = []
          for _ in range(n_atoms):
              coords.append([float(value) for value in handle.readline().split()])

      species_tensor = torch.tensor([species], dtype=torch.int64)
      coordinates_tensor = torch.tensor([coords], dtype=torch.double)
      ani = torch.jit.load(torch_script)
      gradient = ani(species_tensor, coordinates_tensor)
      if hasattr(gradient, "detach"):
          gradient = gradient.detach().cpu()
      if hasattr(gradient, "tolist"):
          gradient = gradient.tolist()
      flat = flatten(gradient)
      energy = flat[-1]
      gradients = flat[:-1]

      with open(output_path, "w", encoding="utf-8") as handle:
          handle.write(f"{energy}\\n")
          handle.write(" ".join(str(value) for value in gradients))
          handle.write("\\n")
          handle.flush()

      # Torch may leave background activity that upsets GraalPy shutdown on exit.
      os._exit(0)
      """;

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  public ANIExternal() {
    super();
  }

  public ANIExternal(FFXBinding binding) {
    super(binding);
  }

  public ANIExternal(String[] args) {
    super(args);
  }

  @Override
  public ANIExternal run() {
    if (!init()) {
      return this;
    }

    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    filename = activeAssembly.getFile().getAbsolutePath();
    logger.info("\n Running External ANI Energy on " + filename);

    String ffxHome = System.getProperty("basedir");
    Path defaultGraalPyRoot = Paths.get(ffxHome, "python-resources");
    String graalpyProperty = System.getProperty("graalpy", defaultGraalPyRoot.toString());
    Path graalpyPath = Paths.get(graalpyProperty);
    Path graalpyExecutable = resolvePythonExecutable(graalpyPath);
    logger.info(" graalpy (-Dgraalpy=path.to.graalpy):             " + graalpyPath);
    logger.info(" executable:                                      " + graalpyExecutable);

    String torchScriptProperty = System.getProperty("torchscript", "ANI2x.pt");
    Path torchScript = resolveTorchScript(ffxHome, torchScriptProperty);
    logger.info(" torchscript (-Dtorchscript=path.to.torchscript): " + torchScript);

    if (!Files.isExecutable(graalpyExecutable)) {
      logger.severe(" External GraalPy executable was not found: " + graalpyExecutable);
      return this;
    }
    if (!Files.isRegularFile(torchScript)) {
      logger.severe(" TorchScript model was not found: " + torchScript);
      return this;
    }

    Atom[] atoms = activeAssembly.getAtomArray();
    int nAtoms = atoms.length;
    int[] species = new int[nAtoms];
    double[][] coords = new double[nAtoms][3];
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      species[i] = atom.getAtomicNumber();
      coords[i][0] = atom.getX();
      coords[i][1] = atom.getY();
      coords[i][2] = atom.getZ();
    }

    Path inputFile = null;
    Path outputFile = null;
    try {
      Path tempDir = baseDir == null ? null : baseDir.toPath();
      inputFile = tempDir == null
          ? Files.createTempFile("ani-external-input", ".txt")
          : Files.createTempFile(tempDir, "ani-external-input", ".txt");
      outputFile = tempDir == null
          ? Files.createTempFile("ani-external-output", ".txt")
          : Files.createTempFile(tempDir, "ani-external-output", ".txt");

      writeInputFile(inputFile, torchScript, species, coords);
      runSubprocess(graalpyExecutable, inputFile, outputFile, Paths.get(ffxHome, "python"));
      readAndLogResults(outputFile, nAtoms);
    } catch (Exception e) {
      logger.severe(" External ANI subprocess failed: " + e);
    } finally {
      deleteIfExists(inputFile);
      deleteIfExists(outputFile);
    }

    return this;
  }

  private Path resolvePythonExecutable(Path graalpyPath) {
    if (Files.isRegularFile(graalpyPath)) {
      return graalpyPath;
    }

    boolean isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");
    String exeDir = isWindows ? "Scripts" : "bin";
    String[] candidateNames = isWindows
        ? new String[]{"python.exe", "graalpy.exe"}
        : new String[]{"python", "graalpy"};

    List<Path> candidates = new ArrayList<>();
    for (String candidateName : candidateNames) {
      candidates.add(graalpyPath.resolve("venv").resolve(exeDir).resolve(candidateName));
      candidates.add(graalpyPath.resolve(exeDir).resolve(candidateName));
    }

    for (Path candidate : candidates) {
      if (Files.isExecutable(candidate)) {
        return candidate;
      }
    }

    return candidates.getFirst();
  }

  private Path resolveTorchScript(String ffxHome, String torchScriptProperty) {
    Path torchScript = Paths.get(torchScriptProperty);
    if (torchScript.isAbsolute() || Files.exists(torchScript)) {
      return torchScript.toAbsolutePath();
    }

    Path pythonDirTorchScript = Paths.get(ffxHome, "python", torchScriptProperty);
    if (Files.exists(pythonDirTorchScript)) {
      return pythonDirTorchScript.toAbsolutePath();
    }

    return torchScript.toAbsolutePath();
  }

  private void writeInputFile(Path inputFile, Path torchScript, int[] species, double[][] coords)
      throws IOException {
    List<String> lines = new ArrayList<>();
    lines.add(torchScript.toAbsolutePath().toString());
    lines.add(Integer.toString(species.length));

    StringBuilder speciesLine = new StringBuilder();
    for (int i = 0; i < species.length; i++) {
      if (i > 0) {
        speciesLine.append(' ');
      }
      speciesLine.append(species[i]);
    }
    lines.add(speciesLine.toString());

    for (double[] coord : coords) {
      lines.add(coord[0] + " " + coord[1] + " " + coord[2]);
    }

    Files.write(inputFile, lines, StandardCharsets.UTF_8);
  }

  private void runSubprocess(Path graalpyExecutable, Path inputFile, Path outputFile, Path workingDir)
      throws IOException, InterruptedException {
    ProcessBuilder processBuilder = new ProcessBuilder(
        graalpyExecutable.toAbsolutePath().toString(),
        "-c",
        PYTHON_SCRIPT,
        inputFile.toAbsolutePath().toString(),
        outputFile.toAbsolutePath().toString());
    processBuilder.directory(workingDir.toFile());
    processBuilder.redirectErrorStream(true);

    Process process = processBuilder.start();
    try (BufferedReader reader = process.inputReader(StandardCharsets.UTF_8)) {
      String line;
      while ((line = reader.readLine()) != null) {
        logger.info(line);
      }
    }

    int exitCode = process.waitFor();
    if (exitCode != 0) {
      throw new IOException(" External ANI subprocess exited with code " + exitCode);
    }
  }

  private void readAndLogResults(Path outputFile, int nAtoms) throws IOException {
    List<String> lines = Files.readAllLines(outputFile, StandardCharsets.UTF_8);
    if (lines.isEmpty()) {
      throw new IOException(" External ANI subprocess produced no output.");
    }

    double energy = parseDouble(lines.getFirst().trim());
    if (lines.size() > 1 && !lines.get(1).isBlank()) {
      String[] gradientTokens = lines.get(1).trim().split("\\s+");
      if (gradientTokens.length != nAtoms * 3) {
        logger.warning(
            " Expected " + (nAtoms * 3) + " gradient elements but received " + gradientTokens.length);
      }
    }
    logger.info(" ANI-2x Energy (Hartree): " + energy);
  }

  private void deleteIfExists(Path path) {
    if (path == null) {
      return;
    }
    try {
      Files.deleteIfExists(path);
    } catch (IOException e) {
      logger.fine(" Could not delete temporary file " + path + ": " + e.getMessage());
    }
  }
}
