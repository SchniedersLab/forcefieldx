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
package ffx.xray.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.refine.RefinedParameter;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;
import ffx.xray.XRayEnergy;
import ffx.xray.cli.XrayOptions;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.ArrayList;
import java.util.List;

import static ffx.utilities.StringUtils.parseAtomRanges;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The Gradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Gradient [options] &lt;filename&gt;
 */
@Command(description = " Test the potential energy gradient.", name = "test.Gradient")
public class Gradient extends AlgorithmsCommand {

  @Mixin
  AtomSelectionOptions atomSelectionOptions;

  @Mixin
  XrayOptions xrayOptions;

  /**
   * -d or --dx Finite-difference step size.
   */
  @Option(names = {"-d", "--dx"}, defaultValue = "1.0e-5", paramLabel = "1.0e-5",
      description = "Finite-difference step size.")
  public double dx = 1.0e-5;

  /**
   * --tol or --tolerance Gradient error tolerance.
   */
  @Option(names = {"--tol", "--tolerance"}, defaultValue = "1.0e-3",
      paramLabel = "1.0e-3 kcal/mol/param", description = "Gradient error tolerance.")
  public double tolerance = 1.0e-3;

  /**
   * --params or --gradientAtoms Ranges of degrees of freedom to test [ALL, NONE, Range(s): 1-3,6-N].
   */
  @Option(names = {"--params", "--gradientParams"}, paramLabel = "ALL", defaultValue = "ALL",
      description = "Ranges of degrees of freedom to test [ALL, NONE, Range(s): 1-3,6-N].")
  public String gradientParams = "ALL";


  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames;

  private MolecularAssembly[] molecularAssemblies;
  private DiffractionData diffractionData;
  private XRayEnergy xrayEnergy;
  public int nFailures = 0;

  /**
   * Gradient constructor.
   */
  public Gradient() {
    super();
  }

  /**
   * Gradient constructor.
   *
   * @param binding The Binding to use.
   */
  public Gradient(FFXBinding binding) {
    super(binding);
  }

  /**
   * Gradient constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Gradient(String[] args) {
    super(args);
  }

  /**
   * Execute the script.
   */
  @Override
  public Gradient run() {

    if (!init()) {
      return this;
    }

    xrayOptions.init();

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.getFirst());
      activeAssembly = molecularAssemblies[0];
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
    }

    // Update the active filename
    filename = activeAssembly.getFile().getAbsolutePath();

    // Set active atoms for each assembly.
    for (MolecularAssembly assembly : molecularAssemblies) {
      atomSelectionOptions.setActiveAtoms(assembly);
    }

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties();
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();
    xrayEnergy = new XRayEnergy(diffractionData);

    RefinementModel refinementModel = diffractionData.getRefinementModel();
    RefinementMode refinementMode = refinementModel.getRefinementMode();

    logger.info("\n Testing the gradient of " + filename);
    logger.info(refinementMode.toString());
    logger.info(refinementModel.toString());

    int n = refinementModel.getNumParameters();
    double[] x = new double[n];
    double[] g = new double[n];
    // Collect the analytic gradient.
    refinementModel.getParameters(x);
    xrayEnergy.energyAndGradient(x, g);

    List<RefinedParameter> refinedParameters = refinementModel.getRefinedParameters();
    int numParameters = refinedParameters.size();

    List<Integer> parametersToTest;
    if (gradientParams.equalsIgnoreCase("NONE")) {
      logger.info(" The gradient of no parameters will be evaluated.");
      return this;
    } else if (gradientParams.equalsIgnoreCase("ALL")) {
      logger.info(" Checking the gradient for all parameters.\n");
      parametersToTest = new ArrayList<>();
      for (int i = 0; i < numParameters; i++) {
        parametersToTest.add(i);
      }
    } else {
      parametersToTest = parseAtomRanges(" Parameters", gradientParams, numParameters);
      logger.info(" Checking the gradient for parameters in the range: " + gradientParams + "\n");
    }

    for (Integer parameter : parametersToTest) {
      RefinedParameter refinedParameter = refinedParameters.get(parameter);
      int index = refinedParameter.getIndex();
      int numParams = refinedParameter.getNumberOfParameters();
      double[] grad = new double[numParams];
      double error = 0.0;
      // Collect the gradient
      for (int i = 0; i < numParams; i++) {
        double orig = x[index + i];
        x[index + i] = orig + dx;
        double eplus = xrayEnergy.energy(x);
        x[index + i] = orig - dx;
        double eminus = xrayEnergy.energy(x);
        x[index + i] = orig;
        grad[i] = (eplus - eminus) / (2.0 * dx);
        double diff = grad[i] - g[index + i];
        error += diff * diff;
      }
      error = sqrt(error);
      Atom atom = refinedParameter.getAtom();
      if (numParams == 1) {
        logger.info(format("\n Refinement parameter %d with index %d: %s", parameter + 1, index, atom));
      } else {
        logger.info(format("\n Refinement parameter %d with indices %d - %s: %s",
            parameter + 1, index, index + numParams - 1, atom));
      }
      logger.info(format(" %s", refinedParameter));
      if (error > tolerance) {
        logger.info(format("  Failed: %10.6f", error));
        nFailures++;
      } else {
        logger.info(format("  Passed: %10.6f", error));
      }
      StringBuilder numeric = new StringBuilder("  Numeric:  ");
      StringBuilder analytic = new StringBuilder("  Analytic: ");
      for (int i = 0; i < numParams; i++) {
        numeric.append(format(" %10.6f", grad[i]));
        analytic.append(format(" %10.6f", g[index + i]));
      }
      logger.info(numeric.toString());
      logger.info(analytic.toString());
    }

    logger.info("\n %d / %d failures.".formatted(nFailures, numParameters));

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials = getPotentialsFromAssemblies(molecularAssemblies);
    if (xrayEnergy != null) {
      potentials.add(xrayEnergy);
    }
    return potentials;
  }

}