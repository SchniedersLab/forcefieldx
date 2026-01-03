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
package ffx.xray.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.TimerOptions;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementEnergy;
import ffx.xray.cli.XrayOptions;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Parameters;

import java.util.List;

import static java.lang.String.format;

/**
 * The X-ray Timer script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Timer [options] &lt;filename&gt;
 */
@Command(description = " Time calculation of the X-ray target.", name = "xray.Timer")
public class Timer extends AlgorithmsCommand {

  @Mixin
  private TimerOptions timerOptions;

  @Mixin
  private XrayOptions xrayOptions;

  @Mixin
  AtomSelectionOptions atomSelectionOptions;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames;

  private DiffractionData diffractionData;
  private MolecularAssembly[] molecularAssemblies;

  private RefinementEnergy refinementEnergy;

  /**
   * Timer constructor.
   */
  public Timer() {
    super();
  }

  /**
   * Timer constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Timer(String[] args) {
    super(args);
  }

  /**
   * Timer constructor.
   *
   * @param binding The Binding to use.
   */
  public Timer(FFXBinding binding) {
    super(binding);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Timer run() {

    if (!init()) {
      return this;
    }

    xrayOptions.init();

    // Set the number of threads (needs to be done before opening the files).
    if (timerOptions.getThreads() > 0) {
      System.setProperty("pj.nt", Integer.toString(timerOptions.getThreads()));
    }

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0));
      activeAssembly = molecularAssemblies[0];
      filename = filenames.get(0);
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
      filename = activeAssembly.getFile().getAbsolutePath();
    }

    logger.info("\n Running xray.Timer on " + filename);

    // Apply active atom flags.
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      atomSelectionOptions.setActiveAtoms(molecularAssembly);
    }

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties();
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);
    refinementEnergy = xrayOptions.toXrayEnergy(diffractionData);

    // Print the initial energy of each conformer.
    algorithmFunctions.energy(molecularAssemblies);

    int n = refinementEnergy.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    refinementEnergy.getCoordinates(x);
    Potential energy = refinementEnergy.getDataEnergy();

    logger.info("\n Beginning Timings\n");

    int nEvals = timerOptions.getIterations();
    boolean gradient = !timerOptions.getNoGradient();
    long minTime = Long.MAX_VALUE;
    double sumTime2 = 0.0;
    int halfnEvals = (nEvals % 2 == 1) ? (nEvals / 2) : (nEvals / 2) - 1; // Halfway point
    for (int i = 0; i < nEvals; i++) {
      long time = -System.nanoTime();
      double e;
      if (gradient) {
        e = energy.energyAndGradient(x, g);
      } else {
        e = energy.energy(x);
      }
      time += System.nanoTime();
      if (gradient) {
        logger.info(format(" Energy & Gradient: %12.6f in %6.3f (sec)", e, time * 1.0E-9));
      } else {
        logger.info(format(" Energy: %12.6f in %6.3f (sec)", e, time * 1.0E-9));
      }

      minTime = time < minTime ? time : minTime;
      if (i >= nEvals / 2) {
        double time2 = time * 1.0E-9;
        sumTime2 += (time2 * time2);
      }
    }

    ++halfnEvals;
    double rmsTime = Math.sqrt(sumTime2 / halfnEvals);
    logger.info(format("\n Minimum time:           %6.3f (sec)", minTime * 1.0E-9));
    logger.info(format(" RMS time (latter half): %6.3f (sec)", rmsTime));

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }

  @Override
  public boolean destroyPotentials() {
    return diffractionData == null ? true : diffractionData.destroy();
  }
}
