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

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.TimerOptions
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.sqrt

/**
 * The Timer script evaluates the wall clock time for energy and forces.
 * <br>
 * Usage:
 * <br>
 * ffxc Timer [options] &lt;filename&gt;
 */
@Command(description = " Time energy evaluations.", name = "ffxc Timer")
class Timer extends PotentialScript {

  /**
   * Mix in Timing Options.
   */
  @Mixin
  private TimerOptions timer

  /**
   * The final argument is an atomic coordinate file in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  private ForceFieldEnergy energy = null

  /**
   * Timer Constructor.
   */
  Timer() {
    this(new Binding())
  }

  /**
   * Timer Constructor.
   * @param binding Groovy Binding to use.
   */
  Timer(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Timer run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Set the number of threads.
    if (timer.threads > 0) {
      System.setProperty("pj.nt", Integer.toString(timer.threads))
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    if (timer.noGradient) {
      logger.info("\n Timing energy for " + filename)
    } else {
      logger.info("\n Timing energy and gradient for " + filename)
    }

    // The number of iterations.
    int nEvals = timer.iterations

    energy = activeAssembly.getPotentialEnergy()

    long minTime = Long.MAX_VALUE
    double sumTime2 = 0.0
    // Halfway point
    int halfnEvals = (int) ((nEvals % 2 == 1) ? (nEvals / 2) : (nEvals / 2) - 1)
    for (int i = 0; i < nEvals; i++) {
      long time = -System.nanoTime()
      double e = energy.energy(!timer.noGradient, timer.getVerbose())
      time += System.nanoTime()
      if (!timer.getVerbose()) {
        logger.info(format(" Energy %16.8f in %6.3f (sec)", e, time * 1.0E-9))
      }
      minTime = time < minTime ? time : minTime
      if (i >= (int) (nEvals / 2)) {
        double time2 = time * 1.0E-9
        sumTime2 += (time2 * time2)
      }
    }

    ++halfnEvals
    double rmsTime = sqrt(sumTime2 / halfnEvals)
    logger.info(format("\n Minimum time:           %6.3f (sec)", minTime * 1.0E-9))
    logger.info(format(" RMS time (latter half): %6.3f (sec)", rmsTime))

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (energy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(energy)
    }
    return potentials
  }

}

