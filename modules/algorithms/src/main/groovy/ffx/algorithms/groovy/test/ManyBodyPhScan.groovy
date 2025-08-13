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
package ffx.algorithms.groovy.test

import edu.rit.pj.Comm
import ffx.potential.cli.PotentialScript
import ffx.utilities.FFXScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Unmatched

import static java.lang.String.format

@Command(description = " Run a pH Scan with ManyBody.", name = "test.ManyBodyPHScan")
class ManyBodyPhScan extends PotentialScript {

  /**
   * --spH --startpH Lower end of the pH range to be evaluated.
   */
  @Option(names = ['--spH', "--startpH"], paramLabel = "0.0", defaultValue = "0.0",
      description = 'Lower end of the pH range to be evaluated.')
  double start = 0.0

  /**
   * --epH --endpH Upper end of the pH range to be evaluated.
   */
  @Option(names = ['--epH', "--endpH"], paramLabel = "14.0", defaultValue = "14.0",
      description = 'Upper end of the pH range to be evaluated.')
  double end = 14.0

  /**
   * --ns --nSteps Number of steps for the given pH range.
   */
  @Option(names = ['--ns', "--nSteps"], paramLabel = "2.0", defaultValue = "2.0",
      description = 'Number of steps for a given pH range.')
  double nSteps = 2.0

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Unmatched
  List<String> unmatched = null

  /**
   * ManyBodyPHScan Constructor.
   */
  ManyBodyPhScan() {
    this(new groovy.lang.Binding())
  }

  /**
   * ManyBodyPHScan Constructor.
   * @param binding The Groovy Binding to use.
   */
  ManyBodyPhScan(Binding binding) {
    super(binding)
  }

  @Override
  ManyBodyPhScan run() {

    if (!init()) {
      return this
    }

    // Set a flag to avoid double use of MPI in downstream commands.
    System.setProperty("pj.use.mpi", "false")

    Class<? extends FFXScript> script = getScript(unmatched.get(0))
    Comm world = Comm.world()
    int numProc = world.size()
    int rank = world.rank()

    if (numProc > 1) {
      logger.info(format(" Number of processes:  %d", numProc))
      logger.info(format(" Rank of this process: %d", rank))
    }

    // Remove ManyBodyPHScan command.
    unmatched.remove(0)

    double stepSize = (end - start) / (nSteps - 1)

    int pHIndex = unmatched.indexOf("0.0")
    for (int i = 0; i < nSteps; i++) {
      double pHValue = start + (stepSize * i)
      List<String> commandArgs = new ArrayList<>()
      for (String arg : unmatched) {
        unmatched.set(pHIndex, pHValue.toString())
        commandArgs.add(arg)
      }

      // Create a Binding for command line arguments.
      Binding binding = new Binding()
      binding.setVariable("args", commandArgs)

      Script groovyScript = script.getDeclaredConstructor().newInstance()
      groovyScript.setBinding(binding)

      try {
        groovyScript.run()
      } catch (Exception e) {
        logger.info(format(" Exception for pH value: %s", pHValue.toString()))
      }
    }
    // Clear the pj.use.mpi flag.
    System.clearProperty("pj.use.mpi")

    return this

  }
}
