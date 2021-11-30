//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.algorithms.groovy

import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.utilities.FFXScript
import picocli.CommandLine.Command
import picocli.CommandLine.Unmatched

import static groovy.io.FileType.DIRECTORIES
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*
import static org.apache.commons.lang3.math.NumberUtils.isParsable

/**
 * Run an FFX command in a series of sub-directories sequentially or in parallel.
 * <br>
 * Usage:
 * <br>
 * ffxc ForEachDir [command] [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Run an FFX command in a series of sub-directories.", name = "ffxc ForEachDir")
class ForEachDir extends AlgorithmsScript {

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Unmatched
  List<String> unmatched = null

  /**
   * Minimize Constructor.
   */
  ForEachDir() {
    this(new Binding())
  }

  /**
   * Minimize Constructor.
   * @param binding The Groovy Binding to use.
   */
  ForEachDir(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  ForEachDir run() {

    if (!init()) {
      return this
    }

    // Set a flag to avoid double use of MPI in downstream commands.
    System.setProperty("pj.use.mpi", "false")

    Class<? extends FFXScript> script = getScript(unmatched.get(0))
    if (script != null) {
      logger.info(format(" The %s will be run in each subdirectory.", script))
    } else {
      logger.info(format(" %s was not recognized.", unmatched.get(0)))
      return this
    }

    Comm world= Comm.world()
    int numProc = world.size()
    int rank = world.rank()

    if (numProc > 1) {
      logger.info(format(" Number of processes:  %d", numProc))
      logger.info(format(" Rank of this process: %d", rank))
    }

    // Remove the ForEachDir command.
    unmatched.remove(0)

    File cwd = new File(".")
    List<File> directories = []
    cwd.traverse(type: DIRECTORIES, maxDepth: 0) {
      directories.add(it)
    }

    Collections.sort(directories, new Comparator<File>() {
      @Override
      int compare(File o1, File o2) {
        String s1 = o1.getName()
        String s2 = o2.getName()
        if (isParsable(s1) && isParsable(s2)) {
          // Numeric comparison
          return Double.valueOf(s1).compareTo(Double.valueOf(s2))
        }
        // Fall back to String comparison
        return s1.compareTo(s2)
      }
    })

    int numDir = directories.size()

    for (int i=0; i<numDir; i++) {
      if (i % numProc == rank) {
        File dir = directories.get(i)
        if (dir.exists()) {
          String path = normalize(dir.getAbsolutePath())
          logger.info(format(" Current Directory: %s", path))
          List<String> dirParameters = new ArrayList<>()
          // Pass along the unmatched parameters
          for (String arg : unmatched) {
            if (arg.containsIgnoreCase("SUBDIR")) {
              arg = arg.replaceFirst("SUBDIR", "")
              arg = concat(path, getName(arg))
            }
            dirParameters.add(arg)
          }

          // Create a Binding for command line arguments.
          Binding binding = new Binding()
          binding.setVariable("args", dirParameters)

          // Create a new instance of the script and run it.
          Script groovyScript = script.getDeclaredConstructor().newInstance()
          groovyScript.setBinding(binding)
          groovyScript.run()
        }
      }
    }

    // Clear the pj.use.mpi flag.
    System.clearProperty("pj.use.mpi")

    return this
  }

  @Override
  List<Potential> getPotentials() {
    return Collections.emptyList()
  }

}
