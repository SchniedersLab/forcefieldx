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
package ffx.algorithms.groovy

import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.utilities.FFXScript
import picocli.CommandLine.Command
import picocli.CommandLine.Unmatched

import static groovy.io.FileType.FILES
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*

/**
 * Run an FFX command in on a series of files.
 * <br>
 * Usage:
 * <br>
 * ffxc ForEachFile [command] [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Run an FFX command on a series of files in the current directory.", name = "ffxc ForEachFile")
class ForEachFile extends AlgorithmsScript {

  /**
   * --ext or --filetype Only apply the command to files with given filetype (ARC, XYZ, PDB, etc).
   */
  // @Option(names = ['--ext', '--filetype'], paramLabel = 'ext',
  //    description = 'Only apply the command to files with given filetype (ARC, XYZ, PDB, etc).')
  // String ext = null

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Unmatched
  List<String> unmatched = null

  /**
   * Minimize Constructor.
   */
  ForEachFile() {
    this(new Binding())
  }

  /**
   * Minimize Constructor.
   * @param binding The Groovy Binding to use.
   */
  ForEachFile(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  ForEachFile run() {

    if (!init()) {
      return this
    }

    // Set a flag to avoid double use of MPI in downstream commands.
    System.setProperty("pj.use.mpi", "false")

    Class<? extends FFXScript> script = getScript(unmatched.get(0))
    if (script != null) {
      logger.info(format(" The %s will be run on each file.", script))
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

    // Remove the ForEachFile command.
    unmatched.remove(0)

    File cwd = new File(".")
    List<File> files = []
    cwd.traverse(type: FILES, maxDepth: 0) {
      files.add(it)
    }

    int numFiles = files.size()
    for (int i=0; i<numFiles; i++) {
      if (i % numProc == rank) {
        File file = files.get(i)
        if (file.exists()) {
          String path = normalize(file.getAbsolutePath())
          logger.info(format(" Current File: %s", path))
          List<String> dirParameters = new ArrayList<>()
          // Pass along the unmatched parameters
          for (String arg : unmatched) {
//            if (arg.containsIgnoreCase("SUBDIR")) {
//              arg = arg.replaceFirst("SUBDIR", "")
//              arg = concat(path, getName(arg))
//            }
            dirParameters.add(arg)
          }
          // Add the current file.
          dirParameters.add(path)

          // Create a Binding for command line arguments.
          Binding binding = new Binding()
          binding.setVariable("args", dirParameters)

          // Create a new instance of the script and run it.
          Script groovyScript = script.getDeclaredConstructor().newInstance()
          groovyScript.setBinding(binding)

          try {
            groovyScript.run()
          } catch (Exception e) {
            logger.info(format(" Exception for file: %s", path))
          }
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
