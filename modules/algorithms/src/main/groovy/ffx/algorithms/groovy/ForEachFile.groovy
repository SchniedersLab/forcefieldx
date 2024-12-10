//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import edu.rit.pj.*
import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.Utilities
import ffx.utilities.FFXScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Unmatched

import static groovy.io.FileType.FILES
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.normalize

/**
 * Run an FFX command on a series of files. Parallel Java across nodes is supported.
 *
 * Recursion through the directory structure is supported to a supplied level using the
 * --recurse flag (0 only includes the current directory).
 *
 * Files can be selected using a regular expression -- the default matches all
 * files (".*").
 *
 * To control placement of the variable file on the command line, the "FILE" string can be used. If
 * present, it is replaced with the current file. If absent, the current file is last argument
 * to the FFX command.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc ForEachFile [ForEachFile options] Command [Command options] &lt;FILE&gt; [&lt;FILE2&gt;]
 */
@Command(description = " Run an FFX command on a series of files.", name = "ForEachFile")
class ForEachFile extends AlgorithmsScript {

  /**
   * --recurse Maximum recursion level (0 only includes the current directory).
   */
  @Option(names = ['--recurse'], defaultValue = "0", paramLabel = "0",
      description = 'Maximum recursion level (0 only includes the current directory).')
  int recurse

  /**
   * --regex --fileSelectionRegex Evaluate files that match a Regular expression (.* includes all files).
   */
  @Option(names = ['--regex', "--fileSelectionRegex"], paramLabel = ".*", defaultValue = ".*",
      description = 'Locate files that match a regular expression (\'.*\' matches all files).')
  String regex

  /**
   * --regex2 --fileSelectionRegex2 Evaluate files that match a Regular expression (.* includes all files).
   * The second regular expression can be used to place a second file on each command line (e.g., for dual-topology simulations).
   */
  @Option(names = ['--regex2', "--fileSelectionRegex2"], paramLabel = "", defaultValue = "",
      description = 'Locate files that match a 2nd regular expression (\'.*\' matches all files).')
  String regex2

  /**
   * --schedule Load balancing will use a [Dynamic, Fixed, or Guided] schedule.
   */
  @Option(names = ['--schedule'], defaultValue = "dynamic", paramLabel = "dynamic",
      description = 'Load balancing will use a [Dynamic, Fixed, or Guided] schedule.')
  String schedule

  /**
   * -v --verbose Decide whether to print additional logging.
   */
  @Option(names = ['-v', '--verbose'], defaultValue = "false", paramLabel = "false",
      description = 'Print additional logging for errors.')
  boolean verbose

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Unmatched
  List<String> unmatched = null

  /**
   * FFX Script to run in each process.
   */
  Class<? extends FFXScript> script

  /**
   * List of files.
   */
  List<File> files

  /**
   * Place two files on each command line.
   */
  boolean twoFilesPerCommand = false

  /**
   * List of second files to place on the command line.
   */
  List<File> files2

  /**
   * Parallel Java Schedule.
   */
  IntegerSchedule integerSchedule

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

    script = getScript(unmatched.get(0))
    if (script != null) {
      logger.info(format(" The %s will be run on each file.", script))
    } else {
      logger.info(format(" %s was not recognized.", unmatched.get(0)))
      return this
    }

    Comm world = Comm.world()
    int numProc = world.size()
    int rank = world.rank()

    if (numProc > 1) {
      logger.info(format(" Number of processes:  %d", numProc))
      logger.info(format(" Rank of this process: %d", rank))
    }

    // Remove the ForEachFile command.
    unmatched.remove(0)

    // Collect the files.
    File cwd = new File(".")
    files = []
    cwd.traverse(type: FILES, maxDepth: recurse, nameFilter: ~/$regex/) {
      files.add(it)
    }

    if (!regex2.isEmpty()) {
      twoFilesPerCommand = true
      files2 = []
      cwd.traverse(type: FILES, maxDepth: recurse, nameFilter: ~/$regex2/) {
        files2.add(it)
      }
    }

    // Sort the files.
    Collections.sort(files)
    if (twoFilesPerCommand) {
      if (files.size() != files2.size()) {
        logger.info(" The number of files matched by the two regular expressions do not agree.")
        logger.info(" The number of files matched by the first regular expression: " + files.size())
        for (File file : files) {
          logger.info("  File: " + file.getAbsolutePath())
        }
        logger.info(" The number of files matched by the second regular expression: " + files2.size())
        for (File file : files2) {
          logger.info("  File: " + file.getAbsolutePath())
        }
        return this
      }
      Collections.sort(files2)
    }

    // Create the Parallel Java execution Schedule.
    try {
      integerSchedule = IntegerSchedule.parse(schedule.toLowerCase())
      logger.info(" Parallel Schedule: " + schedule)
    } catch (Exception e) {
      integerSchedule = IntegerSchedule.dynamic()
      logger.info(" Parallel Schedule: Dynamic")
    }

    // Create a HybridTeam and then execute the ForEachFileRegion
    WorkerTeam workerTeam = new WorkerTeam(world)
    workerTeam.execute(new ForEachFileRegion())

    // Clear the pj.use.mpi flag.
    System.clearProperty("pj.use.mpi")

    return this
  }

  /**
   * ForEachFileRegion delegates work to ForEachFileLoop instances in Parallel Java processes.
   */
  private class ForEachFileRegion extends WorkerRegion {

    @Override
    void run() throws Exception {
      int numFiles = files.size()
      execute(0, numFiles - 1, new ForEachFileLoop())
    }

  }

  /**
   * ForEachFileLoop is executed in a Parallel Java process.
   */
  private class ForEachFileLoop extends WorkerIntegerForLoop {


    @Override
    IntegerSchedule schedule() {
      return integerSchedule
    }

    @Override
    void run(int lb, int ub) throws Exception {
      for (int i = lb; i <= ub; i++) {
        File file = files.get(i)
        if (!file.exists()) {
          logger.info(format(" Ignoring file that does not exist: %s", file.getAbsolutePath()))
          continue
        }
        File dualTopologyFile = null
        if (twoFilesPerCommand) {
          dualTopologyFile = files2.get(i)
          if (!dualTopologyFile.exists()) {
            logger.info(format(" Ignoring dual topology file that does not exist: %s", dualTopologyFile.getAbsolutePath()))
            continue
          }
        }

        String path = normalize(file.getAbsolutePath())
        logger.info(format(" Current File: %s", path))

        String dualTopologyPath = null
        if (twoFilesPerCommand) {
          dualTopologyPath = normalize(dualTopologyFile.getAbsolutePath())
          logger.info(format(" Current Dual Topology File: %s", dualTopologyPath))
        }

        List<String> commandArgs = new ArrayList<>()
        // Pass along the unmatched parameters
        boolean foundFile = false
        boolean found2File = false
        for (String arg : unmatched) {
          if (arg.equalsIgnoreCase("FILE")) {
            // Replace FILE with the current file path.
            commandArgs.add(path)
            foundFile = true
          } else if (twoFilesPerCommand && arg.equalsIgnoreCase("FILE2")) {
            commandArgs.add(dualTopologyPath)
            found2File = true
          } else {
            commandArgs.add(arg)
          }
        }

        if (!foundFile) {
          // Add the current file as the final argument.
          commandArgs.add(path)
        }
        if (twoFilesPerCommand && !found2File) {
          // Add the dual topology file as the final argument to a dual topology simulation.
          commandArgs.add(dualTopologyPath)
        }

        // Create a Binding for command line arguments.
        Binding binding = new Binding()
        binding.setVariable("args", commandArgs)

        // Create a new instance of the script and run it.
        Script groovyScript = script.getDeclaredConstructor().newInstance()
        groovyScript.setBinding(binding)

        try {
          groovyScript.run()
        } catch (Exception e) {
          logger.info(format(" Exception for file: %s", path))
          if (twoFilesPerCommand) {
            logger.info(format(" Dual topology file: %s", dualTopologyPath))
          }
          if (verbose) {
            logger.info(e.toString());
            logger.info(Utilities.stackTraceToString(e));
          }
        }
      }
    }
  }

  @Override
  List<Potential> getPotentials() {
    return Collections.emptyList()
  }

}
