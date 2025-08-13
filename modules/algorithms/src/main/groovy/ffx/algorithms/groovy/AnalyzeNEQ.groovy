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

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.Utilities
import ffx.potential.parsers.BARFilter
import ffx.utilities.FFXScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.regex.Matcher
import java.util.regex.Pattern

import static groovy.io.FileType.FILES
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.normalize

/**
 * Performs analysis of Non-equilibrium work simulations. Takes directory paths to the forward and reverse work logs and
 * uses them to create BAR files and calculate the BAR free energy difference.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt &lt;structures2&gt;
 */
@Command(description = "Performs analysis of Non-equilibrium work simulations.", name = "AnalyzeNEQ")
class AnalyzeNEQ extends AlgorithmsScript {

  /**
   * --reFile --fileSelectionRegex Locate files that match a Regular expression (.* includes all files).
   */
  @Option(names = ['--reFile', "--fileSelectionRegex"], paramLabel = "work.log", defaultValue = "work.log",
          description = 'Locate files that match a regular expression.')
  String reFile

  /**
   * --reSearch --fileSearchRegex Search files for a Regular expression.
   */
  @Option(names = ['--reSearch', "--fileSearchRegex"], paramLabel = "Boole", defaultValue = "Boole",
          description = 'Locate this regular expression in log files.')
  String reSearch

  /**
   * --reSearch --fileSearchRegex Maximum iterations for BAR calculation.
   */
  @Option(names = ['--bi', "--barIterations"], paramLabel = "100", defaultValue = "100",
          description = 'Maximum iterations for BAR calculation.')
  int barIterations

  /**
   * The final argument(s) should be filenames for lambda windows in order.
   */
  @Parameters(arity = "2", paramLabel = "path",
          description = 'Two paths to directories. The first to the forward work directory and second to the reverse work directory.')
  List<String> directories = null

  /**
   * BAR Constructor.
   */
  AnalyzeNEQ() {
    this(new groovy.lang.Binding())
  }

  /**
   * BAR Constructor.
   * @param binding The Groovy Binding to use.
   */
  AnalyzeNEQ(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  AnalyzeNEQ run() {
    // Begin boilerplate code.
    if (!init()) {
      return this
    }

    int recurse = 1 // todo have as input option

    String forwardDir = directories[0]
    String reverseDir = directories[1]

    // Collect the forward works.
    File fdir = new File(forwardDir)
    List<File> ffiles = []
    fdir.traverse(type: FILES, maxDepth: recurse, nameFilter: ~/$reFile/) {
      ffiles.add(it)
    }
    double[] fworks = grabWorks(ffiles)

    // Collect the reverse works.
    File rdir = new File(reverseDir)
    List<File> rfiles = []
    rdir.traverse(type: FILES, maxDepth: recurse, nameFilter: ~/$reFile/) {
      rfiles.add(it)
    }
    double[] rworks = grabWorks(rfiles)

    // Create BAR file
    String outputName = fdir.getBaseName() + "-" + rdir.getBaseName() + ".bar" // todo could be specified

    File barFile = new File(".", "this.xyz")
    double temp = 300.0 // todo could be specified
    BARFilter barFilter = new BARFilter(barFile, new double[fworks.length], fworks, rworks, new double[rworks.length],
        new double[3], new double[3], temp, temp)
    barFilter.writeFile(outputName, false, false)

    // Create a Binding for command line arguments.
    List<String> commandArgs = new ArrayList<>()

    // Options
    commandArgs.add("--ns")
    commandArgs.add("2")
    commandArgs.add("--ni")
    commandArgs.add(barIterations as String)
    commandArgs.add("--useTinker")
    commandArgs.add("--bf")
    commandArgs.add(outputName)

    // Parameters
    // commandArgs.add("test.pdb") // todo fix - need coord file

    Binding binding = new Binding()
    binding.setVariable("args", commandArgs)

    Class<? extends FFXScript> script
    script = getScript("BAR")

    // Create a new instance of the script and run it.
    Script groovyScript = script.getDeclaredConstructor().newInstance()
    groovyScript.setBinding(binding)

    try {
      groovyScript.run()
    } catch (Exception e) {
      logger.info(" Exception running BAR.")
      logger.info(e.toString())
      logger.info(Utilities.stackTraceToString(e))
    }

    return this
  }

  /**
   * Search the files for the regular expression and separate matching lines to get work values.
   * @param files - list of found files
   * @param reverseNeg - boolean that controls whether to multiply the reverse works by -1
   * @return works - double array with work values from files
   */
  double[] grabWorks(List<File> files) {
    // Sort the files.
    Collections.sort(files)

    int numFiles = files.size()
    List<Double> works = new ArrayList<Double>();
    for (int i = 0; i < numFiles; i++) {
      File file = files.get(i)
      if (!file.exists()) {
        logger.info(format(" Ignoring file that does not exist: %s", file.getAbsolutePath()))
        continue
      }

      String path = normalize(file.getAbsolutePath())

      String workLine = ""

      try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
        Pattern re = Pattern.compile(reSearch)
        String line

        // todo above - have ability to request single files with all works or multiple files and grab last work

        while ((line = reader.readLine()) != null) {
          Matcher matcher = re.matcher(line)
          if (matcher.find()) {
            String[] workSplit = line.split()
            if (!(workSplit.length == 4 || workSplit.length == 5)) {
              logger.warning(format("%s line is NOT length four or five: \"%s\"", reSearch, workLine))
              continue
            }
            works.add(workSplit[workSplit.length - 1].toDouble())
          }
        }

      } catch (IOException e) {
        System.err.println("Error reading file: " + e.getMessage())
      }
    }

    return works
  }

/**
 * {@inheritDoc}
 */
  @Override
  List<Potential> getPotentials() {
    return Collections.emptyList()
  }
}
