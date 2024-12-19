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

import edu.rit.mp.DoubleBuf
import edu.rit.mp.IntegerBuf
import edu.rit.pj.Comm
import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.Crystal
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.EstimateBootstrapper
import ffx.numerics.estimator.SequentialEstimator
import ffx.numerics.math.BootStrapStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.Utilities
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.BARFilter
import ffx.potential.parsers.SystemFilter
import ffx.utilities.FFXScript
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.regex.Matcher
import java.util.regex.Pattern

import static ffx.utilities.Constants.NS2SEC
import static groovy.io.FileType.FILES
import static java.lang.Integer.MIN_VALUE
import static java.lang.Integer.parseInt
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.normalize
import static org.apache.commons.math3.util.FastMath.min

/**
 * This script ...
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt &lt;structures2&gt;
 */
@Command(description = "Performs analysis of for NEQ simulations with work logs in directories ./1-N.", name = "AnalyzeNEQ")
class AnalyzeNEQ extends AlgorithmsScript {

  /**
   * --fDir --forwardDirPath Path to forward work directory.
   */
  @Option(names = ['--fDir', "--forwardDirPath"], paramLabel = ".", defaultValue = ".",
          description = 'Path to forward work directory.')
  String forwardDir // Todo set as required instead of option

  /**
   * --fDir --forwardDirPath Path to reverse work directory.
   */
  @Option(names = ['--rDir', "--reverseDirPath"], paramLabel = ".", defaultValue = ".",
          description = 'Path to reverse work directory.')
  String reverseDir // todo set as required instead of option

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

//  /**
//   * The final argument(s) should be filenames for lambda windows in order.
//   */
//  @Parameters(arity = "1..*", paramLabel = "files",
//      description = 'A single PDB/XYZ when windows are auto-determined (or two for dual topology). Two trajectory files for BAR between two ensembles (or four for dual topology).')
//  List<String> filenames = null

//  /** List of files */
//  List<File> files


  /**
   * BAR Constructor.
   */
    AnalyzeNEQ() {
    this(new Binding())
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
    double temp = 300.0
    BARFilter barFilter = new BARFilter(barFile, new double[fworks.length], fworks, rworks, new double[rworks.length], new double[3], new double[3], temp)
    barFilter.writeFile(outputName, false, false)

    // ffxc BAR --nw 2 --ni 10000 --useTinker end.pdb
    // Create a Binding for command line arguments.
    // String[] args = {"-I", "5", getResourcePath("5awl.pdb")};
    List<String> commandArgs = new ArrayList<>()

    // Options
    commandArgs.add("--nw")
    commandArgs.add("2")
    commandArgs.add("--ni")
    commandArgs.add("10000")
    commandArgs.add("--useTinker")
    commandArgs.add("--tinkerDir")
    commandArgs.add(outputName)

    // Parameters
    commandArgs.add("test.pdb") // /Users/jakemiller/forcefieldx/examples/1EY0.pdb todo fix - coord file

//    ClassLoader classLoader = getClass().getClassLoader()
//    URL url = getClass().location
//    URL url = classLoader.getResource("4icb_ca_a.xyz")
//    commandArgs.add(url.getPath())

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
    for (int i=0; i < numFiles; i++) {
    File file = files.get(i)
    if (!file.exists()) {
      logger.info(format(" Ignoring file that does not exist: %s", file.getAbsolutePath()))
      continue
    }

    String path = normalize(file.getAbsolutePath())
//      logger.info(format(" Current File: %s", path))

    String workLine = ""

    try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
      Pattern re = Pattern.compile(reSearch)
      String line

      while ((line = reader.readLine()) != null) {
        Matcher matcher = re.matcher(line)
        if (matcher.find()) {
//            System.out.println(line)
          workLine = line
        }
      }
    } catch (IOException e) {
      System.err.println("Error reading file: " + e.getMessage())
    }
    String[] workSplit = workLine.split()

    if (workSplit.length != 4) {
      logger.warning(format("%s line is NOT length four: \"%s\"", reSearch, workLine))
      continue
    }
      works.add(workSplit[3].toDouble())
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
