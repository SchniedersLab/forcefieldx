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
package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsScript;
import ffx.numerics.Potential;
import ffx.potential.Utilities;
import ffx.potential.parsers.BARFilter;
import ffx.utilities.FFXScript;
import groovy.lang.Binding;
import groovy.lang.Script;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static ffx.utilities.FileUtils.traverseFiles;
import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.getBaseName;
import static org.apache.commons.io.FilenameUtils.normalize;

/**
 * Performs analysis of Non-equilibrium work simulations. Takes directory paths to the forward and reverse work logs and
 * uses them to create BAR files and calculate the BAR free energy difference.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt; &lt;structures2&gt;
 */
@Command(description = "Performs analysis of Non-equilibrium work simulations.", name = "AnalyzeNEQ")
public class AnalyzeNEQ extends AlgorithmsScript {

  /**
   * --reFile --fileSelectionRegex Locate files that match a Regular expression (.* includes all files).
   */
  @Option(names = {"--reFile", "--fileSelectionRegex"}, paramLabel = "work.log", defaultValue = "work.log",
      description = "Locate files that match a regular expression.")
  private String reFile;

  /**
   * --reSearch --fileSearchRegex Search files for a Regular expression.
   */
  @Option(names = {"--reSearch", "--fileSearchRegex"}, paramLabel = "Boole", defaultValue = "Boole",
      description = "Locate this regular expression in log files.")
  private String reSearch;

  /**
   * --bi --barIterations Maximum iterations for BAR calculation.
   */
  @Option(names = {"--bi", "--barIterations"}, paramLabel = "100", defaultValue = "100",
      description = "Maximum iterations for BAR calculation.")
  private int barIterations;

  /**
   * The final argument(s) should be filenames for lambda windows in order.
   */
  @Parameters(arity = "2", paramLabel = "path",
      description = "Two paths to directories. The first to the forward work directory and second to the reverse work directory.")
  private List<String> directories = null;

  /**
   * AnalyzeNEQ Constructor.
   */
  public AnalyzeNEQ() {
    super();
  }

  /**
   * AnalyzeNEQ Constructor.
   *
   * @param binding The Groovy Binding to use.
   */
  public AnalyzeNEQ(Binding binding) {
    super(binding);
  }

  /**
   * AnalyzeNEQ constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public AnalyzeNEQ(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public AnalyzeNEQ run() {
    // Begin boilerplate code.
    if (!init()) {
      return this;
    }

    int recurse = 1; // todo have as input option

    String forwardDir = directories.get(0);
    String reverseDir = directories.get(1);

    // Collect the forward works.
    File fdir = new File(forwardDir);
    List<File> ffiles = traverseFiles(fdir, recurse, reFile);
    double[] fworks = grabWorks(ffiles);

    // Collect the reverse works.
    File rdir = new File(reverseDir);
    List<File> rfiles = traverseFiles(rdir, recurse, reFile);
    double[] rworks = grabWorks(rfiles);

    // Create BAR file
    String outputName = getBaseName(fdir.getName()) + "-" + getBaseName(rdir.getName()) + ".bar"; // todo could be specified

    File barFile = new File(".", "this.xyz");
    double temp = 300.0; // todo could be specified
    BARFilter barFilter = new BARFilter(barFile, new double[fworks.length], fworks, rworks, new double[rworks.length],
        new double[3], new double[3], temp, temp);
    barFilter.writeFile(outputName, false, false);

    // Create a Binding for command line arguments.
    List<String> commandArgs = new ArrayList<>();

    // Options
    commandArgs.add("--ns");
    commandArgs.add("2");
    commandArgs.add("--ni");
    commandArgs.add(Integer.toString(barIterations));
    commandArgs.add("--useTinker");
    commandArgs.add("--bf");
    commandArgs.add(outputName);

    // Parameters
    // commandArgs.add("test.pdb"); // todo fix - need coord file

    Binding binding = new Binding();
    binding.setVariable("args", commandArgs);

    Class<? extends FFXScript> script;
    script = getScript("BAR");

    // Create a new instance of the script and run it.
    try {
      Script groovyScript = script.getDeclaredConstructor().newInstance();
      groovyScript.setBinding(binding);
      groovyScript.run();
    } catch (Exception e) {
      logger.info(" Exception running BAR.");
      logger.info(e.toString());
      logger.info(Utilities.stackTraceToString(e));
    }

    return this;
  }

  /**
   * Search the files for the regular expression and separate matching lines to get work values.
   *
   * @param files - list of found files
   * @return works - double array with work values from files
   */
  private double[] grabWorks(List<File> files) {
    // Sort the files.
    Collections.sort(files);

    int numFiles = files.size();
    List<Double> works = new ArrayList<>();

    // Compile the regular expression pattern.
    Pattern workRegEx = Pattern.compile(reSearch);

    for (int i = 0; i < numFiles; i++) {
      File file = files.get(i);
      if (!file.exists()) {
        logger.info(format(" Ignoring file that does not exist: %s", file.getAbsolutePath()));
        continue;
      }

      String path = normalize(file.getAbsolutePath());
      try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
        String line;
        // todo above - have ability to request single files with all works or multiple files and grab last work
        while ((line = reader.readLine()) != null) {
          Matcher matcher = workRegEx.matcher(line);
          if (matcher.find()) {
            String[] workSplit = line.trim().split("\\s+");
            if (!(workSplit.length == 4 || workSplit.length == 5)) {
              logger.warning(format("%s line is NOT length four or five: \"%s\"", reSearch, line));
              continue;
            }
            works.add(Double.parseDouble(workSplit[workSplit.length - 1]));
          }
        }

      } catch (IOException e) {
        System.err.println("Error reading file: " + e.getMessage());
      }
    }

    // Convert List<Double> to double[]
    double[] result = new double[works.size()];
    for (int i = 0; i < works.size(); i++) {
      result[i] = works.get(i);
    }
    return result;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public List<Potential> getPotentials() {
    return Collections.emptyList();
  }
}
