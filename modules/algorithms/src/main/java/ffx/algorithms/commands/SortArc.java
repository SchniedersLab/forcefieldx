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
package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.String.format;

/**
 * The SortArc script sort Monte Carlo archive files by lambda value. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc SortArc [options] &lt;structures1&gt; &lt;structures2&gt;
 */

@Command(description = " Unwind .ARC files for nWindows", name = "SortArc")
public class SortArc extends AlgorithmsCommand {

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Mixin
  private TopologyOptions topologyOptions;

  @Option(names = {"--nw", "--nWindows"}, paramLabel = "-1",
      description = "If set, auto-determine lambda values and subdirectories (overrides other flags).")
  private int nWindows = -1;

  @Option(names = {"--bT", "--sortByTemp"}, paramLabel = "false",
      description = "If set, sort archive files by temperature values")
  private boolean sortTemp = false;

  @Option(names = {"--sT", "--startTemp"}, paramLabel = "298.15",
      defaultValue = "298.15",
      description = "Sets the starting temperature for the exponential temperature ladder if sorting by temperature.")
  private double lowTemperature = 298.15;

  @Option(names = {"--ex", "--exponent"}, paramLabel = "0.5",
      defaultValue = "0.5",
      description = "Sets the exponent for the exponential temperature ladder if sorting by temperature.")
  private double exponent = 0.05;

  /**
   * The final argument(s) should be filenames for lambda windows in order..
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "Trajectory files for the first end of the window, followed by trajectories for the other end")
  List<String> filenames = null;

  private double[] lambdaValues;
  private double[] temperatureValues;
  private SystemFilter[] openers;
  private SystemFilter[][] writers;
  private String[] files;
  private CompositeConfiguration additionalProperties;
  private List<String> windowFiles = new ArrayList<>();
  MolecularAssembly[] topologies;
  MolecularAssembly ma;

  /**
   * Sets an optional Configuration with additional properties.
   *
   * @param additionalProps
   */
  public void setProperties(CompositeConfiguration additionalProps) {
    this.additionalProperties = additionalProps;
  }

  /**
   * SortArc Constructor.
   */
  public SortArc() {
    super();
  }

  /**
   * SortArc Constructor.
   *
   * @param binding The Binding to use.
   */
  public SortArc(FFXBinding binding) {
    super(binding);
  }

  /**
   * SortArc constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public SortArc(String[] args) {
    super(args);
  }

  @Override
  public SortArc run() {
    if (!init()) {
      return this;
    }

    // Determine the number of topologies to be read and allocate the array.
    int numTopologies = topologyOptions.getNumberOfTopologies(filenames);
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];

    // Turn on computation of lambda derivatives if softcore atoms exist.
    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    files = new String[numTopologies];
    for (int i = 0; i < numTopologies; i++) {
      files[i] = filenames.get(i);
    }

    if (nWindows != -1) {
      for (int i = 0; i < nWindows; i++) {
        for (int j = 0; j < numTopologies; j++) {
          String fullPathToFile = FilenameUtils.getFullPath(files[j]);
          String directoryFullPath = fullPathToFile.replace(files[j], "") + i;
          windowFiles.add(directoryFullPath + File.separator + i);
        }

      }

      lambdaValues = new double[nWindows];
      temperatureValues = new double[nWindows];
      for (int i = 0; i < nWindows; i++) {
        if (sortTemp) {
          temperatureValues[i] = lowTemperature * Math.exp(exponent * i);
        } else {
          lambdaValues[i] = alchemicalOptions.getInitialLambda(nWindows, i, false);
        }
      }
    }

    if (filenames == null) {
      return this;
    }

    String[][] archiveFullPaths = new String[nWindows][numTopologies];
    File file = new File(files[0]);
    String directoryPath = file.getAbsoluteFile().getParent() + File.separator;
    String[][] archiveNewPath = new String[nWindows][numTopologies];
    File[][] saveFile = new File[nWindows][numTopologies];
    File[][] arcFiles = new File[nWindows][numTopologies];

    for (int j = 0; j < numTopologies; j++) {
      String archiveName = FilenameUtils.getBaseName(files[j]) + ".arc";
      for (int i = 0; i < nWindows; i++) {
        archiveFullPaths[i][j] = directoryPath + i + File.separator + archiveName;
        File arcFile = new File(archiveFullPaths[i][j]);
        arcFiles[i][j] = arcFile;
        archiveNewPath[i][j] = directoryPath + i + File.separator + FilenameUtils.getBaseName(files[j]) + "_E" + i + ".arc";
        saveFile[i][j] = new File(archiveNewPath[i][j]);
      }
    }

    openers = new XYZFilter[numTopologies];
    writers = new XYZFilter[nWindows][numTopologies];

    for (int j = 0; j < numTopologies; j++) {
      if (filenames.get(j).contains(".pdb")) {
        ma = alchemicalOptions.openFile(algorithmFunctions, topologyOptions,
            threadsPerTopology, archiveFullPaths[0][j], j);
      } else {
        ma = alchemicalOptions.openFile(algorithmFunctions, topologyOptions,
            threadsPerTopology, filenames.get(j), j);
      }
      topologies[j] = ma;
      openers[j] = algorithmFunctions.getFilter();

      for (int i = 0; i < nWindows; i++) {
        File arc = saveFile[i][j];
        writers[i][j] = new XYZFilter(arc, topologies[j], topologies[j].getForceField(), additionalProperties);
      }
    }

    double tolerance;
    if (sortTemp) {
      tolerance = 1.0e-2;
    } else {
      tolerance = 1.0e-4;
    }

    for (int j = 0; j < numTopologies; j++) {
      for (int i = 0; i < nWindows; i++) {
        logger.info(format(" Initializing %d topologies for each end", numTopologies));
        openers[j].setFile(arcFiles[i][j]);
        topologies[j].setFile(arcFiles[i][j]);
        logger.info("Set file to:" + arcFiles[i][j].toString());

        int snapshots = openers[j].countNumModels();
        logger.info(String.valueOf(snapshots));

        for (int n = 0; n < snapshots; n++) {
          boolean resetPosition = (n == 0);
          openers[j].readNext(resetPosition, false);
          String remarkLine = openers[j].getRemarkLines()[0];

          double lambda = 0;
          double temp = 0;
          if (remarkLine.contains(" Lambda: ")) {
            String[] tokens = remarkLine.split(" +");
            for (int p = 0; p < tokens.length; p++) {
              if (tokens[p].startsWith("Lambda")) {
                lambda = Double.parseDouble(tokens[p + 1]);
              }
              if (tokens[p].startsWith("Temp")) {
                temp = Double.parseDouble(tokens[p + 1]);
              }
            }

          }

          double diff;
          for (int k = 0; k < nWindows; k++) {
            if (sortTemp) {
              diff = Math.abs(temperatureValues[k] - temp);
            } else {
              diff = Math.abs(lambdaValues[k] - lambda);
            }

            if (diff < tolerance) {
              writers[k][j].writeFile(saveFile[k][j], true, new String[]{remarkLine});
              //set topology back to archive being read in
              topologies[j].setFile(arcFiles[i][j]);
              break;
            }
          }
        }
      }
    }
    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return Collections.emptyList();
  }
}
