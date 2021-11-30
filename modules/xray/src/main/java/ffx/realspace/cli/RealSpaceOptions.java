// ******************************************************************************
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
// ******************************************************************************
package ffx.realspace.cli;

import ffx.algorithms.AlgorithmFunctions;
import ffx.potential.MolecularAssembly;
import ffx.realspace.RealSpaceData;
import ffx.realspace.parsers.RealSpaceFile;
import ffx.xray.RefinementEnergy;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.cli.DataRefinementOptions;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize a real-space density map.
 *
 * @author Michael J. Schnieders
 * @author Soham Ali
 * @since 1.0
 */
public class RealSpaceOptions extends DataRefinementOptions {

  private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

  /** The refinement mode to use. */
  public RefinementMode refinementMode = RefinementMode.COORDINATES;

  /**
   * -X or --data Specify input data filename, weight applied to the data (wA) and if the data is
   * from a neutron experiment.
   */
  @Option(
      names = {"-X", "--data"},
      arity = "2",
      description = "Specify input data filename and its weight (wA) (e.g. -X filename 1.0).")
  String[] data = null;

  /**
   * Process input to collect Real Space Files.
   *
   * @param filenames Input filenames (first filename is ignored).
   * @param molecularAssembly Currently open molecularAssembly.
   * @return a list of Real Space File instances.
   */
  public List<RealSpaceFile> processData(
      List<String> filenames, MolecularAssembly molecularAssembly) {

    logger.info("\n");

    // set up real space map data (can be multiple files)
    List<RealSpaceFile> mapfiles = new ArrayList<>();
    if (filenames.size() > 1) {
      RealSpaceFile realspacefile = new RealSpaceFile(filenames.get(1), wA);
      mapfiles.add(realspacefile);
    }
    if (data != null) {
      for (int i = 0; i < data.length; i += 2) {
        double w = wA;
        if (data.length > i + 1) {
          try {
            w = Double.parseDouble(data[i + 1]);
          } catch (Exception e) {
            //
          }
        }
        RealSpaceFile realspacefile = new RealSpaceFile(data[i], w);
        mapfiles.add(realspacefile);
      }
    }

    if (mapfiles.size() == 0) {
      RealSpaceFile realspacefile = new RealSpaceFile(molecularAssembly, wA);
      mapfiles.add(realspacefile);
    }

    // TODO - write out map files if needed.
    //        DiffractionFile diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
    //        DiffractionData diffractionData = new DiffractionData(molecularAssembly,
    // molecularAssembly[0].getProperties(),
    //                SolventModel.POLYNOMIAL, diffractionFile);
    //        diffractionData.scaleBulkFit();
    //        diffractionData.printStats();
    //        String mapFileName = String.format("%s_ffx_%d",
    // FilenameUtils.removeExtension(dataFileName), ++nDiffractionData);
    //        diffractionData.writeMaps(mapFileName);
    //        mapFiles.add(new RealSpaceFile(mapFileName + "_2fofc.map", 1.0));

    return mapfiles;
  }

  /**
   * Process input from opened molecular assemblies to a RefinementEnergy
   *
   * @param filenames All filenames included in the real-space data.
   * @param molecularAssemblies Array of MolecularAssembly instances.
   * @param algorithmFunctions An AlgorithmFunctions object.
   * @return An assembled RefinementEnergy with real-space energy.
   */
  public RefinementEnergy toRealSpaceEnergy(
      List<String> filenames,
      MolecularAssembly[] molecularAssemblies,
      AlgorithmFunctions algorithmFunctions) {
    RealSpaceFile[] mapFiles =
        processData(filenames, molecularAssemblies[0]).toArray(new RealSpaceFile[0]);
    RealSpaceData realspacedata =
        new RealSpaceData(
            molecularAssemblies,
            molecularAssemblies[0].getProperties(),
            molecularAssemblies[0].getParallelTeam(),
            mapFiles);
    algorithmFunctions.energy(molecularAssemblies[0]);
    return new RefinementEnergy(realspacedata, RefinementMode.COORDINATES);
  }
}
