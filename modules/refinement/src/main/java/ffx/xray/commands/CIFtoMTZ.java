//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.xray.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionRefinementData;
import ffx.xray.parsers.CIFFilter;
import ffx.xray.parsers.MTZWriter;
import ffx.xray.parsers.MTZWriter.MTZType;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.List;

import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * The CIF2MTZ script saves a CIF file to MTZ format.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.CIF2MTZ &lt;filename&gt;
 */
@Command(description = " Convert a CIF file to MTZ format.", name = "xray.CIFtoMTZ")
public class CIFtoMTZ extends AlgorithmsCommand {

  /**
   * A CIF filename.
   */
  @Parameters(arity = "2", paramLabel = "file", description = "A PDB file and a CIF diffraction file.")
  private List<String> filenames = null;

  private MolecularAssembly[] molecularAssemblies;
  private DiffractionRefinementData refinementData;

  /**
   * CIF2MTZ constructor.
   */
  public CIFtoMTZ() {
    super();
  }

  /**
   * CIF2MTZ constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public CIFtoMTZ(String[] args) {
    super(args);
  }

  /**
   * CIF2MTZ constructor.
   *
   * @param binding The Binding to use.
   */
  public CIFtoMTZ(FFXBinding binding) {
    super(binding);
  }

  /**
   * Execute the script.
   */
  @Override
  public CIFtoMTZ run() {

    if (!init()) {
      return this;
    }

    String pdb = filenames.get(0);
    String cif = filenames.get(1);

    logger.info("\n Running CIF2MTZ on " + cif);
    molecularAssemblies = algorithmFunctions.openAll(pdb);

    CIFFilter cifFilter = new CIFFilter();
    ReflectionList reflectionlist =
        cifFilter.getReflectionList(new File(cif), molecularAssemblies[0].getProperties());

    if (reflectionlist == null) {
      logger.info(" Using crystal information from the PDB file to generate MTZ file.");

      Crystal crystal = molecularAssemblies[0].getCrystal().getUnitCell();
      double res = cifFilter.getResolution(new File(cif), crystal);
      if (res < 0.0) {
        logger.info(" The resolution could not be determined from the PDB and CIF files.");
        return this;
      }

      Resolution resolution = new Resolution(res);
      reflectionlist = new ReflectionList(crystal, resolution,
          molecularAssemblies[0].getProperties());
    }

    refinementData = new DiffractionRefinementData(molecularAssemblies[0].getProperties(), reflectionlist);
    cifFilter.readFile(new File(cif), reflectionlist, refinementData, molecularAssemblies[0].getProperties());

    MTZWriter mtzwriter = new MTZWriter(reflectionlist, refinementData,
        removeExtension(cif) + ".mtz", MTZType.DATAONLY);
    mtzwriter.write();

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }
}
