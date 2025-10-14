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
package ffx.xray.commands;

import ffx.algorithms.cli.AlgorithmsScript;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.xray.DiffractionData;
import ffx.xray.cli.XrayOptions;
import groovy.lang.Binding;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Parameters;

import java.util.List;

import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * The SaveConformerPDB script saves alternate conformers with proper alt loc labels.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.SaveConformerPDB [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "xray.SaveConformerPDB")
public class SaveConformerPDB extends AlgorithmsScript {

  @Mixin
  private XrayOptions xrayOptions;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames;
  private MolecularAssembly[] molecularAssemblies;
  private DiffractionData diffractionData;

  /**
   * SaveConformerPDB constructor.
   */
  public SaveConformerPDB() {
    super();
  }

  /**
   * SaveConformerPDB constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public SaveConformerPDB(String[] args) {
    super(args);
  }

  /**
   * SaveConformerPDB constructor.
   * @param binding The Groovy Binding to use.
   */
  public SaveConformerPDB(Binding binding) {
    super(binding);
  }

  @Override
  public SaveConformerPDB run() {
    if (!init()) {
      return this;
    }

    xrayOptions.init();

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0));
      activeAssembly = molecularAssemblies[0];
      filename = filenames.get(0);
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
      filename = activeAssembly.getFile().getAbsolutePath();
    }

    if (molecularAssemblies.length == 1) {
      logger.info("No alternate conformers");
      return this;
    }

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties();
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);
    diffractionData.scaleBulkFit();
    diffractionData.printStats();
    algorithmFunctions.energy(molecularAssemblies);

    List<Residue> residuesA = molecularAssemblies[0].getResidueList();
    List<Residue> residuesB = molecularAssemblies[1].getResidueList();

    for (int j = 0; j < residuesA.size(); j++) {
      List<Atom> atomsA = residuesA.get(j).getAtomList();
      List<Atom> atomsB = residuesB.get(j).getAtomList();
      for (int i = 0; i < atomsA.size(); i++) {
        Atom atomA = atomsA.get(i);
        Atom atomB = atomsB.get(i);

        if (atomA.getAltLoc() == null || atomA.getAltLoc() == ' ') {
          atomA.setAltLoc('A');
        }
        if (atomB.getAltLoc() == null || atomB.getAltLoc() == ' ') {
          atomB.setAltLoc('B');
        }
      }
    }

    logger.info(" ");
    diffractionData.writeModel(removeExtension(filename) + ".pdb");
    diffractionData.writeData(removeExtension(filename) + ".mtz");

    return this;
  }
}
