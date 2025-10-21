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
package ffx.potential.commands;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.parameters.BioType;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * Print out Biotype records for the atoms in an XYZ file.
 *
 * Usage:
 *   ffxc Biotype &lt;filename&gt;
 */
@Command(name = "Biotype", description = " Print out Biotype records for the atoms in an XYZ file.")
public class Biotype extends PotentialCommand {

  @Option(names = {"--name", "--moleculeName"}, paramLabel = "MOL", defaultValue = "MOL",
      description = "The molecule name to use for the Biotype records.")
  private String molName;

  @Option(names = {"-a", "--useAtomNames"}, paramLabel = "false", defaultValue = "false",
      description = "Use the atom names in the XYZ file.")
  private boolean useAtomNames;

  @Option(names = {"-c", "--writeCONECT"}, paramLabel = "false", defaultValue = "false",
      description = "Write out CONECT records to append to the PDB file.")
  private boolean writeCONNECT;

  @Option(names = {"-w", "--writePDB"}, paramLabel = "false", defaultValue = "false",
      description = "Write out a PDB file with the updated atom and molecule names.")
  private boolean writePDB;

  @Parameters(arity = "1", paramLabel = "files",
      description = "An XYZ coordinate file.")
  private String filename = null;

  /** List of created Biotype records. */
  private List<BioType> bioTypes = null;

  public Biotype() {
    super();
  }

  public Biotype(FFXBinding binding) {
    super(binding);
  }

  public Biotype(String[] args) {
    super(args);
  }

  /** Get the list of created Biotype records. */
  public List<BioType> getBioTypes() {
    return bioTypes;
  }

  @Override
  public Biotype run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    logger.info("\n Running Biotype on " + filename);

    Molecule[] molecules = activeAssembly.getMoleculeArray();
    Atom[] atoms = activeAssembly.getAtomArray();

    if (writeCONNECT) {
      for (Atom a : atoms) {
        // PDB CONECT format columns reference (7-31 are bond target serials)
        StringBuilder sb = new StringBuilder(String.format("CONECT%5s", Integer.toString(a.getXyzIndex())));
        List<Bond> bonds = a.getBonds();
        for (Bond b : bonds) {
          sb.append(String.format("%5s", Integer.toString(b.get1_2(a).getXyzIndex())));
        }
        logger.info(sb.toString());
      }
    }

    if (molecules.length > 1) {
      logger.info(" Biotype is intended for a system with one molecule.");
      return this;
    }

    // Update the molecule name.
    molecules[0].setName(molName);

    // Create a List of biotype entries.
    bioTypes = new ArrayList<>();

    // Update atom names if requested.
    if (!useAtomNames) {
      Atom.ElementSymbol[] symbols = Atom.ElementSymbol.values();
      int[] elementCounts = new int[symbols.length + 1]; // atomic numbers start at 1
      for (Atom atom : atoms) {
        int element = atom.getAtomicNumber();
        int n = elementCounts[element];
        String name = symbols[element - 1].name().toUpperCase() + n;
        atom.setName(name);
        elementCounts[element] = n + 1;
      }
    }

    int index = 1;
    for (Atom atom : atoms) {
      // Update the residue/molecule name.
      atom.setResName(molName);

      // Collect bond neighbor atom names.
      List<Bond> bonds = atom.getBonds();
      String[] bondString = null;
      if (bonds != null) {
        bondString = new String[bonds.size()];
        int i = 0;
        for (Bond bond : bonds) {
          bondString[i++] = bond.get1_2(atom).getName();
        }
      }

      // Create a Biotype entry and log it.
      BioType biotype = new BioType(index++, atom.getName(), molName,
          atom.getAtomType().type, bondString);

      bioTypes.add(biotype);
      logger.info(biotype.toString());
    }

    // Optionally, save a PDB with updated atom/molecule names.
    if (writePDB) {
      File pdbFile = createOutputFile(filename, "pdb");
      logger.info("\n Saving PDB file: " + pdbFile);
      potentialFunctions.saveAsPDB(activeAssembly, pdbFile);
    }

    // Return the bioTypes via the Binding for compatibility with callers.
    binding.setVariable("bioTypes", bioTypes);

    return this;
  }
}
