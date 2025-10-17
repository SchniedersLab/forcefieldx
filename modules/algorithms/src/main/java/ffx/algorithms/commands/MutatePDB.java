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
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXBinding;
import ffx.utilities.Keyword;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static ffx.utilities.StringUtils.parseAtomRanges;

/**
 * The MutatePDB script mutates residue(s) of a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc MutatePDB [options] &lt;pdb&gt;
 */
@Command(description = " Mutate PDB residue(s).", name = "MutatePDB")
public class MutatePDB extends AlgorithmsCommand {

  /**
   * -r or --resid Residue numbers.
   */
  @Option(names = {"--resid", "-r"}, paramLabel = "1", defaultValue = "1",
      description = "Residue number(s).")
  private String resIDs;

  /**
   * -n or --resname New residue names.
   */
  @Option(names = {"--resname", "-n"}, paramLabel = "ALA", defaultValue = "ALA",
      description = "New residue name(s).")
  private String resNameString;

  /**
   * -ch or --chain Single character chain name (default is ' '). If only one chain exists, that chain will be mutated.
   */
  @Option(names = {"--chain", "--ch"}, paramLabel = " ", defaultValue = " ",
      description = "Single character chain name (default is ' ').")
  private String chainString;

  /**
   * -R or --rotamer Rotamer number to apply.
   */
  @Option(names = {"--rotamer", "-R"}, paramLabel = "-1", defaultValue = "-1",
      description = "Rotamer number to apply.")
  private int rotamer;

  /**
   * --allChains  Mutate all copies of a chain in a multimeric protein.
   */
  @Option(names = {"--allChains"}, paramLabel = "false", defaultValue = "false",
      description = "Mutate all copies of a chains in a multimeric protein.")
  private boolean allChains;

  /**
   * A PDB filename.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "A PDB input file.")
  private String filename;

  private ForceFieldEnergy forceFieldEnergy;

  /**
   * MutatePDB Constructor.
   */
  public MutatePDB() {
    super();
  }

  /**
   * MutatePDB Constructor.
   *
   * @param binding The Binding to use.
   */
  public MutatePDB(FFXBinding binding) {
    super(binding);
  }

  /**
   * MutatePDB constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public MutatePDB(String[] args) {
    super(args);
  }

  /**
   * Execute the script.
   */
  @Override
  public MutatePDB run() {

    if (!init()) {
      return this;
    }

    // The "false" assembly provides access to the chainIDs without compromising the mutated molecular assembly.
    // Used if --allChains is true.

    // Load the MolecularAssembly.
    MolecularAssembly falseAssembly = getActiveAssembly(filename);
    if (falseAssembly == null) {
      logger.info(helpString());
      return this;
    }

    List<Integer> resIDint = parseAtomRanges("residueID", resIDs, falseAssembly.getAtomList().size());

    String[] resNameArr = Arrays.stream(resNameString.split("\\.|,|;")).map(String::trim).toArray(String[]::new);

    String[] chainStringArr = Arrays.stream(chainString.split("\\.|,|;")).map(String::trim).toArray(String[]::new);
    String s = "";
    for (String sub : chainStringArr) {
      if (sub.length() >= 2) {
        logger.warning("Chain ID's have to be single characters separated by commas!");
        return this;
      } else {
        s += sub;
      }
    }
    char[] chainArr = s.toCharArray();

    // For every chain, mutate the residue.
    Polymer[] chains = falseAssembly.getChains();

    if (chains.length == 1 && chainArr.length == 0) {
      chainArr = new char[]{chains[0].getChainID()};
    }

    if (resIDint.size() != resNameArr.length) { // || resIDint.size() != chainArr.length) {
      logger.warning("The number of chains, residue names, and residue ids must be the same!");
      return this;
    }

    int destRotamer = 0;
    if (rotamer > -1) {
      destRotamer = rotamer;
    }

    // Read in command line.
    File structureFile = new File(filename);
    int index = filename.lastIndexOf(".");
    String name = filename.substring(0, index);
    MolecularAssembly molecularAssembly = new MolecularAssembly(name);
    molecularAssembly.setFile(structureFile);

    CompositeConfiguration properties = Keyword.loadProperties(structureFile);
    ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField = forceFieldFilter.parse();
    molecularAssembly.setForceField(forceField);

    PDBFilter pdbFilter = new PDBFilter(structureFile, molecularAssembly, forceField, properties);

    List<PDBFilter.Mutation> mutations = new ArrayList<>();
    for (int i = 0; i <= resIDint.size() - 1; i++) {
      // need to add one to get the correct resID because parseAtomRanges removes one
      if (allChains) {
        for (Polymer currentChain : chains) {
          logger.info("\n Mutating residue number " + (resIDint.get(i) + 1) + " of chain " + currentChain.getChainID() + " to " + resNameArr[i]);
          mutations.add(new PDBFilter.Mutation(resIDint.get(i) + 1, currentChain.getChainID(), resNameArr[i]));
        }
      } else {
        logger.info("\n Mutating residue number " + (resIDint.get(i) + 1) + " of chain " + chainArr[i] + " to " + resNameArr[i]);
        mutations.add(new PDBFilter.Mutation(resIDint.get(i) + 1, chainArr[i], resNameArr[i]));
      }
    }
    pdbFilter.mutate(mutations);

    pdbFilter.readFile();
    pdbFilter.applyAtomProperties();
    molecularAssembly.finalize(true, forceField);

    if (destRotamer > -1) {
      if (allChains) {
        chains = molecularAssembly.getChains();
        for (Polymer currentChain : chains) {
          for (int i = 0; i <= resIDint.size() - 1; i++) {
            Residue residue = currentChain.getResidue(resIDint.get(i) + 1);
            Rotamer[] rotamers = residue.getRotamers();
            if (rotamers != null && rotamers.length > 0) {
              RotamerLibrary.applyRotamer(residue, rotamers[destRotamer]);
            } else {
              logger.info(" No rotamer to apply.");
            }
          }
        }
      } else {
        for (int i = 0; i <= resIDint.size() - 1; i++) {
          Polymer polymer = molecularAssembly.getChain(String.valueOf(chainArr[i]));
          Residue residue = polymer.getResidue(resIDint.get(i) + 1);
          Rotamer[] rotamers = residue.getRotamers();
          if (rotamers != null && rotamers.length > destRotamer) {
            RotamerLibrary.applyRotamer(residue, rotamers[destRotamer]);
          } else {
            logger.info(" No rotamer to apply.");
          }
        }
      }
    }
    pdbFilter.writeFile(structureFile, false);

    for (PDBFilter.Mutation mut : mutations) {
      mut.calculateTorsion();
    }

    forceFieldEnergy = molecularAssembly.getPotentialEnergy();

    return this;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials;
    if (forceFieldEnergy == null) {
      potentials = Collections.emptyList();
    } else {
      potentials = Collections.singletonList(forceFieldEnergy);
    }
    return potentials;
  }

}
