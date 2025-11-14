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
package ffx.realspace.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.ManyBodyOptions;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.algorithms.optimize.TitrationManyBody;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.realspace.cli.RealSpaceOptions;
import ffx.utilities.FFXBinding;
import ffx.xray.RefinementEnergy;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard;
import static java.lang.String.format;

/**
 * Run GenZ for ensemble generation.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.GenZ [options] &lt;filename&gt;
 */
@Command(description = " Run GenZ for ensemble generation.", name = "ffxc realspace.GenZ")
public class GenZ extends AlgorithmsCommand {

  @Mixin
  private ManyBodyOptions manyBodyOptions;

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Mixin
  private RealSpaceOptions realSpaceOptions;

  @Option(names = {"--rEE", "--ro-ensembleEnergy"}, paramLabel = "0.0",
      description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
  private String ensembleEnergy = "0.0";

  @Option(names = {"--pF", "--printFiles"}, paramLabel = "false",
      description = "Write to an energy restart file and ensemble file.")
  private boolean printFiles = false;

  @Option(names = {"--pKa"}, paramLabel = "false",
      description = "Calculating protonation populations for pKa shift.")
  private boolean pKa = false;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames;

  private RefinementEnergy refinementEnergy;

  ForceFieldEnergy potentialEnergy;
  private MolecularAssembly[] molecularAssemblies;
  private ForceField forceField;
  TitrationManyBody titrationManyBody;
  List<Residue> selectedResidues;
  MolecularAssembly[] conformerAssemblies = new MolecularAssembly[3];

  /**
   * GenZ Constructor.
   */
  public GenZ() {
    super();
  }

  /**
   * GenZ constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public GenZ(String[] args) {
    super(args);
  }

  /**
   * GenZ Constructor.
   * @param binding The Binding to use.
   */
  public GenZ(FFXBinding binding) {
    super(binding);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public GenZ run() {
    if (!init()) {
      return this;
    }

    // Atomic clashes are expected and will be handled using direct induced dipoles.
    System.setProperty("sor-scf-fallback", "false");
    System.setProperty("direct-scf-fallback", "true");

    double titrationPH = manyBodyOptions.getTitrationPH();
    double inclusionCutoff = manyBodyOptions.getInclusionCutoff();
    int mutatingResidue = manyBodyOptions.getInterestedResidue();
    boolean onlyTitration = manyBodyOptions.getOnlyTitration();
    double pHRestraint = manyBodyOptions.getPHRestraint();
    if (manyBodyOptions.getTitration()) {
      System.setProperty("manybody-titration", "true");
    }

    // Many-Body expansion of the X-ray target converges much more quickly with the NEA.
    String nea = System.getProperty("native-environment-approximation", "true");
    System.setProperty("native-environment-approximation", nea);

    boolean lambdaTerm = alchemicalOptions.hasSoftcore();
    if (lambdaTerm) {
      // Turn on softcore van der Waals
      System.setProperty("lambdaterm", "true");
      // Turn of alchemical electrostatics
      System.setProperty("elec-lambdaterm", "false");
      // Turn on intra-molecular softcore
      System.setProperty("intramolecular-softcore", "true");
    }
    System.setProperty("ro-ensembleEnergy", ensembleEnergy);

    String modelFilename;
    if (filenames != null && !filenames.isEmpty()) {
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0));
      activeAssembly = molecularAssemblies[0];
      modelFilename = filenames.get(0);
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
      modelFilename = activeAssembly.getFile().getAbsolutePath();
    }

    CompositeConfiguration properties = activeAssembly.getProperties();
    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);
    potentialEnergy = activeAssembly.getPotentialEnergy();
    forceField = activeAssembly.getForceField();
    if (forceField == null) {
      logger.info("This force field is null");
    }

    String filename = filenames.get(0);

    List<Residue> titrateResidues = new ArrayList<>();

    //Prepare variables for saving out the highest population rotamers (optimal rotamers)
    Set<Atom> excludeAtoms = new HashSet<>();
    boolean isTitrating = false;

    String[] titratableResidues = new String[]{"HIS", "HIE", "HID", "GLU", "GLH", "ASP", "ASH", "LYS", "LYD", "CYS", "CYD"};
    List<String> titratableResiudesList = Arrays.asList(titratableResidues);
    double boltzmannWeights;
    double offsets;
    double[][] populationArray = new double[1][55];
    double[][] titrateBoltzmann;
    double[] protonationBoltzmannSums;
    double totalBoltzmann = 0;
    List<Residue> residueList = activeAssembly.getResidueList();

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly);
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.");
      return this;
    }

    // Application of rotamers uses side-chain atom naming from the PDB.
    if (properties.getBoolean("standardizeAtomNames", false)) {
      renameAtomsToPDBStandard(activeAssembly);
    }

    // Handle rotamer optimization with titration.
    if (manyBodyOptions.getTitration()) {
      logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n");

      // Collect residue numbers.
      List<Integer> resNumberList = new ArrayList<>();
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber());
      }

      // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
      titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
          resNumberList, titrationPH, manyBodyOptions);
      activeAssembly = titrationManyBody.getProtonatedAssembly();
      potentialEnergy = activeAssembly.getPotentialEnergy();
    }
    molecularAssemblies = new MolecularAssembly[]{activeAssembly};
    refinementEnergy = realSpaceOptions.toRealSpaceEnergy(filenames, molecularAssemblies);

    if (lambdaTerm) {
      alchemicalOptions.setFirstSystemAlchemistry(activeAssembly);
      LambdaInterface lambdaInterface = potentialEnergy;
      double lambda = alchemicalOptions.getInitialLambda();
      logger.info(format(" Setting ManyBody softcore lambda to: %5.3f", lambda));
      lambdaInterface.setLambda(lambda);
    }

    RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, refinementEnergy, algorithmListener);
    rotamerOptimization.setPrintFiles(printFiles);
    rotamerOptimization.setWriteEnergyRestart(printFiles);
    rotamerOptimization.setPHRestraint(pHRestraint);
    rotamerOptimization.setpH(titrationPH);

    manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly);

    double[] x = new double[refinementEnergy.getNumberOfVariables()];
    x = refinementEnergy.getCoordinates(x);
    double e = refinementEnergy.energy(x, true);
    logger.info(format("\n Initial target energy: %16.8f ", e));

    selectedResidues = rotamerOptimization.getResidues();
    rotamerOptimization.initFraction(selectedResidues);

    RotamerLibrary.measureRotamers(residueList, false);

    rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()));

    int[] currentRotamers = new int[selectedResidues.size()];

    //Calculate possible permutations for assembly
    try {
      rotamerOptimization.getFractions(selectedResidues.toArray(new Residue[0]), 0, currentRotamers);
      if (pKa) {
        rotamerOptimization.getProtonationPopulations(selectedResidues.toArray(new Residue[0]));
      }
    } catch (Exception ex) {
      logger.severe(" Error calculating rotamer fractions: " + ex.getMessage());
      return this;
    }

    //Collect the Bolztmann weights and calculated offset of each assembly
    boltzmannWeights = rotamerOptimization.getTotalBoltzmann();
    offsets = rotamerOptimization.getRefEnergy();

    //Calculate the populations for the all residue rotamers
    populationArray = rotamerOptimization.getFraction();

    try (FileWriter fileWriter = new FileWriter("populations.txt")) {
      int residueIndex = 0;
      for (Residue residue : selectedResidues) {
        fileWriter.write("\n");
        Rotamer[] rotamers = residue.getRotamers();
        for (int rotIndex = 0; rotIndex < rotamers.length; rotIndex++) {
          String rotPop = format("%.6f", populationArray[residueIndex][rotIndex]);
          fileWriter.write(residue.getName() + residue.getResidueNumber() + "\t" +
              rotamers[rotIndex].toString() + "\t" + rotPop + "\n");
        }
        residueIndex += 1;
      }
      System.out.println("\n Successfully wrote to the populations file.");
    } catch (IOException ioException) {
      logger.severe(" Error writing populations file: " + ioException.getMessage());
    }

    int[][] conformers;
    try {
      conformers = rotamerOptimization.getConformers();
    } catch (Exception ex) {
      logger.severe(" Error getting conformers: " + ex.getMessage());
      return this;
    }

    char[] altLocs = new char[]{'C', 'B', 'A'};
    List<String> residueChainNum = new ArrayList<>();
    for (Residue residue : activeAssembly.getResidueList()) {
      char chain = residue.getChainID();
      int resNum = residue.getResidueNumber();
      residueChainNum.add(String.valueOf(chain) + String.valueOf(resNum));
    }

    File structureFile = new File(filename);
    int[] optimalRotamers = new int[selectedResidues.size()];
    int assemblyIndex = 0;
    String[] rotNames = new String[selectedResidues.size()];
    for (int confIndex = 2; confIndex > -1; confIndex--) {
      List<Residue> conformerResidueList = new ArrayList<>();
      MolecularAssembly conformerAssembly = algorithmFunctions.open(filename);
      if (manyBodyOptions.getTitration()) {
        logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n");

        // Collect residue numbers.
        List<Integer> resNumberList = new ArrayList<>();
        for (Residue residue : residues) {
          resNumberList.add(residue.getResidueNumber());
        }

        // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
        titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
            resNumberList, titrationPH, manyBodyOptions);
        conformerAssembly = titrationManyBody.getProtonatedAssembly();
        potentialEnergy = conformerAssembly.getPotentialEnergy();
      }
      for (int resIndex = 0; resIndex < selectedResidues.size(); resIndex++) {
        Residue residueSelect = selectedResidues.get(resIndex);
        String resChainNum = String.valueOf(residueSelect.getChainID()) + String.valueOf(residueSelect.getResidueNumber());
        int index = residueChainNum.indexOf(resChainNum);
        Residue residue = conformerAssembly.getResidueList().get(index);
        conformerResidueList.add(residue);
        residue.setRotamers(manyBodyOptions.getRotamerLibrary(true));
        Rotamer[] rotamers = residue.getRotamers();
        int rotIndex = conformers[resIndex][confIndex];
        if (populationArray[resIndex][rotIndex] != 0) {
          optimalRotamers[resIndex] = rotIndex;
          RotamerLibrary.applyRotamer(residue, rotamers[rotIndex]);
          double occupancy = populationArray[resIndex][rotIndex];
          boolean diffStates = false;
          for (int i = 2; i > -1; i--) {
            int rotamerInd = conformers[resIndex][i];
            String rotName = rotamers[rotamerInd].getName();
            double occupancyTest = populationArray[resIndex][rotamerInd];
            if (i == 2 && rotNames[resIndex] != null) {
              rotNames[resIndex] = rotName;
            } else if (i < 2 && occupancyTest != 0 && !rotNames[resIndex].contains(rotName)) {
              diffStates = true;
              String newString = rotNames[resIndex] + rotName;
              logger.info(newString);
              rotNames[resIndex] = newString;
            } else if (i == 2) {
              rotNames[resIndex] = rotName;
            }
          }
          for (Atom atom : residue.getAtomList()) {
            if (!residue.getBackboneAtoms().contains(atom) || diffStates) {
              if (occupancy != 1) {
                atom.setAltLoc(altLocs[confIndex]);
                atom.setOccupancy(occupancy);
              } else {
                atom.setOccupancy(occupancy);
                atom.setAltLoc(' ');
              }
            } else {
              atom.setOccupancy(1.0);
              atom.setAltLoc(' ');
            }
          }
        }
      }

      conformerAssemblies[assemblyIndex] = conformerAssembly;
      titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, conformerAssemblies[assemblyIndex], conformerResidueList);
      excludeAtoms.addAll(titrationManyBody.getExcludeAtoms());
      assemblyIndex++;
    }
    PDBFilter pdbFilter = new PDBFilter(structureFile, Arrays.asList(conformerAssemblies), forceField, properties);
    pdbFilter.writeFile(structureFile, false, excludeAtoms, true, true);
    System.setProperty("standardizeAtomNames", "false");

    return this;
  }
}
