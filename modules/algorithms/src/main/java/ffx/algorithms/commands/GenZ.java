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
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard;
import static java.lang.String.format;

/**
 * The ReductionPartition script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Run GenZ function for free energy change.", name = "GenZ")
public class GenZ extends AlgorithmsCommand {

  @Mixin
  private ManyBodyOptions manyBodyOptions;

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @CommandLine.Option(names = {"--resC", "--residueChain"}, paramLabel = "A",
      description = "The chain that is mutating.")
  private String mutatingChain = "A";

  @CommandLine.Option(names = {"-n", "--residueName"}, paramLabel = "ALA",
      description = "Mutant residue.")
  private String resName;

  @CommandLine.Option(names = {"--rEE", "--ro-ensembleEnergy"}, paramLabel = "0.0",
      description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
  private String ensembleEnergy = "0.0";

  @CommandLine.Option(names = {"--un", "--unfolded"}, paramLabel = "false",
      description = "Run the unfolded state tripeptide.")
  private boolean unfolded = false;

  @CommandLine.Option(names = {"--pKa"}, paramLabel = "false",
      description = "Calculating protonation populations for pKa shift.")
  private boolean pKa = false;

  @CommandLine.Option(names = {"--pB", "--printBoltzmann"}, paramLabel = "false",
      description = "Save the Boltzmann weights of protonated residue and total Boltzmann weights.")
  private boolean printBoltzmann = false;

  @CommandLine.Option(names = {"--pF", "--printFiles"}, paramLabel = "false",
      description = "Write to an energy restart file and ensemble file.")
  private boolean printFiles = false;

  @CommandLine.Option(names = {"--rCS", "--recomputeSelf"}, paramLabel = "false",
      description = "Recompute the self energies after loading a restart file.")
  private boolean recomputeSelf = false;

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private List<String> filenames = null;


  ForceFieldEnergy potentialEnergy;
  TitrationManyBody titrationManyBody;
  /**
   * Assembly with mutations after running mutate pdb
   */
  MolecularAssembly mutatedAssembly;
  /**
   * Binding to run mutate pdb
   */
  FFXBinding mutatorBinding;
  /**
   * All the residues
   */
  List<Residue> residues;
  /**
   * Residues included in the partition function
   */
  List<Residue> selectedResidues;
  /**
   * File to save the unfolded state of a protein for free energy prediction
   */
  private String unfoldedFileName;
  /**
   * Population array for residue rotamers
   */
  private double[][] populationArray;

  /**
   * ManyBody Constructor.
   */
  public GenZ() {
    super();
  }

  /**
   * ManyBody Constructor.
   *
   * @param binding The Binding to use.
   */
  public GenZ(FFXBinding binding) {
    super(binding);
  }

  /**
   * GenZ constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public GenZ(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public GenZ run() {
    if (!init()) {
      return this;
    }

    // Get all the important flags from the manybody options
    double titrationPH = manyBodyOptions.getTitrationPH();
    double inclusionCutoff = manyBodyOptions.getInclusionCutoff();
    int mutatingResidue = manyBodyOptions.getInterestedResidue();
    boolean onlyTitration = manyBodyOptions.getOnlyTitration();
    double pHRestraint = manyBodyOptions.getPHRestraint();
    // Set system property to propagate titration
    if (manyBodyOptions.getTitration()) {
      System.setProperty("manybody-titration", "true");
    }

    // If soft coring
    boolean lambdaTerm = alchemicalOptions.hasSoftcore();
    if (lambdaTerm) {
      // Turn on softcore van der Waals
      System.setProperty("lambdaterm", "true");
      // Turn of alchemical electrostatics
      System.setProperty("elec-lambdaterm", "false");
      // Turn on intra-molecular softcore
      System.setProperty("intramolecular-softcore", "true");
    }
    // Set the energy cutoff for permutations to include in the ensemble
    System.setProperty("ro-ensembleEnergy", ensembleEnergy);
    activeAssembly = getActiveAssembly(filenames.get(0));

    //Make an unfolded state assembly when predicting folding free energy difference
    if (unfolded) {
      unfoldedFileName = "wt" + mutatingResidue + ".pdb";
      List<Atom> atoms = activeAssembly.getAtomList();
      Set<Atom> excludeAtoms = new HashSet<>();
      for (Atom atom : atoms) {
        if (atom.getResidueNumber() < mutatingResidue - 1 || atom.getResidueNumber() > mutatingResidue + 1) {
          excludeAtoms.add(atom);
        } else if (atom.getResidueNumber() == mutatingResidue - 1 && "H".equals(atom.getName())) {
          excludeAtoms.add(atom);
        }
      }
      File file = new File(unfoldedFileName);
      PDBFilter pdbFilter = new PDBFilter(file, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties());
      pdbFilter.writeFile(file, false, excludeAtoms, true, true);
      setActiveAssembly(getActiveAssembly(unfoldedFileName));
    }

    //Allocate arrays for different values coming out of the partition function
    double[] boltzmannWeights = new double[2];
    double[] offsets = new double[2];
    double[][] populationArray = new double[1][55];
    double[][] titrateBoltzmann = null;
    double[] protonationBoltzmannSums = null;
    double totalBoltzmann = 0;
    List<Residue> residueList = activeAssembly.getResidueList();

    String mutatedFileName = "";
    //Call the MutatePDB script and mutate the residue of interest
    if (filenames.size() == 1 && mutatingResidue != -1) {
      mutatorBinding = new FFXBinding();
      if (unfolded) {
        // mutatorBinding = new Binding('-r', mutatingResidue.toString(), '-n', resName, unfoldedFileName)
        String[] args = {"-r", String.valueOf(mutatingResidue), "-n", resName, unfoldedFileName};
        mutatorBinding.setVariable("args", args);
      } else {
        // mutatorBinding = new Binding('-r', mutatingResidue.toString(), '-n', resName, filenames.get(0), '--ch', mutatingChain)
        String[] args = {"-r", String.valueOf(mutatingResidue), "-n", resName, "--ch", mutatingChain, unfoldedFileName};
        mutatorBinding.setVariable("args", args);
      }
      MutatePDB mutatePDB = new MutatePDB(mutatorBinding);
      mutatePDB.run();
      mutatedFileName = (String) mutatorBinding.getVariable("versionFileName");
    }

    String listResidues = "";
    // Select residues with alpha carbons within the inclusion cutoff or
    // Select only the titrating residues or the titrating residues and those within the inclusion cutoff
    if ((mutatingResidue != -1 && inclusionCutoff != -1) || onlyTitration) {
      listResidues = manyBodyOptions.selectInclusionResidues(residueList, mutatingResidue, onlyTitration, inclusionCutoff);
    }

    String filename = filenames.get(0);

    //Set the number of assemblies the partition function will be calculated for
    int numLoop = 1;
    if (mutatingResidue != -1) {
      numLoop = 2;
    }

    //Prepare variables for saving out the highest population rotamers (optimal rotamers)
    int[] optimalRotamers;
    Set<Atom> excludeAtoms = new HashSet<>();

    //Calculate all possible permutations for the number of assembles
    for (int j = 0; j < numLoop; j++) {

      // Load the MolecularAssembly second molecular assembly if applicable.
      if (j > 0) {
        if (filenames.size() == 1) {
          mutatedAssembly = getActiveAssembly(mutatedFileName);
          setActiveAssembly(mutatedAssembly);
          logger.info(activeAssembly.getResidueList().toString());
          activeAssembly.getPotentialEnergy().energy();
          filename = mutatedFileName;
        } else {
          setActiveAssembly(getActiveAssembly(filenames.get(j)));
          filename = filenames.get(j);
        }
      }

      if (activeAssembly == null) {
        logger.info(helpString());
        return this;
      }

      CompositeConfiguration properties = activeAssembly.getProperties();

      // Application of rotamers uses side-chain atom naming from the PDB.
      if (properties.getBoolean("standardizeAtomNames", false)) {
        renameAtomsToPDBStandard(activeAssembly);
      }

      // Update the potential energy to match current assembly
      activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);
      potentialEnergy = activeAssembly.getPotentialEnergy();

      // Selecting residues
      if (!pKa || onlyTitration) {
        manyBodyOptions.setListResidues(listResidues);
      }

      // Collect residues to optimize.
      residues = manyBodyOptions.collectResidues(activeAssembly);
      if (residues == null || residues.isEmpty()) {
        logger.info(" There are no residues in the active system to optimize.");
        return this;
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
        MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly();
        setActiveAssembly(protonatedAssembly);
        potentialEnergy = protonatedAssembly.getPotentialEnergy();
      }

      // Turn on softcoring lambda
      if (lambdaTerm) {
        alchemicalOptions.setFirstSystemAlchemistry(activeAssembly);
        LambdaInterface lambdaInterface = (LambdaInterface) potentialEnergy;
        double lambda = alchemicalOptions.getInitialLambda();
        logger.info(format(" Setting ManyBody softcore lambda to: %5.3f", lambda));
        lambdaInterface.setLambda(lambda);
      }

      //Run rotamer optimization with specified parameters
      RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly,
          potentialEnergy, algorithmListener);
      rotamerOptimization.setPrintFiles(printFiles);
      rotamerOptimization.setWriteEnergyRestart(printFiles);
      rotamerOptimization.setPHRestraint(pHRestraint);
      rotamerOptimization.setRecomputeSelf(recomputeSelf);
      rotamerOptimization.setpH(titrationPH);

      manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly);

      // Initialize fractions for selected residues
      selectedResidues = rotamerOptimization.getResidues();
      rotamerOptimization.initFraction(selectedResidues);

      logger.info("\n Initial Potential Energy:");
      potentialEnergy.energy(false, true);

      logger.info("\n Initial Rotamer Torsion Angles:");
      RotamerLibrary.measureRotamers(selectedResidues, false);

      // Run the optimization.
      rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(selectedResidues.size()));

      int[] currentRotamers = new int[selectedResidues.size()];

      //Calculate possible permutations for assembly
      try {
        rotamerOptimization.getFractions(selectedResidues.toArray(new Residue[0]), 0, currentRotamers);
      } catch (Exception e) {
        logger.severe(" Error calculating fractions: " + e.getMessage());
        return this;
      }

      //Collect the Bolztmann weights and calculated offset of each assembly
      boltzmannWeights[j] = rotamerOptimization.getTotalBoltzmann();
      offsets[j] = rotamerOptimization.getRefEnergy();

      //Calculate the populations for the residue rotamers
      populationArray = rotamerOptimization.getFraction();
      if (printBoltzmann) {
        titrateBoltzmann = rotamerOptimization.getPopulationBoltzmann();
        totalBoltzmann = rotamerOptimization.getTotalBoltzmann();
      }

      // Collect the most populous rotamers
      optimalRotamers = rotamerOptimization.getOptimumRotamers();
      if (manyBodyOptions.getTitration()) {
        // Remove excess atoms from titratable residues
        titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, selectedResidues);
      }

      // Calculate the protonation populations
      if (pKa) {
        rotamerOptimization.getProtonationPopulations(selectedResidues.toArray(new Residue[0]));
      }


    }

    //Print information from the fraction protonated calculations
    try (FileWriter fileWriter = new FileWriter("populations.txt")) {
      int residueIndex = 0;
      for (Residue residue : selectedResidues) {
        fileWriter.write("\n");
        protonationBoltzmannSums = new double[selectedResidues.size()];
        // Set sums for to protonated, deprotonated, and tautomer states of titratable residues
        Rotamer[] rotamers = residue.getRotamers();
        for (int rotIndex = 0; rotIndex < rotamers.length; rotIndex++) {
          String rotPop = format("%.6f", populationArray[residueIndex][rotIndex]);
          fileWriter.write(residue.getName() + residue.getResidueNumber() + "\t" +
              rotamers[rotIndex].toString() + "\t" + rotPop + "\n");
          if (pKa) {
            switch (rotamers[rotIndex].getName()) {
              case "HIS":
              case "LYS":
              case "GLH":
              case "ASH":
              case "CYS":
                if (printBoltzmann) {
                  protonationBoltzmannSums[residueIndex] += titrateBoltzmann[residueIndex][rotIndex];
                }
                break;
              default:
                break;
            }
          }

        }
        // Print protonated and total boltzmann values
        if (printBoltzmann) {
          logger.info("\n Residue " + residue.getName() + residue.getResidueNumber() + " Protonated Boltzmann: " +
              protonationBoltzmannSums[residueIndex]);
          logger.info("\n Total Boltzmann: " + totalBoltzmann);
        }
        residueIndex += 1;
      }
      System.out.println("\n Successfully wrote to the populations file.");
    } catch (Exception e) {
      logger.severe("Error writing populations file: " + e.getMessage());
    }


    // Save the pdb file with the most popular rotamers for all residues included in the partition function
    System.setProperty("standardizeAtomNames", "false");
    File modelFile = saveDirFile(activeAssembly.getFile());
    PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
        activeAssembly.getProperties());
    if (manyBodyOptions.getTitration()) {
      String remark = format("Titration pH: %6.3f", titrationPH);
      if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true, new String[]{remark})) {
        logger.info(format(" Save failed for %s", activeAssembly));
      }
    } else {
      if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
        logger.info(format(" Save failed for %s", activeAssembly));
      }
    }
    if (mutatingResidue != -1) {
      //Calculate Gibbs free energy change of mutating residues
      double gibbs = -(0.6) * (Math.log(boltzmannWeights[1] / boltzmannWeights[0]));
      logger.info("\n Gibbs Free Energy Change: " + gibbs);
    }


    return this;
  }


  /**
   * Returns the potential energy of the active assembly. Used during testing assertions.
   *
   * @return potentialEnergy Potential energy of the active assembly.
   */
  public ForceFieldEnergy getPotential() {
    return potentialEnergy;
  }

  public double[][] getPopulationArray() {
    return populationArray;
  }
}
