package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.ManyBodyOptions;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.algorithms.optimize.TitrationManyBody;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
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

@Command(description = "Run pKa calculations for protonation populations.", name = "GenZpKa")
public class GenZpKa extends AlgorithmsCommand {

  @Mixin
  private ManyBodyOptions manyBodyOptions;

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Option(names = {"--rEE", "--ro-ensembleEnergy"}, paramLabel = "0.0", defaultValue = "0.0",
          description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
  private String ensembleEnergy;

  @Option(names = {"--pB", "--printBoltzmann"}, paramLabel = "false", defaultValue = "false",
          description = "Save the Boltzmann weights of protonated residue and total Boltzmann weights.")
  private boolean printBoltzmann;

  @Option(names = {"--pF", "--printFiles"}, paramLabel = "false", defaultValue = "false",
          description = "Write to an energy restart file and ensemble file.")
  private boolean printFiles;

  @Option(names = {"--rCS", "--recomputeSelf"}, paramLabel = "false", defaultValue = "false",
          description = "Recompute the self energies after loading a restart file.")
  private boolean recomputeSelf;

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file", defaultValue = "",
          description = "XYZ or PDB input file.")
  private String filename;

  ForceFieldEnergy potentialEnergy;
  /**
   * Populations for each rotamer of each residue.
   */
  private double[][] populationArray;

  /**
   * ManyBody Constructor.
   */
  public GenZpKa() {
    super();
  }

  /**
   * ManyBody Constructor.
   *
   * @param binding The Binding to use.
   */
  public GenZpKa(FFXBinding binding) {
    super(binding);
  }

  /**
   * GenZ constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public GenZpKa(String[] args) {
    super(args);
  }


  /**
   * {@inheritDoc}
   */
  @Override
  public GenZpKa run() {
    if (!init()) {
      return this;
    }

    // Get all the important flags from the manybody options
    double titrationPH = manyBodyOptions.getTitrationPH();
    double inclusionCutoff = manyBodyOptions.getInclusionCutoff();
    boolean onlyTitration = manyBodyOptions.getOnlyTitration();
    double pHRestraint = manyBodyOptions.getPHRestraint();
    // Set system property to propagate titration
    System.setProperty("manybody-titration", "true");

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
    System.setProperty("standardizeAtomNames", "false");

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    // Allocate arrays for different values coming out of the partition function
    double boltzmannWeights = 0;
    double[][] titrateBoltzmann = null;
    double totalBoltzmann = 0;
    List<Residue> residueList = activeAssembly.getResidueList();

    String listResidues = "";
    // Select residues with alpha carbons within the inclusion cutoff or
    // Select only the titrating residues or the titrating residues and those within the inclusion cutoff
    if (onlyTitration) {
      listResidues = manyBodyOptions.selectInclusionResidues(residueList, -1, onlyTitration, inclusionCutoff);
    }

    //Prepare variables for saving out the highest population rotamers (optimal rotamers)
    int[] optimalRotamers;
    Set<Atom> excludeAtoms = new HashSet<>();

    // Calculate all possible permutations for the number of assembles

    // Load the MolecularAssembly second molecular assembly if applicable.

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
    if (onlyTitration) {
      manyBodyOptions.setListResidues(listResidues);
    }

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly);
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.");
      return this;
    }

    // Handle rotamer optimization with titration.
    TitrationManyBody titrationManyBody = null;
    logger.info("\n Adding titration hydrogen to : " + filename + "\n");

    // Collect residue numbers.
    List<Integer> resNumberList = new ArrayList<>();
    for (Residue residue : residues) {
      resNumberList.add(residue.getResidueNumber());
    }

    // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
    titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
            resNumberList, titrationPH, manyBodyOptions);
    MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly();
    activeAssembly = protonatedAssembly;
    potentialEnergy = protonatedAssembly.getPotentialEnergy();


    // Turn on softcoring lambda
    if (lambdaTerm) {
      alchemicalOptions.setFirstSystemAlchemistry(activeAssembly);
      LambdaInterface lambdaInterface = potentialEnergy;
      double lambda = alchemicalOptions.getInitialLambda();
      logger.info(format(" Setting ManyBody softcore lambda to: %5.3f", lambda));
      lambdaInterface.setLambda(lambda);
    }

    //Run rotamer optimization with specified parameters
    RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, potentialEnergy, algorithmListener);
    rotamerOptimization.setPrintFiles(printFiles);
    rotamerOptimization.setWriteEnergyRestart(printFiles);
    rotamerOptimization.setPHRestraint(pHRestraint);
    rotamerOptimization.setRecomputeSelf(recomputeSelf);
    rotamerOptimization.setpH(titrationPH);

    manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly);

    // Initialize fractions for selected residues
    List<Residue> selectedResidues = rotamerOptimization.getResidues();
    rotamerOptimization.initFraction(selectedResidues);

    logger.info("\n Initial Potential Energy:");
    potentialEnergy.energy(false, true);

    logger.info("\n Initial Rotamer Torsion Angles:");
    RotamerLibrary.measureRotamers(selectedResidues, false);

    // Run the optimization.
    rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(selectedResidues.size()));

    int[] currentRotamers = new int[selectedResidues.size()];

    // Calculate possible permutations for assembly
    try {
      rotamerOptimization.getFractions(selectedResidues.toArray(new Residue[0]), 0, currentRotamers);
    } catch (Exception e) {
      logger.severe(" Error calculating fractions: " + e.getMessage());
      return this;
    }

    // Collect the Boltzmann weights and calculated offset of each assembly
    boltzmannWeights = rotamerOptimization.getTotalBoltzmann();

    // Calculate the populations for the residue rotamers
    populationArray = rotamerOptimization.getFraction();
    if (printBoltzmann) {
      titrateBoltzmann = rotamerOptimization.getPopulationBoltzmann();
      totalBoltzmann = rotamerOptimization.getTotalBoltzmann();
    }

    // Collect the most populous rotamers
    optimalRotamers = rotamerOptimization.getOptimumRotamers();
    // Remove excess atoms from titratable residues
    titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, selectedResidues);


    // Calculate the protonation populations
    rotamerOptimization.getProtonationPopulations(selectedResidues.toArray(new Residue[0]));


    // Print information from the fraction protonated calculations
    String populationFilename = "populations.txt";

    populationFilename = "populations" + ".txt";
    try (FileWriter fileWriter = new FileWriter(populationFilename)) {
      int residueIndex = 0;
      for (Residue residue : selectedResidues) {
        fileWriter.write("\n");
        double protonationBoltzmannSum = 0.0;
        // Set sums for to protonated, deprotonated, and tautomer states of titratable residues
        Rotamer[] rotamers = residue.getRotamers();
        for (int rotIndex = 0; rotIndex < rotamers.length; rotIndex++) {
          String rotPop = format("%.6f", populationArray[residueIndex][rotIndex]);
          fileWriter.write(residue.getName() + residue.getResidueNumber() + "\t" +
                  rotamers[rotIndex].toString() + "\t" + rotPop + "\n");

          switch (rotamers[rotIndex].getName()) {
            case "HIS":
            case "LYS":
            case "GLH":
            case "ASH":
            case "CYS":
              if (printBoltzmann) {
                protonationBoltzmannSum += titrateBoltzmann[residueIndex][rotIndex];
              }
              break;
            default:
              break;

          }

        }
        // Print protonated and total boltzmann values
        if (printBoltzmann) {
          logger.info("\n Residue " + residue.getName() + residue.getResidueNumber()
                  + " Protonated Boltzmann: " + protonationBoltzmannSum);
        }
        residueIndex += 1;
      }
      logger.info("\n Total Boltzmann: " + totalBoltzmann);
      logger.info("\n Successfully wrote to the populations file: " + populationFilename);
    } catch (Exception e) {
      logger.severe("Error writing populations file: " + e.getMessage());
    }

    // Save the pdb file with the most popular rotamers for all residues included in the partition function

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

    return this;
  }

  /**
   * The population for each rotamer of each residue.
   *
   * @return The population array.
   */
  public double[][] getPopulationArray() {
    return populationArray;
  }

  /**
   * Returns the potential energy of the active assembly. Used during testing assertions.
   *
   * @return potentialEnergy Potential energy of the active assembly.
   */
  public ForceFieldEnergy getPotential() {
    return potentialEnergy;
  }
}