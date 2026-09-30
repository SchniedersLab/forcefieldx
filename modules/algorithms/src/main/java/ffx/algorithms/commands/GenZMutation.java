package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.ManyBodyOptions;
import ffx.algorithms.commands.MutatePDB;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.algorithms.optimize.TitrationManyBody;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Constants;
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
import static java.lang.Double.parseDouble;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.log;

@Command(description = "Run mutation calculations for residues.", name = "GenZMutation")
public class GenZMutation extends AlgorithmsCommand {

  @Mixin
  private ManyBodyOptions manyBodyOptions;

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Option(names = {"--resC", "--residueChain"}, paramLabel = "A", defaultValue = "A",
          description = "The chain that is mutating.")
  private String mutatingChain;

  @Option(names = {"-n", "--residueName"}, paramLabel = "ALA", defaultValue = "",
          description = "Mutant residue.")
  private String resName;

  @Option(names = {"--rEE", "--ro-ensembleEnergy"}, paramLabel = "0.0", defaultValue = "0.0",
          description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
  private String ensembleEnergy;

  @Option(names = {"--un", "--unfolded"}, paramLabel = "false", defaultValue = "false",
          description = "Run the unfolded state tripeptide.")
  private boolean unfolded;

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
  public GenZMutation() {
    super();
  }

  /**
   * ManyBody Constructor.
   *
   * @param binding The Binding to use.
   */
  public GenZMutation(FFXBinding binding) {
    super(binding);
  }

  /**
   * GenZ constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public GenZMutation(String[] args) {
    super(args);
  }


  /**
     * {@inheritDoc}
     */
    @Override
    public GenZMutation run() {
        if (!init()) {
            return this;
        }

      // Get all the important flags from the manybody options
      double inclusionCutoff = manyBodyOptions.getInclusionCutoff();
      int mutatingResidue = manyBodyOptions.getInterestedResidue();
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

      // Load the MolecularAssembly.
      activeAssembly = getActiveAssembly(filename);
      if (activeAssembly == null) {
        logger.info(helpString());
        return this;
      }

      // Set the filename.
      filename = activeAssembly.getFile().getAbsolutePath();

        if (unfolded) {
            logger.info("Running unfolded state calculations...");
            // Add unfolded state-specific logic here.
        } else {
            logger.info(format("Mutating chain %s to residue %s.", mutatingChain, resName));
            // Add mutation-specific logic here.
        }

      // Make an unfolded state assembly when predicting folding free energy difference
      String unfoldedFileName = "";
      if (unfolded) {
        // File to save the unfolded state of a protein for free energy prediction
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

        // Load the unfolded state system.
        activeAssembly = getActiveAssembly(unfoldedFileName);
        filename = activeAssembly.getFile().getAbsolutePath();
      }

      // Allocate arrays for different values coming out of the partition function
      double[] boltzmannWeights = new double[2];
      double[][] titrateBoltzmann = null;
      double totalBoltzmann = 0;
      List<Residue> residueList = activeAssembly.getResidueList();

      String mutatedFileName = "";
      // Call the MutatePDB script and mutate the residue of interest
      if (mutatingResidue != -1) {
        FFXBinding mutatorBinding = new FFXBinding();
        if (unfolded) {
          String[] args = {"-r", String.valueOf(mutatingResidue), "-n", resName, unfoldedFileName};
          mutatorBinding.setVariable("args", args);
        } else {
          String[] args = {"-r", String.valueOf(mutatingResidue), "-n", resName, "--ch", mutatingChain, filename};
          mutatorBinding.setVariable("args", args);
        }
        MutatePDB mutatePDB = new MutatePDB(mutatorBinding);
        mutatePDB.run();
        mutatedFileName = (String) mutatorBinding.getVariable("versionFileName");
      }


      // Set the number of assemblies the partition function will be calculated for
      int numLoop = 1;
      if (mutatingResidue != -1) {
        numLoop = 2;
      }

      //Prepare variables for saving out the highest population rotamers (optimal rotamers)
      int[] optimalRotamers;
      Set<Atom> excludeAtoms = new HashSet<>();

      // Calculate all possible permutations for the number of assembles
      for (int j = 0; j < numLoop; j++) {

        // Load the MolecularAssembly second molecular assembly if applicable.
        if (j > 0) {
          activeAssembly = getActiveAssembly(mutatedFileName);
          activeAssembly.getPotentialEnergy().energy();
          filename = activeAssembly.getFile().getAbsolutePath();
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

        String listResidues = "";
        // Select residues with alpha carbons within the inclusion cutoff
        listResidues = manyBodyOptions.selectInclusionResidues(residueList, mutatingResidue, false, inclusionCutoff);

        manyBodyOptions.setListResidues(listResidues);

        // Collect residues to optimize.
        List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly);
        if (residues == null || residues.isEmpty()) {
          logger.info(" There are no residues in the active system to optimize.");
          return this;
        }

        // Handle rotamer optimization with titration.
        TitrationManyBody titrationManyBody = null;
        if (manyBodyOptions.getTitration()) {
          logger.info("\n Adding titration hydrogen to : " + filename + "\n");

          // Collect residue numbers.
          List<Integer> resNumberList = new ArrayList<>();
          for (Residue residue : residues) {
            resNumberList.add(residue.getResidueNumber());
          }

          // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
          titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
                  resNumberList, manyBodyOptions.getTitrationPH(), manyBodyOptions);
          MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly();
          activeAssembly = protonatedAssembly;
          potentialEnergy = protonatedAssembly.getPotentialEnergy();
        }

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
        rotamerOptimization.setPHRestraint(manyBodyOptions.getPHRestraint());
        rotamerOptimization.setRecomputeSelf(recomputeSelf);
        rotamerOptimization.setpH(manyBodyOptions.getTitrationPH());

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
        boltzmannWeights[j] = rotamerOptimization.getTotalBoltzmann();

        // Calculate the populations for the residue rotamers
        populationArray = rotamerOptimization.getFraction();

        // Collect the most populous rotamers
        optimalRotamers = rotamerOptimization.getOptimumRotamers();
        if (manyBodyOptions.getTitration()) {
          // Remove excess atoms from titratable residues
          titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, selectedResidues);
        }


        // Print information from the fraction protonated calculations
        String populationFilename = "populations.txt";
        if (j > 0) {
          populationFilename = "populations."  + j + ".txt";
        }
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

            }

            residueIndex += 1;
          }
          logger.info("\n Total Boltzmann: " + totalBoltzmann);
          logger.info("\n Successfully wrote to the populations file: " + populationFilename);
        } catch (Exception e) {
          logger.severe("Error writing populations file: " + e.getMessage());
        }
      }

      // Save the pdb file with the most popular rotamers for all residues included in the partition function
      System.setProperty("standardizeAtomNames", "false");
      File modelFile = saveDirFile(activeAssembly.getFile());
      PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
              activeAssembly.getProperties());
      if (manyBodyOptions.getTitration()) {
        String remark = format("Titration pH: %6.3f", manyBodyOptions.getTitrationPH());
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true, new String[]{remark})) {
          logger.info(format(" Save failed for %s", activeAssembly));
        }
      } else {
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
          logger.info(format(" Save failed for %s", activeAssembly));
        }
      }

      if (mutatingResidue != -1) {
        CompositeConfiguration properties = activeAssembly.getProperties();
        String temp = properties.getString("temperature", "298.15");
        double temperature = 298.15;
        if (temp != null) {
          temperature = parseDouble(temp);
        }
        double kT = temperature * Constants.R;
        // Calculate free energy difference for mutating a residue.
        double dG = -kT * log(boltzmannWeights[1] / boltzmannWeights[0]);
        logger.info(format("\n Mutation Free Energy Difference: %12.8f (kcal/mol)", dG));
      }


      return this;
    }
}