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
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.numerics.math.RunningStatistics;
import ffx.numerics.math.SummaryStatistics;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.XPHFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.lang3.ArrayUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static java.lang.String.format;

/**
 * Use the Rao-Blackwell Estimator to estimate a free energy difference of protonation for a CpHMD system.
 * <br>
 * Usage: Umod calculation for model compounds
 * <br>
 * ffxc test.RaoBlackwellEstimator [options] &lt;filename&gt; [file2...];
 */
@Command(description = " Use the Rao-Blackwell estimator to get a free energy difference for residues in a CpHMD system.", name = "test.RaoBlackwellEstimator")
public class RaoBlackwellEstimator extends AlgorithmsCommand {

  @Option(names = {"--aFi", "--arcFile"}, paramLabel = "traj",
          description = "A file containing the the PDB from which to build the ExtendedSystem. There is currently no default.")
  private String arcFileName = null;

  @Option(names = {"--numSnaps"}, paramLabel = "-1", defaultValue = "-1",
          description = "Number of snapshots to use from an archive file. Default is all.")
  private int numSnaps;

  @Option(names = {"--specifiedResidues", "--sR"}, paramLabel = "<selection>", defaultValue = "",
          description = "Specified residues to do analysis.")
  private String specified;

  @Option(names = {"--startSnap"}, paramLabel = "-1", defaultValue = "-1",
          description = "Start energy evaluations at a snap other than 2.")
  private int startSnap;

  @Option(names = {"--bootstrapIter"}, paramLabel = "100000", defaultValue = "100000",
          description = "Number of bootstrap iterations. Set -1 for no bootstrapping.")
  private int bootstrapIter;

  @Option(names = {"--skip"}, paramLabel = "-1", defaultValue = "-1",
          description = "Calculate energies on snaps with this interval.")
  private int skip;

  @Option(names = {"--writeFrequency"}, paramLabel = "100", defaultValue = "100",
          description = "Calculate the RBE and print at this snapshot read frequency.")
  private int writeFrequency;

  @Parameters(arity = "1..*", paramLabel = "files",
          description = "PDB input file in the same directory as the ARC file.")
  private String filename;

  private Potential forceFieldEnergy;
  private ArrayList<Double>[][] oneZeroDeltaLists;
  private ArrayList<Double>[][] tautomerOneZeroDeltaList;
  private int numESVs;
  private int numTautomerESVs;

  /**
   * RaoBlackwellEstimator Constructor.
   */
  public RaoBlackwellEstimator() {
    super();
  }

  /**
   * RaoBlackwellEstimator Constructor.
   * @param binding The Binding to use.
   */
  public RaoBlackwellEstimator(FFXBinding binding) {
    super(binding);
  }

  /**
   * RaoBlackwellEstimator constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public RaoBlackwellEstimator(String[] args) {
    super(args);
  }

  @Override
  @SuppressWarnings("unchecked")
  public RaoBlackwellEstimator run() {

    if (!init()) {
      return this;
    }

    // See if the ARC file exists
    File arcFile = new File(arcFileName);
    if(!arcFile.exists()){
      logger.severe(format(" ARC file %s does not exist.", arcFile));
    }
    else{
      logger.info(format("Using ARC file %s.", arcFile));
    }

    boolean bootstrap = false;
    if(bootstrapIter >= 50)
    {
      bootstrap = true;
    } else if (bootstrapIter != -1){
      logger.severe("Too few bootstrap iterations specified. Must be at least 50.");
    }

    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }
    forceFieldEnergy = activeAssembly.getPotentialEnergy();

    // Set the filename.
    String filename = activeAssembly.getFile().getAbsolutePath();

    // Initialize and attach extended system first.
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, 7.0, null);

    // Set up a special residue if one is specified.
    Residue specialResidue = null;
    int numberOfStates = 1; // Regular Rao-Blackwell Estimator if this does not change (special residue does not change until j != 0 in main loop)
    int[][] states = null;
    if(esvSystem.getSpecialResidueList().size() > 1){
      logger.severe(" Multiple special residues were identified in the key file. " +
              "Only one can be specified with this algorithm.");
    } else if (esvSystem.getSpecialResidueList().size() == 1) {
      int specialResidueNumber = esvSystem.getSpecialResidueList().get(0).intValue();
      for (Residue residue : esvSystem.getTitratingResidueList()) {
        if (residue.getResidueNumber() == specialResidueNumber) {
          specialResidue = residue;
        }
      }
      if(specialResidue != null){
        numberOfStates = !esvSystem.isTautomer(specialResidue) ? 3 : 4;
        switch (specialResidue.getName()) {
        // How this array is used later in the code --> states[:][0] = titration   states[:][1] = tautomer
          case "ASD":
          case "GLD":
            states = new int[3][2];
            states[0][0] = 0;
            states[0][1] = 0;
            states[1][0] = 1;
            states[1][1] = 0;
            states[2][0] = 1;
            states[2][1] = 1;
            break;

          case "HIS":
            states = new int[3][2];
            states[0][0] = 0;
            states[0][1] = 0;
            states[1][0] = 0;
            states[1][1] = 1;
            states[2][0] = 1;
            states[2][1] = 0;
            break;

          case "LYS":
          case "CYS":
            states = new int[2][2];
            states[0][0] = 0;
            states[0][1] = 0; // Ignored
            states[1][0] = 1;
            states[1][1] = 0; // Ignored
            break;
        }
        // Specifies the different states that the special residue will be evaluated in
      } else {
        logger.severe(" The special residue specified in the key file was not found in the titrating residue list.");
      }
    }

    // Look for specified residues
    // Convert string to int array
    int[] specifiedResidues = null;
    if(specified != null && !specified.isEmpty()){
      String[] specifiedResiduesString = specified.split(",");
      specifiedResidues = new int[specifiedResiduesString.length];
      for (int i = 0; i < specifiedResiduesString.length; i++) {
        specifiedResidues[i] = Integer.parseInt(specifiedResiduesString[i].trim());
      }
    }
    ArrayList<Residue> onlyResidues = new ArrayList<>();
    ArrayList<Integer> onlyResidueIndices = new ArrayList<>();
    if(specifiedResidues != null){
      for (int i = 0; i < esvSystem.getTitratingResidueList().size(); i++) {
        Residue residue = esvSystem.getTitratingResidueList().get(i);
        if (ArrayUtils.contains(specifiedResidues, residue.getResidueNumber())) {
          onlyResidues.add(residue);
          onlyResidueIndices.add(i);
        }
      }
      if(onlyResidues.size() != specifiedResidues.length){
        logger.severe("Could not find all residues from --specifiedResidues input.");
      }
    }
    else{
      for (int i = 0; i < esvSystem.getTitratingResidueList().size(); i++) {
        onlyResidueIndices.add(i);
      }
    }

    // Create the oneZeroDeltaLists and tautomerOneZeroDeltaList arrays
    // Make a list for ESV's energy differences (energy evals are done at tautomer = 0 for these arrays)
    numESVs = esvSystem.getTitratingResidueList().size();
    oneZeroDeltaLists = new ArrayList[numESVs][numberOfStates + 1];
    for (int i = 0; i < numESVs; i++) {
      for (int j = 0; j < numberOfStates + 1; j++) {
        oneZeroDeltaLists[i][j] = new ArrayList<>();
      }
    }
    // Make a list for tautomerizing ESV's energy differences (energy evals are done at tautomer = 1 for these arrays)
    numTautomerESVs = esvSystem.getTautomerizingResidueList().size();
    tautomerOneZeroDeltaList = new ArrayList[numTautomerESVs][numberOfStates + 1];
    for (int i = 0; i < numTautomerESVs; i++) {
      for (int j = 0; j < numberOfStates + 1; j++) {
        tautomerOneZeroDeltaList[i][j] = new ArrayList<>();
      }
    }

    // Set up the XPHFilter
    activeAssembly.setFile(arcFile);
    XPHFilter xphFilter = new XPHFilter(
            arcFile,
            activeAssembly,
            activeAssembly.getForceField(),
            activeAssembly.getProperties(),
            esvSystem);
    xphFilter.readFile();

    // Set up the force field energy and extended system
    esvSystem.setFixedTitrationState(true);
    esvSystem.setFixedTautomerState(true);
    ((ForceFieldEnergy) forceFieldEnergy).attachExtendedSystem(esvSystem);
    logger.info(format(" Attached extended system with %d residues.", numESVs));

    // Read in first energy snapshot
    double[] x = new double[forceFieldEnergy.getNumberOfVariables()];
    forceFieldEnergy.getCoordinates(x);
    forceFieldEnergy.energy(x, true);

    // Get pH from ARC file
    double pH = 0.0;
    String[] parts = xphFilter.getRemarkLines()[0].split(" ");
    for(int i = 0; i < parts.length; i++) {
      if (parts[i].contains("pH")) {
        pH = Double.parseDouble(parts[i+1]);
      }
    }
    logger.info("\n Setting constant pH to " + pH + ".");
    esvSystem.setConstantPh(pH);

    // Get and log the number of snapshots to use.
    int evals = 0;
    if(numSnaps != -1) {
      logger.info(format(" Using %d snapshots.", numSnaps));
    }
    else {
      logger.info(format(" Using all %d snapshots.", xphFilter.countNumModels()));
    }

    // Read through the ARC file.
    while(xphFilter.readNext()) {
      // Read through snaps we aren't interested in
      if(startSnap != -1 && startSnap > 2 && evals == 0){
        for(int i = 0; i < startSnap - 2; i++){
          xphFilter.readNext();
        }
      }

      // Get coordinates/energies with each new snap and calculate energy differences
      forceFieldEnergy.getCoordinates(x);
      for (int i : onlyResidueIndices) {
        double titrationState = 0;
        double tautomerState = 0;
        if(specialResidue != null){
          titrationState = esvSystem.getTitrationLambda(specialResidue);
          tautomerState = esvSystem.getTautomerLambda(specialResidue);
        }

        // Loop through the special residues different possible states (tautomers/titrations)
        for (int j = 0; j < numberOfStates; j++) {
          if (j != 0) {
            esvSystem.setTitrationLambda(specialResidue, states[j-1][0], false);
            if(numberOfStates != 3) {
              esvSystem.setTautomerLambda(specialResidue, states[j-1][1], false);
            }
          }

          // Evaluate energy (regular Rao-Blackwell if no special residue is found)
          ArrayList<Double> results = getZeroOneDeltas(i, esvSystem, (ForceFieldEnergy) forceFieldEnergy, x);

          Residue res = esvSystem.getTitratingResidueList().get(i);
          if(esvSystem.isTautomer(res)){
            tautomerOneZeroDeltaList[esvSystem.getTautomerizingResidueList().indexOf(res)][j].add(results.get(0));
            oneZeroDeltaLists[i][j].add(results.get(1));
          } else{
            oneZeroDeltaLists[i][j].add(results.get(0));
          }
        }

        if(specialResidue == esvSystem.getExtendedResidueList().get(i)) {
          break;
        }

        // Reset special residue for next residue
        if(specialResidue != null) {
          esvSystem.setTitrationLambda(specialResidue, titrationState, false);
          esvSystem.setTautomerLambda(specialResidue, tautomerState, false);
        }
      }

      evals++;
      if(evals % writeFrequency == 0 || evals == numSnaps) {
        // Calculate the Rao-Blackwell estimator for each residue.
        int tautomerCount = 0;
        double[][] energyLists = new double[numESVs][numberOfStates];
        double[][] energyStdLists = new double[numESVs][numberOfStates];
        double[][] tautomerEnergyLists = new double[numTautomerESVs][numberOfStates];
        double[][] tautomerEnergyStdLists = new double[numTautomerESVs][numberOfStates];
        for (int i : onlyResidueIndices) {
          Residue res = esvSystem.getExtendedResidueList().get(i);
          logger.info("\n Performing Rao-Blackwell Estimator on " + res.getAminoAcid3() + ".");
          if (bootstrap) {
            logger.info("  Performing bootstrap with " + bootstrapIter + " iterations.");
          } else {
            logger.info("  Performing RBE without bootstrap. Ignore standard deviation values.");
          }
          for (int j = 0; j < numberOfStates; j++) {
            double[] bootstrapMeanStd = RBE(oneZeroDeltaLists[i][j], bootstrap, bootstrapIter);
            energyLists[i][j] = bootstrapMeanStd[0];
            if (bootstrap) {
              energyStdLists[i][j] = bootstrapMeanStd[1];
            }

            if (esvSystem.getTautomerizingResidueList().contains(res)) {
              bootstrapMeanStd = RBE(tautomerOneZeroDeltaList[esvSystem.getTautomerizingResidueList().indexOf(res)][j],
                      bootstrap, bootstrapIter);
              tautomerEnergyLists[tautomerCount][j] = bootstrapMeanStd[0];
              if (bootstrap) {
                tautomerEnergyStdLists[tautomerCount][j] = bootstrapMeanStd[1];
              }
            }
          }
          if (esvSystem.isTautomer(res)) {
            tautomerCount++;
          }
          if (specialResidue == res) {
            break;
          }
        }

        // Print the results.
        printResults(specialResidue, esvSystem, energyLists, energyStdLists, tautomerEnergyLists, tautomerEnergyStdLists,
                states, numberOfStates, numESVs, onlyResidueIndices);
      }
      if (numSnaps != -1 && evals >= numSnaps) {
        break;
      }

      if(skip != -1){
        for(int i = 0; i < skip-1; i++){
          xphFilter.readNext();
        }
      }
    }
    return this;
  }


  // Calculate the RBE
  private static double[] RBE(ArrayList<Double> deltaUList, boolean bootstrap, int bootstrapIter){
    ArrayList<Double> deltaU = deltaUList;
    double temperature = 298.0;
    double boltzmann = 0.001985875;
    double beta = 1.0 / (temperature * boltzmann);

    ArrayList<Double> deltaExp = exp(mult(-beta, deltaU));
    ArrayList<Double> numerator = div(mult(beta, mult(deltaU, deltaExp)), subtract(1.0, deltaExp));
    ArrayList<Double> denominator = div(mult(beta, deltaU), subtract(1.0, deltaExp));
    // Calculate averages, std devs., and uncertainties of the numerator and denominator distributions using bootstrap.
    double[] deltaGRBE = bootstrap ?
            bootStrap(numerator, denominator, bootstrapIter) :
            new double[] {-(1.0 / beta) * Math.log(average(numerator) / average(denominator))};
    return deltaGRBE;
  }

  // Iterate calculating the RBE
  private static double[] bootStrap(ArrayList<Double> numerator, ArrayList<Double> denominator, int iter) {
    RunningStatistics estimates = new RunningStatistics();
    for (int k = 0; k < iter; k++) {
      Random rng = new Random();
      int[] trial = getBootstrapIndices(numerator.size(), rng);
      double estimate = estimateDg(numerator, denominator, trial);
      estimates.addValue(estimate);
    }
    SummaryStatistics stats = new SummaryStatistics(estimates);
    return new double[] {stats.mean, stats.getSd()};
  }

  // Calculate the RBE based on given indicies
  private static double estimateDg(ArrayList<Double> num, ArrayList<Double> denom, int[] index) {
    double temperature = 298.0;
    double boltzmann = 0.001985875;
    double beta = 1.0 / (temperature * boltzmann);
    ArrayList<Double> numerator = new ArrayList<>();
    numerator.ensureCapacity(index.length);
    ArrayList<Double> denominator = new ArrayList<>();
    denominator.ensureCapacity(index.length);

    for(int i = 0; i < index.length; i++) {
      numerator.add(num.get(index[i]));
      denominator.add(denom.get(index[i]));
    }

    return -(1.0 / beta) * Math.log(average(numerator) / average(denominator));
  }

  // Helper/Printing Methods
  private static void printResults(Residue specialResidue, ExtendedSystem esvSystem, double[][] energyLists, double[][] energyStdLists,
                              double[][] tautomerEnergyLists, double[][]tautomerStdLists, int[][] states,
                              int numberOfStates, int numESVs, ArrayList<Integer> onlyResidueIndex)
  {
    logger.info("\n Rao-Blackwell Estimator Results: ");
    ArrayList<String> line = new ArrayList<>();
    if(specialResidue != null){
      logger.info(" Special Residue: " + specialResidue.toString());
      if(esvSystem.isTautomer(specialResidue)){
        logger.info(format("  %-10s %-10s %-23s %-28s %-28s %-28s", "Residue", "Tautomer", "DeltaGTitr", "DeltaG-SpecialRes=(" + states[0][0] + "," + states[0][1] + ")", "DeltaG-SpecialRes=(" + states[1][0] + "," + states[1][1] + ")", "DeltaG-SpecialRes=(" + states[2][0] + "," + states[2][1] + ")"));
      } else{
        logger.info(format("  %-10s %-10s %-23s %-28s %-28s", "Residue", "Tautomer", "DeltaGTitr","DeltaG-SpecialRes=(" + states[0][0] + "," + states[0][1] + ")", "DeltaG-SpecialRes=(" + states[1][0] + "," + states[1][1] + ")"));
      }
    } else
    {
      logger.info(format("  %-10s %-10s %-23s", "Residue", "Tautomer", "DeltaGTitr"));
    }
    int tautomerCount = 0;
    for(int i : onlyResidueIndex){
      Residue res = esvSystem.getExtendedResidueList().get(i);
      line.add(res.toString());
      line.add("0");
      line.add(Double.toString(energyLists[i][0]));
      line.add(Double.toString(energyStdLists[i][0]));
      for(int j = 1; j < numberOfStates; j++){
        line.add(Double.toString(energyLists[i][j]));
        line.add(Double.toString(energyStdLists[i][j]));
      }
      if(specialResidue != null && esvSystem.isTautomer(specialResidue)) {
        logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f    %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f",
                line.get(0),
                line.get(1),
                Double.parseDouble(line.get(2)),
                Double.parseDouble(line.get(3)),
                Double.parseDouble(line.get(4)),
                Double.parseDouble(line.get(5)),
                Double.parseDouble(line.get(6)),
                Double.parseDouble(line.get(7)),
                Double.parseDouble(line.get(8)),
                Double.parseDouble(line.get(9))));
      } else if (specialResidue != null) {
        logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f    %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f",
                line.get(0),
                line.get(1),
                Double.parseDouble(line.get(2)),
                Double.parseDouble(line.get(3)),
                Double.parseDouble(line.get(4)),
                Double.parseDouble(line.get(5)),
                Double.parseDouble(line.get(6)),
                Double.parseDouble(line.get(7))));
      } else {
        logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f",
                line.get(0),
                line.get(1),
                Double.parseDouble(line.get(2)),
                Double.parseDouble(line.get(3))));
      }
      line.clear();

      if(esvSystem.isTautomer(res)) {
        line.add(res.toString());
        line.add("1");
        line.add(Double.toString(tautomerEnergyLists[tautomerCount][0]));
        line.add(Double.toString(tautomerStdLists[tautomerCount][0]));
        for(int j = 1; j < numberOfStates; j++){
          line.add(Double.toString(tautomerEnergyLists[tautomerCount][j]));
            line.add(Double.toString(tautomerStdLists[tautomerCount][j]));
        }
        tautomerCount++;
        if(specialResidue != null && esvSystem.isTautomer(specialResidue)) {
          logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f    %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f",
                  line.get(0),
                  line.get(1),
                  Double.parseDouble(line.get(2)),
                  Double.parseDouble(line.get(3)),
                  Double.parseDouble(line.get(4)),
                  Double.parseDouble(line.get(5)),
                  Double.parseDouble(line.get(6)),
                  Double.parseDouble(line.get(7)),
                  Double.parseDouble(line.get(8)),
                  Double.parseDouble(line.get(9))));
        } else if (specialResidue != null) {
          logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f    %-10.5f +/- %-5.3f         %-10.5f +/- %-5.3f",
                  line.get(0),
                  line.get(1),
                  Double.parseDouble(line.get(2)),
                  Double.parseDouble(line.get(3)),
                  Double.parseDouble(line.get(4)),
                  Double.parseDouble(line.get(5)),
                  Double.parseDouble(line.get(6)),
                  Double.parseDouble(line.get(7))));
        } else {
          logger.info(format("  %-10s %-10s %-10.5f +/- %-5.3f",
                  line.get(0),
                  line.get(1),
                  Double.parseDouble(line.get(2)),
                  Double.parseDouble(line.get(3))));
        }
        line.clear();
      }
    }
  }

  private static ArrayList<Double> getZeroOneDeltas(int i, ExtendedSystem esv,
                                                    ForceFieldEnergy forceFieldEnergy, double[] x)
  {
      ArrayList<Double> deltaU = new ArrayList<Double>();
      Residue res = esv.getExtendedResidueList().get(i);
      double titrationState = esv.getTitrationLambda(res);
      double tautomerState = esv.getTautomerLambda(res);

      if (esv.getTautomerizingResidueList().contains(res)) {
        esv.setTautomerLambda(res, 1, false);

        esv.setTitrationLambda(res, 0, false);
        double zeroEnergy = forceFieldEnergy.energy(x, false);

        esv.setTitrationLambda(res, 1, false);
        double oneEnergy = forceFieldEnergy.energy(x, false);

        esv.setTitrationLambda(res, titrationState, false);

        deltaU.add(oneEnergy - zeroEnergy);
        esv.setTautomerLambda(res, 0, false);
      }

      esv.setTitrationLambda(res, 0, false);
      double zeroEnergy = forceFieldEnergy.energy(x, false);

      esv.setTitrationLambda(res, 1, false);
      double oneEnergy = forceFieldEnergy.energy(x, false);

      esv.setTitrationLambda(res, titrationState, false);

      deltaU.add(oneEnergy - zeroEnergy);
      if (esv.getTautomerizingResidueList().contains(res)) {
        esv.setTautomerLambda(res, tautomerState, false);
      }

      return deltaU;
    }

  private static double average(ArrayList<Double> list) {
    double sum = 0.0;
    for (Double d : list) {
      sum += d;
    }
    return sum / list.size();
  }

  private static ArrayList<Double> mult(double a, ArrayList<Double> u) {
        ArrayList<Double> result = new ArrayList<Double>();
        for (Double d : u) {
            result.add(a * d);
        }
      return result;
    }

  private static ArrayList<Double> mult(ArrayList<Double> v, ArrayList<Double> u) {
    if (v.size() != u.size()) {
      throw new IllegalArgumentException("Vector sizes must be equal.");
    }
    ArrayList<Double> result = new ArrayList<Double>();
    for (int i = 0; i < v.size(); i++) {
      result.add(v.get(i) * u.get(i));
    }
    return result;
  }

  private static ArrayList<Double> subtract(double a, ArrayList<Double> u) {
    ArrayList<Double> result = new ArrayList<Double>();
    for (Double d : u) {
      result.add(a - d);
    }
    return result;
  }

  private static ArrayList<Double> exp(ArrayList<Double> u) {
    ArrayList<Double> result = new ArrayList<Double>();
    for (Double d : u) {
      result.add(Math.exp(d));
    }
    return result;
  }

  private static ArrayList<Double> div(ArrayList<Double> a, ArrayList<Double> b){
    if (a.size() != b.size()) {
      throw new IllegalArgumentException("Vector sizes must be equal.");
    }
    ArrayList<Double> result = new ArrayList<Double>();
    for (int i = 0; i < a.size(); i++) {
      result.add(a.get(i) / b.get(i));
    }
    return result;
  }
}
