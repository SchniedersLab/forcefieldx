package ffx.potential.commands;

import static java.lang.String.format;

import com.google.common.collect.MinMaxPriorityQueue;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.utils.StructureMetrics;
import ffx.utilities.FFXContext;
import java.io.File;
import java.util.Collections;
import java.util.List;
import org.apache.commons.io.FilenameUtils;
//import org.codehaus.groovy.runtime.ResourceGroovyMethods;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc Energy")
public class Energy extends PotentialCommand {

  @Mixin
  private AtomSelectionOptions atomSelectionOptions;
  /**
   * -m or --moments print out electrostatic moments.
   */
  @Option(names = {"-m",
      "--moments"}, paramLabel = "false", defaultValue = "false", description = "Print out electrostatic moments.")
  private boolean moments = false;
  /**
   * --rg or --gyrate Print out the radius of gyration.
   */
  @Option(names = {"--rg",
      "--gyrate"}, paramLabel = "false", defaultValue = "false", description = "Print out the radius of gyration.")
  private boolean gyrate = false;
  /**
   * --in or --inertia Print out the moments of inertia.
   */
  @Option(names = {"--in",
      "--inertia"}, paramLabel = "false", defaultValue = "false", description = "Print out the moments of inertia.")
  private boolean inertia = false;
  /**
   * -g or --gradient to print out gradients.
   */
  @Option(names = {"-g",
      "--gradient"}, paramLabel = "false", defaultValue = "false", description = "Compute the atomic gradient as well as energy.")
  private boolean gradient = false;
  /**
   * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
   */
  @Option(names = {"--fl",
      "--findLowest"}, paramLabel = "0", defaultValue = "0", description = "Return the n lowest energies from an ARC/PDB file.")
  private int fl = 0;
  /**
   * -v or --verbose enables printing out all energy components for multi-snapshot files ( the first
   * snapshot is always printed verbosely).
   */
  @Option(names = {"-v",
      "--verbose"}, paramLabel = "false", defaultValue = "false", description = "Print out all energy components for each snapshot.")
  private boolean verbose = false;
  /**
   * --dc or --densityCutoff Collect structures above a specified density.
   */
  @Option(names = {"--dc",
      "--densityCutoff"}, paramLabel = "0.0", defaultValue = "0.0", description = "Create ARC file of structures above a specified density.")
  private double dCutoff = 0.0;
  /**
   * --ec or --energyCutOff Collect structures below a specified energy range from the minimum
   * energy.
   */
  @Option(names = {"--ec",
      "--energyCutoff"}, paramLabel = "0.0", defaultValue = "0.0", description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
  private double eCutoff = 0.0;
  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file", description = "The atomic coordinate file in PDB or XYZ format.")
  private String filename = null;
  public double energy = 0.0;
  public ForceFieldEnergy forceFieldEnergy = null;
  private AssemblyState assemblyState = null;

  /**
   * Energy constructor.
   */
  public Energy() {
    this(new FFXContext());
  }

  /**
   * Energy constructor.
   *
   * @param ffxContext The FFXContext to use.
   */
  public Energy(FFXContext ffxContext) {
    super(ffxContext);
  }

  /**
   * Execute the script.
   */
  public Energy run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    setActiveAssembly(getActiveAssembly(filename));
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    logger.info("\n Running Energy on " + filename);

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly);

    forceFieldEnergy = activeAssembly.getPotentialEnergy();
    int nVars = forceFieldEnergy.getNumberOfVariables();
    double[] x = new double[nVars];
    forceFieldEnergy.getCoordinates(x);

    if (gradient) {
      double[] g = new double[nVars];
      int nAts = (nVars / 3);
      energy = forceFieldEnergy.energyAndGradient(x, g, true);
      logger.info("    Atom       X, Y and Z Gradient Components (kcal/mol/A)");
      for (int i = 0; i < nAts; i++) {
        int i3 = 3 * i;
        logger.info(String.format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]));
      }
    } else {
      energy = forceFieldEnergy.energy(x, true);
    }

    if (moments) {
      forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false);
    }

    if (gyrate) {
      double rg = StructureMetrics.radiusOfGyration(activeAssembly.getActiveAtomArray());
      logger.info(String.format(" Radius of gyration:           %10.5f A", rg));
    }

    if (inertia) {
      double[][] inertiaValue = StructureMetrics.momentsOfInertia(
          activeAssembly.getActiveAtomArray(), false, true, true);
    }

    SystemFilter systemFilter = potentialFunctions.getFilter();
    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {

      int numSnaps = fl;
      double lowestEnergy = energy;
      assemblyState = new AssemblyState(activeAssembly);
      int index = 1;

      // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
      MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = null;
      if (fl > 0) {
        lowestEnergyQueue = MinMaxPriorityQueue.maximumSize(numSnaps).create();
        lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy));
      }

      int numModels = systemFilter.countNumModels();
      //Store densities in ordered encountered (used in density cutoff implementation).
      double[] densities = new double[numModels];
      //Store energies in ordered encountered (used in energy cutoff implementation).
      double[] energies = new double[numModels];
      energies[0] = energy;

      while (systemFilter.readNext()) {
        index = index++;
        Crystal crystal = activeAssembly.getCrystal();
        densities[index - 1] = crystal.getDensity(activeAssembly.getMass());
        forceFieldEnergy.setCrystal(crystal);
        forceFieldEnergy.getCoordinates(x);
        if (verbose) {
          logger.info(String.format(" Snapshot %4d", index));
          if (!crystal.aperiodic()) {
            logger.info(String.format("\n Density:                                %6.3f (g/cc)",
                crystal.getDensity(activeAssembly.getMass())));
          }

          energy = forceFieldEnergy.energy(x, true);
        } else {
          energy = forceFieldEnergy.energy(x, false);
          logger.info(String.format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy));
        }

        energies[index - 1] = energy;

        //Update lowest encountered energy
        if (energy < lowestEnergy) {
          lowestEnergy = energy;
        }

        if (fl > 0) {
          lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), energy));
        }

        if (moments) {
          forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false);
        }

        if (gyrate) {
          double rg = StructureMetrics.radiusOfGyration(activeAssembly.getActiveAtomArray());
          logger.info(String.format(" Radius of gyration:          %10.5f A", rg));
        }

        if (inertia) {
          double[][] inertiaValue = StructureMetrics.momentsOfInertia(
              activeAssembly.getActiveAtomArray(), false, true, true);
        }

      }

      // Use the current base directory, or update if necessary based on the given filename.
      String dirString = getBaseDirString(filename);

      // If cutoffs have been selected create an ARC or PDB to store structures that satisfy cutoff.
      try {
        if ((eCutoff > 0.0 || dCutoff > 0.0) && numModels > 1) {
          systemFilter.readNext(true);
          setActiveAssembly(systemFilter.getActiveMolecularSystem());

          String name = FilenameUtils.getName(filename);
          String ext = FilenameUtils.getExtension(name);
          name = FilenameUtils.removeExtension(name);

          File saveFile;
          SystemFilter writeFilter;
          if (ext.toUpperCase().contains("XYZ")) {
            saveFile = new File(dirString + name + ".xyz");
            writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                activeAssembly.getProperties());
            potentialFunctions.saveAsXYZ(activeAssembly, saveFile);
          } else if (ext.toUpperCase().contains("ARC")) {
            saveFile = new File(dirString + name + ".arc");
            saveFile = potentialFunctions.versionFile(saveFile);
            writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                activeAssembly.getProperties());
            logger.info(" Saving to " + saveFile.getAbsolutePath());
            saveFile.createNewFile();
          } else {
            saveFile = new File(dirString + name + ".pdb");
            saveFile = potentialFunctions.versionFile(saveFile);
            writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                activeAssembly.getProperties());
            ((PDBFilter) writeFilter).setModelNumbering(0);
            saveFile.createNewFile();
          }

          // Determine if each structure meets the cutoff condition
          if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
            int structNum = 0;
            if (eCutoff > 0.0) {
              logger.info(
                  format("Lowest Energy of: %16.4f\n Saving structures with energy below: %16.4f",
                      lowestEnergy, lowestEnergy + eCutoff));
            }

            do {
              if (dCutoff > 0.0 && eCutoff > 0.0) {
                if (energies[structNum] < lowestEnergy + eCutoff && densities[structNum] > dCutoff) {
                  if (systemFilter instanceof PDBFilter) {
                    PDBFilter pdbFilter = (PDBFilter) systemFilter;
                    pdbFilter.writeFile(saveFile, true, false, false);
                    // ResourceGroovyMethods.append(saveFile, "ENDMDL\n");
                  } else if (systemFilter instanceof XYZFilter) {
                    writeFilter.writeFile(saveFile, true);
                  }
                }
              } else if (dCutoff > 0.0) {
                if (densities[structNum] > dCutoff) {
                  if (systemFilter instanceof PDBFilter) {
                    PDBFilter pdbFilter = (PDBFilter) systemFilter;
                    pdbFilter.writeFile(saveFile, true, false, false);
                    // ResourceGroovyMethods.append(saveFile, "ENDMDL\n");
                  } else if (systemFilter instanceof XYZFilter) {
                    writeFilter.writeFile(saveFile, true);
                  }
                }
              } else if (eCutoff > 0.0) {
                if (energies[structNum] < lowestEnergy + eCutoff) {
                  if (systemFilter instanceof PDBFilter) {
                    PDBFilter pdbFilter = (PDBFilter) systemFilter;
                    pdbFilter.writeFile(saveFile, true, false, false);
                    // ResourceGroovyMethods.append(saveFile, "ENDMDL\n");
                  } else if (systemFilter instanceof XYZFilter) {
                    writeFilter.writeFile(saveFile, true);
                  }
                }
              }
              structNum++;
            } while (systemFilter.readNext());
            if (systemFilter instanceof PDBFilter) {
              // ResourceGroovyMethods.append(saveFile, "END\n");
            }
          }
        }
      } catch (Exception e) {
        logger.info(" Exception evaluating cutoffs:\n " + e);
      }

      if (fl > 0) {
        if (numSnaps > index) {
          logger.warning(String.format(
              " Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported",
              numSnaps, filename, index, index));
          numSnaps = index;
        }

        String name = FilenameUtils.getName(filename);

        for (int i = 0; i < numSnaps - 1; i++) {
          StateContainer savedState = lowestEnergyQueue.removeLast();
          AssemblyState finalAssembly = savedState.getState();
          finalAssembly.revertState();
          double finalEnergy = savedState.getEnergy();
          logger.info(
              String.format(" The potential energy found is %16.8f (kcal/mol)", finalEnergy));
          File saveFile = potentialFunctions.versionFile(new File(dirString + name));
          MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly();
          potentialFunctions.saveAsPDB(molecularAssembly, saveFile);
        }

        StateContainer savedState = lowestEnergyQueue.removeLast();
        AssemblyState lowestAssembly = savedState.getState();
        lowestEnergy = savedState.getEnergy();

        assemblyState.revertState();
        logger.info(
            String.format(" The lowest potential energy found is %16.8f (kcal/mol)", lowestEnergy));

        // Prints our final energy (which will be the lowest energy
        File saveFile = potentialFunctions.versionFile(new File(dirString + name));
        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly();
        potentialFunctions.saveAsPDB(molecularAssembly, saveFile);
      }

    }

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials;
    if (forceFieldEnergy == null) {
      potentials = Collections.emptyList();
    } else {
      potentials = Collections.singletonList((Potential) forceFieldEnergy);
    }

    return potentials;
  }

  /**
   * This entry point is being used to test GraalVM ahead-of-time compilation.
   *
   * @param args Command line arguments.
   */
  public static void main(String... args) {
    System.setProperty("java.awt.headless", "true");
    System.setProperty("j3d.rend", "noop");
    FFXContext binding = new FFXContext(args);
    Energy energyScript = new Energy(binding);
    energyScript.run();
    System.exit(0);
  }

  public AtomSelectionOptions getAtomSelectionOptions() {
    return atomSelectionOptions;
  }

  public void setAtomSelectionOptions(AtomSelectionOptions atomSelectionOptions) {
    this.atomSelectionOptions = atomSelectionOptions;
  }

  private class StateContainer implements Comparable<StateContainer> {

    public StateContainer(AssemblyState state, double e) {
      this.state = state;
      this.e = e;
    }

    public AssemblyState getState() {
      return state;
    }

    public double getEnergy() {
      return e;
    }

    @Override
    public int compareTo(StateContainer o) {
      return Double.compare(e, o.getEnergy());
    }

    private final AssemblyState state;
    private final double e;
  }
}
