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
package ffx.potential.commands;

import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import java.util.logging.Level;

import static ffx.potential.utils.StructureMetrics.momentsOfInertia;
import static ffx.potential.utils.StructureMetrics.radiusOfGyration;
import static ffx.utilities.StringUtils.parseAtomRanges;
import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.getExtension;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * Compute the force field potential energy.
 *
 * Usage:
 *   ffxc Energy &lt;filename&gt;
 */
@Command(name = "Energy", description = " Compute the force field potential energy.")
public class Energy extends PotentialCommand {

  @Mixin
  private AtomSelectionOptions atomSelectionOptions = new AtomSelectionOptions();

  @Mixin
  private AlchemicalOptions alchemicalOptions = new AlchemicalOptions();

  @Mixin
  private TopologyOptions topologyOptions = new TopologyOptions();

  @Option(names = {"-m", "--moments"}, paramLabel = "false", defaultValue = "false",
      description = "Print out electrostatic moments.")
  private boolean moments = false;

  @Option(names = {"--rg", "--gyrate"}, paramLabel = "false", defaultValue = "false",
      description = "Print out the radius of gyration.")
  private boolean gyrate = false;

  @Option(names = {"--in", "--inertia"}, paramLabel = "false", defaultValue = "false",
      description = "Print out the moments of inertia.")
  private boolean inertia = false;

  @Option(names = {"-g", "--gradient"}, paramLabel = "false", defaultValue = "false",
      description = "Compute the atomic gradient as well as energy.")
  private boolean gradient = false;

  @Option(names = {"--fl", "--findLowest"}, paramLabel = "0", defaultValue = "0",
      description = "Return the n lowest energies from an ARC/PDB file.")
  private int fl = 0;

  @Option(names = {"-v", "--verbose"}, paramLabel = "false", defaultValue = "false",
      description = "Print out all energy components for each snapshot.")
  private boolean verbose = false;

  @Option(names = {"--dc", "--densityCutoff"}, paramLabel = "0.0", defaultValue = "0.0",
      description = "Create ARC file of structures above a specified density.")
  private double dCutoff = 0.0;

  @Option(names = {"--ec", "--energyCutoff"}, paramLabel = "0.0", defaultValue = "0.0",
      description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
  private double eCutoff = 0.0;

  @Option(names = {"--pd", "--printDihedral"}, paramLabel = "", defaultValue = "",
      description = "Atom indices to print dihedral angle values.")
  private String dihedralAtoms = "";

  @Option(names = {"--pb", "--printBondedTerms"}, paramLabel = "", defaultValue = "false",
      description = "Print all bonded energy terms.")
  private boolean printBondedTerms = false;

  @Parameters(arity = "1..4", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format.")
  private List<String> filenames = null;

  public double energy = 0.0;
  public ForceFieldEnergy forceFieldEnergy = null;
  private AssemblyState assemblyState = null;

  private Potential potential;
  MolecularAssembly[] topologies;

  public Energy() {
    super();
  }

  public Energy(FFXBinding binding) {
    super(binding);
  }

  public Energy(String[] args) {
    super(args);
  }

  @Override
  public Energy run() {
    if (!init()) {
      return this;
    }

    int numTopologies = topologyOptions.getNumberOfTopologies(filenames);
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];

    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    if (filenames == null || filenames.isEmpty()) {
      activeAssembly = getActiveAssembly(filenames);
      if (activeAssembly == null) {
        logger.info(helpString());
        return this;
      }
      filenames = new ArrayList<>();
      filenames.add(activeAssembly.getFile().getName());
      topologies[0] = alchemicalOptions.processFile(topologyOptions, activeAssembly, 0);
    } else {
      logger.info(format(" Initializing %d topologies...", numTopologies));
      for (int i = 0; i < numTopologies; i++) {
        topologies[i] = alchemicalOptions.openFile(potentialFunctions, topologyOptions, threadsPerTopology, filenames.get(i), i);
      }
      activeAssembly = topologies[0];
    }

    if (topologies.length == 1) {
      atomSelectionOptions.setActiveAtoms(topologies[0]);
    }

    StringBuilder sb = new StringBuilder("\n Calculating energy of ");
    potential = topologyOptions.assemblePotential(topologies, sb);
    logger.info(sb.toString());

    LambdaInterface linter = (potential instanceof LambdaInterface) ? (LambdaInterface) potential : null;
    if (linter != null) {
      linter.setLambda(alchemicalOptions.getInitialLambda());
    }

    String filename = activeAssembly.getFile().getAbsolutePath();
    logger.info("\n Running Energy on " + filename);

    forceFieldEnergy = activeAssembly.getPotentialEnergy();
    int nVars = potential.getNumberOfVariables();
    double[] x = new double[nVars];
    potential.getCoordinates(x);

    if (gradient) {
      double[] g = new double[nVars];
      int nAts = nVars / 3;
      energy = potential.energyAndGradient(x, g, true);
      logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"));
      for (int i = 0; i < nAts; i++) {
        int i3 = 3 * i;
        logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]));
      }
    } else {
      energy = potential.energy(x, true);
    }

    if (moments) {
      forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false);
    }

    if (gyrate) {
      double rg = radiusOfGyration(activeAssembly.getActiveAtomArray());
      logger.info(format(" Radius of gyration:           %10.5f A", rg));
    }

    if (inertia) {
      momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true);
    }

    ArrayList<Integer> unique = null;
    if (dihedralAtoms != null && !dihedralAtoms.isEmpty()) {
      unique = new ArrayList<>(parseAtomRanges("Dihedral Atoms", dihedralAtoms, activeAssembly.getAtomList().size()));
      printDihedral(activeAssembly, unique);
    }

    if (printBondedTerms) {
      forceFieldEnergy.logBondedTermsAndRestraints();
    }

    SystemFilter systemFilter = potentialFunctions.getFilter();
    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
      int numSnaps = fl;
      double lowestEnergy = energy;
      assemblyState = new AssemblyState(activeAssembly);
      int index = 1;
      PriorityQueue<StateContainer> lowestEnergyQueue = null;
      if (fl > 0) {
        lowestEnergyQueue = new PriorityQueue<>(numSnaps, Collections.reverseOrder());
        lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy));
      }

      int numModels = systemFilter.countNumModels();
      double[] densities = new double[numModels];
      double[] energies = new double[numModels];
      energies[0] = energy;

      while (systemFilter.readNext()) {
        index++;
        Crystal crystal = activeAssembly.getCrystal();
        densities[index - 1] = crystal.getDensity(activeAssembly.getMass());
        forceFieldEnergy.setCrystal(crystal);
        forceFieldEnergy.getCoordinates(x);
        if (verbose) {
          logger.info(format(" Snapshot %4d", index));
          if (!crystal.aperiodic()) {
            logger.info(format("\n Density:                                %6.3f (g/cc)", crystal.getDensity(activeAssembly.getMass())));
            if (logger.isLoggable(Level.FINE)) {
              logger.fine(crystal.toString());
            }
          }
          energy = forceFieldEnergy.energy(x, true);
        } else {
          energy = forceFieldEnergy.energy(x, false);
          logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy));
        }

        energies[index - 1] = energy;
        if (energy < lowestEnergy) {
          lowestEnergy = energy;
        }

        if (fl > 0) {
          StateContainer sc = new StateContainer(new AssemblyState(activeAssembly), energy);
          if (lowestEnergyQueue.size() < numSnaps) {
            lowestEnergyQueue.add(sc);
          } else {
            StateContainer worst = lowestEnergyQueue.peek();
            if (worst != null && energy < worst.getEnergy()) {
              lowestEnergyQueue.poll();
              lowestEnergyQueue.add(sc);
            }
          }
        }

        if (moments) {
          forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false);
        }
        if (gyrate) {
          double rg = radiusOfGyration(activeAssembly.getActiveAtomArray());
          logger.info(format(" Radius of gyration:          %10.5f A", rg));
        }
        if (inertia) {
          momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true);
        }
        if (dihedralAtoms != null && !dihedralAtoms.isEmpty()) {
          printDihedral(activeAssembly, unique);
        }
        if (printBondedTerms) {
          forceFieldEnergy.logBondedTermsAndRestraints();
        }
      }

      String dirString = getBaseDirString(filename);

      if ((eCutoff > 0.0 || dCutoff > 0.0) && numModels > 1) {
        systemFilter.readNext(true);
        activeAssembly = systemFilter.getActiveMolecularSystem();

        String name = getName(filename);
        String ext = getExtension(name);
        name = removeExtension(name);

        File saveFile;
        SystemFilter writeFilter;
        if (ext.toUpperCase().contains("XYZ")) {
          saveFile = new File(dirString + name + ".xyz");
          writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties());
          potentialFunctions.saveAsXYZ(activeAssembly, saveFile);
        } else if (ext.toUpperCase().contains("ARC")) {
          saveFile = new File(dirString + name + ".arc");
          saveFile = potentialFunctions.versionFile(saveFile);
          writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties());
          logger.info(" Saving to " + saveFile.getAbsolutePath());
          try { saveFile.createNewFile(); } catch (Exception e) { /* ignore */ }
        } else {
          saveFile = new File(dirString + name + ".pdb");
          saveFile = potentialFunctions.versionFile(saveFile);
          writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties());
          if (numModels > 1 && writeFilter instanceof PDBFilter) {
            ((PDBFilter) writeFilter).setModelNumbering(0);
          }
          try { saveFile.createNewFile(); } catch (Exception e) { /* ignore */ }
        }

        int structNum = 0;
        do {
          if (dCutoff > 0.0 && eCutoff > 0.0) {
            if (energies[structNum] < lowestEnergy + eCutoff && densities[structNum] > dCutoff) {
              if (systemFilter instanceof PDBFilter pdbFilter) {
                pdbFilter.writeFile(saveFile, true, false, false);
                try { Files.writeString(saveFile.toPath(), "ENDMDL\n", StandardOpenOption.APPEND); } catch (Exception e) { /* ignore */ }
              } else if (systemFilter instanceof XYZFilter) {
                writeFilter.writeFile(saveFile, true);
              }
            }
          } else if (dCutoff > 0.0) {
            if (densities[structNum] > dCutoff) {
              if (systemFilter instanceof PDBFilter pdbFilter) {
                pdbFilter.writeFile(saveFile, true, false, false);
                try { Files.writeString(saveFile.toPath(), "ENDMDL\n", StandardOpenOption.APPEND); } catch (Exception e) { /* ignore */ }
              } else if (systemFilter instanceof XYZFilter) {
                writeFilter.writeFile(saveFile, true);
              }
            }
          } else if (eCutoff > 0.0) {
            if (energies[structNum] < lowestEnergy + eCutoff) {
              if (systemFilter instanceof PDBFilter pdbFilter) {
                pdbFilter.writeFile(saveFile, true, false, false);
                try { Files.writeString(saveFile.toPath(), "ENDMDL\n", StandardOpenOption.APPEND); } catch (Exception e) { /* ignore */ }
              } else if (systemFilter instanceof XYZFilter) {
                writeFilter.writeFile(saveFile, true);
              }
            }
          }
          structNum++;
        } while (systemFilter.readNext());
        if (systemFilter instanceof PDBFilter) {
          try { Files.writeString(saveFile.toPath(), "END\n", StandardOpenOption.APPEND); } catch (Exception e) { /* ignore */ }
        }
      }

      if (fl > 0) {
        if (numSnaps > index) {
          logger.warning(format(" Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported", numSnaps, filenames, index, index));
          numSnaps = index;
        }

        String name = getName(filenames.get(0));
        logger.info(" Saving the " + numSnaps + " lowest energy snapshots to " + dirString + name);
        SystemFilter saveFilter = potentialFunctions.getFilter();
        File saveFile = potentialFunctions.versionFile(new File(dirString + name));
        int toSave = Math.min(numSnaps, lowestEnergyQueue.size());
        for (int i = 0; i < toSave; i++) {
          StateContainer savedState = lowestEnergyQueue.poll();
          AssemblyState finalAssembly = savedState.getState();
          finalAssembly.revertState();
          String remark = format("Potential Energy: %16.8f (kcal/mol)", savedState.getEnergy());
          logger.info(format(" Snapshot %d %s", i + 1, remark));
          saveFilter.writeFile(saveFile, true, new String[]{remark});
        }
      }
    }

    return this;
  }

  private static void printDihedral(MolecularAssembly molecularAssembly, ArrayList<Integer> indices) {
    if (indices != null && indices.size() == 4) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      Atom atom0 = atoms[indices.get(0)];
      Atom atom1 = atoms[indices.get(1)];
      Atom atom2 = atoms[indices.get(2)];
      Atom atom3 = atoms[indices.get(3)];
      try {
        for (Torsion torsion : atom0.getTorsions()) {
          if (torsion.compare(atom0, atom1, atom2, atom3)) {
            logger.info("\n Torsion: " + torsion);
            return;
          }
        }
      } catch (Exception e) {
        logger.info(" Exception during dihedral print " + e);
      }
      logger.info("\n No torsion between atoms:\n  " + atom0 + "\n  " + atom1 + "\n  " + atom2 + "\n  " + atom3);
    }
  }

  @Override
  public List<Potential> getPotentials() {
    if (forceFieldEnergy == null) {
      return Collections.emptyList();
    } else {
      return Collections.singletonList(forceFieldEnergy);
    }
  }

  public ForceFieldEnergy getForceFieldEnergy() {
    return forceFieldEnergy;
  }

  public double getEnergy() {
    return energy;
  }

  private static class StateContainer implements Comparable<StateContainer> {
    private final AssemblyState state;
    private final double e;
    StateContainer(AssemblyState state, double e) { this.state = state; this.e = e; }
    AssemblyState getState() { return state; }
    double getEnergy() { return e; }
    @Override
    public int compareTo(StateContainer o) { return Double.compare(e, o.getEnergy()); }
  }
}