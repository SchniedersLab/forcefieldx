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
package ffx.xray.commands;

import edu.rit.pj.Comm;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.ManyBodyOptions;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.algorithms.optimize.TitrationManyBody;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementEnergy;
import ffx.xray.cli.XrayOptions;
import ffx.xray.refine.RefinementMode;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static java.lang.String.format;

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "xray.ManyBody")
public class ManyBody extends AlgorithmsCommand {

  @Mixin
  private XrayOptions xrayOptions;

  @Mixin
  private ManyBodyOptions manyBodyOptions;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames;
  private MolecularAssembly[] molecularAssemblies;
  private DiffractionData diffractionData;
  private double initialTargetEnergy;
  private double finalTargetEnergy;

  /**
   * ManyBody constructor.
   */
  public ManyBody() {
    super();
  }

  /**
   * ManyBody constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public ManyBody(String[] args) {
    super(args);
  }

  /**
   * ManyBody constructor.
   *
   * @param binding The Binding to use.
   */
  public ManyBody(FFXBinding binding) {
    super(binding);
  }

  @Override
  public ManyBody run() {

    if (!init()) {
      return this;
    }

    xrayOptions.init();

    // Atomic clashes are expected and will be handled using direct induced dipoles.
    System.setProperty("sor-scf-fallback", "false");
    System.setProperty("direct-scf-fallback", "true");

    // This flag is for ForceFieldEnergyOpenMM and must be set before reading files.
    // It enforces that all torsions include a Fourier series with 6 terms.
    // Otherwise, during titration the number of terms for each torsion may change and
    // cause updateParametersInContext to throw an exception.
    // Note that OpenMM is not usually used for crystals (it doesn't handle space groups).
    double titrationPH = manyBodyOptions.getTitrationPH();
    if (titrationPH > 0) {
      System.setProperty("manybody-titration", "true");
    }

    // Many-Body expansion of the X-ray target converges much more quickly with the NEA.
    String nea = System.getProperty("native-environment-approximation", "true");
    System.setProperty("native-environment-approximation", nea);

    String filename;
    if (filenames != null && !filenames.isEmpty()) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.getFirst());
      activeAssembly = molecularAssemblies[0];
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      molecularAssemblies = new MolecularAssembly[]{activeAssembly};
    }

    // Update the active filename
    filename = activeAssembly.getFile().getAbsolutePath();

    CompositeConfiguration properties = activeAssembly.getProperties();
    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);

    // The refinement mode must be coordinates.
    if (xrayOptions.refinementMode != RefinementMode.COORDINATES) {
      logger.info(" Refinement mode set to COORDINATES.");
      xrayOptions.refinementMode = RefinementMode.COORDINATES;
    }

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly);
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.");
      return this;
    }

    // Handle rotamer optimization with titration.
    TitrationManyBody titrationManyBody = null;
    if (titrationPH > 0) {
      logger.info(format("\n Adding titration hydrogen to: %s\n", filenames.getFirst()));
      List<Integer> resNumberList = new ArrayList<>();
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber());
      }

      // Create a new MolecularAssembly with additional protons and update the ForceFieldEnergy
      titrationManyBody = new TitrationManyBody(filenames.getFirst(), activeAssembly.getForceField(),
          resNumberList, titrationPH, manyBodyOptions);
      molecularAssemblies = titrationManyBody.getProtonatedAssemblies();
      activeAssembly = molecularAssemblies[0];
    }

    // Combine script flags (in parseResult) with properties.
    xrayOptions.setProperties(parseResult, properties);

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties);
    RefinementEnergy refinementEnergy = xrayOptions.toXrayEnergy(diffractionData);
    refinementEnergy.setScaling(null);

    boolean isTitrating = false;
    Set<Atom> excludeAtoms = new HashSet<>();

    for (int i = 0; i < molecularAssemblies.length; i++) {

      // Only optimization of the root location is currently tested.
      if (i > 0) {
        break;
      }
      activeAssembly = molecularAssemblies[i];
      activeAssembly.setFile(new File(filenames.getFirst()));


      // Save current B-factors
      Atom[] atoms = activeAssembly.getAtomArray();
      double[] bfactors = new double[atoms.length];
      double averageBFactor = 0;
      for (int j = 0; j < atoms.length; j++) {
        double bfactor = atoms[j].getTempFactor();
        bfactors[j] = bfactor;
        averageBFactor += bfactor;
      }
      // Set each atom to use the average b-factor to avoid a bias toward the initial rotamer.
      averageBFactor /= atoms.length;
      for (Atom atom : atoms) {
        atom.setTempFactor(averageBFactor);
      }

      diffractionData.scaleBulkFit();
      diffractionData.printStats();

      RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, refinementEnergy, algorithmListener);
      manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly);

      double[] x = new double[refinementEnergy.getNumberOfVariables()];
      x = refinementEnergy.getCoordinates(x);
      initialTargetEnergy = refinementEnergy.energy(x, true);
      logger.info(format("\n Initial target energy: %16.8f ", initialTargetEnergy));

      List<Residue> residueList = rotamerOptimization.getResidues();

      // For molecular assemblies other than the root assembly, only optimize alternate location residues.
      if (i > 0) {
        Character altLoc = activeAssembly.getAlternateLocation();
        List<Residue> altLocResidues = new ArrayList<>();
        for (Residue r : residues) {
          if (r.conatainsAltLoc(altLoc)) {
            altLocResidues.add(r);
          }
        }
        residueList = altLocResidues;
        rotamerOptimization.setResidues(altLocResidues);
      }

      RotamerLibrary.measureRotamers(residueList, false);
      rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()));

      int[] optimalRotamers = rotamerOptimization.getOptimumRotamers();

      if (titrationPH > 0) {
        isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, residueList);
      }
      logger.info(" Final Minimum Energy");
      // Get final parameters and compute the target function.
      x = refinementEnergy.getCoordinates(x);
      finalTargetEnergy = refinementEnergy.energy(x, true);

      if (isTitrating) {
        double phBias = rotamerOptimization.getEnergyExpansion().getTotalRotamerPhBias(residueList, optimalRotamers, titrationPH, manyBodyOptions.getPHRestraint());
        logger.info(format("\n  Rotamer pH Bias      %16.8f", phBias));
        logger.info(format("  Xray Target with Bias%16.8f\n", phBias + finalTargetEnergy));
      } else {
        logger.info(format("\n  Final Target Energy  %16.8f\n", finalTargetEnergy));
      }

      // Revert to the saved B-factors.
      for (int j = 0; j < atoms.length; j++) {
        atoms[j].setTempFactor(bfactors[j]);
      }
      diffractionData.scaleBulkFit();
      diffractionData.printStats();
    }

    if (Comm.world().rank() == 0) {
      properties.setProperty("standardizeAtomNames", "false");
      File modelFile = saveDirFile(activeAssembly.getFile());
      PDBFilter pdbFilter = new PDBFilter(modelFile, List.of(molecularAssemblies), activeAssembly.getForceField(), properties);
      if (titrationPH > 0) {
        String remark = format("Titration pH: %6.3f", titrationPH);
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true, new String[]{remark})) {
          logger.info(format(" Save failed for %s", activeAssembly));
        }
      } else {
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
          logger.info(format(" Save failed for %s", activeAssembly));
        }
      }
    }

    return this;
  }

  public double getInitialTargetEnergy() {
    return initialTargetEnergy;
  }

  public double getFinalTargetEnergy() {
    return finalTargetEnergy;
  }


  /**
   * Get the ManyBodyOptions.
   *
   * @return The ManyBodyOptions.
   */
  public ManyBodyOptions getManyBodyOptions() {
    return manyBodyOptions;
  }

  @Override
  public List<Potential> getPotentials() {
    return getPotentialsFromAssemblies(molecularAssemblies);
  }

  @Override
  public boolean destroyPotentials() {
    return diffractionData == null ? true : diffractionData.destroy();
  }
}
