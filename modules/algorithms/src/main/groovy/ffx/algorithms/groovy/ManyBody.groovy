//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.algorithms.groovy

import edu.rit.pj.Comm
import ffx.algorithms.optimize.TitrationManyBody
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.LambdaInterface
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.parsers.PDBFilter
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard
import static java.lang.String.format

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Run ManyBody algorithm on a system.", name = "ManyBody")
class ManyBody extends AlgorithmsScript {

  @Mixin
  ManyBodyOptions manyBodyOptions

  @Mixin
  AlchemicalOptions alchemicalOptions

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private String filename

  ForceFieldEnergy potentialEnergy
  TitrationManyBody titrationManyBody

  /**
   * ManyBody Constructor.
   */
  ManyBody() {
    this(new Binding())
  }

  /**
   * ManyBody Constructor.
   * @param binding The Groovy Binding to use.
   */
  ManyBody(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  ManyBody run() {

    if (!init()) {
      return this
    }

    // This flag is for ForceFieldEnergyOpenMM and must be set before reading files.
    // It enforces that all torsions include a Fourier series with 6 terms.
    // Otherwise, during titration the number of terms for each torsion can change and
    // cause updateParametersInContext to throw an exception.
    double titrationPH = manyBodyOptions.getTitrationPH()

    if (manyBodyOptions.getTitration()) {
      System.setProperty("manybody-titration", "true")
    }

    // Turn on computation of lambda derivatives if softcore atoms exist.
    boolean lambdaTerm = alchemicalOptions.hasSoftcore()
    if (lambdaTerm) {
      // Turn on softcore van der Waals
      System.setProperty("lambdaterm", "true")
      // Turn off alchemical electrostatics
      System.setProperty("elec-lambdaterm", "false")
      // Turn on intra-molecular softcore
      System.setProperty("intramolecular-softcore", "true");
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)


    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    String listResidues = "";
    if(manyBodyOptions.getOnlyTitration() || manyBodyOptions.getOnlyProtons() ||
            manyBodyOptions.getInterestedResidue() != -1 && manyBodyOptions.getDistanceCutoff() != -1){
      listResidues = manyBodyOptions.selectDistanceResidues(activeAssembly.getResidueList(),
              manyBodyOptions.getInterestedResidue(),manyBodyOptions.getOnlyTitration(),
              manyBodyOptions.getOnlyProtons(),manyBodyOptions.getDistanceCutoff())
      manyBodyOptions.setListResidues(listResidues)
    }

    CompositeConfiguration properties = activeAssembly.getProperties()

    // Application of rotamers uses side-chain atom naming from the PDB.
    if (properties.getBoolean("standardizeAtomNames", false)) {
      renameAtomsToPDBStandard(activeAssembly)
    }

    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
    potentialEnergy = activeAssembly.getPotentialEnergy()

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly)
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.")
      return this
    }

    // Handle rotamer optimization with titration.
    if (manyBodyOptions.getTitration()) {
      logger.info("\n Adding titration hydrogen to : " + filename + "\n")

      // Collect residue numbers.
      List<Integer> resNumberList = new ArrayList<>()
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber())
      }

      // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
      titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
          resNumberList, titrationPH)
      MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
      setActiveAssembly(protonatedAssembly)
      potentialEnergy = protonatedAssembly.getPotentialEnergy()
    }

    if (lambdaTerm) {
      alchemicalOptions.setFirstSystemAlchemistry(activeAssembly)
      LambdaInterface lambdaInterface = (LambdaInterface) potentialEnergy
      double lambda = alchemicalOptions.getInitialLambda()
      logger.info(format(" Setting ManyBody softcore lambda to: %5.3f", lambda))
      lambdaInterface.setLambda(lambda)
    }

    RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly,
        potentialEnergy, algorithmListener)
    rotamerOptimization.setPHRestraint(manyBodyOptions.getPHRestraint())
    rotamerOptimization.setOnlyProtons(manyBodyOptions.getOnlyProtons())
    rotamerOptimization.setpH(titrationPH)

    manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)
    List<Residue> residueList = rotamerOptimization.getResidues()

    logger.info("\n Initial Potential Energy:")
    potentialEnergy.energy(false, true)

    logger.info("\n Initial Rotamer Torsion Angles:")
    RotamerLibrary.measureRotamers(residueList, false)

    // Run the optimization.
    rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()))

    boolean isTitrating = false
    Set<Atom> excludeAtoms = new HashSet<>()
    int[] optimalRotamers = rotamerOptimization.getOptimumRotamers()
    if (manyBodyOptions.getTitration()) {
      isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, residueList)
    }

    // Log the final result on rank 0.
    int rank = Comm.world().rank()
    if (rank == 0) {
      logger.info(" Final Minimum Energy")
      double energy = potentialEnergy.energy(false, true)
      if (isTitrating) {
        double phBias = rotamerOptimization.getEnergyExpansion().getTotalRotamerPhBias(residueList,
            optimalRotamers, titrationPH, manyBodyOptions.getPHRestraint())
        logger.info(format("\n  Rotamer pH Bias    %16.8f", phBias))
        logger.info(format("  Potential with Bias%16.8f\n", phBias + energy))
      }

      // Prevent residues from being renamed based on the existence of hydrogen
      // atoms (i.e. hydrogen that are excluded from being written out).
      properties.setProperty("standardizeAtomNames", "false")
      File modelFile = saveDirFile(activeAssembly.getFile())
      PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
          properties)
      if (manyBodyOptions.getTitration()) {
        String remark = format("Titration pH: %6.3f", titrationPH)
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true, remark)) {
          logger.info(format(" Save failed for %s", activeAssembly))
        }
      } else {
        if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
          logger.info(format(" Save failed for %s", activeAssembly))
        }
      }
    }

    if (manyBodyOptions.getTitration()) {
      System.clearProperty("manybody-titration")
    }

    return this
  }

  /**
   * Returns the potential energy of the active assembly. Used during testing assertions.
   * @return potentialEnergy Potential energy of the active assembly.
   */
  ForceFieldEnergy getPotential() {
    return potentialEnergy
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potentialEnergy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList((Potential) potentialEnergy)
    }
    return potentials
  }

}
