//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3
import ffx.potential.bonded.*
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.utils.PotentialsUtils
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
@Command(description = " Run ManyBody algorithm on a system.", name = "ffxc ManyBody")
class ManyBody extends AlgorithmsScript {

  @Mixin
  ManyBodyOptions manyBody

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private String filename

  ForceFieldEnergy potentialEnergy
  boolean testing = false
  boolean monteCarloTesting = false

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

    String priorGKwarn = System.getProperty("gk-suppressWarnings")
    if (priorGKwarn == null || priorGKwarn.isEmpty()) {
      System.setProperty("gk-suppressWarnings", "true")
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)

    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    CompositeConfiguration properties = activeAssembly.getProperties()
    if (properties.getBoolean("standardizeAtomNames", false)) {
      renameAtomsToPDBStandard(activeAssembly)
    }

    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)

    // If rotamer optimization with titration, create new molecular assembly with additional protons
    ForceField forceField = activeAssembly.getForceField()
    potentialEnergy = activeAssembly.getPotentialEnergy()
    potentialEnergy.energy(false, true)

    List<String> resNumberList = new ArrayList<>()
    // List<Character> chainList = new ArrayList<>()
    List<Residue> residues
    if (manyBody.residueGroup.start > -1 || manyBody.residueGroup.all > -1) {
      residues = manyBody.getResidues(activeAssembly)
    } else {
      residues = activeAssembly.getResidueList()
    }

    if (residues.isEmpty()){
      logger.info("Residue list is empty")
    }

    for (Residue residue : residues) {
      resNumberList.add(String.valueOf(residue.getResidueNumber()))
    }

    activeAssembly = new MolecularAssembly(filename)

    activeAssembly.setForceField(forceField)

    File structureFile = new File(filename)
    PDBFilter protFilter = new PDBFilter(structureFile, activeAssembly, forceField, forceField.getProperties(), resNumberList)
    if(manyBody.group.titrationPH != 0){
      logger.info("\n Adding rotamer optimization with titration protons to : " + filename + "\n")
      protFilter.setRotamerTitration(true)
    }
    protFilter.readFile()
    protFilter.applyAtomProperties()
    activeAssembly.finalize(true, forceField)
    potentialEnergy = ForceFieldEnergy.energyFactory(activeAssembly)
    potentialEnergy.energy(false, true)
    activeAssembly.setFile(structureFile)

    RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, potentialEnergy, algorithmListener)

    testing = getTesting()
    if (testing) {
      rotamerOptimization.turnRotamerSingleEliminationOff()
      rotamerOptimization.turnRotamerPairEliminationOff()
    }

    if (monteCarloTesting) {
      rotamerOptimization.setMonteCarloTesting(true)
    }

    manyBody.initRotamerOptimization(rotamerOptimization, activeAssembly)

    List<Residue> residueList = rotamerOptimization.getResidues()

    boolean master = true
    if (Comm.world().size() > 1) {
      int rank = Comm.world().rank()
      if (rank != 0) {
        master = false
      }
    }

    algorithmFunctions.energy(activeAssembly)
    RotamerLibrary.measureRotamers(residueList, false)

    RotamerOptimization.Algorithm algorithm
    switch (manyBody.getAlgorithmNumber()) {
      case 1:
        algorithm = RotamerOptimization.Algorithm.INDEPENDENT
        break
      case 2:
        algorithm = RotamerOptimization.Algorithm.ALL
        break
      case 3:
        algorithm = RotamerOptimization.Algorithm.BRUTE_FORCE
        break
      case 4:
        algorithm = RotamerOptimization.Algorithm.WINDOW
        break
      case 5:
        algorithm = RotamerOptimization.Algorithm.BOX
        break
      default:
        throw new IllegalArgumentException(
            format(" Algorithm choice was %d, not in range 1-5!", manyBody.getAlgorithmNumber()))
    }
    rotamerOptimization.optimize(algorithm)

    boolean isTitrating = false
    Set<Atom> excludeAtoms = new HashSet<>()
    int[] optimalRotamers = rotamerOptimization.getOptimumRotamers()

    int i = 0
    for (Residue residue : residueList) {
      Rotamer rotamer = residue.getRotamers()[optimalRotamers[i++]]
      RotamerLibrary.applyRotamer(residue, rotamer)
      if (rotamer.isTitrating) {
        isTitrating = true
        AminoAcid3 aa3 = rotamer.aminoAcid3
        residue.setName(aa3.name())
        switch (aa3) {
          case AminoAcid3.HID:
            // No HE2
            Atom HE2 = residue.getAtomByName("HE2", true)
            excludeAtoms.add(HE2)
            break
          case AminoAcid3.HIE:
            // No HD1
            Atom HD1 = residue.getAtomByName("HD1", true)
            excludeAtoms.add(HD1)
            break
          case AminoAcid3.ASP:
            // No HD2
            Atom HD2 = residue.getAtomByName("HD2", true)
            excludeAtoms.add(HD2)
            break
          case AminoAcid3.GLU:
            // No HE2
            Atom HE2 = residue.getAtomByName("HE2", true)
            excludeAtoms.add(HE2)
            break
          case AminoAcid3.LYD:
            // No HZ3
            Atom HZ3 = residue.getAtomByName("HZ3", true)
            excludeAtoms.add(HZ3)
            break
          default:
            // Do nothing.
            break
        }
      }
    }

    if (master) {
      logger.info(" Final Minimum Energy")
      double energy = potentialEnergy.energy(false, true)
      if (isTitrating) {
        double phBias = rotamerOptimization.getEnergyExpansion().getTotalRotamerPhBias(residueList, optimalRotamers)
        logger.info(format("\n  Rotamer pH Bias    %16.8f", phBias))
        logger.info(format("  Potential with Bias%16.8f\n", phBias + energy))
      }

      // Prevent residues from being renamed based on the existence of hydrogen
      // atoms (i.e. hydrogen that excluded from being written out).
      properties.setProperty("standardizeAtomNames", "false")
      File modelFile = saveDirFile(activeAssembly.getFile())
      PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(), properties)
      if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
        logger.info(format(" Save failed for %s", activeAssembly))
      }
    }

    if (priorGKwarn == null) {
      System.clearProperty("gk-suppressWarnings")
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

  /**
   * Set method for the testing boolean. When true, the testing boolean will shut off all elimination criteria forcing either a monte carlo or brute force search over all permutations.
   * @param testing A boolean flag that turns off elimination criteria for testing purposes.
   */
  void setTesting(boolean testing) {
    this.testing = testing
  }

  /**
   * Get method for the testing boolean. When true, the testing boolean will shut off all elimination criteria forcing either a monte carlo or brute force search over all permutations.
   * @return testing A boolean flag that turns off elimination criteria for testing purposes.
   */
  boolean getTesting() {
    return testing
  }

  /**
   * Set to true when testing the monte carlo rotamer optimization algorithm. True will trigger the "set seed"
   * functionality of the pseudo-random number generator in the RotamerOptimization.java class to create a deterministic monte carlo algorithm.
   * @param bool True ONLY when a deterministic monte carlo approach is desired. False in all other cases.
   */
  void setMonteCarloTesting(boolean bool) {
    this.monteCarloTesting = bool
  }
}
