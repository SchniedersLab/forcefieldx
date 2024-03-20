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
package ffx.realspace.groovy

import edu.rit.pj.Comm
import ffx.algorithms.optimize.TitrationManyBody
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.parsers.PDBFilter
import ffx.realspace.cli.RealSpaceOptions
import ffx.xray.RefinementEnergy
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "realspace.ManyBody")
class ManyBody extends AlgorithmsScript {

  @Mixin
  private RealSpaceOptions realSpaceOptions

  @Mixin
  ManyBodyOptions manyBodyOptions

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames
  private RefinementEnergy refinementEnergy

  ForceFieldEnergy potentialEnergy
  TitrationManyBody titrationManyBody

  /**
   * ManyBody constructor.
   */
  ManyBody() {
    this(new Binding())
  }

  /**
   * ManyBody constructor.
   * @param binding The Groovy Binding to use.
   */
  ManyBody(Binding binding) {
    super(binding)
  }

  @Override
  ManyBody run() {

    if (!init()) {
      return this
    }

    // This flag is for ForceFieldEnergyOpenMM and must be set before reading files.
    // It enforces that all torsions include a Fourier series with 6 terms.
    // Otherwise, during titration the number of terms for each torsion may change and
    // causing updateParametersInContext to throw an exception.
    // Note that OpenMM is not usually used for crystals (it doesn't handle space groups).
    double titrationPH = manyBodyOptions.getTitrationPH()
    if (titrationPH > 0) {
      System.setProperty("manybody-titration", "true")
    }

    String filename
    if (filenames != null && filenames.size() > 0) {
      activeAssembly = algorithmFunctions.open(filenames.get(0))
      filename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      filename = activeAssembly.getFile().getAbsolutePath()
    }
    MolecularAssembly[] molecularAssemblies = [activeAssembly] as MolecularAssembly[]

    CompositeConfiguration properties = activeAssembly.getProperties()
    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
    potentialEnergy = activeAssembly.getPotentialEnergy()

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly)
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.")
      return this
    }

    // Handle rotamer optimization with titration.
    if (titrationPH > 0) {
      logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")

      List<Integer> resNumberList = new ArrayList<>()
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber())
      }

      // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
      titrationManyBody = new TitrationManyBody(filenames.get(0), activeAssembly.getForceField(),
          resNumberList, titrationPH)
      MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
      setActiveAssembly(protonatedAssembly)
      potentialEnergy = protonatedAssembly.getPotentialEnergy()
      molecularAssemblies = [activeAssembly] as MolecularAssembly[]
    }

    refinementEnergy = realSpaceOptions.toRealSpaceEnergy(filenames, molecularAssemblies)

    RotamerOptimization rotamerOptimization = new RotamerOptimization(
        activeAssembly, refinementEnergy, algorithmListener)
    manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)

    double[] x = new double[refinementEnergy.getNumberOfVariables()]
    x = refinementEnergy.getCoordinates(x)
    refinementEnergy.energy(x, true)

    List<Residue> residueList = rotamerOptimization.getResidues()
    RotamerLibrary.measureRotamers(residueList, false)

    rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()))

    boolean isTitrating = false
    Set<Atom> excludeAtoms = new HashSet<>()
    int[] optimalRotamers = rotamerOptimization.getOptimumRotamers()

    if (titrationPH > 0) {
      isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, residueList)
    }

    if (Comm.world().rank() == 0) {
      logger.info(" Final Minimum Energy")
      // Get final parameters and compute the target function.
      x = refinementEnergy.getCoordinates(x)
      double energy = refinementEnergy.energy(x, true)

      if (isTitrating) {
        double phBias = rotamerOptimization.getEnergyExpansion().getTotalRotamerPhBias(residueList,
            optimalRotamers)
        logger.info(format("\n  Rotamer pH Bias    %16.8f", phBias))
        logger.info(format("  Potential with Bias%16.8f\n", phBias + energy))
      } else {
        logger.info(format("\n  Real Space Target  %16.8f\n", energy))
      }
      String ext = FilenameUtils.getExtension(filename)
      filename = FilenameUtils.removeExtension(filename)
      if (ext.toUpperCase().contains("XYZ")) {
        algorithmFunctions.saveAsXYZ(molecularAssemblies[0], new File(filename + ".xyz"))
      } else {
        properties.setProperty("standardizeAtomNames", "false")
        File modelFile = saveDirFile(activeAssembly.getFile())
        PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly,
            activeAssembly.getForceField(),
            properties)
        if (titrationPH > 0) {
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
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    return refinementEnergy == null ? Collections.emptyList() :
        Collections.singletonList((Potential) refinementEnergy)
  }
}
