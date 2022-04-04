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
package ffx.realspace.groovy

import edu.rit.pj.Comm
import ffx.algorithms.TitrationManyBody
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFilter
import ffx.realspace.cli.RealSpaceOptions
import ffx.xray.RefinementEnergy
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static java.lang.String.format

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "ffxc realspace.ManyBody")
class ManyBody extends AlgorithmsScript {

  @Mixin
  private RealSpaceOptions realSpace

  @Mixin
  ManyBodyOptions manyBody

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

    if (manyBody.group.titrationPH != 0) {
      System.setProperty("manybody-titration", "true")
    }

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      activeAssembly = algorithmFunctions.open(filenames.get(0))
      modelFilename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }
    MolecularAssembly[] assemblies = [activeAssembly] as MolecularAssembly[]

    CompositeConfiguration properties = activeAssembly.getProperties()
    if (!properties.containsKey("gk-suppressWarnings")) {
      properties.setProperty("gk-suppressWarnings", "true")
    }
    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
    ForceField forceField = activeAssembly.getForceField()
    potentialEnergy = activeAssembly.getPotentialEnergy()

    // Collect residues to optimize.
    List<Integer> resNumberList = new ArrayList<>()
    List<Residue> residues
    if (manyBody.residueGroup.start > -1 || manyBody.residueGroup.all > -1) {
      residues = manyBody.getResidues(activeAssembly)
    } else {
      residues = activeAssembly.getResidueList()
    }

    if (residues.isEmpty()) {
      logger.info("Residue list is empty")
      return
    }

    for (Residue residue : residues) {
      resNumberList.add(residue.getResidueNumber())
    }

    // If rotamer optimization with titration, create new molecular assembly with additional protons
    if (manyBody.group.titrationPH != 0) {
      logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")
      titrationManyBody = new TitrationManyBody(filenames.get(0),forceField,resNumberList,manyBody.group.titrationPH)
      MolecularAssembly prot = titrationManyBody.getProtonatedAssembly()
      setActiveAssembly(prot)
      potentialEnergy = prot.getPotentialEnergy()
      assemblies = [activeAssembly] as MolecularAssembly[]
    }

    refinementEnergy = realSpace.toRealSpaceEnergy(filenames, assemblies, algorithmFunctions)
    RotamerOptimization rotamerOptimization = new RotamerOptimization(
        activeAssembly, refinementEnergy, algorithmListener)

    manyBody.initRotamerOptimization(rotamerOptimization, activeAssembly)

    double[] x = new double[refinementEnergy.getNumberOfVariables()]
    x = refinementEnergy.getCoordinates(x)
    refinementEnergy.energy(x, true)

    List<Residue> residueList = rotamerOptimization.getResidues()
    RotamerLibrary.measureRotamers(residueList, false)

    if (manyBody.algorithm == 1) {
      rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT)
    } else if (manyBody.algorithm == 2) {
      rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL)
    } else if (manyBody.algorithm == 3) {
      rotamerOptimization.optimize(RotamerOptimization.Algorithm.BRUTE_FORCE)
    } else if (manyBody.algorithm == 4) {
      rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW)
    } else if (manyBody.algorithm == 5) {
      rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX)
    }

    boolean master = true
    if (Comm.world().size() > 1) {
      int rank = Comm.world().rank()
      if (rank != 0) {
        master = false
      }
    }

    boolean isTitrating = false
    Set<Atom> excludeAtoms = new HashSet<>()
    int[] optimalRotamers = rotamerOptimization.getOptimumRotamers()
    excludeAtoms = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, residueList)

    if (master) {
      logger.info(" Final Minimum Energy")
      algorithmFunctions.energy(activeAssembly)
      double energy = potentialEnergy.energy(false, true)
      if (isTitrating) {
        double phBias = rotamerOptimization.getEnergyExpansion().getTotalRotamerPhBias(residueList,
                optimalRotamers)
        logger.info(format("\n  Rotamer pH Bias    %16.8f", phBias))
        logger.info(format("  Potential with Bias%16.8f\n", phBias + energy))
      }
      String ext = FilenameUtils.getExtension(modelFilename)
      modelFilename = FilenameUtils.removeExtension(modelFilename)
      if (ext.toUpperCase().contains("XYZ")) {
        algorithmFunctions.saveAsXYZ(assemblies[0], new File(modelFilename + ".xyz"))
      } else {
        //algorithmFunctions.saveAsPDB(assemblies, new File(modelFilename + ".pdb"))
        properties.setProperty("standardizeAtomNames", "false")
        File modelFile = saveDirFile(activeAssembly.getFile())
        PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
                properties)
        if (manyBody.group.titrationPH != 0){
          String remark = "Titration pH:   " + manyBody.group.titrationPH.toString()
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
