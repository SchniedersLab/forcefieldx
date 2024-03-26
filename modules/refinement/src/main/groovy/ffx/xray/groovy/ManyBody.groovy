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
package ffx.xray.groovy

import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.algorithms.optimize.manybody.EnergyExpansion
import ffx.algorithms.optimize.TitrationManyBody
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.parsers.PDBFilter
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.cli.XrayOptions
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "xray.ManyBody")
class ManyBody extends AlgorithmsScript {

  @Mixin
  XrayOptions xrayOptions

  @Mixin
  ManyBodyOptions manyBodyOptions

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames
  private RefinementEnergy refinementEnergy

  ForceFieldEnergy potentialEnergy
  private MolecularAssembly[] molecularAssemblies
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

    xrayOptions.init()

    // This flag is for ForceFieldEnergyOpenMM and must be set before reading files.
    // It enforces that all torsions include a Fourier series with 6 terms.
    // Otherwise, during titration the number of terms for each torsion may change and
    // causing updateParametersInContext to throw an exception.
    // Note that OpenMM is not usually used for crystals (it doesn't handle space groups).
    double titrationPH = manyBodyOptions.getTitrationPH()
    if (titrationPH > 0) {
      System.setProperty("manybody-titration", "true")
    }

    // Many-Body expansion of the X-ray target converges much more quickly with the NEA.
    String nea = System.getProperty("native-environment-approximation", "true")
    System.setProperty("native-environment-approximation", nea)

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = molecularAssemblies[0]
      logger.info(molecularAssemblies.length.toString())
      modelFilename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      molecularAssemblies = [activeAssembly]
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }

    CompositeConfiguration properties = activeAssembly.getProperties()
    activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
    potentialEnergy = activeAssembly.getPotentialEnergy()

    // The refinement mode must be coordinates.
    if (xrayOptions.refinementMode != RefinementMode.COORDINATES) {
      logger.info(" Refinement mode set to COORDINATES.")
      xrayOptions.refinementMode = RefinementMode.COORDINATES
    }

    // Collect residues to optimize.
    List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly)
    if (residues == null || residues.isEmpty()) {
      logger.info(" There are no residues in the active system to optimize.")
      return this
    }

    // Handle rotamer optimization with titration.
    if (titrationPH > 0) {
      logger.info(format("\n Adding titration hydrogen to: %s\n", filenames.get(0)))
      List<Integer> resNumberList = new ArrayList<>()
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber())
      }

      // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
      titrationManyBody = new TitrationManyBody(filenames.get(0), activeAssembly.getForceField(),
          resNumberList, titrationPH)
      MolecularAssembly[] protonatedAssemblies = titrationManyBody.getProtonatedAssemblies()
      setActiveAssembly(protonatedAssemblies[0])
      potentialEnergy = protonatedAssemblies[0].getPotentialEnergy()
      molecularAssemblies = protonatedAssemblies
    }

    // Load parsed X-ray properties.
    xrayOptions.setProperties(parseResult, properties)

    // Set up the diffraction data, which could be multiple files.
    DiffractionData diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties)
    refinementEnergy = xrayOptions.toXrayEnergy(diffractionData)
    refinementEnergy.setScaling(null)

    boolean isTitrating = false
    Set<Atom> excludeAtoms = new HashSet<>()
    for (MolecularAssembly currentMolecularAssembly : molecularAssemblies) {
      setActiveAssembly(currentMolecularAssembly)
      currentMolecularAssembly.setFile(filenames[0] as File)
      if (currentMolecularAssembly.getAtomList().get(0).getAltLoc() == 'A' && molecularAssemblies.length > 1) {
        for (int i = 0; i < molecularAssemblies[0].getAtomList().size(); i++) {
          molecularAssemblies[0].getAtomList().get(i).setOccupancy(1.0)
          molecularAssemblies[1].getAtomList().get(i).setOccupancy(0.0)
        }
        logger.info(" Occupancy of 1st Molecular Assembly Atoms: " + molecularAssemblies[0].getAtomList().get(0).getOccupancy())
        logger.info(" Occupancy of 2nd Molecular Assembly Atoms: " + molecularAssemblies[1].getAtomList().get(0).getOccupancy())

      } else if (currentMolecularAssembly.getAtomList().get(0).getAltLoc() == 'B' && molecularAssemblies.length > 1) {
        for (int i = 0; i < molecularAssemblies[0].getAtomList().size(); i++) {
          molecularAssemblies[0].getAtomList().get(i).setOccupancy(0.5)
          molecularAssemblies[1].getAtomList().get(i).setOccupancy(0.5)
        }
        logger.info(" Occupancy of 1st Molecular Assembly Atoms: " + molecularAssemblies[0].getAtomList().get(0).getOccupancy())
        logger.info(" Occupancy of 2nd Molecular Assembly Atoms: " + molecularAssemblies[1].getAtomList().get(0).getOccupancy())
      }

      RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, refinementEnergy, algorithmListener)
      manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)

      double[] x = new double[refinementEnergy.getNumberOfVariables()]
      x = refinementEnergy.getCoordinates(x)
      double e = refinementEnergy.energy(x, true)
      logger.info(format("\n Initial target energy: %16.8f ", e))

      List<Residue> residueList = rotamerOptimization.getResidues()
      RotamerLibrary.measureRotamers(residueList, false)

      rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()))

      int[] optimalRotamers = rotamerOptimization.getOptimumRotamers()


      if (titrationPH > 0) {
        isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, residueList)
      }
      logger.info(" Final Minimum Energy")
      // Get final parameters and compute the target function.
      x = refinementEnergy.getCoordinates(x)
      e = refinementEnergy.energy(x, true)

      if (isTitrating) {
        double phBias = EnergyExpansion.getTotalRotamerPhBias(residueList, optimalRotamers)
        logger.info(format("\n  Rotamer pH Bias      %16.8f", phBias))
        logger.info(format("  Xray Target with Bias%16.8f\n", phBias + e))
      } else {
        logger.info(format("\n  Xray Target          %16.8f\n", e))
      }
      diffractionData.scaleBulkFit()
      diffractionData.printStats()
    }

    if (molecularAssemblies.length > 1) {
      List<Residue> residueListA = molecularAssemblies[0].getResidueList()
      List<Residue> residueListB = molecularAssemblies[1].getResidueList()
      int firstRes = residueListA.get(0).getResidueNumber()
      for (Residue residue : residueListA) {
        int resNum = residue.getResidueNumber()
        List<Atom> atomList = residue.getAtomList()
        for (int i = 0; i < residue.getAtomList().size(); i++) {
          Atom atom = atomList.get(i)
          String resNameA = atom.getResidueName()
          //logger.info("Residue A: " + resNameA)
          double coorAX = atom.getX()
          double coorAY = atom.getY()
          double coorAZ = atom.getZ()
          Residue residueB = residueListB.get(resNum - firstRes)
          List<Atom> atomListB = residueB.getAtomList()
          Atom atomB = atomListB.get(i)
          String resNameB = atomB.getResidueName()
          //logger.info("Residue B: " + resNameB)
          double coorBX = atomB.getX()
          double coorBY = atomB.getY()
          double coorBZ = atomB.getZ()
          if (coorAX == coorBX && coorAY == coorBY && coorAZ == coorBZ && resNameA == resNameB) {
            atom.setAltLoc(' ' as Character)
            atomB.setAltLoc(' ' as Character)
            atom.setOccupancy(1.0)
          }
        }
      }
      if (titrationPH > 0) {
        diffractionData.writeModel(removeExtension(filenames[0]) + ".pdb", excludeAtoms, titrationPH)
      } else {
        diffractionData.writeModel(removeExtension(filenames[0]) + ".pdb")
      }
      diffractionData.writeData(removeExtension(filenames[0]) + ".mtz")
    } else if (Comm.world().rank() == 0) {
      String ext = FilenameUtils.getExtension(modelFilename)
      modelFilename = FilenameUtils.removeExtension(modelFilename)
      if (ext.toUpperCase().contains("XYZ")) {
        algorithmFunctions.saveAsXYZ(activeAssembly, new File(modelFilename + ".xyz"))
      } else {
        properties.setProperty("standardizeAtomNames", "false")
        File modelFile = saveDirFile(activeAssembly.getFile())
        PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly,
            activeAssembly.getForceField(), properties)
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


