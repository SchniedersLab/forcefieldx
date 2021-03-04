//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.groovy.test

import ffx.crystal.Crystal
import ffx.crystal.SymOp
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level
import java.util.stream.Collectors

import static ffx.potential.utils.CrystalSuperposeFunctions.generateBaseSphere
import static ffx.potential.utils.CrystalSuperposeFunctions.mimicSphere
import static ffx.potential.utils.Superpose.*
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.PI
import static org.apache.commons.math3.util.FastMath.cbrt

/**
 * The PACCOM script calculates intracrystollographic distances.
 *
 * @author Okimasa OKADA
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler, Kaleb Bierstedt, and Micheal Schnieders 2020
 * <br>
 * Usage:
 * <br>
 * ffxc test.CrystalSuperpose &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Compare crystal packings based on intermolecular distances.", name = "ffxc test.CrystalSuperpose")
class CrystalSuperpose extends PotentialScript {

  /**
   * --nm or --numberMolecules Number of molecules to include from each crystal in RMSD comparison.
   */
  @Option(names = ['--nm', '--numberMolecules'], paramLabel = '50',
      description = 'Determines crystal sphere size for comparison.')
  int nMolecules = 50

  /**
   * --im or --inflatedMolecules Number of molecules in the inflated sphere.
   */
  @Option(names = ['--im', '--inflatedMolecules'], paramLabel = '500',
      description = 'Determines number of molecules in inflated sphere.')
  int inflatedMolecules = 500

  // TODO: implement allVsAll
  /**
   * --all or --allVsAll Compute RMSD of between each file.
   */
  @Option(names = ['--all', '--allVsAll'], paramLabel = "false", defaultValue = "false",
      description = 'Calculate the RMSD between each sphere.')
  private boolean allVsAll = false

  /**
   * --op or --originalPACCOM Use Algorithm supplied by original PACCOM.
   */
  @Option(names = ['--op', '--originalPACCOM'], paramLabel = "false", defaultValue = "false",
      description = 'Utilize a version of PACCOM that mirrors original.')
  private static boolean original = false

  /**
   * -s or --save Save out individual XYZ/PDB files for each created sphere.
   */
  @Option(names = ['-s', '--save'], paramLabel = "false", defaultValue = "false",
      description = 'Save each sphere as an aperiodic system to a PDB file.')
  private boolean saveFiles = false

  // TODO: implement closestDistance
  /**
   * --cd or --closestDistance Neighbor molecules calculated via closest atom (not currently implemented).
   */
//    @Option(names = ['--cd', '--closestDistance'], paramLabel = "false", defaultValue = "false",
//            description = 'Neighbors calculated via closest atom.')
//    private static boolean closestDistance = false

  /**
   * --nh or --noHydrogens Perform comparison without hydrogens.
   */
  @Option(names = ['--nh', '--noHydrogens'], paramLabel = "false", defaultValue = "false",
      description = 'Crystal RMSD calculated without hydrogen atoms.')
  private static boolean noHydrogens = false

  /**
   * Select individual atoms to match and calculate RMSD rather than using all atoms.
   */
  @Option(names = ['--ca', '--centerAtoms'], arity = "1..*", paramLabel = "atom integers",
      description = 'Specify atoms to match and calculate RMSD. Otherwise, all atom.')
  private int[] centerAtoms = null

  /**
   * The final argument(s) should be two or more filenames.
   */
  @Parameters(arity = "2..*", paramLabel = "files",
      description = 'Atomic coordinate files to compare in XYZ format.')
  List<String> filenames = null

  private MolecularAssembly baseAssembly

  private MolecularAssembly[] assemblies

  public double[][] coordinates = null

  /**
   * CrystalSuperpose Constructor.
   */
  CrystalSuperpose() {
    this(new Binding())
  }

  /**
   * CrystalSuperpose Constructor.
   * @param binding Groovy Binding to use.
   */
  CrystalSuperpose(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  CrystalSuperpose run() {
    System.setProperty("vdwterm", "false")
    // Read in structures
    if (!init()) {
      return
    }

    // Ensure file exists/create active molecular assembly
    if (filenames != null && filenames.size() > 1) {
      baseAssembly = potentialFunctions.open(filenames.get(0))
      assemblies = potentialFunctions.openAll(filenames.get(1))
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return
    } else {
      assemblies = [activeAssembly]
    }

    List<Atom> exampleAtoms = baseAssembly.getAtomList()
    int nAtoms = exampleAtoms.size()
    int numHydrogens = 0

    // If comparison atoms are not specified, use all atoms.
    int[] comparisonAtoms
    if (centerAtoms != null) {
      comparisonAtoms = new int[centerAtoms.length]
      for (int i = 0; i < centerAtoms.size(); i++) {
        comparisonAtoms[i] = centerAtoms[i] - 1
        if (comparisonAtoms[i] < 0 || comparisonAtoms[i] >= nAtoms) {
          logger.severe(" Selected atoms are outside of molecular size.")
        }
      }
    } else {
      if (noHydrogens) {
        int n = 0
        for (int i = 0; i < nAtoms; i++) {
          if (!exampleAtoms.get(i).isHydrogen()) {
            n++
          }
        }
        comparisonAtoms = new int[n]
        n = 0
        for (int i = 0; i < nAtoms; i++) {
          if (!exampleAtoms.get(i).isHydrogen()) {
            comparisonAtoms[n++] = i
          }
        }
      } else {
        comparisonAtoms = new int[nAtoms]
        for (int i = 0; i < nAtoms; i++) {
          comparisonAtoms[i] = i
        }
      }
    }

    int compareAtomsSize = comparisonAtoms.size()
    logger.info(format(" Number of atoms being compared: %d of %d",
        compareAtomsSize, nAtoms))

    if (logger.isLoggable(Level.FINE)) {
      for (int integerVal : comparisonAtoms) {
        logger.fine(format(" %d", integerVal + 1))
      }
    }

    // Here we will use the unit cell, to create a new replicates crystal that may be
    // a different size (i.e. larger).
    Crystal unitCell = baseAssembly.getCrystal().getUnitCell()
    double asymmetricUnitVolume = unitCell.volume / unitCell.getNumSymOps()
    // Add wiggle room for boundary cutoffs
    double inflationFactor = 1.0
    // Estimate a radius that will include "nExpanded".
    double radius = cbrt((3.0 / (4.0 * PI) * inflatedMolecules * asymmetricUnitVolume)) + inflationFactor

    logger.info(format(" Copies in target sphere:     %16d", inflatedMolecules))
    logger.info(format(" Estimated spherical radius:  %16.2f", radius))
    logger.info(format(" Number of copies to compare: %16d", nMolecules))

    // Generate the base sphere
    logger.info(" Expanding system 1.")
    // Save the SymOps used to expand System 1.
    List<SymOp> symOps = new ArrayList<>()
    MolecularAssembly baseSphere =
        generateBaseSphere(baseAssembly, nMolecules, original, radius, symOps)

    // This is ReplicatesCrystal that was used to generate a sphere of molecules.
    Crystal baseCrystal = baseAssembly.getCrystal()

    // Loop over systems that should be superposed on the first system.
    for (MolecularAssembly molecularAssembly : assemblies) {
      // Collect asymmetric unit atomic coordinates.
      Atom[] atoms = molecularAssembly.getAtomArray()
      double[] x = new double[nAtoms]
      double[] y = new double[nAtoms]
      double[] z = new double[nAtoms]
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atoms[i]
        x[i] = atom.getX()
        y[i] = atom.getY()
        z[i] = atom.getZ()
      }

      // Allocate space for coordinates after application of a SymOp.
      double[] xS = new double[nAtoms]
      double[] yS = new double[nAtoms]
      double[] zS = new double[nAtoms]

      // Init the RMSD to its max value.
      double rmsdValue = Double.MAX_VALUE

      // Loop over the unit cell SymOps
      Crystal unitCellCrystal = molecularAssembly.getCrystal().getUnitCell()
      int numSymOps = unitCellCrystal.getNumSymOps()

      logger.info(format("\n Looping over SymOps"))
      for (int iSym = 0; iSym < numSymOps; iSym++) {

        // Apply SymOp to the asymmetric unit atoms.
        SymOp symOp = unitCellCrystal.spaceGroup.getSymOp(iSym)
        logger.finer(format(" SymOp: %d \n%s", iSym, symOp.toString()))
        unitCellCrystal.applySymOp(nAtoms, x, y, z, xS, yS, zS, symOp)

        // Update atomic positions
        for (int i = 0; i < nAtoms; i++) {
          atoms[i].moveTo(xS[i], yS[i], zS[i])
        }

        // Compute the RMSD
        double tempRMSD = crystalComparison(baseSphere, baseCrystal, symOps, molecularAssembly, radius)
        logger.info(format(" SymOp %2d RMSD: %16.8f A", iSym, tempRMSD))
        if (tempRMSD < rmsdValue) {
          rmsdValue = tempRMSD
        }
      }

      logger.info(format("\n Best RMSD: %16.8f A", rmsdValue))
    }

    return this
  }

  double crystalComparison(MolecularAssembly baseSphere, Crystal baseCrystal,
      List<SymOp> symOps, MolecularAssembly assembly2, double radius) {

    // Collect the atomic positions for both systems.
    Atom[] baseAtoms = baseSphere.getAtomList()
    int nExpanded = baseAtoms.size()
    double[] x1 = new double[nExpanded * 3]
    for (int i = 0; i < nExpanded; i++) {
      // TODO: Pull out a subset of atoms (i.e. no hydrogen; only selected heavy atoms)
      Atom a = baseAtoms[i]
      x1[i * 3] = a.getX()
      x1[i * 3 + 1] = a.getY()
      x1[i * 3 + 2] = a.getZ()
    }

    // Apply the baseSphere symOps to the second assembly.
    logger.fine(" Generating Mimic Sphere:")
    double[] x2 = mimicSphere(baseCrystal, assembly2, symOps, radius)

    double[] mass = new double[nExpanded]
    Arrays.fill(mass, 1.0)

    logger.fine(format("\n Superposing %d atoms.", nExpanded))
    // RMSD before alignment.
    double originalRMSD = rmsd(x1, x2, mass)

    // Calculate the translation on only the used subset, but apply it to the entire structure.
    applyTranslation(x1, calculateTranslation(x1, mass))
    applyTranslation(x2, calculateTranslation(x2, mass))
    double translatedRMSD = rmsd(x1, x2, mass)

    // Calculate the rotation on only the used subset, but apply it to the entire structure.
    applyRotation(x2, calculateRotation(x1, x2, mass))
    double rotatedRMSD = rmsd(x1, x2, mass)

    logger.fine(
        format(" RMSD Original %16.8f, Translate %16.8f, Rotate %16.8f",
            originalRMSD, translatedRMSD, rotatedRMSD))

    // Compute the RMSD.
    return rotatedRMSD
  }

  @Override
  List<Potential> getPotentials() {
    if (assemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(assemblies).
          filter {a -> a != null
          }.
          map {a -> a.getPotentialEnergy()
          }.
          filter {e -> e != null
          }.
          collect(Collectors.toList())
    }
  }
}