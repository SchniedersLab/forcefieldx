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
package ffx.potential.groovy

import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The Cart2Frac script converts Cartesian coordinates to Fractional.
 * <br>
 * Usage:
 * <br>
 * ffxc Cart2Frac &lt;filename&gt;
 */
@Command(description = " Convert Cartesian coordinates to fractional.", name = "ffxc Cart2Frac")
class Cart2Frac extends PotentialScript {

  /**
   * The final argument should be a file in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  /**
   * Save a reference to the MolecularAssembly instances to destroy their potentials.
   */
  private MolecularAssembly[] molecularAssemblies

  /**
   * Return Cartesian Coordinate input.
   * @return
   */
  double[][][] getCart() {
    return cartCoordinates
  }

  /**
   * Return Fractional Coordinate output.
   * @return
   */
  double[][][] getFrac() {
    return fracCoordinates
  }

  /**
   * Cartesian coordinate input.
   */
  private double[][][] cartCoordinates = null

  /**
   * Fractional coordinate output.
   */
  private double[][][] fracCoordinates = null

  /**
   * Cart2Frac constructor.
   */
  Cart2Frac() {
    this(new Binding())
  }

  /**
   * Cart2Frac constructor.
   * @param binding Binding to use.
   */
  Cart2Frac(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Cart2Frac run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load one or more MolecularAssembly instances.
    molecularAssemblies = getActiveAssemblies(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    int num = molecularAssemblies.length
    cartCoordinates = new double[num][][]
    fracCoordinates = new double[num][][]

    // Loop over each system.
    for (int i = 0; i < num; i++) {
      def molecularAssembly = molecularAssemblies[i]
      logger.info("\n Converting from Cartesian to fractional coordinates for " +
          molecularAssembly.toString())
      Crystal crystal = molecularAssembly.getCrystal().getUnitCell()

      Atom[] atoms = molecularAssembly.getAtomArray()
      int nAtoms = atoms.length;
      fracCoordinates[i] = new double[nAtoms][3]
      cartCoordinates[i] = new double[nAtoms][3]

      double[] frac = new double[3]
      double[] cart = new double[3]

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index]
        atom.getXYZ(cart)
        crystal.toFractionalCoordinates(cart, frac)
        atom.moveTo(frac)

        cartCoordinates[i][index][0] = cart[0]
        cartCoordinates[i][index][1] = cart[1]
        cartCoordinates[i][index][2] = cart[2]

        fracCoordinates[i][index][0] = frac[0]
        fracCoordinates[i][index][1] = frac[1]
        fracCoordinates[i][index][2] = frac[2]
      }
    }

    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    String name = FilenameUtils.getName(filename)
    String ext = FilenameUtils.getExtension(name)
    name = FilenameUtils.removeExtension(name)

    String dirName = FilenameUtils.getFullPath(saveDir.getAbsolutePath())

    if (ext.toUpperCase().contains("XYZ")) {
      potentialFunctions.saveAsXYZ(molecularAssemblies[0], new File(dirName + name + ".xyz"))
    } else {
      potentialFunctions.saveAsPDB(molecularAssemblies, new File(dirName + name + ".pdb"))
    }

    binding.setVariable("cart", cartCoordinates)
    binding.setVariable("frac", fracCoordinates)

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (molecularAssemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(molecularAssemblies).filter {a -> a != null
      }.map {a -> a.getPotentialEnergy()
      }.filter {e -> e != null
      }.collect(Collectors.toList())
    }
  }
}
