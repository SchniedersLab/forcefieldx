//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
 * The Frac2Cart script converts from Fractional to Cartesian coordinates.
 * <br>
 * Usage:
 * <br>
 * ffxc Frac2Cart &lt;filename&gt;
 */
@Command(description = " Convert fractional to Cartesian coordinates.", name = "ffxc Frac2Cart")
class Frac2Cart extends PotentialScript {

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  private MolecularAssembly[] molecularAssemblies

  public double[][] cartCoordinates = null
  public double[][] fracCoordinates = null

  /**
   * Frac2Cart Constructor.
   */
  Frac2Cart() {
    this(new Binding())
  }

  /**
   * Frac2Cart Constructor.
   * @param binding Groovy Binding to use.
   */
  Frac2Cart(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Frac2Cart run() {

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

    logger.info("\n Converting from fractional to Cartesian coordinates for " + filename)

    // Loop over each system.
    for (int i = 0; i < molecularAssemblies.length; i++) {
      def molecularAssembly = molecularAssemblies[i]
      Crystal crystal = molecularAssembly.getCrystal().getUnitCell()

      Atom[] atoms = molecularAssembly.getAtomArray()
      int nAtoms = atoms.length;
      fracCoordinates = new double[nAtoms][3]
      cartCoordinates = new double[nAtoms][3]

      double[] frac = new double[3]
      double[] cart = new double[3]

      for (int index = 0; index < nAtoms; index++) {
        Atom atom = atoms[index]
        atom.getXYZ(frac)
        crystal.toCartesianCoordinates(frac, cart)
        atom.moveTo(cart)

        cartCoordinates[index][0] = cart[0]
        cartCoordinates[index][1] = cart[1]
        cartCoordinates[index][2] = cart[2]

        fracCoordinates[index][0] = frac[0]
        fracCoordinates[index][1] = frac[1]
        fracCoordinates[index][2] = frac[2]
      }
    }

    // Configure the base directory if it has not been set.
    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = saveDir.toString() + File.separator
    String name = FilenameUtils.getName(filename)
    String ext = FilenameUtils.getExtension(name)
    name = FilenameUtils.removeExtension(name)

    if (ext.toUpperCase().contains("XYZ")) {
      potentialFunctions.saveAsXYZ(molecularAssemblies[0], new File(dirName + name + ".xyz"))
    } else {
      potentialFunctions.saveAsPDB(molecularAssemblies, new File(dirName + name + ".pdb"))
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (molecularAssemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(molecularAssemblies).filter {
        a -> a != null
      }.map {a -> a.getPotentialEnergy()
      }.filter {e -> e != null
      }.collect(Collectors.toList())
    }
  }
}
