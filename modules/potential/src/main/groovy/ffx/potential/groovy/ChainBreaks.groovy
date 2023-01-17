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

import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Residue.ResidueType
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.PotentialsUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static ffx.numerics.math.DoubleMath.dist
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getExtension
import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension
import static org.apache.commons.math3.util.FastMath.abs

@Command(description = " Fix chain breaks in a pdb file.", name = "ChainBreaks")
class ChainBreaks extends PotentialScript {

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  private String filename = null

  List<String> chainBreaks = new ArrayList<>()
  List<double[]> newCoordinates
  PotentialsUtils potentialsUtils = new PotentialsUtils()

  /**
   * ChainBreaks constructor.
   */
  ChainBreaks() {
    this(new Binding())
  }

  /**
   * ChainBreaks constructor.
   * @param binding The Groovy Binding to use.
   */
  ChainBreaks(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  ChainBreaks run() {
    // Init the context and bind variables.
    if (!init()) {
      return null
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return null
    }

    // Get the base name of the file and its extension.
    String name = getName(filename)
    String ext = getExtension(name)
    name = removeExtension(name)

    if (!ext.containsIgnoreCase("pdb")) {
      logger.info(format(" The file extension does not include 'pdb': %s", filename))
      return null
    }

    List<Residue> residues = activeAssembly.getResidueList()
    chainBreaks = findChainBreaks(residues, 3)
    logger.info(format(" Fixing Chain Breaks in %s", filename))

    // Use the current base directory, or update if necessary based on the given filename.
    String dirString = getBaseDirString(filename)
    File editedPDBFile = new File(dirString + name + "_edited.pdb")

    logger.info(format(" Saving New Coordinates to:\n %s", editedPDBFile.toString()))
    potentialsUtils.saveAsPDB(activeAssembly, editedPDBFile, false, false)

    return this
  }

  private List<String> findChainBreaks(List<Residue> residues, double cutoff) {
    List<List<Residue>> subChains = new ArrayList<>()
    List<String> chainBreaks = new ArrayList<>()

    // Chain-start atom: N (amino) / O5* (nucleic)
    // Chain-end atom:   C (amino) / O3* (nucleic)
    ResidueType rType = residues.get(0).getResidueType()
    String startAtName
    String endAtName
    switch (rType) {
      case ResidueType.AA:
        startAtName = "N"
        endAtName = "C"
        break
      case ResidueType.NA:
        boolean namedStar =
            residues.stream()
                .flatMap((Residue r) -> r.getAtomList().stream())
                .anyMatch((Atom a) -> a.getName() == "O5*")
        if (namedStar) {
          startAtName = "O5*"
          endAtName = "O3*"
        } else {
          startAtName = "O5\'"
          endAtName = "O3\'"
        }
        break
      case ResidueType.UNK:
      default:
        logger.fine(
            " Not attempting to find chain breaks for chain with residue "
                + residues.get(0).toString())
        return null
    }

    List<Residue> subChain = null
    Residue previousResidue = null
    Atom priorEndAtom = null

    for (Residue residue : residues) {
      List<Atom> resAtoms = residue.getAtomList()
      if (priorEndAtom == null) {
        // Initialization.
        subChain = new ArrayList<>()
        subChain.add(residue)
        subChains.add(subChain)
      } else {
        // Find the start atom of the current residue.
        Atom startAtom = null
        for (Atom a : resAtoms) {
          if (a.getName().equalsIgnoreCase(startAtName)) {
            startAtom = a
            break
          }
        }
        if (startAtom == null) {
          subChain.add(residue)
          continue
        }
        // Compute the distance between the previous carbonyl carbon and the current nitrogen.
        double r = dist(priorEndAtom.getXYZ(null), startAtom.getXYZ(null))

        if (r > cutoff) {
          // Start a new chain.
          subChain = new ArrayList<>()
          subChain.add(residue)
          subChains.add(subChain)
          chainBreaks.add("C " + previousResidue.toString() + " N " + residue.toString())
          fixChainBreaks(priorEndAtom.getXYZ(null), startAtom.getXYZ(null))
          priorEndAtom.setXYZ(newCoordinates.get(0))
          startAtom.setXYZ(newCoordinates.get(1))
        } else {
          // Continue the current chain.
          subChain.add(residue)
        }
      }

      // Save the carbonyl carbon.
      for (Atom a : resAtoms) {
        if (a.getName().equalsIgnoreCase(endAtName)) {
          priorEndAtom = a
          break
        }
      }
      previousResidue = residue
    }

    return chainBreaks
  }

  private void fixChainBreaks(double[] cCoordinates, double[] nCoordinates) {
    logger.info(" Generating new coordinates.")
    newCoordinates = new ArrayList<>()
    double distance = dist(cCoordinates, nCoordinates)
    while (distance > 3) {
      double dx = abs((cCoordinates[0] - nCoordinates[0]) / 4.0)
      double dy = abs((cCoordinates[1] - nCoordinates[1]) / 4.0)
      double dz = abs((cCoordinates[2] - nCoordinates[2]) / 4.0)

      if (cCoordinates[0] > nCoordinates[0]) {
        cCoordinates[0] = cCoordinates[0] - dx
        nCoordinates[0] = nCoordinates[0] + dx
      } else if (cCoordinates[0] < nCoordinates[0]) {
        cCoordinates[0] = cCoordinates[0] + dx
        nCoordinates[0] = nCoordinates[0] - dx
      }

      if (cCoordinates[1] > nCoordinates[1]) {
        cCoordinates[1] = cCoordinates[1] - dy
        nCoordinates[1] = nCoordinates[1] + dy
      } else if (cCoordinates[1] < nCoordinates[1]) {
        cCoordinates[1] = cCoordinates[1] + dy
        nCoordinates[1] = nCoordinates[1] - dy
      }

      if (cCoordinates[2] > nCoordinates[2]) {
        cCoordinates[2] = cCoordinates[2] - dz
        nCoordinates[2] = nCoordinates[2] + dz
      } else if (cCoordinates[2] < nCoordinates[0]) {
        cCoordinates[2] = cCoordinates[2] + dz
        nCoordinates[2] = nCoordinates[2] - dz
      }

      distance = dist(cCoordinates, nCoordinates)
    }
    newCoordinates.add(cCoordinates)
    newCoordinates.add(nCoordinates)
  }

}


