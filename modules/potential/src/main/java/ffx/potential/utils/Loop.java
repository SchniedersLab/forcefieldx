// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.potential.utils;

import static ffx.potential.bonded.BondedUtils.determineIntxyz;
import static java.lang.String.copyValueOf;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.*;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Loop class.
 *
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class Loop {

  private static final Logger logger = Logger.getLogger(Loop.class.getName());
  private final LoopClosure loopClosure;
  private final MolecularAssembly molecularAssembly;
  private final Random random = new Random();
  private static final int MAX_SOLUTIONS = 16;
  private final double[][] rN = new double[3][3];
  private final double[][] rA = new double[3][3];
  private final double[][] rC = new double[3][3];
  private double[][] altCoords;
  private boolean useAltCoords = false;

  /**
   * Constructor for Loop.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param firstResidue a int.
   * @param endResidue a int.
   */
  public Loop(MolecularAssembly molecularAssembly, int firstResidue, int endResidue) {
    loopClosure = new LoopClosure();
    this.molecularAssembly = molecularAssembly;
    generateLoops(firstResidue, endResidue);
  }

  /**
   * Constructor for Loop.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public Loop(MolecularAssembly molecularAssembly) {
    loopClosure = new LoopClosure();
    this.molecularAssembly = molecularAssembly;
    this.altCoords = new double[molecularAssembly.getAtomArray().length][3];
  }

  /**
   * generateLoops.
   *
   * @param firstResidue a int.
   * @param endResidue a int.
   * @param coordinates an array of {@link double} objects.
   * @return a {@link java.util.List} object.
   */
  public List<double[]> generateLoops(int firstResidue, int endResidue, double[] coordinates) {
    setAltCoordinates(coordinates);
    return generateLoops(firstResidue, endResidue);
  }

  /**
   * generateLoops.
   *
   * @param firstResidue a int.
   * @param endResidue a int.
   * @return a {@link java.util.List} object.
   */
  public List<double[]> generateLoops(int firstResidue, int endResidue) {
    List<Atom> backBoneAtoms = molecularAssembly.getBackBoneAtoms();

    List<double[]> solutions = new ArrayList<>();
    logger.info(format(" First residue:   %d\n", firstResidue));
    logger.info(format(" Ending residue:  %d\n", endResidue));
    for (Atom atom : backBoneAtoms) {
      int resID = atom.getResidueNumber();
      if (resID > endResidue) {
        // terminates the collection of atom coordinates
        break;
      }

      if (resID >= firstResidue) {
        int ir = resID - firstResidue;
        double[] initArray;
        if (useAltCoords) {
          initArray = altCoords[atom.getIndex() - 1];
        } else {
          initArray = atom.getXYZ(null);
        }

        switch (atom.getAtomType().name) {
          case "C" -> {
            // Backbone carbon coordinates are stored in rC[]
            rC[ir][0] = initArray[0];
            rC[ir][1] = initArray[1];
            rC[ir][2] = initArray[2];
          }
          case "CA" -> {
            // Backbone alpha carbon coordinates are stored in rA[]
            rA[ir][0] = initArray[0];
            rA[ir][1] = initArray[1];
            rA[ir][2] = initArray[2];
          }
          case "N" -> {
            // Backbone nitrogen coordinates are stored in rN[]
            rN[ir][0] = initArray[0];
            rN[ir][1] = initArray[1];
            rN[ir][2] = initArray[2];
          }
        }
      }
    }

    // Method that solves the 16th degree, 3 peptide polynomial.
    var rSolnN = new double[MAX_SOLUTIONS][3][3];
    var rSolnA = new double[MAX_SOLUTIONS][3][3];
    var rSolnC = new double[MAX_SOLUTIONS][3][3];

    var numSolutions = loopClosure.solve3PepPoly(rN[0], rA[0], rA[2], rC[2], rSolnN, rSolnA, rSolnC);

    StringBuilder sb = new StringBuilder();
    sb.append(format(" First residue:                %d\n", firstResidue));
    sb.append(format(" Ending residue:               %d\n", endResidue));
    sb.append(format(" Number of solutions:          %d\n", numSolutions));
    logger.info(sb.toString());

    for (int k = 0; k < numSolutions; k++) {
      solutions.add(getSolutionCoordinates(k, rSolnN, rSolnA, rSolnC, firstResidue, endResidue));
    }

    return solutions;
  }

  /**
   * getRA.
   *
   * @return an array of {@link double} objects.
   */
  public double[][] getRA() {
    return rA;
  }

  /**
   * getRC.
   *
   * @return an array of {@link double} objects.
   */
  public double[][] getRC() {
    return rC;
  }

  /**
   * getRN.
   *
   * @return an array of {@link double} objects.
   */
  public double[][] getRN() {
    return rN;
  }

  private double[] getSolutionCoordinates(int k, double[][][] rSolnN, double[][][] rSolnA,
      double[][][] rSolnC, int startResidue, int endResidue) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rN[i][j] = rSolnN[k][i][j];
        rA[i][j] = rSolnA[k][i][j];
        rC[i][j] = rSolnC[k][i][j];
      }
    }

    Polymer[] newChain = molecularAssembly.getChains();

    Atom[] atomArray = molecularAssembly.getAtomArray();
    double[][] coordsArray = new double[atomArray.length][3];
    if (useAltCoords) {
      arraycopy(this.altCoords, 0, coordsArray, 0, coordsArray.length);
    } else {
      for (int i = 0; i < atomArray.length; i++) {
        Atom a = atomArray[i];
        coordsArray[i][0] = a.getX();
        coordsArray[i][1] = a.getY();
        coordsArray[i][2] = a.getZ();
      }
    }

    // Update coordinates of solved loop.
    for (int i = startResidue; i <= endResidue; i++) {
      Residue newResidue = newChain[0].getResidue(i);

      // Determine the new position of backbone atoms and get an offset to apply
      // for other atoms.
      var backBoneAtoms = newResidue.getBackboneAtoms();
      var cOffset = new double[3];
      var nOffset = new double[3];
      var aOffset = new double[3];
      for (Atom backBoneAtom : backBoneAtoms) {
        int index = backBoneAtom.getIndex() - 1;
        switch (backBoneAtom.getAtomType().name) {
          case "C" -> {
            cOffset[0] = rC[i - startResidue][0] - coordsArray[index][0];
            cOffset[1] = rC[i - startResidue][1] - coordsArray[index][1];
            cOffset[2] = rC[i - startResidue][2] - coordsArray[index][2];
            arraycopy(rC[i - startResidue], 0, coordsArray[index], 0, 3);
          }
          case "N" -> {
            nOffset[0] = rN[i - startResidue][0] - coordsArray[index][0];
            nOffset[1] = rN[i - startResidue][1] - coordsArray[index][1];
            nOffset[2] = rN[i - startResidue][2] - coordsArray[index][2];
            arraycopy(rN[i - startResidue], 0, coordsArray[index], 0, 3);
          }
          case "CA" -> {
            aOffset[0] = rA[i - startResidue][0] - coordsArray[index][0];
            aOffset[1] = rA[i - startResidue][1] - coordsArray[index][1];
            aOffset[2] = rA[i - startResidue][2] - coordsArray[index][2];
            arraycopy(rA[i - startResidue], 0, coordsArray[index], 0, 3);
          }
          default -> {
          }
        }
      }

      // Loop through all other backbone atoms and apply the above offsets.
      for (Atom backBoneAtom : backBoneAtoms) {
        int index = backBoneAtom.getIndex() - 1;
        switch (backBoneAtom.getAtomType().name) {
          case "C", "N", "CA" -> {
            // Do nothing. The previous loop already set these values.
          }
          case "O" -> {
            coordsArray[index][0] += cOffset[0];
            coordsArray[index][1] += cOffset[1];
            coordsArray[index][2] += cOffset[2];
          }
          case "H" -> {
            coordsArray[index][0] += nOffset[0];
            coordsArray[index][1] += nOffset[1];
            coordsArray[index][2] += nOffset[2];
          }
          default -> {
            coordsArray[index][0] += aOffset[0];
            coordsArray[index][1] += aOffset[1];
            coordsArray[index][2] += aOffset[2];
          }
        }
      }

      // Update sidechain atoms based on the alpha carbon.
      for (Atom sideChainAtom : newResidue.getSideChainAtoms()) {
        int index = sideChainAtom.getIndex() - 1;
        coordsArray[index][0] += aOffset[0];
        coordsArray[index][1] += aOffset[1];
        coordsArray[index][2] += aOffset[2];
      }
    }

    double[] coordsArray1D = new double[atomArray.length * 3];
    for (int i = 0; i < atomArray.length; i++) {
      arraycopy(coordsArray[i], 0, coordsArray1D, i * 3, 3);
    }
    return coordsArray1D;
  }

  private void setAltCoordinates(double[] coordinates) {
    for (int i = 0; i < coordinates.length / 3; i++) {
      arraycopy(coordinates, i * 3, this.altCoords[i], 0, 3);
    }

    useAltCoords = true;
  }

  private double[][] fillCoordsArray(Atom atom, double[][] coordsArray, double[] determinedXYZ) {
    int XYZIndex = atom.getIndex() - 1;
    coordsArray[XYZIndex][0] = determinedXYZ[0];
    coordsArray[XYZIndex][1] = determinedXYZ[1];
    coordsArray[XYZIndex][2] = determinedXYZ[2];
    return coordsArray;
  }
}
