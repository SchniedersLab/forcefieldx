// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.utils;

import static ffx.potential.bonded.BondedUtils.determineIntxyz;
import static java.lang.System.arraycopy;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;

/**
 * Loop class.
 *
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class Loop {

  private static final Logger logger = Logger.getLogger(Loop.class.getName());
  public final LoopClosure loopClosure;
  private final MolecularAssembly molecularAssembly;
  private final Random random = new Random();
  int maxSolution = 16;
  double[][] rN = new double[3][3];
  double[][] rA = new double[3][3];
  double[][] rC = new double[3][3];
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

    boolean bool1 = true;
    int i = 0;
    List<double[]> solutions = new ArrayList<>();
    logger.info(String.format(" First residue:   %d\n", firstResidue));
    logger.info(String.format(" Ending residue:  %d\n", endResidue));
    while (bool1) {
      Atom atom = backBoneAtoms.get(i);
      int resID = atom.getResidueNumber();
      if (resID > endResidue) {
        // terminates the collection of atom coordinates
        bool1 = false;
      }

      if (resID < firstResidue) {
        i++;
      } else if (resID >= firstResidue && resID <= endResidue) {
        int ir = resID - firstResidue;
        String atmname = atom.getAtomType().name;
        String ca = "CA";
        String c = "C";
        String n = "N";
        // coordinateArray temporarily holds the coordinates of a
        //      specific atom from the pdb file
        double[] initArray;
        if (useAltCoords) {
          initArray = altCoords[atom.getIndex() - 1];
        } else {
          initArray = atom.getXYZ(null);
        }
        if (atmname.contentEquals(n)) {
          // Backbone nitrogen coordinates are stored in rN[]
          rN[ir][0] = initArray[0];
          rN[ir][1] = initArray[1];
          rN[ir][2] = initArray[2];
        } else if (atmname.contentEquals(ca)) {
          // Backbone alpha carbon coordinates are stored in rA[]
          rA[ir][0] = initArray[0];
          rA[ir][1] = initArray[1];
          rA[ir][2] = initArray[2];
        } else if (atmname.contentEquals(c)) {
          // Backbone carbon coordinates are stored in rC[]
          rC[ir][0] = initArray[0];
          rC[ir][1] = initArray[1];
          rC[ir][2] = initArray[2];
        }
        i++;
      }
    }

    /** Method that solves the 16th degree, 3 peptide polynomial. */
    double[][][] rSolnN = new double[maxSolution][3][3];
    double[][][] rSolnA = new double[maxSolution][3][3];
    double[][][] rSolnC = new double[maxSolution][3][3];
    int[] nSoln = new int[1];

    loopClosure.solve3PepPoly(rN[0], rA[0], rA[2], rC[2], rSolnN, rSolnA, rSolnC, nSoln);

    StringBuilder sb = new StringBuilder();
    sb.append(String.format(" First residue:                %d\n", firstResidue));
    sb.append(String.format(" Ending residue:               %d\n", endResidue));
    sb.append(String.format(" Number of solutions:          %d\n", nSoln[0]));
    logger.info(sb.toString());

    for (int k = 0; k < nSoln[0]; k++) {
      double[] coordsArray;
      coordsArray = getSolutionCoordinates(k, rSolnN, rSolnA, rSolnC, firstResidue, endResidue);
      solutions.add(coordsArray);
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

  /**
   * useAltCoordinates.
   *
   * @param bool a boolean.
   */
  public void useAltCoordinates(boolean bool) {
    this.useAltCoords = bool;
  }

  private double[] getSolutionCoordinates(
      int k,
      double[][][] rSolnN,
      double[][][] rSolnA,
      double[][][] rSolnC,
      int startResidue,
      int endResidue) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rN[i][j] = rSolnN[k][i][j];
        rA[i][j] = rSolnA[k][i][j];
        rC[i][j] = rSolnC[k][i][j];
      }
    }

    Polymer[] newChain = molecularAssembly.getChains();
    List<Atom> backBoneAtoms;

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

    // Loop through residues to build backbone C,N,CA
    for (int i = startResidue; i <= endResidue; i++) {
      Residue newResidue = newChain[0].getResidue(i);
      Residue backResidue = newChain[0].getResidue(i - 1);
      backBoneAtoms = newResidue.getBackboneAtoms();
      Atom C = null;
      Atom N = null;
      Atom CA = null;
      double[] c = new double[3];
      double[] ca = new double[3];
      double[] n = new double[3];
      double[] bc = new double[3];
      double[] determinedXYZ = new double[3];

      for (Atom backBoneAtom : backBoneAtoms) {
        //    backBoneAtom.setBuilt(true);
        int backBoneIndex = backBoneAtom.getIndex() - 1;

        switch (backBoneAtom.getAtomType().name) {
          case "C":
            C = backBoneAtom;
            coordsArray[backBoneIndex][0] = rC[i - startResidue][0];
            coordsArray[backBoneIndex][1] = rC[i - startResidue][1];
            coordsArray[backBoneIndex][2] = rC[i - startResidue][2];
            break;
          case "N":
            N = backBoneAtom;
            coordsArray[backBoneIndex][0] = rN[i - startResidue][0];
            coordsArray[backBoneIndex][1] = rN[i - startResidue][1];
            coordsArray[backBoneIndex][2] = rN[i - startResidue][2];
            break;
          case "CA":
            CA = backBoneAtom;
            coordsArray[backBoneIndex][0] = rA[i - startResidue][0];
            coordsArray[backBoneIndex][1] = rA[i - startResidue][1];
            coordsArray[backBoneIndex][2] = rA[i - startResidue][2];
            break;
          default:
            break;
        }
      }
    }

    // Loop through again to build Hydrogens and side-chains with current C,N,CA coordinates
    for (int i = startResidue + 1; i < endResidue - 1; i++) {
      Residue newResidue = newChain[0].getResidue(i);
      Residue backResidue = newChain[0].getResidue(i - 1);
      Residue forwardResidue = newChain[0].getResidue(i + 1);
      backBoneAtoms = newResidue.getBackboneAtoms();

      double[] c = new double[3];
      double[] ca = new double[3];
      double[] n = new double[3];
      double[] bc = new double[3];
      double[] fn = new double[3];
      double[] determinedXYZ = new double[3];

      // Obtaining coordinates for sidechains
      Atom C = (Atom) newResidue.getAtomNode("C");
      arraycopy(coordsArray[C.getIndex() - 1], 0, c, 0, 3);
      Atom CA = (Atom) newResidue.getAtomNode("CA");
      arraycopy(coordsArray[CA.getIndex() - 1], 0, ca, 0, 3);
      Atom N = (Atom) newResidue.getAtomNode("N");
      arraycopy(coordsArray[N.getIndex() - 1], 0, n, 0, 3);
      Atom BC = (Atom) backResidue.getAtomNode("C");
      arraycopy(coordsArray[BC.getIndex() - 1], 0, bc, 0, 3);
      Atom FN = (Atom) forwardResidue.getAtomNode("N");
      arraycopy(coordsArray[FN.getIndex() - 1], 0, fn, 0, 3);
      /*
      for (Atom backBoneAtom : backBoneAtoms) {
      //    backBoneAtom.setBuilt(true);
          logger.info(String.format("getAtomType().name "+backBoneAtom.getAtomType().name));
          switch (backBoneAtom.getAtomType().name) {
                  case "H":
                                      logger.info(String.format("H getAtomType().name "+backBoneAtom.getAtomType().name));
                      determinedXYZ = determineIntxyz(ca, 1.0, n, 109.5, c, 109.5, -1);
                      coordsArray = fillCoordsArray(backBoneAtom,coordsArray, determinedXYZ);
                      break;
                  case "HA":
                                      logger.info(String.format("HA getAtomType().name "+backBoneAtom.getAtomType().name));
                      determinedXYZ = determineIntxyz(n, 1.0, bc, 119.0, ca, 119.0, 1);
                      coordsArray = fillCoordsArray(backBoneAtom,coordsArray, determinedXYZ);
                      break;
                  case "O":
                                      logger.info(String.format("O getAtomType().name "+backBoneAtom.getAtomType().name));
                      determinedXYZ = determineIntxyz(c, 1.2255, ca, 122.4, n, 180, 0);
                      coordsArray = fillCoordsArray(backBoneAtom,coordsArray, determinedXYZ);
                      break;
                  default:
                      break;
          }
      }
              */
      Atom H = (Atom) newResidue.getAtomNode("H");
      double[] h = coordsArray[H.getIndex() - 1];
      arraycopy(determineIntxyz(n, 1.0, bc, 119.0, ca, 119.0, 1), 0, h, 0, 3); // H
      coordsArray = fillCoordsArray(H, coordsArray, h);

      Atom O = (Atom) newResidue.getAtomNode("O");
      double[] o = coordsArray[O.getIndex() - 1];
      arraycopy(determineIntxyz(c, 1.2255, ca, 122.4, fn, 180, 0), 0, o, 0, 3); // O
      coordsArray = fillCoordsArray(O, coordsArray, o);

      AminoAcid3 name = AminoAcidUtils.AminoAcid3.valueOf(newResidue.getName());

      if (name != AminoAcidUtils.AminoAcid3.GLY) {
        Atom HA = (Atom) newResidue.getAtomNode("HA");
        double[] ha = coordsArray[HA.getIndex() - 1];
        arraycopy(determineIntxyz(ca, 1.0, n, 109.5, c, 109.5, -1), 0, ha, 0, 3); // HA
        coordsArray = fillCoordsArray(HA, coordsArray, ha);
      }

      double rotScale = random.nextDouble();
      switch (name) {
        case GLY:
          {
            Atom HA2 = (Atom) newResidue.getAtomNode("HA2");
            double[] ha2 = coordsArray[HA2.getIndex() - 1];

            Atom HA3 = (Atom) newResidue.getAtomNode("HA3");
            double[] ha3 = coordsArray[HA3.getIndex() - 1];

            arraycopy(determineIntxyz(ca, 1.00, n, 109.5, c, 109.5, 0), 0, ha2, 0, 3); // HA2
            coordsArray = fillCoordsArray(HA2, coordsArray, ha2);

            arraycopy(determineIntxyz(ca, 1.00, n, 109.5, ha2, 109.5, -1), 0, ha3, 0, 3); // HA3
            coordsArray = fillCoordsArray(HA3, coordsArray, ha3);

            break;
          }
        case ALA:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom HB1 = (Atom) newResidue.getAtomNode("HB1");
            double[] hb1 = coordsArray[HB1.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];

            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3); // CB
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(determineIntxyz(cb, 1.11, ca, 109.5, n, 180.0, 0), 0, hb1, 0, 3); // HB1
            coordsArray = fillCoordsArray(HB1, coordsArray, hb1);

            arraycopy(determineIntxyz(cb, 1.11, ca, 109.5, hb1, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(determineIntxyz(cb, 1.11, ca, 109.5, hb1, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            break;
          }
        case VAL:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG1 = (Atom) newResidue.getAtomNode("CG1");
            double[] cg1 = coordsArray[CG1.getIndex() - 1];
            Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
            double[] cg2 = coordsArray[CG2.getIndex() - 1];
            Atom HB = (Atom) newResidue.getAtomNode("HB");
            double[] hb = coordsArray[HB.getIndex() - 1];
            Atom HG11 = (Atom) newResidue.getAtomNode("HG11");
            double[] hg11 = coordsArray[HG11.getIndex() - 1];
            Atom HG12 = (Atom) newResidue.getAtomNode("HG12");
            double[] hg12 = coordsArray[HG12.getIndex() - 1];
            Atom HG13 = (Atom) newResidue.getAtomNode("HG13");
            double[] hg13 = coordsArray[HG13.getIndex() - 1];
            Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
            double[] hg21 = coordsArray[HG21.getIndex() - 1];
            Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
            double[] hg22 = coordsArray[HG22.getIndex() - 1];
            Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
            double[] hg23 = coordsArray[HG23.getIndex() - 1];
            Bond CG_CB = CB.getBond(CG1);
            Bond HB_CB = CB.getBond(HB);
            Bond HG_CG = HG11.getBond(CG1);
            double dCG_CB = CG_CB.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            Angle CG_CB_CA = CG1.getAngle(CB, CA);
            Angle HB_CB_CA = HB.getAngle(CB, CA);
            Angle HG_CG_CB = HG11.getAngle(CG1, CB);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];

            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg1,
                0,
                3); // CG1
            coordsArray = fillCoordsArray(CG1, coordsArray, cg1);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, cg1, 109.5, -1), 0, cg2, 0, 3); // CG2
            coordsArray = fillCoordsArray(CG2, coordsArray, cg2);

            arraycopy(determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg1, 109.4, 1), 0, hb, 0, 3); // HB
            coordsArray = fillCoordsArray(HB, coordsArray, hb);

            arraycopy(
                determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, ca, 180.0, 0), 0, hg11, 0, 3); // HG11
            coordsArray = fillCoordsArray(HG11, coordsArray, hg11);

            arraycopy(
                determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, hg11, 109.4, 1), 0, hg12, 0, 3); // HG12
            coordsArray = fillCoordsArray(HG12, coordsArray, hg12);

            arraycopy(
                determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, hg11, 109.4, -1),
                0,
                hg13,
                0,
                3); // HG13
            coordsArray = fillCoordsArray(HG13, coordsArray, hg13);

            arraycopy(
                determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, ca, 180.0, 0), 0, hg21, 0, 3); // HG21
            coordsArray = fillCoordsArray(HG21, coordsArray, hg21);

            arraycopy(
                determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, hg21, 109.4, 1), 0, hg22, 0, 3); // HG22
            coordsArray = fillCoordsArray(HG22, coordsArray, hg22);

            arraycopy(
                determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, hg21, 109.4, -1),
                0,
                hg23,
                0,
                3); // HG23
            coordsArray = fillCoordsArray(HG23, coordsArray, hg23);

            break;
          }
        case LEU:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG = (Atom) newResidue.getAtomNode("HG");
            double[] hg = coordsArray[HG.getIndex() - 1];
            Atom HD11 = (Atom) newResidue.getAtomNode("HD11");
            double[] hd11 = coordsArray[HD11.getIndex() - 1];
            Atom HD12 = (Atom) newResidue.getAtomNode("HD12");
            double[] hd12 = coordsArray[HD12.getIndex() - 1];
            Atom HD13 = (Atom) newResidue.getAtomNode("HD13");
            double[] hd13 = coordsArray[HD13.getIndex() - 1];
            Atom HD21 = (Atom) newResidue.getAtomNode("HD21");
            double[] hd21 = coordsArray[HD21.getIndex() - 1];
            Atom HD22 = (Atom) newResidue.getAtomNode("HD22");
            double[] hd22 = coordsArray[HD22.getIndex() - 1];
            Atom HD23 = (Atom) newResidue.getAtomNode("HD23");
            double[] hd23 = coordsArray[HD23.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD1.getBond(CG);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG.getBond(CG);
            Bond HD_CD = HD11.getBond(CD1);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD1.getAngle(CG, CB);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG.getAngle(CG, CB);
            Angle HD_CD_CG = HD11.getAngle(CD1, CG);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 109.5, -1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd1, 109.4, 1), 0, hg, 0, 3); // HG
            coordsArray = fillCoordsArray(HG, coordsArray, hg);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, cb, 180.0, 0), 0, hd11, 0, 3); // HD11
            coordsArray = fillCoordsArray(HD11, coordsArray, hd11);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, hd11, 109.4, 1), 0, hd12, 0, 3); // HD12
            coordsArray = fillCoordsArray(HD12, coordsArray, hd12);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, hd11, 109.4, -1),
                0,
                hd13,
                0,
                3); // HD13
            coordsArray = fillCoordsArray(HD13, coordsArray, hd13);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, cb, 180.0, 0), 0, hd21, 0, 3); // HD21
            coordsArray = fillCoordsArray(HD21, coordsArray, hd21);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, hd21, 109.4, 1), 0, hd22, 0, 3); // HD22
            coordsArray = fillCoordsArray(HD22, coordsArray, hd22);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, hd21, 109.4, -1),
                0,
                hd23,
                0,
                3); // HD23
            coordsArray = fillCoordsArray(HD23, coordsArray, hd23);
            break;
          }
        case ILE:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG1 = (Atom) newResidue.getAtomNode("CG1");
            double[] cg1 = coordsArray[CG1.getIndex() - 1];
            Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
            double[] cg2 = coordsArray[CG2.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom HB = (Atom) newResidue.getAtomNode("HB");
            double[] hb = coordsArray[HB.getIndex() - 1];
            Atom HG12 = (Atom) newResidue.getAtomNode("HG12");
            double[] hg12 = coordsArray[HG12.getIndex() - 1];
            Atom HG13 = (Atom) newResidue.getAtomNode("HG13");
            double[] hg13 = coordsArray[HG13.getIndex() - 1];
            Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
            double[] hg21 = coordsArray[HG21.getIndex() - 1];
            Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
            double[] hg22 = coordsArray[HG22.getIndex() - 1];
            Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
            double[] hg23 = coordsArray[HG23.getIndex() - 1];
            Atom HD11 = (Atom) newResidue.getAtomNode("HD11");
            double[] hd11 = coordsArray[HD11.getIndex() - 1];
            Atom HD12 = (Atom) newResidue.getAtomNode("HD12");
            double[] hd12 = coordsArray[HD12.getIndex() - 1];
            Atom HD13 = (Atom) newResidue.getAtomNode("HD13");
            double[] hd13 = coordsArray[HD13.getIndex() - 1];
            Bond CG1_CB = CG1.getBond(CB);
            Bond CG2_CB = CG2.getBond(CB);
            Bond CD1_CG1 = CD1.getBond(CG1);
            Bond HB_CB = HB.getBond(CB);
            Bond HG1_CG = HG12.getBond(CG1);
            Bond HG2_CG = HG22.getBond(CG2);
            Bond HD_CD = HD12.getBond(CD1);
            double dCG1_CB = CG1_CB.bondType.distance;
            double dCG2_CB = CG2_CB.bondType.distance;
            double dCD1_CG1 = CD1_CG1.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG1_CG = HG1_CG.bondType.distance;
            double dHG2_CG = HG2_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            Angle CG1_CB_CA = CG1.getAngle(CB, CA);
            Angle CG2_CB_CA = CG2.getAngle(CB, CA);
            Angle CD1_CG1_CB = CD1.getAngle(CG1, CB);
            Angle HB_CB_CA = HB.getAngle(CB, CA);
            Angle HG1_CG_CB = HG12.getAngle(CG1, CB);
            Angle HG2_CG_CB = HG21.getAngle(CG2, CB);
            Angle HD_CD1_CG1 = HD11.getAngle(CD1, CG1);
            double dCG1_CB_CA = CG1_CB_CA.angleType.angle[CG1_CB_CA.nh];
            double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
            double dCD1_CG1_CB = CD1_CG1_CB.angleType.angle[CD1_CG1_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG1_CG_CB = HG1_CG_CB.angleType.angle[HG1_CG_CB.nh];
            double dHG2_CG_CB = HG2_CG_CB.angleType.angle[HG2_CG_CB.nh];
            double dHD_CD1_CG1 = HD_CD1_CG1.angleType.angle[HD_CD1_CG1.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG1_CB, ca, dCG1_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg1,
                0,
                3); // CG1
            coordsArray = fillCoordsArray(CG1, coordsArray, cg1);

            arraycopy(
                determineIntxyz(cb, dCG2_CB, ca, dCG2_CB_CA, cg1, 109.5, 1), 0, cg2, 0, 3); // CG2
            coordsArray = fillCoordsArray(CG2, coordsArray, cg2);

            arraycopy(
                determineIntxyz(cg1, dCD1_CG1, cb, dCD1_CG1_CB, ca, 180, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg2, 109.4, 1), 0, hb, 0, 3); // HB
            coordsArray = fillCoordsArray(HB, coordsArray, hb);

            arraycopy(
                determineIntxyz(cg1, dHG1_CG, cb, dHG1_CG_CB, cd1, 109.4, 1),
                0,
                hg12,
                0,
                3); // HG12
            coordsArray = fillCoordsArray(HG12, coordsArray, hg12);

            arraycopy(
                determineIntxyz(cg1, dHG1_CG, cb, dHG1_CG_CB, cd1, 109.4, -1),
                0,
                hg13,
                0,
                3); // HG13
            coordsArray = fillCoordsArray(HG13, coordsArray, hg13);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, cg1, 180.0, 0),
                0,
                hg21,
                0,
                3); // HG21
            coordsArray = fillCoordsArray(HG21, coordsArray, hg21);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, hg21, 109.0, 1),
                0,
                hg22,
                0,
                3); // HG22
            coordsArray = fillCoordsArray(HG22, coordsArray, hg22);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, hg21, 109.0, -1),
                0,
                hg23,
                0,
                3); // HG23
            coordsArray = fillCoordsArray(HG23, coordsArray, hg23);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, cb, 180.0, 0),
                0,
                hd11,
                0,
                3); // HD11
            coordsArray = fillCoordsArray(HD11, coordsArray, hd11);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, hd11, 109.0, 1),
                0,
                hd12,
                0,
                3); // HD12
            coordsArray = fillCoordsArray(HD12, coordsArray, hd12);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, hd11, 109.0, -1),
                0,
                hd13,
                0,
                3); // HD13
            coordsArray = fillCoordsArray(HD13, coordsArray, hd13);
            break;
          }
        case SER:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom OG = (Atom) newResidue.getAtomNode("OG");
            double[] og = coordsArray[OG.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG = (Atom) newResidue.getAtomNode("HG");
            double[] hg = coordsArray[HG.getIndex() - 1];
            Bond OG_CB = OG.getBond(CB);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_OG = HG.getBond(OG);
            double dOG_CB = OG_CB.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_OG = HG_OG.bondType.distance;
            Angle OG_CB_CA = OG.getAngle(CB, CA);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_OG_CB = HG.getAngle(OG, CB);
            double dOG_CB_CA = OG_CB_CA.angleType.angle[OG_CB_CA.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_OG_CB = HG_OG_CB.angleType.angle[HG_OG_CB.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dOG_CB, ca, dOG_CB_CA, n, rotScale * 180.0, 0),
                0,
                og,
                0,
                3); // OG
            coordsArray = fillCoordsArray(OG, coordsArray, og);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og, 106.7, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og, 106.7, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(determineIntxyz(og, dHG_OG, cb, dHG_OG_CB, ca, 180.0, 0), 0, hg, 0, 3); // HG
            coordsArray = fillCoordsArray(HG, coordsArray, hg);
            break;
          }
        case THR:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom OG1 = (Atom) newResidue.getAtomNode("OG1");
            double[] og1 = coordsArray[OG1.getIndex() - 1];
            Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
            double[] cg2 = coordsArray[CG2.getIndex() - 1];
            Atom HB = (Atom) newResidue.getAtomNode("HB");
            double[] hb = coordsArray[HB.getIndex() - 1];
            Atom HG1 = (Atom) newResidue.getAtomNode("HG1");
            double[] hg1 = coordsArray[HG1.getIndex() - 1];
            Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
            double[] hg21 = coordsArray[HG21.getIndex() - 1];
            Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
            double[] hg22 = coordsArray[HG22.getIndex() - 1];
            Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
            double[] hg23 = coordsArray[HG23.getIndex() - 1];
            Bond OG1_CB = OG1.getBond(CB);
            Bond CG2_CB = CG2.getBond(CB);
            Bond HB_CB = HB.getBond(CB);
            Bond HG1_OG1 = HG1.getBond(OG1);
            Bond HG2_CG2 = HG21.getBond(CG2);
            double dOG1_CB = OG1_CB.bondType.distance;
            double dCG2_CB = CG2_CB.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG1_OG1 = HG1_OG1.bondType.distance;
            double dHG2_CG2 = HG2_CG2.bondType.distance;
            Angle OG1_CB_CA = OG1.getAngle(CB, CA);
            Angle CG2_CB_CA = CG2.getAngle(CB, CA);
            Angle HB_CB_CA = HB.getAngle(CB, CA);
            Angle HG1_OG1_CB = HG1.getAngle(OG1, CB);
            Angle HG2_CG2_CB = HG21.getAngle(CG2, CB);
            double dOG1_CB_CA = OG1_CB_CA.angleType.angle[OG1_CB_CA.nh];
            double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG1_OG1_CB = HG1_OG1_CB.angleType.angle[HG1_OG1_CB.nh];
            double dHG2_CG2_CB = HG2_CG2_CB.angleType.angle[HG2_CG2_CB.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dOG1_CB, ca, dOG1_CB_CA, n, rotScale * 180.0, 0),
                0,
                og1,
                0,
                3); // OG1
            coordsArray = fillCoordsArray(OG1, coordsArray, og1);

            arraycopy(
                determineIntxyz(cb, dCG2_CB, ca, dCG2_CB_CA, og1, 107.7, 1), 0, cg2, 0, 3); // CG2
            coordsArray = fillCoordsArray(CG2, coordsArray, cg2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og1, 106.7, -1), 0, hb, 0, 3); // HB
            coordsArray = fillCoordsArray(HB, coordsArray, hb);

            arraycopy(
                determineIntxyz(og1, dHG1_OG1, cb, dHG1_OG1_CB, ca, 180.0, 0), 0, hg1, 0, 3); // HG1
            coordsArray = fillCoordsArray(HG1, coordsArray, hg1);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, ca, 180.0, 0),
                0,
                hg21,
                0,
                3); // HG21
            coordsArray = fillCoordsArray(HG21, coordsArray, hg21);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, hg21, 109.0, 1),
                0,
                hg22,
                0,
                3); // HG22
            coordsArray = fillCoordsArray(HG22, coordsArray, hg22);

            arraycopy(
                determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, hg21, 109.0, -1),
                0,
                hg23,
                0,
                3); // HG23
            coordsArray = fillCoordsArray(HG23, coordsArray, hg23);
            break;
          }
        case CYS:
        case CYX:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom SG = (Atom) newResidue.getAtomNode("SG");
            double[] sg = coordsArray[SG.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG = (Atom) newResidue.getAtomNode("HG");
            double[] hg = coordsArray[HG.getIndex() - 1];
            Bond SG_CB = SG.getBond(CB);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_SG = HG.getBond(SG);
            double dSG_CB = SG_CB.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_SG = HG_SG.bondType.distance;
            Angle SG_CB_CA = SG.getAngle(CB, CA);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_SG_CB = HG.getAngle(SG, CB);
            double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_SG_CB = HG_SG_CB.angleType.angle[HG_SG_CB.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dSG_CB, ca, dSG_CB_CA, n, rotScale * 180.0, 0),
                0,
                sg,
                0,
                3); // SG
            coordsArray = fillCoordsArray(SG, coordsArray, sg);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(determineIntxyz(sg, dHG_SG, cb, dHG_SG_CB, ca, 180.0, 0), 0, hg, 0, 3); // HG
            coordsArray = fillCoordsArray(HG, coordsArray, hg);
            break;
          }
        case CYD:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom SG = (Atom) newResidue.getAtomNode("SG");
            double[] sg = coordsArray[SG.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Bond SG_CB = SG.getBond(CB);
            Bond HB_CB = HB2.getBond(CB);
            double dSG_CB = SG_CB.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            Angle SG_CB_CA = SG.getAngle(CB, CA);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dSG_CB, ca, dSG_CB_CA, n, rotScale * 180.0, 0),
                0,
                sg,
                0,
                3); // SG
            coordsArray = fillCoordsArray(SG, coordsArray, sg);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);
            break;
          }
        case PHE:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
            double[] ce2 = coordsArray[CE2.getIndex() - 1];
            Atom CZ = (Atom) newResidue.getAtomNode("CZ");
            double[] cz = coordsArray[CZ.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Atom HZ = (Atom) newResidue.getAtomNode("HZ");
            double[] hz = coordsArray[HZ.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD1.getBond(CG);
            Bond CE_CD = CE1.getBond(CD1);
            Bond CZ_CE1 = CZ.getBond(CE1);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD_CD = HD1.getBond(CD1);
            Bond HE_CE = HE1.getBond(CE1);
            Bond HZ_CZ = HZ.getBond(CZ);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dCE_CD = CE_CD.bondType.distance;
            double dCZ_CE1 = CZ_CE1.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            double dHZ_CZ = HZ_CZ.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD1.getAngle(CG, CB);
            Angle CE_CD_CG = CE1.getAngle(CD1, CG);
            Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD_CD1_CG = HD1.getAngle(CD1, CG);
            Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
            Angle HZ_CZ_CE1 = HZ.getAngle(CZ, CE1);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
            double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD_CD1_CG = HD_CD1_CG.angleType.angle[HD_CD1_CG.nh];
            double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
            double dHZ_CZ_CE1 = HZ_CZ_CE1.angleType.angle[HZ_CZ_CE1.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce2, 0, 3); // CE2
            coordsArray = fillCoordsArray(CE2, coordsArray, ce2);

            arraycopy(
                determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0), 0, cz, 0, 3); // CZ
            coordsArray = fillCoordsArray(CZ, coordsArray, cz);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD1_CG, ce1, 120.0, 1), 0, hd1, 0, 3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD1_CG, ce2, 120.0, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1), 0, he1, 0, 3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);

            arraycopy(
                determineIntxyz(cz, dHZ_CZ, ce1, dHZ_CZ_CE1, ce2, 120.0, 1), 0, hz, 0, 3); // HZ
            coordsArray = fillCoordsArray(HZ, coordsArray, hz);
            break;
          }
        case PRO:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
            double[] hd3 = coordsArray[HD3.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HD_CD = HD2.getBond(CD);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HD_CD_CG = HD2.getAngle(CD, CG);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 30.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, n, 109.4, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, n, 109.4, -1), 0, hd3, 0, 3); // HD3
            coordsArray = fillCoordsArray(HD3, coordsArray, hd3);
            break;
          }
        case TYR:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
            double[] ce2 = coordsArray[CE2.getIndex() - 1];
            Atom CZ = (Atom) newResidue.getAtomNode("CZ");
            double[] cz = coordsArray[CZ.getIndex() - 1];
            Atom OH = (Atom) newResidue.getAtomNode("OH");
            double[] oh = coordsArray[OH.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Atom HH = (Atom) newResidue.getAtomNode("HH");
            double[] hh = coordsArray[HH.getIndex() - 1];
            Bond CB_CA = CB.getBond(CA);
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD1.getBond(CG);
            Bond CE_CD = CE1.getBond(CD1);
            Bond CZ_CE1 = CZ.getBond(CE1);
            Bond OH_CZ = OH.getBond(CZ);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD_CD = HD1.getBond(CD1);
            Bond HE_CE = HE1.getBond(CE1);
            Bond HH_OH = HH.getBond(OH);
            double dCB_CA = CB_CA.bondType.distance;
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dCE_CD = CE_CD.bondType.distance;
            double dCZ_CE1 = CZ_CE1.bondType.distance;
            double dOH_CZ = OH_CZ.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            double dHH_OH = HH_OH.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD1.getAngle(CG, CB);
            Angle CE_CD_CG = CE1.getAngle(CD1, CG);
            Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
            Angle OH_CZ_CE2 = OH.getAngle(CZ, CE2);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD_CD_CG = HD1.getAngle(CD1, CG);
            Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
            Angle HH_OH_CZ = HH.getAngle(OH, CZ);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
            double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
            double dOH_CZ_CE2 = OH_CZ_CE2.angleType.angle[OH_CZ_CE2.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
            double dHH_OH_CZ = HH_OH_CZ.angleType.angle[HH_OH_CZ.nh];
            arraycopy(determineIntxyz(ca, dCB_CA, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 90.0, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce2, 0, 3); // CE2
            coordsArray = fillCoordsArray(CE2, coordsArray, ce2);

            arraycopy(
                determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0), 0, cz, 0, 3); // CZ
            coordsArray = fillCoordsArray(CZ, coordsArray, cz);

            arraycopy(
                determineIntxyz(cz, dOH_CZ, ce2, dOH_CZ_CE2, ce1, 120.0, 1), 0, oh, 0, 3); // OH
            coordsArray = fillCoordsArray(OH, coordsArray, oh);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, ce1, 120.0, 1), 0, hd1, 0, 3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, ce2, 120.0, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1), 0, he1, 0, 3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);

            arraycopy(determineIntxyz(oh, dHH_OH, cz, dHH_OH_CZ, ce2, 0.0, 0), 0, hh, 0, 3); // HH
            coordsArray = fillCoordsArray(HH, coordsArray, hh);
            break;
          }
        case TYD:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
            double[] ce2 = coordsArray[CE2.getIndex() - 1];
            Atom CZ = (Atom) newResidue.getAtomNode("CZ");
            double[] cz = coordsArray[CZ.getIndex() - 1];
            Atom OH = (Atom) newResidue.getAtomNode("OH");
            double[] oh = coordsArray[OH.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD1.getBond(CG);
            Bond CE_CD = CE1.getBond(CD1);
            Bond CZ_CE1 = CZ.getBond(CE1);
            Bond OH_CZ = OH.getBond(CZ);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD_CD = HD1.getBond(CD1);
            Bond HE_CE = HE1.getBond(CE1);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dCE_CD = CE_CD.bondType.distance;
            double dCZ_CE1 = CZ_CE1.bondType.distance;
            double dOH_CZ = OH_CZ.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD1.getAngle(CG, CB);
            Angle CE_CD_CG = CE1.getAngle(CD1, CG);
            Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
            Angle OH_CZ_CE2 = OH.getAngle(CZ, CE2);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD_CD_CG = HD1.getAngle(CD1, CG);
            Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
            double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
            double dOH_CZ_CE2 = OH_CZ_CE2.angleType.angle[OH_CZ_CE2.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 90.0, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(
                determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0), 0, ce2, 0, 3); // CE2
            coordsArray = fillCoordsArray(CE2, coordsArray, ce2);

            arraycopy(
                determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0), 0, cz, 0, 3); // CZ
            coordsArray = fillCoordsArray(CZ, coordsArray, cz);

            arraycopy(
                determineIntxyz(cz, dOH_CZ, ce2, dOH_CZ_CE2, ce1, 120.0, 1), 0, oh, 0, 3); // OH
            coordsArray = fillCoordsArray(OH, coordsArray, oh);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, ce1, 120.0, 1), 0, hd1, 0, 3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, ce2, 120.0, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1), 0, he1, 0, 3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);
            break;
          }
        case TRP:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
            double[] cd1 = coordsArray[CD1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom NE1 = (Atom) newResidue.getAtomNode("NE1");
            double[] ne1 = coordsArray[NE1.getIndex() - 1];
            Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
            double[] ce2 = coordsArray[CE2.getIndex() - 1];
            Atom CE3 = (Atom) newResidue.getAtomNode("CE3");
            double[] ce3 = coordsArray[CE3.getIndex() - 1];
            Atom CZ2 = (Atom) newResidue.getAtomNode("CZ2");
            double[] cz2 = coordsArray[CZ2.getIndex() - 1];
            Atom CZ3 = (Atom) newResidue.getAtomNode("CZ3");
            double[] cz3 = coordsArray[CZ3.getIndex() - 1];
            Atom CH2 = (Atom) newResidue.getAtomNode("CH2");
            double[] ch2 = coordsArray[CH2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
            double[] he3 = coordsArray[HE3.getIndex() - 1];
            Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
            double[] hz2 = coordsArray[HZ2.getIndex() - 1];
            Atom HZ3 = (Atom) newResidue.getAtomNode("HZ3");
            double[] hz3 = coordsArray[HZ3.getIndex() - 1];
            Atom HH2 = (Atom) newResidue.getAtomNode("HH2");
            double[] hh2 = coordsArray[HH2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD1_CG = CD1.getBond(CG);
            Bond CD2_CG = CD2.getBond(CG);
            Bond NE1_CD1 = NE1.getBond(CD1);
            Bond CE2_NE1 = CE2.getBond(NE1);
            Bond CE3_CD2 = CE3.getBond(CD2);
            Bond CZ2_CE2 = CZ2.getBond(CE2);
            Bond CZ3_CE3 = CZ3.getBond(CE3);
            Bond CH2_CZ2 = CH2.getBond(CZ2);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD1_CD1 = HD1.getBond(CD1);
            Bond HE1_NE1 = HE1.getBond(NE1);
            Bond HE3_CE3 = HE3.getBond(CE3);
            Bond HZ2_CZ2 = HZ2.getBond(CZ2);
            Bond HZ3_CZ3 = HZ3.getBond(CZ3);
            Bond HH2_CH2 = HH2.getBond(CH2);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD1_CG = CD1_CG.bondType.distance;
            double dCD2_CG = CD2_CG.bondType.distance;
            double dNE1_CD1 = NE1_CD1.bondType.distance;
            double dCE2_NE1 = CE2_NE1.bondType.distance;
            double dCE3_CD2 = CE3_CD2.bondType.distance;
            double dCZ2_CE2 = CZ2_CE2.bondType.distance;
            double dCZ3_CE3 = CZ3_CE3.bondType.distance;
            double dCH2_CZ2 = CH2_CZ2.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD1_CD1 = HD1_CD1.bondType.distance;
            double dHE1_NE1 = HE1_NE1.bondType.distance;
            double dHE3_CE3 = HE3_CE3.bondType.distance;
            double dHZ2_CZ2 = HZ2_CZ2.bondType.distance;
            double dHZ3_CZ3 = HZ3_CZ3.bondType.distance;
            double dHH2_CH2 = HH2_CH2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD1_CG_CB = CD1.getAngle(CG, CB);
            Angle CD2_CG_CB = CD2.getAngle(CG, CB);
            Angle NE1_CD1_CG = NE1.getAngle(CD1, CG);
            Angle CE2_NE1_CD1 = CE2.getAngle(NE1, CD1);
            Angle CE3_CD2_CE2 = CE3.getAngle(CD2, CE2);
            Angle CZ2_CE2_CD2 = CZ2.getAngle(CE2, CD2);
            Angle CZ3_CE3_CD2 = CZ3.getAngle(CE3, CD2);
            Angle CH2_CZ2_CE2 = CH2.getAngle(CZ2, CE2);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD1_CD1_CG = HD1.getAngle(CD1, CG);
            Angle HE1_NE1_CD1 = HE1.getAngle(NE1, CD1);
            Angle HE3_CE3_CD2 = HE3.getAngle(CE3, CD2);
            Angle HZ2_CZ2_CE2 = HZ2.getAngle(CZ2, CE2);
            Angle HZ3_CZ3_CE3 = HZ3.getAngle(CZ3, CH2);
            Angle HH2_CH2_CZ2 = HH2.getAngle(CH2, CZ3);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD1_CG_CB = CD1_CG_CB.angleType.angle[CD1_CG_CB.nh];
            double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
            double dNE1_CD1_CG = NE1_CD1_CG.angleType.angle[NE1_CD1_CG.nh];
            double dCE2_NE1_CD1 = CE2_NE1_CD1.angleType.angle[CE2_NE1_CD1.nh];
            double dCE3_CD2_CE2 = CE3_CD2_CE2.angleType.angle[CE3_CD2_CE2.nh];
            double dCZ2_CE2_CD2 = CZ2_CE2_CD2.angleType.angle[CZ2_CE2_CD2.nh];
            double dCZ3_CE3_CD2 = CZ3_CE3_CD2.angleType.angle[CZ3_CE3_CD2.nh];
            double dCH2_CZ2_CE2 = CH2_CZ2_CE2.angleType.angle[CH2_CZ2_CE2.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD1_CD1_CG = HD1_CD1_CG.angleType.angle[HD1_CD1_CG.nh];
            double dHE1_NE1_CD1 = HE1_NE1_CD1.angleType.angle[HE1_NE1_CD1.nh];
            double dHE3_CE3_CD2 = HE3_CE3_CD2.angleType.angle[HE3_CE3_CD2.nh];
            double dHZ2_CZ2_CE2 = HZ2_CZ2_CE2.angleType.angle[HZ2_CZ2_CE2.nh];
            double dHZ3_CZ3_CE3 = HZ3_CZ3_CE3.angleType.angle[HZ3_CZ3_CE3.nh];
            double dHH2_CH2_CZ2 = HH2_CH2_CZ2.angleType.angle[HH2_CH2_CZ2.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dCD1_CG, cb, dCD1_CG_CB, ca, -90.0, 0), 0, cd1, 0, 3); // CD1
            coordsArray = fillCoordsArray(CD1, coordsArray, cd1);

            arraycopy(
                determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, cd1, 108.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(
                determineIntxyz(cd1, dNE1_CD1, cg, dNE1_CD1_CG, cd2, 0.0, 0), 0, ne1, 0, 3); // NE1
            coordsArray = fillCoordsArray(NE1, coordsArray, ne1);

            arraycopy(
                determineIntxyz(ne1, dCE2_NE1, cd1, dCE2_NE1_CD1, cg, 0.0, 0), 0, ce2, 0, 3); // CE2
            coordsArray = fillCoordsArray(CE2, coordsArray, ce2);

            arraycopy(
                determineIntxyz(cd2, dCE3_CD2, ce2, dCE3_CD2_CE2, ne1, 180.0, 0),
                0,
                ce3,
                0,
                3); // CE3
            coordsArray = fillCoordsArray(CE3, coordsArray, ce3);

            arraycopy(
                determineIntxyz(ce2, dCZ2_CE2, cd2, dCZ2_CE2_CD2, ce3, 0.0, 0),
                0,
                cz2,
                0,
                3); // CZ2
            coordsArray = fillCoordsArray(CZ2, coordsArray, cz2);

            arraycopy(
                determineIntxyz(ce3, dCZ3_CE3, cd2, dCZ3_CE3_CD2, ce2, 0.0, 0),
                0,
                cz3,
                0,
                3); // CZ3
            coordsArray = fillCoordsArray(CZ3, coordsArray, cz3);

            arraycopy(
                determineIntxyz(cz2, dCH2_CZ2, ce2, dCH2_CZ2_CE2, cd2, 0.0, 0),
                0,
                ch2,
                0,
                3); // CH2
            coordsArray = fillCoordsArray(CH2, coordsArray, ch2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cd1, dHD1_CD1, cg, dHD1_CD1_CG, ne1, 126.0, 1),
                0,
                hd1,
                0,
                3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(ne1, dHE1_NE1, cd1, dHE1_NE1_CD1, ce2, 126.0, 1),
                0,
                he1,
                0,
                3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ce3, dHE3_CE3, cd2, dHE3_CE3_CD2, cz3, 120.0, 1),
                0,
                he3,
                0,
                3); // HE3
            coordsArray = fillCoordsArray(HE3, coordsArray, he3);

            arraycopy(
                determineIntxyz(cz2, dHZ2_CZ2, ce2, dHZ2_CZ2_CE2, ch2, 120.0, 1),
                0,
                hz2,
                0,
                3); // HZ2
            coordsArray = fillCoordsArray(HZ2, coordsArray, hz2);

            arraycopy(
                determineIntxyz(cz3, dHZ3_CZ3, ce3, dHZ3_CZ3_CE3, ch2, 120.0, 1),
                0,
                hz3,
                0,
                3); // HZ3
            coordsArray = fillCoordsArray(HZ3, coordsArray, hz3);

            arraycopy(
                determineIntxyz(ch2, dHH2_CH2, cz2, dHH2_CH2_CZ2, cz3, 120.0, 1),
                0,
                hh2,
                0,
                3); // HH2
            coordsArray = fillCoordsArray(HH2, coordsArray, hh2);
            break;
          }
        case HIS:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
            double[] nd1 = coordsArray[ND1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
            double[] ne2 = coordsArray[NE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond ND1_CG = ND1.getBond(CG);
            Bond CD2_CG = CD2.getBond(CG);
            Bond CE1_ND1 = CE1.getBond(ND1);
            Bond NE2_CD2 = NE2.getBond(CD2);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD1_ND1 = HD1.getBond(ND1);
            Bond HD2_CD2 = HD2.getBond(CD2);
            Bond HE1_CE1 = HE1.getBond(CE1);
            Bond HE2_NE2 = HE2.getBond(NE2);
            double dCG_CB = CG_CB.bondType.distance;
            double dND1_CG = ND1_CG.bondType.distance;
            double dCD2_CG = CD2_CG.bondType.distance;
            double dCE1_ND1 = CE1_ND1.bondType.distance;
            double dNE2_CD2 = NE2_CD2.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD1_ND1 = HD1_ND1.bondType.distance;
            double dHD2_CD2 = HD2_CD2.bondType.distance;
            double dHE1_CE1 = HE1_CE1.bondType.distance;
            double dHE2_NE2 = HE2_NE2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle ND1_CG_CB = ND1.getAngle(CG, CB);
            Angle CD2_CG_CB = CD2.getAngle(CG, CB);
            Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
            Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD1_ND1_CG = HD1.getAngle(ND1, CG);
            Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
            Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
            Angle HE2_NE2_CD2 = HE2.getAngle(NE2, CD2);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
            double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
            double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
            double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD1_ND1_CG = HD1_ND1_CG.angleType.angle[HD1_ND1_CG.nh];
            double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
            double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
            double dHE2_NE2_CD2 = HE2_NE2_CD2.angleType.angle[HE2_NE2_CD2.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0), 0, nd1, 0, 3); // ND1
            coordsArray = fillCoordsArray(ND1, coordsArray, nd1);

            arraycopy(
                determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(
                determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(
                determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0), 0, ne2, 0, 3); // NE2
            coordsArray = fillCoordsArray(NE2, coordsArray, ne2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(nd1, dHD1_ND1, cg, dHD1_ND1_CG, cb, 0.0, 0), 0, hd1, 0, 3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1),
                0,
                hd2,
                0,
                3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1),
                0,
                he1,
                0,
                3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ne2, dHE2_NE2, cd2, dHE2_NE2_CD2, ce1, 126.0, 1),
                0,
                he2,
                0,
                3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);
            break;
          }
        case HID:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
            double[] nd1 = coordsArray[ND1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
            double[] ne2 = coordsArray[NE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
            double[] hd1 = coordsArray[HD1.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond ND1_CG = ND1.getBond(CG);
            Bond CD2_CG = CD2.getBond(CG);
            Bond CE1_ND1 = CE1.getBond(ND1);
            Bond NE2_CD2 = NE2.getBond(CD2);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD1_ND1 = HD1.getBond(ND1);
            Bond HD2_CD2 = HD2.getBond(CD2);
            Bond HE1_CE1 = HE1.getBond(CE1);
            double dCG_CB = CG_CB.bondType.distance;
            double dND1_CG = ND1_CG.bondType.distance;
            double dCD2_CG = CD2_CG.bondType.distance;
            double dCE1_ND1 = CE1_ND1.bondType.distance;
            double dNE2_CD2 = NE2_CD2.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD1_ND1 = HD1_ND1.bondType.distance;
            double dHD2_CD2 = HD2_CD2.bondType.distance;
            double dHE1_CE1 = HE1_CE1.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle ND1_CG_CB = ND1.getAngle(CG, CB);
            Angle CD2_CG_CB = CD2.getAngle(CG, CB);
            Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
            Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD1_ND1_CG = HD1.getAngle(ND1, CG);
            Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
            Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
            double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
            double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
            double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD1_ND1_CG = HD1_ND1_CG.angleType.angle[HD1_ND1_CG.nh];
            double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
            double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0), 0, nd1, 0, 3); // ND1
            coordsArray = fillCoordsArray(ND1, coordsArray, nd1);

            arraycopy(
                determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(
                determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(
                determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0), 0, ne2, 0, 3); // NE2
            coordsArray = fillCoordsArray(NE2, coordsArray, ne2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(nd1, dHD1_ND1, cg, dHD1_ND1_CG, cb, 0.0, 0), 0, hd1, 0, 3); // HD1
            coordsArray = fillCoordsArray(HD1, coordsArray, hd1);

            arraycopy(
                determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1),
                0,
                hd2,
                0,
                3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1),
                0,
                he1,
                0,
                3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);
            break;
          }
        case HIE:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
            double[] nd1 = coordsArray[ND1.getIndex() - 1];
            Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
            double[] cd2 = coordsArray[CD2.getIndex() - 1];
            Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
            double[] ce1 = coordsArray[CE1.getIndex() - 1];
            Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
            double[] ne2 = coordsArray[NE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond ND1_CG = ND1.getBond(CG);
            Bond CD2_CG = CD2.getBond(CG);
            Bond CE1_ND1 = CE1.getBond(ND1);
            Bond NE2_CD2 = NE2.getBond(CD2);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD2_CD2 = HD2.getBond(CD2);
            Bond HE1_CE1 = HE1.getBond(CE1);
            Bond HE2_NE2 = HE2.getBond(NE2);
            double dCG_CB = CG_CB.bondType.distance;
            double dND1_CG = ND1_CG.bondType.distance;
            double dCD2_CG = CD2_CG.bondType.distance;
            double dCE1_ND1 = CE1_ND1.bondType.distance;
            double dNE2_CD2 = NE2_CD2.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD2_CD2 = HD2_CD2.bondType.distance;
            double dHE1_CE1 = HE1_CE1.bondType.distance;
            double dHE2_NE2 = HE2_NE2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle ND1_CG_CB = ND1.getAngle(CG, CB);
            Angle CD2_CG_CB = CD2.getAngle(CG, CB);
            Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
            Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
            Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
            Angle HE2_NE2_CD2 = HE2.getAngle(NE2, CD2);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
            double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
            double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
            double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
            double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
            double dHE2_NE2_CD2 = HE2_NE2_CD2.angleType.angle[HE2_NE2_CD2.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0), 0, nd1, 0, 3); // ND1
            coordsArray = fillCoordsArray(ND1, coordsArray, nd1);

            arraycopy(
                determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1), 0, cd2, 0, 3); // CD2
            coordsArray = fillCoordsArray(CD2, coordsArray, cd2);

            arraycopy(
                determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0), 0, ce1, 0, 3); // CE1
            coordsArray = fillCoordsArray(CE1, coordsArray, ce1);

            arraycopy(
                determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0), 0, ne2, 0, 3); // NE2
            coordsArray = fillCoordsArray(NE2, coordsArray, ne2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1),
                0,
                hd2,
                0,
                3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1),
                0,
                he1,
                0,
                3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ne2, dHE2_NE2, cd2, dHE2_NE2_CD2, ce1, 126.0, 1),
                0,
                he2,
                0,
                3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);
            break;
          }
        case ASP:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
            double[] od1 = coordsArray[OD1.getIndex() - 1];
            Atom OD2 = (Atom) newResidue.getAtomNode("OD2");
            double[] od2 = coordsArray[OD2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond OD1_CG = OD1.getBond(CG);
            Bond OD2_CG = OD2.getBond(CG);
            Bond HB_CB = HB2.getBond(CB);
            double dCG_CB = CG_CB.bondType.distance;
            double dOD1_CG = OD1_CG.bondType.distance;
            double dOD2_CG = OD2_CG.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle OD1_CG_CB = OD1.getAngle(CG, CB);
            Angle OD2_CG_CB = OD2.getAngle(CG, CB);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
            double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0), 0, od1, 0, 3); // OD1
            coordsArray = fillCoordsArray(OD1, coordsArray, od1);

            arraycopy(
                determineIntxyz(cg, dOD2_CG, cb, dOD2_CG_CB, od1, 126.0, 1), 0, od2, 0, 3); // OD2
            coordsArray = fillCoordsArray(OD2, coordsArray, od2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);
            break;
          }
        case ASH:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
            double[] od1 = coordsArray[OD1.getIndex() - 1];
            Atom OD2 = (Atom) newResidue.getAtomNode("OD2");
            double[] od2 = coordsArray[OD2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond OD1_CG = OD1.getBond(CG);
            Bond OD2_CG = OD2.getBond(CG);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD2_OD2 = HD2.getBond(OD2);
            double dCG_CB = CG_CB.bondType.distance;
            double dOD1_CG = OD1_CG.bondType.distance;
            double dOD2_CG = OD2_CG.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD2_OD2 = HD2_OD2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle OD1_CG_CB = OD1.getAngle(CG, CB);
            Angle OD2_CG_CB = OD2.getAngle(CG, CB);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD2_OD2_CG = HD2.getAngle(OD2, CG);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
            double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD2_OD2_CG = HD2_OD2_CG.angleType.angle[HD2_OD2_CG.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0), 0, od1, 0, 3); // OD1
            coordsArray = fillCoordsArray(OD1, coordsArray, od1);

            arraycopy(
                determineIntxyz(cg, dOD2_CG, cb, dOD2_CG_CB, od1, 126.0, 1), 0, od2, 0, 3); // OD2
            coordsArray = fillCoordsArray(OD2, coordsArray, od2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(od2, dHD2_OD2, cg, dHD2_OD2_CG, od1, 0.0, 0), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);
            break;
          }
        case ASN:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
            double[] od1 = coordsArray[OD1.getIndex() - 1];
            Atom ND2 = (Atom) newResidue.getAtomNode("ND2");
            double[] nd2 = coordsArray[ND2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HD21 = (Atom) newResidue.getAtomNode("HD21");
            double[] hd21 = coordsArray[HD21.getIndex() - 1];
            Atom HD22 = (Atom) newResidue.getAtomNode("HD22");
            double[] hd22 = coordsArray[HD22.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond OD1_CG = OD1.getBond(CG);
            Bond ND2_CG = ND2.getBond(CG);
            Bond HB_CB = HB2.getBond(CB);
            Bond HD2_ND2 = HD21.getBond(ND2);
            double dCG_CB = CG_CB.bondType.distance;
            double dOD1_CG = OD1_CG.bondType.distance;
            double dND2_CG = ND2_CG.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHD2_ND2 = HD2_ND2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle OD1_CG_CB = OD1.getAngle(CG, CB);
            Angle ND2_CG_CB = ND2.getAngle(CG, CB);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HD2_ND2_CG = HD21.getAngle(ND2, CG);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
            double dND2_CG_CB = ND2_CG_CB.angleType.angle[ND2_CG_CB.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHD2_ND2_CG = HD2_ND2_CG.angleType.angle[HD2_ND2_CG.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(
                determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0), 0, od1, 0, 3); // OD1
            coordsArray = fillCoordsArray(OD1, coordsArray, od1);

            arraycopy(
                determineIntxyz(cg, dND2_CG, cb, dND2_CG_CB, od1, 124.0, 1), 0, nd2, 0, 3); // ND2
            coordsArray = fillCoordsArray(ND2, coordsArray, nd2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(nd2, dHD2_ND2, cg, dHD2_ND2_CG, cb, 0.0, 0), 0, hd21, 0, 3); // HD21
            coordsArray = fillCoordsArray(HD21, coordsArray, hd21);

            arraycopy(
                determineIntxyz(nd2, dHD2_ND2, cg, dHD2_ND2_CG, hd21, 120.0, 1),
                0,
                hd22,
                0,
                3); // HD22
            coordsArray = fillCoordsArray(HD22, coordsArray, hd22);
            break;
          }
        case GLU:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
            double[] oe1 = coordsArray[OE1.getIndex() - 1];
            Atom OE2 = (Atom) newResidue.getAtomNode("OE2");
            double[] oe2 = coordsArray[OE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond OE1_CD = OE1.getBond(CD);
            Bond OE2_CD = OE2.getBond(CD);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dOE1_CD = OE1_CD.bondType.distance;
            double dOE2_CD = OE2_CD.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle OE1_CD_CG = OE1.getAngle(CD, CG);
            Angle OE2_CD_CG = OE2.getAngle(CD, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
            double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(
                determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0), 0, oe1, 0, 3); // OE1
            coordsArray = fillCoordsArray(OE1, coordsArray, oe1);

            arraycopy(
                determineIntxyz(cd, dOE2_CD, cg, dOE2_CD_CG, oe1, 126.0, 1), 0, oe2, 0, 3); // OE2
            coordsArray = fillCoordsArray(OE2, coordsArray, oe2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);
            break;
          }
        case GLH:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
            double[] oe1 = coordsArray[OE1.getIndex() - 1];
            Atom OE2 = (Atom) newResidue.getAtomNode("OE2");
            double[] oe2 = coordsArray[OE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond OE1_CD = OE1.getBond(CD);
            Bond OE2_CD = OE2.getBond(CD);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HE2_OE2 = HE2.getBond(OE2);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dOE1_CD = OE1_CD.bondType.distance;
            double dOE2_CD = OE2_CD.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHE2_OE2 = HE2_OE2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle OE1_CD_CG = OE1.getAngle(CD, CG);
            Angle OE2_CD_CG = OE2.getAngle(CD, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HE2_OE2_CD = HE2.getAngle(OE2, CD);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
            double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHE2_OE2_CD = HE2_OE2_CD.angleType.angle[HE2_OE2_CD.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(
                determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0), 0, oe1, 0, 3); // OE1
            coordsArray = fillCoordsArray(OE1, coordsArray, oe1);

            arraycopy(
                determineIntxyz(cd, dOE2_CD, cg, dOE2_CD_CG, oe1, 126.0, 1), 0, oe2, 0, 3); // OE2
            coordsArray = fillCoordsArray(OE2, coordsArray, oe2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(oe2, dHE2_OE2, cd, dHE2_OE2_CD, oe1, 0.0, 0), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);
            break;
          }
        case GLN:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
            double[] oe1 = coordsArray[OE1.getIndex() - 1];
            Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
            double[] ne2 = coordsArray[NE2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HE21 = (Atom) newResidue.getAtomNode("HE21");
            double[] he21 = coordsArray[HE21.getIndex() - 1];
            Atom HE22 = (Atom) newResidue.getAtomNode("HE22");
            double[] he22 = coordsArray[HE22.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond OE1_CD = OE1.getBond(CD);
            Bond NE2_CD = NE2.getBond(CD);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HE2_NE2 = HE21.getBond(NE2);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dOE1_CD = OE1_CD.bondType.distance;
            double dNE2_CD = NE2_CD.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHE2_NE2 = HE2_NE2.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle OE1_CD_CG = OE1.getAngle(CD, CG);
            Angle NE2_CD_CG = NE2.getAngle(CD, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HE2_NE2_CD = HE21.getAngle(NE2, CD);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
            double dNE2_CD_CG = NE2_CD_CG.angleType.angle[NE2_CD_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHE2_NE2_CD = HE2_NE2_CD.angleType.angle[HE2_NE2_CD.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(
                determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0), 0, oe1, 0, 3); // OE1
            coordsArray = fillCoordsArray(OE1, coordsArray, oe1);

            arraycopy(
                determineIntxyz(cd, dNE2_CD, cg, dNE2_CD_CG, oe1, 124.0, 1), 0, ne2, 0, 3); // NE2
            coordsArray = fillCoordsArray(NE2, coordsArray, ne2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(ne2, dHE2_NE2, cd, dHE2_NE2_CD, cg, 0.0, 0), 0, he21, 0, 3); // HE21
            coordsArray = fillCoordsArray(HE21, coordsArray, he21);

            arraycopy(
                determineIntxyz(ne2, dHE2_NE2, cd, dHE2_NE2_CD, he21, 120.0, 1),
                0,
                he22,
                0,
                3); // HE22
            coordsArray = fillCoordsArray(HE22, coordsArray, he22);
            break;
          }
        case MET:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom SD = (Atom) newResidue.getAtomNode("SD");
            double[] sd = coordsArray[SD.getIndex() - 1];
            Atom CE = (Atom) newResidue.getAtomNode("CE");
            double[] ce = coordsArray[CE.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
            double[] he1 = coordsArray[HE1.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
            double[] he3 = coordsArray[HE3.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond SD_CG = SD.getBond(CG);
            Bond CE_SD = CE.getBond(SD);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HE_CE = HE1.getBond(CE);
            double dCG_CB = CG_CB.bondType.distance;
            double dSD_CG = SD_CG.bondType.distance;
            double dCE_SD = CE_SD.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle SD_CG_CB = SD.getAngle(CG, CB);
            Angle CE_SD_CG = CE.getAngle(SD, CG);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HE_CE_SD = HE1.getAngle(CE, SD);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dSD_CG_CB = SD_CG_CB.angleType.angle[SD_CG_CB.nh];
            double dCE_SD_CG = CE_SD_CG.angleType.angle[CE_SD_CG.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHE_CE_SD = HE_CE_SD.angleType.angle[HE_CE_SD.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dSD_CG, cb, dSD_CG_CB, ca, 180.0, 0), 0, sd, 0, 3); // SD
            coordsArray = fillCoordsArray(SD, coordsArray, sd);

            arraycopy(determineIntxyz(sd, dCE_SD, cg, dCE_SD_CG, cb, 180.0, 0), 0, ce, 0, 3); // CE
            coordsArray = fillCoordsArray(CE, coordsArray, ce);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, sd, 112.0, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, sd, 112.0, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, cg, 180.0, 0), 0, he1, 0, 3); // HE1
            coordsArray = fillCoordsArray(HE1, coordsArray, he1);

            arraycopy(
                determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, he1, 109.4, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);

            arraycopy(
                determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, he1, 109.4, -1), 0, he3, 0, 3); // HE3
            coordsArray = fillCoordsArray(HE3, coordsArray, he3);
            break;
          }
        case LYS:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom CE = (Atom) newResidue.getAtomNode("CE");
            double[] ce = coordsArray[CE.getIndex() - 1];
            Atom NZ = (Atom) newResidue.getAtomNode("NZ");
            double[] nz = coordsArray[NZ.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
            double[] hd3 = coordsArray[HD3.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
            double[] he3 = coordsArray[HE3.getIndex() - 1];
            Atom HZ1 = (Atom) newResidue.getAtomNode("HZ1");
            double[] hz1 = coordsArray[HZ1.getIndex() - 1];
            Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
            double[] hz2 = coordsArray[HZ2.getIndex() - 1];
            Atom HZ3 = (Atom) newResidue.getAtomNode("HZ3");
            double[] hz3 = coordsArray[HZ3.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond CE_CD = CE.getBond(CD);
            Bond NZ_CE = NZ.getBond(CE);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HD_CD = HD2.getBond(CD);
            Bond HE_CE = HE2.getBond(CE);
            Bond HZ_NZ = HZ1.getBond(NZ);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dCE_CD = CE_CD.bondType.distance;
            double dNZ_CE = NZ_CE.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            double dHZ_NZ = HZ_NZ.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle CE_CD_CG = CE.getAngle(CD, CG);
            Angle NZ_CE_CD = NZ.getAngle(CE, CD);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HD_CD_CG = HD2.getAngle(CD, CG);
            Angle HE_CE_CD = HE2.getAngle(CE, CD);
            Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
            double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
            double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(determineIntxyz(cd, dCE_CD, cg, dCE_CD_CG, cb, 180.0, 0), 0, ce, 0, 3); // CE
            coordsArray = fillCoordsArray(CE, coordsArray, ce);

            arraycopy(determineIntxyz(ce, dNZ_CE, cd, dNZ_CE_CD, cg, 180.0, 0), 0, nz, 0, 3); // NZ
            coordsArray = fillCoordsArray(NZ, coordsArray, nz);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, -1), 0, hd3, 0, 3); // HD3
            coordsArray = fillCoordsArray(HD3, coordsArray, hd3);

            arraycopy(
                determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);

            arraycopy(
                determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, -1), 0, he3, 0, 3); // HE3
            coordsArray = fillCoordsArray(HE3, coordsArray, he3);

            arraycopy(
                determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, cd, 180.0, 0), 0, hz1, 0, 3); // HZ1
            coordsArray = fillCoordsArray(HZ1, coordsArray, hz1);

            arraycopy(
                determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, 1), 0, hz2, 0, 3); // HZ2
            coordsArray = fillCoordsArray(HZ2, coordsArray, hz2);

            arraycopy(
                determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, -1), 0, hz3, 0, 3); // HZ3
            coordsArray = fillCoordsArray(HZ3, coordsArray, hz3);
            break;
          }
        case LYD:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom CE = (Atom) newResidue.getAtomNode("CE");
            double[] ce = coordsArray[CE.getIndex() - 1];
            Atom NZ = (Atom) newResidue.getAtomNode("NZ");
            double[] nz = coordsArray[NZ.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HD2.getIndex() - 1];
            Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
            double[] hd3 = coordsArray[HD3.getIndex() - 1];
            Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
            double[] he2 = coordsArray[HE2.getIndex() - 1];
            Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
            double[] he3 = coordsArray[HE3.getIndex() - 1];
            Atom HZ1 = (Atom) newResidue.getAtomNode("HZ1");
            double[] hz1 = coordsArray[HZ1.getIndex() - 1];
            Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
            double[] hz2 = coordsArray[HZ2.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond CE_CD = CE.getBond(CD);
            Bond NZ_CE = NZ.getBond(CE);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HD_CD = HD2.getBond(CD);
            Bond HE_CE = HE2.getBond(CE);
            Bond HZ_NZ = HZ1.getBond(NZ);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dCE_CD = CE_CD.bondType.distance;
            double dNZ_CE = NZ_CE.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_CE = HE_CE.bondType.distance;
            double dHZ_NZ = HZ_NZ.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle CE_CD_CG = CE.getAngle(CD, CG);
            Angle NZ_CE_CD = NZ.getAngle(CE, CD);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HD_CD_CG = HD2.getAngle(CD, CG);
            Angle HE_CE_CD = HE2.getAngle(CE, CD);
            Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
            double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
            double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(determineIntxyz(cd, dCE_CD, cg, dCE_CD_CG, cb, 180.0, 0), 0, ce, 0, 3); // CE
            coordsArray = fillCoordsArray(CE, coordsArray, ce);

            arraycopy(determineIntxyz(ce, dNZ_CE, cd, dNZ_CE_CD, cg, 180.0, 0), 0, nz, 0, 3); // NZ
            coordsArray = fillCoordsArray(NZ, coordsArray, nz);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, -1), 0, hd3, 0, 3); // HD3
            coordsArray = fillCoordsArray(HD3, coordsArray, hd3);

            arraycopy(
                determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, 1), 0, he2, 0, 3); // HE2
            coordsArray = fillCoordsArray(HE2, coordsArray, he2);

            arraycopy(
                determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, -1), 0, he3, 0, 3); // HE3
            coordsArray = fillCoordsArray(HE3, coordsArray, he3);

            arraycopy(
                determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, cd, 180.0, 0), 0, hz1, 0, 3); // HZ1
            coordsArray = fillCoordsArray(HZ1, coordsArray, hz1);

            arraycopy(
                determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, 1), 0, hz2, 0, 3); // HZ2
            coordsArray = fillCoordsArray(HZ2, coordsArray, hz2);
            break;
          }
        case ARG:
          {
            Atom CB = (Atom) newResidue.getAtomNode("CB");
            double[] cb = coordsArray[CB.getIndex() - 1];
            Atom CG = (Atom) newResidue.getAtomNode("CG");
            double[] cg = coordsArray[CG.getIndex() - 1];
            Atom CD = (Atom) newResidue.getAtomNode("CD");
            double[] cd = coordsArray[CD.getIndex() - 1];
            Atom NE = (Atom) newResidue.getAtomNode("NE");
            double[] ne = coordsArray[NE.getIndex() - 1];
            Atom CZ = (Atom) newResidue.getAtomNode("CZ");
            double[] cz = coordsArray[CZ.getIndex() - 1];
            Atom NH1 = (Atom) newResidue.getAtomNode("NH1");
            double[] nh1 = coordsArray[NH1.getIndex() - 1];
            Atom NH2 = (Atom) newResidue.getAtomNode("NH2");
            double[] nh2 = coordsArray[NH2.getIndex() - 1];
            Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
            double[] hb2 = coordsArray[HB2.getIndex() - 1];
            Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
            double[] hb3 = coordsArray[HB3.getIndex() - 1];
            Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
            double[] hg2 = coordsArray[HG2.getIndex() - 1];
            Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
            double[] hg3 = coordsArray[HG3.getIndex() - 1];
            Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
            double[] hd2 = coordsArray[HG2.getIndex() - 1];
            Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
            double[] hd3 = coordsArray[HD3.getIndex() - 1];
            Atom HE = (Atom) newResidue.getAtomNode("HE");
            double[] he = coordsArray[HE.getIndex() - 1];
            Atom HH11 = (Atom) newResidue.getAtomNode("HH11");
            double[] hh11 = coordsArray[HH11.getIndex() - 1];
            Atom HH12 = (Atom) newResidue.getAtomNode("HH12");
            double[] hh12 = coordsArray[HH12.getIndex() - 1];
            Atom HH21 = (Atom) newResidue.getAtomNode("HH21");
            double[] hh21 = coordsArray[HH21.getIndex() - 1];
            Atom HH22 = (Atom) newResidue.getAtomNode("HH22");
            double[] hh22 = coordsArray[HH22.getIndex() - 1];
            Bond CG_CB = CG.getBond(CB);
            Bond CD_CG = CD.getBond(CG);
            Bond NE_CD = NE.getBond(CD);
            Bond CZ_NE = CZ.getBond(NE);
            Bond NH_CZ = NH1.getBond(CZ);
            Bond HB_CB = HB2.getBond(CB);
            Bond HG_CG = HG2.getBond(CG);
            Bond HD_CD = HD2.getBond(CD);
            Bond HE_NE = HE.getBond(NE);
            Bond HH_NH = HH11.getBond(NH1);
            double dCG_CB = CG_CB.bondType.distance;
            double dCD_CG = CD_CG.bondType.distance;
            double dNE_CD = NE_CD.bondType.distance;
            double dCZ_NE = CZ_NE.bondType.distance;
            double dNH_CZ = NH_CZ.bondType.distance;
            double dHB_CB = HB_CB.bondType.distance;
            double dHG_CG = HG_CG.bondType.distance;
            double dHD_CD = HD_CD.bondType.distance;
            double dHE_NE = HE_NE.bondType.distance;
            double dHH_NH = HH_NH.bondType.distance;
            Angle CG_CB_CA = CG.getAngle(CB, CA);
            Angle CD_CG_CB = CD.getAngle(CG, CB);
            Angle NE_CD_CG = NE.getAngle(CD, CG);
            Angle CZ_NE_CD = CZ.getAngle(NE, CD);
            Angle NH_CZ_NE = NH1.getAngle(CZ, NE);
            Angle HB_CB_CA = HB2.getAngle(CB, CA);
            Angle HG_CG_CB = HG2.getAngle(CG, CB);
            Angle HD_CD_CG = HD2.getAngle(CD, CG);
            Angle HE_NE_CD = HE.getAngle(NE, CD);
            Angle HH_NH_CZ = HH11.getAngle(NH1, CZ);
            double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
            double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
            double dNE_CD_CG = NE_CD_CG.angleType.angle[NE_CD_CG.nh];
            double dCZ_NE_CD = CZ_NE_CD.angleType.angle[CZ_NE_CD.nh];
            double dNH_CZ_NE = NH_CZ_NE.angleType.angle[NH_CZ_NE.nh];
            double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
            double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
            double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
            double dHE_NE_CD = HE_NE_CD.angleType.angle[HE_NE_CD.nh];
            double dHH_NH_CZ = HH_NH_CZ.angleType.angle[HH_NH_CZ.nh];
            arraycopy(determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1), 0, cb, 0, 3);
            coordsArray = fillCoordsArray(CB, coordsArray, cb);

            arraycopy(
                determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, rotScale * 180.0, 0),
                0,
                cg,
                0,
                3); // CG
            coordsArray = fillCoordsArray(CG, coordsArray, cg);

            arraycopy(determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0), 0, cd, 0, 3); // CD
            coordsArray = fillCoordsArray(CD, coordsArray, cd);

            arraycopy(determineIntxyz(cd, dNE_CD, cg, dNE_CD_CG, cb, 180.0, 0), 0, ne, 0, 3); // NE
            coordsArray = fillCoordsArray(NE, coordsArray, ne);

            arraycopy(determineIntxyz(ne, dCZ_NE, cd, dCZ_NE_CD, cg, 180.0, 0), 0, cz, 0, 3); // CZ
            coordsArray = fillCoordsArray(CZ, coordsArray, cz);

            arraycopy(determineIntxyz(cz, dNH_CZ, ne, dNH_CZ_NE, cd, 180, 0), 0, nh1, 0, 3); // NH1
            coordsArray = fillCoordsArray(NH1, coordsArray, nh1);

            arraycopy(
                determineIntxyz(cz, dNH_CZ, ne, dNH_CZ_NE, nh1, 120.0, 1), 0, nh2, 0, 3); // NH2
            coordsArray = fillCoordsArray(NH2, coordsArray, nh2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1), 0, hb2, 0, 3); // HB2
            coordsArray = fillCoordsArray(HB2, coordsArray, hb2);

            arraycopy(
                determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1), 0, hb3, 0, 3); // HB3
            coordsArray = fillCoordsArray(HB3, coordsArray, hb3);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1), 0, hg2, 0, 3); // HG2
            coordsArray = fillCoordsArray(HG2, coordsArray, hg2);

            arraycopy(
                determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1), 0, hg3, 0, 3); // HG3
            coordsArray = fillCoordsArray(HG3, coordsArray, hg3);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ne, 109.4, 1), 0, hd2, 0, 3); // HD2
            coordsArray = fillCoordsArray(HD2, coordsArray, hd2);

            arraycopy(
                determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ne, 109.4, -1), 0, hd3, 0, 3); // HD3
            coordsArray = fillCoordsArray(HD3, coordsArray, hd3);

            arraycopy(determineIntxyz(ne, dHE_NE, cd, dHE_NE_CD, cz, 120.0, 1), 0, he, 0, 3); // HE
            coordsArray = fillCoordsArray(HE, coordsArray, he);

            arraycopy(
                determineIntxyz(nh1, dHH_NH, cz, dHH_NH_CZ, ne, 180.0, 0), 0, hh11, 0, 3); // HH11
            coordsArray = fillCoordsArray(HH11, coordsArray, hh11);

            arraycopy(
                determineIntxyz(nh1, dHH_NH, cz, dHH_NH_CZ, hh11, 120.0, 1), 0, hh12, 0, 3); // HH12
            coordsArray = fillCoordsArray(HH12, coordsArray, hh12);

            arraycopy(
                determineIntxyz(nh2, dHH_NH, cz, dHH_NH_CZ, ne, 180.0, 0), 0, hh21, 0, 3); // HH21
            coordsArray = fillCoordsArray(HH21, coordsArray, hh21);

            arraycopy(
                determineIntxyz(nh2, dHH_NH, cz, dHH_NH_CZ, hh21, 120.0, 1), 0, hh22, 0, 3); // HH22
            coordsArray = fillCoordsArray(HH22, coordsArray, hh22);
            break;
          }
        default:
          break;
      }
    }
    double[] coordsArray1D = new double[atomArray.length * 3];
    for (int i = 0; i < atomArray.length; i++) {
      for (int j = 0; j < 3; j++) {
        coordsArray1D[i * 3 + j] = coordsArray[i][j];
      }
    }
    return coordsArray1D;
  }

  private void setAltCoordinates(double[] coordinates) {
    for (int i = 0; i < coordinates.length / 3; i++) {
      for (int j = 0; j < 3; j++) {
        this.altCoords[i][j] = coordinates[i * 3 + j];
      }
    }
    useAltCoordinates(true);
  }

  private double[][] fillCoordsArray(Atom atom, double[][] coordsArray, double[] determinedXYZ) {
    int XYZIndex = atom.getIndex() - 1;
    coordsArray[XYZIndex][0] = determinedXYZ[0];
    coordsArray[XYZIndex][1] = determinedXYZ[1];
    coordsArray[XYZIndex][2] = determinedXYZ[2];
    return coordsArray;
  }
}
