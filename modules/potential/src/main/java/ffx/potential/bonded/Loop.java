/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.bonded;

import static ffx.numerics.VectorMath.diff;
import java.io.File;
import java.util.ArrayList;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.MolecularAssembly;

import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;

import java.util.List;

/**
 * @author Mallory R. Tollefson
 */
public class Loop {

    private static final Logger logger = Logger.getLogger(Loop.class.getName());
    private final MolecularAssembly molAss;
    private final boolean writeFile = false;
    
    int max_soln = 16;
    double[][] r_n = new double[5][3];
    double[][] r_a = new double[5][3];
    double[][] r_c = new double[5][3];
    double[][] xyz_o = new double[5][3];
    public final LoopClosure loopClosure;
    public final SturmMethod sturmMethod;

    public Loop(MolecularAssembly molAss, int firstResidue, int endResidue, boolean writeFile) {

        loopClosure = new LoopClosure();
        sturmMethod = new SturmMethod();
        this.molAss = molAss;
        generateLoops(firstResidue, endResidue);
        
    }

    public Loop(MolecularAssembly molAss) {
        loopClosure = new LoopClosure();
        sturmMethod = new SturmMethod();
        this.molAss = molAss;
    }

    public List<double[]> generateLoops(int firstResidue, int endResidue) {
        ArrayList<Atom> backBoneAtoms = molAss.getBackBoneAtoms();
        
        boolean bool1 = true;
        int i = 0;
        List<double[]> solutions = new ArrayList<>();
        
        while (bool1) {
            Atom atom = backBoneAtoms.get(i);
            int resID = atom.getResidueNumber();

            if (resID > endResidue) {
                //terminates the collection of atom coordinates
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
                //coordinateArray temporarily holds the coordinates of a
                //      specific atom from the pdb file
                double[] coordinateArray = atom.getXYZ(null);

                if (atmname.contentEquals(n)) {
                    //Backbone nitrogen coordinates are stored in r_n[]
                    r_n[ir][0] = coordinateArray[0];
                    r_n[ir][1] = coordinateArray[1];
                    r_n[ir][2] = coordinateArray[2];
                } else if (atmname.contentEquals(ca)) {
                    //Backbone alpha carbon coordinates are stored in r_a[]
                    r_a[ir][0] = coordinateArray[0];
                    r_a[ir][1] = coordinateArray[1];
                    r_a[ir][2] = coordinateArray[2];
                } else if (atmname.contentEquals(c)) {
                    //Backbone carbon coordinates are stored in r_c[]
                    r_c[ir][0] = coordinateArray[0];
                    r_c[ir][1] = coordinateArray[1];
                    r_c[ir][2] = coordinateArray[2];
                }
                i++;
            }
        }

        /**
         * Method that solves the 16th degree, 3 peptide polynomial.
         */
        double[][][] r_soln_n = new double[max_soln][3][3];
        double[][][] r_soln_a = new double[max_soln][3][3];
        double[][][] r_soln_c = new double[max_soln][3][3];
        int[] n_soln = new int[1];
         
        loopClosure.solve3PepPoly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, n_soln);
        
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Starting res.:             %d\n", firstResidue));
        sb.append(String.format(" Ending res.:               %d\n", endResidue));
        sb.append(String.format(" No. of solutions:          %d\n", n_soln[0]));
        logger.info(sb.toString());

        for (int k = 0; k < n_soln[0]; k++) {
            double [] coordsArray;
            coordsArray = getSolutionCoordinates(k, r_soln_n, r_soln_a, r_soln_c, firstResidue, endResidue);
            solutions.add(coordsArray);
        }

        return solutions;
    }
    
    // Only for JUnit testing.
    public double[][] getRN() {
        return r_n;
    }

    public double[][] getRA() {
        return r_a;
    }

    public double[][] getRC() {
        return r_c;
    }

    private double[] getSolutionCoordinates(int k, double[][][] r_soln_n, double[][][] r_soln_a, double[][][] r_soln_c, int stt_res, int end_res) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                r_n[i + 1][j] = r_soln_n[k][i][j];
                r_a[i + 1][j] = r_soln_a[k][i][j];
                r_c[i + 1][j] = r_soln_c[k][i][j];
            }
        }
       // double sum = 0.0;

        /* For RMSD?
         for (int i = 0; i < 3; i++) {
         for (int j = 0; j < 3; j++) {
         dr[j] = r_soln_n[k][i][j] - r_n[i + 1][j];
         }
         sum += dot(dr, dr);
         for (int j = 0; j < 3; j++) {
         dr[j] = r_soln_a[k][i][j] - r_a[i + 1][j];
         }
         sum += dot(dr, dr);
         for (int j = 0; j < 3; j++) {
         dr[j] = r_soln_c[k][i][j] - r_c[i + 1][j];
         }
         sum += dot(dr, dr);
         }
         */

        Polymer[] newChain = molAss.getChains();
        ArrayList<Atom> backBoneAtoms;

        Atom[] atomArray = molAss.getAtomArray();
        int index = 0;
        double[] coordsArray = new double[atomArray.length * 3];
        for (int i = 0; i < atomArray.length; i++) {
            Atom a = atomArray[i];
            coordsArray[index++] = a.getX();
            coordsArray[index++] = a.getY();
            coordsArray[index++] = a.getZ();
        }
        
        for (int i = stt_res + 1; i < end_res; i++) {
            Residue newResidue = newChain[0].getResidue(i);
            Residue backResidue = newChain[0].getResidue(i-1);
            backBoneAtoms = newResidue.getBackboneAtoms();
            Atom C = null;
            Atom backC = null;
            Atom N = null;
            Atom CA = null;
            Atom O = null;
            double[] c = new double[3];
            double[] ca = new double[3];
            double[] n = new double[3];
            double[] bc = new double[3];
            double[] determinedXYZ = new double[3];
            
            for (Atom backBoneAtom : backBoneAtoms) {
                backBoneAtom.setBuilt(true);
                int backBoneIndex = 0;
                switch (backBoneAtom.getAtomType().name) {
                    case "C":
                        backC = C;
                        C = backBoneAtom;
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        coordsArray[backBoneIndex] = r_c[i - stt_res][0];
                        coordsArray[backBoneIndex + 1] = r_c[i - stt_res][1];
                        coordsArray[backBoneIndex + 2] = r_c[i - stt_res][2];
                        break;
                    case "N":
                        N = backBoneAtom;
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        coordsArray[backBoneIndex] = r_n[i - stt_res][0];
                        coordsArray[backBoneIndex + 1] = r_n[i - stt_res][1];
                        coordsArray[backBoneIndex + 2] = r_n[i - stt_res][2];
                        break;
                    case "CA":
                        CA = backBoneAtom;
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        coordsArray[backBoneIndex] = r_a[i - stt_res][0];
                        coordsArray[backBoneIndex + 1] = r_a[i - stt_res][1];
                        coordsArray[backBoneIndex + 2] = r_a[i - stt_res][2];
                        break;
                    case "H":
                        ((Atom) backResidue.getAtomNode("C")).getXYZ(bc);
                        ((Atom) newResidue.getAtomNode("CA")).getXYZ(ca);
                        ((Atom) newResidue.getAtomNode("N")).getXYZ(n);
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        determinedXYZ = BondedUtils.determineIntxyz(n, 1.0, bc, 109.5, ca, 109.5, 0);
                        coordsArray[backBoneIndex] =     determinedXYZ[0];
                        coordsArray[backBoneIndex + 1] = determinedXYZ[1];
                        coordsArray[backBoneIndex + 2] = determinedXYZ[2];
                        break;
                    case "HA":
                        ((Atom) newResidue.getAtomNode("C")).getXYZ(c);
                        ((Atom) newResidue.getAtomNode("CA")).getXYZ(ca);
                        ((Atom) newResidue.getAtomNode("N")).getXYZ(n);
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        determinedXYZ = BondedUtils.determineIntxyz(ca, 1.0, n, 109.5, c, 109.5, 0);
                        coordsArray[backBoneIndex] =     determinedXYZ[0];
                        coordsArray[backBoneIndex + 1] = determinedXYZ[1];
                        coordsArray[backBoneIndex + 2] = determinedXYZ[2];
                        break;
                    case "O":
                        ((Atom) newResidue.getAtomNode("C")).getXYZ(c);
                        ((Atom) newResidue.getAtomNode("CA")).getXYZ(ca);
                        ((Atom) newResidue.getAtomNode("N")).getXYZ(n);
                        backBoneIndex = backBoneAtom.getXYZIndex();
                        determinedXYZ = BondedUtils.determineIntxyz(c, 1.2255, ca, 122.4, n, 180, 0);
                        coordsArray[backBoneIndex] =     determinedXYZ[0];
                        coordsArray[backBoneIndex + 1] = determinedXYZ[1];
                        coordsArray[backBoneIndex + 2] = determinedXYZ[2];
                        logger.info(String.format("coords O. X: %f Y: %f Z: %f",coordsArray[backBoneIndex],coordsArray[backBoneIndex + 1],coordsArray[backBoneIndex + 2]));
                        break;
                }
            }
        
            //Obtaining coordinates for sidechains
            AminoAcid3 name = AminoAcid3.valueOf(newResidue.getName());
            ((Atom) newResidue.getAtomNode("C")).getXYZ(c);
            ((Atom) newResidue.getAtomNode("CA")).getXYZ(ca);
            ((Atom) newResidue.getAtomNode("N")).getXYZ(n);
            switch (name) {
                    case VAL: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG1 = (Atom) newResidue.getAtomNode("CG1");
                    double [] cg1 = CG1.getXYZ(null);
                    Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
                    double [] cg2 = CG2.getXYZ(null);
                    Atom HB = (Atom) newResidue.getAtomNode("HB");
                    double [] hb = HB.getXYZ(null);
                    Atom HG11 = (Atom) newResidue.getAtomNode("HG11");
                    double [] hg11 = HG11.getXYZ(null);
                    Atom HG12 = (Atom) newResidue.getAtomNode("HG12");
                    double [] hg12 = HG12.getXYZ(null);
                    Atom HG13 = (Atom) newResidue.getAtomNode("HG13");
                    double [] hg13 = HG13.getXYZ(null);
                    Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
                    double [] hg21 = HG21.getXYZ(null);
                    Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
                    double [] hg22 = HG22.getXYZ(null);
                    Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
                    double [] hg23 = HG23.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);                    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);         //CG1    
                    coordsArray = fillCoordsArray(CG1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, cg1, 109.5, -1);      //CG2
                    coordsArray = fillCoordsArray(CG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg1, 109.4, 1);        //HB
                    coordsArray = fillCoordsArray(HB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, ca, 180.0, 0);      //HG11
                    coordsArray = fillCoordsArray(HG11,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, hg11, 109.4, 1);    //HG12
                    coordsArray = fillCoordsArray(HG12,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dHG_CG, cb, dHG_CG_CB, hg11, 109.4, -1);   //HG13
                    coordsArray = fillCoordsArray(HG13,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, ca, 180.0, 0);      //HG21
                    coordsArray = fillCoordsArray(HG21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, hg21, 109.4, 1);    //HG22
                    coordsArray = fillCoordsArray(HG22,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG_CG, cb, dHG_CG_CB, hg21, 109.4, -1);   //HG23
                    coordsArray = fillCoordsArray(HG23,coordsArray, determinedXYZ);
                    break;
                }
                case LEU: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG = (Atom) newResidue.getAtomNode("HG");
                    double [] hg = HG.getXYZ(null);
                    Atom HD11 = (Atom) newResidue.getAtomNode("HD11");
                    double [] hd11 = HD11.getXYZ(null);
                    Atom HD12 = (Atom) newResidue.getAtomNode("HD12");
                    double [] hd12 = HD12.getXYZ(null);
                    Atom HD13 = (Atom) newResidue.getAtomNode("HD13");
                    double [] hd13 = HD13.getXYZ(null);
                    Atom HD21 = (Atom) newResidue.getAtomNode("HD21");
                    double [] hd21 = HD21.getXYZ(null);
                    Atom HD22 = (Atom) newResidue.getAtomNode("HD22");
                    double [] hd22 = HD22.getXYZ(null);
                    Atom HD23 = (Atom) newResidue.getAtomNode("HD23");
                    double [] hd23 = HD23.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);     
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);          //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 109.5, -1);        //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd1, 109.4, 1);          //HG
                    coordsArray = fillCoordsArray(HG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, cb, 180.0, 0);        //HD11
                    coordsArray = fillCoordsArray(HD11,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, hd11, 109.4, 1);      //HD12
                    coordsArray = fillCoordsArray(HD12,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, hd11, 109.4, -1);     //HD13
                    coordsArray = fillCoordsArray(HD13,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, cb, 180.0, 0);        //HD21
                    coordsArray = fillCoordsArray(HD21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, hd21, 109.4, 1);      //HD22
                    coordsArray = fillCoordsArray(HD22,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, hd21, 109.4, -1);     //HD23
                    coordsArray = fillCoordsArray(HD23,coordsArray, determinedXYZ);
                    break;
                }
                case ILE: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG1 = (Atom) newResidue.getAtomNode("CG1");
                    double [] cg1 = CG1.getXYZ(null);
                    Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
                    double [] cg2 = CG2.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom HB = (Atom) newResidue.getAtomNode("HB");
                    double [] hb = HB.getXYZ(null);
                    Atom HG12 = (Atom) newResidue.getAtomNode("HG12");
                    double [] hg12 = HG12.getXYZ(null);
                    Atom HG13 = (Atom) newResidue.getAtomNode("HG13");
                    double [] hg13 = HG13.getXYZ(null);
                    Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
                    double [] hg21 = HG21.getXYZ(null);
                    Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
                    double [] hg22 = HG22.getXYZ(null);
                    Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
                    double [] hg23 = HG23.getXYZ(null);
                    Atom HD11 = (Atom) newResidue.getAtomNode("HD11");
                    double [] hd11 = HD11.getXYZ(null);
                    Atom HD12 = (Atom) newResidue.getAtomNode("HD12");
                    double [] hd12 = HD12.getXYZ(null);
                    Atom HD13 = (Atom) newResidue.getAtomNode("HD13");
                    double [] hd13 = HD13.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG1_CB, ca, dCG1_CB_CA, n, 0.0, 0);               //CG1
                    coordsArray = fillCoordsArray(CG1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG2_CB, ca, dCG2_CB_CA, cg1, 109.5, 1);           //CG2
                    coordsArray = fillCoordsArray(CG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dCD1_CG1, cb, dCD1_CG1_CB, ca, 180, 0);           //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg2, 109.4, 1);              //HB
                    coordsArray = fillCoordsArray(HB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dHG1_CG, cb, dHG1_CG_CB, cd1, 109.4, 1);         //HG12
                    coordsArray = fillCoordsArray(HG12,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg1, dHG1_CG, cb, dHG1_CG_CB, cd1, 109.4, -1);        //HG13
                    coordsArray = fillCoordsArray(HG13,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, cg1, 180.0, 0);         //HG21
                    coordsArray = fillCoordsArray(HG21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, hg21, 109.0, 1);        //HG22
                    coordsArray = fillCoordsArray(HG22,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG, cb, dHG2_CG_CB, hg21, 109.0, -1);       //HG23
                    coordsArray = fillCoordsArray(HG23,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, cb, 180.0, 0);         //HD11
                    coordsArray = fillCoordsArray(HD11,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, hd11, 109.0, 1);       //HD12
                    coordsArray = fillCoordsArray(HD12,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg1, dHD_CD1_CG1, hd11, 109.0, -1);      //HD13
                    coordsArray = fillCoordsArray(HD13,coordsArray, determinedXYZ);
                    break;
                }
                case SER: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom OG = (Atom) newResidue.getAtomNode("OG");
                    double [] og = OG.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG = (Atom) newResidue.getAtomNode("HG");
                    double [] hg = HG.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dOG_CB, ca, dOG_CB_CA, n, 180.0, 0);            //OG
                    coordsArray = fillCoordsArray(OG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og, 106.7, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og, 106.7, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(og, dHG_OG, cb, dHG_OG_CB, ca, 180.0, 0);           //HG
                    coordsArray = fillCoordsArray(HG,coordsArray, determinedXYZ);
                    break;
                }
                case THR: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom OG1 = (Atom) newResidue.getAtomNode("OG1");
                    double [] og1 = OG1.getXYZ(null);
                    Atom CG2 = (Atom) newResidue.getAtomNode("CG2");
                    double [] cg2 = CG2.getXYZ(null);
                    Atom HB = (Atom) newResidue.getAtomNode("HB");
                    double [] hb = HB.getXYZ(null);
                    Atom HG1 = (Atom) newResidue.getAtomNode("HG1");
                    double [] hg1 = HG1.getXYZ(null);
                    Atom HG21 = (Atom) newResidue.getAtomNode("HG21");
                    double [] hg21 = HG21.getXYZ(null);
                    Atom HG22 = (Atom) newResidue.getAtomNode("HG22");
                    double [] hg22 = HG22.getXYZ(null);
                    Atom HG23 = (Atom) newResidue.getAtomNode("HG23");
                    double [] hg23 = HG23.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dOG1_CB, ca, dOG1_CB_CA, n, 180.0, 0);                 //OG1
                    coordsArray = fillCoordsArray(OG1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG2_CB, ca, dCG2_CB_CA, og1, 107.7, 1);               //CG2
                    coordsArray = fillCoordsArray(CG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, og1, 106.7, -1);                 //HB
                    coordsArray = fillCoordsArray(HB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(og1, dHG1_OG1, cb, dHG1_OG1_CB, ca, 180.0, 0);             //HG1
                    coordsArray = fillCoordsArray(HG1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, ca, 180.0, 0);            //HG21
                    coordsArray = fillCoordsArray(HG21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, hg21, 109.0, 1);          //HG22
                    coordsArray = fillCoordsArray(HG22,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg2, dHG2_CG2, cb, dHG2_CG2_CB, hg21, 109.0, -1);         //HG23
                    coordsArray = fillCoordsArray(HG23,coordsArray, determinedXYZ);
                    break;
                }
                case CYS:
                case CYX: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom SG = (Atom) newResidue.getAtomNode("SG");
                    double [] sg = SG.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG = (Atom) newResidue.getAtomNode("HG");
                    double [] hg = HG.getXYZ(null);
                    if (CA == null || CB == null || N == null || SG == null || HB2 == null || HB3 == null || HG == null) {
                        break;
                    }
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dSG_CB, ca, dSG_CB_CA, n, 180.0, 0);             //SG
                    coordsArray = fillCoordsArray(SG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, 1);           //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, -1);          //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(sg, dHG_SG, cb, dHG_SG_CB, ca, 180.0, 0);            //HG
                    coordsArray = fillCoordsArray(HG,coordsArray, determinedXYZ);
                    break;
                }
                case CYD: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom SG = (Atom) newResidue.getAtomNode("SG");
                    double [] sg = SG.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Bond SG_CB = SG.getBond(CB);
                    Bond HB_CB = HB2.getBond(CB);
                    double dSG_CB = SG_CB.bondType.distance;
                    double dHB_CB = HB_CB.bondType.distance;
                    Angle SG_CB_CA = SG.getAngle(CB, CA);
                    Angle HB_CB_CA = HB2.getAngle(CB, CA);
                    double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
                    double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dSG_CB, ca, dSG_CB_CA, n, 180.0, 0);            //SG
                    coordsArray = fillCoordsArray(SG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, sg, 112.0, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    break;
                }
                case PHE: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
                    double [] ce2 = CE2.getXYZ(null);
                    Atom CZ = (Atom) newResidue.getAtomNode("CZ");
                    double [] cz = CZ.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
                    Atom HZ = (Atom) newResidue.getAtomNode("HZ");
                    double [] hz = HZ.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);          //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1);         //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);           //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);           //CE2
                    coordsArray = fillCoordsArray(CE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0);        //CZ
                    coordsArray = fillCoordsArray(CZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD1_CG, ce1, 120.0, 1);       //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD1_CG, ce2, 120.0, 1);       //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1);        //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1);        //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz, dHZ_CZ, ce1, dHZ_CZ_CE1, ce2, 120.0, 1);        //HZ
                    coordsArray = fillCoordsArray(HZ,coordsArray, determinedXYZ);
                    break;
                }
                case PRO: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
                    double [] hd3 = HD3.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 30.0, 0);             //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 30.0, 0);            //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, n, 109.4, 1);           //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, n, 109.4, -1);          //HD3
                    coordsArray = fillCoordsArray(HD3,coordsArray, determinedXYZ);
                    break;
                }
                case TYR: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
                    double [] ce2 = CE2.getXYZ(null);
                    Atom CZ = (Atom) newResidue.getAtomNode("CZ");
                    double [] cz = CZ.getXYZ(null);
                    Atom OH = (Atom) newResidue.getAtomNode("OH");
                    double [] oh = OH.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
                    Atom HH = (Atom) newResidue.getAtomNode("HH");
                    double [] hh = HH.getXYZ(null);
                    Bond CG_CB = CG.getBond(CB);
                    Bond CD_CG = CD1.getBond(CG);
                    Bond CE_CD = CE1.getBond(CD1);
                    Bond CZ_CE1 = CZ.getBond(CE1);
                    Bond OH_CZ = OH.getBond(CZ);
                    Bond HB_CB = HB2.getBond(CB);
                    Bond HD_CD = HD1.getBond(CD1);
                    Bond HE_CE = HE1.getBond(CE1);
                    Bond HH_OH = HH.getBond(OH);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);    
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 62.0, 0);                     //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 90.0, 0);                   //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1);                 //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);                   //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);                   //CE2
                    coordsArray = fillCoordsArray(CE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0);                //CZ
                    coordsArray = fillCoordsArray(CZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz, dOH_CZ, ce2, dOH_CZ_CE2, ce1, 120.0, 1);                //OH
                    coordsArray = fillCoordsArray(OH,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);                  //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);                 //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, ce1, 120.0, 1);                //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, ce2, 120.0, 1);                //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1);                //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1);                //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(oh, dHH_OH, cz, dHH_OH_CZ, ce2, 0.0, 0);                    //HH
                    coordsArray = fillCoordsArray(HH,coordsArray, determinedXYZ);
                    break;
                }
                case TYD: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
                    double [] ce2 = CE2.getXYZ(null);
                    Atom CZ = (Atom) newResidue.getAtomNode("CZ");
                    double [] cz = CZ.getXYZ(null);
                    Atom OH = (Atom) newResidue.getAtomNode("OH");
                    double [] oh = OH.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 62.0, 0);             //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 90.0, 0);           //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, cd1, 120.0, 1);         //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);           //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dCE_CD, cg, dCE_CD_CG, cb, 180, 0);           //CE2
                    coordsArray = fillCoordsArray(CE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dCZ_CE1, cd1, dCZ_CE1_CD1, cg, 0.0, 0);        //CZ
                    coordsArray = fillCoordsArray(CZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz, dOH_CZ, ce2, dOH_CZ_CE2, ce1, 120.0, 1);        //OH
                    coordsArray = fillCoordsArray(OH,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD_CD, cg, dHD_CD_CG, ce1, 120.0, 1);        //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD_CD, cg, dHD_CD_CG, ce2, 120.0, 1);        //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE_CE, cd1, dHE_CE_CD, cz, 120.0, 1);        //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce2, dHE_CE, cd2, dHE_CE_CD, cz, 120.0, 1);        //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    break;
                }
                case TRP: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD1 = (Atom) newResidue.getAtomNode("CD1");
                    double [] cd1 = CD1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom NE1 = (Atom) newResidue.getAtomNode("NE1");
                    double [] ne1 = NE1.getXYZ(null);
                    Atom CE2 = (Atom) newResidue.getAtomNode("CE2");
                    double [] ce2 = CE2.getXYZ(null);
                    Atom CE3 = (Atom) newResidue.getAtomNode("CE3");
                    double [] ce3 = CE3.getXYZ(null);
                    Atom CZ2 = (Atom) newResidue.getAtomNode("CZ2");
                    double [] cz2 = CZ2.getXYZ(null);
                    Atom CZ3 = (Atom) newResidue.getAtomNode("CZ3");
                    double [] cz3 = CZ3.getXYZ(null);
                    Atom CH2 = (Atom) newResidue.getAtomNode("CH2");
                    double [] ch2 = CH2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
                    double [] he3 = HE3.getXYZ(null);
                    Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
                    double [] hz2 = HZ2.getXYZ(null);
                    Atom HZ3 = (Atom) newResidue.getAtomNode("HZ3");
                    double [] hz3 = HZ3.getXYZ(null);
                    Atom HH2 = (Atom) newResidue.getAtomNode("HH2");
                    double [] hh2 = HH2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 62.0, 0);                 //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD1_CG, cb, dCD1_CG_CB, ca, -90.0, 0);            //CD1
                    coordsArray = fillCoordsArray(CD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, cd1, 108.0, 1);           //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dNE1_CD1, cg, dNE1_CD1_CG, cd2, 0.0, 0);          //NE1
                    coordsArray = fillCoordsArray(NE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne1, dCE2_NE1, cd1, dCE2_NE1_CD1, cg, 0.0, 0);         //CE2
                    coordsArray = fillCoordsArray(CE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dCE3_CD2, ce2, dCE3_CD2_CE2, ne1, 180.0, 0);      //CE3
                    coordsArray = fillCoordsArray(CE3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce2, dCZ2_CE2, cd2, dCZ2_CE2_CD2, ce3, 0.0, 0);        //CZ2
                    coordsArray = fillCoordsArray(CZ2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce3, dCZ3_CE3, cd2, dCZ3_CE3_CD2, ce2, 0.0, 0);        //CZ3
                    coordsArray = fillCoordsArray(CZ3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz2, dCH2_CZ2, ce2, dCH2_CZ2_CE2, cd2, 0.0, 0);        //CH2
                    coordsArray = fillCoordsArray(CH2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);              //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);             //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd1, dHD1_CD1, cg, dHD1_CD1_CG, ne1, 126.0, 1);        //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne1, dHE1_NE1, cd1, dHE1_NE1_CD1, ce2, 126.0, 1);      //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce3, dHE3_CE3, cd2, dHE3_CE3_CD2, cz3, 120.0, 1);      //HE3
                    coordsArray = fillCoordsArray(HE3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz2, dHZ2_CZ2, ce2, dHZ2_CZ2_CE2, ch2, 120.0, 1);      //HZ2
                    coordsArray = fillCoordsArray(HZ2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz3, dHZ3_CZ3, ce3, dHZ3_CZ3_CE3, ch2, 120.0, 1);      //HZ3
                    coordsArray = fillCoordsArray(HZ3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ch2, dHH2_CH2, cz2, dHH2_CH2_CZ2, cz3, 120.0, 1);      //HH2
                    coordsArray = fillCoordsArray(HH2,coordsArray, determinedXYZ);
                    break;
                }
                case HIS: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
                    double [] nd1 = ND1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
                    double [] ne2 = NE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0);        //ND1
                    coordsArray = fillCoordsArray(ND1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1);       //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0);      //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0);      //NE2
                    coordsArray = fillCoordsArray(NE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd1, dHD1_ND1, cg, dHD1_ND1_CG, cb, 0.0, 0);       //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1);    //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1);  //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne2, dHE2_NE2, cd2, dHE2_NE2_CD2, ce1, 126.0, 1);  //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    break;
                }
                case HID: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
                    double [] nd1 = ND1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
                    double [] ne2 = NE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD1 = (Atom) newResidue.getAtomNode("HD1");
                    double [] hd1 = HD1.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0);        //ND1
                    coordsArray = fillCoordsArray(ND1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1);       //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0);      //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0);      //NE2
                    coordsArray = fillCoordsArray(NE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd1, dHD1_ND1, cg, dHD1_ND1_CG, cb, 0.0, 0);       //HD1
                    coordsArray = fillCoordsArray(HD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1);    //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1);  //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    break;
                }
                case HIE: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom ND1 = (Atom) newResidue.getAtomNode("ND1");
                    double [] nd1 = ND1.getXYZ(null);
                    Atom CD2 = (Atom) newResidue.getAtomNode("CD2");
                    double [] cd2 = CD2.getXYZ(null);
                    Atom CE1 = (Atom) newResidue.getAtomNode("CE1");
                    double [] ce1 = CE1.getXYZ(null);
                    Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
                    double [] ne2 = NE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dND1_CG, cb, dND1_CG_CB, ca, 180.0, 0);        //ND1
                    coordsArray = fillCoordsArray(ND1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD2_CG, cb, dCD2_CG_CB, nd1, 108.0, 1);       //CD2
                    coordsArray = fillCoordsArray(CD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd1, dCE1_ND1, cg, dCE1_ND1_CG, cd2, 0.0, 0);      //CE1
                    coordsArray = fillCoordsArray(CE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dNE2_CD2, cg, dNE2_CD2_CG, nd1, 0.0, 0);      //NE2
                    coordsArray = fillCoordsArray(NE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd2, dHD2_CD2, cg, dHD2_CD2_CG, ne2, 126.0, 1);    //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce1, dHE1_CE1, nd1, dHE1_CE1_ND1, ne2, 126.0, 1);  //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne2, dHE2_NE2, cd2, dHE2_NE2_CD2, ce1, 126.0, 1);  //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    break;
                }
                case ASP: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
                    double [] od1 = OD1.getXYZ(null);
                    Atom OD2 = (Atom) newResidue.getAtomNode("OD2");
                    double [] od2 = OD2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);        //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0);      //OD1
                    coordsArray = fillCoordsArray(OD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dOD2_CG, cb, dOD2_CG_CB, od1, 126.0, 1);   //OD2
                    coordsArray = fillCoordsArray(OD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1);      //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1);     //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    break;
                }
                case ASH: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
                    double [] od1 = OD1.getXYZ(null);
                    Atom OD2 = (Atom) newResidue.getAtomNode("OD2");
                    double [] od2 = OD2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0);          //OD1
                    coordsArray = fillCoordsArray(OD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dOD2_CG, cb, dOD2_CG_CB, od1, 126.0, 1);       //OD2
                    coordsArray = fillCoordsArray(OD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(od2, dHD2_OD2, cg, dHD2_OD2_CG, od1, 0.0, 0);      //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    break;
                }
                case ASN: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom OD1 = (Atom) newResidue.getAtomNode("OD1");
                    double [] od1 = OD1.getXYZ(null);
                    Atom ND2 = (Atom) newResidue.getAtomNode("ND2");
                    double [] nd2 = ND2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HD21 = (Atom) newResidue.getAtomNode("HD21");
                    double [] hd21 = HD21.getXYZ(null);
                    Atom HD22 = (Atom) newResidue.getAtomNode("HD22");
                    double [] hd22 = HD22.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dOD1_CG, cb, dOD1_CG_CB, ca, 0.0, 0);          //OD1
                    coordsArray = fillCoordsArray(OD1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dND2_CG, cb, dND2_CG_CB, od1, 124.0, 1);       //ND2
                    coordsArray = fillCoordsArray(ND2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 107.9, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd2, dHD2_ND2, cg, dHD2_ND2_CG, cb, 0.0, 0);      //HD21
                    coordsArray = fillCoordsArray(HD21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nd2, dHD2_ND2, cg, dHD2_ND2_CG, hd21, 120.0, 1);  //HD22
                    coordsArray = fillCoordsArray(HD22,coordsArray, determinedXYZ);
                    break;
                }
                case GLU: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
                    double [] oe1 = OE1.getXYZ(null);
                    Atom OE2 = (Atom) newResidue.getAtomNode("OE2");
                    double [] oe2 = OE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0);        //OE1
                    coordsArray = fillCoordsArray(OE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dOE2_CD, cg, dOE2_CD_CG, oe1, 126.0, 1);       //OE2
                    coordsArray = fillCoordsArray(OE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    break;
                }
                case GLH: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
                    double [] oe1 = OE1.getXYZ(null);
                    Atom OE2 = (Atom) newResidue.getAtomNode("OE2");
                    double [] oe2 = OE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0);        //OE1
                    coordsArray = fillCoordsArray(OE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dOE2_CD, cg, dOE2_CD_CG, oe1, 126.0, 1);       //OE2
                    coordsArray = fillCoordsArray(OE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(oe2, dHE2_OE2, cd, dHE2_OE2_CD, oe1, 0.0, 0);      //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    break;
                }
                case GLN: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom OE1 = (Atom) newResidue.getAtomNode("OE1");
                    double [] oe1 = OE1.getXYZ(null);
                    Atom NE2 = (Atom) newResidue.getAtomNode("NE2");
                    double [] ne2 = NE2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HE21 = (Atom) newResidue.getAtomNode("HE21");
                    double [] he21 = HE21.getXYZ(null);
                    Atom HE22 = (Atom) newResidue.getAtomNode("HE22");
                    double [] he22 = HE22.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dOE1_CD, cg, dOE1_CD_CG, cb, 180.0, 0);        //OE1
                    coordsArray = fillCoordsArray(OE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dNE2_CD, cg, dNE2_CD_CG, oe1, 124.0, 1);       //NE2
                    coordsArray = fillCoordsArray(NE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 107.9, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne2, dHE2_NE2, cd, dHE2_NE2_CD, cg, 0.0, 0);      //HE21
                    coordsArray = fillCoordsArray(HE21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne2, dHE2_NE2, cd, dHE2_NE2_CD, he21, 120.0, 1);  //HE22
                    coordsArray = fillCoordsArray(HE22,coordsArray, determinedXYZ);
                    break;
                }
                case MET: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom SD = (Atom) newResidue.getAtomNode("SD");
                    double [] sd = SD.getXYZ(null);
                    Atom CE = (Atom) newResidue.getAtomNode("CE");
                    double [] ce = CE.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HE1 = (Atom) newResidue.getAtomNode("HE1");
                    double [] he1 = HE1.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
                    Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
                    double [] he3 = HE3.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dSD_CG, cb, dSD_CG_CB, ca, 180.0, 0);           //SD
                    coordsArray = fillCoordsArray(SD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(sd, dCE_SD, cg, dCE_SD_CG, cb, 180.0, 0);           //CE
                    coordsArray = fillCoordsArray(CE,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, sd, 112.0, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, sd, 112.0, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, cg, 180.0, 0);          //HE1
                    coordsArray = fillCoordsArray(HE1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, he1, 109.4, 1);         //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, sd, dHE_CE_SD, he1, 109.4, -1);        //HE3
                    coordsArray = fillCoordsArray(HE3,coordsArray, determinedXYZ);
                    break;
                }
                case LYS: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom CE = (Atom) newResidue.getAtomNode("CE");
                    double [] ce = CE.getXYZ(null);
                    Atom NZ = (Atom) newResidue.getAtomNode("NZ");
                    double [] nz = NZ.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
                    double [] hd3 = HD3.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
                    Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
                    double [] he3 = HE3.getXYZ(null);
                    Atom HZ1 = (Atom) newResidue.getAtomNode("HZ1");
                    double [] hz1 = HZ1.getXYZ(null);
                    Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
                    double [] hz2 = HZ2.getXYZ(null);
                    Atom HZ3 = (Atom) newResidue.getAtomNode("HZ3");
                    double [] hz3 = HZ3.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dCE_CD, cg, dCE_CD_CG, cb, 180.0, 0);           //CE
                    coordsArray = fillCoordsArray(CE,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dNZ_CE, cd, dNZ_CE_CD, cg, 180.0, 0);           //NZ
                    coordsArray = fillCoordsArray(NZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, 1);          //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, -1);         //HD3
                    coordsArray = fillCoordsArray(HD3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, 1);          //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, -1);         //HE3
                    coordsArray = fillCoordsArray(HE3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, cd, 180.0, 0);          //HZ1
                    coordsArray = fillCoordsArray(HZ1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, 1);         //HZ2
                    coordsArray = fillCoordsArray(HZ2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, -1);        //HZ3
                    coordsArray = fillCoordsArray(HZ3,coordsArray, determinedXYZ);
                    break;
                }
                case LYD: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom CE = (Atom) newResidue.getAtomNode("CE");
                    double [] ce = CE.getXYZ(null);
                    Atom NZ = (Atom) newResidue.getAtomNode("NZ");
                    double [] nz = NZ.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HD2.getXYZ(null);
                    Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
                    double [] hd3 = HD3.getXYZ(null);
                    Atom HE2 = (Atom) newResidue.getAtomNode("HE2");
                    double [] he2 = HE2.getXYZ(null);
                    Atom HE3 = (Atom) newResidue.getAtomNode("HE3");
                    double [] he3 = HE3.getXYZ(null);
                    Atom HZ1 = (Atom) newResidue.getAtomNode("HZ1");
                    double [] hz1 = HZ1.getXYZ(null);
                    Atom HZ2 = (Atom) newResidue.getAtomNode("HZ2");
                    double [] hz2 = HZ2.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dCE_CD, cg, dCE_CD_CG, cb, 180.0, 0);           //CE
                    coordsArray = fillCoordsArray(CE,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dNZ_CE, cd, dNZ_CE_CD, cg, 180.0, 0);           //NZ
                    coordsArray = fillCoordsArray(NZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, 1);          //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ce, 109.4, -1);         //HD3
                    coordsArray = fillCoordsArray(HD3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, 1);          //HE2
                    coordsArray = fillCoordsArray(HE2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ce, dHE_CE, cd, dHE_CE_CD, nz, 108.8, -1);         //HE3
                    coordsArray = fillCoordsArray(HE3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, cd, 180.0, 0);          //HZ1
                    coordsArray = fillCoordsArray(HZ1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nz, dHZ_NZ, ce, dHZ_NZ_CE, hz1, 109.5, 1);         //HZ2
                    coordsArray = fillCoordsArray(HZ2,coordsArray, determinedXYZ);
                    break;
                }
                case ARG: {
                    Atom CB = (Atom) newResidue.getAtomNode("CB");
                    double [] cb = CB.getXYZ(null);
                    Atom CG = (Atom) newResidue.getAtomNode("CG");
                    double [] cg = CG.getXYZ(null);
                    Atom CD = (Atom) newResidue.getAtomNode("CD");
                    double [] cd = CD.getXYZ(null);
                    Atom NE = (Atom) newResidue.getAtomNode("NE");
                    double [] ne = NE.getXYZ(null);
                    Atom CZ = (Atom) newResidue.getAtomNode("CZ");
                    double [] cz = CZ.getXYZ(null);
                    Atom NH1 = (Atom) newResidue.getAtomNode("NH1");
                    double [] nh1 = NH1.getXYZ(null);
                    Atom NH2 = (Atom) newResidue.getAtomNode("NH2");
                    double [] nh2 = NH2.getXYZ(null);
                    Atom HB2 = (Atom) newResidue.getAtomNode("HB2");
                    double [] hb2 = HB2.getXYZ(null);
                    Atom HB3 = (Atom) newResidue.getAtomNode("HB3");
                    double [] hb3 = HB3.getXYZ(null);
                    Atom HG2 = (Atom) newResidue.getAtomNode("HG2");
                    double [] hg2 = HG2.getXYZ(null);
                    Atom HG3 = (Atom) newResidue.getAtomNode("HG3");
                    double [] hg3 = HG3.getXYZ(null);
                    Atom HD2 = (Atom) newResidue.getAtomNode("HD2");
                    double [] hd2 = HG2.getXYZ(null);
                    Atom HD3 = (Atom) newResidue.getAtomNode("HD3");
                    double [] hd3 = HD3.getXYZ(null);
                    Atom HE = (Atom) newResidue.getAtomNode("HE");
                    double [] he = HE.getXYZ(null);
                    Atom HH11 = (Atom) newResidue.getAtomNode("HH11");
                    double [] hh11 = HH11.getXYZ(null);
                    Atom HH12 = (Atom) newResidue.getAtomNode("HH12");
                    double [] hh12 = HH12.getXYZ(null);
                    Atom HH21 = (Atom) newResidue.getAtomNode("HH21");
                    double [] hh21 = HH21.getXYZ(null);
                    Atom HH22 = (Atom) newResidue.getAtomNode("HH22");
                    double [] hh22 = HH22.getXYZ(null);
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
                    determinedXYZ = BondedUtils.determineIntxyz(ca, 1.54, n, 109.5, c, 107.8, 1);
                    coordsArray = fillCoordsArray(CB,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dCG_CB, ca, dCG_CB_CA, n, 180.0, 0);            //CG
                    coordsArray = fillCoordsArray(CG,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dCD_CG, cb, dCD_CG_CB, ca, 180.0, 0);           //CD
                    coordsArray = fillCoordsArray(CD,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dNE_CD, cg, dNE_CD_CG, cb, 180.0, 0);           //NE
                    coordsArray = fillCoordsArray(NE,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne, dCZ_NE, cd, dCZ_NE_CD, cg, 180.0, 0);           //CZ
                    coordsArray = fillCoordsArray(CZ,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz, dNH_CZ, ne, dNH_CZ_NE, cd, 180, 0);            //NH1
                    coordsArray = fillCoordsArray(NH1,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cz, dNH_CZ, ne, dNH_CZ_NE, nh1, 120.0, 1);         //NH2
                    coordsArray = fillCoordsArray(NH2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, 1);          //HB2
                    coordsArray = fillCoordsArray(HB2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cb, dHB_CB, ca, dHB_CB_CA, cg, 109.4, -1);         //HB3
                    coordsArray = fillCoordsArray(HB3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, 1);          //HG2
                    coordsArray = fillCoordsArray(HG2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cg, dHG_CG, cb, dHG_CG_CB, cd, 109.4, -1);         //HG3
                    coordsArray = fillCoordsArray(HG3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ne, 109.4, 1);          //HD2
                    coordsArray = fillCoordsArray(HD2,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(cd, dHD_CD, cg, dHD_CD_CG, ne, 109.4, -1);         //HD3
                    coordsArray = fillCoordsArray(HD3,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(ne, dHE_NE, cd, dHE_NE_CD, cz, 120.0, 1);           //HE
                    coordsArray = fillCoordsArray(HE,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nh1, dHH_NH, cz, dHH_NH_CZ, ne, 180.0, 0);        //HH11
                    coordsArray = fillCoordsArray(HH11,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nh1, dHH_NH, cz, dHH_NH_CZ, hh11, 120.0, 1);      //HH12
                    coordsArray = fillCoordsArray(HH12,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nh2, dHH_NH, cz, dHH_NH_CZ, ne, 180.0, 0);        //HH21
                    coordsArray = fillCoordsArray(HH21,coordsArray, determinedXYZ);
                    
                    determinedXYZ = BondedUtils.determineIntxyz(nh2, dHH_NH, cz, dHH_NH_CZ, hh21, 120.0, 1);      //HH22
                    coordsArray = fillCoordsArray(HH22,coordsArray, determinedXYZ);
                    break;
                }
                default:
                    break;
            }     
        }
        return coordsArray; 
    }

    private double[] fillCoordsArray(Atom atom, double[] coordsArray, double[] determinedXYZ) {
        int XYZIndex = atom.getXYZIndex();
        coordsArray[XYZIndex] = determinedXYZ[0];
        coordsArray[XYZIndex + 1] = determinedXYZ[1];
        coordsArray[XYZIndex + 2] = determinedXYZ[2];
        return coordsArray;
    }
}
