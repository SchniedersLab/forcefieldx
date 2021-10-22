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
package ffx.potential.utils;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.nonbonded.VanDerWaals;

import java.util.*;

public class GetProteinFeatures {
    private static final HashMap<AminoAcid3, String> polarityMap = new HashMap<>();
    private static final HashMap<AminoAcid3, String> acidityMap = new HashMap<>();
    private static final NavigableMap<Double, String> phiToStructure = new TreeMap<>();
    private static final NavigableMap<Double, String> psiToStructure = new TreeMap<>();
    private static final HashMap<AminoAcid3, Double> standardSurfaceArea = new HashMap<>();
    private double phi;
    private double psi;
    private double omega;
    private double totalSurfaceArea = 0.0;


    public GetProteinFeatures() {

    }

    static {
        //Map residue to its polarity
        polarityMap.put(AminoAcid3.ARG, "polar");
        polarityMap.put(AminoAcid3.ASN, "polar");
        polarityMap.put(AminoAcid3.ASP, "polar");
        polarityMap.put(AminoAcid3.ASH, "polar");
        polarityMap.put(AminoAcid3.CYS, "polar");
        polarityMap.put(AminoAcid3.GLN, "polar");
        polarityMap.put(AminoAcid3.GLU, "polar");
        polarityMap.put(AminoAcid3.GLH, "polar");
        polarityMap.put(AminoAcid3.HIS, "polar");
        polarityMap.put(AminoAcid3.HIE, "polar");
        polarityMap.put(AminoAcid3.HID, "polar");
        polarityMap.put(AminoAcid3.LYS, "polar");
        polarityMap.put(AminoAcid3.LYD, "polar");
        polarityMap.put(AminoAcid3.SER, "polar");
        polarityMap.put(AminoAcid3.THR, "polar");
        polarityMap.put(AminoAcid3.TYR, "polar");
        polarityMap.put(AminoAcid3.ALA, "nonpolar");
        polarityMap.put(AminoAcid3.GLY, "nonpolar");
        polarityMap.put(AminoAcid3.ILE, "nonpolar");
        polarityMap.put(AminoAcid3.LEU, "nonpolar");
        polarityMap.put(AminoAcid3.MET, "nonpolar");
        polarityMap.put(AminoAcid3.PHE, "nonpolar");
        polarityMap.put(AminoAcid3.PRO, "nonpolar");
        polarityMap.put(AminoAcid3.TRP, "nonpolar");
        polarityMap.put(AminoAcid3.VAL, "nonpolar");

        //Map residue to its acidity
        acidityMap.put(AminoAcid3.ASP, "acidic");
        acidityMap.put(AminoAcid3.ASH, "acidic");
        acidityMap.put(AminoAcid3.ASN, "acidic");
        acidityMap.put(AminoAcid3.GLU, "acidic");
        acidityMap.put(AminoAcid3.GLH, "acidic");
        acidityMap.put(AminoAcid3.ARG, "basic");
        acidityMap.put(AminoAcid3.HIS, "basic");
        acidityMap.put(AminoAcid3.HIE, "basic");
        acidityMap.put(AminoAcid3.HID, "basic");
        acidityMap.put(AminoAcid3.LYS, "basic");
        acidityMap.put(AminoAcid3.LYD, "basic");
        acidityMap.put(AminoAcid3.GLN, "neutral");
        acidityMap.put(AminoAcid3.GLY, "neutral");
        acidityMap.put(AminoAcid3.ALA, "neutral");
        acidityMap.put(AminoAcid3.VAL, "neutral");
        acidityMap.put(AminoAcid3.LEU, "neutral");
        acidityMap.put(AminoAcid3.ILE, "neutral");
        acidityMap.put(AminoAcid3.SER, "neutral");
        acidityMap.put(AminoAcid3.CYS, "neutral");
        acidityMap.put(AminoAcid3.THR, "neutral");
        acidityMap.put(AminoAcid3.MET, "neutral");
        acidityMap.put(AminoAcid3.PRO, "neutral");
        acidityMap.put(AminoAcid3.PHE, "neutral");
        acidityMap.put(AminoAcid3.TYR, "neutral");
        acidityMap.put(AminoAcid3.TRP, "neutral");

        //Map phi to the Alpha Helix (L) and (R) and Beta Sheet Range
        phiToStructure.put(-180.0, "Extended");
        phiToStructure.put(-150.0, "Structure");
        phiToStructure.put(-50.0, "Structure");
        phiToStructure.put(0.0, "Extended");
        phiToStructure.put(50.0, "Structure");
        phiToStructure.put(70.0, "Structure");
        phiToStructure.put(180.0, "Extended");

        //Map psi to the Alpha Helix (L) and (R) and Beta Sheet Range
        psiToStructure.put(-180.0, "Extended");
        psiToStructure.put(-70.0, "Alpha Helix");
        psiToStructure.put(-50.0, "Alpha Helix");
        psiToStructure.put(0.0, "Extended");
        //psiToStructure.put(30.0, "Alpha Helix (L)");
        //psiToStructure.put(70.0, "Alpha Helix (L)");
        psiToStructure.put(85.0, "Extended");
        psiToStructure.put(100.0, "Beta Sheet");
        psiToStructure.put(150.0, "Beta Sheet");
        //psiToStructure.put(150.0, "Anti-Parallel Beta Sheet");
        psiToStructure.put(180.0, "Extended");

        //Map of exposed surface area to standard from model compounds
        standardSurfaceArea.put(AminoAcid3.ALA, 127.15871);
        standardSurfaceArea.put(AminoAcid3.ARG, 269.42558);
        standardSurfaceArea.put(AminoAcid3.ASH, 158.67385);
        standardSurfaceArea.put(AminoAcid3.ASP, 159.97984);
        standardSurfaceArea.put(AminoAcid3.ASN, 159.82709);
        standardSurfaceArea.put(AminoAcid3.CYS, 154.52801);
        standardSurfaceArea.put(AminoAcid3.GLH, 195.30608);
        standardSurfaceArea.put(AminoAcid3.GLN, 197.75170);
        standardSurfaceArea.put(AminoAcid3.GLU, 195.27081);
        standardSurfaceArea.put(AminoAcid3.GLY, 95.346188);
        standardSurfaceArea.put(AminoAcid3.HID, 202.31302);
        standardSurfaceArea.put(AminoAcid3.HIE, 203.22195);
        standardSurfaceArea.put(AminoAcid3.HIS, 205.40232);
        standardSurfaceArea.put(AminoAcid3.ILE, 196.15872);
        standardSurfaceArea.put(AminoAcid3.LEU, 192.85123);
        standardSurfaceArea.put(AminoAcid3.LYD, 235.49182);
        standardSurfaceArea.put(AminoAcid3.LYS, 236.71473);
        standardSurfaceArea.put(AminoAcid3.MET, 216.53318);
        standardSurfaceArea.put(AminoAcid3.PHE, 229.75038);
        standardSurfaceArea.put(AminoAcid3.PRO, 157.30011);
        standardSurfaceArea.put(AminoAcid3.SER, 137.90720);
        standardSurfaceArea.put(AminoAcid3.THR, 157.33759);
        standardSurfaceArea.put(AminoAcid3.TRP, 262.32819);
        standardSurfaceArea.put(AminoAcid3.TYR, 239.91172);
        standardSurfaceArea.put(AminoAcid3.VAL, 171.89211);
    }

    public String[] saveFeatures(Residue residue, double surfaceArea) {
        String[] features = new String[11];
        String name = residue.getName();
        AminoAcid3 aa3 = residue.getAminoAcid3();
        String acid = acidityMap.getOrDefault(aa3, null);
        String structure = "";
        String phiString = "";
        String psiString = "";
        String omegaString = "";

        if (residue.getNextResidue() == null) {
            //Since phi, psi angles are determined between two residues, the first and last residue will not have values
            //and are default labeled as Extended Secondary Structure
            structure = "Extended";
            getPhi(residue);
            phiString = String.valueOf(phi);
            psiString = null;
            omegaString = null;
        } else if (residue.getPreviousResidue() == null) {
            structure = "Extended";
            getPsi(residue);
            getOmega(residue);
            psiString = String.valueOf(psi);
            omegaString = String.valueOf(omega);
            phiString = null;
        } else {
            getPhi(residue);
            getPsi(residue);
            getOmega(residue);
            phiString = String.valueOf(phi);
            psiString = String.valueOf(psi);
            omegaString = String.valueOf(omega);
            structure = getSecondaryStructure();
        }

        totalSurfaceArea += surfaceArea;
        String surfaceAreaString = String.valueOf(surfaceArea);

        double standSurfaceArea = standardSurfaceArea.getOrDefault(residue.getAminoAcid3(), 0.0);
        String normalizedSA = "";
        if(standSurfaceArea != 0.0){
            normalizedSA = String.valueOf(surfaceArea / standSurfaceArea);
        }
        String confidence =  String.valueOf(getConfidenceScore(residue));

        features[0] = name;
        features[1] = String.valueOf(residue.getResidueNumber());
        features[2] = polarityMap.getOrDefault(aa3, "null");
        features[3] = acid;
        features[4] = structure;
        features[5] = phiString;
        features[6] = psiString;
        features[7] = omegaString;
        features[8] = surfaceAreaString;
        features[9] = normalizedSA;
        features[10] = confidence;

        return features;
    }

    public void getPhi(Residue currentRes) {
        Residue previousRes = currentRes.getPreviousResidue();
        double[] cCoor = new double[3];
        double[] nCoor = new double[3];
        double[] caCoor = new double[3];
        double[] c2Coor = new double[3];
        phi = (DoubleMath.dihedralAngle(previousRes.getAtomByName("C", true).getXYZ(cCoor),
                currentRes.getAtomByName("N", true).getXYZ(nCoor),
                currentRes.getAtomByName("CA", true).getXYZ(caCoor),
                currentRes.getAtomByName("C", true).getXYZ(c2Coor))) * 180 / Math.PI;
    }

    public void getPsi(Residue currentRes) {
        //res[0] is always current Res
        Residue nextRes = currentRes.getNextResidue();
        double[] nCoor = new double[3];
        double[] caCoor = new double[3];
        double[] cCoor = new double[3];
        double[] n2Coor = new double[3];
        psi = (DoubleMath.dihedralAngle(currentRes.getAtomByName("N", true).getXYZ(nCoor),
                currentRes.getAtomByName("CA", true).getXYZ(caCoor),
                currentRes.getAtomByName("C", true).getXYZ(cCoor),
                nextRes.getAtomByName("N", true).getXYZ(n2Coor))) * 180 / Math.PI;
    }

    public void getOmega(Residue currentRes) {
        //res[0] is always current Res
        Residue nextRes = currentRes.getNextResidue();
        double[] ca1Coor = new double[3];
        double[] cCoor = new double[3];
        double[] nCoor = new double[3];
        double[] ca2Coor = new double[3];
        omega = (DoubleMath.dihedralAngle(currentRes.getAtomByName("CA", true).getXYZ(ca1Coor),
                currentRes.getAtomByName("C", true).getXYZ(cCoor),
                nextRes.getAtomByName("N", true).getXYZ(nCoor),
                nextRes.getAtomByName("CA", true).getXYZ(ca2Coor))) * 180 / Math.PI;
    }

    //Data on phi,psi and secondary structure correlations is from Tamar Schlick's
    // Molecular Modeling and Simulation Textbook. Approx. page 108.
    //Alpha Helix: phi,psi of approximately -60,-50. Tighter helices are around -50,-25.
    //Looser helices are around -60,-70.

    //Collagen Helix: phi,psi of approximately -60,+125.

    //Parallel Beta Sheet: phi,psi of approximately -120,+115
    //Antiparallel Beta Sheet: phi,psi of approximately -140,+135
    //Get surface area region to get surface area
    public String getSecondaryStructure() {
        String secondaryStructure = "";
        //Use phi, psi ranges to determine which is the most likely secondary structure based on ramachandran values
        Double lowPhiKey = phiToStructure.floorKey(phi);
        String lowPhiStruct = phiToStructure.get(lowPhiKey);
        Double highPhiKey = phiToStructure.ceilingKey(phi);
        String highPhiStruct = phiToStructure.get(highPhiKey);
        if (lowPhiStruct.equals("Extended") || highPhiStruct.equals("Extended")) {
            secondaryStructure = "Extended";
        } else {
            Double lowPsiKey = psiToStructure.floorKey(psi);
            String lowPsiStruct = psiToStructure.get(lowPsiKey);
            Double highPsiKey = psiToStructure.ceilingKey(psi);
            String highPsiStruct = psiToStructure.get(highPsiKey);
            if (lowPsiStruct.equals("Extended") || highPsiStruct.equals("Extended")) {
                secondaryStructure = "Extended";
            } else {
                secondaryStructure = lowPsiStruct;
            }
        }
        return secondaryStructure;
    }

    /*public double getSurfaceArea(Residue residue) {
        List<Atom> atoms = residue.getAtomList();
        int nAtoms = atoms.size();
        int endIndex = surfaceAreaIndex + nAtoms;
        double sumResidueArea = 0.0;

        for (int i = surfaceAreaIndex; i < endIndex; i++) {
            sumResidueArea += surfaceAreaArray[i];
        }
        surfaceAreaIndex += nAtoms;
        totalSurfaceArea += sumResidueArea;

        return sumResidueArea;
    }*/

    public double getTotalSurfaceArea(){
        return totalSurfaceArea;
    }

    public double getConfidenceScore(Residue currentRes){
        double bFactor = 0;
        bFactor = currentRes.getAtomByName("CA", true).getTempFactor();
        return bFactor;
    }


}
