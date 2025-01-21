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
package ffx.potential.utils;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.NavigableMap;
import java.util.TreeMap;

public class GetProteinFeatures {

  private static final HashMap<AminoAcid3, String> polarityMap = new HashMap<>();
  private static final HashMap<AminoAcid3, String> acidityMap = new HashMap<>();
  private static final NavigableMap<Double, String> phiToStructure = new TreeMap<>();
  private static final NavigableMap<Double, String> psiToStructure = new TreeMap<>();
  private static final HashMap<AminoAcid3, Double> standardSurfaceArea = new HashMap<>();
  private static final HashMap<String, AminoAcid3> aminoAcidCodes = new HashMap<>();
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
    acidityMap.put(AminoAcid3.LEU, "neutral");
    acidityMap.put(AminoAcid3.GLN, "neutral");
    acidityMap.put(AminoAcid3.GLY, "neutral");
    acidityMap.put(AminoAcid3.ALA, "neutral");
    acidityMap.put(AminoAcid3.VAL, "neutral");
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

    // Map amino acid codes from 1 letter to 3 letter AA3
    aminoAcidCodes.put("A", AminoAcid3.ALA);
    aminoAcidCodes.put("R", AminoAcid3.ARG);
    aminoAcidCodes.put("N", AminoAcid3.ASN);
    aminoAcidCodes.put("D", AminoAcid3.ASP);
    aminoAcidCodes.put("C", AminoAcid3.CYS);
    aminoAcidCodes.put("E", AminoAcid3.GLU);
    aminoAcidCodes.put("Q", AminoAcid3.GLN);
    aminoAcidCodes.put("G", AminoAcid3.GLY);
    aminoAcidCodes.put("H", AminoAcid3.HIS);
    aminoAcidCodes.put("I", AminoAcid3.ILE);
    aminoAcidCodes.put("L", AminoAcid3.LEU);
    aminoAcidCodes.put("K", AminoAcid3.LYS);
    aminoAcidCodes.put("M", AminoAcid3.MET);
    aminoAcidCodes.put("F", AminoAcid3.PHE);
    aminoAcidCodes.put("P", AminoAcid3.PRO);
    aminoAcidCodes.put("S", AminoAcid3.SER);
    aminoAcidCodes.put("T", AminoAcid3.THR);
    aminoAcidCodes.put("W", AminoAcid3.TRP);
    aminoAcidCodes.put("Y", AminoAcid3.TYR);
    aminoAcidCodes.put("V", AminoAcid3.VAL);
  }

  /**
   * Make a string array of surface area and additional selected features (phi,psi,omega,and
   * structure annotations)
   *
   * @param residue Residue
   * @param surfaceArea residue surface area
   * @param includeAngles select angles
   * @param includeStructure select structure annotation
   * @return String array of features
   */
  public String[] saveFeatures(Residue residue, double surfaceArea, boolean includeAngles,
      boolean includeStructure) {
    int nFeat = 3;
    if (includeAngles) {
      nFeat += 3;
    }
    if (includeStructure) {
      nFeat += 1;
    }

    String[] features = new String[nFeat];
    String structure;
    String phiString;
    String psiString;
    String omegaString;

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
    if (standSurfaceArea != 0.0) {
      normalizedSA = String.valueOf(surfaceArea / standSurfaceArea);
    }
    String confidence = String.valueOf(getConfidenceScore(residue));

    features[0] = surfaceAreaString;
    features[1] = normalizedSA;
    features[2] = confidence;
    if (includeAngles) {
      features[3] = phiString;
      features[4] = psiString;
      features[5] = omegaString;
      if (includeStructure) {
        features[6] = structure;
      }
    } else if (includeStructure) {
      features[3] = structure;
    }
    return features;
  }

  /**
   * Get the phi angle of a residue
   *
   * @param currentRes current residue
   */
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

  /**
   * Get the psi angle of a residue
   *
   * @param currentRes current residue
   */
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

  /**
   * Get the omega angle of a residue
   *
   * @param currentRes current residue
   */
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

  /**
   * Get the secondary structure annotation from the ramachandran angle map
   *
   * @return string of secondary structure
   */
  public String getSecondaryStructure() {
    String secondaryStructure;
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

  /**
   * Get the total surface area for the protein
   *
   * @return The total surface area.
   */
  public double getTotalSurfaceArea() {
    return totalSurfaceArea;
  }

  /**
   * Get the alphafold confidence score or b-factor from an X-ray model.
   *
   * @param currentRes current residue
   * @return confidence/b-factor value
   */
  public double getConfidenceScore(Residue currentRes) {
    return currentRes.getAtomByName("CA", true).getTempFactor();
  }

  /**
   * Use the ddgun output file to get the amino acid changes
   *
   * @param ddgun List of lines from ddGun output file
   * @return List of NP Changes
   */
  public List<String> ddgunToNPChange(List<String> ddgun) {
    List<String> npChanges = new ArrayList<>();
    for (String s : ddgun) {
      String[] splits = s.split("\t");
      String currentNP = splits[2];
      String wt = String.valueOf(currentNP.charAt(0));
      String mut = String.valueOf(currentNP.charAt(currentNP.length() - 1));
      String pos = currentNP.substring(1, currentNP.length() - 1);
      String wt3Letter = aminoAcidCodes.get(wt).toString();
      String mut3Letter = aminoAcidCodes.get(mut).toString();
      String wildType = wt3Letter.charAt(0) + wt3Letter.substring(1, 3).toLowerCase();
      String mutant = mut3Letter.charAt(0) + mut3Letter.substring(1, 3).toLowerCase();
      String npChange = "p." + wildType + pos + mutant;
      npChanges.add(npChange);
    }
    return npChanges;
  }

  /**
   * Get ddgun values from ddgun file
   *
   * @param ddgun List of lines from ddGun output file
   * @return List of ddGun values (raw and abs value)
   */
  public List<Double[]> getDDGunValues(List<String> ddgun) {
    List<Double[]> values = new ArrayList<>();
    for (String s : ddgun) {
      String[] splits = s.split("\t");
      Double[] value = new Double[2];
      value[0] = Double.parseDouble(splits[3]) * -1.0;
      value[1] = Math.abs(Double.parseDouble(splits[3]));
      values.add(value);
    }
    return values;
  }

  /**
   * Get the polarity and acidity changes
   *
   * @param npChanges list of protein changes
   * @param includePolarity select polarity
   * @param includeAcidity select acidity
   * @return list of polarity and acidity changes
   */
  public List<String[]> getPolarityAndAcidityChange(List<String> npChanges, boolean includePolarity,
      boolean includeAcidity) {
    List<String[]> polarityAndAcidity = new ArrayList<>();
    for (String npChange : npChanges) {
      String change = npChange.split("p\\.")[1].toUpperCase(Locale.ROOT);
      String[] value = new String[2];
      AminoAcid3 wt = AminoAcid3.valueOf(change.substring(0, 3));
      AminoAcid3 mut = AminoAcid3.valueOf(change.substring(change.length() - 3));
      if (includeAcidity) {
        if (acidityMap.get(wt).equals("basic") && acidityMap.get(mut).equals("neutral")) {
          value[0] = "bn";
        } else if (acidityMap.get(wt).equals("neutral") && acidityMap.get(mut).equals("basic")) {
          value[0] = "nb";
        } else if (acidityMap.get(wt).equals("acidic") && acidityMap.get(mut).equals("neutral")) {
          value[0] = "an";
        } else if (acidityMap.get(wt).equals("neutral") && acidityMap.get(mut).equals("acidic")) {
          value[0] = "na";
        } else if (acidityMap.get(wt).equals("basic") && acidityMap.get(mut).equals("acidic")) {
          value[0] = "ba";
        } else if (acidityMap.get(wt).equals("acidic") && acidityMap.get(mut).equals("basic")) {
          value[0] = "ab";
        } else if (acidityMap.get(wt).equals(acidityMap.get(mut))) {
          value[0] = "=";
        }
      } else {
        value[0] = null;
      }

      if (includePolarity) {
        if (polarityMap.get(wt).equals("polar") && polarityMap.get(mut).equals("nonpolar")) {
          value[1] = "-";
        } else if (polarityMap.get(wt).equals("nonpolar") && polarityMap.get(mut).equals("polar")) {
          value[1] = "+";
        } else {
          value[1] = "=";
        }
      } else {
        value[1] = null;
      }

      polarityAndAcidity.add(value);
    }
    return polarityAndAcidity;
  }


}
