// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.bonded;

import static ffx.potential.bonded.AminoAcidUtils.AA_CB;
import static ffx.potential.bonded.BondedUtils.findAtomType;

import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.AminoAcidUtils.HisStates;
import ffx.potential.bonded.AminoAcidUtils.HistidineAtomNames;
import ffx.potential.bonded.AminoAcidUtils.LysStates;
import ffx.potential.bonded.AminoAcidUtils.LysineAtomNames;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;

/**
 * Utilities for interpolating between Amino Acid protonation and tautomer states.
 *
 * @author Michael Schnieders
 * @author Andrew Thiel
 * @since 1.0
 */
public class ConstantPhUtils {

  private static final MultipoleType zeroMultipoleType =
      new MultipoleType(MultipoleType.zeroM, new int[] {0, 0, 0},
          MultipoleFrameDefinition.ZTHENX, false);

  private static final PolarizeType zeroPolarizeType =
      new PolarizeType(0, 0.0, 0.39, new int[] {0});

  enum AspStates {
    ASP, ASH1, ASH2
  }

  enum GluStates {
    GLU, GLH1, GLH2
  }

  /** Constant <code>AspartateAtomNames</code> */
  private enum AspartateAtomNames {
    CB(0, 0, 0),
    HB2(1, 1, 1),
    HB3(1, 1, 1),
    CG(2, 2, 2),
    OD1(3, 4, 3),
    OD2(3, 3, 4),
    HD1(-1, 5, -1),
    HD2(-1, -1, 5);

    /**
     * Biotype offset relative to the CB biotype for charged aspartate (ASP).
     */
    private int offsetASP;

    /**
     * Biotype offset relative to the CB biotype for neutral aspartic acid protonated on OD1 (ASH1).
     * <p>
     * This is set to negative -1 for the OD2 hydrogen.
     */
    private int offsetASH1;

    /**
     * Biotype offset relative to the CB biotype for neutral aspartic acid protonated on OD2 (ASH2).
     * <p>
     * This is set to negative -1 for the OD1 hydrogen.
     */
    private int offsetASH2;

    public int getOffset(AspStates state) {
      if (state == AspStates.ASP) {
        return offsetASP;
      } else if (state == AspStates.ASH1) {
        return offsetASH1;
      } else {
        return offsetASH2;
      }
    }

    /**
     * Init the Histidine atom names.
     *
     * @param offsetASP Biotype relative to the CB biotype for ASP.
     * @param offsetASH1 Biotype relative to the CB biotype for ASH.
     * @param offsetASH2 Biotype relative to the CB biotype for ASH.
     */
    AspartateAtomNames(int offsetASP, int offsetASH1, int offsetASH2) {
      this.offsetASP = offsetASP;
      this.offsetASH1 = offsetASH1;
      this.offsetASH2 = offsetASH2;
    }
  }

  /** Constant <code>GlutamateAtomNames</code> */
  private enum GlutamateAtomNames {
    CB(0, 0, 0),
    HB2(1, 1, 1),
    HB3(1, 1, 1),
    CG(2, 2, 2),
    HG2(3, 3, 3),
    HG3(3, 3, 3),
    CD(4, 4, 4),
    OE1(5, 6, 5),
    OE2(5, 5, 6),
    HE1(-1, 7, -1),
    HE2(-1, -1, 7);

    /**
     * Biotype offset relative to the CB biotype for charged Glutamate (GLU).
     */
    private int offsetGLU;

    /**
     * Biotype offset relative to the CB biotype for neutral Glutamate acid protonated on OE1
     * (GLU1).
     * <p>
     * This is set to negative -1 for the OE2 hydrogen.
     */
    private int offsetGLH1;

    /**
     * Biotype offset relative to the CB biotype for neutral Glutamate acid protonated on OE2
     * (GLU2).
     * <p>
     * This is set to negative -1 for the OE1 hydrogen.
     */
    private int offsetGLH2;

    public int getOffset(GluStates state) {
      if (state == GluStates.GLU) {
        return offsetGLU;
      } else if (state == GluStates.GLH1) {
        return offsetGLH1;
      } else {
        return offsetGLH2;
      }
    }

    /**
     * Init the Glutamate atom names.
     *
     * @param offsetGLU Biotype relative to the CB biotype for GLU.
     * @param offsetGLH1 Biotype relative to the CB biotype for GLH.
     * @param offsetGLH2 Biotype relative to the CB biotype for GLH.
     */
    GlutamateAtomNames(int offsetGLU, int offsetGLH1, int offsetGLH2) {
      this.offsetGLU = offsetGLU;
      this.offsetGLH1 = offsetGLH1;
      this.offsetGLH2 = offsetGLH2;
    }
  }

  /**
   * Lysine atom types.
   */
  private final int lysTypes = LysineAtomNames.values().length;
  private final AtomType[][] lysAtomTypes = new AtomType[2][lysTypes];
  private final MultipoleType[][] lysMultipoleTypes = new MultipoleType[2][lysTypes];
  private final PolarizeType[][] lysPolarizeTypes = new PolarizeType[2][lysTypes];

  /**
   * Histidine atom types.
   */
  private final int hisTypes = HistidineAtomNames.values().length;
  private final AtomType[][] hisAtomTypes = new AtomType[3][hisTypes];
  private final MultipoleType[][] hisMultipoleTypes = new MultipoleType[3][hisTypes];
  private final PolarizeType[][] hisPolarizeTypes = new PolarizeType[3][hisTypes];

  /**
   * Aspartic acid atom types.
   */
  private final int aspTypes = AspartateAtomNames.values().length;
  private final AtomType[][] aspAtomTypes = new AtomType[3][aspTypes];
  private final MultipoleType[][] aspMultipoleTypes = new MultipoleType[3][aspTypes];
  private final PolarizeType[][] aspPolarizeTypes = new PolarizeType[3][aspTypes];

  /**
   * Glutamic acid atom types.
   */
  private final int gluTypes = GlutamateAtomNames.values().length;
  private final AtomType[][] gluAtomTypes = new AtomType[3][gluTypes];
  private final MultipoleType[][] gluMultipoleTypes = new MultipoleType[3][gluTypes];
  private final PolarizeType[][] gluPolarizeTypes = new PolarizeType[3][gluTypes];

  private final ForceField forceField;

  public ConstantPhUtils(ForceField forceField) {
    this.forceField = forceField;

    // Populate the Lysine types.
    constructLYSState(AA_CB[AminoAcid3.LYS.ordinal()], LysStates.LYS);
    constructLYSState(AA_CB[AminoAcid3.LYD.ordinal()], LysStates.LYD);

    // Populate the Histidine types.
    constructHISState(AA_CB[AminoAcid3.HIS.ordinal()], HisStates.HIS);
    constructHISState(AA_CB[AminoAcid3.HID.ordinal()], HisStates.HID);
    constructHISState(AA_CB[AminoAcid3.HIE.ordinal()], HisStates.HIE);

    // Populate the Aspartic acid types.
    constructASPState(AA_CB[AminoAcid3.ASP.ordinal()], AspStates.ASP);
    constructASPState(AA_CB[AminoAcid3.ASH.ordinal()], AspStates.ASH1);
    constructASPState(AA_CB[AminoAcid3.ASH.ordinal()], AspStates.ASH2);

    // Populate the Glutamic acid types.
    constructGLUState(AA_CB[AminoAcid3.GLU.ordinal()], GluStates.GLU);
    constructGLUState(AA_CB[AminoAcid3.GLH.ordinal()], GluStates.GLH1);
    constructGLUState(AA_CB[AminoAcid3.GLH.ordinal()], GluStates.GLH2);
  }

  public double[] getMultipole(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda, double[] multipole) {
    switch (aminoAcid3) {
      case LYS:
        double[] lys = lysMultipoleTypes[LysStates.LYS.ordinal()][atomIndex].getMultipole();
        double[] lyd = lysMultipoleTypes[LysStates.LYD.ordinal()][atomIndex].getMultipole();
        for (int i=0; i<multipole.length; i++) {
          multipole[i] = titrationLambda * lys[i] + (1.0 - titrationLambda) * lyd[i];
        }
        break;
      case HIS:
        double[] his = hisMultipoleTypes[HisStates.HIS.ordinal()][atomIndex].getMultipole();
        double[] hid = hisMultipoleTypes[HisStates.HID.ordinal()][atomIndex].getMultipole();
        double[] hie = hisMultipoleTypes[HisStates.HIE.ordinal()][atomIndex].getMultipole();
        for (int i=0; i<multipole.length; i++) {
          //multipole[i] = titrationLambda * lys[i] + (1.0 - titrationLambda) * lyd[i];
        }
        break;
      case ASD:
        break;
      case GLD:
        break;
      default:
        return multipole;
    }
    return multipole;
  }

  public double[] getMultipoleTitrationDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda, double[] multipole) {
    switch (aminoAcid3) {
      case LYS:
        double[] lys = lysMultipoleTypes[LysStates.LYS.ordinal()][atomIndex].getMultipole();
        double[] lyd = lysMultipoleTypes[LysStates.LYD.ordinal()][atomIndex].getMultipole();
        for (int i=0; i<multipole.length; i++) {
          multipole[i] = lys[i] - lyd[i];
        }
        break;
      case HIS:
        break;
      case ASD:
        break;
      case GLD:
        break;
      default:
        return multipole;
    }
    return multipole;
  }

  public double[] getMultipoleTautomerDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda, double[] multipole) {
    switch (aminoAcid3) {
      case HIS:
        break;
      case ASD:
        break;
      case GLD:
        break;
      case LYS: // No tautomers for LYS.
      default:
        return multipole;
    }
    return multipole;
  }

  public double getPolarizability(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda, double defaultPolarizability) {
    switch (aminoAcid3) {
      case LYS:
        double lys = lysPolarizeTypes[LysStates.LYS.ordinal()][atomIndex].polarizability;
        double lyd = lysPolarizeTypes[LysStates.LYD.ordinal()][atomIndex].polarizability;
        return titrationLambda * lys + (1.0 - titrationLambda) * lyd;
      case HIS:
        break;
      case ASD:
        break;
      case GLD:
        break;
      default:
        return defaultPolarizability;
    }
    return defaultPolarizability;
  }

  public double getPolarizabilityTitrationDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda) {
    switch (aminoAcid3) {
      case LYS:
        double lys = lysPolarizeTypes[LysStates.LYS.ordinal()][atomIndex].polarizability;
        double lyd = lysPolarizeTypes[LysStates.LYD.ordinal()][atomIndex].polarizability;
        return lys - lyd;
      case HIS:
        break;
      case ASD:
        break;
      case GLD:
        break;
      default:
        return 0.0;
    }
    return 0.0;
  }

  public double getPolarizabilityTautomerDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda) {
    switch (aminoAcid3) {
      case HIS:
        break;
      case ASD:
        break;
      case GLD:
        break;
      case LYS: // No tautomers for LYS.
      default:
        return 0.0;
    }
    return 0.0;
  }

  private void constructHISState(int biotypeCB, HisStates hisState) {
    int state = hisState.ordinal();
    for (HistidineAtomNames atomName : HistidineAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffsetHIS(hisState);
      if (offset < 0) {
        hisAtomTypes[state][index] = null;
        // Zero out the MultipoleType and Polarizetype.
        hisMultipoleTypes[state][index] = zeroMultipoleType;
        hisPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        hisAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = hisAtomTypes[state][index].getKey();
        hisMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        hisPolarizeTypes[state][index] = forceField.getPolarizeType(key);
      }
    }
  }

  private void constructLYSState(int biotypeCB, LysStates lysState) {
    int state = lysState.ordinal();
    for (LysineAtomNames atomName : LysineAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffsetLYS(lysState);
      if (offset < 0) {
        // Set the AtomType to null.
        lysAtomTypes[state][index] = null;
        // Zero out the MultipoleType and Polarizetype.
        lysMultipoleTypes[state][index] = zeroMultipoleType;
        lysPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        lysAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = lysAtomTypes[state][index].getKey();
        lysMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        lysPolarizeTypes[state][index] = forceField.getPolarizeType(key);
      }
    }
  }

  private void constructASPState(int biotypeCB, AspStates aspState) {
    int state = aspState.ordinal();
    for (AspartateAtomNames atomName : AspartateAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffset(aspState);
      if (offset < 0) {
        // Set the AtomType to null.
        aspAtomTypes[state][index] = null;
        // Zero out the MultipoleType and Polarizetype.
        aspMultipoleTypes[state][index] = zeroMultipoleType;
        aspPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        aspAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = lysAtomTypes[state][index].getKey();
        aspMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        aspPolarizeTypes[state][index] = forceField.getPolarizeType(key);
      }
    }
  }

  private void constructGLUState(int biotypeCB, GluStates gluState) {
    int state = gluState.ordinal();
    for (GlutamateAtomNames atomName : GlutamateAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffset(gluState);
      if (offset < 0) {
        // Set the AtomType to null.
        gluAtomTypes[state][index] = null;
        // Zero out the MultipoleType and Polarizetype.
        gluMultipoleTypes[state][index] = zeroMultipoleType;
        gluPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        gluAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = lysAtomTypes[state][index].getKey();
        gluMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        gluPolarizeTypes[state][index] = forceField.getPolarizeType(key);
      }
    }
  }

}
