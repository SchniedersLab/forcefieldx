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
import static java.lang.String.format;

import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utilities for interpolating between Amino Acid protonation and tautomer states.
 *
 * @author Michael Schnieders
 * @author Andrew Thiel
 * @since 1.0
 */
public class ConstantPhUtils {

  private static final Logger logger = Logger.getLogger(ConstantPhUtils.class.getName());

  private static final MultipoleType zeroMultipoleType =
      new MultipoleType(MultipoleType.zeroM, new int[] {0, 0, 0},
          MultipoleFrameDefinition.ZTHENX, false);

  private static final PolarizeType zeroPolarizeType =
      new PolarizeType(0, 0.0, 0.39, new int[] {0});

  private static final AtomType deprotonatedAtomType = new AtomType(0, 0,
      "H", "\"Deprotonated Hydrogen\"", 1, 1.0080, 1);

  enum AspStates {
    ASP, ASH1, ASH2
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

  enum GluStates {
    GLU, GLH1, GLH2
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

  enum LysStates {
    LYD, LYS
  }

  /** Constant <code>lysineAtoms</code> */
  public enum LysineAtomNames {
    CB(0, 0), HB2(1, 1), HB3(1, 1),
    CG(2, 2), HG2(3, 3), HG3(3, 3),
    CD(4, 4), HD2(5, 5), HD3(5, 5),
    CE(6, 6), HE2(7, 7), HE3(7, 7),
    NZ(8, 8), HZ1(9, 9), HZ2(9, 9),
    HZ3(9, -1);

    /**
     * Biotype offset relative to the CB biotype for LYS.
     */
    private int offsetLYS;

    /**
     * Biotype offset relative to the CB biotype for LYD.
     */
    private int offsetLYD;

    public int getOffsetLYS(LysStates state) {
      if (state == LysStates.LYS) {
        return offsetLYS;
      } else {
        return offsetLYD;
      }
    }

    /**
     * Init the Lysine atom names.
     *
     * @param offsetLYS Biotype offset relative to the CB biotype for LYS.
     * @param offsetLYD Biotype offset relative to the CB biotype for LYD.
     */
    LysineAtomNames(int offsetLYS, int offsetLYD) {
      this.offsetLYS = offsetLYS;
      this.offsetLYD = offsetLYD;
    }
  }

  enum HisStates {
    HIS, HID, HIE
  }

  /** Constant <code>HistidineAtoms</code> */
  public enum HistidineAtomNames {
    // HIS, HID, HIE
    CB(0, 0, 0),
    HB2(1, 1, 1),
    HB3(1, 1, 1),
    CG(2, 2, 2),
    ND1(3, 3, 3),
    // No HD1 proton for HIE; HIE HD1 offset is -1.
    HD1(4, 4, -1),
    CD2(5, 5, 4),
    HD2(6, 6, 5),
    CE1(7, 7, 6),
    HE1(8, 8, 7),
    NE2(9, 9, 8),
    // No HE2 proton for HID; HID HE2 offset is -1
    HE2(10, -1, 9);

    /**
     * Biotype offset relative to the CB biotype for charged histidine (HIS).
     */
    private int offsetHIS;

    /**
     * Biotype offset relative to the CB biotype for neutral histidine protonated on the delta
     * nitrogren (HID).
     * <p>
     * This is set to negative -1 for the epsilon hydrogen.
     */
    private int offsetHID;

    /**
     * Biotype offset relative to the CB biotype for neutral histidine protonated the epsilon
     * nitrogen (HIE).
     * <p>
     * This is set to negative -1 for the delta hydrogen.
     */
    private int offsetHIE;

    public int getOffsetHIS(HisStates state) {
      if (state == HisStates.HIS) {
        return offsetHIS;
      } else if (state == HisStates.HID) {
        return offsetHID;
      } else {
        return offsetHIE;
      }
    }

    /**
     * Init the Histidine atom names.
     *
     * @param offsetHIS Biotype relative to the CB biotype for HIS.
     * @param offsetHID Biotype relative to the CB biotype for HID.
     * @param offsetHIE Biotype relative to the CB biotype for HIE.
     */
    HistidineAtomNames(int offsetHIS, int offsetHID, int offsetHIE) {
      this.offsetHIS = offsetHIS;
      this.offsetHID = offsetHID;
      this.offsetHIE = offsetHIE;
    }
  }

  /**
   * Lysine atom types.
   */
  private final int nLysTypes = LysineAtomNames.values().length;
  private final int nLysStates = LysStates.values().length;
  private final AtomType[][] lysAtomTypes = new AtomType[nLysStates][nLysTypes];
  private final MultipoleType[][] lysMultipoleTypes = new MultipoleType[nLysStates][nLysTypes];
  private final PolarizeType[][] lysPolarizeTypes = new PolarizeType[nLysStates][nLysTypes];

  /**
   * Histidine atom types.
   */
  private final int nHisTypes = HistidineAtomNames.values().length;
  private final int nHisStates = HisStates.values().length;
  private final AtomType[][] hisAtomTypes = new AtomType[nHisStates][nHisTypes];
  private final MultipoleType[][] hisMultipoleTypes = new MultipoleType[nHisStates][nHisTypes];
  private final PolarizeType[][] hisPolarizeTypes = new PolarizeType[nHisStates][nHisTypes];

  /**
   * Aspartic acid atom types.
   */
  private final int nAspTypes = AspartateAtomNames.values().length;
  private final int nAspStates = AspStates.values().length;
  private final AtomType[][] aspAtomTypes = new AtomType[nAspStates][nAspTypes];
  private final MultipoleType[][] aspMultipoleTypes = new MultipoleType[nAspStates][nAspTypes];
  private final PolarizeType[][] aspPolarizeTypes = new PolarizeType[nAspStates][nAspTypes];

  /**
   * Glutamic acid atom types.
   */
  private final int nGluTypes = GlutamateAtomNames.values().length;
  private final int nGluStates = GluStates.values().length;
  private final AtomType[][] gluAtomTypes = new AtomType[nGluStates][nGluTypes];
  private final MultipoleType[][] gluMultipoleTypes = new MultipoleType[nGluStates][nGluTypes];
  private final PolarizeType[][] gluPolarizeTypes = new PolarizeType[nGluStates][nGluTypes];

  private final ForceField forceField;

  public ConstantPhUtils(ForceField forceField) {
    this.forceField = forceField;

    // Populate the Lysine types.
    constructLYSState(AA_CB[AminoAcid3.LYS.ordinal()], LysStates.LYS);
    constructLYSState(AA_CB[AminoAcid3.LYD.ordinal()], LysStates.LYD);
    checkMultipoleFrames("LYS", lysAtomTypes, lysPolarizeTypes, lysMultipoleTypes);

    // Populate the Histidine types.
    constructHISState(AA_CB[AminoAcid3.HIS.ordinal()], HisStates.HIS);
    constructHISState(AA_CB[AminoAcid3.HID.ordinal()], HisStates.HID);
    constructHISState(AA_CB[AminoAcid3.HIE.ordinal()], HisStates.HIE);
    checkMultipoleFrames("HIS", hisAtomTypes, hisPolarizeTypes, hisMultipoleTypes);

    // Populate the Aspartic acid types.
    constructASPState(AA_CB[AminoAcid3.ASP.ordinal()], AspStates.ASP);
    constructASPState(AA_CB[AminoAcid3.ASH.ordinal()], AspStates.ASH1); // First ASH Tautomer
    constructASPState(AA_CB[AminoAcid3.ASH.ordinal()], AspStates.ASH2); // Second ASH Tautomer
    checkMultipoleFrames("ASP", aspAtomTypes, aspPolarizeTypes, aspMultipoleTypes);

    // Populate the Glutamic acid types.
    constructGLUState(AA_CB[AminoAcid3.GLU.ordinal()], GluStates.GLU);
    constructGLUState(AA_CB[AminoAcid3.GLH.ordinal()], GluStates.GLH1); // First GLH Tautomer
    constructGLUState(AA_CB[AminoAcid3.GLH.ordinal()], GluStates.GLH2); // Second GLH Tautomer
    checkMultipoleFrames("GLU", gluAtomTypes, gluPolarizeTypes, gluMultipoleTypes);
  }

  public double[] getMultipole(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda, double[] multipole) {
    switch (aminoAcid3) {
      case LYS:
        double[] lys = lysMultipoleTypes[LysStates.LYS.ordinal()][atomIndex].getMultipole();
        double[] lyd = lysMultipoleTypes[LysStates.LYD.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * lys[i] + (1.0 - titrationLambda) * lyd[i];
        }
        break;
      case HIS:
        double[] his = hisMultipoleTypes[HisStates.HIS.ordinal()][atomIndex].getMultipole();
        double[] hid = hisMultipoleTypes[HisStates.HID.ordinal()][atomIndex].getMultipole();
        double[] hie = hisMultipoleTypes[HisStates.HIE.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * his[i] + (1.0 - titrationLambda) * (tautomerLambda * hid[i] + (1 - tautomerLambda) * hie[i]);
        }
        break;
      case ASD:
        double[] asp =  aspMultipoleTypes[AspStates.ASP.ordinal()][atomIndex].getMultipole();
        double[] ash1 = aspMultipoleTypes[AspStates.ASH1.ordinal()][atomIndex].getMultipole();
        double[] ash2 = aspMultipoleTypes[AspStates.ASH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (tautomerLambda * ash1[i] + (1 - tautomerLambda) * ash2[i]) + (1.0 - titrationLambda) * asp[i];
        }
        break;
      case GLD:
        double[] glu =  gluMultipoleTypes[GluStates.GLU.ordinal()][atomIndex].getMultipole();
        double[] glh1 = gluMultipoleTypes[GluStates.GLH1.ordinal()][atomIndex].getMultipole();
        double[] glh2 = gluMultipoleTypes[GluStates.GLH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (tautomerLambda * glh1[i] + (1 - tautomerLambda) * glh2[i]) + (1.0 - titrationLambda) * glu[i];
        }
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
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = lys[i] - lyd[i];
        }
        break;
      case HIS:
        double[] his = hisMultipoleTypes[HisStates.HIS.ordinal()][atomIndex].getMultipole();
        double[] hid = hisMultipoleTypes[HisStates.HID.ordinal()][atomIndex].getMultipole();
        double[] hie = hisMultipoleTypes[HisStates.HIE.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = his[i] - (tautomerLambda * hid[i] + (1 - tautomerLambda) * hie[i]);
        }
        break;
      case ASD:
        double[] asp =  aspMultipoleTypes[AspStates.ASP.ordinal()][atomIndex].getMultipole();
        double[] ash1 = aspMultipoleTypes[AspStates.ASH1.ordinal()][atomIndex].getMultipole();
        double[] ash2 = aspMultipoleTypes[AspStates.ASH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = (tautomerLambda * ash1[i] + (1 - tautomerLambda) * ash2[i]) - asp[i];
        }
        break;
      case GLD:
        double[] glu =  gluMultipoleTypes[GluStates.GLU.ordinal()][atomIndex].getMultipole();
        double[] glh1 = gluMultipoleTypes[GluStates.GLH1.ordinal()][atomIndex].getMultipole();
        double[] glh2 = gluMultipoleTypes[GluStates.GLH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = (tautomerLambda * glh1[i] + (1 - tautomerLambda) * glh2[i]) - glu[i];
        }
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
        double[] his = hisMultipoleTypes[HisStates.HIS.ordinal()][atomIndex].getMultipole();
        double[] hid = hisMultipoleTypes[HisStates.HID.ordinal()][atomIndex].getMultipole();
        double[] hie = hisMultipoleTypes[HisStates.HIE.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * his[i] + (1.0 - titrationLambda) * (hid[i] - hie[i]);
        }
        break;
      case ASD:
        double[] asp =  aspMultipoleTypes[AspStates.ASP.ordinal()][atomIndex].getMultipole();
        double[] ash1 = aspMultipoleTypes[AspStates.ASH1.ordinal()][atomIndex].getMultipole();
        double[] ash2 = aspMultipoleTypes[AspStates.ASH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (ash1[i] - ash2[i]) + (1.0 - titrationLambda) * asp[i];
        }
        break;
      case GLD:
        double[] glu =  gluMultipoleTypes[GluStates.GLU.ordinal()][atomIndex].getMultipole();
        double[] glh1 = gluMultipoleTypes[GluStates.GLH1.ordinal()][atomIndex].getMultipole();
        double[] glh2 = gluMultipoleTypes[GluStates.GLH2.ordinal()][atomIndex].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (glh1[i] - glh2[i]) + (1.0 - titrationLambda) * glu[i];
        }
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
        double his = hisPolarizeTypes[HisStates.HIS.ordinal()][atomIndex].polarizability;
        double hid = hisPolarizeTypes[HisStates.HID.ordinal()][atomIndex].polarizability;
        double hie = hisPolarizeTypes[HisStates.HIE.ordinal()][atomIndex].polarizability;
        return titrationLambda * his + (1.0 - titrationLambda) * (tautomerLambda * hid + (1 - tautomerLambda) * hie);
      case ASD:
        double asp =  aspPolarizeTypes[AspStates.ASP.ordinal()][atomIndex].polarizability;
        double ash1 = aspPolarizeTypes[AspStates.ASH1.ordinal()][atomIndex].polarizability;
        double ash2 = aspPolarizeTypes[AspStates.ASH2.ordinal()][atomIndex].polarizability;
        return titrationLambda * (tautomerLambda * ash1 + (1 - tautomerLambda) * ash2) + (1.0 - titrationLambda) * asp;
      case GLD:
        double glu =  gluPolarizeTypes[GluStates.GLU.ordinal()][atomIndex].polarizability;
        double glh1 = gluPolarizeTypes[GluStates.GLH1.ordinal()][atomIndex].polarizability;
        double glh2 = gluPolarizeTypes[GluStates.GLH2.ordinal()][atomIndex].polarizability;
        return titrationLambda * (tautomerLambda * glh1 + (1 - tautomerLambda) * glh2) + (1.0 - titrationLambda) * glu;
      default:
        return defaultPolarizability;
    }
  }

  public double getPolarizabilityTitrationDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda) {
    switch (aminoAcid3) {
      case LYS:
        double lys = lysPolarizeTypes[LysStates.LYS.ordinal()][atomIndex].polarizability;
        double lyd = lysPolarizeTypes[LysStates.LYD.ordinal()][atomIndex].polarizability;
        return lys - lyd;
      case HIS:
        double his = hisPolarizeTypes[HisStates.HIS.ordinal()][atomIndex].polarizability;
        double hid = hisPolarizeTypes[HisStates.HID.ordinal()][atomIndex].polarizability;
        double hie = hisPolarizeTypes[HisStates.HIE.ordinal()][atomIndex].polarizability;
        return his - (tautomerLambda * hid + (1 - tautomerLambda) * hie);
      case ASD:
        double asp =  aspPolarizeTypes[AspStates.ASP.ordinal()][atomIndex].polarizability;
        double ash1 = aspPolarizeTypes[AspStates.ASH1.ordinal()][atomIndex].polarizability;
        double ash2 = aspPolarizeTypes[AspStates.ASH2.ordinal()][atomIndex].polarizability;
        return (tautomerLambda * ash1 + (1 - tautomerLambda) * ash2) - asp;
      case GLD:
        double glu =  gluPolarizeTypes[GluStates.GLU.ordinal()][atomIndex].polarizability;
        double glh1 = gluPolarizeTypes[GluStates.GLH1.ordinal()][atomIndex].polarizability;
        double glh2 = gluPolarizeTypes[GluStates.GLH2.ordinal()][atomIndex].polarizability;
        return (tautomerLambda * glh1 + (1 - tautomerLambda) * glh2) - glu;
      default:
        return 0.0;
    }
  }

  public double getPolarizabilityTautomerDeriv(AminoAcid3 aminoAcid3, int atomIndex,
      double titrationLambda, double tautomerLambda) {
    switch (aminoAcid3) {
      case HIS:
        double his = hisPolarizeTypes[HisStates.HIS.ordinal()][atomIndex].polarizability;
        double hid = hisPolarizeTypes[HisStates.HID.ordinal()][atomIndex].polarizability;
        double hie = hisPolarizeTypes[HisStates.HIE.ordinal()][atomIndex].polarizability;
        return titrationLambda * his + (1.0 - titrationLambda) * (hid - hie);
      case ASD:
        double asp =  aspPolarizeTypes[AspStates.ASP.ordinal()][atomIndex].polarizability;
        double ash1 = aspPolarizeTypes[AspStates.ASH1.ordinal()][atomIndex].polarizability;
        double ash2 = aspPolarizeTypes[AspStates.ASH2.ordinal()][atomIndex].polarizability;
        return titrationLambda * (ash1 - ash2) + (1.0 - titrationLambda) * asp;
      case GLD:
        double glu =  gluPolarizeTypes[GluStates.GLU.ordinal()][atomIndex].polarizability;
        double glh1 = gluPolarizeTypes[GluStates.GLH1.ordinal()][atomIndex].polarizability;
        double glh2 = gluPolarizeTypes[GluStates.GLH2.ordinal()][atomIndex].polarizability;
        return titrationLambda * (glh1 - glh2) + (1.0 - titrationLambda) * glu;
      case LYS: // No tautomers for LYS.
      default:
        return 0.0;
    }
  }

  private void constructHISState(int biotypeCB, HisStates hisState) {
    int state = hisState.ordinal();
    for (HistidineAtomNames atomName : HistidineAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffsetHIS(hisState);
      if (offset < 0) {
        hisAtomTypes[state][index] = deprotonatedAtomType;
        // Zero out the MultipoleType and Polarizetype.
        hisMultipoleTypes[state][index] = zeroMultipoleType;
        hisPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        hisAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = hisAtomTypes[state][index].getKey();
        hisMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        hisPolarizeTypes[state][index] = forceField.getPolarizeType(key);
        if (hisMultipoleTypes[state][index] == null || hisPolarizeTypes[state][index] == null) {
          logger.severe(format(" A multipole could not be assigned for Lys atom %s.\n %s\n",
              atomName, hisAtomTypes[state][index]));
        }
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
        lysAtomTypes[state][index] = deprotonatedAtomType;
        // Zero out the MultipoleType and Polarizetype.
        lysMultipoleTypes[state][index] = zeroMultipoleType;
        lysPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        lysAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = lysAtomTypes[state][index].getKey();
        lysMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        lysPolarizeTypes[state][index] = forceField.getPolarizeType(key);
        if (lysMultipoleTypes[state][index] == null || lysPolarizeTypes[state][index] == null) {
          logger.severe(format(" A multipole could not be assigned for Lys atom %s.\n %s\n",
              atomName, lysAtomTypes[state][index]));
        }
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
        aspAtomTypes[state][index] = deprotonatedAtomType;
        // Zero out the MultipoleType and Polarizetype.
        aspMultipoleTypes[state][index] = zeroMultipoleType;
        aspPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        aspAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = aspAtomTypes[state][index].getKey();
        aspMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        aspPolarizeTypes[state][index] = forceField.getPolarizeType(key);
        if (aspMultipoleTypes[state][index] == null || aspPolarizeTypes[state][index] == null) {
          logger.severe(format(" A multipole could not be assigned for Asp atom %s.\n %s\n",
              atomName, aspAtomTypes[state][index]));
        }
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
        gluAtomTypes[state][index] = deprotonatedAtomType;
        // Zero out the MultipoleType and Polarizetype.
        gluMultipoleTypes[state][index] = zeroMultipoleType;
        gluPolarizeTypes[state][index] = zeroPolarizeType;
      } else {
        int biotype = biotypeCB + offset;
        gluAtomTypes[state][index] = findAtomType(biotype, forceField);
        String key = gluAtomTypes[state][index].getKey();
        gluMultipoleTypes[state][index] = forceField.getMultipoleTypeBeginsWith(key);
        gluPolarizeTypes[state][index] = forceField.getPolarizeType(key);
        if (gluMultipoleTypes[state][index] == null || gluPolarizeTypes[state][index] == null) {
          logger.severe(format(" A multipole could not be assigned for Glu atom %s.\n %s\n",
              atomName, gluAtomTypes[state][index]));
        }
      }
    }
  }

  private void checkMultipoleFrames(String label, AtomType[][] atomTypes,
      PolarizeType[][] polarizeTypes, MultipoleType[][] multipoleTypes) {
    int states = multipoleTypes.length;
    int types = multipoleTypes[0].length;
    StringBuilder sb = new StringBuilder();
    for (int t = 0; t < types; t++) {
      MultipoleFrameDefinition frame0 = multipoleTypes[0][t].frameDefinition;
      sb.append(format("\n %s Type %d\n", label, t));
      sb.append(
          format(" %s\n  %s\n  %s\n", atomTypes[0][t], polarizeTypes[0][t], multipoleTypes[0][t]));
      for (int s = 1; s < states; s++) {
        sb.append(
            format(" %s\n  %s\n  %s\n", atomTypes[s][t], polarizeTypes[s][t], multipoleTypes[s][t]));
        MultipoleFrameDefinition frame = multipoleTypes[s][t].frameDefinition;
        if (!frame0.equals(frame)) {
          StringBuilder sb2 = new StringBuilder("\n Incompatible atom types:\n");
          sb2.append(format(" %s\n  %s\n  %s\n", atomTypes[0][t], polarizeTypes[0][t],
              multipoleTypes[0][t]));
          sb2.append(format(" %s\n  %s\n  %s\n", atomTypes[s][t], polarizeTypes[s][t],
              multipoleTypes[s][t]));
          logger.severe(sb2.toString());
        }
      }
    }

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(sb.toString());
    }
  }

}
