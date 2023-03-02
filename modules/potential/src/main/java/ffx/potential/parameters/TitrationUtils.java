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
package ffx.potential.parameters;

import static ffx.potential.bonded.AminoAcidUtils.AA_CB;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.ASH;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.ASP;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.CYD;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.CYS;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.GLH;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.GLU;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.HID;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.HIE;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.HIS;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.LYD;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3.LYS;
import static ffx.potential.bonded.BondedUtils.findAtomType;
import static ffx.potential.parameters.MultipoleType.assignAxisAtoms;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.log;

import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.SoluteType.SOLUTE_RADII_TYPE;
import ffx.utilities.Constants;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utilities for interpolating between Amino Acid protonation and tautomer states.
 *
 * @author Michael Schnieders
 * @author Andrew Thiel
 * @since 1.0
 */
public class TitrationUtils {

  private static final Logger logger = Logger.getLogger(TitrationUtils.class.getName());

  private static final double LOG10 = log(10.0);

  private static final MultipoleType aspZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 140, 139}, MultipoleFrameDefinition.ZTHENX, false);
  private static final MultipoleType ashZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 144, 143}, MultipoleFrameDefinition.ZTHENX, false);

  private static final MultipoleType gluZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 158, 157}, MultipoleFrameDefinition.ZTHENX, false);
  private static final MultipoleType glhZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 164, 163}, MultipoleFrameDefinition.ZTHENX, false);

  private static final MultipoleType hieZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 130, 129}, MultipoleFrameDefinition.ZTHENX, false);
  private static final MultipoleType hidZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 126, 124}, MultipoleFrameDefinition.ZTHENX, false);

  private static final MultipoleType lydZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 200, 198}, MultipoleFrameDefinition.ZTHENX, false);

  private static final MultipoleType cydZeroMultipoleType = new MultipoleType(MultipoleType.zeroM,
      new int[] {0, 49, 43}, MultipoleFrameDefinition.ZTHENX, false);

  private static final PolarizeType zeroPolarizeType = new PolarizeType(0, 0.0, 0.39, 0.0,
      new int[] {0});

  private static final SoluteType zeroSoluteType = new SoluteType(0, 1.0);

  private static final AtomType dummyHydrogenAtomType = new AtomType(0, 0, "H", "\"Dummy Hydrogen\"",
      1, 1.0080, 1);

  private static final BondType zeroBondType = new BondType(new int[] {0, 0}, 0.0, 1.0);
  private static final AngleType zeroAngleType = new AngleType(new int[] {0, 0, 0}, 0.0,
      new double[] {0.0});
  private static final StretchBendType zeroStretchBendType = new StretchBendType(new int[] {0, 0, 0},
      new double[] {0.0, 0.0});
  private static final OutOfPlaneBendType zeroOutOfPlaneBendType = new OutOfPlaneBendType(
      new int[] {0, 0, 0, 0}, 0.0);
  private static final TorsionType zeroTorsionType = new TorsionType(new int[] {0, 0, 0, 0},
      new double[] {0.0}, new double[] {0.0}, new int[] {0});
  private static final PiOrbitalTorsionType zeroPiOrbitalTorsionType = new PiOrbitalTorsionType(
      new int[] {0, 0}, 0.0);

  enum AspStates {
    ASP, ASH1, ASH2
  }

  /** Constant <code>AspartateAtomNames</code> */
  private enum AspartateAtomNames {
    CB(0, 0, 0, 0), HB2(1, 1, 1, 0), HB3(1, 1, 1, 0), CG(2, 2, 2, 0), OD1(3, 4, 3, 0), OD2(3, 3, 4,
        0), HD1(-1, 5, -1, 1), HD2(-1, -1, 5, -1);

    /**
     * Biotype offset relative to the CB biotype for charged aspartate (ASP).
     */
    private final int offsetASP;

    /**
     * Biotype offset relative to the CB biotype for neutral aspartic acid protonated on OD1 (ASH1).
     * <p>
     * This is set to negative -1 for the OD2 hydrogen.
     */
    private final int offsetASH1;

    /**
     * Biotype offset relative to the CB biotype for neutral aspartic acid protonated on OD2 (ASH2).
     * <p>
     * This is set to negative -1 for the OD1 hydrogen.
     */
    private final int offsetASH2;

    private final int tautomerDirection;

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
    AspartateAtomNames(int offsetASP, int offsetASH1, int offsetASH2, int tautomerDirection) {
      this.offsetASP = offsetASP;
      this.offsetASH1 = offsetASH1;
      this.offsetASH2 = offsetASH2;
      this.tautomerDirection = tautomerDirection;
    }
  }

  enum GluStates {
    GLU, GLH1, GLH2
  }

  /** Constant <code>GlutamateAtomNames</code> */
  private enum GlutamateAtomNames {
    CB(0, 0, 0, 0), HB2(1, 1, 1, 0), HB3(1, 1, 1, 0), CG(2, 2, 2, 0), HG2(3, 3, 3, 0), HG3(3, 3, 3,
        0), CD(4, 4, 4, 0), OE1(5, 6, 5, 0), OE2(5, 5, 6, 0), HE1(-1, 7, -1, 1), HE2(-1, -1, 7, -1);

    /**
     * Biotype offset relative to the CB biotype for charged Glutamate (GLU).
     */
    private final int offsetGLU;

    /**
     * Biotype offset relative to the CB biotype for neutral Glutamate acid protonated on OE1
     * (GLU1).
     * <p>
     * This is set to negative -1 for the OE2 hydrogen.
     */
    private final int offsetGLH1;

    /**
     * Biotype offset relative to the CB biotype for neutral Glutamate acid protonated on OE2
     * (GLU2).
     * <p>
     * This is set to negative -1 for the OE1 hydrogen.
     */
    private final int offsetGLH2;

    private final int tautomerDirection;

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
    GlutamateAtomNames(int offsetGLU, int offsetGLH1, int offsetGLH2, int tautomerDirection) {
      this.offsetGLU = offsetGLU;
      this.offsetGLH1 = offsetGLH1;
      this.offsetGLH2 = offsetGLH2;
      this.tautomerDirection = tautomerDirection;
    }
  }

  enum LysStates {
    LYD, LYS
  }

  /** Constant <code>lysineAtoms</code> */
  public enum LysineAtomNames {
    CB(0, 0), HB2(1, 1), HB3(1, 1), CG(2, 2), HG2(3, 3), HG3(3, 3), CD(4, 4), HD2(5, 5), HD3(5,
        5), CE(6, 6), HE2(7, 7), HE3(7, 7), NZ(8, 8), HZ1(9, 9), HZ2(9, 9), HZ3(9, -1);

    /**
     * Biotype offset relative to the CB biotype for LYS.
     */
    private final int offsetLYS;

    /**
     * Biotype offset relative to the CB biotype for LYD.
     */
    private final int offsetLYD;

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
    CB(0, 0, 0, 0), HB2(1, 1, 1, 0), HB3(1, 1, 1, 0), CG(2, 2, 2, 0), ND1(3, 3, 3,
        0), // No HD1 proton for HIE; HIE HD1 offset is -1.
    HD1(4, 4, -1, -1), CD2(5, 5, 4, 0), HD2(6, 6, 5, 0), CE1(7, 7, 6, 0), HE1(8, 8, 7, 0), NE2(9, 9,
        8, 0), // No HE2 proton for HID; HID HE2 offset is -1
    HE2(10, -1, 9, 1);

    /**
     * Biotype offset relative to the CB biotype for charged histidine (HIS).
     */
    private final int offsetHIS;

    /**
     * Biotype offset relative to the CB biotype for neutral histidine protonated on the delta
     * nitrogren (HID).
     * <p>
     * This is set to negative -1 for the epsilon hydrogen.
     */
    private final int offsetHID;

    /**
     * Biotype offset relative to the CB biotype for neutral histidine protonated the epsilon
     * nitrogen (HIE).
     * <p>
     * This is set to negative -1 for the delta hydrogen.
     */
    private final int offsetHIE;

    private final int tautomerDirection;

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
    HistidineAtomNames(int offsetHIS, int offsetHID, int offsetHIE, int tautomerDirection) {
      this.offsetHIS = offsetHIS;
      this.offsetHID = offsetHID;
      this.offsetHIE = offsetHIE;
      this.tautomerDirection = tautomerDirection;
    }
  }

  enum CysStates {
    CYS, CYD
  }

  /** Constant <code>CysteineAtoms</code> */
  public enum CysteineAtomNames {
    CB(0, 0), HB2(1, 1), HB3(1, 1), SG(2, 2), HG(3, -1);

    /**
     * Biotype offset relative to the CB biotype for neutral cysteine (CYS).
     */
    private final int offsetCYS;

    /**
     * Biotype offset relative to the CB biotype for negatively charged cysteine (CYD).
     * <p>
     * This is set to negative -1 for the gamma hydrogen.
     */
    private final int offsetCYD;

    public int getOffsetCYS(CysStates state) {
      if (state == CysStates.CYS) {
        return offsetCYS;
      } else {
        return offsetCYD;
      }
    }

    /**
     * Init the Cysteine atom names.
     *
     * @param offsetCYS Biotype relative to the CB biotype for CYS.
     * @param offsetCYD Biotype relative to the CB biotype for CYD.
     */
    CysteineAtomNames(int offsetCYS, int offsetCYD) {
      this.offsetCYS = offsetCYS;
      this.offsetCYD = offsetCYD;
    }
  }

  /**
   * Lysine atom types.
   */
  private final int nLysAtomNames = LysineAtomNames.values().length;
  private final int nLysStates = LysStates.values().length;
  private final AtomType[][] lysAtomTypes = new AtomType[nLysAtomNames][nLysStates];
  private final MultipoleType[][] lysMultipoleTypes = new MultipoleType[nLysAtomNames][nLysStates];
  private final PolarizeType[][] lysPolarizeTypes = new PolarizeType[nLysAtomNames][nLysStates];
  private final VDWType[][] lysVDWTypes = new VDWType[nLysAtomNames][nLysStates];
  private final SoluteType[][] lysSoluteTypes = new SoluteType[nLysAtomNames][nLysStates];

  /**
   * Histidine atom types.
   */
  private final int nHisAtomNames = HistidineAtomNames.values().length;
  private final int nHisStates = HisStates.values().length;
  private final AtomType[][] hisAtomTypes = new AtomType[nHisAtomNames][nHisStates];
  private final MultipoleType[][] hisMultipoleTypes = new MultipoleType[nHisAtomNames][nHisStates];
  private final PolarizeType[][] hisPolarizeTypes = new PolarizeType[nHisAtomNames][nHisStates];
  private final VDWType[][] hisVDWTypes = new VDWType[nHisAtomNames][nHisStates];
  private final SoluteType[][] hisSoluteTypes = new SoluteType[nHisAtomNames][nHisStates];

  /**
   * Aspartic acid atom types.
   */
  private final int nAspAtomNames = AspartateAtomNames.values().length;
  private final int nAspStates = AspStates.values().length;
  private final AtomType[][] aspAtomTypes = new AtomType[nAspAtomNames][nAspStates];
  private final MultipoleType[][] aspMultipoleTypes = new MultipoleType[nAspAtomNames][nAspStates];
  private final PolarizeType[][] aspPolarizeTypes = new PolarizeType[nAspAtomNames][nAspStates];
  private final VDWType[][] aspVDWTypes = new VDWType[nAspAtomNames][nAspStates];
  private final SoluteType[][] aspSoluteTypes = new SoluteType[nAspAtomNames][nAspStates];

  /**
   * Glutamic acid atom types.
   */
  private final int nGluAtomNames = GlutamateAtomNames.values().length;
  private final int nGluStates = GluStates.values().length;
  private final AtomType[][] gluAtomTypes = new AtomType[nGluAtomNames][nGluStates];
  private final MultipoleType[][] gluMultipoleTypes = new MultipoleType[nGluAtomNames][nGluStates];
  private final PolarizeType[][] gluPolarizeTypes = new PolarizeType[nGluAtomNames][nGluStates];
  private final VDWType[][] gluVDWTypes = new VDWType[nGluAtomNames][nGluStates];
  private final SoluteType[][] gluSoluteTypes = new SoluteType[nGluAtomNames][nGluStates];

  /**
   * Cystine atom types.
   */
  private final int nCysAtomNames = CysteineAtomNames.values().length;
  private final int nCysStates = CysStates.values().length;
  private final AtomType[][] cysAtomTypes = new AtomType[nCysAtomNames][nCysStates];
  private final MultipoleType[][] cysMultipoleTypes = new MultipoleType[nCysAtomNames][nCysStates];
  private final PolarizeType[][] cysPolarizeTypes = new PolarizeType[nCysAtomNames][nCysStates];
  private final VDWType[][] cysVDWTypes = new VDWType[nCysAtomNames][nCysStates];
  private final SoluteType[][] cysSoluteTypes = new SoluteType[nCysAtomNames][nCysStates];

  private final ForceField forceField;
  private final SOLUTE_RADII_TYPE soluteRadiiType;
  private final boolean updateBondedTerms;

  private final HashMap<AminoAcid3, Double> rotamerPhBiasMap = new HashMap<>();

  public TitrationUtils(ForceField forceField) {
    this.forceField = forceField;

    String gkRadius = forceField.getString("GK_RADIUS", "SOLUTE");
    SOLUTE_RADII_TYPE tempType;
    try {
      tempType = SOLUTE_RADII_TYPE.valueOf(gkRadius.trim().toUpperCase());
    } catch (Exception e) {
      tempType = SOLUTE_RADII_TYPE.SOLUTE;
    }
    soluteRadiiType = tempType;

    updateBondedTerms = forceField.getBoolean("TITRATION_UPDATE_BONDED_TERMS", true);

    // Populate the Lysine types.
    constructLYSState(AA_CB[LYS.ordinal()], LysStates.LYS);
    constructLYSState(AA_CB[LYD.ordinal()], LysStates.LYD);
    checkParameterTypes("LYS", lysAtomTypes, lysPolarizeTypes, lysMultipoleTypes, lysVDWTypes);

    // Populate the Histidine types.
    constructHISState(AA_CB[HIS.ordinal()], HisStates.HIS);
    constructHISState(AA_CB[HID.ordinal()], HisStates.HID);
    constructHISState(AA_CB[HIE.ordinal()], HisStates.HIE);
    checkParameterTypes("HIS", hisAtomTypes, hisPolarizeTypes, hisMultipoleTypes, hisVDWTypes);

    // Populate the Aspartic acid types.
    constructASPState(AA_CB[ASP.ordinal()], AspStates.ASP);
    constructASPState(AA_CB[ASH.ordinal()], AspStates.ASH1); // First ASH Tautomer
    constructASPState(AA_CB[ASH.ordinal()], AspStates.ASH2); // Second ASH Tautomer
    checkParameterTypes("ASP", aspAtomTypes, aspPolarizeTypes, aspMultipoleTypes, aspVDWTypes);

    // Populate the Glutamic acid types.
    constructGLUState(AA_CB[GLU.ordinal()], GluStates.GLU);
    constructGLUState(AA_CB[GLH.ordinal()], GluStates.GLH1); // First GLH Tautomer
    constructGLUState(AA_CB[GLH.ordinal()], GluStates.GLH2); // Second GLH Tautomer
    checkParameterTypes("GLU", gluAtomTypes, gluPolarizeTypes, gluMultipoleTypes, gluVDWTypes);

    // Populate the Cystine types.
    constructCYSState(AA_CB[CYS.ordinal()], CysStates.CYS);
    constructCYSState(AA_CB[CYD.ordinal()], CysStates.CYD);
    checkParameterTypes("CYS", cysAtomTypes, cysPolarizeTypes, cysMultipoleTypes, cysVDWTypes);
  }

  public boolean testResidueTypes(Residue residue) {

    boolean testPassed = true;
    int nStates = 1;
    AminoAcid3 aminoAcid3 = residue.getAminoAcid3();
    switch (aminoAcid3) {
      case ASP, ASH, ASD, GLU, GLH, GLD, HIS, HID, HIE -> nStates = 3;
      case CYS, CYD, LYS, LYD -> nStates = 2;
      default -> logger.info(format(" Only one state for atom %s.", aminoAcid3));
    }

    List<Atom> atomList = residue.getSideChainAtoms();
    int nAtoms = atomList.size();
    int[][][] axisAtomIndices = new int[nAtoms][nStates][];
    AtomType[][] atomTypes = new AtomType[nAtoms][nStates];
    MultipoleType[][] multipoleTypes = new MultipoleType[nAtoms][nStates];

    AtomType[] initialAtomTypes = new AtomType[nAtoms];
    MultipoleType[] initialMultipoleTypes = new MultipoleType[nAtoms];
    // Store initial state
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atomList.get(i);
      initialAtomTypes[i] = atom.getAtomType();
      initialMultipoleTypes[i] = atom.getMultipoleType();
    }

    // Load information for each state.
    for (int state = 0; state < nStates; state++) {
      // Load AtomType and MultipoleType instances for each atom for this state.
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atomList.get(i);
        String atomName = atom.getName();
        switch (aminoAcid3) {
          case ASP:
          case ASH:
          case ASD:
            int index = AspartateAtomNames.valueOf(atomName).ordinal();
            atom.setAtomType(aspAtomTypes[index][state]);
            atom.setMultipoleType(aspMultipoleTypes[index][state]);
            break;
          case CYS:
          case CYD:
            index = CysteineAtomNames.valueOf(atomName).ordinal();
            atom.setAtomType(cysAtomTypes[index][state]);
            atom.setMultipoleType(cysMultipoleTypes[index][state]);
            break;
          case GLU:
          case GLH:
          case GLD:
            index = GlutamateAtomNames.valueOf(atomName).ordinal();
            atom.setAtomType(gluAtomTypes[index][state]);
            atom.setMultipoleType(gluMultipoleTypes[index][state]);
            break;
          case HIS:
          case HID:
          case HIE:
            index = HistidineAtomNames.valueOf(atomName).ordinal();
            atom.setAtomType(hisAtomTypes[index][state]);
            atom.setMultipoleType(hisMultipoleTypes[index][state]);
            break;
          case LYS:
          case LYD:
            index = LysineAtomNames.valueOf(atomName).ordinal();
            atom.setAtomType(lysAtomTypes[index][state]);
            atom.setMultipoleType(lysMultipoleTypes[index][state]);
            break;
          default:
            logger.info(format(" Only one state for atom %s.", atom));
        }
        atomTypes[i][state] = atom.getAtomType();
        multipoleTypes[i][state] = atom.getMultipoleType();
      }
      // Assign axis atoms for each atom for this state.
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atomList.get(i);
        assignAxisAtoms(atom);
        axisAtomIndices[i][state] = atom.getAxisAtomIndices();
      }
    }

    // Check the local multipole frames.
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atomList.get(i);
      int[] referenceIndices = axisAtomIndices[i][0];
      AtomType referenceAtomType = atomTypes[i][0];
      MultipoleType referenceMultipoleType = multipoleTypes[i][0];
      for (int state = 1; state < nStates; state++) {
        int[] stateIndices = axisAtomIndices[i][state];
        AtomType stateAtomType = atomTypes[i][state];
        MultipoleType stateMultipoleType = multipoleTypes[i][state];
        if (referenceMultipoleType.frameDefinition != stateMultipoleType.frameDefinition) {
          logger.warning(format(" Local frame definition is inconsistent for atom %s", atom));
          logger.warning(format("  %s\n  %s", referenceAtomType, referenceMultipoleType));
          logger.warning(format("  %s\n  %s", stateAtomType, stateMultipoleType));
          testPassed = false;
          continue;
        }
        if (Arrays.compare(referenceIndices, stateIndices) != 0) {
          // Atom order does not matter for BISECTOR.
          if (referenceMultipoleType.frameDefinition == MultipoleFrameDefinition.BISECTOR) {
            if (referenceIndices[0] == stateIndices[1] && referenceIndices[1] == stateIndices[0]) {
              continue;
            }
          }
          logger.warning(format(" Local frame atom indices are inconsistent for atom %s", atom));
          logger.warning(
              format("  %s %s\n  %s", referenceAtomType, Arrays.toString(referenceIndices),
                  referenceMultipoleType));
          logger.warning(format("  %s %s\n  %s", stateAtomType, Arrays.toString(stateIndices),
              stateMultipoleType));
          testPassed = false;
        }
      }
    }

    // Revert initial state
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atomList.get(i);
      atom.setAtomType(initialAtomTypes[i]);
      atom.setMultipoleType(initialMultipoleTypes[i]);
    }

    return testPassed;
  }

  /**
   * Update force field parameters for the side-chain atoms of the given residue based on the rotamer
   * amino acid type.
   *
   * @param residue Residue to update.
   * @param rotamer Rotamer that contains the amino acid residue identity.
   */
  public void updateResidueParameters(Residue residue, Rotamer rotamer) {
    if (!rotamer.isTitrating) {
      return;
    }

    AminoAcid3 aminoAcid3 = residue.getAminoAcid3();
    switch (aminoAcid3) {
      case ASH, ASP -> {
        // Assume ASP types
        int aspIndex = AspStates.ASP.ordinal();
        if (rotamer.aminoAcid3 == ASH) {
          // Use ASH2 types
          aspIndex = AspStates.ASH2.ordinal();
        }
        for (AspartateAtomNames atomName : AspartateAtomNames.values()) {
          if (atomName.name().equals("HD1")) {
            // Skip the HD1 atom name (used only for ASD during constant pH).
            // This atom should not be present in the residue for ASH/ASP rot opt.
            continue;
          }
          int atomIndex = atomName.ordinal();
          Atom atom = (Atom) residue.getAtomNode(atomName.name());
          if (atom == null) {
            logger.severe(" Atom is null for " + atomName);
            return;
          }
          atom.setAtomType(aspAtomTypes[atomIndex][aspIndex]);
          atom.setMultipoleType(aspMultipoleTypes[atomIndex][aspIndex]);
          atom.setPolarizeType(aspPolarizeTypes[atomIndex][aspIndex]);
          atom.setVDWType(aspVDWTypes[atomIndex][aspIndex]);
          atom.setSoluteType(aspSoluteTypes[atomIndex][aspIndex]);
        }
      }
      case GLU, GLH -> {
        // Assume GLU types
        int gluIndex = GluStates.GLU.ordinal();
        if (rotamer.aminoAcid3 == GLH) {
          // Use GLH2 types
          gluIndex = GluStates.GLH2.ordinal();
        }
        for (GlutamateAtomNames atomName : GlutamateAtomNames.values()) {
          if (atomName.name().equals("HE1")) {
            // Skip the HE1 atom name (used only for GLD during constant pH).
            // This atom should not be present in the residue for GLH/GLU rot opt.
            continue;
          }
          int atomIndex = atomName.ordinal();
          Atom atom = (Atom) residue.getAtomNode(atomName.name());
          if (atom == null) {
            logger.severe(" Atom is null for " + atomName);
            return;
          }
          atom.setAtomType(gluAtomTypes[atomIndex][gluIndex]);
          atom.setMultipoleType(gluMultipoleTypes[atomIndex][gluIndex]);
          atom.setPolarizeType(gluPolarizeTypes[atomIndex][gluIndex]);
          atom.setVDWType(gluVDWTypes[atomIndex][gluIndex]);
          atom.setSoluteType(gluSoluteTypes[atomIndex][gluIndex]);
        }
      }
      case LYS, LYD -> {
        // Assume LYS types
        int lysIndex = LysStates.LYS.ordinal();
        if (rotamer.aminoAcid3 == LYD) {
          // Use LYD types
          lysIndex = LysStates.LYD.ordinal();
        }
        for (LysineAtomNames atomName : LysineAtomNames.values()) {
          int atomIndex = atomName.ordinal();
          Atom atom = (Atom) residue.getAtomNode(atomName.name());
          if (atom == null) {
            logger.severe(" Atom is null for " + atomName);
            return;
          }
          atom.setAtomType(lysAtomTypes[atomIndex][lysIndex]);
          atom.setMultipoleType(lysMultipoleTypes[atomIndex][lysIndex]);
          atom.setPolarizeType(lysPolarizeTypes[atomIndex][lysIndex]);
          atom.setVDWType(lysVDWTypes[atomIndex][lysIndex]);
          atom.setSoluteType(lysSoluteTypes[atomIndex][lysIndex]);
        }
      }
      case CYS, CYD -> {
        // Assume CYS types
        int cysIndex = CysStates.CYS.ordinal();
        if (rotamer.aminoAcid3 == CYD) {
          // Use CYD types
          cysIndex = CysStates.CYD.ordinal();
        }
        for (CysteineAtomNames atomName : CysteineAtomNames.values()) {
          int atomIndex = atomName.ordinal();
          Atom atom = (Atom) residue.getAtomNode(atomName.name());
          if (atom == null) {
            logger.severe(" Atom is null for " + atomName);
            return;
          }
          atom.setAtomType(cysAtomTypes[atomIndex][cysIndex]);
          atom.setMultipoleType(cysMultipoleTypes[atomIndex][cysIndex]);
          atom.setPolarizeType(cysPolarizeTypes[atomIndex][cysIndex]);
          atom.setVDWType(cysVDWTypes[atomIndex][cysIndex]);
          atom.setSoluteType(cysSoluteTypes[atomIndex][cysIndex]);
        }
      }
      case HIS, HIE, HID -> {
        // Assume HIS types.
        int hisIndex = switch (rotamer.aminoAcid3) {
          case HIE -> HisStates.HIE.ordinal();
          case HID -> HisStates.HID.ordinal();
          default -> HisStates.HIS.ordinal();
        };
        for (HistidineAtomNames atomName : HistidineAtomNames.values()) {
          int atomIndex = atomName.ordinal();
          Atom atom = (Atom) residue.getAtomNode(atomName.name());
          if (atom == null) {
            logger.severe(" Atom is null for " + atomName);
            return;
          }
          atom.setAtomType(hisAtomTypes[atomIndex][hisIndex]);
          atom.setMultipoleType(hisMultipoleTypes[atomIndex][hisIndex]);
          atom.setPolarizeType(hisPolarizeTypes[atomIndex][hisIndex]);
          atom.setVDWType(hisVDWTypes[atomIndex][hisIndex]);
          atom.setSoluteType(hisSoluteTypes[atomIndex][hisIndex]);
        }
      }
      default -> logger.severe(
          format(" No support for titrating residue %s with rotamer %s.", residue, rotamer));
    }

    // Update local frame defining atoms now that the AtomType and MultipoleType values are set.
    for (Atom atom : residue.getSideChainAtoms()) {
      assignAxisAtoms(atom);
    }

    // Should bonded terms be updated.
    if (!updateBondedTerms) {
      return;
    }

    // Update Bond force field terms.
    for (Bond bond : residue.getBondList()) {
      AtomType a1 = bond.getAtom(0).getAtomType();
      AtomType a2 = bond.getAtom(1).getAtomType();
      BondType bondType = forceField.getBondType(a1, a2);
      if (bondType == null) {
        bondType = zeroBondType;
      }
      bond.setBondType(bondType);
    }

    // Update Angle force field terms.
    for (Angle angle : residue.getAngleList()) {
      AtomType a1 = angle.getAtom(0).getAtomType();
      AtomType a2 = angle.getAtom(1).getAtomType();
      AtomType a3 = angle.getAtom(2).getAtomType();
      AngleType angleType = forceField.getAngleType(a1, a2, a3);
      if (angleType == null) {
        angleType = zeroAngleType;
      }
      angle.setAngleType(angleType);
    }

    // Update Stretch-Bend force field terms.
    for (StretchBend stretchBend : residue.getStretchBendList()) {
      AtomType a1 = stretchBend.getAtom(0).getAtomType();
      AtomType a2 = stretchBend.getAtom(1).getAtomType();
      AtomType a3 = stretchBend.getAtom(2).getAtomType();
      StretchBendType stretchBendType = forceField.getStretchBendType(a1, a2, a3);
      if (stretchBendType == null) {
        stretchBendType = zeroStretchBendType;
      }
      stretchBend.setStretchBendType(stretchBendType);
    }

    // Update Out-of-Plane Bend force field terms.
    for (OutOfPlaneBend outOfPlaneBend : residue.getOutOfPlaneBendList()) {
      AtomType a4 = outOfPlaneBend.getFourthAtom().getAtomType();
      AtomType a0 = outOfPlaneBend.getFirstAngleAtom().getAtomType();
      AtomType a1 = outOfPlaneBend.getTrigonalAtom().getAtomType();
      AtomType a2 = outOfPlaneBend.getLastAngleAtom().getAtomType();
      OutOfPlaneBendType outOfPlaneBendType = forceField.getOutOfPlaneBendType(a4, a0, a1, a2);
      if (outOfPlaneBendType == null) {
        outOfPlaneBendType = zeroOutOfPlaneBendType;
      }
      outOfPlaneBend.setOutOfPlaneBendType(outOfPlaneBendType);
    }

    // Update torsion force field terms.
    for (Torsion torsion : residue.getTorsionList()) {
      AtomType a1 = torsion.getAtom(0).getAtomType();
      AtomType a2 = torsion.getAtom(1).getAtomType();
      AtomType a3 = torsion.getAtom(2).getAtomType();
      AtomType a4 = torsion.getAtom(3).getAtomType();
      TorsionType torsionType = forceField.getTorsionType(a1, a2, a3, a4);
      if (torsionType == null) {
        torsionType = zeroTorsionType;
      }
      torsion.setTorsionType(torsionType);
    }

    // Update Pi-Orbital Torsion force field terms.
    for (PiOrbitalTorsion piOrbitalTorsion : residue.getPiOrbitalTorsionList()) {
      Bond middleBond = piOrbitalTorsion.getMiddleBond();
      AtomType a1 = middleBond.getAtom(0).getAtomType();
      AtomType a2 = middleBond.getAtom(1).getAtomType();
      PiOrbitalTorsionType piOrbitalTorsionType = forceField.getPiOrbitalTorsionType(a1, a2);
      if (piOrbitalTorsionType == null) {
        piOrbitalTorsionType = zeroPiOrbitalTorsionType;
      }
      piOrbitalTorsion.setPiOrbitalTorsionType(piOrbitalTorsionType);
    }

    // The following terms are not supported yet.
    List<ImproperTorsion> improperTorsions = residue.getImproperTorsionList();
    if (improperTorsions != null && improperTorsions.size() > 0) {
      logger.severe(
          " Improper torsions are not supported yet for pH-dependent rotamer optimization.");
    }

    List<StretchTorsion> stretchTorsions = residue.getStretchTorsionList();
    if (stretchTorsions != null && stretchTorsions.size() > 0) {
      logger.severe(
          " Stretch-torsions are not supported yet for pH-dependent rotamer optimization.");
    }

    List<AngleTorsion> angleTorsions = residue.getAngleTorsionList();
    if (angleTorsions != null && angleTorsions.size() > 0) {
      logger.severe(" Angle-torsions are not supported yet for pH-dependent rotamer optimization.");
    }

    List<TorsionTorsion> torsionTorsions = residue.getTorsionTorsionList();
    if (torsionTorsions != null && torsionTorsions.size() > 0) {
      logger.severe(
          " Torsion-torsions are not supported yet for pH-dependent rotamer optimization.");
    }

    List<UreyBradley> ureyBradleys = residue.getUreyBradleyList();
    if (ureyBradleys != null && ureyBradleys.size() > 0) {
      logger.severe(" Urey-Bradleys are not supported yet for pH-dependent rotamer optimization.");
    }

  }

  public double[] getMultipole(Atom atom, double titrationLambda, double tautomerLambda,
      double[] multipole) {
    /*
    Step 1: retrieve the atomName from atom instance.
    Step 2: retrieve the oridnal from the atom instance + residueType
     */

    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return multipole;
    }
    String atomName = atom.getName();

    switch (aminoAcid3) {
      case LYS:
        int atomIndex = LysineAtomNames.valueOf(atomName).ordinal();
        MultipoleType lysM = lysMultipoleTypes[atomIndex][LysStates.LYS.ordinal()];
        MultipoleType lydM = lysMultipoleTypes[atomIndex][LysStates.LYD.ordinal()];
        double[] lys = lysM.getMultipole();
        double[] lyd = lydM.getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * lys[i] + (1.0 - titrationLambda) * lyd[i];
        }
        break;
      case CYS:
        atomIndex = CysteineAtomNames.valueOf(atomName).ordinal();
        MultipoleType cysM = cysMultipoleTypes[atomIndex][CysStates.CYS.ordinal()];
        MultipoleType cydM = cysMultipoleTypes[atomIndex][CysStates.CYD.ordinal()];
        double[] cys = cysM.getMultipole();
        double[] cyd = cydM.getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * cys[i] + (1.0 - titrationLambda) * cyd[i];
        }
        break;
      case HIS:
        atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        MultipoleType hisM = hisMultipoleTypes[atomIndex][HisStates.HIS.ordinal()];
        MultipoleType hidM = hisMultipoleTypes[atomIndex][HisStates.HID.ordinal()];
        MultipoleType hieM = hisMultipoleTypes[atomIndex][HisStates.HIE.ordinal()];
        double[] his = hisM.getMultipole();
        double[] hid = hidM.getMultipole();
        double[] hie = hieM.getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] =
              titrationLambda * his[i] + (1.0 - titrationLambda) * (tautomerLambda * hie[i]
                  + (1 - tautomerLambda) * hid[i]);
        }
        break;
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        MultipoleType aspM = aspMultipoleTypes[atomIndex][AspStates.ASP.ordinal()];
        MultipoleType ash1M = aspMultipoleTypes[atomIndex][AspStates.ASH1.ordinal()];
        MultipoleType ash2M = aspMultipoleTypes[atomIndex][AspStates.ASH2.ordinal()];
        double[] asp = aspM.getMultipole();
        double[] ash1 = ash1M.getMultipole();
        double[] ash2 = ash2M.getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] =
              titrationLambda * (tautomerLambda * ash1[i] + (1 - tautomerLambda) * ash2[i])
                  + (1.0 - titrationLambda) * asp[i];
        }
        break;
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        MultipoleType gluM = gluMultipoleTypes[atomIndex][GluStates.GLU.ordinal()];
        MultipoleType glh1M = gluMultipoleTypes[atomIndex][GluStates.GLH1.ordinal()];
        MultipoleType glh2M = gluMultipoleTypes[atomIndex][GluStates.GLH2.ordinal()];
        double[] glu = gluM.getMultipole();
        double[] glh1 = glh1M.getMultipole();
        double[] glh2 = glh2M.getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] =
              titrationLambda * (tautomerLambda * glh1[i] + (1 - tautomerLambda) * glh2[i])
                  + (1.0 - titrationLambda) * glu[i];
        }
        break;
      default:
        return multipole;
    }
    return multipole;
  }

  public double[] getMultipoleTitrationDeriv(Atom atom, double titrationLambda,
      double tautomerLambda, double[] multipole) {
    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return multipole;
    }
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case LYS:
        int atomIndex = LysineAtomNames.valueOf(atomName).ordinal();
        double[] lys = lysMultipoleTypes[atomIndex][LysStates.LYS.ordinal()].getMultipole();
        double[] lyd = lysMultipoleTypes[atomIndex][LysStates.LYD.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = lys[i] - lyd[i];
        }
        break;
      case CYS:
        atomIndex = CysteineAtomNames.valueOf(atomName).ordinal();
        double[] cys = cysMultipoleTypes[atomIndex][CysStates.CYS.ordinal()].getMultipole();
        double[] cyd = cysMultipoleTypes[atomIndex][CysStates.CYD.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = cys[i] - cyd[i];
        }
        break;
      case HIS:
        atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        double[] his = hisMultipoleTypes[atomIndex][HisStates.HIS.ordinal()].getMultipole();
        double[] hid = hisMultipoleTypes[atomIndex][HisStates.HID.ordinal()].getMultipole();
        double[] hie = hisMultipoleTypes[atomIndex][HisStates.HIE.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = his[i] - (tautomerLambda * hie[i] + (1 - tautomerLambda) * hid[i]);
        }
        break;
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        double[] asp = aspMultipoleTypes[atomIndex][AspStates.ASP.ordinal()].getMultipole();
        double[] ash1 = aspMultipoleTypes[atomIndex][AspStates.ASH1.ordinal()].getMultipole();
        double[] ash2 = aspMultipoleTypes[atomIndex][AspStates.ASH2.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = (tautomerLambda * ash1[i] + (1 - tautomerLambda) * ash2[i]) - asp[i];
        }
        break;
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        double[] glu = gluMultipoleTypes[atomIndex][GluStates.GLU.ordinal()].getMultipole();
        double[] glh1 = gluMultipoleTypes[atomIndex][GluStates.GLH1.ordinal()].getMultipole();
        double[] glh2 = gluMultipoleTypes[atomIndex][GluStates.GLH2.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = (tautomerLambda * glh1[i] + (1 - tautomerLambda) * glh2[i]) - glu[i];
        }
        break;
      default:
        return multipole;
    }
    return multipole;
  }

  public double[] getMultipoleTautomerDeriv(Atom atom, double titrationLambda, double tautomerLambda,
      double[] multipole) {
    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return multipole;
    }
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case HIS:
        int atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        double[] his = hisMultipoleTypes[atomIndex][HisStates.HIS.ordinal()].getMultipole();
        double[] hid = hisMultipoleTypes[atomIndex][HisStates.HID.ordinal()].getMultipole();
        double[] hie = hisMultipoleTypes[atomIndex][HisStates.HIE.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = (1.0 - titrationLambda) * (hie[i] - hid[i]);
        }
        break;
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        double[] asp = aspMultipoleTypes[atomIndex][AspStates.ASP.ordinal()].getMultipole();
        double[] ash1 = aspMultipoleTypes[atomIndex][AspStates.ASH1.ordinal()].getMultipole();
        double[] ash2 = aspMultipoleTypes[atomIndex][AspStates.ASH2.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (ash1[i] - ash2[i]);
        }
        break;
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        double[] glu = gluMultipoleTypes[atomIndex][GluStates.GLU.ordinal()].getMultipole();
        double[] glh1 = gluMultipoleTypes[atomIndex][GluStates.GLH1.ordinal()].getMultipole();
        double[] glh2 = gluMultipoleTypes[atomIndex][GluStates.GLH2.ordinal()].getMultipole();
        for (int i = 0; i < multipole.length; i++) {
          multipole[i] = titrationLambda * (glh1[i] - glh2[i]);
        }
        break;
      case LYS: // No tautomers for LYS.
      case CYS: // No tautomers for CYS.
      default:
        return multipole;
    }
    return multipole;
  }

  public double getPolarizability(Atom atom, double titrationLambda, double tautomerLambda,
      double defaultPolarizability) {
    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return defaultPolarizability;
    }
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case LYS:
        int atomIndex = LysineAtomNames.valueOf(atomName).ordinal();
        double lys = lysPolarizeTypes[atomIndex][LysStates.LYS.ordinal()].polarizability;
        double lyd = lysPolarizeTypes[atomIndex][LysStates.LYD.ordinal()].polarizability;
        return titrationLambda * lys + (1.0 - titrationLambda) * lyd;
      case CYS:
        atomIndex = CysteineAtomNames.valueOf(atomName).ordinal();
        double cys = cysPolarizeTypes[atomIndex][CysStates.CYS.ordinal()].polarizability;
        double cyd = cysPolarizeTypes[atomIndex][CysStates.CYD.ordinal()].polarizability;
        return titrationLambda * cys + (1.0 - titrationLambda) * cyd;
      case HIS:
        atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        double his = hisPolarizeTypes[atomIndex][HisStates.HIS.ordinal()].polarizability;
        double hid = hisPolarizeTypes[atomIndex][HisStates.HID.ordinal()].polarizability;
        double hie = hisPolarizeTypes[atomIndex][HisStates.HIE.ordinal()].polarizability;
        return titrationLambda * his + (1.0 - titrationLambda) * (tautomerLambda * hie
            + (1.0 - tautomerLambda) * hid);
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        double asp = aspPolarizeTypes[atomIndex][AspStates.ASP.ordinal()].polarizability;
        double ash1 = aspPolarizeTypes[atomIndex][AspStates.ASH1.ordinal()].polarizability;
        double ash2 = aspPolarizeTypes[atomIndex][AspStates.ASH2.ordinal()].polarizability;
        return titrationLambda * (tautomerLambda * ash1 + (1.0 - tautomerLambda) * ash2)
            + (1.0 - titrationLambda) * asp;
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        double glu = gluPolarizeTypes[atomIndex][GluStates.GLU.ordinal()].polarizability;
        double glh1 = gluPolarizeTypes[atomIndex][GluStates.GLH1.ordinal()].polarizability;
        double glh2 = gluPolarizeTypes[atomIndex][GluStates.GLH2.ordinal()].polarizability;
        return titrationLambda * (tautomerLambda * glh1 + (1.0 - tautomerLambda) * glh2)
            + (1.0 - titrationLambda) * glu;
      default:
        return defaultPolarizability;
    }
  }

  public double getPolarizabilityTitrationDeriv(Atom atom, double titrationLambda,
      double tautomerLambda) {
    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return 0.0;
    }
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case LYS:
        int atomIndex = LysineAtomNames.valueOf(atomName).ordinal();
        double lys = lysPolarizeTypes[atomIndex][LysStates.LYS.ordinal()].polarizability;
        double lyd = lysPolarizeTypes[atomIndex][LysStates.LYD.ordinal()].polarizability;
        return lys - lyd;
      case CYS:
        atomIndex = CysteineAtomNames.valueOf(atomName).ordinal();
        double cys = cysPolarizeTypes[atomIndex][CysStates.CYS.ordinal()].polarizability;
        double cyd = cysPolarizeTypes[atomIndex][CysStates.CYD.ordinal()].polarizability;
        return cys - cyd;
      case HIS:
        atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        double his = hisPolarizeTypes[atomIndex][HisStates.HIS.ordinal()].polarizability;
        double hid = hisPolarizeTypes[atomIndex][HisStates.HID.ordinal()].polarizability;
        double hie = hisPolarizeTypes[atomIndex][HisStates.HIE.ordinal()].polarizability;
        return his - (tautomerLambda * hie + (1.0 - tautomerLambda) * hid);
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        double asp = aspPolarizeTypes[atomIndex][AspStates.ASP.ordinal()].polarizability;
        double ash1 = aspPolarizeTypes[atomIndex][AspStates.ASH1.ordinal()].polarizability;
        double ash2 = aspPolarizeTypes[atomIndex][AspStates.ASH2.ordinal()].polarizability;
        return (tautomerLambda * ash1 + (1.0 - tautomerLambda) * ash2) - asp;
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        double glu = gluPolarizeTypes[atomIndex][GluStates.GLU.ordinal()].polarizability;
        double glh1 = gluPolarizeTypes[atomIndex][GluStates.GLH1.ordinal()].polarizability;
        double glh2 = gluPolarizeTypes[atomIndex][GluStates.GLH2.ordinal()].polarizability;
        return (tautomerLambda * glh1 + (1.0 - tautomerLambda) * glh2) - glu;
      default:
        return 0.0;
    }
  }

  public double getPolarizabilityTautomerDeriv(Atom atom, double titrationLambda,
      double tautomerLambda) {
    AminoAcid3 aminoAcid3;
    try {
      aminoAcid3 = atom.getMSNode(Residue.class).getAminoAcid3();
    } catch (Exception exception) {
      return 0.0;
    }
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case HIS:
        int atomIndex = HistidineAtomNames.valueOf(atomName).ordinal();
        double his = hisPolarizeTypes[atomIndex][HisStates.HIS.ordinal()].polarizability;
        double hid = hisPolarizeTypes[atomIndex][HisStates.HID.ordinal()].polarizability;
        double hie = hisPolarizeTypes[atomIndex][HisStates.HIE.ordinal()].polarizability;
        return (1.0 - titrationLambda) * (hie - hid);
      case ASD:
        atomIndex = AspartateAtomNames.valueOf(atomName).ordinal();
        double asp = aspPolarizeTypes[atomIndex][AspStates.ASP.ordinal()].polarizability;
        double ash1 = aspPolarizeTypes[atomIndex][AspStates.ASH1.ordinal()].polarizability;
        double ash2 = aspPolarizeTypes[atomIndex][AspStates.ASH2.ordinal()].polarizability;
        return titrationLambda * (ash1 - ash2);
      case GLD:
        atomIndex = GlutamateAtomNames.valueOf(atomName).ordinal();
        double glu = gluPolarizeTypes[atomIndex][GluStates.GLU.ordinal()].polarizability;
        double glh1 = gluPolarizeTypes[atomIndex][GluStates.GLH1.ordinal()].polarizability;
        double glh2 = gluPolarizeTypes[atomIndex][GluStates.GLH2.ordinal()].polarizability;
        return titrationLambda * (glh1 - glh2);
      case LYS: // No tautomers for LYS.
      case CYS: // No tautomers for LYS.
      default:
        return 0.0;
    }
  }

  public static boolean isTitratingHydrogen(AminoAcid3 aminoAcid3, Atom atom) {
    boolean isTitratingHydrogen = false;
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case ASD -> {
        if (atomName.equals(AspartateAtomNames.HD1.name()) || atomName.equals(
            AspartateAtomNames.HD2.name())) {
          isTitratingHydrogen = true;
        }
      }
      case GLD -> {
        if (atomName.equals(GlutamateAtomNames.HE1.name()) || atomName.equals(
            GlutamateAtomNames.HE2.name())) {
          isTitratingHydrogen = true;
        }
      }
      case HIS -> {
        if (atomName.equals(HistidineAtomNames.HD1.name()) || atomName.equals(
            HistidineAtomNames.HE2.name())) {
          isTitratingHydrogen = true;
        }
      }
      case LYS -> {
        if (atomName.equals(LysineAtomNames.HZ3.name())) {
          isTitratingHydrogen = true;
        }
      }
      case CYS -> {
        if (atomName.equals(CysteineAtomNames.HG.name())) {
          isTitratingHydrogen = true;
        }
      }
    }
    return isTitratingHydrogen;
  }

  /**
   * Used to keep track of heavy atoms with changing polarizability. Only affects carboxylic oxygen
   * and sulfur.
   *
   * @param aminoAcid3 The amino acid type.
   * @param atom The atom to check.
   * @return True if the atom is a heavy atom with changing polarizability.
   */

  public static boolean isTitratingHeavy(AminoAcid3 aminoAcid3, Atom atom) {
    boolean isTitratingHeavy = false;
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case ASD -> {
        if (atomName.equals(AspartateAtomNames.OD1.name()) || atomName.equals(
            AspartateAtomNames.OD2.name())) {
          isTitratingHeavy = true;
        }
      }
      case GLD -> {
        if (atomName.equals(GlutamateAtomNames.OE1.name()) || atomName.equals(
            GlutamateAtomNames.OE2.name())) {
          isTitratingHeavy = true;
        }
      }
      case CYS -> {
        if (atomName.equals(CysteineAtomNames.SG.name())) {
          isTitratingHeavy = true;
        }
      }
    }
    return isTitratingHeavy;
  }

  public static int getTitratingHydrogenDirection(AminoAcid3 aminoAcid3, Atom atom) {
    int tautomerDirection = 0;
    String atomName = atom.getName();
    switch (aminoAcid3) {
      case ASD -> {
        if (atomName.equals(AspartateAtomNames.HD1.name())) {
          tautomerDirection = AspartateAtomNames.HD1.tautomerDirection;
        } else if (atomName.equals(AspartateAtomNames.HD2.name())) {
          tautomerDirection = AspartateAtomNames.HD2.tautomerDirection;
        }
      }
      case GLD -> {
        if (atomName.equals(GlutamateAtomNames.HE1.name())) {
          tautomerDirection = GlutamateAtomNames.HE1.tautomerDirection;
        } else if (atomName.equals(GlutamateAtomNames.HE2.name())) {
          tautomerDirection = GlutamateAtomNames.HE2.tautomerDirection;
        }
      }
      case HIS -> {
        if (atomName.equals(HistidineAtomNames.HD1.name())) {
          tautomerDirection = HistidineAtomNames.HD1.tautomerDirection;
        } else if (atomName.equals(HistidineAtomNames.HE2.name())) {
          tautomerDirection = HistidineAtomNames.HE2.tautomerDirection;
        }
      }
    }
    return tautomerDirection;
  }

  private void constructHISState(int biotypeCB, HisStates hisState) {
    int state = hisState.ordinal();
    for (HistidineAtomNames atomName : HistidineAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffsetHIS(hisState);
      if (offset < 0) {
        hisAtomTypes[index][state] = dummyHydrogenAtomType;
        // Zero out the MultipoleType and PolarizeType.
        if (hisState == HisStates.HID) {
          hisMultipoleTypes[index][state] = hidZeroMultipoleType;
        } else if (hisState == HisStates.HIE) {
          hisMultipoleTypes[index][state] = hieZeroMultipoleType;
        } else {
          logger.severe(" Error constructing HIS states.");
        }
        hisPolarizeTypes[index][state] = zeroPolarizeType;
        hisVDWTypes[index][state] = forceField.getVDWType(Integer.toString(0));
        hisSoluteTypes[index][state] = zeroSoluteType;
      } else {
        int biotype = biotypeCB + offset;
        hisAtomTypes[index][state] = findAtomType(biotype, forceField);
        String key = hisAtomTypes[index][state].getKey();
        hisMultipoleTypes[index][state] = forceField.getMultipoleTypeBeginsWith(key);
        hisPolarizeTypes[index][state] = forceField.getPolarizeType(key);
        int atomClass = hisAtomTypes[index][state].atomClass;
        hisVDWTypes[index][state] = forceField.getVDWType("" + atomClass);
        hisSoluteTypes[index][state] = getSoluteType(forceField, hisAtomTypes[index][state],
            hisVDWTypes[index][state]);
        if (hisMultipoleTypes[index][state] == null || hisPolarizeTypes[index][state] == null
            || hisSoluteTypes[index][state] == null) {
          logger.severe(format(" Titration parameters could not be assigned for His atom %s.\n %s\n",
              atomName, hisAtomTypes[index][state]));
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
        lysAtomTypes[index][state] = dummyHydrogenAtomType;
        // Zero out the MultipoleType and PolarizeType.
        lysMultipoleTypes[index][state] = lydZeroMultipoleType;
        lysPolarizeTypes[index][state] = zeroPolarizeType;
        lysVDWTypes[index][state] = forceField.getVDWType(Integer.toString(0));
        lysSoluteTypes[index][state] = zeroSoluteType;
      } else {
        int biotype = biotypeCB + offset;
        lysAtomTypes[index][state] = findAtomType(biotype, forceField);
        String key = lysAtomTypes[index][state].getKey();
        lysMultipoleTypes[index][state] = forceField.getMultipoleTypeBeginsWith(key);
        lysPolarizeTypes[index][state] = forceField.getPolarizeType(key);
        int atomClass = lysAtomTypes[index][state].atomClass;
        lysVDWTypes[index][state] = forceField.getVDWType("" + atomClass);
        lysSoluteTypes[index][state] = getSoluteType(forceField, lysAtomTypes[index][state],
            lysVDWTypes[index][state]);
        if (lysMultipoleTypes[index][state] == null || lysPolarizeTypes[index][state] == null
            || lysSoluteTypes[index][state] == null) {
          logger.severe(
              format(" Titration parameters could not be assigned for Lys atom %s.\n %s\n", atomName,
                  lysAtomTypes[index][state]));
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
        aspAtomTypes[index][state] = dummyHydrogenAtomType;
        // Zero out the MultipoleType and PolarizeType.
        if (aspState == AspStates.ASP) {
          aspMultipoleTypes[index][state] = aspZeroMultipoleType;
        } else {
          aspMultipoleTypes[index][state] = ashZeroMultipoleType;
        }
        aspPolarizeTypes[index][state] = zeroPolarizeType;
        aspVDWTypes[index][state] = forceField.getVDWType(Integer.toString(0));
        aspSoluteTypes[index][state] = zeroSoluteType;
      } else {
        int biotype = biotypeCB + offset;
        aspAtomTypes[index][state] = findAtomType(biotype, forceField);
        String key = aspAtomTypes[index][state].getKey();
        aspMultipoleTypes[index][state] = forceField.getMultipoleTypeBeginsWith(key);
        aspPolarizeTypes[index][state] = forceField.getPolarizeType(key);
        int atomClass = aspAtomTypes[index][state].atomClass;
        aspVDWTypes[index][state] = forceField.getVDWType("" + atomClass);
        aspSoluteTypes[index][state] = getSoluteType(forceField, aspAtomTypes[index][state],
            aspVDWTypes[index][state]);
        if (aspMultipoleTypes[index][state] == null || aspPolarizeTypes[index][state] == null
            || aspSoluteTypes[index][state] == null) {
          logger.severe(
              format(" Titration parameters could not be assigned for Asp atom %s.\n %s\n", atomName,
                  aspAtomTypes[index][state]));
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
        gluAtomTypes[index][state] = dummyHydrogenAtomType;
        // Zero out the MultipoleType and PolarizeType.
        if (gluState == GluStates.GLU) {
          gluMultipoleTypes[index][state] = gluZeroMultipoleType;
        } else {
          gluMultipoleTypes[index][state] = glhZeroMultipoleType;
        }
        gluPolarizeTypes[index][state] = zeroPolarizeType;
        gluVDWTypes[index][state] = forceField.getVDWType(Integer.toString(0));
        gluSoluteTypes[index][state] = zeroSoluteType;
      } else {
        int biotype = biotypeCB + offset;
        gluAtomTypes[index][state] = findAtomType(biotype, forceField);
        String key = gluAtomTypes[index][state].getKey();
        gluMultipoleTypes[index][state] = forceField.getMultipoleTypeBeginsWith(key);
        gluPolarizeTypes[index][state] = forceField.getPolarizeType(key);
        int atomClass = gluAtomTypes[index][state].atomClass;
        gluVDWTypes[index][state] = forceField.getVDWType("" + atomClass);
        gluSoluteTypes[index][state] = getSoluteType(forceField, gluAtomTypes[index][state],
            gluVDWTypes[index][state]);
        if (gluMultipoleTypes[index][state] == null || gluPolarizeTypes[index][state] == null
            || gluSoluteTypes[index][state] == null) {
          logger.severe(
              format(" Titration parameters could not be assigned for Glu atom %s.\n %s\n", atomName,
                  gluAtomTypes[index][state]));
        }
      }
    }
  }

  private void constructCYSState(int biotypeCB, CysStates cysState) {
    int state = cysState.ordinal();
    for (CysteineAtomNames atomName : CysteineAtomNames.values()) {
      int index = atomName.ordinal();
      int offset = atomName.getOffsetCYS(cysState);
      if (offset < 0) {
        // Set the AtomType to null.
        cysAtomTypes[index][state] = dummyHydrogenAtomType;
        // Zero out the MultipoleType and PolarizeType.
        cysMultipoleTypes[index][state] = cydZeroMultipoleType;
        cysPolarizeTypes[index][state] = zeroPolarizeType;
        cysVDWTypes[index][state] = forceField.getVDWType(Integer.toString(0));
        cysSoluteTypes[index][state] = zeroSoluteType;
      } else {
        int biotype = biotypeCB + offset;
        cysAtomTypes[index][state] = findAtomType(biotype, forceField);
        String key = cysAtomTypes[index][state].getKey();
        cysMultipoleTypes[index][state] = forceField.getMultipoleTypeBeginsWith(key);
        // This is an edge case since the CB/HB atom types have more than 1 matching multipole
        if (cysMultipoleTypes[index][state] == null) {
          if (cysState == CysStates.CYS) {
            if (atomName == CysteineAtomNames.CB) {
              cysMultipoleTypes[index][state] = forceField.getMultipoleType(key + " 8 45");
            } else {
              // HB2 & HB3
              cysMultipoleTypes[index][state] = forceField.getMultipoleType(key + " 43 8");
            }
          } else {
            if (atomName == CysteineAtomNames.CB) {
              cysMultipoleTypes[index][state] = forceField.getMultipoleType(key + " 8 49");
            } else {
              // HB2 & HB3
              cysMultipoleTypes[index][state] = forceField.getMultipoleType(key + " 43 8");
            }
          }
        }
        cysPolarizeTypes[index][state] = forceField.getPolarizeType(key);
        int atomClass = cysAtomTypes[index][state].atomClass;
        cysVDWTypes[index][state] = forceField.getVDWType("" + atomClass);
        cysSoluteTypes[index][state] = getSoluteType(forceField, cysAtomTypes[index][state],
            cysVDWTypes[index][state]);
        if (cysMultipoleTypes[index][state] == null || cysPolarizeTypes[index][state] == null
            || cysSoluteTypes[index][state] == null) {
          logger.severe(
              format(" Titration parameters could not be assigned for Cys atom %s.\n %s\n", atomName,
                  cysAtomTypes[index][state]));
        }
      }
    }
  }

  private void checkParameterTypes(String label, AtomType[][] atomTypes,
      PolarizeType[][] polarizeTypes, MultipoleType[][] multipoleTypes, VDWType[][] vdwTypes) {
    int states = multipoleTypes.length;
    int types = multipoleTypes[0].length;
    StringBuilder sb = new StringBuilder();
    for (int t = 0; t < types; t++) {
      MultipoleFrameDefinition frame0 = multipoleTypes[0][t].frameDefinition;
      double eps0 = vdwTypes[0][t].wellDepth;
      double rad0 = vdwTypes[0][t].radius;
      sb.append(format("\n %s Type %d\n", label, t));
      sb.append(format(" %s\n  %s\n  %s\n  %s\n", atomTypes[0][t], polarizeTypes[0][t],
          multipoleTypes[0][t], vdwTypes[0][t]));
      for (int s = 1; s < states; s++) {
        sb.append(format(" %s\n  %s\n  %s\n  %s\n", atomTypes[s][t], polarizeTypes[s][t],
            multipoleTypes[s][t], vdwTypes[s][t]));
        MultipoleFrameDefinition frame = multipoleTypes[s][t].frameDefinition;

        if (!frame0.equals(frame)) {
          sb.append("\n Incompatible multipole frames:\n");
          sb.append(format(" %s\n  %s\n  %s\n", atomTypes[0][t], polarizeTypes[0][t],
              multipoleTypes[0][t]));
          sb.append(format(" %s\n  %s\n  %s\n", atomTypes[s][t], polarizeTypes[s][t],
              multipoleTypes[s][t]));
        }

        if (atomTypes[0][t].atomicNumber != 1) {
          double epsS = vdwTypes[s][t].wellDepth;
          double radS = vdwTypes[s][t].radius;
          if (epsS != eps0 || radS != rad0) {
            sb.append("\n Incompatible vdW types:\n");
            sb.append(format(" %s\n  %s\n", atomTypes[0][t], vdwTypes[0][t]));
            sb.append(format(" %s\n  %s\n", atomTypes[s][t], vdwTypes[s][t]));
          }
        }
      }
    }

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(sb.toString());
    }
  }

  private SoluteType getSoluteType(ForceField forceField, AtomType atomType, VDWType vdwType) {
    SoluteType soluteType = SoluteType.getCensusSoluteType(atomType.atomicNumber);
    switch (soluteRadiiType) {
      case SOLUTE:
        SoluteType type = SoluteType.getFitSoluteType(forceField, atomType.type);
        if (type != null) {
          soluteType = type;
        }
        break;
      case VDW:
        type = SoluteType.getVDWSoluteType(vdwType);
        if (type != null) {
          soluteType = type;
        }
        break;
    }
    if (soluteType == null) {
      logger.severe(
          format(" No solute type (%s) for %d:\n  \"%s\"\n  %s", soluteRadiiType, atomType.type,
              atomType, vdwType));
    }
    return soluteType;
  }

  public void setRotamerPhBias(double temperature, double pH) {
    /*
     * Set ASH pH bias as sum of Fmod and acidostat energy
     */
    rotamerPhBiasMap.put(ASH, 0.0);

    /*
     * Set ASP pH bias as sum of Fmod and acidostat energy
     */
    double acidostat = LOG10 * Constants.R * temperature * (Titration.ASHtoASP.pKa - pH);
    double fMod = Titration.ASHtoASP.freeEnergyDiff;
    rotamerPhBiasMap.put(ASP, acidostat - fMod);

    /*
     * Set ASH pH bias as sum of Fmod and acidostat energy
     */
    rotamerPhBiasMap.put(GLH, 0.0);

    /*
     * Set GLU pH bias as sum of Fmod and acidostat energy
     */
    acidostat = LOG10 * Constants.R * temperature * (Titration.GLHtoGLU.pKa - pH);
    fMod = Titration.GLHtoGLU.freeEnergyDiff;
    rotamerPhBiasMap.put(GLU, acidostat - fMod);


    /*
     * Set LYS pH bias as sum of Fmod and acidostat energy
     */
    rotamerPhBiasMap.put(LYS, 0.0);

    /*
     * Set LYD pH bias as sum of Fmod and acidostat energy
     */
    acidostat = LOG10 * Constants.R * temperature * (Titration.LYStoLYD.pKa - pH);
    fMod = Titration.LYStoLYD.freeEnergyDiff;
    rotamerPhBiasMap.put(LYD, acidostat - fMod);

    /*
     * Set LYS pH bias as sum of Fmod and acidostat energy
     */
    rotamerPhBiasMap.put(CYS, 0.0);

    /*
     * Set LYD pH bias as sum of Fmod and acidostat energy
     */
    acidostat = LOG10 * Constants.R * temperature * (Titration.CYStoCYD.pKa - pH);
    fMod = Titration.CYStoCYD.freeEnergyDiff;
    rotamerPhBiasMap.put(CYD, acidostat - fMod);

    /*
     * Set HIS pH bias as sum of Fmod and acidostat energy
     */
    rotamerPhBiasMap.put(HIS, 0.0);

    /*
     * Set HID pH bias as sum of Fmod and acidostat energy
     */
    acidostat = LOG10 * Constants.R * temperature * (Titration.HIStoHID.pKa - pH);
    fMod = Titration.HIStoHID.freeEnergyDiff;
    rotamerPhBiasMap.put(HID, acidostat - fMod);

    /*
     * Set HIE pH bias as sum of Fmod and acidostat energy
     */
    acidostat = LOG10 * Constants.R * temperature * (Titration.HIStoHIE.pKa - pH);
    fMod = Titration.HIStoHIE.freeEnergyDiff;
    rotamerPhBiasMap.put(HIE, acidostat - fMod);
  }

  public double getRotamerPhBias(AminoAcid3 AA3) {
    return rotamerPhBiasMap.getOrDefault(AA3, 0.0);
  }

  public double getTotalRotamerPhBias(Rotamer[] rotamers) {
    double total = 0.0;
    for (Rotamer r : rotamers) {
      if (r.isTitrating) {
        total += getRotamerPhBias(r.aminoAcid3);
      }
    }
    return total;
  }

  /**
   * Amino acid protonation reactions. Constructors below specify intrinsic pKa and reference free
   * energy of protonation, obtained via BAR on capped monomers. pKa values from Thurlkill, Richard
   * L., et al. "pK values of the ionizable groups of proteins." Protein science 15.5 (2006):
   * 1214-1218.
   * <p>
   * HIS to HID/HIE pKa values from Bashford, Donald, et al. "Electrostatic calculations of
   * side-chain pKa values in myoglobin and comparison with NMR data for histidines." Biochemistry
   * 32.31 (1993): 8045-8056.
   * <p>
   * -(quadratic * lambda^2 + linear * lambda)
   */
  public enum Titration {


    ASHtoASP(3.67, -71.10, 13.817, -105.690, 162.780, AminoAcid3.ASH, AminoAcid3.ASP),
    GLHtoGLU(4.25, -83.40, 26.619, -128.530, 187.210, AminoAcid3.GLH, AminoAcid3.GLU),
    LYStoLYD(10.40, 41.77, 0.0, -69.29, 24.17778, AminoAcid3.LYS, AminoAcid3.LYD),
    //TYRtoTYD(10.07, 34.961, 0.0, AminoAcidUtils.AminoAcid3.TYR, AminoAcidUtils.AminoAcid3.TYD),
    CYStoCYD(8.55, -66.2, 39.039, -170.920, 216.620, AminoAcid3.CYS, AminoAcid3.CYD),
    //HE2 is the proton that is lost
    HIStoHID(7.00, 40.20, 0.0, -64.317, 30.35, AminoAcid3.HIS, AminoAcid3.HID),
    //HD1 is the proton that is lost
    HIStoHIE(6.60, 37.44, 0.0, -62.931, 32.00, AminoAcid3.HIS, AminoAcid3.HIE),
    HIDtoHIE(Double.NaN, 0.00, 0.0, -36.83, 34.325, AminoAcid3.HID, AminoAcid3.HIE);

    //TerminalNH3toNH2(8.23, 0.0, 00.00, AminoAcidUtils.AminoAcid3.UNK, AminoAcidUtils.AminoAcid3.UNK),
    //TerminalCOOHtoCOO(3.55, 0.0, 00.00, AminoAcidUtils.AminoAcid3.UNK, AminoAcidUtils.AminoAcid3.UNK);


    public final double pKa;
    // Free energy differences used in rotamer optimization
    public final double freeEnergyDiff;
    public final double cubic;
    public final double quadratic;
    public final double linear;
    public final AminoAcid3 protForm;
    public final AminoAcid3 deprotForm;

    /** Invoked by Enum; use the factory method to obtain instances. */

    Titration(double pKa, double freeEnergyDiff, double cubic, double quadratic, double linear,
        AminoAcid3 protForm, AminoAcid3 deprotForm) {
      this.pKa = pKa;
      this.freeEnergyDiff = freeEnergyDiff;
      this.cubic = cubic;
      this.quadratic = quadratic;
      this.linear = linear;
      this.protForm = protForm;
      this.deprotForm = deprotForm;
    }
  }
}
