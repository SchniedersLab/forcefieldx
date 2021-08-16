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

import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static ffx.potential.bonded.AminoAcidUtils.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.AminoAcidUtils.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.AminoAcidUtils.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.bonded.BondedUtils.buildBond;
import static ffx.potential.bonded.BondedUtils.buildH;
import static ffx.potential.bonded.BondedUtils.buildHeavy;
import static ffx.potential.bonded.BondedUtils.buildHydrogenAtom;
import static ffx.potential.bonded.BondedUtils.findAtomType;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.math.Double3;
import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utilities for creating Amino Acid residues.
 *
 * @author Michael Schnieders
 * @since 1.0
 */
public class AminoAcidUtils {

  private static final Logger logger = Logger.getLogger(AminoAcidUtils.class.getName());

  /** Repeating atomic numbers of an amino acid chain. */
  public static final int[] AAPATTERN = {7, 6, 6};
  /** Constant <code>AA3toAA1</code> */
  static final HashMap<AminoAcid3, AminoAcid1> AA3toAA1 = new HashMap<>();
  /** Constant <code>AA1toAA3</code> */
  public static final HashMap<AminoAcid1, AminoAcid3> AA1toAA3 = new HashMap<>();
  /** List of values from the AminoAcid1 enum. */
  public static final List<AminoAcid1> aminoAcid1List = Arrays.asList(AminoAcid1.values());
  /** Constant <code>aminoAcidList</code> */
  public static final List<AminoAcid3> aminoAcidList = Arrays.asList(AminoAcid3.values());

  /** Constant <code>AminoAcidBackboneAtoms</code> */
  public enum AminoAcidBackboneAtoms {N, H, CA, HA, C, O}

  /** Constant <code>GlycineBackboneAtoms</code> */
  public enum GlycineBackboneAtoms {N, H, CA, HA2, HA3, C, O}

  /** Constant <code>ProlineBackboneAtoms</code> */
  public enum ProlineBackboneAtoms {N, CA, HA, C, O}

  /**
   * Turn an Enum into String array.
   *
   * @param e The enum class.
   * @return A String array with all the enum values.
   */
  public static String[] getNames(Class<? extends Enum<?>> e) {
    return Arrays.stream(e.getEnumConstants()).map(Enum::name).toArray(String[]::new);
  }

  /**
   * Stoichiometry of side chains can be used for identification, accept for a couple cases: 1.)
   * Proline & Valine 2.) (Iso)Leucine 3.) DNA Gaunine/RNA Adenine. This Hashtable returns the
   * 3-letter name for amino acids, a single character for nucleic acids, or an integer indicating a
   * special case.
   */
  public static final HashMap<String, String> sidechainStoichiometry = new HashMap<>();

  static {
    // Amino Acid Side Chains
    sidechainStoichiometry.put("S1C3", "MET");
    sidechainStoichiometry.put("S1C1", "CYS");
    sidechainStoichiometry.put("O1C1", "SER");
    sidechainStoichiometry.put("O1C2", "THR");
    sidechainStoichiometry.put("O1C7", "TYR");
    sidechainStoichiometry.put("O2C2", "ASP");
    sidechainStoichiometry.put("O2C3", "GLU");
    sidechainStoichiometry.put("O1N1C2", "ASN");
    sidechainStoichiometry.put("O1N1C3", "GLN");
    sidechainStoichiometry.put("N3C4", "ARG");
    sidechainStoichiometry.put("N2C4", "HIS");
    sidechainStoichiometry.put("N1C9", "TRP");
    sidechainStoichiometry.put("N1C4", "LYS");
    sidechainStoichiometry.put("C7", "PHE");
    sidechainStoichiometry.put("H", "GLY");
    sidechainStoichiometry.put("C1", "ALA");
    // DNA
    sidechainStoichiometry.put("O2N3C6", "DC");
    sidechainStoichiometry.put("O1N5C7", "DA");
    sidechainStoichiometry.put("O3N2C7", "DT");
    // RNA
    sidechainStoichiometry.put("O3N5C7", "G");
    sidechainStoichiometry.put("O3N3C6", "C");
    sidechainStoichiometry.put("O4N2C6", "U");
    // SPECIAL CASES
    sidechainStoichiometry.put("C3", "1"); // Proline / Valine
    sidechainStoichiometry.put("C4", "2"); // (Iso)Leucine
    sidechainStoichiometry.put("O2N5C7", "3"); // DNA Gaunine / RNA Adenine

    AminoAcid1[] aa1 = AminoAcid1.values();
    AminoAcid3[] aa3 = AminoAcid3.values();
    for (int i = 0; i < AminoAcid1.values().length; i++) {
      AA1toAA3.put(aa1[i], aa3[i]);
      AA3toAA1.put(aa3[i], aa1[i]);
    }
  }

  /** Private constructor. */
  private AminoAcidUtils() {
  }

  /**
   * Assign atom types to an amino acid polymer.
   *
   * @param residues The residues to assign atom types to.
   * @param forceField The ForceField to apply.
   * @param bondList Created bonds are added to this List.
   * @throws MissingHeavyAtomException A needed heavy atom was not found.
   * @throws MissingAtomTypeException An atom type could not be found.
   * @since 1.0
   */
  public static void assignAminoAcidAtomTypes(
      List<Residue> residues, ForceField forceField, List<Bond> bondList)
      throws MissingHeavyAtomException, MissingAtomTypeException {
    // Loop over amino acid residues.
    int numberOfResidues = residues.size();
    for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
      Residue residue = residues.get(residueNumber);
      Residue previousResidue = null;
      Residue nextResidue = null;
      if (residueNumber > 0) {
        previousResidue = residues.get(residueNumber - 1);
      }
      if (residueNumber < numberOfResidues - 1) {
        nextResidue = residues.get(residueNumber + 1);
      }
      AminoAcidUtils.assignAminoAcidAtomTypes(
          residue, previousResidue, nextResidue, forceField, bondList);
    }
  }

  /**
   * assignAminoAcidAtomTypes.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param previousResidue a {@link ffx.potential.bonded.Residue} object.
   * @param nextResidue a {@link ffx.potential.bonded.Residue} object.
   * @param forceField a {@link ForceField} object.
   * @param bondList a {@link java.util.List} object.
   * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
   * @throws ffx.potential.bonded.BondedUtils.MissingAtomTypeException if any.
   */
  public static void assignAminoAcidAtomTypes(
      Residue residue,
      Residue previousResidue,
      Residue nextResidue,
      ForceField forceField,
      List<Bond> bondList)
      throws MissingHeavyAtomException, MissingAtomTypeException {

    String residueName = residue.getName().toUpperCase();

    int j = 1;
    ResiduePosition position = MIDDLE_RESIDUE;
    if (previousResidue == null) {
      j = 0;
      position = FIRST_RESIDUE;
    } else if (nextResidue == null) {
      j = 2;
      position = LAST_RESIDUE;
      // If the last residue only contains a nitrogen turn it into an NH2 group.
      Atom N = (Atom) residue.getAtomNode("N");
      if (residue.getAtomNodeList().size() == 1 && N != null) {
        residueName = "NH2";
        residue.setName(residueName);
      }
    }

    AminoAcid3 aminoAcid = getAminoAcid(residueName);
    int aminoAcidNumber = getAminoAcidNumber(residueName);
    // Non-standard Amino Acid; use ALA backbone types.
    boolean nonStandard = false;
    if (aminoAcid == AminoAcid3.UNK) {
      aminoAcidNumber = getAminoAcidNumber("ALA");
      nonStandard = true;
    }

    // Only the last residue in a chain should have an OXT/OT2 atom.
    if (nextResidue != null) {
      removeOXT_OT2(residue);
    }

    // Only the first nitrogen should have H1, H2 and H3 atoms, unless it's an NME cap.
    if (previousResidue != null) {
      removeH1_H2_H3(aminoAcid, residue);
    }

    // Check for missing heavy atoms. This check ignores special terminating groups like FOR, NH2,
    // etc.
    if (!nonStandard) {
      try {
        checkForMissingHeavyAtoms(aminoAcid, residue);
      } catch (BondedUtils.MissingHeavyAtomException e) {
        logger.log(Level.INFO, " {0} could not be parsed.", residue.toString());
        throw e;
      }
    }

    Atom pC = null;
    Atom pCA = null;
    if (previousResidue != null) {
      pC = (Atom) previousResidue.getAtomNode("C");
      pCA = (Atom) previousResidue.getAtomNode("CA");
    }

    // Backbone heavy atoms.
    Atom N = (Atom) residue.getAtomNode("N");
    if (N != null) {
      N.setAtomType(findAtomType(AA_N[j][aminoAcidNumber], forceField));
      if (position != FIRST_RESIDUE) {
        buildBond(pC, N, forceField, bondList);
      }
    }

    Atom CA = null;
    Atom C = null;
    Atom O = null;
    if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NH2)) {
      if (aminoAcid == AminoAcid3.ACE || aminoAcid == AminoAcid3.NME) {
        CA = buildHeavy(residue, "CH3", N, AA_CA[j][aminoAcidNumber], forceField, bondList);
      } else {
        CA = buildHeavy(residue, "CA", N, AA_CA[j][aminoAcidNumber], forceField, bondList);
      }
      if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NME)) {
        C = buildHeavy(residue, "C", CA, AA_C[j][aminoAcidNumber], forceField, bondList);
        O = (Atom) residue.getAtomNode("O");
        if (O == null) {
          O = (Atom) residue.getAtomNode("OT1");
        }
        if (O == null) {
          // Build the carbonyl oxygen; assume the final residue has beta-sheet psi
          double psi = 135.0;
          if (nextResidue != null && N != null) {
            Atom nextN = (Atom) nextResidue.getAtomNode("N");
            psi =
                toDegrees(
                    dihedralAngle(
                        N.getXYZ(null), CA.getXYZ(null), C.getXYZ(null), nextN.getXYZ(null)));
          }
          O =
              buildHeavy(
                  residue,
                  "O",
                  C,
                  1.25,
                  CA,
                  117.0,
                  N,
                  psi - 180.0,
                  0,
                  AA_O[j][aminoAcidNumber],
                  forceField,
                  bondList);
        } else {
          AtomType atomType = findAtomType(AA_O[j][aminoAcidNumber], forceField);
          O.setAtomType(atomType);
          buildBond(C, O, forceField, bondList);
        }
      }
    }

    // Nitrogen hydrogen atoms.
    AtomType atomType = findAtomType(AA_HN[j][aminoAcidNumber], forceField);
    switch (position) {
      case FIRST_RESIDUE:
        switch (aminoAcid) {
          case PRO:
            // Check if, for some reason, the structure has an H1 but not an H2 or H3.
            List<Atom> resAtoms = residue.getAtomList();
            Optional<Atom> H1 =
                resAtoms.stream().filter((Atom a) -> a.getName().equals("H1")).findAny();
            Optional<Atom> H2 =
                resAtoms.stream().filter((Atom a) -> a.getName().equals("H2")).findAny();
            Optional<Atom> H3 =
                resAtoms.stream().filter((Atom a) -> a.getName().equals("H3")).findAny();
            if (H1.isPresent() && H2.isPresent() && H3.isPresent()) {
              logger.severe(
                  String.format(
                      " Proline residue %s somehow has three N-terminal hydrogens!", residue));
            }
            if (H1.isPresent() && H2.isPresent()) {
              H1.get().setName("H3");
            } else if (H1.isPresent() && H3.isPresent()) {
              H1.get().setName("H2");
            } // Else, all is hunky-dory!
            buildHydrogenAtom(
                residue, "H2", N, 1.02, CA, 109.5, C, 0.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H3", N, 1.02, CA, 109.5, C, -120.0, 0, atomType, forceField, bondList);
            break;
          case PCA:
            buildHydrogenAtom(
                residue, "H", N, 1.02, CA, 109.5, C, -60.0, 0, atomType, forceField, bondList);
            break;
          case ACE:
            break;
          default:
            buildHydrogenAtom(
                residue, "H1", N, 1.02, CA, 109.5, C, 180.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H2", N, 1.02, CA, 109.5, C, 60.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H3", N, 1.02, CA, 109.5, C, -60.0, 0, atomType, forceField, bondList);
        }
        break;
      case LAST_RESIDUE:
        switch (aminoAcid) {
          case NH2:
            buildHydrogenAtom(
                residue, "H1", N, 1.02, pC, 119.0, pCA, 0.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H2", N, 1.02, pC, 119.0, pCA, 180.0, 0, atomType, forceField, bondList);
            break;
          case NME:
            buildHydrogenAtom(
                residue, "H", N, 1.02, pC, 118.0, CA, 121.0, 1, atomType, forceField, bondList);
            break;
          default:
            buildHydrogenAtom(
                residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1, atomType, forceField, bondList);
        }
        break;
      default:
        // Mid-chain nitrogen hydrogen.
        buildHydrogenAtom(
            residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1, atomType, forceField, bondList);
    }

    // C-alpha hydrogen atoms.
    String haName = "HA";
    if (aminoAcid == AminoAcid3.GLY) {
      haName = "HA2";
    }
    atomType = findAtomType(AA_HA[j][aminoAcidNumber], forceField);
    switch (position) {
      case FIRST_RESIDUE:
        switch (aminoAcid) {
          case FOR:
            buildHydrogenAtom(
                residue, "H", C, 1.12, O, 0.0, null, 0.0, 0, atomType, forceField, bondList);
            break;
          case ACE:
            buildHydrogenAtom(
                residue, "H1", CA, 1.10, C, 109.5, O, 180.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H2", CA, 1.10, C, 109.5, O, 60.0, 0, atomType, forceField, bondList);
            buildHydrogenAtom(
                residue, "H3", CA, 1.10, C, 109.5, O, -60.0, 0, atomType, forceField, bondList);
            break;
          default:
            buildHydrogenAtom(
                residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1, atomType, forceField, bondList);
            break;
        }
        break;
      case LAST_RESIDUE:
        if (aminoAcid == AminoAcid3.NME) {
          buildHydrogenAtom(
              residue, "H1", CA, 1.10, N, 109.5, pC, 180.0, 0, atomType, forceField, bondList);
          buildHydrogenAtom(
              residue, "H2", CA, 1.10, N, 109.5, pC, 60.0, 0, atomType, forceField, bondList);
          buildHydrogenAtom(
              residue, "H3", CA, 1.10, N, 109.5, pC, -60.0, 0, atomType, forceField, bondList);
        } else {
          buildHydrogenAtom(
              residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1, atomType, forceField, bondList);
        }
        break;
      default:
        buildHydrogenAtom(
            residue, haName, CA, 1.10, N, 109.5, C, 109.0, -1, atomType, forceField, bondList);
    }

    // Build the amino acid side chain.
    assignAminoAcidSideChain(position, aminoAcid, residue, CA, N, C, forceField, bondList);

    // Build the terminal oxygen if the residue is not NH2 or NME.
    if (position == LAST_RESIDUE && !(aminoAcid == AminoAcid3.NH2 || aminoAcid == AminoAcid3.NME)) {
      atomType = findAtomType(AA_O[2][aminoAcidNumber], forceField);
      Atom OXT = (Atom) residue.getAtomNode("OXT");
      if (OXT == null) {
        OXT = (Atom) residue.getAtomNode("OT2");
        if (OXT != null) {
          OXT.setName("OXT");
        }
      }
      if (OXT == null) {
        String resName = C.getResidueName();
        int resSeq = C.getResidueNumber();
        Character chainID = C.getChainID();
        Character altLoc = C.getAltLoc();
        String segID = C.getSegID();
        double occupancy = C.getOccupancy();
        double tempFactor = C.getTempFactor();
        OXT =
            new Atom(
                0,
                "OXT",
                altLoc,
                new double[3],
                resName,
                resSeq,
                chainID,
                occupancy,
                tempFactor,
                segID);
        OXT.setAtomType(atomType);
        residue.addMSNode(OXT);
        intxyz(OXT, C, 1.25, CA, 117.0, O, 126.0, 1);
      } else {
        OXT.setAtomType(atomType);
      }
      buildBond(C, OXT, forceField, bondList);
    }

    // Do some checks on the current residue to make sure all atoms have been assigned an atom type.
    List<Atom> resAtoms = residue.getAtomList();
    for (Atom atom : resAtoms) {
      atomType = atom.getAtomType();
      if (atomType == null) {
        /*
         * Sometimes, with deuterons, a proton has been constructed in
         * its place, so we have a "dummy" deuteron still hanging
         * around.
         */
        String protonEq = atom.getName().replaceFirst("D", "H");
        Atom correspH = (Atom) residue.getAtomNode(protonEq);
        if (correspH == null || correspH.getAtomType() == null) {
          MissingAtomTypeException missingAtomTypeException =
              new MissingAtomTypeException(residue, atom);
          logger.warning("MissingAtomTypeException incoming from 393.");
          throw missingAtomTypeException;
        } else {
          correspH.setName(atom.getName());
          atom.removeFromParent();
          atom = correspH;
          atomType = atom.getAtomType();
        }
      }
      int numberOfBonds = atom.getNumBonds();
      if (numberOfBonds != atomType.valence) {
        if (atom == C && numberOfBonds == atomType.valence - 1 && position != LAST_RESIDUE) {
          continue;
        }
        logger.warning(
            format(
                " An atom for residue %s has the wrong number of bonds:\n %s",
                residueName, atom));
        logger.info(format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
        for (Bond bond : atom.getBonds()) {
          logger.info(" " + bond.toString());
        }
      }
    }
  }

  /**
   * This interface is used by the "Build" routines.
   */
  public interface SideChainType {

    /**
     * This is already implemented by all enum instances.
     */
    String name();

    /**
     * Returns Biotype of this atom.
     *
     * @return The integer biotype.
     */
    int getType();

  }

  /** Constant <code>AIB</code> */
  public enum AIB implements SideChainType {
    CB1(0), CB2(0),
    HB11(1), HB12(1), HB13(1),
    HB21(1), HB22(1), HB23(1);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    AIB(int offset) {
      biotype = offset + AA_CB[AminoAcid3.AIB.ordinal()];
    }
  }

  /** Constant <code>ALA</code> */
  public enum ALA implements SideChainType {
    CB(0),
    HB1(1),
    HB2(1),
    HB3(1);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ALA(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ALA.ordinal()];
    }
  }

  /** Constant <code>ARG</code> */
  public enum ARG implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4), HD2(5), HD3(5),
    NE(6), HE(7),
    CZ(8),
    NH1(9), HH11(10), HH12(10),
    NH2(9), HH21(10), HH22(10);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ARG(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ARG.ordinal()];
    }
  }

  /** Constant <code>ASN</code> */
  public enum ASN implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    OD1(3), ND2(4), HD21(5), HD22(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ASN(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ASN.ordinal()];
    }
  }

  /** Constant <code>ASP</code> */
  public enum ASP implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    OD1(3), OD2(3);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ASP(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ASP.ordinal()];
    }
  }

  /** Constant <code>ASH</code> */
  public enum ASH implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    OD1(3), OD2(4), HD2(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ASH(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ASH.ordinal()];
    }
  }

  /** Constant <code>ASD</code> */
  public enum ASD implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    OD1(3), OD2(3), HD1(4), HD2(4);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ASD(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ASD.ordinal()];
    }
  }

  /** Constant <code>CYS</code> */
  public enum CYS implements SideChainType {
    CB(0), HB2(1), HB3(1),
    SG(2), HG(3);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    CYS(int offset) {
      biotype = offset + AA_CB[AminoAcid3.CYS.ordinal()];
    }
  }

  /** Constant <code>CYS</code> */
  public enum CYX implements SideChainType {
    CB(0), HB2(1), HB3(1),
    SG(2);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    CYX(int offset) {
      biotype = offset + AA_CB[AminoAcid3.CYX.ordinal()];
    }
  }

  /** Constant <code>CYD</code> */
  public enum CYD implements SideChainType {
    CB(0), HB2(1), HB3(1),
    SG(2);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    CYD(int offset) {
      biotype = offset + AA_CB[AminoAcid3.CYD.ordinal()];
    }
  }

  /** Constant <code>GLU</code> */
  public enum GLU implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4),
    OE1(5), OE2(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    GLU(int offset) {
      biotype = offset + AA_CB[AminoAcid3.GLU.ordinal()];
    }
  }

  /** Constant <code>GLH</code> */
  public enum GLH implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4),
    OE1(5), OE2(6), HE2(7);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    GLH(int offset) {
      biotype = offset + AA_CB[AminoAcid3.GLH.ordinal()];
    }
  }

  /** Constant <code>GLD</code> */
  public enum GLD implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4),
    OE1(5), OE2(5), HE1(6), HE2(6);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    GLD(int offset) {
      biotype = offset + AA_CB[AminoAcid3.GLD.ordinal()];
    }
  }

  /** Constant <code>GlutamineAtomNames</code> */
  public enum GLN implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4),
    OE1(5), NE2(6), HE21(7), HE22(7);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    GLN(int offset) {
      biotype = offset + AA_CB[AminoAcid3.GLN.ordinal()];
    }
  }

  /** Constant <code>GLY</code> */
  public enum GLY {HA2}

  /** Constant <code>HistidineAtoms</code> */
  public enum HIS implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    ND1(3), HD1(4),
    CD2(5), HD2(6),
    CE1(7), HE1(8),
    NE2(9), HE2(10);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    HIS(int offset) {
      biotype = offset + AA_CB[AminoAcid3.HIS.ordinal()];
    }
  }

  /** Constant <code>HID</code> */
  public enum HID implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    ND1(3), HD1(4),
    CD2(5), HD2(6),
    CE1(7), HE1(8),
    NE2(9);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    HID(int offset) {
      biotype = offset + AA_CB[AminoAcid3.HID.ordinal()];
    }
  }

  /** Constant <code>HIE</code> */
  public enum HIE implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    ND1(3),
    CD2(4), HD2(5),
    CE1(6), HE1(7),
    NE2(8), HE2(9);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    HIE(int offset) {
      biotype = offset + AA_CB[AminoAcid3.HIE.ordinal()];
    }
  }

  /** Constant <code>ILE</code> */
  public enum ILE implements SideChainType {
    CB(0), HB(1),
    CG1(2), HG12(3), HG13(3),
    CG2(4), HG21(5), HG22(5), HG23(5),
    CD1(6), HD11(7), HD12(7), HD13(7);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ILE(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ILE.ordinal()];
    }
  }

  /** Constant <code>LEU</code> */
  public enum LEU implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG(3),
    CD1(4), HD11(5), HD12(5), HD13(5),
    CD2(6), HD21(7), HD22(7), HD23(7);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    LEU(int offset) {
      biotype = offset + AA_CB[AminoAcid3.LEU.ordinal()];
    }
  }

  /** Constant <code>LYS</code> */
  public enum LYS implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4), HD2(5), HD3(5),
    CE(6), HE2(7), HE3(7),
    NZ(8), HZ1(9), HZ2(9), HZ3(9);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    LYS(int offset) {
      biotype = offset + AA_CB[AminoAcid3.LYS.ordinal()];
    }
  }

  /** Constant <code>LYD</code> */
  public enum LYD implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4), HD2(5), HD3(5),
    CE(6), HE2(7), HE3(7),
    NZ(8), HZ1(9), HZ2(9);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    LYD(int offset) {
      biotype = offset + AA_CB[AminoAcid3.LYD.ordinal()];
    }
  }

  /** Constant <code>MethionineAtomNames</code> */
  public enum MET implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    SD(4),
    CE(5), HE1(6), HE2(6), HE3(6);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    MET(int offset) {
      biotype = offset + AA_CB[AminoAcid3.MET.ordinal()];
    }
  }

  /** Constant <code>ORN</code> */
  public enum ORN implements SideChainType {
    CB(0), HB2(1), HB3(2),
    CG(2), HG2(3), HG3(3),
    CD(4), HD2(5), HD3(5),
    NE(6), HE1(7), HE2(7), HE3(7);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    ORN(int offset) {
      biotype = offset + AA_CB[AminoAcid3.ORN.ordinal()];
    }
  }

  /** Constant <code>PCA</code> */
  public enum PCA implements SideChainType {
    CB(0), HB2(1), HB3(2),
    CG(2), HG2(3), HG3(3),
    CD(4),
    OE(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    PCA(int offset) {
      biotype = offset + AA_CB[AminoAcid3.PCA.ordinal()];
    }
  }

  /** Constant <code>PHE</code> */
  public enum PHE implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    CD1(3), HD1(4),
    CD2(3), HD2(4),
    CE1(5), HE1(6),
    CE2(5), HE2(6),
    CZ(7), HZ(8);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    PHE(int offset) {
      biotype = offset + AA_CB[AminoAcid3.PHE.ordinal()];
    }
  }

  /** Constant <code>PRO</code> */
  public enum PRO implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2), HG2(3), HG3(3),
    CD(4), HD2(5), HD3(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    PRO(int offset) {
      biotype = offset + AA_CB[AminoAcid3.PRO.ordinal()];
    }
  }

  /** Constant <code>SER</code> */
  public enum SER implements SideChainType {
    CB(0), HB2(1), HB3(1),
    OG(2), HG(3);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    SER(int offset) {
      biotype = offset + AA_CB[AminoAcid3.SER.ordinal()];
    }
  }

  /** Constant <code>THR</code> */
  public enum THR implements SideChainType {
    CB(0), HB(1),
    OG1(2), HG1(3),
    CG2(4), HG21(5), HG22(5), HG23(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    THR(int offset) {
      biotype = offset + AA_CB[AminoAcid3.THR.ordinal()];
    }
  }

  /** Constant <code>TRP</code> */
  public enum TRP implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    CD1(3), HD1(4),
    CD2(5),
    NE1(6), HE1(7),
    CE2(8),
    CE3(9), HE3(10),
    CZ2(11), HZ2(12),
    CZ3(13), HZ3(14),
    CH2(15), HH2(16);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    TRP(int offset) {
      biotype = offset + AA_CB[AminoAcid3.TRP.ordinal()];
    }
  }

  /** Constant <code>TYR</code> */
  public enum TYR implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    CD1(3), HD1(4),
    CD2(3), HD2(4),
    CE1(5), HE1(6),
    CE2(5), HE2(6),
    CZ(7),
    OH(8), HH(9);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    TYR(int offset) {
      biotype = offset + AA_CB[AminoAcid3.TYR.ordinal()];
    }
  }

  /** Constant <code>TYD</code> */
  public enum TYD implements SideChainType {
    CB(0), HB2(1), HB3(1),
    CG(2),
    CD1(3), HD1(4),
    CD2(3), HD2(4),
    CE1(5), HE1(6),
    CE2(5), HE2(6),
    CZ(7),
    OH(8);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    TYD(int offset) {
      biotype = offset + AA_CB[AminoAcid3.TYD.ordinal()];
    }
  }

  /** Constant <code>VAL</code> */
  public enum VAL implements SideChainType {
    CB(0), HB(1),
    CG1(2), HG11(3), HG12(3), HG13(3),
    CG2(4), HG21(5), HG22(5), HG23(5);

    /**
     * Biotype for this atom.
     */
    private final int biotype;

    /**
     * {@inheritDoc}
     */
    public int getType() {
      return biotype;
    }

    /**
     * Init the Atom.
     *
     * @param offset Biotype offset relative to the CB biotype.
     */
    VAL(int offset) {
      biotype = offset + AA_CB[AminoAcid3.VAL.ordinal()];
    }
  }

  /**
   * buildAIB.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildAIB(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB1 = buildHeavy(res, AIB.CB1, CA, 1.54, N, 109.5, C, 107.8, -1, ff, bonds);
    Atom CB2 = buildHeavy(res, AIB.CB2, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom HB11 = buildH(res, AIB.HB11, CB1, 1.11, CA, 109.4, N, 180.0, 0, ff, bonds);
    buildH(res, AIB.HB12, CB1, 1.11, CA, 109.4, HB11, 109.4, 1, ff, bonds);
    buildH(res, AIB.HB13, CB1, 1.11, CA, 109.4, HB11, 109.4, -1, ff, bonds);
    Atom HB21 = buildH(res, AIB.HB21, CB2, 1.11, CA, 109.4, N, 180.0, 0, ff, bonds);
    buildH(res, AIB.HB22, CB2, 1.11, CA, 109.4, HB21, 109.4, 1, ff, bonds);
    buildH(res, AIB.HB23, CB2, 1.11, CA, 109.4, HB21, 109.4, -1, ff, bonds);
    return res;
  }

  /**
   * buildGlycine.
   *
   * @param res a {@link ffx.potential.bonded.Residue} object.
   * @param CA a {@link ffx.potential.bonded.Atom} object.
   * @param N a {@link ffx.potential.bonded.Atom} object.
   * @param C a {@link ffx.potential.bonded.Atom} object.
   * @param position {@link ffx.potential.bonded.AminoAcidUtils.ResiduePosition} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link java.util.List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildGlycine(Residue res, Atom CA, Atom N, Atom C, ResiduePosition position,
      ForceField ff,
      List<Bond> bonds) {
    int k = AA_CB[AminoAcid3.GLY.ordinal()];
    switch (position) {
      case FIRST_RESIDUE:
        k = AA_HA[0][k];
        break;
      case LAST_RESIDUE:
        k = AA_HA[2][k];
        break;
      default:
        k = AA_HA[1][k];
    }
    buildH(res, "HA3", CA, 1.10, N, 109.5, C, 109.5, 1, k, ff, bonds);
    return res;
  }

  /**
   * buildAlanine.
   *
   * @param res a {@link ffx.potential.bonded.Residue} object.
   * @param CA a {@link ffx.potential.bonded.Atom} object.
   * @param N a {@link ffx.potential.bonded.Atom} object.
   * @param C a {@link ffx.potential.bonded.Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link java.util.List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildAlanine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ALA.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom HB1 = buildH(res, ALA.HB1, CB, 1.11, CA, 109.4, N, 180.0, 0, ff, bonds);
    buildH(res, ALA.HB2, CB, 1.11, CA, 109.4, HB1, 109.4, 1, ff, bonds);
    buildH(res, ALA.HB3, CB, 1.11, CA, 109.4, HB1, 109.4, -1, ff, bonds);
    return res;
  }

  /**
   * buildArginine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildArginine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ARG.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ARG.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, ARG.CD, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    Atom NE = buildHeavy(res, ARG.NE, CD, 1.45, CG, 109.5, CB, 180, 0, ff, bonds);
    Atom CZ = buildHeavy(res, ARG.CZ, NE, 1.35, CD, 120.0, CG, 180, 0, ff, bonds);
    Atom NH1 = buildHeavy(res, ARG.NH1, CZ, 1.35, NE, 120.0, CD, 180, 0, ff, bonds);
    Atom NH2 = buildHeavy(res, ARG.NH2, CZ, 1.35, NE, 120.0, NH1, 120.0, 1, ff, bonds);
    buildH(res, ARG.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, ARG.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, ARG.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, ARG.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    buildH(res, ARG.HD2, CD, 1.11, CG, 109.4, NE, 109.4, 1, ff, bonds);
    buildH(res, ARG.HD3, CD, 1.11, CG, 109.4, NE, 109.4, -1, ff, bonds);
    buildH(res, ARG.HE, NE, 1.02, CD, 120.0, CZ, 120.0, 1, ff, bonds);
    Atom HH11 = buildH(res, ARG.HH11, NH1, 1.02, CZ, 120.0, NE, 180.0, 0, ff, bonds);
    buildH(res, ARG.HH12, NH1, 1.02, CZ, 120.0, HH11, 120.0, 1, ff, bonds);
    Atom HH21 = buildH(res, ARG.HH21, NH2, 1.02, CZ, 120.0, NE, 180.0, 0, ff, bonds);
    buildH(res, ARG.HH22, NH2, 1.02, CZ, 120.0, HH21, 120.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildAsparagine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildAsparagine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ASN.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ASN.CG, CB, 1.51, CA, 107.8, N, 180, 0, ff, bonds);
    Atom OD1 = buildHeavy(res, ASN.OD1, CG, 1.22, CB, 122.5, CA, 180, 0, ff, bonds);
    Atom ND2 = buildHeavy(res, ASN.ND2, CG, 1.34, CB, 112.7, OD1, 124.0, 1, ff, bonds);
    buildH(res, ASN.HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1, ff, bonds);
    buildH(res, ASN.HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1, ff, bonds);
    Atom HD21 = buildH(res, ASN.HD21, ND2, 1.02, CG, 119.0, CB, 0.0, 0, ff, bonds);
    buildH(res, ASN.HD22, ND2, 1.02, CG, 119.0, HD21, 120.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildAspartate.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildAspartate(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ASP.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ASP.CG, CB, 1.51, CA, 107.8, N, 180, 0, ff, bonds);
    Atom OD1 = buildHeavy(res, ASP.OD1, CG, 1.25, CB, 117.0, CA, 0.0, 0, ff, bonds);
    buildHeavy(res, ASP.OD2, CG, 1.25, CB, 117.0, OD1, 126.0, 1, ff, bonds);
    buildH(res, ASP.HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1, ff, bonds);
    buildH(res, ASP.HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1, ff, bonds);
    return res;
  }

  /**
   * buildCysteine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildCysteine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, CYS.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom SG = buildHeavy(res, CYS.SG, CB, 1.82, CA, 109.0, N, 180, 0, ff, bonds);
    buildH(res, CYS.HB2, CB, 1.11, CA, 109.4, SG, 112.0, 1, ff, bonds);
    buildH(res, CYS.HB3, CB, 1.11, CA, 109.4, SG, 112.0, -1, ff, bonds);
    buildH(res, CYS.HG, SG, 1.34, CB, 96.0, CA, 180.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildCystine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildCystine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, CYX.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom SG = buildHeavy(res, CYX.SG, CB, 1.82, CA, 109.0, N, 180, 0, ff, bonds);
    buildH(res, CYX.HB2, CB, 1.11, CA, 109.4, SG, 112.0, 1, ff, bonds);
    buildH(res, CYX.HB3, CB, 1.11, CA, 109.4, SG, 112.0, -1, ff, bonds);
    List<Atom> resAtoms = res.getAtomList();
    for (Atom atom : resAtoms) {
      atom.setResName("CYS");
    }
    res.setName("CYS");
    return res;
  }

  /**
   * buildDeprotonatedCysteine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildDeprotonatedCysteine(Residue res, Atom CA, Atom N, Atom C,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, CYD.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom SG = buildHeavy(res, CYD.SG, CB, 1.82, CA, 109.0, N, 180, 0, ff, bonds);
    buildH(res, CYD.HB2, CB, 1.11, CA, 109.4, SG, 112.0, 1, ff, bonds);
    buildH(res, CYD.HB3, CB, 1.11, CA, 109.4, SG, 112.0, -1, ff, bonds);
    return res;
  }

  /**
   * buildDeprotonatedLysine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildDeprotonatedLysine(Residue res, Atom CA, Atom N, Atom C,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, LYD.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, LYD.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, LYD.CD, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    Atom CE = buildHeavy(res, LYD.CE, CD, 1.54, CG, 109.5, CB, 180, 0, ff, bonds);
    Atom NZ = buildHeavy(res, LYD.NZ, CE, 1.50, CD, 109.5, CG, 180, 0, ff, bonds);
    buildH(res, LYD.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, LYD.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, LYD.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, LYD.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    buildH(res, LYD.HD2, CD, 1.11, CG, 109.4, CE, 109.4, 1, ff, bonds);
    buildH(res, LYD.HD3, CD, 1.11, CG, 109.4, CE, 109.4, -1, ff, bonds);
    buildH(res, LYD.HE2, CE, 1.11, CD, 109.4, NZ, 108.8, 1, ff, bonds);
    buildH(res, LYD.HE3, CE, 1.11, CD, 109.4, NZ, 108.8, -1, ff, bonds);
    Atom HZ1 = buildH(res, LYD.HZ1, NZ, 1.02, CE, 109.5, CD, 180.0, 0, ff, bonds);
    buildH(res, LYD.HZ2, NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, ff, bonds);
    return res;
  }

  /**
   * buildDeprotonatedTyrosine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildDeprotonatedTyrosine(Residue res, Atom CA, Atom N, Atom C,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, TYD.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, TYD.CG, CB, 1.50, CA, 109.5, N, 62, 0, ff, bonds);
    Atom CD1 = buildHeavy(res, TYD.CD1, CG, 1.39, CB, 120.0, CA, 90, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, TYD.CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, TYD.CE1, CD1, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CE2 = buildHeavy(res, TYD.CE2, CD2, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CZ = buildHeavy(res, TYD.CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0, ff, bonds);
    buildBond(CE2, CZ, ff, bonds);
    buildHeavy(res, TYD.OH, CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, ff, bonds);
    buildH(res, TYD.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, TYD.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, TYD.HD1, CD1, 1.10, CG, 120.0, CE1, 120.0, 1, ff, bonds);
    buildH(res, TYD.HD2, CD2, 1.10, CG, 120.0, CE2, 120.0, 1, ff, bonds);
    buildH(res, TYD.HE1, CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, ff, bonds);
    buildH(res, TYD.HE2, CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildGlutamate.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildGlutamate(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, GLU.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, GLU.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, GLU.CD, CG, 1.51, CB, 107.8, CA, 180, 0, ff, bonds);
    Atom OE1 = buildHeavy(res, GLU.OE1, CD, 1.25, CG, 117.0, CB, 180, 0, ff, bonds);
    buildHeavy(res, GLU.OE2, CD, 1.25, CG, 117.0, OE1, 126.0, 1, ff, bonds);
    buildH(res, GLU.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, GLU.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, GLU.HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1, ff, bonds);
    buildH(res, GLU.HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1, ff, bonds);
    return res;
  }

  /**
   * buildGlutamine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildGlutamine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, GLN.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, GLN.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, GLN.CD, CG, 1.51, CB, 107.8, CA, 180, 0, ff, bonds);
    Atom OE1 = buildHeavy(res, GLN.OE1, CD, 1.22, CG, 122.5, CB, 180, 0, ff, bonds);
    Atom NE2 = buildHeavy(res, GLN.NE2, CD, 1.34, CG, 112.7, OE1, 124.0, 1, ff, bonds);
    buildH(res, GLN.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, GLN.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, GLN.HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1, ff, bonds);
    buildH(res, GLN.HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1, ff, bonds);
    Atom HE21 = buildH(res, GLN.HE21, NE2, 1.02, CD, 119.0, CG, 0.0, 0, ff, bonds);
    buildH(res, GLN.HE22, NE2, 1.02, CD, 119.0, HE21, 120.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildHistidine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildHistidine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, HIS.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bonds);
    Atom CG = buildHeavy(res, HIS.CG, CB, 1.50, CA, 109.5, N, 180, 0, ff, bonds);
    Atom ND1 = buildHeavy(res, HIS.ND1, CG, 1.35, CB, 126.0, CA, 180, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, HIS.CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, HIS.CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0, ff, bonds);
    Atom NE2 = buildHeavy(res, HIS.NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0, ff, bonds);
    buildBond(NE2, CE1, ff, bonds);
    buildH(res, HIS.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, HIS.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, HIS.HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0, ff, bonds);
    buildH(res, HIS.HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1, ff, bonds);
    buildH(res, HIS.HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, ff, bonds);
    buildH(res, HIS.HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildIsoleucine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildIsoleucine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ILE.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bonds);
    Atom CG1 = buildHeavy(res, ILE.CG1, CB, 1.54, CA, 109.5, N, 0, 0, ff, bonds);
    Atom CG2 = buildHeavy(res, ILE.CG2, CB, 1.54, CA, 109.5, CG1, 109.5, 1, ff, bonds);
    Atom CD1 = buildHeavy(res, ILE.CD1, CG1, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    buildH(res, ILE.HB, CB, 1.11, CA, 109.4, CG2, 109.4, 1, ff, bonds);
    buildH(res, ILE.HG12, CG1, 1.11, CB, 109.4, CD1, 109.4, 1, ff, bonds);
    buildH(res, ILE.HG13, CG1, 1.11, CB, 109.4, CD1, 109.4, -1, ff, bonds);
    Atom HG21 = buildH(res, ILE.HG21, CG2, 1.11, CB, 110.0, CG1, 180.0, 0, ff, bonds);
    buildH(res, ILE.HG22, CG2, 1.11, CB, 110.0, HG21, 109.0, 1, ff, bonds);
    buildH(res, ILE.HG23, CG2, 1.11, CB, 110.0, HG21, 109.0, -1, ff, bonds);
    Atom HD11 = buildH(res, ILE.HD11, CD1, 1.11, CG1, 110.0, CB, 180.0, 0, ff, bonds);
    buildH(res, ILE.HD12, CD1, 1.11, CG1, 110.0, HD11, 109.0, 1, ff, bonds);
    buildH(res, ILE.HD13, CD1, 1.11, CG1, 110.0, HD11, 109.0, -1, ff, bonds);
    return res;
  }

  /**
   * buildLeucine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildLeucine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, LEU.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, LEU.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD1 = buildHeavy(res, LEU.CD1, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, LEU.CD2, CG, 1.54, CB, 109.5, CD1, 109.5, -1, ff, bonds);
    buildH(res, LEU.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, LEU.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, LEU.HG, CG, 1.11, CB, 109.4, CD1, 109.4, 1, ff, bonds);
    Atom HD11 = buildH(res, LEU.HD11, CD1, 1.11, CG, 109.4, CB, 180.0, 0, ff, bonds);
    buildH(res, LEU.HD12, CD1, 1.11, CG, 109.4, HD11, 109.4, 1, ff, bonds);
    buildH(res, LEU.HD13, CD1, 1.11, CG, 109.4, HD11, 109.4, -1, ff, bonds);
    Atom HD21 = buildH(res, LEU.HD21, CD2, 1.11, CG, 109.4, CB, 180.0, 0, ff, bonds);
    buildH(res, LEU.HD22, CD2, 1.11, CG, 109.4, HD21, 109.4, 1, ff, bonds);
    buildH(res, LEU.HD23, CD2, 1.11, CG, 109.4, HD21, 109.4, -1, ff, bonds);
    return res;
  }

  /**
   * buildLysine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildLysine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, LYS.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, LYS.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, LYS.CD, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    Atom CE = buildHeavy(res, LYS.CE, CD, 1.54, CG, 109.5, CB, 180, 0, ff, bonds);
    Atom NZ = buildHeavy(res, LYS.NZ, CE, 1.50, CD, 109.5, CG, 180, 0, ff, bonds);
    buildH(res, LYS.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, LYS.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, LYS.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, LYS.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    buildH(res, LYS.HD2, CD, 1.11, CG, 109.4, CE, 109.4, 1, ff, bonds);
    buildH(res, LYS.HD3, CD, 1.11, CG, 109.4, CE, 109.4, -1, ff, bonds);
    buildH(res, LYS.HE2, CE, 1.11, CD, 109.4, NZ, 108.8, 1, ff, bonds);
    buildH(res, LYS.HE3, CE, 1.11, CD, 109.4, NZ, 108.8, -1, ff, bonds);
    Atom HZ1 = buildH(res, LYS.HZ1, NZ, 1.02, CE, 109.5, CD, 180.0, 0, ff, bonds);
    buildH(res, LYS.HZ2, NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, ff, bonds);
    buildH(res, LYS.HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1, ff, bonds);
    return res;
  }

  /**
   * buildMethionine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bond a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildMethionine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bond) {
    Atom CB = buildHeavy(res, MET.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bond);
    Atom CG = buildHeavy(res, MET.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bond);
    Atom SD = buildHeavy(res, MET.SD, CG, 1.82, CB, 109.0, CA, 180, 0, ff, bond);
    Atom CE = buildHeavy(res, MET.CE, SD, 1.82, CG, 96.3, CB, 180, 0, ff, bond);
    buildH(res, MET.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bond);
    buildH(res, MET.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bond);
    buildH(res, MET.HG2, CG, 1.11, CB, 109.4, SD, 112.0, 1, ff, bond);
    buildH(res, MET.HG3, CG, 1.11, CB, 109.4, SD, 112.0, -1, ff, bond);
    Atom HE1 = buildH(res, MET.HE1, CE, 1.11, SD, 112.0, CG, 180.0, 0, ff, bond);
    buildH(res, MET.HE2, CE, 1.11, SD, 112.0, HE1, 109.4, 1, ff, bond);
    buildH(res, MET.HE3, CE, 1.11, SD, 112.0, HE1, 109.4, -1, ff, bond);
    return res;
  }

  /**
   * buildNeutralAsparticAcid.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildNeutralAsparticAcid(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ASH.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ASH.CG, CB, 1.51, CA, 107.8, N, 180, 0, ff, bonds);
    Atom OD1 = buildHeavy(res, ASH.OD1, CG, 1.25, CB, 117.0, CA, 0.0, 0, ff, bonds);
    Atom OD2 = buildHeavy(res, ASH.OD2, CG, 1.25, CB, 117.0, OD1, 126.0, 1, ff, bonds);
    Atom HB2 = buildH(res, ASH.HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1, ff, bonds);
    Atom HB3 = buildH(res, ASH.HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1, ff, bonds);
    Atom HD2 = buildH(res, ASH.HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildTwoProtonAsparticAcid.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildTwoProtonAsparticAcid(Residue res, Atom CA, Atom N, Atom C,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, ASD.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ASD.CG, CB, 1.51, CA, 107.8, N, 180, 0, ff, bonds);
    Atom OD1 = buildHeavy(res, ASD.OD1, CG, 1.25, CB, 117.0, CA, 0.0, 0, ff, bonds);
    Atom OD2 = buildHeavy(res, ASD.OD2, CG, 1.25, CB, 117.0, OD1, 126.0, 1, ff, bonds);
    buildH(res, ASD.HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1, ff, bonds);
    buildH(res, ASD.HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1, ff, bonds);
    buildH(res, ASD.HD1, OD1, 0.98, CG, 108.7, OD2, 0.0, 0, ff, bonds);
    buildH(res, ASD.HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildNeutralGlutamicAcid.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildNeutralGlutamicAcid(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, GLH.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, GLH.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, GLH.CD, CG, 1.51, CB, 107.8, CA, 180, 0, ff, bonds);
    Atom OE1 = buildHeavy(res, GLH.OE1, CD, 1.25, CG, 117.0, CB, 180, 0, ff, bonds);
    Atom OE2 = buildHeavy(res, GLH.OE2, CD, 1.25, CG, 117.0, OE1, 126.0, 1, ff, bonds);
    buildH(res, GLH.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, GLH.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, GLH.HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1, ff, bonds);
    buildH(res, GLH.HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1, ff, bonds);
    buildH(res, GLH.HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildTwoProtonGlutamicAcid.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildTwoProtonGlutamicAcid(Residue res, Atom CA, Atom N, Atom C,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, GLD.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, GLD.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, GLD.CD, CG, 1.51, CB, 107.8, CA, 180, 0, ff, bonds);
    Atom OE1 = buildHeavy(res, GLD.OE1, CD, 1.25, CG, 117.0, CB, 180, 0, ff, bonds);
    Atom OE2 = buildHeavy(res, GLD.OE2, CD, 1.25, CG, 117.0, OE1, 126.0, 1, ff, bonds);
    buildH(res, GLD.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, GLD.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, GLD.HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1, ff, bonds);
    buildH(res, GLD.HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1, ff, bonds);
    buildH(res, GLD.HE1, OE1, 0.98, CD, 108.7, OE2, 0.0, 0, ff, bonds);
    buildH(res, GLD.HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildNeutralHistidineD.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildNeutralHistidineD(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, HID.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bonds);
    Atom CG = buildHeavy(res, HID.CG, CB, 1.50, CA, 109.5, N, 180, 0, ff, bonds);
    Atom ND1 = buildHeavy(res, HID.ND1, CG, 1.35, CB, 126.0, CA, 180, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, HID.CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, HID.CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0, ff, bonds);
    Atom NE2 = buildHeavy(res, HID.NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0, ff, bonds);
    buildBond(NE2, CE1, ff, bonds);
    buildH(res, HID.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, HID.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, HID.HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0, ff, bonds);
    buildH(res, HID.HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1, ff, bonds);
    buildH(res, HID.HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildNeutralHistidineE.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildNeutralHistidineE(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, HIE.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bonds);
    Atom CG = buildHeavy(res, HIE.CG, CB, 1.50, CA, 109.5, N, 180, 0, ff, bonds);
    Atom ND1 = buildHeavy(res, HIE.ND1, CG, 1.35, CB, 126.0, CA, 180, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, HIE.CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, HIE.CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0, ff, bonds);
    Atom NE2 = buildHeavy(res, HIE.NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0, ff, bonds);
    buildBond(NE2, CE1, ff, bonds);
    buildH(res, HIE.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, HIE.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, HIE.HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1, ff, bonds);
    buildH(res, HIE.HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, ff, bonds);
    buildH(res, HIE.HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildOrnithine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildOrnithine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, ORN.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, ORN.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, ORN.CD, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    Atom NE = buildHeavy(res, ORN.NE, CD, 1.50, CG, 109.5, CB, 180, 0, ff, bonds);
    buildH(res, ORN.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, ORN.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, ORN.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, ORN.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    buildH(res, ORN.HD2, CD, 1.11, CG, 109.4, NE, 109.4, 1, ff, bonds);
    buildH(res, ORN.HD3, CD, 1.11, CG, 109.4, NE, 109.4, -1, ff, bonds);
    Atom HE1 = buildH(res, ORN.HE1, NE, 1.02, CD, 109.5, CG, 180.0, 0, ff, bonds);
    buildH(res, ORN.HE2, NE, 1.02, CD, 109.5, HE1, 109.5, 1, ff, bonds);
    buildH(res, ORN.HE3, NE, 1.02, CD, 109.5, HE1, 109.5, -1, ff, bonds);
    return res;
  }

  /**
   * buildPCA.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildPCA(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, PCA.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, PCA.CG, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD = buildHeavy(res, PCA.CD, CG, 1.54, CB, 109.5, CA, 180, 0, ff, bonds);
    buildHeavy(res, PCA.OE, CD, 1.22, N, 126.0, CG, 126.0, 1, ff, bonds);
    buildH(res, PCA.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, PCA.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, PCA.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, PCA.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    return res;
  }

  /**
   * buildPhenylalanine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildPhenylalanine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, PHE.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, PHE.CG, CB, 1.50, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CD1 = buildHeavy(res, PHE.CD1, CG, 1.39, CB, 120.0, CA, 180, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, PHE.CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, PHE.CE1, CD1, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CE2 = buildHeavy(res, PHE.CE2, CD2, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CZ = buildHeavy(res, PHE.CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0, ff, bonds);
    buildBond(CE2, CZ, ff, bonds);
    buildH(res, PHE.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, PHE.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, PHE.HD1, CD1, 1.11, CG, 120.0, CE1, 120.0, 1, ff, bonds);
    buildH(res, PHE.HD2, CD2, 1.11, CG, 120.0, CE2, 120.0, 1, ff, bonds);
    buildH(res, PHE.HE1, CE1, 1.11, CD1, 120.0, CZ, 120.0, 1, ff, bonds);
    buildH(res, PHE.HE2, CE2, 1.11, CD2, 120.0, CZ, 120.0, 1, ff, bonds);
    buildH(res, PHE.HZ, CZ, 1.11, CE1, 120.0, CE2, 120.0, 1, ff, bonds);
    return res;
  }

  /**
   * buildProline.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param position a {@link ResiduePosition} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildProline(Residue res, Atom CA, Atom N, Atom C, ResiduePosition position,
      ForceField ff, List<Bond> bonds) {
    Atom CB = buildHeavy(res, PRO.CB, CA, 1.5247, N, 104.0, C, 109.5, 1, ff, bonds);
    Atom CG = buildHeavy(res, PRO.CG, CB, 1.5247, CA, 104.0, N, 30.0, 0, ff, bonds);
    int cdKey = position == FIRST_RESIDUE ? 469 : PRO.CD.getType();
    Atom CD = buildHeavy(res, PRO.CD.name(), N, 1.5247, CA, 104.0, CB, 0, 0, cdKey, ff, bonds);
    buildBond(CD, CG, ff, bonds);
    buildH(res, PRO.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, PRO.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, PRO.HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1, ff, bonds);
    buildH(res, PRO.HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1, ff, bonds);
    if (position == FIRST_RESIDUE) {
      buildH(res, PRO.HD2.name(), CD, 1.11, CG, 109.4, N, 109.4, 1, 470, ff, bonds);
      buildH(res, PRO.HD3.name(), CD, 1.11, CG, 109.4, N, 109.4, -1, 470, ff, bonds);
    } else {
      buildH(res, PRO.HD2, CD, 1.11, CG, 109.4, N, 109.4, 1, ff, bonds);
      buildH(res, PRO.HD3, CD, 1.11, CG, 109.4, N, 109.4, -1, ff, bonds);
    }
    return res;
  }

  /**
   * buildSerine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildSerine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, SER.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom OG = buildHeavy(res, SER.OG, CB, 1.41, CA, 107.5, N, 180, 0, ff, bonds);
    buildH(res, SER.HB2, CB, 1.11, CA, 109.4, OG, 106.7, 1, ff, bonds);
    buildH(res, SER.HB3, CB, 1.11, CA, 109.4, OG, 106.7, -1, ff, bonds);
    buildH(res, SER.HG, OG, 0.94, CB, 106.9, CA, 180.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildThreonine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildThreonine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, THR.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bonds);
    Atom OG1 = buildHeavy(res, THR.OG1, CB, 1.41, CA, 107.5, N, 180, 0, ff, bonds);
    Atom CG2 = buildHeavy(res, THR.CG2, CB, 1.54, CA, 109.5, OG1, 107.7, 1, ff, bonds);
    buildH(res, THR.HB, CB, 1.11, CA, 109.4, OG1, 106.7, -1, ff, bonds);
    buildH(res, THR.HG1, OG1, 0.94, CB, 106.9, CA, 180.0, 0, ff, bonds);
    Atom HG21 = buildH(res, THR.HG21, CG2, 1.11, CB, 110.0, CA, 180.0, 0, ff, bonds);
    buildH(res, THR.HG22, CG2, 1.11, CB, 110.0, HG21, 109.0, 1, ff, bonds);
    buildH(res, THR.HG23, CG2, 1.11, CB, 110.0, HG21, 109.0, -1, ff, bonds);
    return res;
  }

  /**
   * buildTryptophan.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bond a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildTryptophan(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bond) {
    Atom CB = buildHeavy(res, TRP.CB, CA, 1.54, N, 109.5, C, 109.5, 1, ff, bond);
    Atom CG = buildHeavy(res, TRP.CG, CB, 1.50, CA, 109.5, N, 62, 0, ff, bond);
    Atom CD1 = buildHeavy(res, TRP.CD1, CG, 1.35, CB, 126.0, CA, -90, 0, ff, bond);
    Atom CD2 = buildHeavy(res, TRP.CD2, CG, 1.35, CB, 126.0, CD1, 108.0, 1, ff, bond);
    Atom NE1 = buildHeavy(res, TRP.NE1, CD1, 1.35, CG, 108.0, CD2, 0.0, 0, ff, bond);
    Atom CE2 = buildHeavy(res, TRP.CE2, NE1, 1.35, CD1, 108.0, CG, 0.0, 0, ff, bond);
    buildBond(CE2, CD2, ff, bond);
    Atom CE3 = buildHeavy(res, TRP.CE3, CD2, 1.35, CE2, 120.0, NE1, 180.0, 0, ff, bond);
    Atom CZ2 = buildHeavy(res, TRP.CZ2, CE2, 1.35, CD2, 120.0, CE3, 0.0, 0, ff, bond);
    Atom CZ3 = buildHeavy(res, TRP.CZ3, CE3, 1.35, CD2, 120.0, CE2, 0.0, 0, ff, bond);
    Atom CH2 = buildHeavy(res, TRP.CH2, CZ2, 1.35, CE2, 120.0, CD2, 0.0, 0, ff, bond);
    buildBond(CH2, CZ3, ff, bond);
    buildH(res, TRP.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bond);
    buildH(res, TRP.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bond);
    buildH(res, TRP.HD1, CD1, 1.10, CG, 126.0, NE1, 126.0, 1, ff, bond);
    buildH(res, TRP.HE1, NE1, 1.05, CD1, 126.0, CE2, 126.0, 1, ff, bond);
    buildH(res, TRP.HE3, CE3, 1.10, CD1, 120.0, CZ3, 120.0, 1, ff, bond);
    buildH(res, TRP.HZ2, CZ2, 1.10, CE2, 120.0, CH2, 120.0, 1, ff, bond);
    buildH(res, TRP.HZ3, CZ3, 1.10, CE3, 120.0, CH2, 120.0, 1, ff, bond);
    buildH(res, TRP.HH2, CH2, 1.10, CZ2, 120.0, CZ3, 120.0, 1, ff, bond);
    return res;
  }

  /**
   * buildTyrosine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildTyrosine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, TYR.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG = buildHeavy(res, TYR.CG, CB, 1.50, CA, 109.5, N, 62, 0, ff, bonds);
    Atom CD1 = buildHeavy(res, TYR.CD1, CG, 1.39, CB, 120.0, CA, 90, 0, ff, bonds);
    Atom CD2 = buildHeavy(res, TYR.CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1, ff, bonds);
    Atom CE1 = buildHeavy(res, TYR.CE1, CD1, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CE2 = buildHeavy(res, TYR.CE2, CD2, 1.39, CG, 120.0, CB, 180, 0, ff, bonds);
    Atom CZ = buildHeavy(res, TYR.CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0, ff, bonds);
    buildBond(CE2, CZ, ff, bonds);
    Atom OH = buildHeavy(res, TYR.OH, CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, ff, bonds);
    buildH(res, TYR.HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1, ff, bonds);
    buildH(res, TYR.HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1, ff, bonds);
    buildH(res, TYR.HD1, CD1, 1.10, CG, 120.0, CE1, 120.0, 1, ff, bonds);
    buildH(res, TYR.HD2, CD2, 1.10, CG, 120.0, CE2, 120.0, 1, ff, bonds);
    buildH(res, TYR.HE1, CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, ff, bonds);
    buildH(res, TYR.HE2, CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, ff, bonds);
    buildH(res, TYR.HH, OH, 0.97, CZ, 108.0, CE2, 0.0, 0, ff, bonds);
    return res;
  }

  /**
   * buildValine.
   *
   * @param res a {@link Residue} object.
   * @param CA a {@link Atom} object.
   * @param N a {@link Atom} object.
   * @param C a {@link Atom} object.
   * @param ff a {@link ForceField} object.
   * @param bonds a {@link List} object.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public static Residue buildValine(Residue res, Atom CA, Atom N, Atom C, ForceField ff,
      List<Bond> bonds) {
    Atom CB = buildHeavy(res, VAL.CB, CA, 1.54, N, 109.5, C, 107.8, 1, ff, bonds);
    Atom CG1 = buildHeavy(res, VAL.CG1, CB, 1.54, CA, 109.5, N, 180, 0, ff, bonds);
    Atom CG2 = buildHeavy(res, VAL.CG2, CB, 1.54, CA, 109.5, CG1, 109.5, -1, ff, bonds);
    buildH(res, VAL.HB, CB, 1.11, CA, 109.4, CG1, 109.4, 1, ff, bonds);
    Atom HG11 = buildH(res, VAL.HG11, CG1, 1.11, CB, 109.4, CA, 180.0, 0, ff, bonds);
    buildH(res, VAL.HG12, CG1, 1.11, CB, 109.4, HG11, 109.4, 1, ff, bonds);
    buildH(res, VAL.HG13, CG1, 1.11, CB, 109.4, HG11, 109.4, -1, ff, bonds);
    Atom HG21 = buildH(res, VAL.HG21, CG2, 1.11, CB, 109.4, CA, 180.0, 0, ff, bonds);
    buildH(res, VAL.HG22, CG2, 1.11, CB, 109.4, HG21, 109.4, 1, ff, bonds);
    buildH(res, VAL.HG23, CG2, 1.11, CB, 109.4, HG21, 109.4, -1, ff, bonds);
    return res;
  }

  /**
   * copyResidue.
   *
   * @param fromResidue a {@link ffx.potential.bonded.Residue} object.
   * @param toResidue a {@link ffx.potential.bonded.Residue} object.
   */
  public static void copyResidue(Residue fromResidue, Residue toResidue) {
    String resName = fromResidue.getName();
    AminoAcid3 res = AminoAcid3.valueOf(resName);
    List<String> atomNames = new ArrayList<>();
    switch (res) {
      case ALA:
        atomNames.addAll(Arrays.asList(getNames(ALA.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ASD:
        atomNames.addAll(Arrays.asList(getNames(ASD.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ASH:
        atomNames.addAll(Arrays.asList(getNames(ASH.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ASN:
        atomNames.addAll(Arrays.asList(getNames(ASN.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ASP:
        atomNames.addAll(Arrays.asList(getNames(ASP.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ARG:
        atomNames.addAll(Arrays.asList(getNames(ARG.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case CYS:
        atomNames.addAll(Arrays.asList(getNames(CYS.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case CYD:
        atomNames.addAll(Arrays.asList(getNames(CYD.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case CYX:
        atomNames.addAll(Arrays.asList(getNames(CYX.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case GLD:
        atomNames.addAll(Arrays.asList(getNames(GLD.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case GLH:
        atomNames.addAll(Arrays.asList(getNames(GLH.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case GLN:
        atomNames.addAll(Arrays.asList(getNames(GLN.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case GLU:
        atomNames.addAll(Arrays.asList(getNames(GLU.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case GLY:
        atomNames.addAll(Arrays.asList(getNames(GLY.class)));
        atomNames.addAll(Arrays.asList(getNames(GlycineBackboneAtoms.class)));
        break;
      case HID:
        atomNames.addAll(Arrays.asList(getNames(HID.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case HIE:
        atomNames.addAll(Arrays.asList(getNames(HIE.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case HIS:
        atomNames.addAll(Arrays.asList(getNames(HIS.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case ILE:
        atomNames.addAll(Arrays.asList(getNames(ILE.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case LEU:
        atomNames.addAll(Arrays.asList(getNames(LEU.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case LYD:
        atomNames.addAll(Arrays.asList(getNames(LYD.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case LYS:
        atomNames.addAll(Arrays.asList(getNames(LYS.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case MET:
        atomNames.addAll(Arrays.asList(getNames(MET.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case PHE:
        atomNames.addAll(Arrays.asList(getNames(PHE.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case PRO:
        atomNames.addAll(Arrays.asList(getNames(PRO.class)));
        atomNames.addAll(Arrays.asList(getNames(ProlineBackboneAtoms.class)));
        break;
      case SER:
        atomNames.addAll(Arrays.asList(getNames(SER.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case THR:
        atomNames.addAll(Arrays.asList(getNames(THR.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case TRP:
        atomNames.addAll(Arrays.asList(getNames(TRP.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case TYD:
        atomNames.addAll(Arrays.asList(getNames(TYD.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case TYR:
        atomNames.addAll(Arrays.asList(getNames(TYR.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      case VAL:
        atomNames.addAll(Arrays.asList(getNames(VAL.class)));
        atomNames.addAll(Arrays.asList(getNames(AminoAcidBackboneAtoms.class)));
        break;
      default:
        atomNames = null;
    }
    for (String atomName : atomNames) {
      copyCoordinates(fromResidue, toResidue, atomName);
    }
  }

  /**
   * Only the first nitrogen should have H1, H2 and H3 atoms, unless it's an NME cap.
   *
   * @param aminoAcid 3-letter amino acid name.
   * @param residue the amino acid Residue.
   */
  public static void removeH1_H2_H3(AminoAcid3 aminoAcid, Residue residue) {
    if (aminoAcid != AminoAcid3.NME) {
      if (aminoAcid != AminoAcid3.NH2) {
        Atom H1 = (Atom) residue.getAtomNode("H1");
        if (H1 != null) {
          residue.deleteAtom(H1);
        }
        Atom H2 = (Atom) residue.getAtomNode("H2");
        if (H2 != null) {
          residue.deleteAtom(H2);
        }
      }
      Atom H3 = (Atom) residue.getAtomNode("H3");
      if (H3 != null) {
        residue.deleteAtom(H3);
      }
    }
  }

  /**
   * Only the last residue in a chain should have an OXT/OT2 atom.
   *
   * @param residue the amino acid residue.
   */
  public static void removeOXT_OT2(Residue residue) {
    Atom OXT = (Atom) residue.getAtomNode("OXT");
    if (OXT != null) {
      residue.deleteAtom(OXT);
    }
    Atom OT2 = (Atom) residue.getAtomNode("OT2");
    if (OT2 != null) {
      residue.deleteAtom(OT2);
    }
  }

  /**
   * Assign atom types to a single amino acid side chain.
   *
   * @param position The position of this amino acid in the chain.
   * @param aminoAcid The amino acid to use.
   * @param residue The residue node.
   * @param CA The C-alpha carbon of this residue.
   * @param N The peptide nitrogen of this residue.
   * @param C The peptide carbonyl carbon.
   * @param forceField a {@link ForceField} object.
   * @param bondList a {@link java.util.List} object.
   */
  private static void assignAminoAcidSideChain(
      ResiduePosition position,
      AminoAcid3 aminoAcid,
      Residue residue,
      Atom CA,
      Atom N,
      Atom C,
      ForceField forceField,
      List<Bond> bondList) {
    switch (aminoAcid) {
      case ALA:
        buildAlanine(residue, CA, N, C, forceField, bondList);
        break;
      case GLY:
        buildGlycine(residue, CA, N, C, position, forceField, bondList);
        break;
      case VAL:
        buildValine(residue, CA, N, C, forceField, bondList);
        break;
      case LEU:
        buildLeucine(residue, CA, N, C, forceField, bondList);
        break;
      case ILE:
        buildIsoleucine(residue, CA, N, C, forceField, bondList);
        break;
      case SER:
        buildSerine(residue, CA, N, C, forceField, bondList);
        break;
      case THR:
        buildThreonine(residue, CA, N, C, forceField, bondList);
        break;
      case CYS:
        buildCysteine(residue, CA, N, C, forceField, bondList);
        break;
      case CYX:
        buildCystine(residue, CA, N, C, forceField, bondList);
        break;
      case CYD:
        buildDeprotonatedCysteine(residue, CA, N, C, forceField, bondList);
        break;
      case PRO:
        buildProline(residue, CA, N, C, position, forceField, bondList);
        break;
      case PHE:
        buildPhenylalanine(residue, CA, N, C, forceField, bondList);
        break;
      case TYR:
        buildTyrosine(residue, CA, N, C, forceField, bondList);
        break;
      case TYD:
        buildDeprotonatedTyrosine(residue, CA, N, C, forceField, bondList);
        break;
      case TRP:
        buildTryptophan(residue, CA, N, C, forceField, bondList);
        break;
      case HIS:
        buildHistidine(residue, CA, N, C, forceField, bondList);
        break;
      case HID:
        buildNeutralHistidineD(residue, CA, N, C, forceField, bondList);
        break;
      case HIE:
        buildNeutralHistidineE(residue, CA, N, C, forceField, bondList);
        break;
      case ASP:
        buildAspartate(residue, CA, N, C, forceField, bondList);
        break;
      case ASH:
        buildNeutralAsparticAcid(residue, CA, N, C, forceField, bondList);
        break;
      case ASD:
        buildTwoProtonAsparticAcid(residue, CA, N, C, forceField, bondList);
        break;
      case ASN:
        buildAsparagine(residue, CA, N, C, forceField, bondList);
        break;
      case GLU:
        buildGlutamate(residue, CA, N, C, forceField, bondList);
        break;
      case GLH:
        buildNeutralGlutamicAcid(residue, CA, N, C, forceField, bondList);
        break;
      case GLD:
        buildTwoProtonGlutamicAcid(residue, CA, N, C, forceField, bondList);
        break;
      case GLN:
        buildGlutamine(residue, CA, N, C, forceField, bondList);
        break;
      case MET:
        buildMethionine(residue, CA, N, C, forceField, bondList);
        break;
      case LYS:
        buildLysine(residue, CA, N, C, forceField, bondList);
        break;
      case LYD:
        buildDeprotonatedLysine(residue, CA, N, C, forceField, bondList);
        break;
      case ARG:
        buildArginine(residue, CA, N, C, forceField, bondList);
        break;
      case ORN:
        buildOrnithine(residue, CA, N, C, forceField, bondList);
        break;
      case AIB:
        buildAIB(residue, CA, N, C, forceField, bondList);
        break;
      case PCA:
        buildPCA(residue, CA, N, C, forceField, bondList);
        break;
      case UNK:
        String residueName = residue.getName();
        logger.log(Level.INFO, " Patching side-chain {0}", residueName);
        HashMap<String, AtomType> types = forceField.getAtomTypes(residueName);
        if (!types.isEmpty()) {
          boolean patched = true;
          List<Atom> residueAtoms = residue.getAtomList();
          // Assign atom types for side-chain atoms.
          for (Atom atom : residueAtoms) {
            String atomName = atom.getName().toUpperCase();
            AtomType type = atom.getAtomType();
            if (type == null) {
              type = types.get(atomName);
              atom.setAtomType(type);
              types.remove(atomName);
            }
          }
          // Create bonds between known atoms.
          for (Atom atom : residueAtoms) {
            String atomName = atom.getName();
            String[] bonds = forceField.getBonds(residueName, atomName);
            if (bonds != null) {
              for (String name : bonds) {
                Atom atom2 = (Atom) residue.getAtomNode(name);
                if (atom2 != null && !atom.isBonded(atom2)) {
                  buildBond(atom, atom2, forceField, bondList);
                }
              }
            }
          }

          // Create missing hydrogen atoms.
          if (!types.isEmpty()) {
            // Create a hashmap of the molecule's atoms
            HashMap<String, Atom> atomMap = new HashMap<>();
            for (Atom atom : residueAtoms) {
              atomMap.put(atom.getName().toUpperCase(), atom);
            }
            for (String atomName : types.keySet()) {
              AtomType type = types.get(atomName);
              String[] bonds = forceField.getBonds(residueName, atomName.toUpperCase());
              if (bonds == null || bonds.length != 1) {
                patched = false;
                logger.log(Level.INFO, " Check biotype for hydrogen {0}.", type.name);
                break;
              }
              // Get the heavy atom the hydrogen is bonded to.
              Atom ia = atomMap.get(bonds[0].toUpperCase());
              Atom hydrogen =
                  new Atom(
                      0,
                      atomName,
                      ia.getAltLoc(),
                      new double[3],
                      ia.getResidueName(),
                      ia.getResidueNumber(),
                      ia.getChainID(),
                      ia.getOccupancy(),
                      ia.getTempFactor(),
                      ia.getSegID());
              logger.log(Level.FINE, " Created hydrogen {0}.", atomName);
              hydrogen.setAtomType(type);
              hydrogen.setHetero(true);
              residue.addMSNode(hydrogen);
              int valence = ia.getAtomType().valence;
              List<Bond> aBonds = ia.getBonds();
              int numBonds = aBonds.size();
              // Try to find the following configuration: ib-ia-ic
              Atom ib = null;
              Atom ic = null;
              Atom id = null;
              if (numBonds > 0) {
                Bond bond = aBonds.get(0);
                ib = bond.get1_2(ia);
              }
              if (numBonds > 1) {
                Bond bond = aBonds.get(1);
                ic = bond.get1_2(ia);
              }
              if (numBonds > 2) {
                Bond bond = aBonds.get(2);
                id = bond.get1_2(ia);
              }

              // Building the hydrogens depends on hybridization and the locations of other bonded
              // atoms.
              logger.log(
                  Level.FINE,
                  " Bonding {0} to {1} ({2} of {3}).",
                  new Object[] {atomName, ia.getName(), numBonds, valence});
              switch (valence) {
                case 4:
                  switch (numBonds) {
                    case 3:
                      // Find the average coordinates of atoms ib, ic and id.
                      Double3 b = ib.getXYZ();
                      Double3 c = ic.getXYZ();
                      Double3 d = id.getXYZ();
                      Double3 a = b.add(c).addI(d).scaleI(1.0 / 3.0);

                      // Place the hydrogen at chiral position #1.
                      intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 0);
                      Double3 e1 = hydrogen.getXYZ();
                      Double3 ret = a.sub(e1);
                      double l1 = ret.length();

                      // Place the hydrogen at chiral position #2.
                      intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 1);
                      Double3 e2 = hydrogen.getXYZ();
                      ret = a.sub(e2);
                      double l2 = ret.length();

                      // Revert to #1 if it is farther from the average.
                      if (l1 > l2) {
                        hydrogen.setXYZ(e1.get());
                      }
                      break;
                    case 2:
                      intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 0);
                      break;
                    case 1:
                      intxyz(hydrogen, ia, 1.0, ib, 109.5, null, 0.0, 0);
                      break;
                    case 0:
                      intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                      break;
                    default:
                      logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                      patched = false;
                  }
                  break;
                case 3:
                  switch (numBonds) {
                    case 2:
                      intxyz(hydrogen, ia, 1.0, ib, 120.0, ic, 180.0, 0);
                      break;
                    case 1:
                      intxyz(hydrogen, ia, 1.0, ib, 120.0, null, 0.0, 0);
                      break;
                    case 0:
                      intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                      break;
                    default:
                      logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                      patched = false;
                  }
                  break;
                case 2:
                  switch (numBonds) {
                    case 1:
                      intxyz(hydrogen, ia, 1.0, ib, 120.0, null, 0.0, 0);
                      break;
                    case 0:
                      intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                      break;
                    default:
                      logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                      patched = false;
                  }
                  break;
                case 1:
                  if (numBonds == 0) {
                    intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                  } else {
                    logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                    patched = false;
                  }
                  break;
                default:
                  logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                  patched = false;
              }
              if (!patched) {
                break;
              } else {
                buildBond(ia, hydrogen, forceField, bondList);
              }
            }
          }
          if (!patched) {
            logger.log(Level.SEVERE, format(" Could not patch %s.", residueName));
          } else {
            logger.log(Level.INFO, " Patch for {0} succeeded.", residueName);
            residueAtoms = residue.getAtomList();
            // Assign atom types for side-chain atoms.
            double charge = 0.0;
            for (Atom atom : residueAtoms) {
              logger.info(atom.toString() + " -> " + atom.getAtomType().toString());
            }
          }
        } else {
          switch (position) {
            case FIRST_RESIDUE:
              buildH(
                  residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 355, forceField, bondList);
              break;
            case LAST_RESIDUE:
              buildH(
                  residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 506, forceField, bondList);
              break;
            default:
              buildH(
                  residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 6, forceField, bondList);
          }
        }
        break;
    }
  }

  /**
   * Check for missing heavy atoms. This check ignores special terminating groups like FOR, NH2,
   * etc.
   *
   * @param aminoAcid a {@link AminoAcid3} object.
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
   */
  private static void checkForMissingHeavyAtoms(AminoAcid3 aminoAcid, Residue residue)
      throws MissingHeavyAtomException {
    int expected = aminoAcid.heavyAtoms;
    if (aminoAcid != AminoAcid3.GLY && expected >= 4) {
      int actual = 0;
      List<Atom> resAtoms = residue.getAtomList();
      for (Atom atom : resAtoms) {
        String label = atom.getName().toUpperCase();
        if (!(label.equalsIgnoreCase("OXT") || label.equalsIgnoreCase("OT2"))) {
          if (!label.startsWith("H") && !label.startsWith("D")) {
            actual++;
          }
        }
      }
      if (actual != expected) {
        Atom N = (Atom) residue.getAtomNode("N");
        if (N == null) {
          MissingHeavyAtomException e = new MissingHeavyAtomException("N", null, null);
          logger.warning(
              format(
                  " Residue %c-%s is missing its N-terminal amide nitrogen",
                  residue.getChainID(), residue));
          throw e;
        }
        Atom CA = (Atom) residue.getAtomNode("CA");
        if (CA == null) {
          MissingHeavyAtomException e = new MissingHeavyAtomException("CA", null, null);
          logger.warning(
              format(" Residue %c-%s is missing its alpha carbon", residue.getChainID(), residue));
          throw e;
        }
        Atom C = (Atom) residue.getAtomNode("C");
        if (C == null) {
          MissingHeavyAtomException e = new MissingHeavyAtomException("C", null, null);
          logger.warning(
              format(
                  " Residue %c-%s is missing its C-terminal carboxyl carbon",
                  residue.getChainID(), residue));
          throw e;
        }
      }
    }
  }

  /**
   * Copy coordinates of an atom from one residue to another.
   *
   * @param fromResidue Use the coordinates from this residue.
   * @param toResidue Send the coordinates to this residue.
   * @param atomName The name of the atom whose coordinates will be updated.
   */
  private static void copyCoordinates(Residue fromResidue, Residue toResidue, String atomName) {
    Atom fromAtom;
    if (fromResidue.getAtomNode(atomName) != null) {
      fromAtom = (Atom) fromResidue.getAtomNode(atomName);
    } else {
      fromAtom = (Atom) fromResidue.getAtomNode("H1");
    }

    Atom toAtom = (Atom) toResidue.getAtomNode(atomName);
    toAtom.setXYZ(fromAtom.getXYZ(null));
  }

  /**
   * getAminoAcid.
   *
   * @param residueName a {@link String} object.
   * @return a {@link AminoAcid3} object.
   */
  public static AminoAcid3 getAminoAcid(String residueName) {
    for (AminoAcid3 aminoAcid : aminoAcidList) {
      if (aminoAcid.toString().equalsIgnoreCase(residueName)) {
        return aminoAcid;
      }
    }
    return AminoAcid3.UNK;
  }

  /**
   * This method takes a one letter amino acid code and converts it to a three letter amino acid
   * code. This method relies on the AminoAcid1 and AminoAcid3 enums having amino acids in exactly
   * the same order.
   *
   * @param residueName The one letter amino acid code.
   * @return The three letter amino acid code.
   */
  public static AminoAcid3 getAminoAcid3From1(String residueName) {
    for (AminoAcid1 aminoAcid : aminoAcid1List) {
      if (aminoAcid.toString().equalsIgnoreCase(residueName)) {
        int position = AminoAcid1.valueOf(residueName).ordinal();
        AminoAcid3 aminoAcid3 = AminoAcid3.values()[position];
        return aminoAcid3;
      }
    }
    return AminoAcid3.UNK;
  }

  /**
   * getAminoAcidNumber.
   *
   * @param residueName a {@link String} object.
   * @return a int.
   */
  public static int getAminoAcidNumber(String residueName) {
    int aminoAcidNumber = -1;
    for (AminoAcid3 aminoAcid : aminoAcidList) {
      aminoAcidNumber++;
      if (aminoAcid.toString().equalsIgnoreCase(residueName)) {
        break;
      }
    }
    return aminoAcidNumber;
  }

  /**
   * The 20 standard amino acids.
   */
  public enum AA {
    GLYCINE,
    ALANINE,
    VALINE,
    LEUCINE,
    ISOLEUCINE,
    SERINE,
    THREONINE,
    CYSTEINE,
    PROLINE,
    PHENYLALANINE,
    TYROSINE,
    TRYPTOPHAN,
    ASPARTATE,
    ASPARAGINE,
    GLUTAMATE,
    GLUTAMINE,
    METHIONINE,
    LYSINE,
    ARGININE,
    HISTIDINE
  }

  /** Single letter amino acid codes (need to */
  public enum AminoAcid1 {
    G,
    A,
    V,
    L,
    I,
    S,
    T,
    C,
    X,
    c,
    P,
    F,
    Y,
    y,
    W,
    H,
    U,
    Z,
    D, // ASP
    d, // ASH
    p, // ASD (double protonation)
    N,
    E, // GLU
    e, // GLH
    q, // GLD (double protonation)
    Q,
    M,
    K,
    k,
    R,
    O,
    B,
    J,
    t,
    f,
    a,
    o,
    n,
    m,
    x
  }

  public enum AminoAcid3 {
    // TODO: Check GLY, SER, THR, CYS and CYD multipole assignment during MultiResidue use.
    GLY(4),
    ALA(5, true),
    VAL(7, true),
    LEU(8, true),
    ILE(8, true),
    SER(6),
    THR(7),
    CYS(6, false, false, false),
    CYX(6),
    CYD(6, false, false, true),
    PRO(7),
    PHE(11, true),
    TYR(12, true, false, false),
    TYD(12, true, false, true),
    TRP(14, true),
    HIS(10, true, true, false),
    HID(10, true, true, true),
    HIE(10, true, true, true),
    ASP(8, true, true, false),
    ASH(8, true, true, true),
    ASD(8, true, true, true),
    ASN(8, true),
    GLU(9, true, true, false),
    GLH(9, true, true, true),
    GLD(9, true, true, true),
    GLN(9, true),
    MET(8, true),
    LYS(9, true, true, false),
    LYD(9, true, true, true),
    ARG(11, true),
    ORN(8),
    AIB(6),
    PCA(8),
    H2N(0),
    FOR(0),
    ACE(0),
    COH(0),
    NH2(0),
    NME(0),
    UNK(0);

    public final int heavyAtoms;
    public final boolean useWithMultiResidue;
    public final boolean isTitratable;
    public final boolean nonStandardProtonation;

    AminoAcid3(int heavyAtoms) {
      this.heavyAtoms = heavyAtoms;
      useWithMultiResidue = false;
      isTitratable = false;
      nonStandardProtonation = false;
    }

    AminoAcid3(int heavyAtoms, boolean useWithMultiResidue) {
      this.heavyAtoms = heavyAtoms;
      this.useWithMultiResidue = useWithMultiResidue;
      isTitratable = false;
      nonStandardProtonation = false;
    }

    AminoAcid3(int heavyAtoms, boolean useWithMultiResidue, boolean isTitratable,
        boolean nonStandardProtonation) {
      this.heavyAtoms = heavyAtoms;
      this.useWithMultiResidue = useWithMultiResidue;
      this.isTitratable = isTitratable;
      this.nonStandardProtonation = nonStandardProtonation;
      if (nonStandardProtonation && !this.isTitratable) {
        throw new IllegalArgumentException(
            format(
                " Amino acid class %s cannot be both nonstandard and non-titratable!", this));
      }
    }
  }

  /** The location of a residue within a chain. */
  public enum ResiduePosition {
    FIRST_RESIDUE,
    MIDDLE_RESIDUE,
    LAST_RESIDUE
  }

  /**
   * Biotype keys for amino acid backbone atom types. These are consistent with parameter files from
   * TINKER v. 6.1 (June 2012). <br> xType[0][..] are for N-terminal residues. <br> xType[1][..] are
   * mid-chain residues. <br> xType[2][..] are for C-terminal residues.
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_N = {
      {
          403, 409, 415, 421, 427, 433, 439, 445,
          451, 457, 463, 471, 477, 483, 489, 495,
          501, 507, 513, 519, 519, 525, 531, 537,
          537, 543, 549, 555, 561, 567, 573, 579,
          391, 762, 0, 0, 0, 0, 0, 403
      },
      {
          1, 7, 15, 27, 41, 55, 65, 77,
          87, 96, 105, 116, 131, 147, 162, 185,
          202, 218, 234, 244, 244, 256, 268, 280,
          280, 294, 308, 321, 337, 353, 370, 384,
          391, 0, 0, 0, 0, 0, 0, 1
      },
      {
          584, 590, 596, 602, 608, 614, 620, 626,
          632, 638, 644, 649, 655, 661, 667, 673,
          679, 685, 691, 697, 697, 703, 709, 715,
          715, 721, 727, 733, 739, 745, 751, 757,
          0, 0, 0, 0, 773, 775, 777, 584
      }
  };

  /**
   * Constant <code>AA_CA</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_CA = {
      {
          404, 410, 416, 422, 428, 434, 440, 446,
          452, 458, 464, 472, 478, 484, 490, 496,
          502, 508, 514, 520, 520, 526, 532, 538,
          538, 544, 550, 556, 562, 568, 574, 580,
          392, 0, 0, 767, 0, 0, 0, 404
      },
      {
          2, 8, 16, 28, 42, 56, 66, 78,
          88, 97, 106, 117, 132, 148, 163, 186,
          203, 219, 235, 245, 245, 257, 269, 281,
          281, 295, 309, 322, 338, 354, 371, 385,
          392, 0, 0, 0, 0, 0, 0, 2
      },
      {
          585, 591, 597, 603, 609, 615, 621, 627,
          633, 639, 645, 650, 656, 662, 668, 674,
          680, 686, 692, 698, 698, 704, 710, 716,
          716, 722, 728, 734, 740, 746, 752, 758,
          0, 0, 0, 0, 0, 0, 779, 585
      }
  };

  /**
   * Constant <code>AA_C</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_C = {
      {
          405, 411, 417, 423, 429, 435, 441, 447,
          453, 459, 465, 473, 479, 485, 491, 497,
          503, 509, 515, 521, 521, 527, 533, 539,
          539, 545, 551, 557, 563, 569, 575, 581,
          393, 0, 764, 769, 0, 0, 0, 405
      },
      {
          3, 9, 17, 29, 43, 57, 67, 79,
          89, 98, 107, 118, 133, 149, 164, 187,
          204, 220, 236, 246, 246, 258, 270, 282,
          282, 296, 310, 323, 339, 355, 372, 386,
          393, 0, 0, 0, 0, 0, 0, 3
      },
      {
          586, 592, 598, 604, 610, 616, 622, 628,
          634, 640, 646, 651, 657, 663, 669, 675,
          681, 687, 693, 699, 699, 705, 711, 717,
          717, 723, 729, 735, 741, 747, 753, 759,
          0, 0, 0, 0, 771, 0, 0, 586
      }
  };
  /**
   * Constant <code>AA_HN</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_HN = {
      {
          406, 412, 418, 424, 430, 436, 442, 448,
          454, 460, 466, 474, 480, 486, 492, 498,
          504, 510, 516, 522, 522, 528, 534, 540,
          540, 546, 552, 558, 564, 570, 576, 582,
          394, 763, 0, 0, 0, 0, 0, 406
      },
      {
          4, 10, 18, 30, 44, 58, 68, 80,
          90, 99, 0, 119, 134, 150, 165, 188,
          205, 221, 237, 247, 247, 259, 271, 283,
          283, 297, 311, 324, 340, 356, 373, 387,
          394, 0, 0, 0, 0, 0, 0, 4
      },
      {
          587, 593, 599, 605, 611, 617, 623, 629,
          635, 641, 0, 652, 658, 664, 670, 676,
          682, 688, 694, 700, 700, 706, 712, 718,
          718, 724, 730, 736, 742, 748, 754, 760,
          0, 0, 0, 0, 774, 776, 778, 587
      }
  };

  /**
   * Constant <code>AA_O</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_O = {
      {
          407, 413, 419, 425, 431, 437, 443, 449,
          455, 461, 467, 475, 481, 487, 493, 499,
          505, 511, 517, 523, 523, 529, 535, 541,
          541, 547, 553, 559, 565, 571, 577, 583,
          395, 0, 766, 770, 0, 0, 0, 407
      },
      {
          5, 11, 19, 31, 45, 59, 69, 81,
          91, 100, 108, 120, 135, 151, 166, 189,
          206, 222, 238, 248, 248, 260, 272, 284,
          284, 298, 312, 325, 341, 357, 374, 388,
          395, 0, 0, 0, 0, 0, 0, 5
      },
      {
          588, 594, 600, 606, 612, 618, 624, 630,
          636, 642, 647, 653, 659, 665, 671, 677,
          683, 689, 695, 701, 701, 707, 713, 719,
          719, 725, 731, 737, 743, 749, 755, 761,
          0, 0, 0, 0, 772, 0, 0, 588
      }
  };
  /**
   * Constant <code>AA_HA</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS CYX CYD PRO PHE TYR TYD TRP HIS HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[][] AA_HA = {
      {
          408, 414, 420, 426, 432, 438, 444, 450,
          456, 462, 468, 476, 482, 488, 494, 500,
          506, 512, 518, 524, 524, 530, 536, 542,
          542, 548, 554, 560, 566, 572, 578, 0,
          396, 0, 765, 768, 0, 0, 0, 408
      },
      {
          6, 12, 20, 32, 46, 60, 70, 82,
          92, 101, 109, 121, 136, 152, 167, 190,
          207, 223, 239, 249, 249, 261, 273, 285,
          285, 299, 313, 326, 342, 358, 375, 0,
          396, 0, 0, 0, 0, 0, 0, 6
      },
      {
          589, 595, 601, 607, 613, 619, 625, 631,
          637, 643, 648, 654, 660, 666, 672, 678,
          684, 690, 696, 702, 702, 708, 714, 720,
          720, 726, 732, 738, 744, 750, 756, 0,
          0, 0, 0, 0, 0, 0, 780, 589
      }
  };
  /**
   * Constant <code>AA_CB</code>
   * <p>
   * GLY ALA VAL LEU ILE SER THR CYS
   * CYX CYD PRO PHE TYR TYD TRP HIS
   * HID HIE ASP ASH ASD ASN GLU GLH
   * GLD GLN MET LYS LYD ARG ORN AIB
   * PCA H2N FOR ACE COH NH2 NME UNK
   */
  public static final int[] AA_CB = {
      0, 13, 21, 33, 47, 61, 71, 83,
      93, 102, 110, 122, 137, 153, 168, 191,
      208, 224, 240, 250, 806, 262, 274, 286,
      817, 300, 314, 327, 343, 359, 376, 389,
      397, 0, 0, 0, 0, 0, 0, 0
  };
}
