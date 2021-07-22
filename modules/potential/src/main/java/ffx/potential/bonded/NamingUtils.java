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

import static ffx.potential.bonded.BondedUtils.findAtomsOfElement;
import static ffx.potential.bonded.BondedUtils.findBondedAtoms;
import static ffx.potential.bonded.BondedUtils.findNitrogenAtom;
import static ffx.potential.bonded.BondedUtils.findNucleotideO4s;
import static ffx.potential.bonded.BondedUtils.getAlphaCarbon;
import static ffx.potential.bonded.BondedUtils.hasAttachedAtom;
import static ffx.potential.bonded.BondedUtils.sortAtomsByDistance;
import static ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import static ffx.potential.bonded.NucleicAcidUtils.NucleicAcid3;
import static ffx.potential.bonded.AminoAcidUtils.getAminoAcid;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import ffx.numerics.math.DoubleMath;
import ffx.numerics.math.ScalarMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.parsers.PDBFilter.PDBFileStandard;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Utilities for importing atoms from PDB files and checking their names.
 *
 * @author Jacob Litman
 * @author Michael Schnieders
 * @since 1.0
 */
public class NamingUtils {

  private static final Logger logger = Logger.getLogger(NamingUtils.class.getName());

  /**
   * Ensures proper naming of hydrogens according to latest PDB format. Presently mostly guesses at
   * which hydrogens to re-assign, which may cause chirality errors for prochiral hydrogens. If
   * necessary, we will implement more specific mapping.
   *
   * @param residue Residue to examine.
   * @param fileStandard PDB File Standard to use.
   */
  public static void checkHydrogenAtomNames(Residue residue, PDBFileStandard fileStandard) {
    switch (fileStandard) {
      case VERSION3_3:
        return;
      case VERSION3_2:
      default:
        break;
    }
    // May have to get position.
    String residueType = residue.getName().toUpperCase();
    List<Atom> resAtoms = residue.getAtomList();
    for (Atom atom : resAtoms) {
      if (atom == null) {
        continue;
      }
      String atomName = atom.getName().toUpperCase();
      // Handles situations such as 1H where it should be NA_H1, etc.
      if (atomName.contains("H")) {
        try {
          String firstChar = atomName.substring(0, 1);
          parseInt(firstChar);
          atomName = atomName.substring(1);
          atomName = atomName.concat(firstChar);
          atom.setName(atomName);
        } catch (NumberFormatException e) {
          // Do nothing.
        }
      }
    }
    // Ensures proper hydrogen assignment; for example, Gln should have HB2,
    // HB3 instead of HB1, HB2.
    List<Atom> betas;
    List<Atom> gammas;
    List<Atom> deltas;
    List<Atom> epsilons;
    List<Atom> zetas;
    String atomName;
    Atom OH;
    Atom HH;
    Atom HG;
    Atom HD2;
    switch (getAminoAcid(residueType)) {
      case GLY:
        List<Atom> alphas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          if (atom.getName().toUpperCase().contains("HA")) {
            alphas.add(atom);
          }
        }
        renameGlycineAlphaHydrogens(residue, alphas);
        break;
      case ALA:
        // No known errors with alanine
        break;
      case VAL:
        // No known errors with valine
        break;
      case LEU:
      case SER:
      case CYD:
      case ASP:
        betas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          if (atom.getName().toUpperCase().contains("HB")) {
            betas.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        break;
      case ILE:
        List<Atom> ileAtoms = new ArrayList<>();
        for (Atom atom : resAtoms) {
          if (atom.getName().toUpperCase().contains("HG1")) {
            ileAtoms.add(atom);
          }
        }
        renameIsoleucineHydrogens(residue, ileAtoms);
        break;
      case THR:
        Atom HG1 = (Atom) residue.getAtomNode("HG1");
        if (HG1 == null) {
          for (Atom atom : resAtoms) {
            atomName = atom.getName().toUpperCase();
            // Gets first HG-containing name of length < 4
            // Length < 4 avoids bringing in HG21, HG22, or HG23.
            if (atomName.length() < 4 && atomName.contains("HG")) {
              atom.setName("HG1");
              break;
            }
          }
        }
        break;
      case CYS:
        betas = new ArrayList<>();
        HG = (Atom) residue.getAtomNode("HG");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (HG == null && atomName.contains("HG")) {
            HG = atom;
            HG.setName("HG");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        break;
      case CYX:
        // I pray this is never important, because I don't have an example CYX to work from.
        break;
      case PRO:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        deltas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameDeltaHydrogens(residue, deltas, 23);
        break;
      case PHE:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        Atom HZ = (Atom) residue.getAtomNode("HZ");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (HZ == null && atomName.contains("HZ")) {
            HZ = atom;
            HZ.setName("HZ");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case TYR:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        HH = (Atom) residue.getAtomNode("HH");
        OH = (Atom) residue.getAtomNode("OH");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (HH == null && atomName.contains("HH")) {
            HH = atom;
            HH.setName("HH");
          } else if (OH == null && atomName.contains("O") && atomName.contains("H")) {
            OH = atom;
            OH.setName("OH");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case TYD:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        OH = (Atom) residue.getAtomNode("OH");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (OH == null && atomName.contains("O") && atomName.contains("H")) {
            OH = atom;
            OH.setName("OH");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case TRP:
        betas = new ArrayList<>();
        epsilons = new ArrayList<>();
        zetas = new ArrayList<>();
        Atom HD1 = (Atom) residue.getAtomNode("HD1");
        Atom HH2 = (Atom) residue.getAtomNode("HH2");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (atomName.contains("HZ")) {
            zetas.add(atom);
          } else if (HD1 == null && atomName.contains("HD")) {
            HD1 = atom;
            HD1.setName("HD1");
          } else if (HH2 == null && atomName.contains("HH")) {
            HH2 = atom;
            HH2.setName("HH2");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameEpsilonHydrogens(residue, epsilons, 13);
        renameZetaHydrogens(residue, zetas, 23);
        break;
      case HIS:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case HID:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        Atom HE1 = (Atom) residue.getAtomNode("HE1");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (HE1 == null && atomName.contains("HE")) {
            HE1 = atom;
            HE1.setName("HE1");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        break;
      case HIE:
        betas = new ArrayList<>();
        epsilons = new ArrayList<>();
        HD2 = (Atom) residue.getAtomNode("HD2");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (HD2 == null && atomName.contains("HD")) {
            HD2 = atom;
            HD2.setName("HD2");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case ASH:
        betas = new ArrayList<>();
        HD2 = (Atom) residue.getAtomNode("HD2");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (HD2 == null && atomName.contains("HD")) {
            HD2 = atom;
            HD2.setName("HD2");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        break;
      case ASD:
        betas = new ArrayList<>();
        deltas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameDeltaHydrogens(residue, deltas, 12);
        break;
      case ASN:
        betas = new ArrayList<>();
        List<Atom> HD2s = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HD")) {
            HD2s.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameAsparagineHydrogens(residue, HD2s);
        break;
      case GLU:
      case MET:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        break;
      case GLH:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        Atom HE2 = (Atom) residue.getAtomNode("HE2");
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (HE2 == null && atomName.contains("HE")) {
            HE2 = atom;
            HE2.setName("HE2");
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        break;
      case GLD:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        epsilons = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameEpsilonHydrogens(residue, epsilons, 12);
        break;
      case GLN:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        epsilons = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameGlutamineHydrogens(residue, epsilons);
        break;
      // Epsilons should not break, as they are 1-3.
      case LYS:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        // Zetas are 1-3, should not break.
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameDeltaHydrogens(residue, deltas, 23);
        renameEpsilonHydrogens(residue, epsilons, 23);
        break;
      case LYD:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        deltas = new ArrayList<>();
        epsilons = new ArrayList<>();
        zetas = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (atomName.contains("HE")) {
            epsilons.add(atom);
          } else if (atomName.contains("HZ")) {
            zetas.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameDeltaHydrogens(residue, deltas, 23);
        renameEpsilonHydrogens(residue, epsilons, 23);
        renameZetaHydrogens(residue, zetas, 12);
        break;
      case ARG:
        betas = new ArrayList<>();
        gammas = new ArrayList<>();
        deltas = new ArrayList<>();
        Atom HE = (Atom) residue.getAtomNode("HE");
        List<Atom> HHn = new ArrayList<>();
        for (Atom atom : resAtoms) {
          atomName = atom.getName().toUpperCase();
          if (atomName.contains("HB")) {
            betas.add(atom);
          } else if (atomName.contains("HG")) {
            gammas.add(atom);
          } else if (atomName.contains("HD")) {
            deltas.add(atom);
          } else if (HE == null && atomName.contains("HE")) {
            HE = atom;
            HE.setName("HE");
          } else if (atomName.contains("HH")) {
            HHn.add(atom);
          }
        }
        renameBetaHydrogens(residue, betas, 23);
        renameGammaHydrogens(residue, gammas, 23);
        renameDeltaHydrogens(residue, deltas, 23);
        renameArginineHydrogens(residue, HHn);
        break;
      case ORN:
      case AIB:
      case PCA:
      case UNK:
      default:
        // I am currently unaware of how these amino acids are typically
        // labeled under older PDB standards.
        break;
    }
  }

  /**
   * Names the atoms in an N-terminal acetyl ACE capping group.
   *
   * @param residue Residue containing an acetyl cap.
   * @param aceC The acetyl group's C atom.
   */
  public static void nameAcetylCap(Residue residue, Atom aceC) {
    logger.warning(
        format(
            " Probable ACE cap attached to residue %s; duplicate atom names may result.", residue));
    aceC.setName("C");
    findBondedAtoms(aceC, 8).get(0).setName("O");
    Atom CH3 = findBondedAtoms(aceC, 6).get(0);
    CH3.setName("CH3");
    List<Atom> ntermHs = findBondedAtoms(CH3, 1);
    for (int i = 0; i < 3; i++) {
      ntermHs.get(i).setName(format("H%d", (i + 1)));
    }
  }

  /**
   * Renames an atom, its bonded hydrogens, and returns the next atom in the chain.
   *
   * <p>If applied to an atom that is not a carbon, it will be misnamed as a carbon, so fix that
   * afterwards.
   *
   * @param carbon Alkyl carbon to rename.
   * @param priorAtom Prior atom in the chain.
   * @param protonOffset Number of the first hydrogen (such as 2 for HB2-3).
   * @param posName Name of the position (such as B for CB).
   * @return Next atom in the chain if present.
   */
  public static Optional<Atom> renameAlkyl(
      Atom carbon, Atom priorAtom, int protonOffset, char posName) {
    carbon.setName(format("C%c", posName));
    List<Atom> hydrogens = findBondedAtoms(carbon, 1);
    int numH = hydrogens.size();
    if (numH == 1) {
      hydrogens.get(0).setName(format("H%c", posName));
    } else {
      for (int i = 0; i < numH; i++) {
        hydrogens.get(i).setName(format("H%c%d", posName, i + protonOffset));
      }
    }

    return carbon.getBonds().stream()
        .map((Bond b) -> b.get1_2(carbon))
        .filter((Atom a) -> a != priorAtom)
        .filter((Atom a) -> !hydrogens.contains(a))
        .findAny();
  }

  /**
   * Renames the Atoms in an amino acid to PDB standard.
   *
   * @param residue Residue to fix atom names of.
   */
  public static void renameAminoAcidToPDBStandard(Residue residue) {
    if (residue.getChainID() == null) {
      residue.setChainID('Z');
    }
    final Atom N = findNitrogenAtom(residue);
    AminoAcid3 aa3 = residue.getAminoAcid3();
    if (N != null) {
      N.setName("N");

      Atom CA = getAlphaCarbon(residue, N);
      CA.setName("CA");

      List<Atom> has = findBondedAtoms(CA, 1);
      switch (aa3) {
        case NME:
          // Do all renaming here then return out of the method.
          findBondedAtoms(N, 1).get(0).setName("H");
          CA.setName("CH3");
          for (int i = 1; i <= 3; i++) {
            has.get(i - 1).setName(format("H%d", i));
          }
          return;
        case GLY:
          has.get(0).setName("HA2");
          has.get(1).setName("HA3");
          break;
        default:
          has.get(0).setName("HA");
          break;
      }

      Atom C = null;
      Atom CB = null;
      for (Atom carb : findBondedAtoms(CA, 6)) {
        // Second check is because of serine/threonine OG bonded straight to CB.
        if (hasAttachedAtom(carb, 8) && !hasAttachedAtom(carb, 1)) {
          C = carb;
          C.setName("C");
        } else {
          CB = carb;
          CB.setName("CB");
        }
      }
      if (C == null) {
        throw new IllegalArgumentException(
            format(" Could not find carbonyl carbon for residue %s!", residue));
      }
      if (CB == null && aa3 != AminoAcidUtils.AminoAcid3.GLY) {
        throw new IllegalArgumentException(
            format(" Could not find beta carbon for residue %s!", residue));
      }

      List<Atom> ctermOxygens = findBondedAtoms(C, 8);
      switch (ctermOxygens.size()) {
        case 1:
          ctermOxygens.get(0).setName("O");
          break;
        case 2:
          Atom O = null;
          for (Atom oxy : ctermOxygens) {
            if (oxy.getBonds().size() == 2) {
              O = oxy;
              O.setName("O");
              findBondedAtoms(O, 1).get(0).setName("HO");
            }
          }
          if (O == null) {
            ctermOxygens.get(0).setName("O");
            ctermOxygens.get(1).setName("OXT");
          }
      }

      List<Atom> amideProtons = findBondedAtoms(N, 1);
      if (amideProtons.size() == 1) {
        amideProtons.get(0).setName("H");
      } else {// Should catch both N-termini and proline.
        for (int i = 1; i <= amideProtons.size(); i++) {
          amideProtons.get(i - 1).setName(format("H%d", i));
        }
      }

      // All common atoms are now named: N, H[1-3], CA, HA[2-3], CB, C, O[XT], [HO]
      renameCommonAminoAcids(residue, aa3, CA, CB);
    } else {
      if (aa3 == AminoAcid3.ACE) {
        Atom O = findAtomsOfElement(residue, 8).get(0);
        O.setName("O");
        Atom C = findBondedAtoms(O, 6).get(0);
        C.setName("C");
        Atom CH3 = findBondedAtoms(C, 6).get(0);
        CH3.setName("CH3");
        List<Atom> hydrogens = findBondedAtoms(CH3, 1);
        for (int i = 1; i <= 3; i++) {
          hydrogens.get(i - 1).setName(format("H%d", i));
        }
      } else {
        throw new IllegalArgumentException(
            format(" Could not find nitrogen atom for residue %s!", residue));
      }
    }
  }

  /**
   * renameArginineHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   */
  public static void renameArginineHydrogens(Residue residue, List<Atom> resAtoms) {
    Atom HH11 = (Atom) residue.getAtomNode("HH11");
    Atom HH12 = (Atom) residue.getAtomNode("HH12");
    Atom HH21 = (Atom) residue.getAtomNode("HH21");
    Atom HH22 = (Atom) residue.getAtomNode("HH22");
    if (HH11 != null) {
      resAtoms.remove(HH11);
    }
    if (HH12 != null) {
      resAtoms.remove(HH12);
    }
    if (HH21 != null) {
      resAtoms.remove(HH21);
    }
    if (HH22 != null) {
      resAtoms.remove(HH22);
    }
    if (!resAtoms.isEmpty() && HH11 == null) {
      resAtoms.get(0).setName("HH11");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HH12 == null) {
      resAtoms.get(0).setName("HH12");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HH21 == null) {
      resAtoms.get(0).setName("HH21");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HH22 == null) {
      resAtoms.get(0).setName("HH22");
      resAtoms.remove(0);
    }
  }

  /**
   * renameAsparagineHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   */
  public static void renameAsparagineHydrogens(Residue residue, List<Atom> resAtoms) {
    Atom HD21 = (Atom) residue.getAtomNode("HD21");
    Atom HD22 = (Atom) residue.getAtomNode("HD22");
    if (HD21 != null) {
      resAtoms.remove(HD21);
    }
    if (HD22 != null) {
      resAtoms.remove(HD22);
    }
    if (!resAtoms.isEmpty() && HD21 == null) {
      resAtoms.get(0).setName("HD21");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HD22 == null) {
      resAtoms.get(0).setName("HD21");
    }
  }

  /**
   * Renames Atoms to PDB standard using bonding patterns, atomic numbers, and residue types.
   *
   * <p>Will not work if a force field definition botches its atomic numbers.
   *
   * <p>Only implemented for amino acids and nucleic acids at this time.
   *
   * @param molecularAssembly MolecularAssembly to fix.
   */
  public static void renameAtomsToPDBStandard(MolecularAssembly molecularAssembly) {
    Polymer[] polymers = molecularAssembly.getChains();
    if (polymers != null && polymers.length > 0) {
      for (Polymer polymer : polymers) {
        for (Residue residue : polymer.getResidues()) {
          switch (residue.getResidueType()) {
            case AA:
              renameAminoAcidToPDBStandard(residue);
              break;
            case NA:
              renameNucleicAcidToPDBStandard(residue);
              break;
            case UNK:
            default:
              break;
          }
        }
      }
    }
  }

  /**
   * renameBetaHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   * @param indexes a int.
   */
  public static void renameBetaHydrogens(Residue residue, List<Atom> resAtoms, int indexes) {
    Atom[] HBn = new Atom[3];
    switch (indexes) {
      case 12:
        HBn[0] = (Atom) residue.getAtomNode("HB1");
        HBn[1] = (Atom) residue.getAtomNode("HB2");
        break;
      case 13:
        HBn[0] = (Atom) residue.getAtomNode("HB1");
        HBn[2] = (Atom) residue.getAtomNode("HB3");
        break;
      case 23:
        HBn[1] = (Atom) residue.getAtomNode("HB2");
        HBn[2] = (Atom) residue.getAtomNode("HB3");
        break;
      default:
        return;
    }
    for (Atom HBatom : HBn) {
      resAtoms.remove(HBatom);
    }
    if (!resAtoms.isEmpty() && HBn[0] == null && (indexes == 12 || indexes == 13)) {
      resAtoms.get(0).setName("HB1");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HBn[1] == null && (indexes == 12 || indexes == 23)) {
      resAtoms.get(0).setName("HB2");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HBn[2] == null && (indexes == 13 || indexes == 23)) {
      resAtoms.get(0).setName("HB3");
      resAtoms.remove(0);
    }
  }

  /**
   * Renames a numbered carbon, its bonded hydrogens, and returns the next atom in the chain.
   *
   * <p>If applied to an atom that is not a carbon, it will be misnamed as a carbon, so fix that
   * afterwards.
   *
   * <p>This is for carbons like PHE CD1 and CD2.
   *
   * @param carbon Alkyl carbon to rename.
   * @param priorAtom Prior atom in the chain.
   * @param protonOffset Number of the first hydrogen (such as 2 for HB2-3).
   * @param branchNum Index of the branch.
   * @param posName Name of the position (such as B for CB).
   * @return Next atom in the chain if present.
   */
  public static Optional<Atom> renameBranchedAlkyl(
      Atom carbon, Atom priorAtom, int protonOffset, int branchNum, char posName) {
    carbon.setName(format("C%c%d", posName, branchNum));
    List<Atom> hydrogens = findBondedAtoms(carbon, 1);
    int numH = hydrogens.size();
    if (numH == 1) {
      hydrogens.get(0).setName(format("H%c%d", posName, branchNum));
    } else {
      for (int i = 0; i < numH; i++) {
        hydrogens.get(i).setName(format("H%c%d%d", posName, branchNum, i + protonOffset));
      }
    }

    return carbon.getBonds().stream()
        .map((Bond b) -> b.get1_2(carbon))
        .filter((Atom a) -> a != priorAtom)
        .filter((Atom a) -> !hydrogens.contains(a))
        .findAny();
  }

  /**
   * Renames atoms in common amino acids to PDB standard.
   *
   * @param residue Residue to perform renaming for.
   * @param aa3 Its AA3 code.
   * @param CA Its alpha carbon.
   * @param CB Its beta carbon.
   */
  public static void renameCommonAminoAcids(Residue residue, AminoAcid3 aa3, Atom CA, Atom CB) {
    switch (aa3) {
      case ALA:
        {
          renameAlkyl(CB, CA, 1, 'B');
        }
        break;
      case CYS:
      case CYD:
        {
          Atom SG = renameAlkyl(CB, CA, 2, 'B').get();
          SG.setName("SG");
          if (hasAttachedAtom(SG, 1)) {
            assert aa3 == AminoAcidUtils.AminoAcid3.CYS;
            findBondedAtoms(SG, 1).get(0).setName("HG");
          } else if (hasAttachedAtom(SG, 16)) {
            logger.finer(format(" SG atom %s likely part of a disulfide bond.", SG));
          } else {
            residue.setName("CYD");
          }
        }
        break;
      case ASP:
      case ASH:
      case ASD:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          List<Atom> ODs = findBondedAtoms(CG, 8);

          int protonatedOD = -1; // -1: Deprotonated ASP. 0/1: Index of protonated oxygen (ASH).
          for (int i = 0; i < 2; i++) {
            if (hasAttachedAtom(ODs.get(i), 1)) {
              protonatedOD = i;
              break;
            }
          }

          // Check for double protonation for constant pH.
          if (hasAttachedAtom(ODs.get(0),1) && hasAttachedAtom(ODs.get(1), 1)) {
            protonatedOD = 2;
          }

          switch (protonatedOD) {
            case -1:
              ODs.get(0).setName("OD1");
              ODs.get(1).setName("OD2");
              break;
            case 0:
              if (aa3 != AminoAcidUtils.AminoAcid3.ASH) {
                residue.setName("ASH");
              }
              ODs.get(0).setName("OD2");
              findBondedAtoms(ODs.get(0), 1).get(0).setName("HD2");
              ODs.get(1).setName("OD1");
              break;
            case 1:
              if (aa3 != AminoAcidUtils.AminoAcid3.ASH) {
                residue.setName("ASH");
              }
              ODs.get(1).setName("OD2");
              findBondedAtoms(ODs.get(1), 1).get(0).setName("HD2");
              ODs.get(0).setName("OD1");
              break;
            case 2:
              if (aa3 != AminoAcidUtils.AminoAcid3.ASD) {
                residue.setName("ASD");
              }
              ODs.get(0).setName("OD1");
              findBondedAtoms(ODs.get(0), 1).get(0).setName("HD1");
              ODs.get(1).setName("OD2");
              findBondedAtoms(ODs.get(1), 1).get(0).setName("HD2");
              break;
          }
        }
        break;
      case GLU:
      case GLH:
      case GLD:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
          CD.setName("CD");
          List<Atom> OEs = findBondedAtoms(CD, 8);

          int protonatedOE = -1; // If it remains -1, deprotonated ASP, else ASH.
          for (int i = 0; i < 2; i++) {
            if (hasAttachedAtom(OEs.get(i), 1)) {
              protonatedOE = i;
              break;
            }
          }

          // Check for double protonation for constant pH.
          if (hasAttachedAtom(OEs.get(0),1) && hasAttachedAtom(OEs.get(1), 1)) {
            protonatedOE = 2;
          }

          switch (protonatedOE) {
            case -1:
              OEs.get(0).setName("OE1");
              OEs.get(1).setName("OE2");
              break;
            case 0:
              if (aa3 != AminoAcidUtils.AminoAcid3.GLH) {
                residue.setName("GLH");
              }
              OEs.get(0).setName("OE2");
              findBondedAtoms(OEs.get(0), 1).get(0).setName("HE2");
              OEs.get(1).setName("OE1");
              break;
            case 1:
              if (aa3 != AminoAcidUtils.AminoAcid3.GLH) {
                residue.setName("GLH");
              }
              OEs.get(1).setName("OE2");
              findBondedAtoms(OEs.get(1), 1).get(0).setName("HE2");
              OEs.get(0).setName("OE1");
              break;
            case 2:
              if (aa3 != AminoAcidUtils.AminoAcid3.GLD) {
                residue.setName("GLD");
              }
              OEs.get(0).setName("OE1");
              findBondedAtoms(OEs.get(0), 1).get(0).setName("HE1");
              OEs.get(1).setName("OE2");
              findBondedAtoms(OEs.get(1), 1).get(0).setName("HE2");
              break;
          }
        }
        break;
      case PHE:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          List<Atom> CDs = findBondedAtoms(CG, CB, 6);

          Atom CZ = null;
          for (int i = 1; i <= 2; i++) {
            Atom CD = CDs.get(i - 1);
            Atom CE = renameBranchedAlkyl(CD, CG, 0, i, 'D').get();
            CZ = renameBranchedAlkyl(CE, CD, 0, i, 'E').get();
          }
          CZ.setName("CZ");
          findBondedAtoms(CZ, 1).get(0).setName("HZ");
        }
        break;
      case GLY:
        break;
      case HIS:
      case HIE:
      case HID:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");

          Atom CD2 = findBondedAtoms(CG, 6).stream().filter((Atom a) -> a != CB).findAny().get();
          CD2.setName("CD2");
          findBondedAtoms(CD2, 1).get(0).setName("HD2");

          Atom NE2 = findBondedAtoms(CD2, 7).get(0);
          NE2.setName("NE2");
          List<Atom> HE2 = findBondedAtoms(NE2, 1);
          boolean epsProtonated = (HE2 != null && !HE2.isEmpty());
          if (epsProtonated) {
            HE2.get(0).setName("HE2");
          }

          Atom CE1 = findBondedAtoms(NE2, CD2, 6).get(0);
          CE1.setName("CE1");
          findBondedAtoms(CE1, 1).get(0).setName("HE1");

          Atom ND1 = findBondedAtoms(CG, 7).get(0);
          ND1.setName("ND1");
          List<Atom> HD1 = findBondedAtoms(ND1, 1);
          boolean deltaProtonated = (HD1 != null && !HD1.isEmpty());
          if (deltaProtonated) {
            HD1.get(0).setName("HD1");
          }

          // All constant atoms found: now check protonation state.
          if (epsProtonated && deltaProtonated) {
            assert aa3 == AminoAcidUtils.AminoAcid3.HIS;
          } else if (epsProtonated) {
            residue.setName("HIE");
          } else if (deltaProtonated) {
            residue.setName("HID");
          } else {
            throw new IllegalArgumentException(
                format(" Histidine residue %s is doubly deprotonated!", residue));
          }
        }
        break;
      case ILE:
        {
          findBondedAtoms(CB, 1).get(0).setName("HB");
          List<Atom> CGs = findBondedAtoms(CB, CA, 6);

          for (Atom CG : CGs) {
            List<Atom> HGs = findBondedAtoms(CG, 1);
            int numHGs = HGs.size();
            if (numHGs == 3) {
              renameBranchedAlkyl(CG, CB, 1, 2, 'G');
            } else if (numHGs == 2) {
              Atom CD1 = renameBranchedAlkyl(CG, CB, 2, 1, 'G').get();
              renameBranchedAlkyl(CD1, CG, 1, 1, 'D');
            } else {
              throw new IllegalArgumentException(
                  format(
                      " Isoleucine residue %s had %d gamma hydrogens, expecting 2-3!",
                      residue, numHGs));
            }
          }
        }
        break;
      case LYS:
      case LYD:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
          Atom CE = renameAlkyl(CD, CG, 2, 'D').get();
          Atom NZ = renameAlkyl(CE, CD, 2, 'E').get();
          // For a very brief period, NZ will be named CZ.
          renameAlkyl(NZ, CE, 1, 'Z');
          NZ.setName("NZ");
          int numH = findBondedAtoms(NZ, 1).size();
          switch (numH) {
            case 2:
              residue.setName("LYD");
              break;
            case 3:
              assert aa3 == AminoAcidUtils.AminoAcid3.LYS;
              break;
            default:
              throw new IllegalArgumentException(
                  format(" Lysine residue %s had %d amine protons, expecting 2-3!", residue, numH));
          }
        }
        break;
      case LEU:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          findBondedAtoms(CG, 1).get(0).setName("HG");
          List<Atom> CDs = findBondedAtoms(CG, CB, 6);

          for (int i = 0; i < 2; i++) {
            renameBranchedAlkyl(CDs.get(i), CG, 1, (i + 1), 'D');
          }
        }
        break;
      case MET:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom SD = renameAlkyl(CG, CB, 2, 'G').get();
          Atom CE = renameAlkyl(SD, CG, 0, 'D').get();
          // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
          SD.setName("SD");
          renameAlkyl(CE, SD, 1, 'E');
        }
        break;
      case ASN:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          findBondedAtoms(CG, 8).get(0).setName("OD1");
          Atom ND2 = findBondedAtoms(CG, 7).get(0);
          renameBranchedAlkyl(ND2, CG, 1, 2, 'D');
          // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
          ND2.setName("ND2");
        }
        break;
      case PRO:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
          Atom N = renameAlkyl(CD, CG, 2, 'D').get();
          assert N.getName().equals("N");
        }
        break;
      case GLN:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
          CD.setName("CD");

          findBondedAtoms(CD, 8).get(0).setName("OE1");
          Atom NE2 = findBondedAtoms(CD, 7).get(0);
          renameBranchedAlkyl(NE2, CD, 1, 2, 'E');
          // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
          NE2.setName("NE2");
        }
        break;
      case ARG:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
          Atom NE = renameAlkyl(CD, CG, 2, 'D').get();
          Atom CZ = renameAlkyl(NE, CD, 0, 'E').get();
          NE.setName("NE");
          CZ.setName("CZ");

          List<Atom> NHs = findBondedAtoms(CZ, NE, 7);
          assert NHs.size() == 2;
          for (int i = 0; i < 2; i++) {
            Atom NHx = NHs.get(i);
            renameBranchedAlkyl(NHx, CZ, 1, (i + 1), 'H');
            NHx.setName(format("NH%d", (i + 1)));
          }
        }
        break;
      case SER:
        {
          Atom OG = renameAlkyl(CB, CA, 2, 'B').get();
          renameAlkyl(OG, CB, 0, 'G');
          OG.setName("OG");
        }
        break;
      case THR:
        {
          CB.setName("CB"); // Should be unnecessary.
          findBondedAtoms(CB, 1).get(0).setName("HB");

          Atom OG1 = findBondedAtoms(CB, 8).get(0);
          OG1.setName("OG1");
          findBondedAtoms(OG1, 1).get(0).setName("HG1");

          Atom CG2 = findBondedAtoms(CB, CA, 6).get(0);
          renameBranchedAlkyl(CG2, CB, 1, 2, 'G');
        }
        break;
      case VAL:
        {
          CB.setName("CB"); // Should be unnecessary.
          findBondedAtoms(CB, 1).get(0).setName("HB");

          List<Atom> CGs = findBondedAtoms(CB, CA, 6);

          assert CGs.size() == 2;
          for (int i = 0; i < 2; i++) {
            Atom CGx = CGs.get(i);
            renameBranchedAlkyl(CGx, CB, 1, (i + 1), 'G');
          }
        }
        break;
      case TRP:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          List<Atom> CDs = findBondedAtoms(CG, CB, 6);
          Atom CD1 = null;
          Atom CD2 = null;

          for (Atom CDx : CDs) {
            if (hasAttachedAtom(CDx, 1)) {
              CD1 = CDx;
            } else {
              CD2 = CDx;
              CD2.setName("CD2");
            }
          }
          Atom NE1 = renameBranchedAlkyl(CD1, CG, 0, 1, 'D').get();
          Atom CE2 = renameBranchedAlkyl(NE1, CD1, 0, 1, 'E').get();
          NE1.setName("NE1");
          CE2.setName("CE2");

          Atom CZ2 = findBondedAtoms(CE2, CD2, 6).get(0);
          Atom CH2 = renameBranchedAlkyl(CZ2, CE2, 0, 2, 'Z').get();
          Atom CZ3 = renameBranchedAlkyl(CH2, CZ2, 0, 2, 'H').get();
          Atom CE3 = renameBranchedAlkyl(CZ3, CH2, 0, 3, 'Z').get();
          if (CD2 != renameBranchedAlkyl(CE3, CZ3, 0, 3, 'E').get()) {
            throw new IllegalArgumentException(
                format(" Error in cyclizing tryptophan %s!", residue));
          }
        }
        break;
      case TYR:
      case TYD:
        {
          Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
          CG.setName("CG");
          List<Atom> CDs = findBondedAtoms(CG, CB, 6);
          Atom CZ = null;

          assert CDs.size() == 2;
          for (int i = 1; i <= 2; i++) {
            Atom CDx = CDs.get(i - 1);
            Atom CEx = renameBranchedAlkyl(CDx, CG, 0, i, 'D').get();
            CZ = renameBranchedAlkyl(CEx, CDx, 0, i, 'E').get();
          }

          CZ.setName("CZ");
          Atom OH = findBondedAtoms(CZ, 8).get(0);
          OH.setName("OH");
          if (hasAttachedAtom(OH, 1)) {
            assert aa3 == AminoAcidUtils.AminoAcid3.TYR;
            findBondedAtoms(OH, 1).get(0).setName("HH");
          } else {
            residue.setName("TYD");
          }
        }
        break;
      default:
        throw new IllegalArgumentException((format(" Amino acid %s (%s) not recognized!", residue, aa3)));
    }
  }

  /**
   * Renames atoms in common nucleic acids to PDB standard.
   *
   * @param residue Residue to perform renaming for.
   * @param na3 Its NA3 code.
   */
  public static void renameCommonNucleicAcid(Residue residue, NucleicAcid3 na3) {
    Optional<Atom> optO4s = findNucleotideO4s(residue);
    if (optO4s.isPresent()) {
      // Name O4', which is the unique ether oxygen.
      Atom O4s = optO4s.get();
      O4s.setName("O4'");

      // C1' is bonded to a nitrogen (at least for non-abasic sites), C4' isn't
      List<Atom> bondedC = findBondedAtoms(O4s, 6);
      Atom C4s = null;
      Atom C1s = null;
      // Will need the first base nitrogen (N1/N9), and H1' later.
      Atom N19 = null;
      Atom H1s = null;
      for (Atom c : bondedC) {
        if (hasAttachedAtom(c, 7)) {
          C1s = c;
          C1s.setName("C1'");
          H1s = findBondedAtoms(C1s, 1).get(0);
          H1s.setName("H1'");
          N19 = findBondedAtoms(C1s, 7).get(0);
        } else {
          C4s = c;
          C4s.setName("C4'");
          findBondedAtoms(C4s, 1).get(0).setName("H4'");
        }
      }
      assert C4s != null && C1s != null;

      Atom C2s = findBondedAtoms(C1s, 6).get(0);
      C2s.setName("C2'");

      bondedC = findBondedAtoms(C4s, 6);
      Atom C5s = null;
      Atom C3s;
      Atom O3s;
      for (Atom c : bondedC) {
        if (c.getBonds().stream().anyMatch(b -> b.get1_2(c) == C2s)) {
          C3s = c;
          C3s.setName("C3'");
          O3s = findBondedAtoms(C3s, 8).get(0);
          O3s.setName("O3'");
          findBondedAtoms(C3s, 1).get(0).setName("H3'");
          if (hasAttachedAtom(O3s, 1)) {
            findBondedAtoms(O3s, 1).get(0).setName("HO3'");
          } // Else, handle the possibility of 3'-P cap later.
        } else {
          C5s = c;
          C5s.setName("C5'");
          List<Atom> allH5List = findBondedAtoms(C5s, 1);
          Atom[] allH5s = allH5List.toArray(new Atom[0]);
          sortAtomsByDistance(O4s, allH5s);
          allH5s[0].setName("H5'");
          allH5s[1].setName("H5''");
        }
      }

      if (hasAttachedAtom(C2s, 8)) {
        Atom O2s = findBondedAtoms(C2s, 8).get(0);
        O2s.setName("O2'");
        findBondedAtoms(O2s, 1).get(0).setName("HO2'");
        findBondedAtoms(C2s, 1).get(0).setName("H2'");
      } else {
        List<Atom> bothH2List = findBondedAtoms(C2s, 1);
        Atom[] bothH2 = bothH2List.toArray(new Atom[0]);
        sortAtomsByDistance(H1s, bothH2);
        // Best-guess assignment, but is sometimes the other way around.
        bothH2[0].setName("H2''");
        bothH2[1].setName("H2'");
      }

      // logger.info(format(" C5\' null: %b", C5s == null));
      Atom O5s = findBondedAtoms(C5s, 8).get(0);
      O5s.setName("O5'");

      if (hasAttachedAtom(O5s, 1)) {
        findBondedAtoms(O5s, 1).get(0).setName("HO5'");
      } else if (hasAttachedAtom(O5s, 15)) {
        Atom P = findBondedAtoms(O5s, 15).get(0);
        P.setName("P");
        List<Atom> bondedO = findBondedAtoms(P, O5s, 8);
        List<Atom> thisResO =
            bondedO.stream()
                .filter(o -> residue.getAtomList().contains(o))
                .collect(Collectors.toList());
        int nBonded = bondedO.size();
        int nRes = thisResO.size();

        if (nBonded == 0) {
          // Do nothing.
        } else if (nBonded == nRes) {
          Atom OP1 = bondedO.get(0);
          OP1.setName("OP1");

          // OP2 is approximately +120 degrees from OP1, OP3 is -120 degrees.
          final double[] xyzC5s = C5s.getXYZ(new double[3]);
          final double[] xyzO5s = O5s.getXYZ(new double[3]);
          final double[] xyzP = P.getXYZ(new double[3]);
          final double[] xyzOP1 = OP1.getXYZ(new double[3]);
          double dihedral = DoubleMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzOP1);
          double twoPiOver3 = 2.0 * Math.PI / 3.0;
          double target = ScalarMath.modToRange(dihedral + twoPiOver3, -Math.PI, Math.PI);
          List<Atom> otherO =
              bondedO.stream()
                  .filter(o -> o != OP1)
                  .sorted(
                      Comparator.comparingDouble(
                          (Atom o) -> {
                            double[] xyzO = o.getXYZ(new double[3]);
                            double dihedO = DoubleMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzO);
                            double diff = dihedO - target;
                            double twoPi = 2 * Math.PI;
                            diff = ScalarMath.modToRange(diff, 0, twoPi);
                            diff = diff < Math.PI ? diff : twoPi - diff;
                            return diff;
                          }))
                  .collect(Collectors.toList());
          for (int i = 0; i < otherO.size(); i++) {
            otherO.get(i).setName(format("OP%d", i + 2));
          }
        } else {
          Atom nextO3s =
              bondedO.stream().filter(o -> !residue.getAtomList().contains(o)).findAny().get();

          // OP1 is approximately +120 degrees from next O3', OP2 is -120 degrees.
          final double[] xyzC5s = C5s.getXYZ(new double[3]);
          final double[] xyzO5s = O5s.getXYZ(new double[3]);
          final double[] xyzP = P.getXYZ(new double[3]);
          final double[] xyzNextO3s = nextO3s.getXYZ(new double[3]);
          double dihedral = DoubleMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzNextO3s);
          double twoPiOver3 = 2.0 * Math.PI / 3.0;
          double target = ScalarMath.modToRange(dihedral + twoPiOver3, -Math.PI, Math.PI);
          List<Atom> otherO =
              bondedO.stream()
                  .filter(o -> o != nextO3s)
                  .sorted(
                      Comparator.comparingDouble(
                          (Atom o) -> {
                            double[] xyzO = o.getXYZ(new double[3]);
                            double dihedO = DoubleMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzO);
                            double diff = dihedO - target;
                            double twoPi = 2 * Math.PI;
                            diff = ScalarMath.modToRange(diff, 0, twoPi);
                            diff = diff < Math.PI ? diff : twoPi - diff;
                            return diff;
                          }))
                  .collect(Collectors.toList());
          for (int i = 0; i < otherO.size(); i++) {
            otherO.get(i).setName(format("OP%d", i + 1));
          }
        }

        for (Atom op : bondedO) {
          if (hasAttachedAtom(op, 1)) {
            findBondedAtoms(op, 1).get(0).setName("H" + op.getName());
          }
        }
      }
      renameCommonNucleobase(N19, C1s, na3);
    } else {
      logger.warning(" Could not find O4' for residue " + residue);
    }
  }

  /**
   * Renames the atoms of the common nucleobases (A, C, G, T, U, and deoxy variants).
   *
   * @param N19 N1 of pyrimidines, N9 of purines.
   * @param C1s C1' of the ribose sugar.
   * @param na3 Identity of the nucleic acid.
   */
  public static void renameCommonNucleobase(
      Atom N19, Atom C1s, NucleicAcid3 na3) {
    switch (na3) {
      case ADE:
      case DAD:
        {
          Map<String, Atom> purineBase = renameCommonPurine(N19, C1s);
          // Unique to A: H2, N6, H6[12]
          findBondedAtoms(purineBase.get("C2"), 1).get(0).setName("H2");
          Atom C6 = purineBase.get("C6");
          Atom N1 = purineBase.get("N1");
          Atom N6 = findBondedAtoms(C6, N1, 7).get(0);
          N6.setName("N6");
          List<Atom> allH6List = findBondedAtoms(N6, 1);
          Atom[] allH6 = sortAtomsByDistance(N1, allH6List);
          allH6[0].setName("H61");
          allH6[1].setName("H62");
        }
        break;
      case CYT:
      case DCY:
        {
          Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
          // Unique to C: N4, H4[12]
          Atom C4 = pyrimidineBase.get("C4");
          Atom N3 = pyrimidineBase.get("N3");
          Atom N4 = findBondedAtoms(C4, N3, 7).get(0);
          N4.setName("N4");
          Atom[] allH4 = sortAtomsByDistance(N3, findBondedAtoms(N4, 1));
          allH4[0].setName("H41");
          allH4[1].setName("H42");
        }
        break;
      case GUA:
      case DGU:
        {
          Map<String, Atom> purineBase = renameCommonPurine(N19, C1s);
          // Unique to G: H1, N2, H2[12], O6
          Atom N1 = purineBase.get("N1");
          Atom C2 = purineBase.get("C2");
          Atom C6 = purineBase.get("C6");
          findBondedAtoms(N1, 1).get(0).setName("H1");
          Atom N2 =
              findBondedAtoms(C2, N1, 7).stream()
                  .filter(n -> hasAttachedAtom(n, 1))
                  .findAny()
                  .get();
          N2.setName("N2");
          Atom[] allH2 = sortAtomsByDistance(N1, findBondedAtoms(N2, 1));
          allH2[0].setName("H21");
          allH2[1].setName("H22");
          findBondedAtoms(C6, 8).get(0).setName("O6");
        }
        break;
      case URI:
        {
          Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
          // Unique to U: H3, O4
          findBondedAtoms(pyrimidineBase.get("N3"), 1).get(0).setName("H3");
          findBondedAtoms(pyrimidineBase.get("C4"), 8).get(0).setName("O4");
        }
        break;
      case THY:
      case DTY:
        {
          Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
          // Unique to T: H3, O4, C7
          findBondedAtoms(pyrimidineBase.get("N3"), 1).get(0).setName("H3");
          findBondedAtoms(pyrimidineBase.get("C4"), 8).get(0).setName("O4");
          Atom C5 = pyrimidineBase.get("C5");
          for (Atom c : findBondedAtoms(C5, 6)) {
            List<Atom> bondedH = findBondedAtoms(c, 1);
            if (bondedH != null && bondedH.size() == 3) {
              c.setName("C7");
              for (int i = 0; i < 3; i++) {
                bondedH.get(i).setName(format("H7%d", i + 1));
              }
              break;
            }
          }
        }
        break;
    }
  }

  /**
   * Renames atoms common to all standard purines (A, G)
   *
   * @param N9 The N9 atom.
   * @param C1s The C1' atom.
   * @return A Map containing Atoms important to finding and naming base-unique atoms.
   */
  public static Map<String, Atom> renameCommonPurine(Atom N9, Atom C1s) {
    Map<String, Atom> keyAtoms = new HashMap<>(10);
    N9.setName("N9");
    for (Atom c : findBondedAtoms(N9, C1s, 6)) {
      if (hasAttachedAtom(c, 1)) {
        Atom C8 = c;
        C8.setName("C8");
        findBondedAtoms(C8, 1).get(0).setName("H8");
        Atom N7 = findBondedAtoms(C8, N9, 7).get(0);
        N7.setName("N7");
        Atom C5 = findBondedAtoms(N7, C8, 6).get(0);
        C5.setName("C5");
      } else {
        Atom C4 = c;
        C4.setName("C4");
        Atom N3 = findBondedAtoms(C4, N9, 7).get(0);
        N3.setName("N3");
        Atom C2 = findBondedAtoms(N3, C4, 6).get(0);
        C2.setName("C2");
        keyAtoms.put("C2", C2);
        Atom N1 = findBondedAtoms(C2, N3, 7).get(0);
        N1.setName("N1"); // And not, say, "largest non-nuclear explosion ever".
        keyAtoms.put("N1", N1);
        Atom C6 = findBondedAtoms(N1, C2, 6).get(0);
        C6.setName("C6");
        keyAtoms.put("C6", C6);
      }
    }
    /* Common atoms: N1, C2, N3, C4, C5, C6, N7, C8, H8, N9. */
    return keyAtoms;
  }

  /**
   * Renames atoms common to all standard pyrimidines (C, T, U)
   *
   * @param N1 The N1 atom.
   * @param C1s The C1' atom.
   * @return A Map containing Atoms important to finding and naming base-unique atoms.
   */
  public static Map<String, Atom> renameCommonPyrimidine(Atom N1, Atom C1s) {
    Map<String, Atom> keyAtoms = new HashMap<>();
    N1.setName("N1");
    for (Atom c : findBondedAtoms(N1, C1s, 6)) {
      if (hasAttachedAtom(c, 8)) {
        Atom C2 = c;
        C2.setName("C2");
        findBondedAtoms(C2, 8).get(0).setName("O2");
        Atom N3 = findBondedAtoms(C2, N1, 7).get(0);
        N3.setName("N3");
        keyAtoms.put("N3", N3);
        Atom C4 = findBondedAtoms(N3, C2, 6).get(0);
        C4.setName("C4");
        keyAtoms.put("C4", C4);
        Atom C5 = findBondedAtoms(C4, 6).get(0);
        C5.setName("C5");
        keyAtoms.put("C5", C5);
        if (hasAttachedAtom(C5, 1)) {
          findBondedAtoms(C5, 1).get(0).setName("H5");
        }
      } else {
        Atom C6 = c;
        C6.setName("C6");
        findBondedAtoms(C6, 1).get(0).setName("H6");
      }
    }

    // Common atoms: N1, C2, O2, N3, C4, C5, C6, H6
    return keyAtoms;
  }

  /**
   * renameDeltaHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   * @param indexes a int.
   */
  public static void renameDeltaHydrogens(Residue residue, List<Atom> resAtoms, int indexes) {
    Atom[] HDn = new Atom[3];
    switch (indexes) {
      case 12:
        HDn[0] = (Atom) residue.getAtomNode("HD1");
        HDn[1] = (Atom) residue.getAtomNode("HD2");
        break;
      case 13:
        HDn[0] = (Atom) residue.getAtomNode("HD1");
        HDn[2] = (Atom) residue.getAtomNode("HD3");
        break;
      case 23:
        HDn[1] = (Atom) residue.getAtomNode("HD2");
        HDn[2] = (Atom) residue.getAtomNode("HD3");
        break;
      default:
        return;
    }
    for (Atom HDatom : HDn) {
      resAtoms.remove(HDatom);
    }
    if (!resAtoms.isEmpty() && HDn[0] == null && (indexes == 12 || indexes == 13)) {
      resAtoms.get(0).setName("HD1");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HDn[1] == null && (indexes == 12 || indexes == 23)) {
      resAtoms.get(0).setName("HD2");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HDn[2] == null && (indexes == 13 || indexes == 23)) {
      resAtoms.get(0).setName("HD3");
      resAtoms.remove(0);
    }
  }

  /**
   * renameEpsilonHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   * @param indexes a int.
   */
  public static void renameEpsilonHydrogens(Residue residue, List<Atom> resAtoms, int indexes) {
    Atom[] HEn = new Atom[3];
    switch (indexes) {
      case 12:
        HEn[0] = (Atom) residue.getAtomNode("HE1");
        HEn[1] = (Atom) residue.getAtomNode("HE2");
        break;
      case 13:
        HEn[0] = (Atom) residue.getAtomNode("HE1");
        HEn[2] = (Atom) residue.getAtomNode("HE3");
        break;
      case 23:
        HEn[1] = (Atom) residue.getAtomNode("HE2");
        HEn[2] = (Atom) residue.getAtomNode("HE3");
        break;
      default:
        return;
    }
    for (Atom HEatom : HEn) {
      resAtoms.remove(HEatom);
    }
    if (!resAtoms.isEmpty() && HEn[0] == null && (indexes == 12 || indexes == 13)) {
      resAtoms.get(0).setName("HE1");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HEn[1] == null && (indexes == 12 || indexes == 23)) {
      resAtoms.get(0).setName("HE2");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HEn[2] == null && (indexes == 13 || indexes == 23)) {
      resAtoms.get(0).setName("HE3");
      resAtoms.remove(0);
    }
  }

  /**
   * renameGammaHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   * @param indexes a int.
   */
  public static void renameGammaHydrogens(Residue residue, List<Atom> resAtoms, int indexes) {
    Atom[] HGn = new Atom[3];
    switch (indexes) {
      case 12:
        HGn[0] = (Atom) residue.getAtomNode("HG1");
        HGn[1] = (Atom) residue.getAtomNode("HG2");
        break;
      case 13:
        HGn[0] = (Atom) residue.getAtomNode("HG1");
        HGn[2] = (Atom) residue.getAtomNode("HG3");
        break;
      case 23:
        HGn[1] = (Atom) residue.getAtomNode("HG2");
        HGn[2] = (Atom) residue.getAtomNode("HG3");
        break;
      default:
        return;
    }
    for (Atom HGatom : HGn) {
      resAtoms.remove(HGatom);
    }
    if (!resAtoms.isEmpty() && HGn[0] == null && (indexes == 12 || indexes == 13)) {
      resAtoms.get(0).setName("HG1");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HGn[1] == null && (indexes == 12 || indexes == 23)) {
      resAtoms.get(0).setName("HG2");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HGn[2] == null && (indexes == 13 || indexes == 23)) {
      resAtoms.get(0).setName("HG3");
      resAtoms.remove(0);
    }
  }

  /**
   * renameGlutamineHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   */
  public static void renameGlutamineHydrogens(Residue residue, List<Atom> resAtoms) {
    Atom HE21 = (Atom) residue.getAtomNode("HE21");
    Atom HE22 = (Atom) residue.getAtomNode("HE22");
    if (HE21 != null) {
      resAtoms.remove(HE21);
    }
    if (HE22 != null) {
      resAtoms.remove(HE22);
    }
    if (!resAtoms.isEmpty() && HE21 == null) {
      resAtoms.get(0).setName("HE21");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HE22 == null) {
      resAtoms.get(0).setName("HE21");
    }
  }

  /**
   * renameGlycineAlphaHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   */
  public static void renameGlycineAlphaHydrogens(Residue residue, List<Atom> resAtoms) {
    Atom HA2 = (Atom) residue.getAtomNode("HA2");
    Atom HA3 = (Atom) residue.getAtomNode("HA3");
    if (HA2 != null) {
      resAtoms.remove(HA2);
    }
    if (HA3 != null) {
      resAtoms.remove(HA3);
    }
    if (HA2 == null && !resAtoms.isEmpty()) {
      resAtoms.get(0).setName("HA2");
      resAtoms.remove(0);
    }
    if (HA3 == null && !resAtoms.isEmpty()) {
      resAtoms.get(0).setName("HA3");
    }
  }

  /**
   * renameIsoleucineHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   */
  public static void renameIsoleucineHydrogens(Residue residue, List<Atom> resAtoms) {
    Atom HG12 = (Atom) residue.getAtomNode("HG12");
    Atom HG13 = (Atom) residue.getAtomNode("HG13");
    if (HG12 != null) {
      resAtoms.remove(HG12);
    }
    if (HG13 != null) {
      resAtoms.remove(HG13);
    }
    if (HG12 == null && !resAtoms.isEmpty()) {
      resAtoms.get(0).setName("HG12");
      resAtoms.remove(0);
    }
    if (HG13 == null && !resAtoms.isEmpty()) {
      resAtoms.get(0).setName("HG13");
    }
  }

  /**
   * renameNTerminusHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   */
  public static void renameNTerminusHydrogens(Residue residue) {
    Atom[] h = new Atom[3];
    h[0] = (Atom) residue.getAtomNode("H1");
    h[1] = (Atom) residue.getAtomNode("H2");
    h[2] = (Atom) residue.getAtomNode("H3");
    int numAtoms = 0;
    for (Atom atom : h) {
      numAtoms += (atom == null ? 0 : 1);
    }
    if (numAtoms == 3) {
      return;
    }
    List<Atom> resAtoms = residue.getAtomList();
    for (Atom resAtom : resAtoms) {
      // Check if already contained in h[].
      boolean doContinue = false;
      for (Atom hAtom : h) {
        if (resAtom.equals(hAtom)) {
          doContinue = true;
          break;
        }
      }
      if (doContinue) {
        continue;
      }

      // If the hydrogen matches H or H[1-3], assign to first null h entity.
      String atomName = resAtom.getName().toUpperCase();
      if (atomName.equals("H") || atomName.matches("H[1-3]") || atomName.matches("[1-3]H")) {
        ++numAtoms;
        for (int i = 0; i < h.length; i++) {
          if (h[i] == null) {
            resAtom.setName("H" + (i + 1));
            h[i] = resAtom;
            break;
          }
        }
        if (numAtoms == 3) {
          return;
        }
      }
    }
  }

  /**
   * Renames the Atoms in a nucleic acid to PDB standard.
   *
   * @param residue Residue to fix atom names of.
   */
  public static void renameNucleicAcidToPDBStandard(Residue residue) {
    if (residue.getChainID() == null) {
      residue.setChainID('Z');
    }
    assert residue.getResidueType() == Residue.ResidueType.NA;
    NucleicAcid3 na3 = residue.getNucleicAcid3(true);
    residue.setName(na3.toString());
    switch (na3) {
      case ADE:
      case DAD:
      case CYT:
      case DCY:
      case GUA:
      case DGU:
      case THY:
      case DTY:
      case URI:
        renameCommonNucleicAcid(residue, na3);
        break;
      default:
        logger.info(" Could not rename atoms for nonstandard nucleic acid " + na3);
        break;
    }
  }

  /**
   * renameZetaHydrogens.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param resAtoms a {@link java.util.List} object.
   * @param indexes a int.
   */
  public static void renameZetaHydrogens(Residue residue, List<Atom> resAtoms, int indexes) {
    Atom[] HZn = new Atom[3];
    switch (indexes) {
      case 12:
        HZn[0] = (Atom) residue.getAtomNode("HZ1");
        HZn[1] = (Atom) residue.getAtomNode("HZ2");
        break;
      case 13:
        HZn[0] = (Atom) residue.getAtomNode("HZ1");
        HZn[2] = (Atom) residue.getAtomNode("HZ3");
        break;
      case 23:
        HZn[1] = (Atom) residue.getAtomNode("HZ2");
        HZn[2] = (Atom) residue.getAtomNode("HZ3");
        break;
      default:
        return;
    }
    for (Atom HZatom : HZn) {
      resAtoms.remove(HZatom);
    }
    if (!resAtoms.isEmpty() && HZn[0] == null && (indexes == 12 || indexes == 13)) {
      resAtoms.get(0).setName("HZ1");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HZn[1] == null && (indexes == 12 || indexes == 23)) {
      resAtoms.get(0).setName("HZ2");
      resAtoms.remove(0);
    }
    if (!resAtoms.isEmpty() && HZn[2] == null && (indexes == 13 || indexes == 23)) {
      resAtoms.get(0).setName("HZ3");
      resAtoms.remove(0);
    }
  }

  /** Common HETATOM labels for water and ions. */
  public enum HetAtoms {
    BR,
    CA,
    CA2,
    CL,
    K,
    MG,
    MG2,
    NA,
    HOH,
    H2O,
    WAT,
    ZN,
    ZN2;

    /**
     * Slightly more robust parsing function that ignores case and trailing numbers, -, and +
     *
     * @param str String to parse
     * @return Corresponding HetAtoms value.
     */
    public static HetAtoms parse(String str) {
      String hName = str.toUpperCase().replaceFirst("[0-9+\\-]+$", "");
      return valueOf(hName);
    }
  }
}
