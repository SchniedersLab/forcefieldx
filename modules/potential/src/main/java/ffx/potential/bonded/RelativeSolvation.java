// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.RelativeSolvationType;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * A relative solvation term for chemical perturbations.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class RelativeSolvation {

  private static final Logger logger = Logger.getLogger(RelativeSolvation.class.getName());
  /** Look-up of non-standard energies. */
  private final Map<String, Double> nonStdEnergies;
  /** Solvation library in use. */
  private final SolvationLibrary solvationLibrary;

  /**
   * Constructor for RelativeSolvation.
   *
   * @param solvationLibrary a {@link ffx.potential.bonded.RelativeSolvation.SolvationLibrary}
   *     object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public RelativeSolvation(SolvationLibrary solvationLibrary, ForceField forceField) {
    this.solvationLibrary = solvationLibrary;
    nonStdEnergies = new HashMap<>();
    for (RelativeSolvationType rsType : forceField.getRelativeSolvationTypes().values()) {
      String resName = rsType.getResName();
      double e = rsType.getSolvEnergy();
      if (nonStdEnergies.put(resName, e) != null) {
        logger.warning(format(" Repeat relative solvation for %s", resName));
      }
    }
  }

  /**
   * Gets the solvation energy (de-solvation penalty) for a given residue, allowing for sequence
   * optimization to include an estimate of energy relative to the unfolded state.
   *
   * @param residue Residue to check
   * @param checkZeroes Throws an error if not in solvation energy library
   * @return Solvation energy
   * @throws java.lang.IllegalArgumentException if any.
   * @throws java.lang.IllegalArgumentException if any.
   */
  public double getSolvationEnergy(Residue residue, boolean checkZeroes)
      throws IllegalArgumentException {
    String resName = "";
    double energy;
    Residue theRes =
        (residue instanceof MultiResidue) ? ((MultiResidue) residue).getActive() : residue;
    switch (theRes.getResidueType()) {
      case AA -> {
        if (theRes instanceof MultiResidue) {
          resName = ((MultiResidue) theRes).getActive().getName();
        } else {
          resName = theRes.getName();
        }
        energy = getAASolvationEnergy(theRes);
      }
      case NA -> {
        if (theRes instanceof MultiResidue) {
          resName = ((MultiResidue) theRes).getActive().getName();
        } else {
          resName = theRes.getName();
        }
        energy = getNASolvationEnergy(theRes);
      }
      default -> energy = 0;
    }
    if (checkZeroes && energy == 0) {
      throw new IllegalArgumentException(
          format(" Zero de-solvation energy for residue %s: likely not in solvation library.",
              resName));
    }
    return energy;
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return "Relative solvation library: " + solvationLibrary.toString();
  }

  /**
   * Returns amino acid solvation energy based on solvation library.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Solvation energy
   */
  private double getAASolvationEnergy(Residue residue) {
    return switch (solvationLibrary) {
      case WOLFENDEN -> getWolfendenSolvationEnergy(residue);
      case CABANI -> getCabaniSolvationEnergy(residue);
      case EXPLICIT -> getExplicitSolvationEnergy(residue);
      case GK -> getGKSolvationEnergy(residue);
      case MACCALLUM_SPC -> getMacCallumSPCSolvationEnergy(residue);
      case MACCALLUM_TIP4P -> getMacCallumTIP4PSolvationEnergy(residue);
      default -> 0;
    };
  }

  /**
   * Will return relative solvation energies for nucleic acids; currently returns 0.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Relative solvation energy
   */
  private double getNASolvationEnergy(Residue residue) {
    return 0;
  }

  /**
   * Will return solvation energies relative to glycine for capped monomers in GK solvent; currently
   * wraps getExplicitSolvationEnergy.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Relative solvation energy
   */
  private double getGKSolvationEnergy(Residue residue) {
    return getExplicitSolvationEnergy(residue);
  }

  /**
   * Will return solvation energies relative to glycine for capped monomers in AMOEBA solvent;
   * currently approximates this with charging energies in AMOEBA solvent.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Relative solvation energy
   */
  private double getExplicitSolvationEnergy(Residue residue) {
    AminoAcid3 name = residue.getAminoAcid3();
    return switch (name) {
      case ALA -> 0.58;
      case CYS -> -0.85;
      case CYD -> -69.82;
      case ASP -> -69.45;
      case ASH -> -4.00;
      case GLU -> -71.40;
      case GLH -> -3.61;
      case PHE -> -1.59;
      case GLY -> 0.67;
      case HIS -> -45.36;
      case HID -> -8.42;
      case HIE -> -7.53;
      case ILE -> 0.14;
      case LYS -> -43.98;
      case LYD -> +0.35;
      case MET -> -3.48;
      case ASN -> -5.89;
      case PRO -> 7.82;
      case GLN -> -6.89;
      case ARG -> -42.57;
      case SER -> -2.14;
      case THR -> 0.58;
      case VAL -> 0.10;
      case TRP -> -4.64;
      case TYR -> 1.76;
      case TYD -> -41.71;
      case UNK -> nonStdEnergies.getOrDefault(residue.getName().toUpperCase(), 0.0);
      default -> 0;
    };
  }

  /**
   * Returns absolute solvation energies for side chain analogs as calculated by MacCallum for OPLS
   * in TIP4P solvent.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Solvation energy
   */
  private double getMacCallumTIP4PSolvationEnergy(Residue residue) {
    AminoAcid3 name = residue.getAminoAcid3();
    return switch (name) {
      case ALA -> 9.8;
      case CYS -> -0.5;
      case ASP -> -30.5;
      case GLU -> -19.0;
      case PHE -> -1.2;
      case HIS -> -28.0;
      case ILE -> 12.2;
      case LYS -> -13.6;
      case LEU -> 13.7;
      case MET -> -7.1;
      case ASN -> -34.5;
      case GLN -> -31.4;
      case ARG -> -43.9;
      case SER -> -20.2;
      case THR -> -20.3;
      case VAL -> 12.0;
      case TRP -> -16.2;
      case TYR -> -18.8;
      case UNK -> nonStdEnergies.getOrDefault(residue.getName().toUpperCase(), 0.0);
      default -> 0;
    };
  }

  /**
   * Returns absolute solvation energies for side chain analogs as calculated by MacCallum for OPLS
   * in SPC solvent.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Solvation energy
   */
  private double getMacCallumSPCSolvationEnergy(Residue residue) {
    AminoAcid3 name = residue.getAminoAcid3();
    return switch (name) {
      case ALA -> 9.3;
      case CYS -> -1.1;
      case ASP -> -30.1;
      case GLU -> -18.8;
      case PHE -> -1.4;
      case HIS -> -27.2;
      case ILE -> 11.9;
      case LYS -> -8.6;
      case LEU -> 12.6;
      case MET -> -5.1;
      case ASN -> -34.3;
      case GLN -> -30.8;
      case ARG -> -46.3;
      case SER -> -18.5;
      case THR -> -19.3;
      case VAL -> 11.3;
      case TRP -> -15.1;
      case TYR -> -18.2;
      case UNK -> nonStdEnergies.getOrDefault(residue.getName().toUpperCase(), 0.0);
      default -> 0;
    };
  }

  /**
   * Returns absolute solvation energies for side chain analogs as experimentally measured by Cabani
   * et al.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Solvation energy
   */
  private double getCabaniSolvationEnergy(Residue residue) {
    AminoAcid3 name = residue.getAminoAcid3();
    switch (name) {
      case ALA:
        return 8.4;
      case CYS:
        return -5.2;
      case ASP:
        return -28.1;
      case GLU:
        return -27.1;
      case PHE:
        return -3.7;
      case HIS:
        return -27.4;
      case ILE:
        return 8.7;
      case LYS:
        return -15.5;
      case LEU:
        return 9.7;
      case MET:
        return 9.0;
      case ASN:
        return -40.6;
      case GLN:
        return -18.7;
      case ARG:
        return -30.1;
      case SER:
        return -21.4;
      case THR:
        return -21.0;
      case VAL:
        return 8.2;
      case TRP:
        return -12.3;
      case TYR:
        return -25.7;
      case UNK:
        return nonStdEnergies.getOrDefault(residue.getName().toUpperCase(), 0.0);
      case GLY:
      case PRO:
      default:
        return 0;
    }
  }

  /**
   * Returns absolute solvation energies for side chain analogs as experimentally measured by
   * Wolfenden et al.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return Solvation energy
   */
  private double getWolfendenSolvationEnergy(Residue residue) {
    AminoAcid3 name = residue.getAminoAcid3();
    switch (name) {
      case ALA:
        return 8.1;
      case CYS:
        return -5.1;
      case ASP:
        return -27.5;
      case GLU:
        return -26.6;
      case PHE:
        return -3.1;
      case HIS:
        return -42.1;
      case ILE:
        return 8.8;
      case LYS:
        return -18.0;
      case LEU:
        return 9.4;
      case MET:
        return -6.1;
      case ASN:
        return -39.9;
      case GLN:
        return -38.7;
      case ARG:
        return -44.8;
      case SER:
        return -20.8;
      case THR:
        return -20.1;
      case VAL:
        return 8.2;
      case TRP:
        return -24.3;
      case TYR:
        return -25.2;
      case UNK:
        return nonStdEnergies.getOrDefault(residue.getName().toUpperCase(), 0.0);
      case GLY:
      case PRO:
      default:
        return 0; // Not listed.
    }
  }

  /**
   * Citations: Wolfenden et al: Wolfenden, R., Andersson, L., Cullis, P. M. and Southgate, C. C. B.
   * (1981) AFFINITIES OF AMINO-ACID SIDE-CHAINS FOR SOLVENT WATER. Biochemistry. 20, 849-855
   *
   * <p>Cabani et al: Cabani, S., Gianni, P., Mollica, V. and Lepori, L. (1981) GROUP
   * CONTRIBUTIONS TO THE THERMODYNAMIC PROPERTIES OF NON-IONIC ORGANIC SOLUTES IN DILUTE
   * AQUEOUS-SOLUTION. Journal of Solution Chemistry. 10, 563-595
   *
   * <p>MacCallum OPLS libraries: Maccallum, J. L. and Tieleman, D. P. (2003) Calculation of the
   * water-cyclohexane transfer free energies of neutral amino acid side-chain analogs using the OPLS
   * all-atom force field. Journal of Computational Chemistry. 24, 1930-1935
   */
  public enum SolvationLibrary {
    WOLFENDEN, CABANI, EXPLICIT, GK, MACCALLUM_SPC, MACCALLUM_TIP4P, OPLS_EXPLICIT, OPLS_GK, AUTO, NONE
  }
}
