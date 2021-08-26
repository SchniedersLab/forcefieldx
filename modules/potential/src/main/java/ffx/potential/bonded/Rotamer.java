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

import static org.apache.commons.math3.util.FastMath.max;

import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.NucleicAcidUtils.NucleicAcid3;
import ffx.potential.parameters.TitrationUtils;

/**
 * The Rotamer Class usually represents one immutable amino acid Rotamer.
 * <p>
 * It is additionally being extended to represent one nucleic acid Rotamer.
 *
 * @author Ava M. Lynn
 * @author Jacob M. Litman
 * @since 1.0
 */
public class Rotamer {

  /** Torsions chi 1-4 are used for amino acids and nucleic acids. */
  public final double chi1;
  public final double chi2;
  public final double chi3;
  public final double chi4;
  /** Torsions chi 5-7 are only currently used for nucleic acids. */
  public final double chi5;
  final double chi6;
  final double chi7;
  /**
   * An array of chi angles for this rotamer.
   */
  public final double[] angles;
  /**
   * An array of sigmas for each chi angle.
   */
  public final double[] sigmas;
  /**
   * Number of chi/sigma values.
   */
  public final int length;
  /**
   * Residue state used to initialize the rotamer.
   */
  public ResidueState originalState;
  /**
   * Flag to indicate the rotamer was initialized from a Residue state.
   */
  public boolean isState;
  /**
   * The A.A. name of this residue (or null for a N.A.).
   */
  public AminoAcid3 aminoAcid3;
  /**
   * The N.A. name of this residue (or null for a A.A.).
   */
  public NucleicAcid3 nucleicAcid3;
  /**
   * If this flag is set, application of a rotamer requires updating force field parameters.
   */
  public boolean isTitrating;
  /**
   * The TitrationUtils handles application of rotamer specific force field parameters.
   */
  private TitrationUtils titrationUtils = null;

  /**
   * Constructor for unknown residue types.
   *
   * @param values a double.
   */
  public Rotamer(double... values) {
    length = values.length / 2;
    angles = new double[max(length, 7)];
    sigmas = new double[max(length, 7)];
    nucleicAcid3 = null;
    aminoAcid3 = null;
    for (int i = 0; i < values.length / 2; i++) {
      int ii = 2 * i;
      angles[i] = values[ii];
      sigmas[i] = values[ii + 1];
    }
    chi1 = angles[0];
    chi2 = angles[1];
    chi3 = angles[2];
    chi4 = angles[3];
    chi5 = angles[4];
    chi6 = angles[5];
    chi7 = angles[6];
    originalState = null;
    isState = false;
    isTitrating = false;
  }

  /**
   * Constructor for Rotamer.
   *
   * @param aminoAcid3 a {@link AminoAcid3} object.
   * @param values a double.
   */
  public Rotamer(AminoAcid3 aminoAcid3, double... values) {
    this(values);
    this.aminoAcid3 = aminoAcid3;
  }

  /**
   * Constructor for Rotamer.
   *
   * @param nucleicAcid3 a {@link NucleicAcid3} object.
   * @param values a double.
   */
  public Rotamer(NucleicAcid3 nucleicAcid3, double... values) {
    this(values);
    this.nucleicAcid3 = nucleicAcid3;
  }

  /**
   * Constructor for Rotamer.
   *
   * @param aminoAcid3 a {@link AminoAcid3} object.
   * @param titrationUtils Use to apply rotamer specific force field parameters.
   * @param values a double.
   */
  public Rotamer(AminoAcid3 aminoAcid3, TitrationUtils titrationUtils, double... values) {
    this(aminoAcid3, values);
    if (titrationUtils != null) {
      this.isTitrating = true;
      this.titrationUtils = titrationUtils;
    } else {
      this.isTitrating = false;
    }
  }

  /**
   * Constructor for unknown residue types.
   *
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(ResidueState residueState, double... values) {
    this(values);
    isState = true;
    originalState = residueState;
  }

  /**
   * Constructor for unknown residue types.
   *
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(ResidueState residueState, TitrationUtils titrationUtils, double... values) {
    this(residueState, values);
    if (titrationUtils != null) {
      this.titrationUtils = titrationUtils;
      isTitrating = true;
    }
  }

  /**
   * Constructor for Rotamer.
   *
   * @param aminoAcid3 a {@link AminoAcid3} object.
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(AminoAcid3 aminoAcid3, ResidueState residueState, double... values) {
    this(residueState, values);
    this.aminoAcid3 = aminoAcid3;
  }

  /**
   * Constructor for Rotamer.
   *
   * @param aminoAcid3 a {@link AminoAcid3} object.
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(AminoAcid3 aminoAcid3, ResidueState residueState, TitrationUtils titrationUtils,
      double... values) {
    this(aminoAcid3, residueState, values);
    if (titrationUtils != null) {
      this.titrationUtils = titrationUtils;
      isTitrating = true;
    }
  }

  /**
   * Constructor for Rotamer.
   *
   * @param nucleicAcid3 a {@link NucleicAcid3} object.
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(NucleicAcid3 nucleicAcid3, ResidueState residueState, double... values) {
    this(residueState, values);
    this.nucleicAcid3 = nucleicAcid3;
  }

  /**
   * Constructor for Rotamer.
   *
   * @param nucleicAcid3 a {@link NucleicAcid3} object.
   * @param residueState a {@link ffx.potential.bonded.ResidueState} object.
   * @param values a double.
   */
  public Rotamer(NucleicAcid3 nucleicAcid3, ResidueState residueState, TitrationUtils titrationUtils,
      double... values) {
    this(nucleicAcid3, residueState, values);
    if (titrationUtils != null) {
      this.titrationUtils = titrationUtils;
      isTitrating = true;
    }
  }

  /**
   * Update force field parameters for force field dependent rotamers.
   *
   * @param residue Residue to update.
   */
  public void updateParameters(Residue residue) {
    titrationUtils.updateResidueParameters(residue, this);
  }

  /**
   * Factory method to construct an original-coordinates Rotamer from a residue.
   *
   * @param residue Residue to construct a default rotamer for.
   * @return Rotamer based on the coordinates of the residue.
   */
  public static Rotamer[] defaultRotamerFactory(Residue residue) {
    return defaultRotamerFactory(residue, null);
  }

  /**
   * Factory method to construct an original-coordinates Rotamer from a residue.
   *
   * @param residue Residue to construct a default rotamer for.
   * @return Rotamer based on the coordinates of the residue.
   */
  public static Rotamer[] defaultRotamerFactory(Residue residue, TitrationUtils titrationUtils) {
    ResidueState resState = residue.storeState();
    double[] chi = RotamerLibrary.measureRotamer(residue, false);

    double[] vals = new double[chi.length * 2];
    for (int i = 0; i < chi.length; i++) {
      int index = i * 2;
      vals[index] = chi[i];
      vals[index + 1] = 0.0;
    }

    switch (residue.getResidueType()) {
      case AA:
        // Only one rotamer for non-titrating cases.
        if (titrationUtils == null) {
          Rotamer[] rotamers = new Rotamer[1];
          rotamers[0] = new Rotamer(residue.getAminoAcid3(), resState, titrationUtils, vals);
          return rotamers;
        }
        switch (residue.getAminoAcid3()) {
          case ASH:
          case ASP:
            Rotamer[] rotamers = new Rotamer[2];
            rotamers[0] = new Rotamer(AminoAcid3.ASP, resState, titrationUtils, vals);
            rotamers[1] = new Rotamer(AminoAcid3.ASH, resState, titrationUtils, vals);
            return rotamers;
          case GLH:
          case GLU:
            rotamers = new Rotamer[2];
            rotamers[0] = new Rotamer(AminoAcid3.GLU, resState, titrationUtils, vals);
            rotamers[1] = new Rotamer(AminoAcid3.GLH, resState, titrationUtils, vals);
            return rotamers;
          case HID:
          case HIE:
          case HIS:
            rotamers = new Rotamer[3];
            rotamers[0] = new Rotamer(AminoAcid3.HIS, resState, titrationUtils, vals);
            rotamers[1] = new Rotamer(AminoAcid3.HID, resState, titrationUtils, vals);
            rotamers[2] = new Rotamer(AminoAcid3.HIE, resState, titrationUtils, vals);
            return rotamers;
          case LYS:
          case LYD:
            rotamers = new Rotamer[2];
            rotamers[0] = new Rotamer(AminoAcid3.LYS, resState, titrationUtils, vals);
            rotamers[1] = new Rotamer(AminoAcid3.LYD, resState, titrationUtils, vals);
            return rotamers;
          default:
            rotamers = new Rotamer[1];
            rotamers[0] = new Rotamer(residue.getAminoAcid3(), resState, titrationUtils, vals);
            return rotamers;
        }
      case NA:
        Rotamer[] rotamers = new Rotamer[1];
        rotamers[0] = new Rotamer(residue.getNucleicAcid3(), resState, titrationUtils, vals);
        return rotamers;
      case UNK:
      default:
        rotamers = new Rotamer[1];
        rotamers[0] = new Rotamer(resState, titrationUtils, vals);
        return rotamers;
    }
  }

  /**
   * toAngleString.
   *
   * @return a {@link java.lang.String} object.
   */
  public String toAngleString() {
    StringBuilder sb = new StringBuilder();
    int n = max(4, length);
    for (int i = 0; i < n; i++) {
      sb.append(String.format(" %6.1f %4.1f", angles[i], sigmas[i]));
    }
    return sb.toString();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder(getName());
    int n = max(4, length);
    for (int i = 0; i < n; i++) {
      sb.append(String.format(" %6.1f %4.1f", angles[i], sigmas[i]));
    }
    return sb.toString();
  }

  public String getName() {
    if (aminoAcid3 != null) {
      return aminoAcid3.toString();
    } else if (nucleicAcid3 != null) {
      return nucleicAcid3.toString();
    } else {
      return "";
    }
  }

  public double getRotamerPhBias() {
    return titrationUtils.getRotamerPhBias(aminoAcid3);
  }
}
