// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.parameters;

import static ffx.potential.parameters.ForceField.ForceFieldType.BIOTYPE;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import java.util.Comparator;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The BioType class maps PDB identifiers to atom types.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class BioType extends BaseType implements Comparator<String> {

  /** A Logger for the BioType class. */
  private static final Logger logger = Logger.getLogger(BioType.class.getName());
  /** The PDB atom name for this BioType. */
  public final String atomName;
  /** The PDB molecule name for this BioType. */
  public final String moleculeName;
  /** Bonds are required to listed atom names. */
  public final String[] bonds;
  /** The index of this BioType. */
  public int index;
  /** The force field atom type to be used for the molecule / atom name combination. */
  public int atomType;

  /**
   * BioType Constructor.
   *
   * @param index int
   * @param atomName String
   * @param moleculeName String
   * @param atomType int
   * @param bonds an array of {@link java.lang.String} objects.
   */
  public BioType(int index, String atomName, String moleculeName, int atomType, String[] bonds) {
    super(BIOTYPE, Integer.toString(index));
    this.index = index;
    this.atomName = atomName;
    if (moleculeName != null) {
      this.moleculeName = moleculeName.replace(',', ' ').replace('"', ' ').trim();
    } else {
      this.moleculeName = null;
    }
    this.atomType = atomType;
    this.bonds = bonds;
  }

  /**
   * Construct an BioType from an input string.
   *
   * @param input The overall input String.
   * @param tokens The input String tokenized.
   * @return an BioType instance.
   */
  public static BioType parse(String input, String[] tokens) {
    if (tokens.length < 5) {
      logger.log(Level.WARNING, "Invalid BIOTYPE type:\n{0}", input);
    } else {
      try {
        int index = parseInt(tokens[1]);
        String atomName = tokens[2];
        // The "residue" string may contain spaces,
        // and is therefore surrounded in quotes located at "first" and
        // "last".
        int first = input.indexOf("\"");
        int last = input.lastIndexOf("\"");
        if (first >= last) {
          logger.log(Level.WARNING, "Invalid BIOTYPE type:\n{0}", input);
          return null;
        }
        // Environment
        String moleculeName = input.substring(first, last + 1).intern();
        // Shrink the tokens array to only include entries
        // after the environment field.
        tokens = input.substring(last + 1).trim().split(" +");
        int atomType = parseInt(tokens[0]);
        int bondCount = tokens.length - 1;
        String[] bonds = null;
        if (bondCount > 0) {
          bonds = new String[bondCount];
          arraycopy(tokens, 1, bonds, 0, bondCount);
        }
        return new BioType(index, atomName, moleculeName, atomType, bonds);
      } catch (NumberFormatException e) {
        String message = "Exception parsing BIOTYPE type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public int compare(String s1, String s2) {
    int t1 = parseInt(s1);
    int t2 = parseInt(s2);
    return Integer.compare(t1, t2);
  }

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    BioType bioType = (BioType) o;
    return bioType.index == this.index;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(index);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted biotype.
   */
  @Override
  public String toString() {
    StringBuilder sb =
        new StringBuilder(
            format("biotype  %5d  %-4s  \"%-23s\"  %5d", index, atomName, moleculeName, atomType));
    if (bonds != null && bonds.length > 0) {
      for (String bond : bonds) {
        sb.append(format("  %-4s", bond));
      }
    }
    return sb.toString();
  }

  /**
   * incrementIndexAndType
   *
   * @param indexIncrement a int.
   * @param typeIncrement a int.
   */
  void incrementIndexAndType(int indexIncrement, int typeIncrement) {
    index += indexIncrement;
    atomType += typeIncrement;
    setKey(Integer.toString(index));
  }
}
