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
package ffx.potential.parameters;

import static ffx.potential.parameters.ForceField.ForceFieldType.POLARIZE;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.pow;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The PolarizeType class defines an isotropic atomic polarizability.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class PolarizeType extends BaseType implements Comparator<String> {

  private static final Logger logger = Logger.getLogger(PolarizeType.class.getName());

  private static final double sixth = 1.0 / 6.0;
  /** Thole damping factor. */
  public final double thole;
  /** Value of polarizability scale factor. */
  public final double pdamp;
  /** Isotropic polarizability in units of Angstroms^3. */
  public final double polarizability;
  /** Atom type number. */
  public int type;
  /** Connected types in the polarization group of each atom. (may be null) */
  public int[] polarizationGroup;

  /**
   * PolarizeType Constructor.
   *
   * @param atomType The atom type.
   * @param polarizability The polarizability.
   * @param thole The Thole dampling constant.
   * @param polarizationGroup The atom types in the polarization group.
   */
  public PolarizeType(int atomType, double polarizability, double thole, int[] polarizationGroup) {
    super(POLARIZE, Integer.toString(atomType));
    this.type = atomType;
    this.thole = thole;
    this.polarizability = polarizability;
    this.polarizationGroup = polarizationGroup;
    if (thole == 0.0) {
      pdamp = 0.0;
    } else {
      pdamp = pow(polarizability, sixth);
    }
  }

  /**
   * assignPolarizationGroups.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param ip11 an array of {@link int} objects.
   * @param ip12 an array of {@link int} objects.
   * @param ip13 an array of {@link int} objects.
   */
  public static void assignPolarizationGroups(
      Atom[] atoms, int[][] ip11, int[][] ip12, int[][] ip13) {

    // Find directly connected group members for each atom.
    List<Integer> group = new ArrayList<>();
    List<Integer> polarizationGroup = new ArrayList<>();
    for (Atom ai : atoms) {
      group.clear();
      polarizationGroup.clear();
      Integer index = ai.getIndex() - 1;
      group.add(index);
      polarizationGroup.add(ai.getType());
      PolarizeType polarizeType = ai.getPolarizeType();
      if (polarizeType != null) {
        if (polarizeType.polarizationGroup != null) {
          for (int i : polarizeType.polarizationGroup) {
            if (!polarizationGroup.contains(i)) {
              polarizationGroup.add(i);
            }
          }
          growGroup(polarizationGroup, group, ai);
          Collections.sort(group);
          ip11[index] = new int[group.size()];
          int j = 0;
          for (int k : group) {
            ip11[index][j++] = k;
          }
        } else {
          ip11[index] = new int[group.size()];
          int j = 0;
          for (int k : group) {
            ip11[index][j++] = k;
          }
        }
      } else {
        String message =
            "The polarize keyword was not found for atom "
                + (index + 1)
                + " with type "
                + ai.getType();
        logger.severe(message);
      }
    }

    // Find 1-2 group relationships.
    int nAtoms = atoms.length;
    int[] mask = new int[nAtoms];
    List<Integer> list = new ArrayList<>();
    List<Integer> keep = new ArrayList<>();
    for (int i = 0; i < nAtoms; i++) {
      mask[i] = -1;
    }
    for (int i = 0; i < nAtoms; i++) {
      list.clear();
      for (int j : ip11[i]) {
        list.add(j);
        mask[j] = i;
      }
      keep.clear();
      for (int j : list) {
        Atom aj = atoms[j];
        List<Bond> bonds = aj.getBonds();
        for (Bond b : bonds) {
          Atom ak = b.get1_2(aj);
          int k = ak.getIndex() - 1;
          if (mask[k] != i) {
            keep.add(k);
          }
        }
      }
      list.clear();
      for (int j : keep) {
        for (int k : ip11[j]) {
          list.add(k);
        }
      }
      Collections.sort(list);
      ip12[i] = new int[list.size()];
      int j = 0;
      for (int k : list) {
        ip12[i][j++] = k;
      }
    }

    // Find 1-3 group relationships.
    for (int i = 0; i < nAtoms; i++) {
      mask[i] = -1;
    }
    for (int i = 0; i < nAtoms; i++) {
      for (int j : ip11[i]) {
        mask[j] = i;
      }
      for (int j : ip12[i]) {
        mask[j] = i;
      }
      list.clear();
      for (int j : ip12[i]) {
        for (int k : ip12[j]) {
          if (mask[k] != i) {
            if (!list.contains(k)) {
              list.add(k);
            }
          }
        }
      }
      ip13[i] = new int[list.size()];
      Collections.sort(list);
      int j = 0;
      for (int k : list) {
        ip13[i][j++] = k;
      }
    }
  }

  /**
   * Average two PolarizeType instances. The atom types to include in the new polarizationGroup must
   * be supplied.
   *
   * @param polarizeType1 a {@link ffx.potential.parameters.PolarizeType} object.
   * @param polarizeType2 a {@link ffx.potential.parameters.PolarizeType} object.
   * @param atomType a int.
   * @param polarizationGroup an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.PolarizeType} object.
   */
  public static PolarizeType average(
      PolarizeType polarizeType1,
      PolarizeType polarizeType2,
      int atomType,
      int[] polarizationGroup) {
    if (polarizeType1 == null || polarizeType2 == null) {
      return null;
    }
    double thole = (polarizeType1.thole + polarizeType2.thole) / 2.0;
    double polarizability = (polarizeType1.polarizability + polarizeType2.polarizability) / 2.0;
    return new PolarizeType(atomType, polarizability, thole, polarizationGroup);
  }

  /**
   * Construct a PolarizeType from an input string.
   *
   * @param input The overall input String.
   * @param tokens The input String tokenized.
   * @return a PolarizeType instance.
   */
  public static PolarizeType parse(String input, String[] tokens) {
    if (tokens.length < 4) {
      logger.log(Level.WARNING, "Invalid POLARIZE type:\n{0}", input);
    } else {
      try {
        int atomType = parseInt(tokens[1]);
        double polarizability = parseDouble(tokens[2]);
        double thole = parseDouble(tokens[3]);
        int entries = tokens.length - 4;
        int[] polarizationGroup = null;
        if (entries > 0) {
          polarizationGroup = new int[entries];
          for (int i = 4; i < tokens.length; i++) {
            polarizationGroup[i - 4] = parseInt(tokens[i]);
          }
        }
        return new PolarizeType(atomType, polarizability, thole, polarizationGroup);
      } catch (NumberFormatException e) {
        String message = "Exception parsing POLARIZE type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * A recursive method that checks all atoms bonded to the seed atom for inclusion in the
   * polarization group. The method is called on each newly found group member.
   *
   * @param polarizationGroup Atom types that should be included in the group.
   * @param group XYZ indices of current group members.
   * @param seed The bonds of the seed atom are queried for inclusion in the group.
   */
  private static void growGroup(List<Integer> polarizationGroup, List<Integer> group, Atom seed) {
    List<Bond> bonds = seed.getBonds();
    for (Bond bi : bonds) {
      Atom aj = bi.get1_2(seed);
      int tj = aj.getType();
      boolean added = false;
      for (int g : polarizationGroup) {
        if (g == tj) {
          Integer index = aj.getIndex() - 1;
          if (!group.contains(index)) {
            group.add(index);
            added = true;
            break;
          }
        }
      }
      if (added) {
        PolarizeType polarizeType = aj.getPolarizeType();
        for (int i : polarizeType.polarizationGroup) {
          if (!polarizationGroup.contains(i)) {
            polarizationGroup.add(i);
          }
        }
        growGroup(polarizationGroup, group, aj);
      }
    }
  }

  /**
   * add
   *
   * @param key a int.
   */
  public void add(int key) {
    for (int i : polarizationGroup) {
      if (key == i) {
        return;
      }
    }
    int len = polarizationGroup.length;
    int[] newGroup = new int[len + 1];
    arraycopy(polarizationGroup, 0, newGroup, 0, len);
    newGroup[len] = key;
    polarizationGroup = newGroup;
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
    PolarizeType polarizeType = (PolarizeType) o;
    return polarizeType.type == this.type;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(type);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted polarization type.
   */
  @Override
  public String toString() {
    StringBuilder polarizeString =
        new StringBuilder(String.format("polarize  %5d  %8.5f %8.5f", type, polarizability, thole));
    if (polarizationGroup != null) {
      for (int a : polarizationGroup) {
        polarizeString.append(String.format("  %5d", a));
      }
    }
    return polarizeString.toString();
  }

  /**
   * incrementType
   *
   * @param increment a int.
   */
  void incrementType(int increment) {
    type += increment;
    setKey(Integer.toString(type));
    if (polarizationGroup != null) {
      for (int i = 0; i < polarizationGroup.length; i++) {
        polarizationGroup[i] += increment;
      }
    }
  }

  /**
   * Add mapped known types to the polarization group of a new patch.
   *
   * @param typeMap a lookup between new atom types and known atom types.
   * @return a boolean.
   */
  boolean patchTypes(HashMap<AtomType, AtomType> typeMap) {
    if (polarizationGroup == null) {
      return false;
    }

    // Append known mapped types.
    int len = polarizationGroup.length;
    int added = 0;
    for (AtomType newType : typeMap.keySet()) {
      for (int i = 1; i < len; i++) {
        if (polarizationGroup[i] == newType.type) {
          AtomType knownType = typeMap.get(newType);
          added++;
          polarizationGroup = Arrays.copyOf(polarizationGroup, len + added);
          polarizationGroup[len + added - 1] = knownType.type;
        }
      }
    }

    return added > 0;
  }
}
