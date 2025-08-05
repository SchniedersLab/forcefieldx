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
package ffx.potential.parameters;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.utilities.FFXProperty;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.ForceFieldType.POLARIZE;
import static ffx.utilities.Constants.ANG_TO_NM;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.String.valueOf;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The PolarizeType class defines an isotropic atomic polarizability.
 * <p>
 * If this modifier is not present, then charge penetration values will be used for polarization damping, as in the HIPPO polarization model.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "polarize", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """ 
    [1 integer, up to 3 reals and up to 8 integers]
    Provides the values for a single atomic dipole polarizability parameter.
    The initial integer modifier, if positive, gives the atom type number for which a polarizability parameter is to be defined.
    If the first integer modifier is negative, then the parameter value to follow applies only to the specific atom whose atom number is the negative of the modifier.
    The first real number modifier gives the value of the dipole polarizability in Ang^3.
    The second real number modifier, if present, gives the Thole damping value.
    A Thole value of zero implies undamped polarization.
    The third real modifier, if present, gives a direct field damping value only used with the AMOEBA+ polarization model.
    The remaining integer modifiers list the atom type numbers of atoms directly bonded to the current atom and which will be considered to be part of the current atom’s polarization group.
    If the parameter is for a specific atom, then the integers defining the polarization group are ignored.
    """)
public final class PolarizeType extends BaseType implements Comparator<String> {

  private static final Logger logger = Logger.getLogger(PolarizeType.class.getName());

  private static final double sixth = 1.0 / 6.0;
  public static final double DEFAULT_DIRECT_11_SCALE = 0.0;
  public static final double DEFAULT_DIRECT_12_SCALE = 1.0;
  public static final double DEFAULT_DIRECT_13_SCALE = 1.0;
  public static final double DEFAULT_DIRECT_14_SCALE = 1.0;
  public static final double DEFAULT_POLAR_12_SCALE = 0.0;
  public static final double DEFAULT_POLAR_13_SCALE = 0.0;
  public static final double DEFAULT_POLAR_14_SCALE = 1.0;
  public static final double DEFAULT_POLAR_15_SCALE = 1.0;
  public static final double DEFAULT_POLAR_12_INTRA = 0.0;
  public static final double DEFAULT_POLAR_13_INTRA = 0.0;
  public static final double DEFAULT_POLAR_14_INTRA = 0.5;
  public static final double DEFAULT_POLAR_15_INTRA = 1.0;
  /**
   * Thole damping factor.
   */
  public final double thole;
  /**
   * Value of polarizability scale factor.
   */
  public double pdamp;
  /**
   * Direct polarization damping.
   */
  public final double ddp;
  /**
   * Isotropic polarizability in units of Angstroms^3.
   */
  public final double polarizability;
  /**
   * Atom type number.
   */
  public int type;
  /**
   * Connected types in the polarization group of each atom (can be null).
   */
  public int[] polarizationGroup;

  /**
   * PolarizeType Constructor.
   *
   * @param atomType          The atom type.
   * @param polarizability    The polarizability.
   * @param thole             The Thole damping constant.
   * @param polarizationGroup The atom types in the polarization group.
   */
  public PolarizeType(int atomType, double polarizability, double thole, double ddp,
                      int[] polarizationGroup) {
    super(POLARIZE, Integer.toString(atomType));
    this.type = atomType;
    this.thole = thole;
    this.polarizability = polarizability;
    this.ddp = 0.0;
    this.polarizationGroup = polarizationGroup;
    if (thole == 0.0) {
      pdamp = 0.0;
    } else {
      pdamp = pow(polarizability, sixth);
    }
  }

  /**
   * Construct a PolarizeType from a reference type and updated polarizability.
   *
   * @param polarizeType   The reference PolarizeType.
   * @param polarizability The updated polarizability.
   */
  public PolarizeType(PolarizeType polarizeType, double polarizability) {
    this(polarizeType.type, polarizability, polarizeType.thole, polarizeType.ddp,
        Arrays.copyOf(polarizeType.polarizationGroup, polarizeType.polarizationGroup.length));
    // Keep pdamp fixed during titration.
    pdamp = polarizeType.pdamp;
  }

  /**
   * Assign polarization groups to atoms based on their connectivity. This method will
   * assign 1-1, 1-2, and 1-3 polarization groups to the provided atoms.
   * The first dimension of the 2D arrays should already be allocated to the number of atoms.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param ip11  the 1-1 polarization groups as 2D int array.
   * @param ip12  the 1-2 polarization groups as 2D int array.
   * @param ip13  the 1-3 polarization groups as 2D int array.
   */
  public static void assignPolarizationGroups(Atom[] atoms, int[][] ip11, int[][] ip12,
                                              int[][] ip13) {

    // Find directly connected group members for each atom.
    List<Integer> group = new ArrayList<>();
    List<Integer> polarizationGroup = new ArrayList<>();
    for (Atom ai : atoms) {
      group.clear();
      polarizationGroup.clear();
      int index = ai.getIndex() - 1;
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
        String message = "The polarize keyword was not found for atom " + (index + 1) + " with type "
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
   * @param polarizeType1     The first PolarizeType.
   * @param polarizeType2     The second PolarizeType.
   * @param atomType          The atom type to use for the new PolarizeType.
   * @param polarizationGroup The atom types to include in the new polarizationGroup.
   * @return The averaged PolarizeType.
   */
  public static PolarizeType average(PolarizeType polarizeType1, PolarizeType polarizeType2,
                                     int atomType, int[] polarizationGroup) {
    if (polarizeType1 == null || polarizeType2 == null) {
      return null;
    }
    double thole = (polarizeType1.thole + polarizeType2.thole) / 2.0;
    double polarizability = (polarizeType1.polarizability + polarizeType2.polarizability) / 2.0;
    double ddp = (polarizeType1.ddp + polarizeType2.ddp) / 2.0;
    return new PolarizeType(atomType, polarizability, thole, ddp, polarizationGroup);
  }

  /**
   * Construct a PolarizeType from an input string.
   *
   * @param input  The overall input String.
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

        // There might be a direct polarization damping value.
        double ddp = 0.0;
        int nextToken = 4;
        try {
          // Try to parse the next token as an Int;
          // If there is an exception then it's a direct damping value.
          parseInt(tokens[nextToken]);
        } catch (NumberFormatException e) {
          ddp = parseDouble(tokens[nextToken]);
          nextToken = 5;
        } catch (ArrayIndexOutOfBoundsException e) {
          // No direct damping value.
        }
        int entries = tokens.length - nextToken;
        int[] polarizationGroup = null;
        if (entries > 0) {
          polarizationGroup = new int[entries];
          for (int i = nextToken; i < tokens.length; i++) {
            polarizationGroup[i - nextToken] = parseInt(tokens[i]);
          }
        }
        return new PolarizeType(atomType, polarizability, thole, ddp, polarizationGroup);
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
   * @param group             XYZ indices of current group members.
   * @param seed              The bonds of the seed atom are queried for inclusion in the group.
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
   * A recursive method that checks all atoms bonded to the seed atom for inclusion in the
   * polarization group. The method is called on each newly found group member.
   *
   * @param group XYZ indices of current group members.
   * @param seed  The bonds of the seed atom are queried for inclusion in the group.
   */
  public static void growGroup(List<Integer> group, Atom seed) {
    List<Bond> bonds = seed.getBonds();
    for (Bond bond : bonds) {
      Atom atom = bond.get1_2(seed);
      int tj = atom.getType();
      boolean added = false;
      PolarizeType polarizeType = seed.getPolarizeType();
      for (int type : polarizeType.polarizationGroup) {
        if (type == tj) {
          Integer index = atom.getIndex() - 1;
          if (!group.contains(index)) {
            group.add(index);
            added = true;
            break;
          }
        }
      }
      if (added) {
        growGroup(group, atom);
      }
    }
  }

  /**
   * Add an atom type to the polarization group.
   *
   * @param atomType The atom type to add.
   */
  public void add(int atomType) {
    for (int i : polarizationGroup) {
      if (atomType == i) {
        return;
      }
    }
    int len = polarizationGroup.length;
    int[] newGroup = new int[len + 1];
    arraycopy(polarizationGroup, 0, newGroup, 0, len);
    newGroup[len] = atomType;
    polarizationGroup = newGroup;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compare(String s1, String s2) {
    int t1 = parseInt(s1);
    int t2 = parseInt(s2);
    return Integer.compare(t1, t2);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    PolarizeType polarizeType = (PolarizeType) o;
    return polarizeType.type == this.type;
  }

  /**
   * {@inheritDoc}
   */
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
    StringBuilder polarizeString = new StringBuilder(
        format("polarize  %5d  %8.5f %8.5f", type, polarizability, thole));
    if (ddp != 0.0) {
      polarizeString.append(format(" %8.5f", ddp));
    }
    if (polarizationGroup != null) {
      for (int a : polarizationGroup) {
        polarizeString.append(format("  %5d", a));
      }
    }
    return polarizeString.toString();
  }

  /**
   * Add constant attributes to the AmoebaMultipoleForce
   *
   * @param node       the AmoebaMultipoleForce element.
   * @param forceField the ForceField for collecting constants.
   */
  public static void addXMLAttributes(Element node, ForceField forceField) {
    node.setAttribute("direct11Scale", valueOf(forceField.getDouble("direct-11-scale", DEFAULT_DIRECT_11_SCALE)));
    node.setAttribute("direct12Scale", valueOf(forceField.getDouble("direct-12-scale", DEFAULT_DIRECT_12_SCALE)));
    node.setAttribute("direct13Scale", valueOf(forceField.getDouble("direct-13-scale", DEFAULT_DIRECT_13_SCALE)));
    node.setAttribute("direct14Scale", valueOf(forceField.getDouble("direct-14-scale", DEFAULT_DIRECT_14_SCALE)));
    node.setAttribute("mutual11Scale", valueOf(forceField.getDouble("mutual-11-scale", 1.0)));
    node.setAttribute("mutual12Scale", valueOf(forceField.getDouble("mutual-12-scale", 1.0)));
    node.setAttribute("mutual13Scale", valueOf(forceField.getDouble("mutual-13-scale", 1.0)));
    node.setAttribute("mutual14Scale", valueOf(forceField.getDouble("mutual-14-scale", 1.0)));
    node.setAttribute("polar12Scale", valueOf(forceField.getDouble("polar-12-scale", DEFAULT_POLAR_12_SCALE)));
    node.setAttribute("polar13Scale", valueOf(forceField.getDouble("polar-13-scale", DEFAULT_POLAR_13_SCALE)));
    node.setAttribute("polar14Scale", valueOf(forceField.getDouble("polar-14-scale", DEFAULT_POLAR_14_SCALE)));
    node.setAttribute("polar15Scale", valueOf(forceField.getDouble("polar-15-scale", DEFAULT_POLAR_15_SCALE)));
    // polar-12-intra and polar-13-intra are zero by default. Other values are not supported by FFX (or OpenMM).
    // polar-15-intra is 1.0 by default. Other values are not supported by FFX or OpenMM.
    node.setAttribute("polar14Intra", valueOf(forceField.getDouble("polar-14-intra", DEFAULT_POLAR_14_INTRA)));
  }

  /**
   * Write PolarizeType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Polarize");
    node.setAttribute("type", format("%d", type));
    // Convert Ang^3 to nm^3
    node.setAttribute("polarizability", format("%f", polarizability * ANG_TO_NM * ANG_TO_NM * ANG_TO_NM));
    node.setAttribute("thole", format("%f", thole));
    // TODO: OMM does not have support for AMOEBA+ and its use of direct polarization damping (ddp)
    int i = 1;
    if (polarizationGroup != null) {
      for (int a : polarizationGroup) {
        node.setAttribute(format("pgrp%d", i), format("%d", a));
        i++;
      }
    }
    return node;
  }

  /**
   * Both the type and polarization group are incremented by the same amount.
   *
   * @param increment The increment to add to the type.
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
