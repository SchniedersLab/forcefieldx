// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
import ffx.utilities.FFXProperty;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.Comparator;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.ForceFieldType.ATOM;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;

/**
 * The AtomType class represents one molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "atom", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [2 integers, name, quoted string, integer, real and integer]
    Provides the values needed to define a single force field atom type.
    The first two integer modifiers denote the atom type and class numbers.
    If the type and class are identical, only a single integer value is required.
    The next modifier is a three-character atom name, followed by an 24-character or less atom description contained in single quotes.
    The next two modifiers are the atomic number and atomic mass.
    The final integer modifier is the "valence" of the atom, defined as the expected number of attached or bonded atoms.
    """)
public final class AtomType extends BaseType implements Comparator<String> {

  /**
   * A Logger for the AngleType class.
   */
  private static final Logger logger = Logger.getLogger(AtomType.class.getName());
  /**
   * Short name (ie CH3/CH2 etc).
   */
  public final String name;
  /**
   * Description of the atom's bonding environment.
   */
  public final String environment;
  /**
   * Atomic Number.
   */
  public final int atomicNumber;
  /**
   * Atomic weight. "An atomic weight (relative atomic weight) of an element from a specified source
   * is the ratio of the average atomicWeight per atom of the element to 1/12 of the atomicWeight of
   * an atom of 12C"
   */
  public final double atomicWeight;
  /**
   * Valence number for this type.
   */
  public final int valence;
  /**
   * Atom type.
   */
  public int type;
  /**
   * Atom class.
   */
  public int atomClass;

  /**
   * AtomType Constructor.
   *
   * @param type         int
   * @param atomClass    int
   * @param name         String
   * @param environment  String
   * @param atomicNumber int
   * @param atomicWeight double
   * @param valence      int
   */
  public AtomType(int type, int atomClass, String name, String environment, int atomicNumber,
                  double atomicWeight, int valence) {
    super(ATOM, Integer.toString(type));
    this.type = type;
    this.atomClass = atomClass;
    this.name = name;
    this.environment = environment;
    this.atomicNumber = atomicNumber;
    this.atomicWeight = atomicWeight;
    this.valence = valence;
  }

  /**
   * Construct an AtomType from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return an AtomType instance.
   */
  public static AtomType parse(String input, String[] tokens) {
    if (tokens.length < 7) {
      logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
    } else {
      try {
        int index = 1;
        // Atom Type
        int type = parseInt(tokens[index++]);
        // Atom Class
        int atomClass;
        // The following try/catch checks for one of the following two cases:
        //
        // NUMBER TYPE CLASS IDENTIFIER ... (example is OPLS-AA)
        // vs.
        // NUMBER TYPE IDENTIFIER ... (example is OPLS-UA)
        //
        // If there is no atom class, the exception will be caught
        // and the atomClass field will remain equal to null.
        try {
          atomClass = parseInt(tokens[index]);
          // If the parseInt succeeds, this force field has atom classes.
          index++;
        } catch (NumberFormatException e) {
          // Some force fields do not use atom classes.
          atomClass = -1;
        }
        // Name
        String name = tokens[index].intern();
        // The "environment" string may contain spaces,
        // and is therefore surrounded in quotes located at "first" and
        // "last".
        int first = input.indexOf("\"");
        int last = input.lastIndexOf("\"");
        if (first >= last) {
          logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
          return null;
        }
        // Environment
        String environment = input.substring(first, last + 1).intern();
        // Shrink the tokens array to only include entries
        // after the environment field.
        tokens = input.substring(last + 1).trim().split(" +");
        index = 0;
        // Atomic Number
        int atomicNumber = parseInt(tokens[index++]);
        // Atomic Mass
        double mass = parseDouble(tokens[index++]);
        // Hybridization
        int hybridization = parseInt(tokens[index]);

        AtomType atomType = new AtomType(type, atomClass, name, environment, atomicNumber, mass,
            hybridization);
        if (!checkAtomicNumberAndMass(atomicNumber, mass)) {
          // Ignore united atom (UA) entries.
          if (!environment.toUpperCase().contains("UA")) {
            logger.warning(" Atomic number and weight do not agree:\n" + atomType);
          }
        }
        return atomType;
      } catch (NumberFormatException e) {
        String message = "Exception parsing AtomType:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
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
    AtomType atomType = (AtomType) o;
    return atomType.type == this.type;
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
   * <p>Nicely formatted atom type string.
   */
  @Override
  public String toString() {
    String s;
    if (atomClass >= 0) {
      s = format("atom  %5d  %5d  %-4s  %-25s  %3d  %8.4f  %d", type, atomClass, name, environment,
          atomicNumber, atomicWeight, valence);
    } else {
      s = format("atom  %5d  %-4s  %-25s  %3d  %8.4f  %d", type, name, environment, atomicNumber,
          atomicWeight, valence);
    }
    return s;
  }

  /**
   * Create an AtomType Element.
   *
   * @param doc        the Document instance.
   * @param forceField the ForceField to grab constants from.
   * @return the AtomType element.
   */
  public static Element getXMLAtomTypes(Document doc, ForceField forceField) {
    Element node = doc.createElement("AtomTypes");
    Map<String, AtomType> types = forceField.getAtomTypes();
    for (AtomType atomType : types.values()) {
      node.appendChild(atomType.toXML(doc));
    }
    return node;
  }

  /**
   * Write AtomType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Type");
    node.setAttribute("name", format("%d", type));
    node.setAttribute("class", format("%d", atomClass));
    if (atomicNumber >= 1) {
      node.setAttribute("element", format("%s", Atom.ElementSymbol.values()[atomicNumber - 1]));
    } else {
      // Handle force fields with dummy atoms that use atomic number 0.
      node.setAttribute("element", "");
    }
    node.setAttribute("mass", format("%.3f", atomicWeight));
    return node;
  }

  /**
   * incrementClassAndType
   *
   * @param classIncrement The value to increment the atom class by.
   * @param typeIncrement  The value to increment the atom type by.
   */
  void incrementClassAndType(int classIncrement, int typeIncrement) {
    atomClass += classIncrement;
    type += typeIncrement;
    setKey(Integer.toString(type));
  }

  /**
   * Check if the supplied atomic mass is within 0.1 AMU of the IUPAC value for the given atomic
   * number.
   * <p>
   * For atomic numbers outside the range 1 to 118, true always is returned.
   *
   * @param atomicNumber The atomic number.
   * @param mass         The atomic mass.
   * @return True if the given mass is within the given tolerance of the IUPAC value for the atomic
   * number.
   */
  public static boolean checkAtomicNumberAndMass(int atomicNumber, double mass) {
    return checkAtomicNumberAndMass(atomicNumber, mass, 0.1);
  }

  /**
   * Check if the supplied atomic mass is within the supplied tolerance (in AMU) of the IUPAC value
   * for the given atomic number.
   * <p>
   * For atomic numbers outside the range 1 to 118, true always is returned.
   *
   * @param atomicNumber The atomic number.
   * @param mass         The atomic mass.
   * @param tolerance    The error tolerance in AMU.
   * @return True if the given mass is within the given tolerance of the IUPAC value for the atomic
   * number.
   */
  public static boolean checkAtomicNumberAndMass(int atomicNumber, double mass, double tolerance) {
    // Ignore atomic numbers outside the range 1 to 118.
    if (atomicNumber == 0 || atomicNumber >= atomicMass.length) {
      return true;
    }

    double expected = atomicMass[atomicNumber - 1];
    return abs(expected - mass) < tolerance;
  }

  /**
   * <a href="https://iupac.qmul.ac.uk/AtWt">IUPAC Commission</a> on Isotopic Abundances and Atomic
   * Weights.
   * Retrieved on 1/24/22.
   */
  public static final double[] atomicMass = { /* H Hydrogen */  1.008,
      /* 2 He Helium */ 4.002,
      /* 3 Li Lithium */ 6.94,
      /* 4 Be Beryllium */ 9.012,
      /* 5 B Boron */ 10.81,
      /* 6 C Carbon */ 12.011,
      /* 7 N Nitrogen */ 14.007,
      /* 8 O Oxygen */ 15.999,
      /* 9 F Fluorine */ 18.998,
      /* 10 Ne Neon */ 20.1797,
      /* 11 Na Sodium */ 22.989,
      /* 12 Mg Magnesium */ 24.305,
      /* 13 Al Aluminium */ 26.981,
      /* 14 Si Silicon */ 28.085,
      /* 15 P Phosphorus */ 30.973,
      /* 16 S Sulfur */ 32.06,
      /* 17 Cl Chlorine */ 35.45,
      /* 18 Ar Argon */ 39.948,
      /* 19 K Potassium */ 39.0983,
      /* 20 Ca Calcium */ 40.078,
      /* 21 Sc Scandium */ 44.955,
      /* 22 Ti Titanium */ 47.867,
      /* 23 V Vanadium */ 50.9415,
      /* 24 Cr Chromium */ 51.9961,
      /* 25 Mn Manganese */ 54.938,
      /* 26 Fe Iron */ 55.845,
      /* 27 Co Cobalt */ 58.933,
      /* 28 Ni Nickel */ 58.6934,
      /* 29 Cu Copper */ 63.546,
      /* 30 Zn Zinc */ 65.38,
      /* 31 Ga Gallium */ 69.723,
      /* 32 Ge Germanium */ 72.630,
      /* 33 As Arsenic */ 74.921,
      /* 34 Se Selenium */ 78.971,
      /* 35 Br Bromine */ 79.904,
      /* 36 Kr Krypton */ 83.798,
      /* 37 Rb Rubidium */ 85.4678,
      /* 38 Sr Strontium */ 87.62,
      /* 39 Y Yttrium */ 88.905,
      /* 40 Zr Zirconium */ 91.224,
      /* 41 Nb Niobium */ 92.906,
      /* 42 Mo Molybdenum */ 95.95,
      /* 43 Tc Technetium */ 97.0,
      /* 44 Ru Ruthenium */ 101.07,
      /* 45 Rh Rhodium */ 102.905,
      /* 46 Pd Palladium */ 106.42,
      /* 47 Ag Silver */ 107.8682,
      /* 48 Cd Cadmium */ 112.414,
      /* 49 In Indium */ 114.818,
      /* 50 Sn Tin */ 118.710,
      /* 51 Sb Antimony */ 121.760,
      /* 52 Te Tellurium */ 127.60,
      /* 53 I Iodine */ 126.904,
      /* 54 Xe Xenon */ 131.293,
      /* 55 Cs Caesium */ 132.905,
      /* 56 Ba Barium */ 137.327,
      /* 57 La Lanthanum */ 138.905,
      /* 58 Ce Cerium */ 140.116,
      /* 59 Pr Praseodymium */ 140.907,
      /* 60 Nd Neodymium */ 144.242,
      /* 61 Pm Promethium */ 145.0,
      /* 62 Sm Samarium */ 150.36,
      /* 63 Eu Europium */ 151.964,
      /* 64 Gd Gadolinium */ 157.25,
      /* 65 Tb Terbium */ 158.925,
      /* 66 Dy Dysprosium */ 162.500,
      /* 67 Ho Holmium */ 164.930,
      /* 68 Er Erbium */ 167.259,
      /* 69 Tm Thulium */ 168.934,
      /* 70 Yb Ytterbium */ 173.045,
      /* 71 Lu Lutetium */ 174.9668,
      /* 72 Hf Hafnium */ 178.486,
      /* 73 Ta Tantalum */ 180.947,
      /* 74 W Tungsten */ 183.84,
      /* 75 Re Rhenium */ 186.207,
      /* 76 Os Osmium */ 190.23,
      /* 77 Ir Iridium */ 192.217,
      /* 78 Pt Platinum */ 195.084,
      /* 79 Au Gold */ 196.966,
      /* 80 Hg Mercury */ 200.592,
      /* 81 Tl Thallium */ 204.38,
      /* 82 Pb Lead */ 207.2,
      /* 83 Bi Bismuth */ 208.980,
      /* 84 Po Polonium */ 209.0,
      /* 85 At Astatine */ 210.0,
      /* 86 Rn Radon */ 222.0,
      /* 87 Fr Francium */ 223.0,
      /* 88 Ra Radium */ 226.0,
      /* 89 Ac Actinium */ 227.0,
      /* 90 Th Thorium */ 232.0377,
      /* 91 Pa Protactinium */ 231.035,
      /* 92 U Uranium */ 238.028,
      /* 93 Np Neptunium */ 237.0,
      /* 94 Pu Plutonium */ 244.0,
      /* 95 Am Americium */ 243.0,
      /* 96 Cm Curium */ 247.0,
      /* 97 Bk Berkelium */ 247.0,
      /* 98 Cf Californium */ 251.0,
      /* 99 Es Einsteinium */ 252.0,
      /* 100 Fm Fermium */ 257.0,
      /* 101 Md Mendelevium */ 258.0,
      /* 102 No Nobelium */ 259.0,
      /* 103 Lr Lawrencium */ 262.0,
      /* 104 Rf Rutherfordium */ 267.0,
      /* 105 Db Dubnium */ 270.0,
      /* 106 Sg Seaborgium */ 269.0,
      /* 107 Bh Bohrium */ 270.0,
      /* 108 Hs Hassium */ 270.0,
      /* 109 Mt Meitnerium */ 278.0,
      /* 110 Ds Darmstadtium */ 281.0,
      /* 111 Rg Roentgenium */ 281.0,
      /* 112 Cn Copernicium */ 285.0,
      /* 113 Nh Nihonium */ 286.0,
      /* 114 Fl Flerovium */ 289.0,
      /* 115 Mc Moscovium */ 289.0,
      /* 116 Lv Livermorium */ 293.0,
      /* 117 Ts Tennessine */ 293.0,
      /* 118 Og Oganesson */ 294.0};
}
