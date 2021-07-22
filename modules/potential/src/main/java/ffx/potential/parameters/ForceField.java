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

import static java.lang.String.format;

import java.net.URL;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * The ForceField class organizes parameters for a molecular mechanics force field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ForceField {

  private static final Logger logger = Logger.getLogger(ForceField.class.getName());
  /** A map between a force field name and its internal parameter file. */
  private static final Map<ForceFieldName, URL> forceFields = new EnumMap<>(ForceFieldName.class);

  static {
    ClassLoader cl = ForceField.class.getClassLoader();
    String prefix = "ffx/potential/parameters/ff/";
    for (ForceFieldName ff : ForceFieldName.values()) {
      forceFields.put(ff, cl.getResource(prefix + ff));
    }
  }

  /** Flag to prevent patch renumbering. */
  private final boolean noRenumbering;
  /** The CompositeConfiguration that contains key=value property pairs from a number of sources. */
  private final CompositeConfiguration properties;

  private final Map<String, AngleType> angleTypes;
  private final Map<String, AngleType> anglepTypes;
  private final Map<String, AtomType> atomTypes;
  private final Map<String, BioType> bioTypes;
  private final Map<String, BondType> bondTypes;
  private final Map<String, ChargeType> chargeTypes;
  private final Map<String, MultipoleType> multipoleTypes;
  private final Map<String, OutOfPlaneBendType> outOfPlaneBendTypes;
  private final Map<String, PolarizeType> polarizeTypes;
  private final Map<String, StretchBendType> stretchBendTypes;
  private final Map<String, StretchTorsionType> stretchTorsionTypes;
  private final Map<String, AngleTorsionType> angleTorsionTypes;
  private final Map<String, PiTorsionType> piTorsionTypes;
  private final Map<String, TorsionType> torsionTypes;
  private final Map<String, TorsionType> improperTypes;
  private final Map<String, ImproperTorsionType> imptorsTypes;
  private final Map<String, SoluteType> soluteTypes;
  private final Map<String, TorsionTorsionType> torsionTorsionTypes;
  private final Map<String, UreyBradleyType> ureyBradleyTypes;
  private final Map<String, VDWType> vanderWaalsTypes;
  private final Map<String, VDWType> vanderWaals14Types;
  private final Map<String, RelativeSolvationType> relativeSolvationTypes;
  private final Map<ForceFieldType, Map<String, ? extends BaseType>> forceFieldTypes;
  /** URL to the force field parameter file. */
  public URL forceFieldURL;
  /**
   * ForceField Constructor.
   *
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public ForceField(CompositeConfiguration properties) {
    this.properties = properties;

    noRenumbering = properties.getBoolean("noPatchRenumbering", false);

    /*
     Each force field "type" implements the "Comparator<String>" interface
     so that passing an "empty" instance of the "type" to its TreeMap
     constructor will keep the types sorted.
    */
    angleTypes = new TreeMap<>(new AngleType(new int[3], 0, new double[1], null));
    anglepTypes = new TreeMap<>(new AngleType(new int[3], 0, new double[1], null, null));
    atomTypes = new TreeMap<>(new AtomType(0, 0, null, null, 0, 0, 0));
    bioTypes = new TreeMap<>(new BioType(0, null, null, 0, null));
    bondTypes = new TreeMap<>(new BondType(new int[2], 0, 0, null));
    chargeTypes = new TreeMap<>(new ChargeType(0, 0));
    soluteTypes = new TreeMap<>(new SoluteType(0, 0.0, 0.0, 0.0));
    multipoleTypes = new TreeMap<>(new MultipoleType(new double[10], null, null, false));
    outOfPlaneBendTypes = new TreeMap<>(new OutOfPlaneBendType(new int[4], 0));
    piTorsionTypes = new TreeMap<>(new PiTorsionType(new int[2], 0));
    polarizeTypes = new TreeMap<>(new PolarizeType(0, 0, 0, new int[1]));
    stretchBendTypes = new TreeMap<>(new StretchBendType(new int[3], new double[1]));
    stretchTorsionTypes = new TreeMap<>(new StretchTorsionType(new int[4], new double[1]));
    angleTorsionTypes = new TreeMap<>(new AngleTorsionType(new int[4], new double[1]));
    torsionTorsionTypes = new TreeMap<>();
    torsionTypes =
        new TreeMap<>(new TorsionType(new int[4], new double[1], new double[1], new int[1]));
    improperTypes =
        new TreeMap<>(new TorsionType(new int[4], new double[1], new double[1], new int[1]));
    imptorsTypes = new TreeMap<>(new ImproperTorsionType(new int[4], 0.0, 0.0, 2));
    ureyBradleyTypes = new TreeMap<>(new UreyBradleyType(new int[3], 0, 0));
    vanderWaalsTypes = new TreeMap<>(new VDWType(0, 0, 0, 0));
    vanderWaals14Types = new TreeMap<>(new VDWType(0, 0, 0, 0));
    relativeSolvationTypes = new TreeMap<>(new RelativeSolvationType("", 0.0));

    forceFieldTypes = new EnumMap<>(ForceFieldType.class);
    forceFieldTypes.put(ForceFieldType.ANGLE, angleTypes);
    forceFieldTypes.put(ForceFieldType.ANGLEP, anglepTypes);
    forceFieldTypes.put(ForceFieldType.ATOM, atomTypes);
    forceFieldTypes.put(ForceFieldType.BOND, bondTypes);
    forceFieldTypes.put(ForceFieldType.BIOTYPE, bioTypes);
    forceFieldTypes.put(ForceFieldType.CHARGE, chargeTypes);
    forceFieldTypes.put(ForceFieldType.SOLUTE, soluteTypes);
    forceFieldTypes.put(ForceFieldType.OPBEND, outOfPlaneBendTypes);
    forceFieldTypes.put(ForceFieldType.MULTIPOLE, multipoleTypes);
    forceFieldTypes.put(ForceFieldType.PITORS, piTorsionTypes);
    forceFieldTypes.put(ForceFieldType.POLARIZE, polarizeTypes);
    forceFieldTypes.put(ForceFieldType.STRBND, stretchBendTypes);
    forceFieldTypes.put(ForceFieldType.STRTORS, stretchTorsionTypes);
    forceFieldTypes.put(ForceFieldType.ANGTORS, angleTorsionTypes);
    forceFieldTypes.put(ForceFieldType.TORSION, torsionTypes);
    forceFieldTypes.put(ForceFieldType.IMPROPER, improperTypes);
    forceFieldTypes.put(ForceFieldType.IMPTORS, imptorsTypes);
    forceFieldTypes.put(ForceFieldType.TORTORS, torsionTorsionTypes);
    forceFieldTypes.put(ForceFieldType.UREYBRAD, ureyBradleyTypes);
    forceFieldTypes.put(ForceFieldType.VDW, vanderWaalsTypes);
    forceFieldTypes.put(ForceFieldType.VDW14, vanderWaals14Types);
    forceFieldTypes.put(ForceFieldType.RELATIVESOLV, relativeSolvationTypes);

    trueImpliedBoolean("ELEC_LAMBDATERM", "GK_LAMBDATERM");
    trueImpliedBoolean("LAMBDATERM", "VDW_LAMBDATERM", "ELEC_LAMBDATERM", "GK_LAMBDATERM");
  }

  /**
   * Get for the URL for the named force field.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField.ForceFieldName} object.
   * @return a {@link java.net.URL} object.
   */
  public static URL getForceFieldURL(ForceFieldName forceField) {
    if (forceField == null) {
      return null;
    }
    return forceFields.get(forceField);
  }

  /**
   * isForceFieldKeyword.
   *
   * @param keyword a {@link java.lang.String} object.
   * @return a boolean.
   */
  public static boolean isForceFieldKeyword(String keyword) {
    keyword = toEnumForm(keyword);
    try {
      ForceFieldType.valueOf(keyword);
      return true;
    } catch (Exception e) {
      // Ignore.
    }
    return false;
  }

  /**
   * Enums are uppercase with underscores, but property files use lower case with dashes.
   *
   * @param key an input keyword
   * @return the keyword in Enum form.
   */
  public static String toEnumForm(String key) {
    if (key == null) {
      return null;
    }
    return key.toUpperCase().replace("-", "_");
  }

  /**
   * Enums are uppercase with underscores, but property files use lower case with dashes.
   *
   * @param s an input Enum string
   * @return the normalized keyword
   */
  public static String toPropertyForm(String s) {
    if (s == null) {
      return null;
    }
    return s.toLowerCase().replace("_", "-");
  }

  /**
   * Add an instance of a force field type. Force Field types are more complicated than simple
   * Strings or doubles, in that they have multiple fields and may occur multiple times.
   *
   * @param <T> ForceFieldType to add that extends BaseType
   * @param type The ForceFieldType to add.
   */
  @SuppressWarnings("unchecked")
  public <T extends BaseType> void addForceFieldType(T type) {
    if (type == null) {
      logger.info(" Null force field type ignored.");
      return;
    }

    Map<String, T> treeMap = (Map<String, T>) forceFieldTypes.get(type.forceFieldType);
    if (treeMap == null) {
      logger.log(Level.INFO, " Unrecognized force field type ignored {0}", type.forceFieldType);
      type.print();
      return;
    }
    if (treeMap.containsKey(type.key)) {
      if (treeMap.get(type.key).toString().equalsIgnoreCase(type.toString())) {
        // Ignore this type if its identical to an existing type.
        return;
      }
      logger.log(
          Level.WARNING,
          " A force field entry of type {0} already exists with the key: {1}\n The (discarded) old entry: {2}\n The new entry            : {3}",
          new Object[] {
            type.forceFieldType, type.key, treeMap.get(type.key).toString(), type.toString()
          });
    }
    treeMap.put(type.key, type);
  }

  /**
   * Add a property from a external parameter file.
   *
   * @param property Property string.
   * @param value double
   */
  public void addProperty(String property, String value) {
    if (property == null) {
      return;
    }
    String key = toPropertyForm(property);
    //        try {
    //            String old = getString(key);
    //            if (old.equalsIgnoreCase(value)) {
    //                return;
    //            }
    //            logger.info(format("  Existing %s  %s", key, old));
    //        } catch (Exception e) {
    //            // Property does not exist yet.
    //        } finally {
    //            properties.addProperty(key, value);
    //            logger.info(format("  Added    %s  %s", key, value));
    //        }
    properties.addProperty(key, value);
  }

  /**
   * Clear a property from the force field instance.
   * @param property Property to clear.
   */
  public void clearProperty(String property) {
    properties.clearProperty(property);
  }

  /**
   * Append a 2nd ForceField "patch" to the current ForceField. Note that only the force field types
   * are appended; properties are ignored.
   *
   * @param patch The force field patch to append.
   */
  public void append(ForceField patch) {
    // Determine the highest current atom class, atom type and biotype index.
    int classOffset = maxClass();
    int typeOffset = maxType();
    int bioTypeOffset = maxBioType();

    int minClass = patch.minClass();
    int minType = patch.minType();
    int minBioType = patch.minBioType();

    classOffset -= (minClass - 1);
    typeOffset -= (minType - 1);
    bioTypeOffset -= (minBioType - 1);

    patch.renumberForceField(classOffset, typeOffset, bioTypeOffset);

    for (AngleType angleType : patch.angleTypes.values()) {
      angleTypes.put(angleType.getKey(), angleType);
    }

    for (AngleType angleType : patch.anglepTypes.values()) {
      angleTypes.put(angleType.getKey(), angleType);
    }

    for (AtomType atomType : patch.atomTypes.values()) {
      atomTypes.put(atomType.getKey(), atomType);
    }

    for (BioType bioType : patch.bioTypes.values()) {
      bioTypes.put(bioType.getKey(), bioType);
    }

    for (BondType bondType : patch.bondTypes.values()) {
      bondTypes.put(bondType.getKey(), bondType);
    }

    for (MultipoleType multipoleType : patch.multipoleTypes.values()) {
      multipoleTypes.put(multipoleType.getKey(), multipoleType);
    }

    for (OutOfPlaneBendType outOfPlaneBendType : patch.outOfPlaneBendTypes.values()) {
      outOfPlaneBendTypes.put(outOfPlaneBendType.getKey(), outOfPlaneBendType);
    }

    for (PiTorsionType piTorsionType : patch.piTorsionTypes.values()) {
      piTorsionTypes.put(piTorsionType.getKey(), piTorsionType);
    }

    for (PolarizeType polarizeType : patch.polarizeTypes.values()) {
      polarizeTypes.put(polarizeType.getKey(), polarizeType);
    }

    for (StretchBendType stretchBendType : patch.stretchBendTypes.values()) {
      stretchBendTypes.put(stretchBendType.getKey(), stretchBendType);
    }

    for (StretchTorsionType stretchTorsionType : patch.stretchTorsionTypes.values()) {
      stretchTorsionTypes.put(stretchTorsionType.getKey(), stretchTorsionType);
    }

    for (AngleTorsionType angleTorsionType : patch.angleTorsionTypes.values()) {
      angleTorsionTypes.put(angleTorsionType.getKey(), angleTorsionType);
    }

    for (TorsionTorsionType torsionTorsionType : patch.torsionTorsionTypes.values()) {
      torsionTorsionTypes.put(torsionTorsionType.getKey(), torsionTorsionType);
    }

    for (TorsionType torsionType : patch.torsionTypes.values()) {
      torsionTypes.put(torsionType.getKey(), torsionType);
    }

    for (TorsionType torsionType : patch.improperTypes.values()) {
      torsionTypes.put(torsionType.getKey(), torsionType);
    }

    for (ImproperTorsionType improperTorsionType : patch.imptorsTypes.values()) {
      imptorsTypes.put(improperTorsionType.getKey(), improperTorsionType);
    }

    for (UreyBradleyType ureyBradleyType : patch.ureyBradleyTypes.values()) {
      ureyBradleyTypes.put(ureyBradleyType.getKey(), ureyBradleyType);
    }

    for (VDWType vdwType : patch.vanderWaalsTypes.values()) {
      vanderWaalsTypes.put(vdwType.getKey(), vdwType);
    }

    for (VDWType vdwType : patch.vanderWaals14Types.values()) {
      vanderWaals14Types.put(vdwType.getKey(), vdwType);
    }

    for (SoluteType soluteType : patch.soluteTypes.values()) {
      soluteTypes.put(soluteType.getKey(), soluteType);
    }

    for (RelativeSolvationType rsType : patch.relativeSolvationTypes.values()) {
      relativeSolvationTypes.put(rsType.getKey(), rsType);
    }

    // Is this a modified residue patch?
    String modres = patch.getString("MODRES", "false");
    if (!modres.equalsIgnoreCase("false")) {
      logger.info(" Adding modified residue patch.");
      modifiedResidue(modres);
    }
  }

  /**
   * getAngleTorsionType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.AngleTorsionType} object.
   */
  public AngleTorsionType getAngleTorsionType(String key) {
    return angleTorsionTypes.get(key);
  }

  /**
   * getAngleType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.AngleType} object.
   */
  public AngleType getAngleType(String key) {
    AngleType angleType = angleTypes.get(key);
    if (angleType == null) {
      angleType = anglepTypes.get(key);
    }
    return angleType;
  }

  /**
   * getAtomType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.AtomType} object.
   */
  public AtomType getAtomType(String key) {
    return atomTypes.get(key);
  }

  /**
   * getAtomType
   *
   * @param moleculeName a {@link java.lang.String} object.
   * @param atomName a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.AtomType} object.
   */
  public AtomType getAtomType(String moleculeName, String atomName) {
    for (BioType bioType : bioTypes.values()) {
      if (bioType.moleculeName.equalsIgnoreCase(moleculeName)
          && bioType.atomName.equalsIgnoreCase(atomName)) {
        String key = Integer.toString(bioType.atomType);
        return atomTypes.get(key);
      }
    }
    return null;
  }

  /**
   * Getter for the field <code>atomTypes</code>.
   *
   * @param moleculeName a {@link java.lang.String} object.
   * @return a {@link java.util.HashMap} object.
   */
  public HashMap<String, AtomType> getAtomTypes(String moleculeName) {
    HashMap<String, AtomType> types = new HashMap<>();
    for (BioType bioType : bioTypes.values()) {
      if (bioType.moleculeName.equalsIgnoreCase(moleculeName)) {
        String key = Integer.toString(bioType.atomType);
        AtomType type = atomTypes.get(key);
        types.put(bioType.atomName.toUpperCase(), type);
      }
    }
    return types;
  }

  /**
   * getBioType.
   *
   * @param moleculeName a {@link java.lang.String} object.
   * @param atomName a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.BioType} object.
   */
  public BioType getBioType(String moleculeName, String atomName) {
    for (BioType bioType : bioTypes.values()) {
      if (bioType.moleculeName.equalsIgnoreCase(moleculeName)
          && bioType.atomName.equalsIgnoreCase(atomName)) {
        return bioType;
      }
    }
    return null;
  }

  /**
   * getBioType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.BioType} object.
   */
  public BioType getBioType(String key) {
    return bioTypes.get(key);
  }

  /**
   * getBioTypeMap.
   *
   * @return a {@link java.util.Map} object.
   */
  public Map<String, BioType> getBioTypeMap() {
    return bioTypes;
  }

  /**
   * getBondType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.BondType} object.
   */
  public BondType getBondType(String key) {
    return bondTypes.get(key);
  }

  /**
   * getBonds
   *
   * @param moleculeName a {@link java.lang.String} object.
   * @param atomName a {@link java.lang.String} object.
   * @return an array of {@link java.lang.String} objects.
   */
  public String[] getBonds(String moleculeName, String atomName) {
    for (BioType bioType : bioTypes.values()) {
      if (bioType.moleculeName.equalsIgnoreCase(moleculeName)
          && bioType.atomName.equalsIgnoreCase(atomName)) {
        return bioType.bonds;
      }
    }
    return null;
  }

  /**
   * getBoolean
   *
   * @param property The property to return.
   * @return a boolean.
   * @throws java.lang.Exception if any.
   */
  public boolean getBoolean(String property) throws Exception {
    if (property == null) {
      throw new Exception("NULL property");
    }
    String key = toPropertyForm(property);
    if (!properties.containsKey(key)) {
      throw new Exception("Undefined property: " + key);
    }
    return properties.getBoolean(key);
  }

  /**
   * getBoolean
   *
   * @param property The property to return.
   * @param defaultBoolean The default to return.
   * @return a boolean.
   */
  public boolean getBoolean(String property, boolean defaultBoolean) {
    try {
      return getBoolean(property);
    } catch (Exception e) {
      return defaultBoolean;
    }
  }

  /**
   * getDouble
   *
   * @param property The property to return.
   * @return an double.
   * @throws java.lang.Exception if any.
   */
  public double getDouble(String property) throws Exception {
    if (property == null) {
      throw new Exception("NULL property");
    }
    String key = toPropertyForm(property);
    if (!properties.containsKey(key)) {
      throw new Exception("Undefined property: " + key);
    }
    return properties.getDouble(key);
  }

  /**
   * getDouble
   *
   * @param property The property to return.
   * @param defaultDouble The default to return.
   * @return an double.
   */
  public double getDouble(String property, Double defaultDouble) {
    try {
      return getDouble(property);
    } catch (Exception e) {
      return defaultDouble;
    }
  }

  /**
   * getForceFieldTypeCount
   *
   * @param type a {@link ffx.potential.parameters.ForceField.ForceFieldType} object.
   * @return a int.
   */
  @SuppressWarnings("unchecked")
  public int getForceFieldTypeCount(ForceFieldType type) {
    TreeMap<String, BaseType> treeMap = (TreeMap<String, BaseType>) forceFieldTypes.get(type);
    if (treeMap == null) {
      logger.log(Level.WARNING, "Unrecognized Force Field Type: {0}", type);
      return 0;
    }
    return treeMap.size();
  }

  /**
   * getImproperType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.TorsionType} object.
   */
  public ImproperTorsionType getImproperType(String key) {
    return imptorsTypes.get(key);
  }

  /**
   * getImproperType
   *
   * @return a {@link ffx.potential.parameters.TorsionType} object.
   */
  public Collection<ImproperTorsionType> getImproperTypes() {
    return imptorsTypes.values();
  }

  /**
   * getInteger
   *
   * @param property The property to return.
   * @return an int.
   * @throws java.lang.Exception if any.
   */
  public int getInteger(String property) throws Exception {
    if (property == null) {
      throw new Exception("NULL property");
    }
    String key = toPropertyForm(property);
    if (!properties.containsKey(key)) {
      throw new Exception("Undefined property: " + key);
    }
    return properties.getInt(key);
  }

  /**
   * getInteger
   *
   * @param property The property to return.
   * @param defaultInteger The default to return.
   * @return an int.
   */
  public int getInteger(String property, Integer defaultInteger) {
    try {
      return getInteger(property);
    } catch (Exception e) {
      return defaultInteger;
    }
  }

  /**
   * getMultipoleType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.MultipoleType} object.
   */
  public MultipoleType getMultipoleType(String key) {
    return multipoleTypes.get(key);
  }

  /**
   * Fine the MultipoleType whose key begins with the supplied String.
   * If there are more than one MultipoleTypes that begin with the key, null is returned.
   * @param key The key to search for.
   * @return The MultipoleType if one and only one match is found.
   */
  public MultipoleType getMultipoleTypeBeginsWith(String key) {
    int count = 0;
    MultipoleType multipoleType = null;
    for (String s : multipoleTypes.keySet()) {
      if (s.startsWith(key + " ")){
        count++;
        multipoleType = multipoleTypes.get(s);
      }
    }

    if (count == 1) {
      return multipoleType;
    }

    return null;
  }

  /**
   * getOutOfPlaneBendType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
   */
  public OutOfPlaneBendType getOutOfPlaneBendType(String key) {
    return outOfPlaneBendTypes.get(key);
  }

  /**
   * getPiTorsionType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.PiTorsionType} object.
   */
  public PiTorsionType getPiTorsionType(String key) {
    return piTorsionTypes.get(key);
  }

  /**
   * getPolarizeType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.PolarizeType} object.
   */
  public PolarizeType getPolarizeType(String key) {
    return polarizeTypes.get(key);
  }

  /**
   * Getter for the field <code>properties</code>.
   *
   * @return a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public CompositeConfiguration getProperties() {
    return properties;
  }

  /**
   * Getter for the field <code>relativeSolvationTypes</code>.
   *
   * @return a {@link java.util.HashMap} object.
   */
  public HashMap<String, RelativeSolvationType> getRelativeSolvationTypes() {
    HashMap<String, RelativeSolvationType> types = new HashMap<>();
    for (String key : relativeSolvationTypes.keySet()) {
      types.put(key, relativeSolvationTypes.get(key));
    }
    return types;
  }

  /**
   * Get a SoluteType.
   *
   * @param key The class of the atom.
   * @return The SoluteType if its present.
   */
  public SoluteType getSoluteType(String key) {
    return soluteTypes.get(key);
  }

  public Map<String, SoluteType> getSoluteTypes() {
    return soluteTypes;
  }

  /**
   * getStretchBendType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.StretchBendType} object.
   */
  public StretchBendType getStretchBendType(String key) {
    return stretchBendTypes.get(key);
  }

  /**
   * getStretchTorsionType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.StretchTorsionType} object.
   */
  public StretchTorsionType getStretchTorsionType(String key) {
    return stretchTorsionTypes.get(key);
  }

  /**
   * getBoolean
   *
   * @param property The property to return.
   * @return a String.
   * @throws java.lang.Exception if any.
   */
  public String getString(String property) throws Exception {
    if (property == null) {
      throw new Exception("NULL property");
    }
    String key = toPropertyForm(property);
    if (!properties.containsKey(key)) {
      throw new Exception("Undefined property: " + key);
    }
    return properties.getString(key);
  }

  /**
   * getBoolean
   *
   * @param property The property to return.
   * @param defaultString The default to return.
   * @return a boolean.
   */
  public String getString(String property, String defaultString) {
    try {
      return getString(property);
    } catch (Exception e) {
      return defaultString;
    }
  }

  /**
   * getTorsionTorsionType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.TorsionTorsionType} object.
   */
  public TorsionTorsionType getTorsionTorsionType(String key) {
    return torsionTorsionTypes.get(key);
  }

  /**
   * getTorsionType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.TorsionType} object.
   */
  public TorsionType getTorsionType(String key) {
    return torsionTypes.get(key);
  }

  /**
   * getUreyBradleyType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.UreyBradleyType} object.
   */
  public UreyBradleyType getUreyBradleyType(String key) {
    return ureyBradleyTypes.get(key);
  }

  /**
   * getVDW14Type
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.VDWType} object.
   */
  public VDWType getVDW14Type(String key) {
    return vanderWaals14Types.get(key);
  }

  /**
   * getVDW14Types
   *
   * @return a {@link java.util.Map} object.
   */
  public Map<String, VDWType> getVDW14Types() {
    return vanderWaals14Types;
  }

  /**
   * getVDWType
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link ffx.potential.parameters.VDWType} object.
   */
  public VDWType getVDWType(String key) {
    return vanderWaalsTypes.get(key);
  }

  /**
   * getVDWTypes
   *
   * @return a {@link java.util.Map} object.
   */
  public Map<String, VDWType> getVDWTypes() {
    return vanderWaalsTypes;
  }

  /**
   * Checks if a property was specified.
   *
   * @param property String to check.
   * @return Ever specified.
   * @throws java.lang.NullPointerException If forceFieldDouble is null.
   */
  public boolean hasProperty(String property) throws NullPointerException {
    if (property == null) {
      throw new NullPointerException("NULL keyword");
    }
    String key = toPropertyForm(property);
    return properties.containsKey(key);
  }

  /** log */
  public void log() {
    for (ForceFieldType s : forceFieldTypes.keySet()) {
      log(s.toString());
    }
  }

  /**
   * Prints any force field keyword to Standard.out.
   *
   * @param key String
   */
  public void log(String key) {
    ForceFieldType type = ForceFieldType.valueOf(key);
    logger.info(toString(type));
  }

  /** print */
  public void print() {
    for (ForceFieldType s : forceFieldTypes.keySet()) {
      print(s.toString());
    }
  }

  /**
   * print
   *
   * @param key a {@link java.lang.String} object.
   */
  public void print(String key) {
    ForceFieldType type = ForceFieldType.valueOf(key);
    logger.info(toString(type));
  }

  /**
   * Renumber ForceField class, type and biotype values.
   *
   * @param classOffset a int.
   * @param typeOffset a int.
   * @param bioTypeOffset a int.
   */
  public void renumberForceField(int classOffset, int typeOffset, int bioTypeOffset) {
    if (noRenumbering) {
      return;
    }
    for (AngleType angleType : angleTypes.values()) {
      angleType.incrementClasses(classOffset);
    }

    for (AngleType angleType : anglepTypes.values()) {
      angleType.incrementClasses(classOffset);
    }

    for (AtomType atomType : atomTypes.values()) {
      atomType.incrementClassAndType(classOffset, typeOffset);
    }

    for (BioType bioType : bioTypes.values()) {
      bioType.incrementIndexAndType(bioTypeOffset, typeOffset);
    }

    for (BondType bondType : bondTypes.values()) {
      bondType.incrementClasses(classOffset);
    }

    for (MultipoleType multipoleType : multipoleTypes.values()) {
      multipoleType.incrementType(typeOffset);
    }

    for (OutOfPlaneBendType outOfPlaneBendType : outOfPlaneBendTypes.values()) {
      outOfPlaneBendType.incrementClasses(classOffset);
    }

    for (PiTorsionType piTorsionType : piTorsionTypes.values()) {
      piTorsionType.incrementClasses(classOffset);
    }

    for (PolarizeType polarizeType : polarizeTypes.values()) {
      polarizeType.incrementType(typeOffset);
    }

    for (StretchBendType stretchBendType : stretchBendTypes.values()) {
      stretchBendType.incrementClasses(classOffset);
    }

    for (StretchTorsionType stretchTorsionType : stretchTorsionTypes.values()) {
      stretchTorsionType.incrementClasses(classOffset);
    }

    for (AngleTorsionType angleTorsionType : angleTorsionTypes.values()) {
      angleTorsionType.incrementClasses(classOffset);
    }

    for (TorsionTorsionType torsionTorsionType : torsionTorsionTypes.values()) {
      torsionTorsionType.incrementClasses(classOffset);
    }

    for (TorsionType torsionType : torsionTypes.values()) {
      torsionType.incrementClasses(classOffset);
    }

    for (TorsionType torsionType : improperTypes.values()) {
      torsionType.incrementClasses(classOffset);
    }

    for (ImproperTorsionType improperTorsionType : imptorsTypes.values()) {
      improperTorsionType.incrementClasses(classOffset);
    }

    for (UreyBradleyType ureyBradleyType : ureyBradleyTypes.values()) {
      ureyBradleyType.incrementClasses(classOffset);
    }

    for (VDWType vanderWaalsType : vanderWaalsTypes.values()) {
      vanderWaalsType.incrementClass(classOffset);
    }

    for (VDWType vanderWaals14Type : vanderWaals14Types.values()) {
      vanderWaals14Type.incrementClass(classOffset);
    }

    for (SoluteType soluteType : soluteTypes.values()) {
      soluteType.incrementType(typeOffset);
    }
  }

  /**
   * The AngleFunction in use by this ForceField.
   *
   * @param angleFunction The AngleFunction to use.
   */
  public void setAngleFunction(AngleType.AngleFunction angleFunction) {
    for (AngleType angleType : anglepTypes.values()) {
      angleType.setAngleFunction(angleFunction);
    }
    for (AngleType angleType : angleTypes.values()) {
      angleType.setAngleFunction(angleFunction);
    }
  }

  /**
   * The BondFunction in use by this ForceField.
   *
   * @param bondFunction The BondFunction to use.
   */
  public void setBondFunction(BondType.BondFunction bondFunction) {
    for (BondType bondType : bondTypes.values()) {
      bondType.setBondFunction(bondFunction);
    }
  }

  /**
   * setTorsionScale.
   *
   * @param scaleFactor a double.
   */
  public void setTorsionScale(double scaleFactor) {
    for (TorsionType type : torsionTypes.values()) {
      type.setScaleFactor(scaleFactor);
    }
    for (PiTorsionType type : piTorsionTypes.values()) {
      type.setScaleFactor(scaleFactor);
    }
  }

  /**
   * Return a String for any Force Field keyword.
   *
   * @param type ForceFieldType
   * @return String
   */
  public String toString(ForceFieldType type) {
    StringBuilder sb = new StringBuilder("\n");
    Map<String, ? extends BaseType> t = forceFieldTypes.get(type);

    if (t.size() == 0) {
      return "";
    }

    for (Object o : t.values()) {
      sb.append(o.toString()).append("\n");
    }
    return sb.toString();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    String forceFieldName;
    try {
      forceFieldName = getString("FORCEFIELD");
    } catch (Exception ex) {
      forceFieldName = "Unknown";
    }
    return forceFieldName;
  }

  /**
   * toString
   *
   * @param key a {@link java.lang.String} object.
   * @return a {@link java.lang.String} object.
   */
  public String toString(String key) {
    if (key == null) {
      return null;
    }

    key = toPropertyForm(key);

    if (properties.containsKey(key)) {
      List<Object> l = properties.getList(key);
      return key + " " + Arrays.toString(l.toArray());
    } else {
      return key + " is not defined.";
    }
  }

  /**
   * toStringBuffer
   *
   * @return Returns a StringBuffer representation of the ForceField.
   */
  public StringBuffer toStringBuffer() {
    StringBuffer sb = new StringBuffer();
    for (ForceFieldType s : forceFieldTypes.keySet()) {
      ForceFieldType type = ForceFieldType.valueOf(s.toString());
      sb.append(toString(type));
    }
    return sb;
  }

  /**
   * All atoms whose atomName begins with name in the passed molecule will be updated to the new
   * type. The case where an atom such as CD is should map to both CD1 and CD2.
   *
   * @param molecule a {@link java.lang.String} object.
   * @param atom a {@link java.lang.String} object.
   * @param newType a int.
   * @return a {@link ffx.potential.parameters.AtomType} object.
   */
  private AtomType updateBioType(String molecule, String atom, int newType) {
    int type = 0;
    for (BioType bioType : bioTypes.values()) {
      if (bioType.moleculeName.equalsIgnoreCase(molecule)) {
        if (atom.length() <= bioType.atomName.length()) {
          if (bioType.atomName.toUpperCase().startsWith(atom.toUpperCase())) {
            type = bioType.atomType;
            bioType.atomType = newType;
          }
        }
      }
    }
    return getAtomType(Integer.toString(type));
  }

  /**
   * Patches that add new atom classes/types that bond to existing atom classes/types require
   * "hybrid" force field types that include a mixture of new and existing types.
   *
   * @param typeMap A look-up from new types to existing types.
   * @param patchTypes a {@link java.util.HashMap} object.
   */
  private void patchClassesAndTypes(
      HashMap<AtomType, AtomType> typeMap, HashMap<String, AtomType> patchTypes) {

    for (BondType bondType : bondTypes.values().toArray(new BondType[0])) {
      BondType newType = bondType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    for (AngleType angleType : angleTypes.values().toArray(new AngleType[0])) {
      AngleType newType = angleType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    for (OutOfPlaneBendType outOfPlaneBendType :
        outOfPlaneBendTypes.values().toArray(new OutOfPlaneBendType[0])) {
      OutOfPlaneBendType newType = outOfPlaneBendType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    for (PiTorsionType piTorsionType : piTorsionTypes.values().toArray(new PiTorsionType[0])) {
      PiTorsionType newType = piTorsionType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    for (StretchBendType stretchBendType :
        stretchBendTypes.values().toArray(new StretchBendType[0])) {
      StretchBendType newType = stretchBendType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    /* for (TorsionTorsionType torsionTorsionType :
     * torsionTorsionTypes.values().toArray(new TorsionTorsionType[0])) {
     * String currentKey = torsionTorsionType.key;
     * torsionTorsionType.patchClasses(typeMap); if
     * (!torsionTorsionType.key.equals(currentKey)) {
     * torsionTorsionTypes.remove(currentKey);
     * addForceFieldType(torsionTorsionType); } }
     */

    for (TorsionType torsionType : torsionTypes.values().toArray(new TorsionType[0])) {
      TorsionType newType = torsionType.patchClasses(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    /*
    for (ImproperTorsionType improperType : imptorsTypes.values().toArray(new ImproperTorsionType[0])) {
        String currentKey = improperType.key;
        improperType.patchClasses(typeMap);
        if (!improperType.key.equals(currentKey)) {
            torsionTypes.remove(currentKey);
            addForceFieldType(improperType);
        }
    }

    for (UreyBradleyType ureyBradleyType : ureyBradleyTypes.values().toArray(new UreyBradleyType[0])) {
        String currentKey = ureyBradleyType.key;
        ureyBradleyType.patchClasses(typeMap);
        if (!ureyBradleyType.key.equals(currentKey)) {
            ureyBradleyTypes.remove(currentKey);
            addForceFieldType(ureyBradleyType);
        }
    } */
    for (MultipoleType multipoleType : multipoleTypes.values().toArray(new MultipoleType[0])) {
      MultipoleType newType = multipoleType.patchTypes(typeMap);
      if (newType != null) {
        logger.info(" " + newType.toString());
        addForceFieldType(newType);
      }
    }

    try {
      for (AtomType atomType : patchTypes.values()) {
        PolarizeType polarizeType = getPolarizeType(atomType.key);
        if (polarizeType != null && polarizeType.patchTypes(typeMap)) {
          logger.info(" " + polarizeType.toString());
        }
      }
    } catch (Exception e) {
      // Inefficient hack. Should actually check if polarizeTypes are necessary.
    }
  }

  /**
   * Returns the minimum atom class value.
   *
   * @return a int.
   */
  private int minClass() {
    int minClass = maxClass();
    for (AtomType type : atomTypes.values()) {
      if (type.atomClass < minClass) {
        minClass = type.atomClass;
      }
    }
    return minClass;
  }

  /**
   * Returns the minimum atom type value.
   *
   * @return a int.
   */
  private int minType() {
    int minType = maxType();
    for (String key : atomTypes.keySet()) {
      int type = Integer.parseInt(key);
      if (type < minType) {
        minType = type;
      }
    }
    return minType;
  }

  /**
   * Returns the minimum Biotype value.
   *
   * @return a int.
   */
  private int minBioType() {
    int minBioType = maxBioType();
    for (String key : bioTypes.keySet()) {
      int type = Integer.parseInt(key);
      if (type < minBioType) {
        minBioType = type;
      }
    }
    return minBioType;
  }

  /**
   * Returns the maximum atom class value.
   *
   * @return a int.
   */
  private int maxClass() {
    int maxClass = 0;
    for (AtomType type : atomTypes.values()) {
      if (type.atomClass > maxClass) {
        maxClass = type.atomClass;
      }
    }
    return maxClass;
  }

  /**
   * Returns the maximum atom type value.
   *
   * @return a int.
   */
  private int maxType() {
    int maxType = 0;
    for (String key : atomTypes.keySet()) {
      int type = Integer.parseInt(key);
      if (type > maxType) {
        maxType = type;
      }
    }
    return maxType;
  }

  /**
   * Returns the maximum Biotype.
   *
   * @return a int.
   */
  private int maxBioType() {
    int maxBioType = 0;
    for (String key : bioTypes.keySet()) {
      int type = Integer.parseInt(key);
      if (type > maxBioType) {
        maxBioType = type;
      }
    }
    return maxBioType;
  }

  private void modifiedResidue(String modres) {
    String[] tokens = modres.trim().split(" +");
    String modResname = tokens[0].toUpperCase();
    String stdName = tokens[1].toUpperCase();
    HashMap<String, AtomType> patchAtomTypes = getAtomTypes(modResname);
    HashMap<String, AtomType> stdAtomTypes = getAtomTypes(stdName);

    HashMap<String, AtomType> patchTypes = new HashMap<>();
    int len = tokens.length;
    for (int i = 2; i < len; i++) {
      String atomName = tokens[i].toUpperCase();
      if (!patchTypes.containsKey(atomName) && patchAtomTypes.containsKey(atomName)) {
        AtomType type = patchAtomTypes.get(atomName);
        patchTypes.put(atomName, type);
      }
    }

    HashMap<AtomType, AtomType> typeMap = new HashMap<>();
    for (String atomName : stdAtomTypes.keySet()) {
      boolean found = false;
      for (int i = 2; i < len; i++) {
        if (atomName.equalsIgnoreCase(tokens[i])) {
          found = true;
          break;
        }
      }
      if (!found) {
        AtomType stdType = stdAtomTypes.get(atomName);
        // Edit new BioType records to point to an existing force field type.
        AtomType patchType = updateBioType(modResname, atomName, stdType.type);
        if (patchType != null) {
          typeMap.put(patchType, stdType);
          logger.info(" " + patchType.toString() + " -> " + stdType.toString());
        }
      }
    }

    patchClassesAndTypes(typeMap, patchTypes);
  }

  /**
   * If some set of other boolean values imply another boolean is true, set that implied boolean to
   * true.
   *
   * <p>First designed for LAMBDATERM, which is implied by any of VDW_LAMBDATERM, ELEC_LAMBDATERM,
   * or GK_LAMBDATERM.
   *
   * @param toSet Property to set true if otherBooleans true.
   * @param otherBooleans Properties that imply toSet is true.
   */
  private void trueImpliedBoolean(String toSet, String... otherBooleans) {
    // Short-circuit if it's already true.
    if (getBoolean(toSet, false)) {
      return;
    }
    // Check all the other booleans that imply toSet.
    for (String otherBool : otherBooleans) {
      if (getBoolean(otherBool, false)) {
        addProperty(toSet, "true");
        logger.info(format(" Setting implied boolean %s true due to boolean %s", toSet, otherBool));
      }
    }
  }

  /** Check for self-consistent polarization groups. */
  private void checkPolarizationTypes() {
    boolean change = false;
    for (String key : polarizeTypes.keySet()) {
      PolarizeType polarizeType = polarizeTypes.get(key);
      int orig = Integer.parseInt(key);
      int[] types = polarizeType.polarizationGroup;
      if (types == null) {
        continue;
      }

      for (int type : types) {
        String key2 = Integer.toString(type);
        PolarizeType polarizeType2 = polarizeTypes.get(key2);
        if (polarizeType2 == null) {
          logger.severe(
              format("Polarize type %s references nonexistant polarize type %s.", key, key2));
          continue;
        }
        int[] types2 = polarizeType2.polarizationGroup;
        if (types2 == null) {
          polarizeType2.add(orig);
          change = true;
          continue;
        }
        boolean found = false;
        for (int type2 : types2) {
          for (int type3 : types) {
            if (type2 == type3) {
              found = true;
              break;
            }
          }
          if (!found) {
            polarizeType.add(type2);
            change = true;
          }
        }
      }
    }
    if (change) {
      checkPolarizationTypes();
    }
  }

  public enum ELEC_FORM {
    PAM,
    FIXED_CHARGE
  }

  /** Available force fields. */
  public enum ForceFieldName {
    AMBER_1994,
    AMBER_1996,
    AMBER_1998,
    AMBER_1999,
    AMBER_1999_SB,
    AMBER_1999_SB_AMOEBA,
    AMBER_1999_SB_TIP3F,
    AMOEBA_2004,
    AMOEBA_2009,
    AMOEBA_BIO_2009,
    AMOEBA_BIO_2018,
    AMOEBA_NUC_2017,
    AMOEBA_PROTEIN_2004,
    AMOEBA_PROTEIN_2013,
    AMOEBA_WATER_2003,
    AMOEBA_WATER_2014,
    CHARMM_22,
    CHARMM_22_CMAP,
    IAMOEBA_WATER,
    OPLS_AA,
    OPLS_AAL
  }

  public enum ForceFieldType {
    KEYWORD,
    ANGLE,
    ANGLEP,
    ANGTORS,
    ATOM,
    BIOTYPE,
    BOND,
    CHARGE,
    IMPROPER,
    IMPTORS,
    MULTIPOLE,
    OPBEND,
    PITORS,
    POLARIZE,
    SOLUTE,
    STRBND,
    STRTORS,
    TORSION,
    TORTORS,
    UREYBRAD,
    VDW,
    VDW14,
    RELATIVESOLV
  }

  /** Enumerates the types of constraints that can be applied. */
  public enum ConstraintTypes {
    // Constrain a Bond.
    BOND,
    // Constrain a 3-atom Angle and its two component Bonds.
    ANGLEBONDS
    // TODO: Support dihedral constraints, lone angle constraints, etc.
  }
}
