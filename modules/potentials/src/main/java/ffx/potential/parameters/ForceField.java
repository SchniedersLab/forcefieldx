/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential.parameters;

import java.io.File;
import java.net.URL;
import java.util.*;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;

/**
 * The ForceField class organizes parameters for a molecular mechanics force
 * field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class ForceField {

    private static final Logger logger = Logger.getLogger(ForceField.class.getName());

    /**
     * A map between a Force_Field and its internal parameter file.
     */
    private static final Map<ForceFieldName, URL> forceFields = new EnumMap<>(ForceFieldName.class);

    static {
        ClassLoader cl = ForceField.class.getClassLoader();
        String prefix = "ffx/potential/parameters/ff/";
        for (ForceFieldName ff : ForceFieldName.values()) {
            forceFields.put(ff, cl.getResource(prefix + ff));
        }
    }

    /**
     * <p>
     * Getter for the field <code>forceFieldURL</code>.</p>
     *
     * @param forceField a
     * {@link ffx.potential.parameters.ForceField.ForceFieldName} object.
     * @return a {@link java.net.URL} object.
     */
    public static URL getForceFieldURL(ForceFieldName forceField) {
        if (forceField == null) {
            return null;
        }
        return forceFields.get(forceField);
    }

    public URL forceFieldURL;
    public File keywordFile;
    private final CompositeConfiguration properties;
    private final Map<String, AngleType> angleTypes;
    private final Map<String, AtomType> atomTypes;
    private final Map<String, BioType> bioTypes;
    private final Map<String, BondType> bondTypes;
    private final Map<String, ChargeType> chargeTypes;
    private final Map<String, MultipoleType> multipoleTypes;
    private final Map<String, OutOfPlaneBendType> outOfPlaneBendTypes;
    private final Map<String, PolarizeType> polarizeTypes;
    private final Map<String, StretchBendType> stretchBendTypes;
    private final Map<String, PiTorsionType> piTorsionTypes;
    private final Map<String, TorsionType> torsionTypes;
    private final Map<String, ImproperTorsionType> imptorsTypes;
    private final Map<String, TorsionTorsionType> torsionTorsionTypes;
    private final Map<String, UreyBradleyType> ureyBradleyTypes;
    private final Map<String, VDWType> vanderWaalsTypes;
    private final Map<ForceFieldType, Map> forceFieldTypes;

    /**
     * ForceField Constructor.
     *
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param parameterFile a {@link java.io.File} object.
     */
    public ForceField(CompositeConfiguration properties, File parameterFile) {
        this(properties);
        this.keywordFile = parameterFile;
    }

    /**
     * ForceField Constructor.
     *
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public ForceField(CompositeConfiguration properties) {
        this.properties = properties;
        /**
         * Each force field "type" implements the "Comparator<String>" interface
         * so that passing an "empty" instance of the "type" to its TreeMap
         * constructor will keep the types sorted.
         */
        angleTypes = new TreeMap<>(new AngleType(new int[3], 0, new double[1], null));
        atomTypes = new TreeMap<>(new AtomType(0, 0, null, null, 0, 0, 0));
        bioTypes = new TreeMap<>(new BioType(0, null, null, 0, null));
        bondTypes = new TreeMap<>(new BondType(new int[2], 0, 0, null));
        chargeTypes = new TreeMap<>(new ChargeType(0, 0));
        multipoleTypes = new TreeMap<>(new MultipoleType(0, new double[3], new double[3][3], null, null));
        outOfPlaneBendTypes = new TreeMap<>(new OutOfPlaneBendType(new int[4], 0));
        piTorsionTypes = new TreeMap<>(new PiTorsionType(new int[2], 0));
        polarizeTypes = new TreeMap<>(new PolarizeType(0, 0, 0, new int[1]));
        stretchBendTypes = new TreeMap<>(new StretchBendType(new int[3], new double[1]));
        torsionTorsionTypes = new TreeMap<>();
        torsionTypes = new TreeMap<>(new TorsionType(new int[4], new double[1], new double[1], new int[1]));
        imptorsTypes = new TreeMap<>(new ImproperTorsionType(new int[4], 0.0, 0.0, 2));
        ureyBradleyTypes = new TreeMap<>(new UreyBradleyType(new int[3], 0, 0));
        vanderWaalsTypes = new TreeMap<>(new VDWType(0, 0, 0, 0));

        forceFieldTypes = new EnumMap<>(ForceFieldType.class);
        forceFieldTypes.put(ForceFieldType.ANGLE, angleTypes);
        forceFieldTypes.put(ForceFieldType.ATOM, atomTypes);
        forceFieldTypes.put(ForceFieldType.BOND, bondTypes);
        forceFieldTypes.put(ForceFieldType.BIOTYPE, bioTypes);
        forceFieldTypes.put(ForceFieldType.CHARGE, chargeTypes);
        forceFieldTypes.put(ForceFieldType.OPBEND, outOfPlaneBendTypes);
        forceFieldTypes.put(ForceFieldType.MULTIPOLE, multipoleTypes);
        forceFieldTypes.put(ForceFieldType.PITORS, piTorsionTypes);
        forceFieldTypes.put(ForceFieldType.POLARIZE, polarizeTypes);
        forceFieldTypes.put(ForceFieldType.STRBND, stretchBendTypes);
        forceFieldTypes.put(ForceFieldType.TORSION, torsionTypes);
        forceFieldTypes.put(ForceFieldType.IMPTORS, imptorsTypes);
        forceFieldTypes.put(ForceFieldType.TORTORS, torsionTorsionTypes);
        forceFieldTypes.put(ForceFieldType.UREYBRAD, ureyBradleyTypes);
        forceFieldTypes.put(ForceFieldType.VDW, vanderWaalsTypes);
    }

    public boolean getBoolean(ForceFieldDouble forceFieldDouble, double d) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public int minClass() {
        int minClass = maxClass();
        for (AtomType type : atomTypes.values()) {
            if (type.atomClass < minClass) {
                minClass = type.atomClass;
            }
        }
        return minClass;
    }

    public int minType() {
        int minType = maxType();
        for (String key : atomTypes.keySet()) {
            int type = Integer.parseInt(key);
            if (type < minType) {
                minType = type;
            }
        }
        return minType;
    }

    public int minBioType() {
        int minBioType = maxBioType();
        for (String key : bioTypes.keySet()) {
            int type = Integer.parseInt(key);
            if (type < minBioType) {
                minBioType = type;
            }
        }
        return minBioType;
    }

    public int maxClass() {
        int maxClass = 0;
        for (AtomType type : atomTypes.values()) {
            if (type.atomClass > maxClass) {
                maxClass = type.atomClass;
            }
        }
        return maxClass;
    }

    public int maxType() {
        int maxType = 0;
        for (String key : atomTypes.keySet()) {
            int type = Integer.parseInt(key);
            if (type > maxType) {
                maxType = type;
            }
        }
        return maxType;
    }

    public int maxBioType() {
        int maxBioType = 0;
        for (String key : bioTypes.keySet()) {
            int type = Integer.parseInt(key);
            if (type > maxBioType) {
                maxBioType = type;
            }
        }
        return maxBioType;
    }

    public void renumberForceField(int classOffset, int typeOffset, int bioTypeOffset) {

        for (AngleType patchType : angleTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (AtomType patchType : atomTypes.values()) {
            patchType.incrementClassAndType(classOffset, typeOffset);
        }

        for (BioType patchType : bioTypes.values()) {
            patchType.incrementIndexAndType(bioTypeOffset, typeOffset);
        }

        for (BondType patchType : bondTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (MultipoleType patchType : multipoleTypes.values()) {
            patchType.incrementType(typeOffset);
        }

        for (OutOfPlaneBendType patchType : outOfPlaneBendTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (PiTorsionType patchType : piTorsionTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (PolarizeType patchType : polarizeTypes.values()) {
            patchType.incrementType(typeOffset);
        }

        for (StretchBendType patchType : stretchBendTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (TorsionTorsionType patchType : torsionTorsionTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (TorsionType patchType : torsionTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (ImproperTorsionType patchType : imptorsTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (UreyBradleyType patchType : ureyBradleyTypes.values()) {
            patchType.incrementClasses(classOffset);
        }

        for (VDWType patchType : vanderWaalsTypes.values()) {
            patchType.incrementClass(classOffset);
        }
    }

    /**
     * Append a 2nd ForceField "patch" to the current ForceField. Note that only
     * the force field types are appended; properties are ignored.
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

        for (AngleType patchType : patch.angleTypes.values()) {
            angleTypes.put(patchType.getKey(), patchType);
        }

        for (AtomType patchType : patch.atomTypes.values()) {
            atomTypes.put(patchType.getKey(), patchType);
        }

        for (BioType patchType : patch.bioTypes.values()) {
            bioTypes.put(patchType.getKey(), patchType);
        }

        for (BondType patchType : patch.bondTypes.values()) {
            bondTypes.put(patchType.getKey(), patchType);
        }

        for (MultipoleType patchType : patch.multipoleTypes.values()) {
            multipoleTypes.put(patchType.getKey(), patchType);
        }

        for (OutOfPlaneBendType patchType : patch.outOfPlaneBendTypes.values()) {
            outOfPlaneBendTypes.put(patchType.getKey(), patchType);
        }

        for (PiTorsionType patchType : patch.piTorsionTypes.values()) {
            piTorsionTypes.put(patchType.getKey(), patchType);
        }

        for (PolarizeType patchType : patch.polarizeTypes.values()) {
            polarizeTypes.put(patchType.getKey(), patchType);
        }

        for (StretchBendType patchType : patch.stretchBendTypes.values()) {
            stretchBendTypes.put(patchType.getKey(), patchType);
        }

        for (TorsionTorsionType patchType : patch.torsionTorsionTypes.values()) {
            torsionTorsionTypes.put(patchType.getKey(), patchType);
        }

        for (TorsionType patchType : patch.torsionTypes.values()) {
            torsionTypes.put(patchType.getKey(), patchType);
        }

        for (ImproperTorsionType patchType : patch.imptorsTypes.values()) {
            imptorsTypes.put(patchType.getKey(), patchType);
        }

        for (UreyBradleyType patchType : patch.ureyBradleyTypes.values()) {
            ureyBradleyTypes.put(patchType.getKey(), patchType);
        }

        for (VDWType patchType : patch.vanderWaalsTypes.values()) {
            vanderWaalsTypes.put(patchType.getKey(), patchType);
        }
    }

    /**
     * Enums are uppercase with underscores, but property files use lower case
     * with dashes.
     *
     * @param s
     * @return
     */
    private String normalizeKey(String s) {
        if (s == null) {
            return null;
        }
        return s.toLowerCase().replace("_", "-");
    }

    /**
     * <p>
     * getDouble</p>
     *
     * @param forceFieldDouble a
     * {@link ffx.potential.parameters.ForceField.ForceFieldDouble} object.
     * @return a double.
     * @throws java.lang.Exception if any.
     */
    public double getDouble(ForceFieldDouble forceFieldDouble)
            throws Exception {
        if (forceFieldDouble == null) {
            throw new Exception("NULL keyword");
        }
        String key = normalizeKey(forceFieldDouble.toString());
        if (!properties.containsKey(key)) {
            throw new Exception("Undefined keyword: " + key);
        }
        return properties.getDouble(key);
    }

    /**
     * <p>
     * getDouble</p>
     *
     * @param forceFieldDouble a
     * {@link ffx.potential.parameters.ForceField.ForceFieldDouble} object.
     * @param defaultValue a double.
     * @return a double.
     */
    public double getDouble(ForceFieldDouble forceFieldDouble,
            double defaultValue) {
        try {
            return getDouble(forceFieldDouble);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    /**
     * <p>
     * getInteger</p>
     *
     * @param forceFieldInteger a
     * {@link ffx.potential.parameters.ForceField.ForceFieldInteger} object.
     * @return a int.
     * @throws java.lang.Exception if any.
     */
    public int getInteger(ForceFieldInteger forceFieldInteger)
            throws Exception {
        if (forceFieldInteger == null) {
            throw new Exception("NULL keyword");
        }
        String key = normalizeKey(forceFieldInteger.toString());
        if (!properties.containsKey(key)) {
            throw new Exception("Undefined keyword: " + key);
        }
        return properties.getInt(key);
    }

    /**
     * <p>
     * getInteger</p>
     *
     * @param forceFieldInteger a
     * {@link ffx.potential.parameters.ForceField.ForceFieldInteger} object.
     * @param defaultValue a int.
     * @return a int.
     */
    public int getInteger(ForceFieldInteger forceFieldInteger,
            int defaultValue) {
        try {
            return getInteger(forceFieldInteger);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    /**
     * <p>
     * getString</p>
     *
     * @param forceFieldString a
     * {@link ffx.potential.parameters.ForceField.ForceFieldString} object.
     * @return a {@link java.lang.String} object.
     * @throws java.lang.Exception if any.
     */
    public String getString(ForceFieldString forceFieldString)
            throws Exception {
        if (forceFieldString == null) {
            throw new Exception("NULL keyword");
        }
        String key = normalizeKey(forceFieldString.toString());
        if (!properties.containsKey(key)) {
            throw new Exception("Undefined keyword: " + key);
        }
        return properties.getString(key);
    }

    /**
     * <p>
     * getString</p>
     *
     * @param forceFieldString a
     * {@link ffx.potential.parameters.ForceField.ForceFieldString} object.
     * @param defaultString a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getString(ForceFieldString forceFieldString,
            String defaultString) {
        try {
            return getString(forceFieldString);
        } catch (Exception e) {
            return defaultString;
        }
    }

    /**
     * <p>
     * getBoolean</p>
     *
     * @param forceFieldBoolean a
     * {@link ffx.potential.parameters.ForceField.ForceFieldBoolean} object.
     * @return a boolean.
     * @throws java.lang.Exception if any.
     */
    public boolean getBoolean(ForceFieldBoolean forceFieldBoolean)
            throws Exception {
        if (forceFieldBoolean == null) {
            throw new Exception("NULL keyword");
        }
        String key = normalizeKey(forceFieldBoolean.toString());
        if (!properties.containsKey(key)) {
            throw new Exception("Undefined keyword: " + key);
        }
        return properties.getBoolean(key);
    }

    /**
     * <p>
     * getBoolean</p>
     *
     * @param forceFieldBoolean a
     * {@link ffx.potential.parameters.ForceField.ForceFieldBoolean} object.
     * @param defaultBoolean a {@link java.lang.Boolean} object.
     * @return a boolean.
     */
    public boolean getBoolean(ForceFieldBoolean forceFieldBoolean,
            Boolean defaultBoolean) {
        try {
            return getBoolean(forceFieldBoolean);
        } catch (Exception e) {
            return defaultBoolean;
        }
    }

    /**
     * Add a force field keyword that is represented by a double.
     *
     * @param forceFieldDouble ForceFieldDouble
     * @param value double
     */
    public void addForceFieldDouble(ForceFieldDouble forceFieldDouble, double value) {
        if (forceFieldDouble == null) {
            return;
        }
        String key = normalizeKey(forceFieldDouble.toString());
        try {
            double old = getDouble(forceFieldDouble);
            properties.clearProperty(key);
            logger.info(String.format(" Removed %s with value %8.3f.", key, old));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" Added %s with value %8.3f.", key, value));
        }
    }

    /**
     * Add a force field keyword that is represented by an int.
     *
     * @param forceFieldInteger ForceFieldInteger
     * @param value int
     */
    public void addForceFieldInteger(ForceFieldInteger forceFieldInteger, int value) {
        if (forceFieldInteger == null) {
            return;
        }
        String key = normalizeKey(forceFieldInteger.toString());
        try {
            int old = getInteger(forceFieldInteger);
            logger.info(String.format(" Clearing %s with value %d.", key, old));
            properties.clearProperty(key);
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            logger.info(String.format(" Adding %s with value %d.", key, value));
            properties.addProperty(key, value);
        }
    }

    /**
     * Store a force field keyword that is represented by a String.
     *
     * @param forceFieldString ForceFieldString
     * @param value String
     */
    public void addForceFieldString(ForceFieldString forceFieldString, String value) {
        if (forceFieldString == null) {
            return;
        }
        String key = normalizeKey(forceFieldString.toString());
        try {
            String old = getString(forceFieldString);
            properties.clearProperty(key);
            logger.info(String.format(" Removed %s with value %s.", key, old));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" Added %s with value %s.", key, value));
        }
    }

    /**
     * Store a force field keyword that is represented by a Boolean.
     *
     * @param forceFieldBoolean ForceFielBoolean
     * @param value Boolean
     */
    public void addForceFieldBoolean(ForceFieldBoolean forceFieldBoolean, Boolean value) {
        if (forceFieldBoolean == null) {
            return;
        }
        String key = normalizeKey(forceFieldBoolean.toString());
        try {
            boolean old = getBoolean(forceFieldBoolean);
            properties.clearProperty(key);
            logger.info(String.format(" Cleared %s with value %s.", key, Boolean.toString(old)));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" Added %s with value %s.", key, Boolean.toString(value)));
        }
    }

    /**
     * Add an instance of a force field type. Force Field types are more
     * complicated than simple Strings or doubles, in that they have multiple
     * fields and may occur multiple times.
     *
     * @param type BaseType
     */
    public void addForceFieldType(BaseType type) {
        if (type == null) {
            logger.info(" Null force field type ignored.");
            return;
        }

        Map treeMap = forceFieldTypes.get(type.forceFieldType);
        if (treeMap == null) {
            logger.info(" Unrecognized force field type ignored " + type.forceFieldType);
            type.print();
            return;
        }
        if (treeMap.containsKey(type.key)) {
            logger.warning(" A force field entry of type " + type.forceFieldType
                    + " already exists with the key: " + type.key
                    + "\n The (discarded) old entry: " + treeMap.get(type.key).
                    toString() + "\n The new entry            : " + type.toString());
        }
        Class baseTypeClass = type.getClass();
        treeMap.put(type.key, baseTypeClass.cast(type));
    }

    /**
     * <p>
     * getAngleType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.AngleType} object.
     */
    public AngleType getAngleType(String key) {
        return angleTypes.get(key);
    }

    /**
     * <p>
     * getAtomType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.AtomType} object.
     */
    public AtomType getAtomType(String key) {
        return atomTypes.get(key);
    }

    /**
     * <p>
     * getAtomType</p>
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
     * <p>
     * Getter for the field <code>atomTypes</code>.</p>
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
     * <p>
     * getBonds</p>
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
     * <p>
     * getBondType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.BondType} object.
     */
    public BondType getBondType(String key) {
        return bondTypes.get(key);
    }

    /**
     * <p>
     * getBioType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.BioType} object.
     */
    public BioType getBioType(String key) {
        return bioTypes.get(key);
    }

    /**
     * <p>
     * getMultipoleType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.MultipoleType} object.
     */
    public MultipoleType getMultipoleType(String key) {
        return multipoleTypes.get(key);
    }

    /**
     * <p>
     * getOutOfPlaneBendType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     */
    public OutOfPlaneBendType getOutOfPlaneBendType(String key) {
        return outOfPlaneBendTypes.get(key);
    }

    /**
     * <p>
     * getPolarizeType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.PolarizeType} object.
     */
    public PolarizeType getPolarizeType(String key) {
        return polarizeTypes.get(key);
    }

    /**
     * <p>
     * getStretchBendType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.StretchBendType} object.
     */
    public StretchBendType getStretchBendType(String key) {
        return stretchBendTypes.get(key);
    }

    /**
     * <p>
     * getPiTorsionType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.PiTorsionType} object.
     */
    public PiTorsionType getPiTorsionType(String key) {
        return piTorsionTypes.get(key);
    }

    /**
     * <p>
     * getTorsionType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.TorsionType} object.
     */
    public TorsionType getTorsionType(String key) {
        return torsionTypes.get(key);
    }

    /**
     * <p>
     * getImproperType</p>
     *
     * @return a {@link ffx.potential.parameters.TorsionType} object.
     */
    public Collection<ImproperTorsionType> getImproperTypes() {
        return imptorsTypes.values();
    }

    /**
     * <p>
     * getImproperType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.TorsionType} object.
     */
    public ImproperTorsionType getImproperType(String key) {
        return imptorsTypes.get(key);
    }

    /**
     * <p>
     * getTorsionTorsionType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.TorsionTorsionType} object.
     */
    public TorsionTorsionType getTorsionTorsionType(String key) {
        return torsionTorsionTypes.get(key);
    }

    /**
     * <p>
     * getUreyBradleyType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.UreyBradleyType} object.
     */
    public UreyBradleyType getUreyBradleyType(String key) {
        return ureyBradleyTypes.get(key);
    }

    /**
     * <p>
     * getVDWType</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link ffx.potential.parameters.VDWType} object.
     */
    public VDWType getVDWType(String key) {
        return vanderWaalsTypes.get(key);
    }

    /**
     * <p>
     * getVDWTypes</p>
     *
     * @return a {@link java.util.Map} object.
     */
    public Map<String, VDWType> getVDWTypes() {
        return vanderWaalsTypes;
    }

    /**
     * <p>
     * getForceFieldTypeCount</p>
     *
     * @param type a {@link ffx.potential.parameters.ForceField.ForceFieldType}
     * object.
     * @return a int.
     */
    public int getForceFieldTypeCount(ForceFieldType type) {
        TreeMap<String, BaseType> treeMap
                = (TreeMap<String, BaseType>) forceFieldTypes.get(type);
        if (treeMap == null) {
            logger.warning("Unrecognized Force Field Type: " + type);
            return 0;
        }
        return treeMap.size();
    }

    /**
     * Check for self-consistent polarization groups.
     */
    public void checkPolarizationTypes() {
        boolean change = false;
        for (String key : polarizeTypes.keySet()) {
            PolarizeType polarizeType = polarizeTypes.get(key);
            int orig = Integer.parseInt(key);
            int types[] = polarizeType.polarizationGroup;
            if (types == null) {
                continue;
            }
            for (int i = 0; i < types.length; i++) {
                int type = types[i];
                String key2 = Integer.toString(type);
                PolarizeType polarizeType2 = polarizeTypes.get(key2);
                if (polarizeType2 == null) {
                    logger.severe(format("Polarize type %s references nonexistant polarize type %s.",
                            key, key2));
                    continue;
                }
                int types2[] = polarizeType2.polarizationGroup;
                if (types2 == null) {
                    polarizeType2.add(orig);
                    change = true;
                    continue;
                }
                boolean found = false;
                for (int j = 0; j < types2.length; j++) {
                    int type2 = types2[j];
                    for (int k = 0; k < types.length; k++) {
                        if (type2 == types[k]) {
                            found = true;
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

    /**
     * <p>
     * log</p>
     */
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

    /**
     * <p>
     * print</p>
     */
    public void print() {
        for (ForceFieldType s : forceFieldTypes.keySet()) {
            print(s.toString());
        }
    }

    /**
     * <p>
     * print</p>
     *
     * @param key a {@link java.lang.String} object.
     */
    public void print(String key) {
        ForceFieldType type = ForceFieldType.valueOf(key);
        System.out.println(toString(type));
    }

    /**
     * Return a String for any Force Field keyword.
     *
     * @param type ForceFieldType
     * @return String
     */
    public String toString(ForceFieldType type) {
        StringBuilder sb = new StringBuilder("\n");
        Map t = forceFieldTypes.get(type);
        if (t.size() > 0) {
            Set set = t.keySet();
            Iterator i = set.iterator();
            while (i.hasNext()) {
                sb.append(t.get(i.next())).append("\n");
            }
            return sb.toString();
        }
        return "";
    }

    @Override
    public String toString() {
        String forceFieldName;
        try {
            forceFieldName = getString(ForceFieldString.FORCEFIELD);
        } catch (Exception ex) {
            forceFieldName = "Unknown";
        }
        return forceFieldName;
    }

    /**
     * <p>
     * toString</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String toString(String key) {
        if (key == null) {
            return null;
        }

        key = normalizeKey(key);

        if (properties.containsKey(key)) {
            List l = properties.getList(key);
            return key + " " + Arrays.toString(l.toArray());
        } else {
            return key + " is not defined.";
        }
    }

    /**
     * Patches that add new atom classes/types that bond to existing atom
     * classes/types require "hybrid" force field types that include a mixture
     * of new and existing types.
     *
     * @param typeMap A look-up from new types to existing types.
     */
    public void patchClassesAndTypes(HashMap<AtomType, AtomType> typeMap) {

        for (AngleType angleType : angleTypes.values().toArray(new AngleType[angleTypes.size()])) {
            String currentKey = angleType.key;
            angleType.patchClasses(typeMap);
            if (!angleType.key.equals(currentKey)) {
                angleTypes.remove(currentKey);
                addForceFieldType(angleType);
            }
        }

        for (BondType bondType : bondTypes.values().toArray(new BondType[bondTypes.size()])) {
            String currentKey = bondType.key;
            bondType.patchClasses(typeMap);
            if (!bondType.key.equals(currentKey)) {
                bondTypes.remove(currentKey);
                addForceFieldType(bondType);
            }
        }

        for (OutOfPlaneBendType outOfPlaneBendType : outOfPlaneBendTypes.values().toArray(new OutOfPlaneBendType[outOfPlaneBendTypes.size()])) {
            String currentKey = outOfPlaneBendType.key;
            outOfPlaneBendType.patchClasses(typeMap);
            if (!outOfPlaneBendType.key.equals(currentKey)) {
                outOfPlaneBendTypes.remove(currentKey);
                addForceFieldType(outOfPlaneBendType);
            }
        }

        for (PiTorsionType piTorsionType : piTorsionTypes.values().toArray(new PiTorsionType[piTorsionTypes.size()])) {
            String currentKey = piTorsionType.key;
            piTorsionType.patchClasses(typeMap);
            if (!piTorsionType.key.equals(currentKey)) {
                piTorsionTypes.remove(currentKey);
                addForceFieldType(piTorsionType);
            }
        }

        for (StretchBendType stretchBendType : stretchBendTypes.values().toArray(new StretchBendType[stretchBendTypes.size()])) {
            String currentKey = stretchBendType.key;
            stretchBendType.patchClasses(typeMap);
            if (!stretchBendType.key.equals(currentKey)) {
                stretchBendTypes.remove(currentKey);
                addForceFieldType(stretchBendType);
            }
        }

        for (TorsionTorsionType torsionTorsionType : torsionTorsionTypes.values().toArray(new TorsionTorsionType[torsionTorsionTypes.size()])) {
            String currentKey = torsionTorsionType.key;
            torsionTorsionType.patchClasses(typeMap);
            if (!torsionTorsionType.key.equals(currentKey)) {
                torsionTorsionTypes.remove(currentKey);
                addForceFieldType(torsionTorsionType);
            }
        }

        for (TorsionType torsionType : torsionTypes.values().toArray(new TorsionType[torsionTypes.size()])) {
            String currentKey = torsionType.key;
            torsionType.patchClasses(typeMap);
            if (!torsionType.key.equals(currentKey)) {
                torsionTypes.remove(currentKey);
                addForceFieldType(torsionType);
            }
        }

        for (ImproperTorsionType improperType : imptorsTypes.values().toArray(new ImproperTorsionType[imptorsTypes.size()])) {
            String currentKey = improperType.key;
            improperType.patchClasses(typeMap);
            if (!improperType.key.equals(currentKey)) {
                torsionTypes.remove(currentKey);
                addForceFieldType(improperType);
            }
        }

        for (UreyBradleyType ureyBradleyType : ureyBradleyTypes.values().toArray(new UreyBradleyType[ureyBradleyTypes.size()])) {
            String currentKey = ureyBradleyType.key;
            ureyBradleyType.patchClasses(typeMap);
            if (!ureyBradleyType.key.equals(currentKey)) {
                ureyBradleyTypes.remove(currentKey);
                addForceFieldType(ureyBradleyType);
            }
        }

        for (MultipoleType multipoleType : multipoleTypes.values().toArray(new MultipoleType[multipoleTypes.size()])) {
            String currentKey = multipoleType.key;
            multipoleType.patchTypes(typeMap);
            if (!multipoleType.key.equals(currentKey)) {
                multipoleTypes.remove(currentKey);
                addForceFieldType(multipoleType);
            }
        }

        for (PolarizeType polarizeType : polarizeTypes.values().toArray(new PolarizeType[polarizeTypes.size()])) {
            String currentKey = polarizeType.key;
            polarizeType.patchTypes(typeMap);
            if (!polarizeType.key.equals(currentKey)) {
                polarizeTypes.remove(currentKey);
                addForceFieldType(polarizeType);
            }
        }
    }

    /**
     * Available force fields.
     */
    public enum ForceFieldName {

        AMOEBA_WATER, AMOEBA_2004, AMOEBA_PROTEIN_2004, AMOEBA_PROTEIN_2004_U1, AMOEBA_2009, AMOEBA_BIO_2009, AMOEBA_BIO_2009_ORIG, AMOEBA_PROTEIN_2013, OPLS_AAL
    }

    public enum ForceFieldString {

        SPACEGROUP, NCSGROUP, FORCEFIELD, POLARIZATION, VDW_SCHEDULE, REAL_SCHEDULE, RECIP_SCHEDULE, SCF_PREDICTOR, SCF_ALGORITHM
    }

    public enum ForceFieldDouble {

        A_AXIS, B_AXIS, C_AXIS, ALPHA, BETA, GAMMA, POLAR_DAMP, POLAR_SOR, POLAR_EPS, POLAR_EPS_PRECISE, EWALD_CUTOFF, EWALD_ALPHA, EWALD_PRECISION, PME_MESH_DENSITY, VDW_CUTOFF, MPOLE_11_SCALE, MPOLE_12_SCALE, MPOLE_13_SCALE, MPOLE_14_SCALE, MPOLE_15_SCALE, POLAR_12_SCALE, POLAR_13_SCALE, DIRECT_11_SCALE, RIGID_SCALE, VDW_LAMBDA_EXPONENT, VDW_LAMBDA_ALPHA, PERMANENT_LAMBDA_EXPONENT, PERMANENT_LAMBDA_ALPHA, POLARIZATION_LAMBDA_START, POLARIZATION_LAMBDA_END, POLARIZATION_LAMBDA_EXPONENT, DUAL_TOPOLOGY_LAMBDA_EXPONENT, CG_PRECONDITIONER_CUTOFF, CG_PRECONDITIONER_EWALD, CG_PRECONDITIONER_SOR,
        RESTRAINT_K, PROBE_RADIUS
    }

    public enum ForceFieldInteger {

        PME_ORDER, PME_REAL_THREADS, PME_GRID_X, PME_GRID_Y, PME_GRID_Z, LIGAND_START, LIGAND_STOP, SCF_CYCLES, SCF_PREDICTOR_ORDER
    }

    public enum ForceFieldBoolean {

        BONDTERM, ANGLETERM, RESTRAINTERM, STRBNDTERM, UREYTERM, OPBENDTERM, TORSIONTERM, PITORSTERM, TORTORTERM, VDWLRTERM, VDWTERM, IMPROPERTERM, MPOLETERM, POLARIZETERM, GKTERM, SCFCACHE, APERIODIC, CUDAFFT, RIGID_HYDROGENS, LAMBDATERM, NCSTERM, USE_CHARGES, USE_DIPOLES, USE_QUADRUPOLES, ROTATE_MULTIPOLES, LIGAND_VAPOR_ELEC, NO_LIGAND_CONDENSED_SCF, USE_SCF_PRECONDITIONER, INTERMOLECULAR_SOFTCORE, LAMBDA_VALENCE_RESTRAINTS
    }

    public enum ForceFieldType {

        KEYWORD, ANGLE, ATOM, BIOTYPE, BOND, CHARGE, IMPTORS, MULTIPOLE, OPBEND, PITORS, POLARIZE, STRBND, TORSION, TORTORS, UREYBRAD, VDW
    }
}
