/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.parameters;

import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
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
     * A map between a force field name and its internal parameter file.
     */
    private static final Map<ForceFieldName, URL> forceFields = new EnumMap<>(ForceFieldName.class);
    private static final boolean noRenumbering = System.getProperty("noPatchRenumbering") != null;

    static {
        ClassLoader cl = ForceField.class.getClassLoader();
        String prefix = "ffx/potential/parameters/ff/";
        for (ForceFieldName ff : ForceFieldName.values()) {
            forceFields.put(ff, cl.getResource(prefix + ff));
        }
    }

    /**
     * <p>
     * Get for the URL for the named force field.</p>
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

    /**
     * URL to the force field parameter file.
     */
    public URL forceFieldURL;
    /**
     * Reference to the keyword file in use.
     */
    public File keywordFile;
    /**
     * The CompositeConfiguration that contains key=value property pairs from a
     * number of sources.
     */
    private final CompositeConfiguration properties;
    private final Map<String, AngleType> angleTypes;
    private final Map<String, AtomType> atomTypes;
    private final Map<String, BioType> bioTypes;
    private final Map<String, BondType> bondTypes;
    private final Map<String, ChargeType> chargeTypes;
    private final Map<String, ISolvRadType> iSolvRadTypes;
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
    private final Map<String, RelativeSolvationType> relativeSolvationTypes;
    private final Map<ForceFieldType, Map<String,? extends BaseType>> forceFieldTypes;

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
        iSolvRadTypes = new TreeMap<>(new ISolvRadType(0, 0.0));
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
        relativeSolvationTypes = new TreeMap<>(new RelativeSolvationType("", 0.0));

        forceFieldTypes = new EnumMap<>(ForceFieldType.class);
        forceFieldTypes.put(ForceFieldType.ANGLE, angleTypes);
        forceFieldTypes.put(ForceFieldType.ATOM, atomTypes);
        forceFieldTypes.put(ForceFieldType.BOND, bondTypes);
        forceFieldTypes.put(ForceFieldType.BIOTYPE, bioTypes);
        forceFieldTypes.put(ForceFieldType.CHARGE, chargeTypes);
        forceFieldTypes.put(ForceFieldType.ISOLVRAD, iSolvRadTypes);
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
        forceFieldTypes.put(ForceFieldType.RELATIVESOLV, relativeSolvationTypes);
    }

    /**
     * Returns the minimum atom class value.
     *
     * @return
     */
    public int minClass() {
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
     * @return
     */
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

    /**
     * Returns the minimum Biotype value.
     *
     * @return
     */
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

    /**
     * Returns the maximum atom class value.
     *
     * @return
     */
    public int maxClass() {
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
     * @return
     */
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

    /**
     * Returns the maximum Biotype.
     *
     * @return
     */
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

    /**
     * Renumber ForceField class, type and biotype values.
     *
     * @param classOffset
     * @param typeOffset
     * @param bioTypeOffset
     */
    public void renumberForceField(int classOffset, int typeOffset, int bioTypeOffset) {
        if (noRenumbering) {
            return;
        }
        for (AngleType angleType : angleTypes.values()) {
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

        for (TorsionTorsionType torsionTorsionType : torsionTorsionTypes.values()) {
            torsionTorsionType.incrementClasses(classOffset);
        }

        for (TorsionType torsionType : torsionTypes.values()) {
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

        for (ISolvRadType iSolvRadType : iSolvRadTypes.values()) {
            iSolvRadType.incrementType(typeOffset);
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

        for (AngleType angleType : patch.angleTypes.values()) {
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

        for (TorsionTorsionType torsionTorsionType : patch.torsionTorsionTypes.values()) {
            torsionTorsionTypes.put(torsionTorsionType.getKey(), torsionTorsionType);
        }

        for (TorsionType torsionType : patch.torsionTypes.values()) {
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

        for (ISolvRadType iSolvRadType : patch.iSolvRadTypes.values()) {
            iSolvRadTypes.put(iSolvRadType.getKey(), iSolvRadType);
//            logger.info(" Adding iSolvRadType to map: " + iSolvRadType.getKey() + "/" + iSolvRadType);
        }

        for (RelativeSolvationType rsType : patch.relativeSolvationTypes.values()) {
            relativeSolvationTypes.put(rsType.getKey(), rsType);
        }

        // Is this a modified residue patch?
        String modres = patch.getString(ForceFieldString.MODRES, "false");
        if (!modres.equalsIgnoreCase("false")) {
            logger.info(" Adding modified residue patch.");
            modifiedResidue(modres);
        }

    }

    private void modifiedResidue(String modres) {
        String tokens[] = modres.trim().split(" +");
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
     * Enums are uppercase with underscores, but property files use lower case
     * with dashes.
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
     * Enums are uppercase with underscores, but property files use lower case
     * with dashes.
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

    public CompositeConfiguration getProperties() {
        return properties;
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
        String key = toPropertyForm(forceFieldDouble.toString());
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
        String key = toPropertyForm(forceFieldInteger.toString());
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
        String key = toPropertyForm(forceFieldString.toString());
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
        String key = toPropertyForm(forceFieldBoolean.toString());
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
        String key = toPropertyForm(forceFieldDouble.toString());
        try {
            double old = getDouble(forceFieldDouble);
            if (old == value) {
                return;
            }
            properties.clearProperty(key);
            logger.info(String.format(" Removed %s %8.3f", key, old));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" %s %8.3f", key, value));
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
        String key = toPropertyForm(forceFieldInteger.toString());
        try {
            int old = getInteger(forceFieldInteger);
            if (old == value) {
                return;
            }
            logger.info(String.format(" Clearing %s %d", key, old));
            properties.clearProperty(key);
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            logger.info(String.format(" %s %d", key, value));
            properties.addProperty(key, value);
        }
    }

    public static boolean isForceFieldKeyword(String keyword) {
        keyword = keyword.toUpperCase();
        try {
            ForceFieldBoolean forceFieldBoolean = ForceFieldBoolean.valueOf(keyword);
            return true;
        } catch (Exception e) {
            // Ignore.
        }
        try {
            ForceFieldDouble forceFieldDouble = ForceFieldDouble.valueOf(keyword);
            return true;
        } catch (Exception e) {
            // Ignore.
        }
        try {
            ForceFieldInteger forceFieldInteger = ForceFieldInteger.valueOf(keyword);
            return true;
        } catch (Exception e) {
            // Ignore.
        }
        try {
            ForceFieldString forceFieldString = ForceFieldString.valueOf(keyword);
            return true;
        } catch (Exception e) {
            // Ignore.
        }
        try {
            ForceFieldType forceFieldType = ForceFieldType.valueOf(keyword);
            return true;
        } catch (Exception e) {
            // Ignore.
        }
        return false;
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
        String key = toPropertyForm(forceFieldString.toString());
        try {
            String old = getString(forceFieldString);
            if (old.equalsIgnoreCase(value)) {
                return;
            }
            properties.clearProperty(key);
            logger.info(String.format(" Removed %s %s.", key, old));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" %s %s", key, value));
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
        String key = toPropertyForm(forceFieldBoolean.toString());
        try {
            boolean old = getBoolean(forceFieldBoolean);
            if (old == value) {
                return;
            }
            properties.clearProperty(key);
            logger.info(String.format(" Cleared %s %s", key, Boolean.toString(old)));
        } catch (Exception e) {
            // Property does not exist yet.
        } finally {
            properties.addProperty(key, value);
            logger.info(String.format(" %s %s", key, Boolean.toString(value)));
        }
    }

    /**
     * Add an instance of a force field type. Force Field types are more
     * complicated than simple Strings or doubles, in that they have multiple
     * fields and may occur multiple times.
     *
     * @param type BaseType
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
            logger.log(Level.WARNING,
                    " A force field entry of type {0} already exists with the key: {1}\n The (discarded) old entry: {2}\n The new entry            : {3}",
                    new Object[]{type.forceFieldType, type.key,
                        treeMap.get(type.key).toString(), type.toString()});
        }
        Class<? extends BaseType> baseTypeClass = type.getClass();
        treeMap.put(type.key, (T) baseTypeClass.cast(type));
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
     * All atoms whose atomName begins with name in the passed molecule will be
     * updated to the new type. The case where an atom such as CD is should map
     * to both CD1 and CD2.
     *
     * @param molecule
     * @param atom
     * @param newType
     * @return
     */
    public AtomType updateBioType(String molecule, String atom, int newType) {
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

    public Map<String, BioType> getBioTypeMap() {
        return bioTypes;
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

    public HashMap<String, RelativeSolvationType> getRelativeSolvationTypes() {
        HashMap<String, RelativeSolvationType> types = new HashMap<>();
        for (String key : relativeSolvationTypes.keySet()) {
            types.put(key, relativeSolvationTypes.get(key));
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

    public ISolvRadType getISolvRadType(String key) {
        return iSolvRadTypes.get(key);
    }

    public Map<String, ISolvRadType> getISolvRadTypes() {
        return iSolvRadTypes;
    }

    /**
     * <p>
     * getForceFieldTypeCount</p>
     *
     * @param type a {@link ffx.potential.parameters.ForceField.ForceFieldType}
     * object.
     * @return a int.
     */
    @SuppressWarnings("unchecked")
    public int getForceFieldTypeCount(ForceFieldType type) {
        TreeMap<String, BaseType> treeMap
                = (TreeMap<String, BaseType>) forceFieldTypes.get(type);
        if (treeMap == null) {
            logger.log(Level.WARNING, "Unrecognized Force Field Type: {0}", type);
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

    public void printMultipoleTypes() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Loaded multipole types: \n"));
        for (String key : multipoleTypes.keySet()) {
            sb.append(String.format(" m %s : \n", key, multipoleTypes.get(key).key));
        }
        logger.info(sb.toString());
    }

    public void printAtomTypes() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Loaded atom types: \n"));
        for (String key : atomTypes.keySet()) {
            sb.append(String.format(" a %s : \n", key, atomTypes.get(key).key));
        }
        logger.info(sb.toString());
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

        key = toPropertyForm(key);

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
    public void patchClassesAndTypes(HashMap<AtomType, AtomType> typeMap, HashMap<String, AtomType> patchTypes) {

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

        for (OutOfPlaneBendType outOfPlaneBendType
                : outOfPlaneBendTypes.values().toArray(new OutOfPlaneBendType[0])) {
            OutOfPlaneBendType newType = outOfPlaneBendType.patchClasses(typeMap);
            if (newType != null) {
                logger.info(" " + newType.toString());
                addForceFieldType(newType);
            }
        }

        for (PiTorsionType piTorsionType
                : piTorsionTypes.values().toArray(new PiTorsionType[0])) {
            PiTorsionType newType = piTorsionType.patchClasses(typeMap);
            if (newType != null) {
                logger.info(" " + newType.toString());
                addForceFieldType(newType);
            }
        }

        for (StretchBendType stretchBendType
                : stretchBendTypes.values().toArray(new StretchBendType[0])) {
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
     * Available force fields.
     */
    public enum ForceFieldName {

        AMBER99SB,
        AMBER99SB_TIP3F,
        AMBER99SB_TIP3AMOEBA,
        AMOEBA_WATER,
        AMOEBA_WATER_2003,
        AMOEBA_WATER_2014,
        AMOEBA_2004,
        AMOEBA_PROTEIN_2004,
        AMOEBA_PROTEIN_2004_U1,
        AMOEBA_2009,
        AMOEBA_BIO_2009,
        AMOEBA_BIO_2009_ORIG,
        AMOEBA_PROTEIN_2013,
        IAMOEBA_WATER,
        OPLS_AAL
    }

    public enum ForceFieldString {

        CAVMODEL,
        ARRAY_REDUCTION,
        EPSILONRULE,
        FFT_METHOD,
        FORCEFIELD,
        NCSGROUP,
        MODRES,
        POLARIZATION,
        RADIUSRULE,
        RADIUSSIZE,
        RADIUSTYPE,
        REAL_SCHEDULE,
        RECIP_SCHEDULE,
        RELATIVE_SOLVATION,
        SCF_PREDICTOR,
        SCF_ALGORITHM,
        SPACEGROUP,
        VDWINDEX,
        VDWTYPE,
        VDW_SCHEDULE,
        GK_RADIIOVERRIDE,
        GK_RADIIBYNUMBER
    }

    public enum ForceFieldDouble {

        /* Unit cell parameters */
        A_AXIS, B_AXIS, C_AXIS, ALPHA, BETA, GAMMA,
        /* Van der Waals' cutoff and softcoring parameters */
        VDW_CUTOFF, VDW_LAMBDA_EXPONENT, VDW_LAMBDA_ALPHA,
        /* Van der Waals masking rules */
        VDW_12_SCALE, VDW_13_SCALE, VDW_14_SCALE, VDW_15_SCALE,
        /* Polarization parameters */
        POLAR_DAMP, POLAR_SOR, POLAR_EPS, POLAR_EPS_PRECISE,
        CG_PRECONDITIONER_CUTOFF, CG_PRECONDITIONER_EWALD, CG_PRECONDITIONER_SOR,
        /* Polarization masking rules */
        POLAR_12_SCALE, POLAR_13_SCALE, DIRECT_11_SCALE, RIGID_SCALE,
        /* Electrostatics parameters */
        EWALD_CUTOFF, EWALD_ALPHA, EWALD_PRECISION, PME_MESH_DENSITY,
        /* Electrostatics masking rules */
        MPOLE_11_SCALE, MPOLE_12_SCALE, MPOLE_13_SCALE, MPOLE_14_SCALE, MPOLE_15_SCALE,
        /* Permanent electrostatics softcoring  */
        PERMANENT_LAMBDA_EXPONENT, PERMANENT_LAMBDA_ALPHA,
        /* Permanent electrostatics lambda window */
        PERMANENT_LAMBDA_START, PERMANENT_LAMBDA_END,
        /* Polarization softcoring */
        POLARIZATION_LAMBDA_EXPONENT, POLARIZATION_LAMBDA_ALPHA,
        /* Polarization lambda window */
        POLARIZATION_LAMBDA_START, POLARIZATION_LAMBDA_END,
        DUAL_TOPOLOGY_LAMBDA_EXPONENT,
        /* Generalized Kirkwood dielectric and debugging */
        GK_EPSILON, GK_OVERLAPSCALE, GK_BONDIOVERRIDE,
        /* Miscellaneous */
        RESTRAINT_K, PROBE_RADIUS, BORNAI, SURFACE_TENSION, TORSIONUNIT, IMPTORUNIT
    }

    public enum ForceFieldInteger {

        FF_THREADS,
        PME_ORDER,
        PME_REAL_THREADS,
        PME_GRID_X,
        PME_GRID_Y,
        PME_GRID_Z,
        LIGAND_START,
        LIGAND_STOP,
        SCF_CYCLES,
        SCF_PREDICTOR_ORDER
    }

    public enum ForceFieldBoolean {

        BONDTERM, ANGLETERM, RESTRAINTERM, COMRESTRAINTERM, STRBNDTERM, UREYTERM,
        OPBENDTERM, TORSIONTERM, PITORSTERM, TORTORTERM, VDWLRTERM, VDWTERM, IMPROPERTERM,
        MPOLETERM, POLARIZETERM, GKTERM, SCFCACHE, APERIODIC, RIGID_HYDROGENS,
        LAMBDATERM, NCSTERM, USE_CHARGES, USE_DIPOLES, USE_QUADRUPOLES, ROTATE_MULTIPOLES,
        LIGAND_VAPOR_ELEC, NO_LIGAND_CONDENSED_SCF, USE_SCF_PRECONDITIONER, INTERMOLECULAR_SOFTCORE,
        INTRAMOLECULAR_SOFTCORE, LAMBDA_VALENCE_RESTRAINTS, LAMBDA_TORSIONS, RECIPTERM, BORN_USE_ALL,
        CHECK_ALL_NODE_CHARGES, GK_USEFITRADII, GK_VERBOSERADII, PRINT_ON_FAILURE
    }

    public enum ForceFieldType {

        KEYWORD,
        ANGLE,
        ATOM,
        BIOTYPE,
        BOND,
        CHARGE,
        ISOLVRAD,
        IMPTORS,
        MULTIPOLE,
        OPBEND,
        PITORS,
        POLARIZE,
        STRBND,
        TORSION,
        TORTORS,
        UREYBRAD,
        VDW,
        RELATIVESOLV
    }

    /**
     * 1) Read in protein with HETATM side-chain. 2) Read in side-chain analog
     * patch and add to ForceField. 3) Define HETATM precedence via 1) existing
     * side-chain or 2) patch. 4) Create needed parameters via averaging. A)
     * Create a map between patch atoms and forcefield atoms. B) Loop over all
     * patch terms and create averaged ones if there is mixed precedence.
     */
}
