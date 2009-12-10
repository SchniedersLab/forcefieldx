/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.parameters;

import java.io.File;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 * The ForceField class organizes parameters for a molecular mechanics force
 * field.
 */
public class ForceField {

    private static final Logger log = Logger.getLogger("ffx");

    public enum ForceFieldString {

        SPACEGROUP, EPSILONRULE, FORCEFIELD, RADIUSRULE, RADIUSSIZE, RADIUSTYPE,
        POLARIZATION, VDWTYPE
    }

    public enum ForceFieldDouble {

        A_AXIS, B_AXIS, C_AXIS, ALPHA, BETA, GAMMA, ANGLE_CUBIC, ANGLE_QUARTIC,
        ANGLE_PENTIC, ANGLE_SEXTIC, BOND_CUBIC, BOND_QUARTIC, OPBENDUNIT,
        TORSIONUNIT, DIELECTRIC, POLAR_DAMP, POLAR_SOR, POLAR_EPS, EWALD_CUTOFF,
        EWALD_ALPHA, PME_SPACING, VDW_13_SCALE, VDW_14_SCALE, VDW_15_SCALE,
        VDW_CUTOFF, MPOLE_11_SCALE, MPOLE_12_SCALE, MPOLE_13_SCALE, MPOLE_14_SCALE,
        MPOLE_15_SCALE, POLAR_11_SCALE, POLAR_12_SCALE, POLAR_13_SCALE,
        POLAR_14_SCALE, POLAR_15_SCALE, DIRECT_11_SCALE, DIRECT_12_SCALE,
        DIRECT_13_SCALE, DIRECT_14_SCALE, MUTUAL_11_SCALE, MUTUAL_12_SCALE,
        MUTUAL_13_SCALE, MUTUAL_14_SCALE
    }

    public enum ForceFieldInteger {

        PME_ORDER
    }

    public enum ForceFieldBoolean {

        BONDTERM, ANGLETERM, STRBNDTERM, UREYTERM, OPBENDTERM,
        TORSIONTERM, PITORSTERM, TORTORTERM, VDWTERM,
        MPOLETERM, POLARIZETERM, SCFCACHE;
    }

    public enum ForceFieldType {

        ATOM, ANGLE, BIOTYPE, BOND, CHARGE, MULTIPOLE, OPBEND, PITORS, POLARIZE,
        STRBND, TORSION, TORTORS, UREYBRAD, VDW, KEYWORD
    }
    public File forceFieldFile;
    public File keywordFile;
    private final TreeMap<String, AngleType> angleTypes;
    private final TreeMap<String, AtomType> atomTypes;
    private final TreeMap<String, BondType> bondTypes;
    private final TreeMap<String, BioType> bioTypes;
    private final TreeMap<String, ChargeType> chargeTypes;
    private final TreeMap<String, MultipoleType> multipoleTypes;
    private final TreeMap<String, OutOfPlaneBendType> outOfPlaneBendTypes;
    private final TreeMap<String, PolarizeType> polarizeTypes;
    private final TreeMap<String, StretchBendType> stretchBendTypes;
    private final TreeMap<String, PiTorsionType> piTorsionTypes;
    private final TreeMap<String, TorsionType> torsionTypes;
    private final TreeMap<String, TorsionTorsionType> torsionTorsionTypes;
    private final TreeMap<String, UreyBradleyType> ureyBradleyTypes;
    private final TreeMap<String, VDWType> vanderWaalsTypes;
    private final TreeMap<String, Keyword> keywordTypes;
    private final TreeMap<ForceFieldType, TreeMap> forceFieldTypes;

    /**
     * ForceField Constructor.
     */
    public ForceField(File forceFieldFile, File keyFile) {
        this.forceFieldFile = forceFieldFile;
        this.keywordFile = keyFile;
        angleTypes = new TreeMap<String, AngleType>();
        atomTypes = new TreeMap<String, AtomType>();
        bondTypes = new TreeMap<String, BondType>();
        bioTypes = new TreeMap<String, BioType>();
        chargeTypes = new TreeMap<String, ChargeType>();
        outOfPlaneBendTypes = new TreeMap<String, OutOfPlaneBendType>();
        multipoleTypes = new TreeMap<String, MultipoleType>();
        piTorsionTypes = new TreeMap<String, PiTorsionType>();
        polarizeTypes = new TreeMap<String, PolarizeType>();
        stretchBendTypes = new TreeMap<String, StretchBendType>();
        torsionTypes = new TreeMap<String, TorsionType>();
        torsionTorsionTypes = new TreeMap<String, TorsionTorsionType>();
        ureyBradleyTypes = new TreeMap<String, UreyBradleyType>();
        vanderWaalsTypes = new TreeMap<String, VDWType>();
        keywordTypes = new TreeMap<String, Keyword>();

        forceFieldTypes = new TreeMap<ForceFieldType, TreeMap>();
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
        forceFieldTypes.put(ForceFieldType.TORTORS, torsionTorsionTypes);
        forceFieldTypes.put(ForceFieldType.UREYBRAD, ureyBradleyTypes);
        forceFieldTypes.put(ForceFieldType.VDW, vanderWaalsTypes);
        forceFieldTypes.put(ForceFieldType.KEYWORD, keywordTypes);
    }

    public double getDouble(ForceFieldDouble forceFieldDouble)
            throws Exception {
        Keyword keyword = keywordTypes.get(forceFieldDouble.toString());
        if (keyword == null) {
            throw new Exception("Undefined keyword: " + forceFieldDouble.toString());
        }
        return Double.parseDouble(keyword.getEntry(0));
    }

    public double getDouble(ForceFieldDouble forceFieldDouble,
            double defaultValue) {
        try {
            return getDouble(forceFieldDouble);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    public int getInteger(ForceFieldInteger forceFieldInteger)
            throws Exception {
        Keyword keyword = keywordTypes.get(forceFieldInteger.toString());
        if (keyword == null) {
            throw new Exception("Undefined keyword: " + forceFieldInteger.toString());
        }
        return Integer.parseInt(keyword.getEntry(0));
    }

    public int getInteger(ForceFieldInteger forceFieldInteger,
            int defaultValue) {
        try {
            return getInteger(forceFieldInteger);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    public String getString(ForceFieldString forceFieldString)
            throws Exception {
        Keyword keyword = keywordTypes.get(forceFieldString.toString());
        if (keyword == null) {
            throw new Exception("Undefined keyword: " + forceFieldString.toString());
        }
        return keyword.getEntry(0);
    }

    public String getString(ForceFieldString forceFieldString,
            String defaultString) {
        try {
            return getString(forceFieldString);
        } catch (Exception e) {
            return defaultString;
        }
    }

    public boolean getBoolean(ForceFieldBoolean forceFieldBoolean)
            throws Exception {
        Keyword keyword = keywordTypes.get(forceFieldBoolean.toString());
        if (keyword == null) {
            throw new Exception("Undefined keyword: " + forceFieldBoolean.toString());
        }
        return Boolean.parseBoolean(keyword.getEntry(0));
    }

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
     * @param forceFieldDouble
     *            ForceFieldDouble
     * @param value
     *            double
     */
    public void addForceFieldDouble(ForceFieldDouble forceFieldDouble, double value) {
        Keyword keyword = keywordTypes.get(forceFieldDouble.toString());
        if (keyword != null) {
            double oldValue = getDouble(forceFieldDouble, value);
            if (Double.compare(oldValue, value) != 0) {
                log.warning(String.format(
                        "\nOld %s value of %6.3e replaced with %6.3e\n", keyword,
                        oldValue, value));
            }
            keyword.clear();
            keyword.append(Double.toString(value));
        } else {
            keyword = new Keyword(forceFieldDouble.toString(), Double.toString(value));
            keywordTypes.put(forceFieldDouble.toString(), keyword);
        }
    }

    /**
     * Add a force field keyword that is represented by an int.
     *
     * @param forceFieldInteger
     *            ForceFieldInteger
     * @param value
     *            int
     */
    public void addForceFieldInteger(ForceFieldInteger forceFieldInteger, int value) {
        Keyword keyword = keywordTypes.get(forceFieldInteger.toString());
        if (keyword != null) {
            int oldValue = getInteger(forceFieldInteger, value);
            if (oldValue != value) {
                log.warning(String.format(
                        "\nOld %s value of %6d replaced with %6d\n", keyword,
                        oldValue, value));
            }
            keyword.clear();
            keyword.append(Integer.toString(value));
        } else {
            keyword = new Keyword(forceFieldInteger.toString(), Integer.toString(value));
            keywordTypes.put(forceFieldInteger.toString(), keyword);
        }
    }

    /**
     * Store a force field keyword that is represented by a String.
     *
     * @param forceFieldString
     *            ForceFieldString
     * @param value
     *            String
     */
    public void addForceFieldString(ForceFieldString forceFieldString, String value) {
        Keyword keyword = keywordTypes.get(forceFieldString.toString());
        if (keyword != null) {
            String oldValue = getString(forceFieldString, value);
            if (!oldValue.equalsIgnoreCase(value)) {
                log.warning(String.format(
                        "\nOld %s value of %6d replaced with %6d\n", forceFieldString,
                        oldValue, value));
            }
            keyword.clear();
            keyword.append(value);
        } else {
            keyword = new Keyword(forceFieldString.toString(), value);
            keywordTypes.put(forceFieldString.toString(), keyword);
        }
    }

    /**
     * Store a force field keyword that is represented by a Boolean.
     *
     * @param forceFieldBoolean ForceFielBoolean
     *
     * @param value Boolean
     */
    public void addForceFieldBoolean(ForceFieldBoolean forceFieldBoolean, Boolean value) {
        Keyword keyword = keywordTypes.get(forceFieldBoolean.toString());
        if (keyword != null) {
            boolean oldValue = getBoolean(forceFieldBoolean, value);
            if (oldValue != value) {
                log.warning(String.format(
                        "\nOld %s value of %l replaced with %l\n", keyword,
                        oldValue, value));
            }
            keyword.clear();
            keyword.append(Boolean.toString(value));
        } else {
            keyword = new Keyword(forceFieldBoolean.toString(), Boolean.toString(value));
            keywordTypes.put(forceFieldBoolean.toString(), keyword);
        }
    }

    /**
     * Add an instance of a force field type. Force Field types are more
     * complicated than simple Strings or doubles, in that they have multiple
     * fields and may occur multiple times.
     *
     * @param type
     *            BaseType
     */
    public void addForceFieldType(BaseType type) {
        if (type == null) {
            log.info("Null force field type ignored.");
            return;
        }
        TreeMap treeMap = forceFieldTypes.get(type.forceFieldType);
        if (treeMap == null) {
            log.info("Unrecognized force field type ignored " + type.forceFieldType);
            type.print();
        }
        if (treeMap.containsKey(type.key)) {
            log.fine("A force field entry of type " + type.forceFieldType
                    + " already exists with the key: " + type.key
                    + "\nThe (discarded) old entry:\n" + treeMap.get(type.key).
                    toString() + "\nThe new entry:\n" + type.toString());
        }
        Class baseTypeClass = type.getClass();
        treeMap.put(type.key, baseTypeClass.cast(type));
    }

    public AngleType getAngleType(String key) {
        return angleTypes.get(key);
    }

    public AtomType getAtomType(String key) {
        return atomTypes.get(key);
    }

    public BondType getBondType(String key) {
        return bondTypes.get(key);
    }

    public BioType getBioType(String key) {
        return bioTypes.get(key);
    }

    public ChargeType getChargeType(String key) {
        return chargeTypes.get(key);
    }

    public MultipoleType getMultipoleType(String key) {
        return multipoleTypes.get(key);
    }

    public OutOfPlaneBendType getOutOfPlaneBendType(String key) {
        return outOfPlaneBendTypes.get(key);
    }

    public PolarizeType getPolarizeType(String key) {
        return polarizeTypes.get(key);
    }

    public StretchBendType getStretchBendType(String key) {
        return stretchBendTypes.get(key);
    }

    public PiTorsionType getPiTorsionType(String key) {
        return piTorsionTypes.get(key);
    }

    public TorsionType getTorsionType(String key) {
        return torsionTypes.get(key);
    }

    public TorsionTorsionType getTorsionTorsionType(String key) {
        return torsionTorsionTypes.get(key);
    }

    public UreyBradleyType getUreyBradleyType(String key) {
        return ureyBradleyTypes.get(key);
    }

    public VDWType getVDWType(String key) {
        return vanderWaalsTypes.get(key);
    }

    public TreeMap<String, VDWType> getVDWTypes() {
        return vanderWaalsTypes;
    }

    public int getForceFieldTypeCount(ForceFieldType type) {
        TreeMap<String, BaseType> treeMap =
                (TreeMap<String, BaseType>) forceFieldTypes.get(type);
        if (treeMap == null) {
            log.warning("Unrecognized Force Field Type: " + type);
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
                    log.severe("Polarize type " + key
                            + "references nonexistant polarize type " + key2);
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

    public void print() {
        for (ForceFieldType s : forceFieldTypes.keySet()) {
            print(s.toString());
        }
    }

    /**
     * Prints any force field keyword to Standard.out.
     *
     * @param key
     *            String
     */
    public void print(String key) {
        ForceFieldType type = ForceFieldType.valueOf(key);
        log.info(toString(type));
    }

    /**
     * Return a String for any Force Field keyword.
     *
     * @param type
     *            ForceFieldType
     * @return String
     */
    public String toString(ForceFieldType type) {
        StringBuffer stringBuffer = new StringBuffer("\n");
        TreeMap t = forceFieldTypes.get(type);
        if (t.size() > 0) {
            for (Object o : t.values()) {
                stringBuffer.append(o + "\n");
            }
            return stringBuffer.toString();
        }
        return null;
    }

    public String toString(String key) {
        if (key == null) return null;
        Keyword keyword = keywordTypes.get(key.toUpperCase());
        if (keyword != null) {
            return keyword.toString();
        } else {
            return null;
        }
    }
}
