//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.potential.extended;

import javax.swing.tree.TreeNode;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.Properties;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.AbstractConfiguration;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.lang3.StringUtils;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import static ffx.potential.extended.ExtConstants.RNG;
import static ffx.utilities.Constants.kB;

/**
 * <p>ExtUtils class.</p>
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public final class ExtUtils {

    /**
     * Static class.
     */
    private ExtUtils() {
    }

    private static final Logger logger = Logger.getLogger(ExtUtils.class.getName());

    /**
     * <p>setProp.</p>
     *
     * @param key a {@link java.lang.String} object.
     * @param set a T object.
     * @param <T> a T object.
     */
    public static <T> void setProp(String key, T set) {
        if (key.startsWith("esv.")) {
            if (Arrays.asList(ExtendedSystem.ExtendedSystemConfig.class.getDeclaredFields())
                    .stream().noneMatch((Field f) -> key.substring(4).equalsIgnoreCase(f.getName()))) {
                logger.warning(format("Setting an unused property: %s", key));
            }
        }
        System.setProperty(key, String.valueOf(set));
    }

    /**
     * <p>setProp.</p>
     *
     * @param set  a boolean.
     * @param keys a {@link java.lang.String} object.
     */
    public static void setProp(boolean set, String... keys) {
        for (String key : keys) {
            setProp(key, set);
        }
    }

    /**
     * Constant <code>debug=prop("dbg.debug", false)</code>
     */
    public static final boolean debug = prop("dbg.debug", false);
    /**
     * Constant <code>verbose=prop("dbg.verbose", false)</code>
     */
    public static final boolean verbose = prop("dbg.verbose", false);
    /**
     * Constant <code>doubleFormat="prop(dbg.debugFormat, %.3f)"</code>
     */
    public static final String doubleFormat = prop("dbg.debugFormat", "%.3f");

    /**
     * Print the Properties defined by the given key prefixes.
     *
     * @param header      a {@link java.lang.String} object.
     * @param properties  a {@link java.util.Properties} object.
     * @param keyPrefixes an array of {@link java.lang.String} objects.
     */
    public static void printConfigSet(String header, Properties properties, String[] keyPrefixes) {
        StringBuilder SB = new StringBuilder();
        properties.keySet().stream()
                .filter(k -> {
                    String key = k.toString().toLowerCase();
                    for (String prefix : keyPrefixes) {
                        if (key.startsWith(prefix)) {
                            return true;
                        }
                    }
                    return false;
                })
                .forEach(key -> SB.append(String.format("   %-30s %s",
                        key.toString() + ":", System.getProperty(key.toString()))));
        if (!SB.toString().isEmpty()) {
            if (header != null) {
                SB.insert(0, format(" %s", header));
            }
            logger.info(SB.toString());
        }
    }

    /**
     * <p>printConfigSet.</p>
     *
     * @param header     a {@link java.lang.String} object.
     * @param properties a {@link java.util.Properties} object.
     * @param prefix     a {@link java.lang.String} object.
     */
    public static void printConfigSet(String header, Properties properties, String prefix) {
        printConfigSet(header, properties, new String[]{prefix});
    }

    /**
     * <p>primitivizeIntegers.</p>
     *
     * @param list a {@link java.util.List} object.
     * @return an array of {@link int} objects.
     */
    public static int[] primitivizeIntegers(List<Integer> list) {
        return list.stream().mapToInt((Integer i) -> i).toArray();
    }

    /**
     * <p>primitivizeDoubles.</p>
     *
     * @param list a {@link java.util.List} object.
     * @return an array of {@link double} objects.
     */
    public static double[] primitivizeDoubles(List<Double> list) {
        return list.stream().mapToDouble((Double i) -> i).toArray();
    }

    /**
     * <p>initializeBackgroundMultipoles.</p>
     *
     * @param backgroundAtoms a {@link java.util.List} object.
     * @param ff              a {@link ffx.potential.parameters.ForceField} object.
     */
    public static void initializeBackgroundMultipoles(List<Atom> backgroundAtoms, ForceField ff) {
        for (int i = 0; i < backgroundAtoms.size(); i++) {
            Atom bg = backgroundAtoms.get(i);
            MultipoleType type = MultipoleType.multipoleTypeFactory(bg, ff);
            if (type == null) {
                logger.severe(format("No multipole could be assigned to atom %s of type %s.",
                        bg.toString(), bg.getAtomType()));
            }
        }
    }

    @SafeVarargs
    @SuppressWarnings("serial")
    private static <T> Collection<T> join(Collection<? extends T>... add) {
        return (new ArrayList<T>() {
            {
                for (Collection<? extends T> c : add) {
                    addAll(c);
                }
            }
        });
    }

    /**
     * <p>join.</p>
     *
     * @param add a {@link java.util.List} object.
     * @param <T> a T object.
     * @return a {@link java.util.List} object.
     */
    @SafeVarargs
    @SuppressWarnings("serial")
    public static <T> List<T> join(List<? extends T>... add) {
        return (new ArrayList<T>() {
            {
                for (Collection<? extends T> c : add) {
                    addAll(c);
                }
            }
        });
    }

    /**
     * <p>intRange.</p>
     *
     * @param start a int.
     * @param end   a int.
     * @return a {@link java.util.List} object.
     */
    public static List<Integer> intRange(int start, int end) {
        List<Integer> intArray = new ArrayList<>();
        for (int i = start; i <= end; i++) {
            intArray.add(i);
        }
        return intArray;
    }

    /**
     * <p>intRangeArray.</p>
     *
     * @param start a int.
     * @param end   a int.
     * @return an array of {@link int} objects.
     */
    public static int[] intRangeArray(int start, int end) {
        int array[] = new int[end - start + 1];
        for (int i = start; i <= end; i++) {
            array[i - start] = i;
        }
        return array;
    }

    /**
     * <p>intRange.</p>
     *
     * @param single a int.
     * @return a {@link java.util.List} object.
     */
    public static List<Integer> intRange(int single) {
        return intRange(single, single);
    }

    /**
     * <p>view.</p>
     *
     * @param c   a {@link java.util.Collection} object.
     * @param <T> a T object.
     * @return a {@link java.util.Collection} object.
     */
    public static <T> Collection<T> view(Collection<? extends T> c) {
        return (c != null) ? Collections.unmodifiableCollection(c) : null;
    }

    /**
     * <p>view.</p>
     *
     * @param c   a {@link java.util.List} object.
     * @param <T> a T object.
     * @return a {@link java.util.List} object.
     */
    public static <T> List<T> view(List<? extends T> c) {
        return (c != null) ? Collections.unmodifiableList(c) : null;
    }

    /**
     * <p>arrayCopy.</p>
     *
     * @param from an array of T[] objects.
     * @param <T>  a T object.
     * @return an array of T[] objects.
     */
    @SuppressWarnings("unchecked")
    public static <T> T[] arrayCopy(T[] from) {
        List<T> to = new ArrayList<>();
        to.addAll(Arrays.asList(from));
        return (T[]) to.toArray();
    }

    /**
     * <p>frmt.</p>
     *
     * @param x an array of {@link double} objects.
     * @return a {@link java.lang.String} object.
     */
    public static String frmt(double[] x) {
        return formatArray(x);
    }

    /**
     * Helper method for logging distance and multipole arrays.
     *
     * @param x an array of {@link double} objects.
     * @return a {@link java.lang.String} object.
     */
    public static String formatArray(double[] x) {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < x.length; i++) {
            sb.append(format(doubleFormat, x[i]));
            if (i + 1 < x.length) {
                sb.append(", ");
            }
        }
        sb.append("]");
        return sb.toString();
    }

    /**
     * Helper method for logging distance and multipole arrays.
     *
     * @param x an array of {@link double} objects.
     * @return a {@link java.lang.String} object.
     */
    public static String formatArray(double[][] x) {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < x.length; i++) {
            sb.append("[");
            for (int j = 0; j < x[i].length; j++) {
                sb.append(format(doubleFormat, x[i][j]));
                if (j + 1 < x[i].length) {
                    sb.append(", ");
                }
            }
            sb.append("]");
            if (i + 1 < x.length) {
                sb.append("; ");
            }
        }
        sb.append("]");
        return sb.toString();
    }

    /**
     * Parse configuration properties of type String, Boolean, Integer, Double,
     * OptionalInt, OptionalDouble, and any list or array type. Default value
     * determines the return type. Provides case-insensitivity in property keys.
     *
     * @param key        a {@link java.lang.String} object.
     * @param defaultVal a T object.
     * @param properties a {@link org.apache.commons.configuration2.AbstractConfiguration} object.
     * @param <T>        a T object.
     * @return a T object.
     * @throws java.util.NoSuchElementException if any.
     * @throws java.lang.NumberFormatException  if any.
     * @throws java.util.NoSuchElementException if any.
     * @throws java.lang.NumberFormatException  if any.
     */
    @SuppressWarnings({"unchecked"})
    public static final <T> T prop(String key, final T defaultVal, final AbstractConfiguration properties)
            throws NoSuchElementException, NumberFormatException {
        Objects.requireNonNull(key);
        Objects.requireNonNull(defaultVal);
        Objects.requireNonNull(properties);
        try {
            // Determine whether the property exists in any valid form.
            if (!properties.containsKey(key)) {
                key = key.toLowerCase();
                if (properties.containsKey(key.replaceAll("_", "-"))) {
                    key = key.replaceAll("_", "-");
                } else if (properties.containsKey(key.replaceAll("_", ".").replaceAll("-", "."))) {
                    key = key.replaceAll("_", ".").replaceAll("-", ".");
                } else if (properties.containsKey(key.toUpperCase())) {
                    key = key.toUpperCase();
                } else {
                    boolean found = false;
                    Iterator<String> keys = properties.getKeys();
                    while (keys.hasNext()) {
                        String search = keys.next();
                        if (search.equalsIgnoreCase(key)) {
                            key = search;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        return defaultVal;
                    }
                }
            }
            if (properties.getString(key) == null) {
                return defaultVal;
            }
            // The property exists and is non-null.
            if (defaultVal instanceof String) {
                return (T) properties.getString(key, (String) defaultVal);
            } else if (defaultVal instanceof Integer) {
                return (T) properties.getInteger(key, (Integer) defaultVal);
            } else if (defaultVal instanceof OptionalInt) {
                return (T) OptionalInt.of(properties.getInt(key, (Integer) defaultVal));
            } else if (defaultVal instanceof Double) {
                return (T) properties.getDouble(key, (Double) defaultVal);
            } else if (defaultVal instanceof OptionalDouble) {
                return (T) OptionalDouble.of(properties.getDouble(key, (Double) defaultVal));
            } else if (defaultVal instanceof Boolean) {
                // Infer true for (only) explicitly set boolean properties, eg. -Dgkterm.
                if (properties.getString(key).isEmpty()) {
                    properties.setProperty(key, "true");
                }
                return (T) properties.getBoolean(key, (Boolean) defaultVal);
            } else if (defaultVal instanceof List<?>) {
                List<?> defaultVals = (List<?>) defaultVal;
                List<Object> objs = properties.getList(key, defaultVals);
                if (defaultVals.get(0) instanceof Double) {
                    List<Double> ret = new ArrayList<>();
                    for (Object obj : objs) {
                        ret.add(Double.parseDouble(obj.toString().replaceAll("[\\[\\],]", "")));
                    }
                    return (T) ret;
                } else if (defaultVals.get(0) instanceof Integer) {
                    List<Integer> ret = new ArrayList<>();
                    for (Object obj : objs) {
                        ret.add(Integer.parseInt(obj.toString().replaceAll("[\\[\\],]", "")));
                    }
                    return (T) ret;
                }
                throw new IllegalArgumentException();
            } else if (defaultVal.getClass().isArray()) {
                return (T) properties.getList(key, Arrays.asList(defaultVal)).toArray(new Object[0]);
            } else {
                throw new IllegalArgumentException();
            }
        } catch (RuntimeException ex) {
            logger.warning(String.format("Error parsing property %s with value %s; the default is %s.",
                    key, properties.getString(key, "null"), defaultVal.toString()));

            logger.info(ex.toString());
            ex.printStackTrace();

            throw ex;
        }
    }

    private static AbstractConfiguration systemCache;

    /**
     * If no Configuration is supplied, default to a converted and cached copy
     * of System.getProperties().
     *
     * @param key        a {@link java.lang.String} object.
     * @param defaultVal a T object.
     * @param <T>        a T object.
     * @return a T object.
     */
    public static final <T> T prop(String key, T defaultVal) {
        if (systemCache == null) {
            systemCache = new CompositeConfiguration();
        }
        AbstractConfiguration systemCache = new CompositeConfiguration();
        for (Object sysKey : System.getProperties().keySet()) {
            systemCache.addProperty((String) sysKey, System.getProperty((String) sysKey));
        }
        return prop(key, defaultVal, systemCache);
    }

    /**
     * Parse property but warn if the default is not retured.
     *
     * @param key        a {@link java.lang.String} object.
     * @param defaultVal a T object.
     * @param warning    a {@link java.lang.String} object.
     * @param <T>        a T object.
     * @return a T object.
     */
    public static final <T> T prop(String key, T defaultVal, String warning) {
        T parsed = prop(key, defaultVal);
        if (!parsed.equals(defaultVal) && warning != null) {
            logger.warning(warning);
        }
        return parsed;
    }

    /**
     * Parse system properties into an Enum class.
     *
     * @param <T>  the Class to operate on.
     * @param type a {@link java.lang.Class} object.
     * @param key  a {@link java.lang.String} object.
     * @param def  a T instance.
     * @return a T object.
     * @throws java.lang.IllegalArgumentException if any.
     */
    public static final <T extends Enum<T>> T prop(Class<T> type, String key, T def)
            throws IllegalArgumentException {
        String value = System.getProperty(key);
        if (value == null) {
            return def;
        }
        try {
            T parsed = T.valueOf(type, value);
            return parsed;
        } catch (IllegalArgumentException ex) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("Invalid property definition: %s = %s for type %s.", key, value, type.getName()));
            sb.append(format("  Allowable values:  %s", Arrays.toString(type.getEnumConstants())));
            sb.append(format("  Returning default: %s", def.toString()));
            logger.warning(sb.toString());
            return def;
        }
    }

    /**
     * <p>prop.</p>
     *
     * @param <T>  the Class to operate on.
     * @param key  a {@link java.lang.String} object.
     * @param type a {@link java.lang.Class} object.
     * @param def  a T instance.
     * @return a T object.
     */
    public static final <T extends Enum<T>> T prop(String key, Class<T> type, T def) {
        return prop(type, key, def);
    }

    /**
     * Return velocities from a Maxwell-Boltzmann distribution of momenta. The
     * variance of each independent momentum component is kT * mass.
     *
     * @param mass        a double.
     * @param temperature a double.
     * @return an array of {@link double} objects.
     */
    public static double[] maxwellVelocity(double mass, double temperature) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = RNG.nextGaussian() * sqrt(kB * temperature / mass);
        }
        return vv;
    }

    /**
     * <p>printTreeFromNode.</p>
     *
     * @param target a {@link ffx.potential.bonded.MSNode} object.
     */
    public static void printTreeFromNode(MSNode target) {
        StringBuilder sb = new StringBuilder();
        sb.append(format(" Tree Hierarchy for %s\n", target.getName()));
        sb.append(format(" ===================%s\n", StringUtils.repeat("=", target.getName().length())));
        TreeNode root = target;
        while (root.getParent() != null) {
            root = root.getParent();
        }
        List<MSNode> toLeaves = target.getDescendants(MSNode.class);
        boolean init = false;
        for (MSNode node : toLeaves) {
            if (!init) {
                sb.append(format(" %s> %s\n", StringUtils.repeat(">", node.getLevel()), node.getName()));
                init = true;
                continue;
            }
            if (!node.getName().isEmpty()) {
                sb.append(format(" └%s> %s\n", StringUtils.repeat("-", node.getLevel()), node.getName()));
            }
        }
        logger.info(sb.toString());
    }

    /**
     * Element-wise sum over a list of 1D double arrays. This implementation
     * benchmarks faster than the equivalent Java8 stream() API.
     */
    private static double[] eleSum1DArrays(List<double[]> terms, int numESVs) {
        double[] termSum = new double[terms.size()];
        for (int iTerm = 0; iTerm < terms.size(); iTerm++) {
            double[] currentTerm = terms.get(iTerm);
            if (currentTerm.length != numESVs) {
                logger.warning(format("iTerm %d length: %d, numESVs: %d", iTerm, terms.get(iTerm).length, numESVs));
                throw new IndexOutOfBoundsException();
            }
            for (int iESV = 0; iESV < numESVs; iESV++) {
                termSum[iESV] += currentTerm[iESV];
            }
        }
        return termSum;
    }

    /**
     * Element-wise sum over a list of 2D double arrays.
     */
    private static double[][] eleSum2DArrays(List<double[][]> terms, int numESVs, int nVars) {
        if (terms == null || terms.isEmpty()) {
            throw new NullPointerException("Summing an empty or null terms list.");
        }
        double[][] termSum = new double[numESVs][nVars];
        for (int iTerm = 0; iTerm < terms.size(); iTerm++) {
            double[][] currentTerm = terms.get(iTerm);
            if (currentTerm.length != numESVs) {
                throw new IndexOutOfBoundsException();
            }
            for (int iESV = 0; iESV < numESVs; iESV++) {
                if (currentTerm[iESV].length != nVars) {
                    throw new IndexOutOfBoundsException();
                }
                for (int iAtom = 0; iAtom < nVars; iAtom++) {
                    termSum[iESV][iAtom] += currentTerm[iESV][iAtom];
                }
            }
        }
        return termSum;
    }

    /**
     * <p>removeFinalNewline.</p>
     *
     * @param sb a {@link java.lang.StringBuilder} object.
     * @return a {@link java.lang.StringBuilder} object.
     */
    public static StringBuilder removeFinalNewline(StringBuilder sb) {
        if (sb.lastIndexOf("\n") + 1 == sb.length()) {
            sb.replace(sb.length() - 1, sb.length(), "");
        }
        return sb;
    }

    /**
     * Returns the result of calling toString() on each element of the array.
     * Returned type is Object[] (rather than String[]) to support use as
     * varargs. Example: format("%s", Arrays.toString( array ); // GOOD
     * format("%s %s %s %s", Arrays.toString( array ); // FAIL
     * <p>
     * format("%s", arrayToStrings( array ); // FAIL format("%s %s %s %s",
     * arrayToStrings( array ); // GOOD
     *
     * @param args an array of {@link java.lang.Object} objects.
     * @return an array of {@link java.lang.Object} objects.
     */
    public static Object[] arrayToStrings(Object[] args) {
        return Arrays.stream(args).map(Object::toString).toArray();
    }

    /**
     * Use this to iterate easily over a collection of lists. Returns the
     * inferred parent list type. Duplicate elements allowed. "View" intentional
     * to avoid vain attempts at modifying the component lists.
     *
     * @param lists a {@link java.util.List} object.
     * @param <T>   a T object.
     * @return a {@link java.util.List} object.
     */
    @SafeVarargs
    public static final <T> List<T> joinedListView(List<? extends T>... lists) {
        List<T> joined = new ArrayList<>();
        for (List<? extends T> list : lists) {
            joined.addAll(list);
        }
        return Collections.unmodifiableList(joined);
    }

    @SuppressWarnings("serial")
    public static class UnderConstructionException extends RuntimeException {
    }

}
