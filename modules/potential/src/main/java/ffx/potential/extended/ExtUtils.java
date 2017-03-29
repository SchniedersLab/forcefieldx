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
package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;

import javax.swing.tree.TreeNode;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.lang.StringUtils;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.nonbonded.Multipole;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ExtConstants.RNG;
import static ffx.potential.extended.ExtConstants.kB;
import static ffx.potential.extended.ExtConstants.roomTemperature;
import static ffx.potential.extended.ExtUtils.DebugHandler.debugFormat;
import static ffx.potential.extended.SBLogger.SB;

/**
 * 
 * @author slucore
 */
public final class ExtUtils {
    private ExtUtils() {}   // static class
    private static final Logger logger = Logger.getLogger(ExtUtils.class.getName());
	
	/**
	 * Print the Properties defined by the given key prefixes.
	 */
	public static void printConfigSet(String header, Properties properties, String[] keyPrefixes) {
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
				.forEach(key -> SB.logfn("   %-20s %s",
						key.toString() + ":", System.getProperty(key.toString())));
		if (!SB.isEmpty()) {
			if (header != null) {
				SB.headern(" %s", header);
			}
			SB.print();
		}
	}
	public static void printConfigSet(String header, Properties properties, String prefix) {
		printConfigSet(header, properties, new String[]{prefix});
	}
    
    public static void initializeBackgroundMultipoles(List<Atom> backgroundAtoms, ForceField ff) {
        for (int i = 0; i < backgroundAtoms.size(); i++) {
            Atom bg = backgroundAtoms.get(i);
            Multipole multipole = Multipole.buildMultipole(bg, ff);
            if (multipole == null) {
                logger.severe(format("No multipole could be assigned to atom %s of type %s.",
                        bg.toString(), bg.getAtomType()));
            }
        }
    }
    
    /**
     * Helper method for logging distance and multipole arrays.
     */
    public static String formatArray(double[] x) {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < x.length; i++) {
            sb.append(format(debugFormat(), x[i]));
            if (i + 1 < x.length) {
                sb.append(", ");
            }
        }
        sb.append("]");
        return sb.toString();
    }
    
    /**
     * Helper method for logging distance and multipole arrays.
     */
    public static String formatArray(double[][] x) {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < x.length; i++) {
            sb.append("[");
            for (int j = 0; j < x[i].length; j++) {
                sb.append(format(debugFormat(), x[i][j]));
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
     * Parse system properties of type String, Boolean, Integer, Double, OptionalInt, and OptionalDouble.
     * Default value determines the return type. Provides case-insensitivity in property keys.
     * TODO: adapt this to use a CompositeConfiguration object rather than the unsafe System.getProperties().
     */
    @SuppressWarnings("unchecked")
    public static final <T> T prop(String key, T defaultVal)
            throws NoSuchElementException, NumberFormatException {
        if (defaultVal == null) {
            throw new IllegalArgumentException();
        }
        T parsed = null;
        try {
            if (System.getProperty(key) == null) {
                if (System.getProperty(key.toLowerCase()) != null) {
                    key = key.toLowerCase();
                } else if (System.getProperty(key.toUpperCase()) != null) {
                    key = key.toUpperCase();
                } else {
                    boolean found = false;
                    for (Object search : System.getProperties().keySet()) {
                        if (search instanceof String && ((String) search).equalsIgnoreCase(key)) {
                            key = (String) search;
                        }
                    }
                    if (!found) {
                        return defaultVal;
                    }
                }
            }
            if (defaultVal instanceof String) {
                parsed = (System.getProperty(key) != null)
                        ? (T) System.getProperty(key) : defaultVal;
            } else if (defaultVal instanceof Integer) {
                parsed = (System.getProperty(key) != null)
                        ? (T) Integer.valueOf(System.getProperty(key)) : defaultVal;
            } else if (defaultVal instanceof OptionalInt) {
                parsed = (System.getProperty(key) != null)
                        ? (T) OptionalInt.of(Integer.parseInt(System.getProperty(key))) : defaultVal;
            } else if (defaultVal instanceof Double) {
                parsed = (System.getProperty(key) != null)
                        ? (T) Double.valueOf(System.getProperty(key)) : defaultVal;
            } else if (defaultVal instanceof OptionalDouble) {
                parsed = (System.getProperty(key) != null)
                        ? (T) OptionalDouble.of(Double.parseDouble(System.getProperty(key))) : defaultVal;
            } else if (defaultVal instanceof Boolean) {
                if (System.getProperty(key) != null) {
                    if (System.getProperty(key).isEmpty()) {
                        System.setProperty(key, "true");
                    }
                    parsed = (T) Boolean.valueOf(System.getProperty(key));
                }
            } else {
                throw new IllegalArgumentException();
            }
        } catch (IllegalArgumentException ex) {
            String value = (System.getProperty(key) != null) ? System.getProperty(key) : "null";
            logger.warning(String.format("Error parsing property %s with value %s; the default is %s.",
                    key, value, defaultVal.toString()));
            throw ex;
        }
        return parsed;
    }
    
    /**
     * Parse property but warn if the default is not used.
     */
    public static final <T> T prop(String key, T defaultVal, String warning) {
        T parsed = prop(key, defaultVal);
        if (!parsed.equals(defaultVal) && warning != null) {
            SB.warning(warning);
        }
        return parsed;
    }
    
    /**
     * Parse system properties into an Enum class.
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
            SB.nlogfn("Invalid property definition: %s = %s for type %s.", key, value, type.getName());
            SB.logfn( "  Allowable values:  %s", Arrays.toString(type.getEnumConstants()));
            SB.logf(  "  Returning default: %s", def.toString());
            SB.warning();
            return def;
        }
    }
    
    /**
     * Return room-temperature velocities from a Maxwell-Boltzmann distribution of momenta.
     */
    public static double[] singleRoomtempMaxwell(double mass) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = RNG.nextGaussian() * sqrt(kB * roomTemperature / mass);
        }
        return vv;
    }

    /**
     * Return velocities from a Maxwell-Boltzmann distribution of momenta.
     * The variance of each independent momentum component is kT * mass.
     */
    public static double[] singleMaxwell(double mass, double temperature) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = RNG.nextGaussian() * sqrt(kB * temperature / mass);
        }
        return vv;
    }
    
    public static void printTreeFromNode(MSNode target) {
        StringBuilder sb = new StringBuilder();
        sb.append(format(" Tree Hierarchy for %s\n", target.getName()));
        sb.append(format(" ===================%s\n", StringUtils.repeat("=", target.getName().length())));
        TreeNode root = target;
        while (root.getParent() != null) {
            root = root.getParent();
        }
        Enumeration<?> fromRoot = target.pathFromAncestorEnumeration(root);
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
     * Element-wise sum over a list of 1D double arrays.
     * This implementation benchmarks faster than the equivalent Java8 stream() API.
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
    
    public static StringBuilder removeFinalNewline(StringBuilder sb) {
        if (sb.lastIndexOf("\n") + 1 == sb.length()) {
            sb.replace(sb.length() - 1, sb.length(), "");
        }
        return sb;
    }
    
    /**
     * Returns the result of calling toString() on each element of the array.
     *  Returned type is Object[] (rather than String[]) to support use as varargs.
     *  Example: 
     *      format("%s",          Arrays.toString( array );  // GOOD
     *      format("%s %s %s %s", Arrays.toString( array );  // FAIL
     * 
     *      format("%s",          arrayToStrings( array );   // FAIL
     *      format("%s %s %s %s", arrayToStrings( array );   // GOOD
     */
    public static Object[] arrayToStrings(Object[] args) {
        return Arrays.stream(args).map(Object::toString).toArray();
    }
    
    /**
     * Consolidates debugging flags/variables and logging.
     * Call buglog() to write statements that go unprinted unless DEBUG
     *  is set or explicit class/hierarchy/method sources are requested.
     */
    public static class DebugHandler {
        private DebugHandler() {}   // utility class
        public static boolean DEBUG() { return DEBUG; }
        public static boolean VERBOSE() { return VERBOSE; }
        public static String debugFormat() { return debugFormat; }
        public static OptionalInt debugIntI() { return debugIntI; }
        public static OptionalInt debugIntK() { return debugIntK; }

        private static final boolean DEBUG = prop("esv-debug", false);
        private static final boolean VERBOSE = prop("esv-verbose", false);
        private static final String debugFormat = prop("esv-debugFormat", "%.4g");
        private static final OptionalInt debugIntI = prop("esv-debugIntI", OptionalInt.empty());
        private static final OptionalInt debugIntK = prop("esv-debugIntK", OptionalInt.empty());
    }
    
    /**
     * Use this to iterate easily over a collection of lists.
     * Returns the inferred parent list type. Duplicate elements allowed.
     * "View" intentional to avoid vain attempts at modifying the component lists.
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
    public static class UnderConstructionException extends RuntimeException {}

}
