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

import java.io.File;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;

import javax.swing.tree.TreeNode;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.lang.StringUtils;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MSNode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

import static ffx.potential.extended.ExtConstants.kB;
import static ffx.potential.extended.ExtConstants.roomTemperature;
import static ffx.potential.extended.ExtUtils.DebugHandler.debugFormat;

/**
 * 
 * @author slucore
 */
public final class ExtUtils {

    // Private constructor implies static class.
    private ExtUtils() {}
    
    private static final Logger logger = Logger.getLogger(ExtUtils.class.getName());    
    private static final Random rng = ThreadLocalRandom.current();
    private static final CallerID cid = new CallerID();
    private static final CompositeConfiguration config = new CompositeConfiguration();
    
    public static void mutateToProtonatedForms(MolecularAssembly mola) {
        if (true) {
            throw new UnderConstructionException();
        }
        File structure = mola.getFile();
        MolecularAssembly molecularAssembly = new MolecularAssembly("AllProtonated");
        molecularAssembly.setFile(structure);
        
        CompositeConfiguration properties = Keyword.loadProperties(structure);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);
        
        PDBFilter pdbFilter = new PDBFilter(structure, molecularAssembly, forceField, properties);
        
        // TODO: 
//        pdbFilter.mutate(chain,resID,resName);
//        pdbFilter.readFile();
//        pdbFilter.applyAtomProperties();
    }
    
    /**
     * ("log-format") Shorthand for the ubiquitous logger.info(String.format("why",42));
     */
    public static void logf(String msg, Object... args) {
        cid.getCallingLogger().info(format(msg, args));
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
                return (System.getProperty(key) != null) 
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
                        System.setProperty(key,"true");
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
        config.addProperty(key, parsed);
        return parsed;
    }
    
    /**
     * Parse system properties into an Enum class.
     */
    public static final <T extends Enum<T>> T prop(Class<T> type, String key, T def)
            throws IllegalArgumentException {
        T parsed = (System.getProperty(key) != null) ? T.valueOf(type, System.getProperty(key)) : def;
        config.addProperty(key, parsed);
//        logf(" ESV Properties Manager: %s = %s", key, parsed.toString());
        return parsed;
    }
        
    public static final CompositeConfiguration getExtConfig() {
        return config;
    }
    
    public static void printExtConfig() {
        Iterator<String> it = config.getKeys();
        SB.logfn("\n Extended Configuration: ");
        while (it.hasNext()) {
            String key = it.next();
            SB.logfn("   %16s = %6s", it, config.getProperty(key));
        }
        SB.nl();
        SB.print();
    }
    
    /**
     * Return room-temperature velocities from a Maxwell-Boltzmann distribution of momenta.
     */
    public static double[] singleRoomtempMaxwell(double mass) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = rng.nextGaussian() * sqrt(kB * roomTemperature / mass);
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
            vv[i] = rng.nextGaussian() * sqrt(kB * temperature / mass);
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
        Enumeration fromRoot = target.pathFromAncestorEnumeration(root);
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
    
    public static class SB {
        private static StringBuilder sb = new StringBuilder();
        private static int queued = 0;
        public static void clear() {
            sb = new StringBuilder();
            queued = 0;
        }
        public static void logf(String msg, Object... args) {
            sb.append(format(msg, args));
            queued++;
        }
        public static void logfn(String msg, Object... args) {
            sb.append(format(msg, args)).append("\n");
            queued++;
        }
        public static void nl() {
            sb.append("\n");
        }
        public static void print() {
            cid.getCallingLogger().info(sb.toString());
            clear();
        }
        public static void printIf(boolean print) {
            if (print) {
                print();
            } else {
                clear();
            }
        }
        public static void printIfPresent(String header, Object... args) {
            if (queued > 0) {
                cid.getCallingLogger().info(format(header, args));
                print();
            } else {
                clear();
            }
        }
    }
    
    /**
     * Consolidates debugging flags/variables and logging.
     * Call buglog() to write statements that go unprinted unless DEBUG
     *  is set or explicit class/hierarchy/method sources are requested.
     */
    public static class DebugHandler {
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
    
    private static class CallerID extends SecurityManager {
        public Class<?> getCallingClass() {
            Class<?>[] callStack = getClassContext();
            for (Class<?> caller : callStack) {
                if (!caller.getName().contains("ExtUtils")) {
                    return caller;
                }
            }
            return ExtUtils.class;  // fallback
        }
        public Logger getCallingLogger() {
            return Logger.getLogger(getCallingClass().getName());
        }
    }
    
    public static class UnderConstructionException extends RuntimeException {}
    
}
