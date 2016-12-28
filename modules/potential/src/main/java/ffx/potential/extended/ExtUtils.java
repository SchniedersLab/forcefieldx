package ffx.potential.extended;

import java.io.File;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
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
import ffx.potential.utils.DebugException;
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
    
    private static Logger logger = Logger.getLogger(ExtUtils.class.getName());    
    private static final Random rng = ThreadLocalRandom.current();
    
//    public void setR(double r[]) {
//        setR(r, 0.0);
//    }
//    
//    public void setR(double r[], double buffer) {
//        switch (coordinates) {
//            case QI:
//                setQIRotationMatrix(r[0], r[1], r[2]);
//                double zl = r[2] + buffer;
//                r2 = r[0] * r[0] + r[1] * r[1] + zl * zl;
//                x = 0.0;
//                y = 0.0;
//                z = sqrt(r2);
//                R = z;
//                break;
//            case GLOBAL:
//                r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//                x = r[0];
//                y = r[1];
//                z = r[2] + buffer;
//                R = sqrt(r2);
//                break;
//        }
//    }
    
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
//        molecularAssembly.finalize(true, forceField);
    }
    
    public static void setLogSource(Logger logger) {
        ExtUtils.logger = logger;
    }
    
    /**
     * ("log-format") Shorthand for the ubiquitous logger.info(String.format("why",42));
     */
    public static void logf(String msg, Object... args) {
        logger.info(format(msg, args));
    }
    
    /**
     * ("log-format") Shorthand for the ubiquitous logger.info(String.format("why",42));
     */
    public static void logfn(String msg, Object... args) {
        logger.info(format(msg + "\n", args));
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
    
    public static void THROW_DEBUG(String message) {
        throw new DebugException(message);
    }

    /**
     * Parse system properties of type String, Boolean, Integer, Double, OptionalInt, and OptionalDouble.
     * Default value determines the return type; null assigned to String.
     */
    public static <T> T prop(String key, T defaultVal)
            throws NoSuchElementException, NumberFormatException {
        if (defaultVal == null) {
            throw new IllegalArgumentException();
        }
        T parsed = defaultVal;
        try {
            if (System.getProperty(key) == null) {
                if (System.getProperty(key.toLowerCase()) != null) {
                    key = key.toLowerCase();
                } else if (System.getProperty(key.toLowerCase()) != null) {
                    key = key.toUpperCase();
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
                    if (System.getProperty(key).equals("")) {
                        System.setProperty(key,"true");
                    }
                    parsed = (T) Boolean.valueOf(System.getProperty(key));
                }
            } else {
                throw new IllegalArgumentException();
            }
        } catch (Exception ex) {
            String value = (System.getProperty(key) != null) ? System.getProperty(key) : "null";
            logger.warning(String.format("Error parsing property %s with value %s; the default is %s.", 
                    key, value, defaultVal.toString()));
            throw ex;
        }
        logfn(" ESV Properties Manager: %s = %s", key, parsed.toString());
        return parsed;
    }
    
    /**
     * Parse system properties into an Enum class.
     */
    public static <T extends Enum<T>> T prop(Class<T> type, String key, T def)
            throws IllegalArgumentException {
        T parsed = (System.getProperty(key) != null) ? T.valueOf(type, System.getProperty(key)) : def;
        logfn(" ESV Properties Manager: %s = %s", key, parsed.toString());
        return parsed;
    }
    
    /**
     * Parse multiple system properties into a single boolean; useful for handling synonyms.
     */
    public static boolean prop(String[] keys, boolean defaultValue) {
        if (defaultValue) { // When default is true, any explicitly false member kills it.
            for (String key : keys) {
                if (!prop(key, defaultValue)) {
                    return false;
                }
            }
        } else {            // When default is false, any explicitly true member activates it.
            for (String key : keys) {
                if (prop(key, defaultValue)) {
                    return true;
                }
            }
        }
        return defaultValue;
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
            logger.info(sb.toString());
            clear();
        }
        public static void printIfPresent(String header, Object... args) {
            if (queued > 0) {
                logger.info(format(header, args));
                print();
            }
        }
    }
    
    /**
     * Consolidates debugging flags/variables and logging.
     * Call buglog() to write statements that go unprinted unless DEBUG
     *  is set or explicit class/hierarchy/method sources are requested.
     */
    public static class DebugHandler {
        
        public static void buglog(String msg, Object... args) {
            buglog.info(format(msg, args));
        }
        public static void buglog(int depr, String msg, Object... args) {
            buglog(msg, args);
        }
        public static void buglog(int depr, StringBuilder sb) {
            buglog(sb.toString());
        }
        public static void buglog(StringBuilder sb) {
            buglog(sb.toString());
        }
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
        private static final Logger buglog = Logger.getLogger(DebugHandler.class.getName());
        private static final List<String> debugClasses = new ArrayList<>();
        private static final List<String> debugMethods = new ArrayList<>();
        
        static {
            buglog.setLevel(Level.ALL);
        }
    }
    
    public static class UnderConstructionException extends RuntimeException {}
    
}
