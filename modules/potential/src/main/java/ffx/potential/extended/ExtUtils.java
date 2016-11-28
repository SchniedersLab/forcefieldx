package ffx.potential.extended;

import java.util.Enumeration;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;

import javax.swing.tree.TreeNode;

import org.apache.commons.lang.StringUtils;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.MSNode;
import ffx.potential.utils.DebugException;

import static ffx.potential.extended.ExtConstants.kB;
import static ffx.potential.extended.ExtConstants.roomTemperature;

/**
 * 
 * @author slucore
 */
public final class ExtUtils {

    // Private constructor implies static class.
    private ExtUtils() {}
    
    private static final Logger logger = Logger.getLogger(ExtUtils.class.getName());    
    private static final Random rng = ThreadLocalRandom.current();
    
    public static void THROW_DEBUG(String message) {
        throw new DebugException(message);
    }
    
    public static boolean prop(String key) {
        return (System.getProperty(key) != null && !System.getProperty(key).equalsIgnoreCase("false"));
    }
    public static <T> T prop(String key, T defaultVal) {
        try {
            if (defaultVal instanceof Integer) {
                return (System.getProperty(key) != null) 
                        ? (T) Integer.valueOf(System.getProperty(key)) : defaultVal;
            } else if (defaultVal instanceof Double) {
                return (System.getProperty(key) != null) 
                        ? (T) Double.valueOf(System.getProperty(key)) : defaultVal;
            } else if (defaultVal instanceof Boolean) {
                Boolean def = (Boolean) defaultVal;
                if (System.getProperty(key) != null) {
                    if (System.getProperty(key).equals("")) {
                        System.setProperty(key,"true");
                    }
                    return (T) Boolean.valueOf(System.getProperty(key));
                } else {
                    return defaultVal;
                }
            } else {    // else assume OptionalDouble
                return (System.getProperty(key) != null) 
                        ? (T) OptionalDouble.of(Double.valueOf(System.getProperty(key))) : (T) OptionalDouble.empty();
            }
        } catch (Exception ex) {
            String value = (System.getProperty(key) != null) ? System.getProperty(key).toString() : "null";
            logger.severe(String.format("Error parsing property %s with value %s; the default is %s.", 
                    key, defaultVal.toString(), value));
            throw ex;
        }
    }
    public static <T extends Enum<T>> T prop(Class<T> type, String key, T def) {
        return (System.getProperty(key) != null) ? T.valueOf(type, System.getProperty(key)) : def;
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
                sb.append(format(" >%s> %s\n", StringUtils.repeat(">>", node.getLevel()), node.getName()));
                init = true;
                continue;
            }
            sb.append(format(" └%s> %s\n", StringUtils.repeat("--", node.getLevel()), node.getName()));
        }
    }
    
}
