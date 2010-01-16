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
package ffx.ui;

import java.io.File;
import java.util.Hashtable;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.apache.commons.configuration.SystemConfiguration;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Keyword;
import java.util.Iterator;
import java.util.logging.Level;
import org.apache.commons.configuration.Configuration;

/**
 * The FFXSystem class contains extensions to the generic
 * ffe.lang.MolecularAssembly class specific to Force Field X interacting
 * with TINKER.
 */
public class FFXSystem extends MolecularAssembly {

    private static final Logger logger = Logger.getLogger(ffx.ui.FFXSystem.class.getName());
    private static final long serialVersionUID = 50L;
    public static final int MultiScaleLevel = 4;
    // Log file being used for modeling commands
    private File logFile;
    // Key file for this system
    private File keyFile;
    private Hashtable<String, Keyword> keywords = new Hashtable<String, Keyword>();
    // Command Description if this System is the result of a TINKER commad
    private String commandDescription = null;
    // Archive
    private Trajectory trajectory = null;
    // Simulation data
    //private double time, temperature, energy;
    private int step;
    // Flag to indicate this System is being closed
    private boolean closing = false;
    private CompositeConfiguration properties = null;

    /**
     * FFXSystem Constructor
     *
     * @param name
     *            String
     */
    public FFXSystem(String name, String description, File file) {
        super(name);
        setFile(file);
        commandDescription = description;
        loadProperties();
    }

    /**
     * This method sets up configuration properties in the following precedence
     * order:
     * 1.) Java system properties
     *     a.) -Dkey=value from the Java command line
     *     b.) System.setProperty("key","value") within Java code.
     *
     * 2.) Structure specific properties (for example pdbname.properties)
     *
     * 3.) User specific properties (~/.ffx/ffx.properties)
     *
     * 4.) System wide properties (file defined by environment variable FFX_PROPERTIES)
     */
    private void loadProperties() {
        /**
         * Command line options take precedences.
         */
        properties = new CompositeConfiguration();
        properties.addConfiguration(new SystemConfiguration());

        /**
         * Structure specific options are 2nd.
         */
        File file = getFile();
        String filename = file.getAbsolutePath();
        filename = org.apache.commons.io.FilenameUtils.removeExtension(filename);
        filename = filename + ".properties";
        File structurePropFile = new File(filename);
        if (structurePropFile.exists() && structurePropFile.canRead()) {
            try {
                properties.addConfiguration(new PropertiesConfiguration(structurePropFile));
            } catch (Exception e) {
                logger.info("Error loading " + filename + ".");
            }
        }

        /**
         * User specific options are 3rd.
         */
        filename = System.getProperty("user.home") + File.separator + ".ffx/ffx.properties";
        File userPropFile = new File(filename);
        if (userPropFile.exists() && userPropFile.canRead()) {
            try {
                properties.addConfiguration(new PropertiesConfiguration(userPropFile));
            } catch (Exception e) {
                logger.info("Error loading " + filename + ".");
            }

        }

        /**
         * System wide options are last.
         */
        filename = System.getenv("FFX_PROPERTIES");
        if (filename != null) {
            File systemPropFile = new File(filename);
            if (systemPropFile.exists() && systemPropFile.canRead()) {
                try {
                    properties.addConfiguration(new PropertiesConfiguration(systemPropFile));
                } catch (Exception e) {
                    logger.info("Error loading " + filename + ".");
                }
            }
        }

        /**
         * Echo the interpolated configuration.
         */
        if (logger.isLoggable(Level.FINE)) {
            Configuration config = properties.interpolatedConfiguration();
            Iterator<String> i = config.getKeys();
            while (i.hasNext()) {
                String s = i.next();
                logger.info("Key: " + s + ", Value: " + config.getString(s));
            }
        }
    }

    public void addKeyword(Keyword k) {
        if (keywords.containsKey(k.getKeyword())) {
            return;
        }
        keywords.put(k.getKeyword(), k);
    }

    @Override
    public boolean destroy() {
        setClosing(true);
        return super.destroy();
    }

    public File getKeyFile() {
        return keyFile;
    }

    public Keyword getKeyword(String k) {
        return keywords.get(k);
    }

    public Hashtable<String, Keyword> getKeywords() {
        return keywords;
    }

    public File getLogFile() {
        if (logFile == null) {
            if (getFile() == null) {
                return null;
            }
            String fileName = getFile().getName();
            int dot = fileName.lastIndexOf(".");
            fileName = fileName.subSequence(0, dot) + ".log";
            logFile = new File(fileName);
        }
        return logFile;
    }

    public String getStepString() {
        return String.format("Step: %12d", step);
    }

    public Trajectory getTrajectory() {
        return trajectory;
    }

    public boolean isClosing() {
        return closing;
    }

    public boolean isStale() {
        for (Atom a : getAtomList()) {
            if (a.isStale()) {
                return true;
            }
        }
        return false;
    }

    public void removeKeyword(Keyword kd) {
        if (keywords.containsKey(kd.getKeyword())) {
            keywords.remove(kd.getKeyword());
        }
    }

    public void setClosing(boolean b) {
        closing = b;
    }

    public void setCommandDescription(String command) {
        commandDescription = command;
    }

    public void setKeyFile(File f) {
        keyFile = f;
    }

    public void setKeywords(Hashtable<String, Keyword> k) {
        keywords = k;
    }

    public void setLogFile(File f) {
        logFile = f;
    }

    public void setStep(int s) {
        step = s;
    }

    public void setTrajectory(Trajectory t) {
        trajectory = t;
    }

    public String toFFString() {
        StringBuffer sb = new StringBuffer(toString());
        if (forceField != null) {
            String ff = forceField.toString("forcefield");
            if (ff != null) {
                ff = ff.substring(10).trim();
                sb.append(" (");
                sb.append(ff);
                sb.append(")");
            }
        }
        return sb.toString();
    }

    public String toFileString() {
        if (getFile() == null) {
            return toFFString();
        }
        StringBuffer sb = new StringBuffer(getFile().getAbsolutePath());
        if (forceField != null) {
            String ff = forceField.toString("forcefield");
            if (ff != null) {
                ff = ff.substring(10).trim();
                sb.append(" (");
                sb.append(ff);
                sb.append(")");
            }
        }
        return sb.toString();
    }

    @Override
    public String toString() {
        if (getFile() != null) {
            if (commandDescription != null) {
                return getFile().getName() + " (" + commandDescription + ")";
            }
            return getFile().getName();
        }
        if (getName() != null) {
            if (commandDescription != null) {
                return getName() + commandDescription;
            }
            return getName();
        }
        return "FFX System";
    }
}
