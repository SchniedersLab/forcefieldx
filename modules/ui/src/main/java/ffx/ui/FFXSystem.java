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

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Keyword;

/**
 * The FFXSystem class contains extensions to the generic
 * ffe.lang.MolecularAssembly class specific to Force Field X interacting
 * with TINKER.
 */
public class FFXSystem extends MolecularAssembly {

    private static final Logger logger = Logger.getLogger(FFXSystem.class.getName());
    private static final long serialVersionUID = 50L;
    public static final int MultiScaleLevel = 4;
    // Log file being used for modeling commands
    private File logFile;
    // Key file for this system
    private File keyFile;

    private Hashtable<String, Keyword> keywords = new Hashtable<String, Keyword>();
    private CompositeConfiguration properties = null;

    private String commandDescription = null;
    // Archive
    private Trajectory trajectory = null;
    // Flag to indicate this System is being closed
    private boolean closing = false;


    /**
     * Constructor.
     * 
     * @param file Coordinate file.
     * @param description Short description of the command that created this system.
     * @param properties Properties controlling operations on this system.
     */
    public FFXSystem(File file, String description, CompositeConfiguration properties) {
        super(FilenameUtils.getBaseName(file.getName()));
        setFile(file);
        commandDescription = description;
        this.properties = properties;
    }

    public CompositeConfiguration getProperties() {
        return properties;
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
