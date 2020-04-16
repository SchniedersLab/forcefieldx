// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
// ******************************************************************************
package ffx.ui;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Keyword;
import java.io.File;
import java.util.Hashtable;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * The FFXSystem class contains extensions to the generic MolecularAssembly class.
 *
 * @author Michael J. Schnieders
 */
public class FFXSystem extends MolecularAssembly {

  /** Constant <code>MultiScaleLevel=4</code> */
  public static final int MultiScaleLevel = 4;
  // Log file being used for modeling commands
  private File logFile;
  // Key file for this system
  private File keyFile;
  private Hashtable<String, Keyword> keywords = new Hashtable<String, Keyword>();
  private CompositeConfiguration properties;
  private String commandDescription;
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

  /**
   * addKeyword
   *
   * @param k a {@link ffx.utilities.Keyword} object.
   */
  public void addKeyword(Keyword k) {
    if (keywords.containsKey(k.getKeyword())) {
      return;
    }
    keywords.put(k.getKeyword(), k);
  }

  /** {@inheritDoc} */
  @Override
  public boolean destroy() {
    setClosing(true);
    return super.destroy();
  }

  /**
   * Getter for the field <code>keyFile</code>.
   *
   * @return a {@link java.io.File} object.
   */
  public File getKeyFile() {
    return keyFile;
  }

  /**
   * Setter for the field <code>keyFile</code>.
   *
   * @param f a {@link java.io.File} object.
   */
  public void setKeyFile(File f) {
    keyFile = f;
  }

  /**
   * getKeyword
   *
   * @param k a {@link java.lang.String} object.
   * @return a {@link ffx.utilities.Keyword} object.
   */
  public Keyword getKeyword(String k) {
    return keywords.get(k);
  }

  /**
   * Getter for the field <code>keywords</code>.
   *
   * @return a {@link java.util.Hashtable} object.
   */
  public Hashtable<String, Keyword> getKeywords() {
    return keywords;
  }

  /**
   * Setter for the field <code>keywords</code>.
   *
   * @param k a {@link java.util.Hashtable} object.
   */
  public void setKeywords(Hashtable<String, Keyword> k) {
    keywords = k;
  }

  /**
   * Getter for the field <code>logFile</code>.
   *
   * @return a {@link java.io.File} object.
   */
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

  /**
   * Setter for the field <code>logFile</code>.
   *
   * @param f a {@link java.io.File} object.
   */
  public void setLogFile(File f) {
    logFile = f;
  }

  /**
   * Getter for the field <code>properties</code>.
   *
   * @return a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  @Override
  public CompositeConfiguration getProperties() {
    return properties;
  }

  /**
   * Getter for the field <code>trajectory</code>.
   *
   * @return a {@link ffx.ui.Trajectory} object.
   */
  public Trajectory getTrajectory() {
    return trajectory;
  }

  /**
   * Setter for the field <code>trajectory</code>.
   *
   * @param t a {@link ffx.ui.Trajectory} object.
   */
  public void setTrajectory(Trajectory t) {
    trajectory = t;
  }

  /**
   * isClosing
   *
   * @return a boolean.
   */
  public boolean isClosing() {
    return closing;
  }

  /**
   * Setter for the field <code>closing</code>.
   *
   * @param b a boolean.
   */
  public void setClosing(boolean b) {
    closing = b;
  }

  /**
   * isStale
   *
   * @return a boolean.
   */
  public boolean isStale() {
    for (Atom a : getAtomList()) {
      if (a.isStale()) {
        return true;
      }
    }
    return false;
  }

  /**
   * removeKeyword
   *
   * @param kd a {@link ffx.utilities.Keyword} object.
   */
  public void removeKeyword(Keyword kd) {
    if (keywords.containsKey(kd.getKeyword())) {
      keywords.remove(kd.getKeyword());
    }
  }

  /**
   * Setter for the field <code>commandDescription</code>.
   *
   * @param command a {@link java.lang.String} object.
   */
  public void setCommandDescription(String command) {
    commandDescription = command;
  }

  /**
   * toFFString
   *
   * @return a {@link java.lang.String} object.
   */
  public String toFFString() {
    StringBuilder sb = new StringBuilder(toString());
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

  /**
   * toFileString
   *
   * @return a {@link java.lang.String} object.
   */
  public String toFileString() {
    if (getFile() == null) {
      return toFFString();
    }
    StringBuilder sb = new StringBuilder(getFile().getAbsolutePath());
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

  /** {@inheritDoc} */
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
