// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.parsers;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import static ffx.utilities.StringUtils.parseAtomRange;
import static java.lang.String.format;

/**
 * The SystemFilter class is the base class for most Force Field X file parsers.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class SystemFilter {

  protected static final Pattern lambdaPattern = Pattern.compile("Lambda: +([01]\\.\\d+)");
  private static final Logger logger = Logger.getLogger(SystemFilter.class.getName());
  private static Versioning vers = Versioning.TINKER;
  private static int absoluteCounter = 0;
  /**
   * Constant <code>dieOnMissingAtom=</code>
   */
  protected final boolean dieOnMissingAtom;
  /**
   * Standardize atom names to PDB standard by default.
   */
  protected boolean standardizeAtomNames;
  /**
   * True if atoms are to be printed to their van der Waals centers instead of nuclear centers
   * (applies primarily to hydrogen).
   */
  protected final boolean vdwH;
  /**
   * The atomList is filled by filters that extend SystemFilter.
   */
  protected List<Atom> atomList = null;
  /**
   * The bondList may be filled by the filters that extend SystemFilter.
   */
  protected List<Bond> bondList = null;
  /**
   * All MolecularAssembly instances defined. More than one MolecularAssembly should be defined for
   * PDB entries with alternate locations.
   */
  protected List<MolecularAssembly> systems = new Vector<>();
  /**
   * Append multiple files into one MolecularAssembly.
   */
  protected List<File> files;
  /**
   * The file format being handled.
   */
  protected FileType fileType = FileType.UNK;
  /**
   * Properties associated with this file.
   */
  protected CompositeConfiguration properties;
  /**
   * The molecular mechanics force field being used.
   */
  protected ForceField forceField;
  /**
   * True after the file has been read successfully.
   */
  protected boolean fileRead = false;
  /**
   * The current MolecularAssembly being populated. Note that more than one MolecularAssembly should
   * be defined for PDB files with alternate locations.
   */
  MolecularAssembly activeMolecularAssembly;
  /**
   * File currently being read.
   */
  File currentFile = null;

  /**
   * Initializations common to the all the constructors.
   *
   * @param forceField The force field.
   * @param properties The CompositeConfiguration properties.
   */
  private SystemFilter(ForceField forceField, CompositeConfiguration properties) {
    this.forceField = forceField;
    this.properties = properties;
    if (properties == null && forceField != null) {
      this.properties = forceField.getProperties();
    }
    if (this.properties != null) {
      vdwH = this.properties.getBoolean("vdwHydrogens", false);
      dieOnMissingAtom = this.properties.getBoolean("trajectory-dieOnMissing", false);
      standardizeAtomNames = this.properties.getBoolean("standardizeAtomNames", true);
    } else {
      logger.info(" SystemFilter: Using default values due to no properties.");
      vdwH = false;
      dieOnMissingAtom = false;
      standardizeAtomNames = true;
    }
  }

  /**
   * Constructor for SystemFilter.
   *
   * @param files             a {@link java.util.List} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
   * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *                          object.
   */
  public SystemFilter(List<File> files, MolecularAssembly molecularAssembly, ForceField forceField,
                      CompositeConfiguration properties) {
    this(forceField, properties);
    this.files = files;
    if (files != null) {
      this.currentFile = files.get(0);
    }
    this.activeMolecularAssembly = molecularAssembly;
  }

  /**
   * Constructor for SystemFilter.
   *
   * @param file              a {@link java.io.File} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
   * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *                          object.
   */
  public SystemFilter(File file, MolecularAssembly molecularAssembly, ForceField forceField,
                      CompositeConfiguration properties) {
    this(forceField, properties);
    files = new ArrayList<>();
    if (file != null) {
      files.add(file);
    }
    this.currentFile = file;
    this.activeMolecularAssembly = molecularAssembly;
  }

  /**
   * Constructor for SystemFilter.
   *
   * @param file                a {@link java.io.File} object.
   * @param molecularAssemblies a {@link java.util.List} object.
   * @param forceField          a {@link ffx.potential.parameters.ForceField} object.
   * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *                            object.
   */
  public SystemFilter(File file, List<MolecularAssembly> molecularAssemblies, ForceField forceField,
                      CompositeConfiguration properties) {
    this(forceField, properties);
    files = new ArrayList<>();
    if (file != null) {
      files.add(file);
    }
    this.currentFile = file;
    this.systems = new ArrayList<>(molecularAssemblies);
    this.activeMolecularAssembly = systems.get(0);
  }

  /**
   * previousVersion
   *
   * @param file a {@link java.io.File} object.
   * @return a {@link java.io.File} object.
   */
  public static File previousVersion(File file) {
    if (file == null) {
      return null;
    }
    String fileName = file.getAbsolutePath();
    int dot = file.getAbsolutePath().lastIndexOf(".");
    int under = file.getAbsolutePath().lastIndexOf("_");
    File newFile = file;
    if (under > dot) {
      String name = fileName.substring(0, under);
      newFile = new File(name);
    }
    File baseFile = newFile;
    File previousFile = null;
    int i = 1;
    while (newFile.exists()) {
      i = i + 1;
      previousFile = newFile;
      newFile = baseFile;
      int thousand = i / 1000;
      int hundred = (i - 1000 * thousand) / 100;
      int tens = (i - 1000 * thousand - 100 * hundred) / 10;
      int ones = i - 1000 * thousand - 100 * hundred - 10 * tens;
      StringBuilder newFileString = new StringBuilder(baseFile.getAbsolutePath());
      if (thousand != 0) {
        newFileString.append('_').append(thousand).append(hundred).append(tens).append(ones);
      } else if (hundred != 0) {
        newFileString.append('_').append(hundred).append(tens).append(ones);
      } else if (tens != 0) {
        newFileString.append('_').append(tens).append(ones);
      } else {
        newFileString.append('_').append(ones);
      }
      newFile = new File(newFileString.toString());
    }
    return previousFile;
  }

  /**
   * Negative: prefix a version number; Positive: postfix; Zero: TINKER-style.
   *
   * @param vers a {@link ffx.potential.parsers.SystemFilter.Versioning} object.
   */
  public static void setVersioning(Versioning vers) {
    SystemFilter.vers = vers;
  }

  /**
   * Use setVersioning() to choose between prefix/postfix.
   *
   * @param file a {@link java.io.File} object.
   * @return a {@link java.io.File} object.
   */
  public static File version(File file) {
    if (vers == Versioning.TINKER) {
      return versionTinker(file);
    } else {
      if (vers == Versioning.PREFIX_ABSOLUTE || vers == Versioning.POSTFIX_ABSOLUTE) {
        return versionAbsolute(file, (vers == Versioning.PREFIX_ABSOLUTE));
      } else {
        return version(file, (vers == Versioning.PREFIX));
      }
    }
  }

  private static File version(File file, boolean prefix) {
    if (file == null || !(file.exists())) {
      return file;
    }
    String fn = file.getAbsolutePath();
    int dot = fn.lastIndexOf(".");
    int under = fn.lastIndexOf("_");
    if (dot < 0) {
      fn += ".unk";
      dot = fn.lastIndexOf(".");
    }
    if (under < 0) {
      under = dot;
    }
    String name = (prefix) ? fn.substring(0, under) : fn.substring(0, dot);
    String extension = (prefix) ? fn.substring(dot + 1) : fn.substring(dot + 1, under);
    int number = 0;
    String newFn = (prefix) ? format("%s_%d.%s", name, number, extension)
        : format("%s.%s_%d", name, extension, number);
    if (prefix && under < dot) {
      try {
        number = Integer.parseInt(fn.substring(under + 1, dot));
      } catch (NumberFormatException ex) {
        // Then we have something like "AKA_dyn.pdb"
        name = fn.substring(0, dot);
        number++;
        newFn = format("%s_%d.%s", name, number, extension);
      }
    } else if (!prefix && under > dot) {
      try {
        number = Integer.parseInt(fn.substring(under + 1));
        number++;
      } catch (NumberFormatException ex) {
        //
      }
    }
    File newFile = new File(newFn);
    while (newFile.exists()) {
      ++number;
      newFn = (prefix) ? format("%s_%d.%s", name, number, extension)
          : format("%s.%s_%d", name, extension, number);
      newFile = new File(newFn);
    }
    return newFile;
  }

  private static synchronized File versionAbsolute(File file, boolean prefix) {
    if (file == null || !(file.exists())) {
      return file;
    }
    String fn = file.getAbsolutePath();
    int dot = fn.lastIndexOf(".");
    int under = fn.lastIndexOf("_");
    if (dot < 0) {
      fn += ".unk";
      dot = fn.lastIndexOf(".");
    }
    if (under < 0) {
      under = dot;
    }
    String name = (prefix) ? fn.substring(0, under) : fn.substring(0, dot);
    String extension = (prefix) ? fn.substring(dot + 1) : fn.substring(dot + 1, under);
    String newFn = (prefix) ? format("%s_%d.%s", name, absoluteCounter, extension)
        : format("%s.%s_%d", name, extension, absoluteCounter);
    File newFile = new File(newFn);
    while (newFile.exists()) {
      absoluteCounter++;
      newFn = (prefix) ? format("%s_%d.%s", name, absoluteCounter, extension)
          : format("%s.%s_%d", name, extension, absoluteCounter);
      newFile = new File(newFn);
    }
    return newFile;
  }

  /**
   * This follows the TINKER file versioning scheme.
   *
   * @param file File to find a version for.
   * @return File Versioned File.
   */
  private static File versionTinker(File file) {
    if (file == null) {
      return null;
    }
    if (!file.exists()) {
      return file;
    }
    String fileName = file.getAbsolutePath();
    int dot = file.getAbsolutePath().lastIndexOf(".");
    int under = file.getAbsolutePath().lastIndexOf("_");
    File newFile = file;
    if (under > dot) {
      String name = fileName.substring(0, under);
      newFile = new File(name);
    }
    File oldFile = newFile;
    int i = 1;
    while (newFile.exists()) {
      i = i + 1;
      String newFileString = format("%s_%d", oldFile.getAbsolutePath(), i);
      newFile = new File(newFileString);
    }
    return newFile;
  }

  /**
   * Converts a list of atom indices to an array of atoms.
   *
   * @param atomList List of atom indices.
   * @param atoms    Array of atoms.
   * @return Array of atoms.
   */
  public static Set<Atom> atomListToSet(List<Integer> atomList, Atom[] atoms) {
    Set<Atom> atomSet = new HashSet<>();
    for (int i = 0; i < atomList.size(); i++) {
      atomSet.add(atoms[atomList.get(i)]);
    }
    return atomSet;
  }

  /**
   * Automatically sets atom-specific flags, particularly nouse and inactive, and apply harmonic
   * restraints. Intended to be called at the end of readFile() implementations.
   *
   * <p>Supported syntax: "(\\d+)-(\\d+)"
   */
  public void applyAtomProperties() {
    /*
      What may be a more elegant implementation is to make readFile() a
      public concrete, but skeletal method, and then have readFile() call a
      protected abstract readFile method for each implementation.
    */
    Atom[] atomArray = activeMolecularAssembly.getAtomArray();
    int nAtoms = atomArray.length;
    String[] nouseKeys = properties.getStringArray("nouse");
    for (String nouseKey : nouseKeys) {
      String[] toks = nouseKey.split("\\s+");
      for (String tok : toks) {
        try {
          List<Integer> nouseRange = parseAtomRange("nouse", tok, nAtoms);
          for (int j : nouseRange) {
            atomArray[j].setUse(false);
          }
        } catch (IllegalArgumentException ex) {
          boolean atomFound = false;
          for (Atom atom : atomArray) {
            if (atom.getName().equalsIgnoreCase(tok)) {
              atomFound = true;
              atom.setUse(false);
            }
          }
          if (atomFound) {
            logger.info(format(" Setting atoms with name %s to not be used", tok));
          } else {
            logger.log(Level.INFO, ex.getLocalizedMessage());
          }
        }
      }
    }

    if (properties.containsKey("active")) {
      for (Atom atom : atomArray) {
        atom.setActive(false);
      }
      String[] activeKeys = properties.getStringArray("active");
      for (String inactiveKey : activeKeys) {
        try {
          List<Integer> inactiveRange = parseAtomRange("inactive", inactiveKey, nAtoms);
          for (int i : inactiveRange) {
            atomArray[i].setActive(false);
          }
        } catch (IllegalArgumentException ex) {
          logger.log(Level.INFO, ex.getLocalizedMessage());
        }
      }
    } else if (properties.containsKey("inactive")) {
      for (Atom atom : atomArray) {
        atom.setActive(true);
      }
      String[] inactiveKeys = properties.getStringArray("inactive");
      for (String inactiveKey : inactiveKeys) {
        try {
          List<Integer> inactiveRange = parseAtomRange("inactive", inactiveKey, nAtoms);
          for (int i : inactiveRange) {
            atomArray[i].setActive(false);
          }
        } catch (IllegalArgumentException ex) {
          logger.log(Level.INFO, ex.getLocalizedMessage());
        }
      }
    }

    String[] noElStrings = properties.getStringArray("noElectro");
    for (String noE : noElStrings) {
      String[] toks = noE.split("\\s+");
      for (String tok : toks) {
        try {
          List<Integer> noERange = parseAtomRange("noElectro", tok, nAtoms);
          for (int i : noERange) {
            atomArray[i].setElectrostatics(false);
          }
        } catch (IllegalArgumentException ex) {
          boolean atomFound = false;
          for (Atom atom : atomArray) {
            if (atom.getName().equalsIgnoreCase(tok)) {
              atomFound = true;
              atom.setElectrostatics(false);
            }
          }
          if (atomFound) {
            logger.info(format(" Disabled electrostatics for atoms with name %s", tok));
          } else {
            logger.log(Level.INFO, format(" No electrostatics input %s could not be parsed as a numerical "
                + "range or atom type present in assembly", tok));
          }
        }
      }
    }
  }

  /**
   * Attempts to close any open resources associated with the underlying file; primarily to be used
   * when finished reading a trajectory.
   */
  public abstract void closeReader();

  public int countNumModels() {
    return -1;
  }

  /**
   * Returns true if the read was successful
   *
   * @return a boolean.
   */
  public boolean fileRead() {
    return fileRead;
  }

  /**
   * Return the MolecularSystem that has been read in
   *
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  public MolecularAssembly getActiveMolecularSystem() {
    return activeMolecularAssembly;
  }

  /**
   * Getter for the field <code>atomList</code>.
   *
   * @return a {@link java.util.List} object.
   */
  public List<Atom> getAtomList() {
    return atomList;
  }

  /**
   * getFile
   *
   * @return a {@link java.io.File} object.
   */
  public File getFile() {
    return currentFile;
  }

  /**
   * setFile
   *
   * @param file a {@link java.io.File} object.
   */
  public void setFile(File file) {
    this.currentFile = file;
    files = new ArrayList<>();
    if (file != null) {
      files.add(file);
    }
  }

  /**
   * Getter for the field <code>files</code>.
   *
   * @return a {@link java.util.List} object.
   */
  public List<File> getFiles() {
    return files;
  }

  /**
   * Setter for the field <code>files</code>.
   *
   * @param files a {@link java.util.List} object.
   */
  public void setFiles(List<File> files) {
    this.files = files;
    if (files != null && !files.isEmpty()) {
      this.currentFile = files.get(0);
    } else {
      this.currentFile = null;
    }
  }

  /**
   * Gets the last read lambda value read by the filter, if any.
   *
   * @return Last lambda value read by this filter.
   */
  public OptionalDouble getLastReadLambda() {
    return OptionalDouble.empty();
  }

  /**
   * Get the MolecularAssembly array.
   *
   * @return an array of {@link ffx.potential.MolecularAssembly} objects.
   */
  public MolecularAssembly[] getMolecularAssemblyArray() {
    if (!systems.isEmpty()) {
      return systems.toArray(new MolecularAssembly[0]);
    } else {
      return new MolecularAssembly[]{activeMolecularAssembly};
    }
  }

  /**
   * Gets all remark lines read by the last readFile or readNext call.
   *
   * @return Array of Strings representing remark lines, if any.
   */
  public String[] getRemarkLines() {
    return new String[0];
  }

  /**
   * Return snapshot number.
   *
   * @return The snapshot number.
   */
  public int getSnapshot() {
    return -1;
  }

  /**
   * getType
   *
   * @return a {@link ffx.potential.Utilities.FileType} object.
   */
  public FileType getType() {
    return fileType;
  }

  /**
   * setType
   *
   * @param fileType a {@link ffx.potential.Utilities.FileType} object.
   */
  public void setType(FileType fileType) {
    this.fileType = fileType;
  }

  /**
   * This method is different for each subclass and must be overridden.
   *
   * @return a boolean.
   */
  public abstract boolean readFile();

  /**
   * Reads the next model if applicable (currently, ARC and PDB files only).
   *
   * @return If next model read.
   */
  public abstract boolean readNext();

  /**
   * Reads the next model if applicable (currently, ARC files only).
   *
   * @param resetPosition Resets to first frame.
   * @return If next model read.
   */
  public abstract boolean readNext(boolean resetPosition);

  /**
   * Reads the next model if applicable (currently, ARC files only).
   *
   * @param resetPosition Resets to first frame.
   * @param print         Flag to print.
   * @return If next model read.
   */
  public abstract boolean readNext(boolean resetPosition, boolean print);

  /**
   * Reads the next model if applicable (currently, ARC files only).
   *
   * @param resetPosition Resets to first frame.
   * @param print         Flag to print.
   * @param parse         Parse data in file. May want to skip structures for parallel jobs.
   * @return If next model read.
   */
  public abstract boolean readNext(boolean resetPosition, boolean print, boolean parse);

  /**
   * Setter for the field <code>forceField</code>.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public void setForceField(ForceField forceField) {
    this.forceField = forceField;
  }

  /**
   * Setter for the field <code>properties</code>.
   *
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *                   object.
   */
  public void setProperties(CompositeConfiguration properties) {
    this.properties = properties;
  }

  /**
   * This method is different for each subclass and must be overridden.
   *
   * <p>If the append flag is true, "saveFile" will be appended to. Otherwise, the default
   * versioning
   * scheme will be applied.
   *
   * @param saveFile a {@link java.io.File} object.
   * @param append   a boolean.
   * @return a boolean.
   */
  public boolean writeFile(File saveFile, boolean append) {
    return writeFile(saveFile, append, null);
  }

  /**
   * This method is different for each subclass and must be overridden.
   *
   * <p>If the append flag is true, "saveFile" will be appended to. Otherwise, the default
   * versioning
   * scheme will be applied.
   *
   * @param saveFile   a {@link java.io.File} object.
   * @param append     a boolean.
   * @param extraLines Additional lines to append to a comments section, or null.
   * @return a boolean.
   */
  public abstract boolean writeFile(File saveFile, boolean append, String[] extraLines);

  /**
   * Setter for the field <code>fileRead</code>.
   *
   * @param fileRead a boolean.
   */
  protected void setFileRead(boolean fileRead) {
    this.fileRead = fileRead;
  }

  /**
   * setMolecularSystem
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  void setMolecularSystem(MolecularAssembly molecularAssembly) {
    activeMolecularAssembly = molecularAssembly;
  }

  public enum Versioning {
    TINKER, PREFIX, POSTFIX, PREFIX_ABSOLUTE, POSTFIX_ABSOLUTE, NONE
  }
}
