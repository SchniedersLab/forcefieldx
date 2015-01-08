/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential.parsers;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;

/**
 * The SystemFilter class is the base class for most Force Field X file parsers.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public abstract class SystemFilter {

    /**
     * This follows the TINKER file versioning scheme.
     *
     * @param file File to find a version for.
     * @return File Versioned File.
     */
    public static File version(File file) {
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
            newFile = oldFile;
            int thousand = i / 1000;
            int hundred = (i - 1000 * thousand) / 100;
            int tens = (i - 1000 * thousand - 100 * hundred) / 10;
            int ones = i - 1000 * thousand - 100 * hundred - 10 * tens;
            StringBuilder newFileString = new StringBuilder(oldFile.getAbsolutePath());
            if (thousand != 0) {
                newFileString.append('_').append(thousand).append(hundred).append(tens).append(ones);
            } else if (hundred != 0) {
                newFileString.append('_').append(hundred).append(tens).append(
                        ones);
            } else if (tens != 0) {
                newFileString.append('_').append(tens).append(ones);
            } else {
                newFileString.append('_').append(ones);
            }
            newFile = new File(newFileString.toString());
        }
        return newFile;
    }

    /**
     * <p>
     * previousVersion</p>
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
                newFileString.append('_').append(hundred).append(tens).append(
                        ones);
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
     * The atomList is filled by filters that extend SystemFilter.
     */
    protected ArrayList<Atom> atomList = null;
    /**
     * The bondList may be filled by the filters that extend SystemFilter.
     */
    protected ArrayList<Bond> bondList = null;
    /**
     * The current MolecularAssembly being populated. Note that more than one
     * MolecularAssembly should be defined for PDB files with alternate
     * locations.
     */
    protected MolecularAssembly activeMolecularAssembly = null;
    /**
     * All MolecularAssembly instances defined. More than one MolecularAssembly
     * should be defined for PDB entries with alternate locations.
     */
    protected Vector<MolecularAssembly> systems = new Vector<MolecularAssembly>();
    /**
     * File currently being read.
     */
    protected File currentFile = null;
    /**
     * Append multiple files into one MolecularAssembly.
     */
    protected List<File> files = null;
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
    protected ForceField forceField = null;
    /**
     * True after the file has been read successfully.
     */
    protected boolean fileRead = false;

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param files a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public SystemFilter(List<File> files, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        this.files = files;
        if (files != null) {
            this.currentFile = files.get(0);
        }
        this.activeMolecularAssembly = molecularAssembly;
        this.forceField = forceField;
        this.properties = properties;
    }

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public SystemFilter(File file, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        files = new ArrayList<File>();
        if (file != null) {
            files.add(file);
        }
        this.currentFile = file;
        this.activeMolecularAssembly = molecularAssembly;
        this.forceField = forceField;
        this.properties = properties;
    }

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public SystemFilter(File file, List<MolecularAssembly> molecularAssemblies,
            ForceField forceField, CompositeConfiguration properties) {
        files = new ArrayList<File>();
        if (file != null) {
            files.add(file);
        }
        this.currentFile = file;
        this.systems = new Vector(molecularAssemblies);
        this.activeMolecularAssembly = systems.firstElement();
        this.forceField = forceField;
        this.properties = properties;
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
     * <p>
     * getAtomCount</p>
     *
     * @return a int.
     */
    public int getAtomCount() {
        if (atomList == null) {
            return 0;
        }
        return atomList.size();
    }

    /**
     * <p>
     * Getter for the field <code>atomList</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Atom> getAtomList() {
        return atomList;
    }

    /**
     * <p>
     * getBondCount</p>
     *
     * @return a int.
     */
    public int getBondCount() {
        if (bondList == null) {
            return 0;
        }
        return bondList.size();
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
     * <p>
     * getMolecularAssemblys</p>
     *
     * @return an array of {@link ffx.potential.MolecularAssembly} objects.
     */
    public MolecularAssembly[] getMolecularAssemblys() {
        if (systems.size() > 0) {
            MolecularAssembly assemblies[] = new MolecularAssembly[systems.size()];
            return systems.toArray(assemblies);
        } else {
            MolecularAssembly assemblies[] = {activeMolecularAssembly};
            return assemblies;
        }
    }

    /**
     * <p>
     * getType</p>
     *
     * @return a {@link ffx.potential.Utilities.FileType} object.
     */
    public FileType getType() {
        return fileType;
    }

    /**
     * <p>
     * getFile</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getFile() {
        return currentFile;
    }

    /**
     * <p>
     * Getter for the field <code>files</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<File> getFiles() {
        return files;
    }

    /**
     * This method is different for each subclass and must be overidden
     *
     * @return a boolean.
     */
    public abstract boolean readFile();

    /**
     * <p>
     * Setter for the field <code>fileRead</code>.</p>
     *
     * @param fileRead a boolean.
     */
    public void setFileRead(boolean fileRead) {
        this.fileRead = fileRead;
    }

    /**
     * <p>
     * Setter for the field <code>forceField</code>.</p>
     *
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     */
    public void setForceField(ForceField forceField) {
        this.forceField = forceField;
    }

    /**
     * <p>
     * Setter for the field <code>properties</code>.</p>
     *
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public void setProperties(CompositeConfiguration properties) {
        this.properties = properties;
    }

    /**
     * <p>
     * setMolecularSystem</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     */
    public void setMolecularSystem(MolecularAssembly molecularAssembly) {
        activeMolecularAssembly = molecularAssembly;
    }

    /**
     * <p>
     * setType</p>
     *
     * @param fileType a {@link ffx.potential.Utilities.FileType} object.
     */
    public void setType(FileType fileType) {
        this.fileType = fileType;
    }

    /**
     * <p>
     * setFile</p>
     *
     * @param file a {@link java.io.File} object.
     */
    public void setFile(File file) {
        this.currentFile = file;
        files = new ArrayList<File>();
        if (file != null) {
            files.add(file);
        }
    }

    /**
     * <p>
     * Setter for the field <code>files</code>.</p>
     *
     * @param files a {@link java.util.List} object.
     */
    public void setFiles(List<File> files) {
        this.files = files;
        if (files != null) {
            this.currentFile = files.get(0);
        } else {
            this.currentFile = null;
        }
    }

    /**
     * This method is different for each subclass and must be overidden.
     *
     * If the append flag is true, "saveFile" will be appended to. Otherwise the
     * default versioning scheme will be applied.
     *
     * @param saveFile a {@link java.io.File} object.
     * @param append a boolean.
     * @return a boolean.
     */
    public abstract boolean writeFile(File saveFile, boolean append);
}
