/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.potential.parsers;

import java.io.File;
import java.util.ArrayList;
import java.util.Vector;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.parameters.ForceField;
import java.util.List;

/**
 * The SystemFilter class is the base class for most Force Field X file parsers.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public abstract class SystemFilter {

    /**
     * This follows the TINKER file versioning scheme.
     *
     * @param file
     *            File to find a version for.
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
     * MolecularAssembly should be defined for PDB files with
     * alternate locations.
     */
    protected MolecularAssembly activeMolecularAssembly = null;
    /**
     * All MolecularAssembly instances defined. More than one
     * MolecularAssembly should be defined for PDB entries with
     * alternate locations.
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
     */
    public boolean fileRead() {
        return fileRead;
    }

    public int getAtomCount() {
        if (atomList == null) {
            return 0;
        }
        return atomList.size();
    }

    public ArrayList<Atom> getAtomList() {
        return atomList;
    }

    public int getBondCount() {
        if (bondList == null) {
            return 0;
        }
        return bondList.size();
    }

    /**
     * Return the MolecularSystem that has been read in
     */
    public MolecularAssembly getActiveMolecularSystem() {
        return activeMolecularAssembly;
    }

    public MolecularAssembly[] getMolecularAssemblys() {
        if (systems.size() > 0) {
            MolecularAssembly assemblies[] = new MolecularAssembly[systems.size()];
            return systems.toArray(assemblies);
        } else {
            MolecularAssembly assemblies[] = {activeMolecularAssembly};
            return assemblies;
        }
    }

    public FileType getType() {
        return fileType;
    }

    public File getFile() {
        return currentFile;
    }

    public List<File> getFiles() {
        return files;
    }

    /**
     * This method is different for each subclass and must be overidden
     */
    public abstract boolean readFile();

    public void setFileRead(boolean fileRead) {
        this.fileRead = fileRead;
    }

    public void setForceField(ForceField forceField) {
        this.forceField = forceField;
    }

    public void setProperties(CompositeConfiguration properties) {
        this.properties = properties;
    }

    public void setMolecularSystem(MolecularAssembly molecularAssembly) {
        activeMolecularAssembly = molecularAssembly;
    }

    public void setType(FileType fileType) {
        this.fileType = fileType;
    }

    public void setFile(File file) {
        this.currentFile = file;
        files = new ArrayList<File>();
        if (file != null) {
            files.add(file);
        }
    }

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
     * If the append flag is true, "saveFile" will be appended to. Otherwise
     * the default versioning scheme will be applied.
     */
    public abstract boolean writeFile(File saveFile, boolean append);
}
