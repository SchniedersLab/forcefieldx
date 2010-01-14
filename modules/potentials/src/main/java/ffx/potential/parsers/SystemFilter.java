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
package ffx.potential.parsers;

import java.io.File;
import java.util.ArrayList;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Keyword;
import java.util.Hashtable;

/**
 * The SystemFilter class is the base class for most file parsers.
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
            StringBuffer newFileString = new StringBuffer(oldFile.getAbsolutePath());
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
            StringBuffer newFileString = new StringBuffer(baseFile.getAbsolutePath());
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
    protected MolecularAssembly molecularAssembly = null;
    protected FileType fileType = FileType.UNK;
    /**
     * The atomList and bondList are filled by the filters that extend this base
     * class.
     */
    protected ArrayList<Atom> atomList = null;
    protected ArrayList<Bond> bondList = null;
    protected ForceField forceField = null;
    protected boolean fileRead = false;
    protected Hashtable<String, Keyword> keywordHash;

    /**
     * Default constructor.
     */
    public SystemFilter() {
    }

    /**
     * SystemFilter constructor.
     *
     * @param f
     *            MolecularAssembly
     */
    public SystemFilter(MolecularAssembly f) {
        molecularAssembly = f;
    }

    public SystemFilter(MolecularAssembly f, ForceField mm) {
        this(f);
        forceField = mm;
    }

    public SystemFilter(MolecularAssembly f, ForceField mm,
            Hashtable<String,Keyword> keywordHash) {
        this(f,mm);
        this.keywordHash = keywordHash;
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
    public MolecularAssembly getMolecularSystem() {
        return molecularAssembly;
    }

    public FileType getType() {
        return fileType;
    }

    /**
     * This method is different for each subclass and must be overidden
     */
    public abstract boolean readFile();

    public void setFileRead(boolean b) {
        fileRead = b;
    }

    public void setForceField(ForceField f) {
        forceField = f;
    }

    public void setKeywordHash(Hashtable<String,Keyword> keywordHash){
        this.keywordHash = keywordHash;
    }

    public void setMolecularSystem(MolecularAssembly f) {
        molecularAssembly = f;
    }

    public void setType(FileType fileType) {
        this.fileType = fileType;
    }

    /**
     * This method is different for each subclass and must be overidden
     */
    public abstract boolean writeFile();
}
