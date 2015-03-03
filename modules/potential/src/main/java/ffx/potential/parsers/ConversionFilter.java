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

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 * The ConversionFilter class is the base class for most Force Field X parsers for
 * non-Force Field X data structures
 *
 * @author Michael J. Schnieders
 * @autor Jacob M. Litman
 * @since 1.0
 *
 */
public abstract class ConversionFilter {
    
    /**
     * The atomList is filled by filters that extend ConversionFilter.
     */
    protected ArrayList<Atom> atomList = null;
    /**
     * The bondList may be filled by the filters that extend ConversionFilter.
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
     * Structure currently being converted.
     */
    protected Object currentStructure = null;
    /**
     * Append multiple data structures into one MolecularAssembly.
     */
    protected List<Object> structures = null;
    /**
     * Destination file when saved.
     */
    protected File file = null;
    /**
     * Format associated with saving the structure.
     */
    protected Utilities.FileType fileType = Utilities.FileType.UNK;
    /**
     * The file format being handled.
     */
    protected Utilities.DataType dataType = Utilities.DataType.UNK;
    /**
     * Properties associated with this file.
     */
    protected CompositeConfiguration properties;
    /**
     * The molecular mechanics force field being used.
     */
    protected ForceField forceField = null;
    /**
     * True after the data structure has been successfully converted to an FFX
     * data structure.
     */
    protected boolean hasConverted = false;
    /**
     * True if atoms are to be printed to their van der Waals centers instead of
     * nuclear centers (applies primarily to hydrogens).
     */
    protected boolean vdwH = false;
    
    
    public ConversionFilter (Object structure, MolecularAssembly molecularAssembly, 
            ForceField forcefield, CompositeConfiguration properties) {
        this.currentStructure = structure;
        this.structures = new ArrayList<>();
        if (structure != null) {
            structures.add(structure);
        }
        this.activeMolecularAssembly = molecularAssembly;
        this.forceField = forcefield;
        this.properties = properties;
        
        String vdwHydrogens = System.getProperty("vdwHydrogens");
        if (vdwHydrogens != null && vdwHydrogens.equalsIgnoreCase("true")) {
            vdwH = true;
        }
    }
    public ConversionFilter (List<Object> structures, MolecularAssembly molecularAssembly, 
            ForceField forcefield, CompositeConfiguration properties) {
        this.structures = structures;
        if (structures != null && !structures.isEmpty()) {
            this.currentStructure = structures.get(0);
        }
        this.activeMolecularAssembly = molecularAssembly;
        this.forceField = forcefield;
        this.properties = properties;
        
        String vdwHydrogens = System.getProperty("vdwHydrogens");
        if (vdwHydrogens != null && vdwHydrogens.equalsIgnoreCase("true")) {
            vdwH = true;
        }
    }
    public ConversionFilter (Object structure, List<MolecularAssembly> molecularAssemblies, 
            ForceField forcefield, CompositeConfiguration properties) {
        this.currentStructure = structure;
        this.structures = new ArrayList<>();
        if (structure != null) {
            structures.add(structure);
        }
        this.systems = new Vector(molecularAssemblies);
        this.activeMolecularAssembly = systems.firstElement();
        this.forceField = forcefield;
        this.properties = properties;
        
        String vdwHydrogens = System.getProperty("vdwHydrogens");
        if (vdwHydrogens != null && vdwHydrogens.equalsIgnoreCase("true")) {
            vdwH = true;
        }
    }
    
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
     * Returns true if the conversion was successful
     *
     * @return a boolean.
     */
    public boolean converted() {
        return hasConverted;
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
     * @return an array of {@link ffx.potential.MolecularAssembly}
     * objects.
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
     * getFileType</p>
     *
     * @return a {@link ffx.potential.Utilities.FileType} object.
     */
    public Utilities.FileType getFileType() {
        return fileType;
    }
    
    /**
     * <p>
     * getDataType</p>
     *
     * @return a {@link ffx.potential.Utilities.DataType} object.
     */
    public Utilities.DataType getDataType() {
        return dataType;
    }

    /**
     * <p>
     * getFile</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getFile() {
        return file;
    }

    /**
     * This method is different for each subclass and must be overidden
     *
     * @return a boolean.
     */
    public abstract boolean convert();

    /**
     * <p>
     * Setter for the field <code>hasConverted</code>.</p>
     *
     * @param converted a boolean.
     */
    public void setConverted(boolean converted) {
        this.hasConverted = converted;
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
     * setFileType</p>
     *
     * @param fileType a {@link ffx.potential.Utilities.FileType} object.
     */
    public void setFileType(Utilities.FileType fileType) {
        this.fileType = fileType;
    }

    /**
     * <p>
     * setDataType</p>
     *
     * @param dataType a {@link ffx.potential.Utilities.FileType} object.
     */
    public void setDataType(Utilities.DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * <p>
     * setFile</p>
     *
     * @param file a {@link java.io.File} object.
     */
    public void setFile(File file) {
        this.file = file;
    }
    
    /**
     * <p>
     * setStructure</p>
     *
     * @param structure a {@link java.lang.Object} object.
     */
    public void setStructure(Object structure) {
        this.currentStructure = structure;
        this.structures = new ArrayList<>();
        if (structure != null) {
            this.structures.add(structure);
        }
    }

    /**
     * <p>
     * Setter for the field <code>structures</code>.</p>
     *
     * @param structures a {@link java.util.List} object.
     */
    public void setStructures(List<Object> structures) {
        this.structures = structures;
        if (structures != null && !structures.isEmpty()) {
            this.currentStructure = structures.get(0);
        } else {
            this.currentStructure = null;
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
