//******************************************************************************
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
//******************************************************************************
package ffx.potential.parsers;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.parameters.ForceField;

/**
 * The ConversionFilter class is the base class for most Force Field X parsers
 * for non-Force Field X data structures
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public abstract class ConversionFilter {

    private static final Logger logger = Logger.getLogger(ConversionFilter.class.getName());
    private static final Pattern intrangePattern = Pattern.compile("(\\d+)-(\\d+)");
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
    MolecularAssembly activeMolecularAssembly;
    /**
     * All MolecularAssembly instances defined. More than one MolecularAssembly
     * should be defined for PDB entries with alternate locations.
     */
    protected List<MolecularAssembly> systems = new Vector<>();
    /**
     * Structure currently being converted.
     */
    private Object currentStructure;
    /**
     * Append multiple data structures into one MolecularAssembly.
     */
    protected List<Object> structures;
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
    protected ForceField forceField;
    /**
     * True after the data structure has been successfully converted to an FFX
     * data structure.
     */
    private boolean hasConverted = false;
    /**
     * True if atoms are to be printed to their van der Waals centers instead of
     * nuclear centers (applies primarily to hydrogens).
     */
    boolean vdwH;
    /**
     * A set of coordinate restraints obtained from the properties.
     */
    private List<CoordRestraint> coordRestraints;

    /**
     * <p>Constructor for ConversionFilter.</p>
     *
     * @param structure         a {@link java.lang.Object} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forcefield        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    ConversionFilter(Object structure, MolecularAssembly molecularAssembly,
                     ForceField forcefield, CompositeConfiguration properties) {
        this.currentStructure = structure;
        this.structures = new ArrayList<>();
        if (structure != null) {
            structures.add(structure);
        }
        this.activeMolecularAssembly = molecularAssembly;
        this.forceField = forcefield;
        this.properties = properties;

        if (properties != null) {
            vdwH = properties.getBoolean("vdwHydrogens", false);
        } else {
            vdwH = false;
        }
    }

    /**
     * <p>Constructor for ConversionFilter.</p>
     *
     * @param structure           a {@link java.lang.Object} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forcefield          a {@link ffx.potential.parameters.ForceField} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    ConversionFilter(Object structure, List<MolecularAssembly> molecularAssemblies,
                     ForceField forcefield, CompositeConfiguration properties) {
        this.currentStructure = structure;
        this.structures = new ArrayList<>();
        if (structure != null) {
            structures.add(structure);
        }
        this.systems = new ArrayList<>(molecularAssemblies);
        this.activeMolecularAssembly = systems.get(0);
        this.forceField = forcefield;
        this.properties = properties;

        if (properties != null) {
            vdwH = properties.getBoolean("vdwHydrogens", false);
        } else {
            vdwH = false;
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
            MolecularAssembly[] assemblies = new MolecularAssembly[systems.size()];
            return systems.toArray(assemblies);
        } else {
            MolecularAssembly[] assemblies = {activeMolecularAssembly};
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
     * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public void setProperties(CompositeConfiguration properties) {
        this.properties = properties;
    }

    /**
     * <p>
     * setMolecularSystem</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
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
     * Automatically sets atom-specific flags, particularly nouse and inactive,
     * and apply harmonic restraints. Intended to be called at the end of
     * convert() implementations.
     */
    public void applyAtomProperties() {
        /*
          What may be a more elegant implementation is to make convert() a
          public concrete, but skeletal method, and then have convert()
          call a protected abstract readFile method for each implementation.
         */

        Atom[] molaAtoms = activeMolecularAssembly.getAtomArray();
        int nmolaAtoms = molaAtoms.length;
        String[] nouseKeys = properties.getStringArray("nouse");
        for (String nouseKey : nouseKeys) {
            String[] toks = nouseKey.split("\\s+");

            for (String tok : toks) {
                try {
                    int[] nouseRange = parseAtNumArg("restraint", tok, nmolaAtoms);
                    logger.info(format(" Setting atoms %d-%d to not be used",
                            nouseRange[0] + 1, nouseRange[1] + 1));
                    for (int j = nouseRange[0]; j <= nouseRange[1]; j++) {
                        molaAtoms[j].setUse(false);
                    }
                } catch (IllegalArgumentException ex) {
                    boolean atomFound = false;
                    for (Atom atom : molaAtoms) {
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

        String[] inactiveKeys = properties.getStringArray("inactive");
        for (String inactiveKey : inactiveKeys) {
            try {
                int[] inactiveRange = parseAtNumArg("inactive", inactiveKey, nmolaAtoms);
                logger.log(Level.INFO, format(" Atoms %d-%d set to be not "
                        + "active", inactiveRange[0] + 1, inactiveRange[1] + 1));
                for (int i = inactiveRange[0]; i <= inactiveRange[1]; i++) {
                    molaAtoms[i].setActive(false);
                }
            } catch (IllegalArgumentException ex) {
                logger.log(Level.INFO, ex.getLocalizedMessage());
            }
        }

        coordRestraints = new ArrayList<>();
        String[] cRestraintStrings = properties.getStringArray("restraint");
        for (String coordRestraint : cRestraintStrings) {
            String[] toks = coordRestraint.split("\\s+");
            double forceconst = Double.parseDouble(toks[0]);
            logger.info(format(" Adding lambda-disabled coordinate restraint "
                    + "with force constant %10.4f", forceconst));
            Set<Atom> restraintAtoms = new HashSet<>();

            for (int i = 1; i < toks.length; i++) {
                try {
                    int[] nouseRange = parseAtNumArg("restraint", toks[i], nmolaAtoms);
                    logger.info(format(" Adding atoms %d-%d to restraint",
                            nouseRange[0], nouseRange[1]));
                    for (int j = nouseRange[0]; j <= nouseRange[1]; j++) {
                        restraintAtoms.add(molaAtoms[j]);
                    }
                } catch (IllegalArgumentException ex) {
                    boolean atomFound = false;
                    for (Atom atom : molaAtoms) {
                        if (atom.getName().equalsIgnoreCase(toks[i])) {
                            atomFound = true;
                            restraintAtoms.add(atom);
                        }
                    }
                    if (atomFound) {
                        logger.info(format(" Added atoms with name %s to restraint", toks[i]));
                    } else {
                        logger.log(Level.INFO, ex.getLocalizedMessage());
                    }
                }
            }
            if (!restraintAtoms.isEmpty()) {
                Atom[] ats = restraintAtoms.toArray(new Atom[0]);
                coordRestraints.add(new CoordRestraint(ats, forceField, false, forceconst));
            } else {
                logger.warning(format(" Empty or unparseable restraint argument %s", coordRestraint));
            }
        }


        String[] lamRestraintStrings = properties.getStringArray("lamrestraint");
        for (String coordRestraint : lamRestraintStrings) {
            String[] toks = coordRestraint.split("\\s+");
            double forceconst = Double.parseDouble(toks[0]);
            logger.info(format(" Adding lambda-enabled coordinate restraint "
                    + "with force constant %10.4f", forceconst));
            Set<Atom> restraintAtoms = new HashSet<>();

            for (int i = 1; i < toks.length; i++) {
                try {
                    int[] nouseRange = parseAtNumArg("restraint", toks[i], nmolaAtoms);
                    logger.info(format(" Adding atoms %d-%d to restraint",
                            nouseRange[0], nouseRange[1]));
                    for (int j = nouseRange[0]; j <= nouseRange[1]; j++) {
                        restraintAtoms.add(molaAtoms[j]);
                    }
                } catch (IllegalArgumentException ex) {
                    boolean atomFound = false;
                    for (Atom atom : molaAtoms) {
                        if (atom.getName().equalsIgnoreCase(toks[i])) {
                            atomFound = true;
                            restraintAtoms.add(atom);
                        }
                    }
                    if (atomFound) {
                        logger.info(format(" Added atoms with name %s to restraint", toks[i]));
                    } else {
                        logger.log(Level.INFO, ex.getLocalizedMessage());
                    }
                }
            }
            if (!restraintAtoms.isEmpty()) {
                Atom[] ats = restraintAtoms.toArray(new Atom[0]);
                coordRestraints.add(new CoordRestraint(ats, forceField, true, forceconst));
            } else {
                logger.warning(format(" Empty or unparseable restraint argument %s", coordRestraint));
            }
        }

        String[] xyzRestStrings = properties.getStringArray("xyzRestraint");

        // Variables to let sequential and otherwise identical xyzRestraint strings to be globbed into one restraint.
        List<Atom> restraintAts = new ArrayList<>();
        List<double[]> coords = new ArrayList<>();
        double lastForceConst = 0;
        boolean lastUseLam = false;

        for (String xR : xyzRestStrings) {
            String[] toks = xR.split("\\s+");
            int nToks = toks.length;
            if (nToks != 6) {
                logger.info(" XYZ restraint rejected: must have force constant, lambda boolean (true/false), 3 coordinates, and an atom number");
                logger.info(" For a coordinate restraint centered on original coordinates, use restraint or lamrestraint keys.");
                logger.info(format(" Rejected restraint string: %s", xR));
            } else {
                try {
                    double forceConst = Double.parseDouble(toks[0]);
                    boolean useLambda = Boolean.parseBoolean(toks[1]);

                    if (forceConst != lastForceConst || useLambda != lastUseLam) {
                        int nAts = restraintAts.size();
                        if (nAts != 0) {
                            double[][] restXYZ = new double[3][nAts];
                            Atom[] atArr = restraintAts.toArray(new Atom[nAts]);
                            for (int i = 0; i < 3; i++) {
                                for (int j = 0; j < nAts; j++) {
                                    restXYZ[i][j] = coords.get(j)[i];
                                }
                            }
                            CoordRestraint thePin = new CoordRestraint(atArr, forceField, lastUseLam, lastForceConst);
                            thePin.setCoordinatePin(restXYZ);
                            thePin.setIgnoreHydrogen(false);
                            coordRestraints.add(thePin);
                        }
                        restraintAts = new ArrayList<>();
                        coords = new ArrayList<>();
                        lastForceConst = forceConst;
                        lastUseLam = useLambda;
                    }

                    double[] atXYZ = new double[3];
                    for (int i = 0; i < 3; i++) {
                        atXYZ[i] = Double.parseDouble(toks[i + 2]);
                    }
                    int atNum = Integer.parseInt(toks[5]) - 1;
                    restraintAts.add(molaAtoms[atNum]);
                    coords.add(atXYZ);
                } catch (Exception ex) {
                    logger.info(format(" Exception parsing xyzRestraint %s: %s", xR, ex.toString()));
                }
            }
        }

        int nAts = restraintAts.size();
        if (nAts != 0) {
            double[][] restXYZ = new double[3][nAts];
            Atom[] atArr = restraintAts.toArray(new Atom[nAts]);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < nAts; j++) {
                    restXYZ[i][j] = coords.get(j)[i];
                }
            }
            CoordRestraint thePin = new CoordRestraint(atArr, forceField, lastUseLam, lastForceConst);
            thePin.setCoordinatePin(restXYZ);
            thePin.setIgnoreHydrogen(false);
            coordRestraints.add(thePin);
        }

        String[] noElStrings = properties.getStringArray("noElectro");
        for (String noE : noElStrings) {
            String[] toks = noE.split("\\s+");
            for (String tok : toks) {
                try {
                    int[] noERange = parseAtNumArg("noElectro", tok, nmolaAtoms);
                    for (int i = noERange[0]; i <= noERange[1]; i++) {
                        molaAtoms[i].setElectrostatics(false);
                    }
                    logger.log(Level.INFO, format(" Disabled electrostatics "
                            + "for atoms %d-%d", noERange[0] + 1, noERange[1] + 1));
                } catch (IllegalArgumentException ex) {
                    boolean atomFound = false;
                    for (Atom atom : molaAtoms) {
                        if (atom.getName().equalsIgnoreCase(tok)) {
                            atomFound = true;
                            atom.setElectrostatics(false);
                        }
                    }
                    if (atomFound) {
                        logger.info(format(" Disabled electrostatics for atoms with name %s", tok));
                    } else {
                        logger.log(Level.INFO, format(" No electrostatics "
                                + "input %s could not be parsed as a numerical "
                                + "range or atom type present in assembly", tok));
                    }
                }
            }
        }
    }

    /**
     * Gets the coordinate restraints parsed by this Filter.
     *
     * @return Coordinate restraints.
     */
    public List<CoordRestraint> getCoordRestraints() {
        if (!coordRestraints.isEmpty()) {
            return new ArrayList<>(coordRestraints);
        } else {
            return null;
        }
    }

    /**
     * Parses a numerical argument for an atom-specific flag. Intended to reduce
     * the amount of repetitive code in applyAtomProperties by parsing and
     * checking for validity, and then returning the appropriate range. Input
     * should be 1-indexed (user end), output 0-indexed.
     *
     * @param keyType Type of key
     * @param st      Input string
     * @param nAtoms  Number of atoms in the MolecularAssembly
     * @return An int[2] with start, end indices (inclusive).
     * @throws IllegalArgumentException if an invalid argument
     */
    private int[] parseAtNumArg(String keyType, String st, int nAtoms) throws IllegalArgumentException {
        Matcher m = intrangePattern.matcher(st);
        if (m.matches()) {
            int start = Integer.parseInt(m.group(1)) - 1;
            int end = Integer.parseInt(m.group(2)) - 1;
            if (start > end) {
                throw new IllegalArgumentException(format(" %s input %s not "
                        + "valid; start > end", keyType, st));
            } else if (start < 0) {
                throw new IllegalArgumentException(format(" %s input %s not "
                        + "valid; atoms should be indexed starting from 1", keyType, st));
            } else if (start >= nAtoms) {
                throw new IllegalArgumentException(format(" %s input %s not "
                        + "valid; atom range is out of bounds for assembly of "
                        + "length %d", keyType, st, nAtoms));
            } else {
                if (end >= nAtoms) {
                    logger.log(Level.INFO, format(" Truncating range %s "
                            + "to end of valid range %d", st, nAtoms));
                    end = nAtoms - 1;
                }
                int[] indices = {start, end};
                return indices;
            }
        } else {
            try {
                int atNum = Integer.parseUnsignedInt(st) - 1;
                if (atNum < 0 || atNum >= nAtoms) {
                    throw new IllegalArgumentException(format(" %s numerical "
                                    + "argument %s out-of-bounds for range 1 to %d", keyType,
                            st, nAtoms));
                }
                int[] indices = {atNum, atNum};
                return indices;
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException(format(" %s input %s "
                        + "could not be parsed as a positive number or range of "
                        + "positive integers", keyType, st));
            }
        }
    }

    /**
     * This method is different for each subclass and must be overidden.
     * <p>
     * If the append flag is true, "saveFile" will be appended to. Otherwise the
     * default versioning scheme will be applied.
     *
     * @param saveFile a {@link java.io.File} object.
     * @param append   a boolean.
     * @return a boolean.
     */
    public abstract boolean writeFile(File saveFile, boolean append);
}
