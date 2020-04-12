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
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Polymer;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.parameters.ForceField;
import static ffx.utilities.StringUtils.parseAtNumArg;

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
    protected final boolean standardizeAtomNames;
    /**
     * True if atoms are to be printed to their van der Waals centers instead of
     * nuclear centers (applies primarily to hydrogens).
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
     * All MolecularAssembly instances defined. More than one MolecularAssembly
     * should be defined for PDB entries with alternate locations.
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
     * The current MolecularAssembly being populated. Note that more than one
     * MolecularAssembly should be defined for PDB files with alternate
     * locations.
     */
    MolecularAssembly activeMolecularAssembly;
    /**
     * File currently being read.
     */
    File currentFile = null;
    /**
     * A set of coordinate restraints obtained from the properties.
     */
    private List<CoordRestraint> coordRestraints;

    /**
     * Initializations common to the all the constructors.
     *
     * @param forceField The force field.
     * @param properties The CompositeConfiguration properties.
     */
    private SystemFilter(ForceField forceField, CompositeConfiguration properties) {
        this.forceField = forceField;
        this.properties = properties;
        if (properties != null) {
            vdwH = properties.getBoolean("vdwHydrogens", false);
            dieOnMissingAtom = properties.getBoolean("trajectory-dieOnMissing", false);
            standardizeAtomNames = properties.getBoolean("standardizeAtomNames", true);
        } else {
            vdwH = false;
            dieOnMissingAtom = false;
            standardizeAtomNames = true;
        }
    }

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param files             a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}  object.
     * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public SystemFilter(List<File> files, MolecularAssembly molecularAssembly,
                        ForceField forceField, CompositeConfiguration properties) {
        this(forceField, properties);
        this.files = files;
        if (files != null) {
            this.currentFile = files.get(0);
        }
        this.activeMolecularAssembly = molecularAssembly;
    }

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param file              a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}  object.
     * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public SystemFilter(File file, MolecularAssembly molecularAssembly,
                        ForceField forceField, CompositeConfiguration properties) {
        this(forceField, properties);
        files = new ArrayList<>();
        if (file != null) {
            files.add(file);
        }
        this.currentFile = file;
        this.activeMolecularAssembly = molecularAssembly;
    }

    /**
     * <p>
     * Constructor for SystemFilter.</p>
     *
     * @param file                a {@link java.io.File} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forceField          a {@link ffx.potential.parameters.ForceField} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public SystemFilter(File file, List<MolecularAssembly> molecularAssemblies,
                        ForceField forceField, CompositeConfiguration properties) {
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
     * Automatically sets atom-specific flags, particularly nouse and inactive,
     * and apply harmonic restraints. Intended to be called at the end of
     * readFile() implementations.
     */
    public void applyAtomProperties() {
        /*
          What may be a more elegant implementation is to make readFile() a
          public concrete, but skeletal method, and then have readFile() call a
          protected abstract readFile method for each implementation.
        */
        Atom[] molaAtoms = activeMolecularAssembly.getAtomArray();
        int nmolaAtoms = molaAtoms.length;
        String[] nouseKeys = properties.getStringArray("nouse");
        for (String nouseKey : nouseKeys) {
            String[] toks = nouseKey.split("\\s+");
            for (String tok : toks) {
                try {
                    List<Integer> nouseRange = parseAtNumArg("nouse", tok, nmolaAtoms);
                    for (int j : nouseRange) {
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
                List<Integer> inactiveRange = parseAtNumArg("inactive", inactiveKey, nmolaAtoms);
                for (int i : inactiveRange) {
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
            double forceconst;
            try {
                forceconst = Double.parseDouble(toks[0]);
            } catch (NumberFormatException ex) {
                logger.log(Level.INFO, " First argument to coordinate restraint must be a positive force constant; discarding coordinate restraint.");
                continue;
            }
            if (forceconst < 0) {
                logger.log(Level.INFO, " Force constants must be positive. Discarding coordinate restraint.");
                continue;
            }
            logger.info(format(" Adding lambda-disabled coordinate restraint "
                    + "with force constant %10.4f kcal/mol/A", forceconst));
            Set<Atom> restraintAtoms = parseAtomicRanges(activeMolecularAssembly, toks, 1);
            if (!(restraintAtoms == null || restraintAtoms.isEmpty())) {
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
                    + "with force constant %10.4f kcal/mol/A", forceconst));
            Set<Atom> restraintAtoms = new HashSet<>();

            for (int i = 1; i < toks.length; i++) {
                try {
                    List<Integer> restrRange = parseAtNumArg("restraint", toks[i], nmolaAtoms);
                    for (int j : restrRange) {
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
                        logger.log(Level.INFO, format(" Restraint input %s "
                                + "could not be parsed as a numerical range or "
                                + "an atom type present in assembly", toks[i]));
                    }
                }
            }
            if (!restraintAtoms.isEmpty()) {
                Atom[] ats = restraintAtoms.toArray(new Atom[restraintAtoms.size()]);
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
                    List<Integer> noERange = parseAtNumArg("noElectro", tok, nmolaAtoms);
                    for (int i : noERange) {
                        molaAtoms[i].setElectrostatics(false);
                    }
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
                        logger.log(Level.INFO, format(" No electrostatics input %s could not be parsed as a numerical "
                                + "range or atom type present in assembly", tok));
                    }
                }
            }
        }
    }

    /**
     * Attempts to close any open resources associated with the underlying file;
     * primarily to be used when finished reading a trajectory.
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
     * <p>
     * Getter for the field <code>atomList</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<Atom> getAtomList() {
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
     * setFile</p>
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
     * <p>
     * Getter for the field <code>files</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<File> getFiles() {
        return files;
    }

    /**
     * <p>
     * Setter for the field <code>files</code>.</p>
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
     * <p>
     * getMolecularAssemblys</p>
     *
     * @return an array of {@link ffx.potential.MolecularAssembly} objects.
     */
    public MolecularAssembly[] getMolecularAssemblys() {
        if (systems.size() > 0) {
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
     * setType</p>
     *
     * @param fileType a {@link ffx.potential.Utilities.FileType} object.
     */
    public void setType(FileType fileType) {
        this.fileType = fileType;
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
     * This method is different for each subclass and must be overidden
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
     *                   {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public void setProperties(CompositeConfiguration properties) {
        this.properties = properties;
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
                return version(file, (vers == Versioning.PREFIX_ABSOLUTE));
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
    public boolean writeFile(File saveFile, boolean append) {
        return writeFile(saveFile, append, null);
    }

    /**
     * This method is different for each subclass and must be overidden.
     * <p>
     * If the append flag is true, "saveFile" will be appended to. Otherwise the
     * default versioning scheme will be applied.
     *
     * @param saveFile   a {@link java.io.File} object.
     * @param append     a boolean.
     * @param extraLines Additional lines to append to a comments section, or null.
     * @return a boolean.
     */
    public abstract boolean writeFile(File saveFile, boolean append, String[] extraLines);

    private static Set<Atom> parseAtomicRanges(MolecularAssembly mola, String[] toks, int tokOffset) {
        return parseAtomicRanges(mola, Arrays.copyOfRange(toks, tokOffset, toks.length));
    }

    private static Set<Atom> parseAtomicRanges(MolecularAssembly mola, String[] toks) {
        Atom[] atoms = mola.getAtomArray();
        return Arrays.stream(toks).
                parallel().
                flatMap((String tok) -> {
                    if (tok.equalsIgnoreCase("HEAVY")) {
                        return Arrays.stream(mola.getChains()).
                                flatMap((Polymer poly) -> poly.getAtomList().stream().filter(Atom::isHeavy));
                    } else if (tok.matches("^[0-9]+$")) {
                        return Stream.of(atoms[Integer.parseInt(tok) - 1]);
                    } else if (tok.matches("^[0-9]+-[0-9]+")) {
                        String[] subtoks = tok.split("-");
                        int first = Integer.parseInt(subtoks[0]) - 1;
                        int last = Integer.parseInt(subtoks[1]) - 1;
                        return IntStream.rangeClosed(first, last).mapToObj((int i) -> atoms[i]);
                    } else {
                        return Arrays.stream(atoms).filter((Atom a) -> a.getName().equals(tok));
                    }
                }).collect(Collectors.toSet());
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
        String newFn = (prefix)
                ? format("%s_%d.%s", name, number, extension)
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
            newFn = (prefix)
                    ? format("%s_%d.%s", name, number, extension)
                    : format("%s.%s_%d", name, extension, number);
            newFile = new File(newFn);
        }
        return newFile;
    }

    private synchronized static File versionAbsolute(File file, boolean prefix) {
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
        String newFn = (prefix)
                ? format("%s_%d.%s", name, absoluteCounter, extension)
                : format("%s.%s_%d", name, extension, absoluteCounter);
        File newFile = new File(newFn);
        while (newFile.exists()) {
            absoluteCounter++;
            newFn = (prefix)
                    ? format("%s_%d.%s", name, absoluteCounter, extension)
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

    public enum Versioning {
        TINKER, PREFIX, POSTFIX, PREFIX_ABSOLUTE, POSTFIX_ABSOLUTE, NONE;
    }

    /**
     * <p>
     * Setter for the field <code>fileRead</code>.</p>
     *
     * @param fileRead a boolean.
     */
    protected void setFileRead(boolean fileRead) {
        this.fileRead = fileRead;
    }

    /**
     * <p>
     * setMolecularSystem</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     *                          object.
     */
    void setMolecularSystem(MolecularAssembly molecularAssembly) {
        activeMolecularAssembly = molecularAssembly;
    }
}
