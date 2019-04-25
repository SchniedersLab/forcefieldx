//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.numerics.math.VectorMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedUtils;
import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.bonded.MSGroup;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Hybrid36;
import ffx.utilities.StringUtils;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.r;
import static ffx.potential.bonded.AminoAcidUtils.renameArginineHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameAsparagineHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameBetaHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameDeltaHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameEpsilonHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameGammaHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameGlutamineHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameGlycineAlphaHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameIsoleucineHydrogens;
import static ffx.potential.bonded.AminoAcidUtils.renameZetaHydrogens;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static ffx.potential.bonded.NucleicAcidUtils.assignNucleicAcidAtomTypes;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidList;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcid;
import static ffx.potential.bonded.ResidueEnumerations.nucleicAcidList;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_2;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_3;
import static ffx.utilities.StringUtils.padLeft;
import static ffx.utilities.StringUtils.padRight;

/**
 * The PDBFilter class parses data from a Protein DataBank (*.PDB) file. The
 * following records are recognized: ANISOU, ATOM, CONECT, CRYST1, END, HELIX,
 * HETATM, LINK, SHEET, SSBOND, REMARK. The rest are currently ignored.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://www.wwpdb.org/documentation/format32/v3.2.html"> PDB
 * format 3.2</a>
 * @since 1.0
 */
public final class PDBFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(PDBFilter.class.getName());

    /**
     * Map of SEQRES entries.
     */
    private final Map<Character, String[]> seqres = new HashMap<>();
    /**
     * Map of DBREF entries.
     */
    private final Map<Character, int[]> dbref = new HashMap<>();
    /**
     * List of altLoc characters seen in the PDB file.
     */
    private final List<Character> altLocs = new ArrayList<>();
    /**
     * The current altLoc - ie. the one we are defining a chemical system for.
     */
    private Character currentAltLoc = 'A';
    /**
     * List of segIDs defined for the PDB file.
     * <p>
     * The expectation is for chain naming from A-Z, then from 0-9. For large
     * systems, chain names are sometimes reused due to limitations in the PDB
     * format.
     * <p>
     * However, we define segIDs to always be unique. For the first A-Z,0-9
     * series chainID == segID. Then, for second A-Z,0-9 series, the segID =
     * 1A-1Z,10-19, and for the third series segID = 2A-2Z,20-29, and so on.
     */
    private final List<String> segIDs = new ArrayList<>();
    private final Map<Character, List<String>> segidMap = new HashMap<>();
    /**
     * Maps a chain to the number of insertion codes encountered in that chain.
     */
    private final Map<Character, Integer> inscodeCount = new HashMap<>();
    /**
     * Maps chainIDResNumInsCode to renumbered chainIDResNum. For example,
     * residue 52A in chain C might be renumbered to residue 53, and mapped as
     * "C52A" to "C53".
     */
    private final Map<String, String> pdbToNewResMap = new HashMap<>();
    /**
     * List of modified residues *
     */
    private final Map<String, String> modres = new HashMap<>();
    /**
     * Character for the current chain ID.
     */
    private Character currentChainID = null;
    /**
     * String for the current SegID.
     */
    private String currentSegID = null;
    /**
     * Flag to indicate a mutation is requested.
     */
    private boolean mutate = false;
    private List<Mutation> mutations = null;
    /**
     * Flag to indicate if missing fields should be printed (i.e. missing
     * B-factors).
     */
    private boolean printMissingFields = true;
    /**
     * Number of symmetry operators in the current crystal.
     */
    private int nSymOp = 0;
    /**
     * Assume current standard.
     */
    private PDBFileStandard fileStandard = VERSION3_3;
    /**
     * If true, output is directed into arrayOutput instead of the file.
     */
    private boolean listMode = false;
    private ArrayList<String> listOutput = new ArrayList<>();
    /**
     * Don't output atoms which fail Atom.getUse().
     */
    private boolean ignoreUnusedAtoms = false;
    /**
     * Keep track of ATOM record serial numbers to match them with ANISOU
     * records.
     */
    private final HashMap<Integer, Atom> atoms = new HashMap<>();
    /**
     * If false, skip logging "Saving file".
     */
    private boolean logWrites = true;
    /**
     * Keep track of the current MODEL in the file.
     */
    private int modelsRead = 1;
    private final Map<MolecularAssembly, BufferedReader> readers = new HashMap<>();
    /**
     * Tracks output MODEL numbers. Unused if below zero.
     */
    private int modelsWritten = -1;
    private boolean noVersioning = false;
    private final static Set<String> backboneNames;

    private File readFile;

    /**
     * Standardize atom names to PDB standard by default.
     */
    private static final boolean standardizeAtomNames = Boolean.parseBoolean(System.getProperty("standardizeAtomNames", "true"));

    static {
        Set<String> bbAts = new HashSet<>();
        //String[] names = {"C", "CA", "N", "O", "H", "H1", "H2", "H3", "OXT", "OT2"};
        String[] names = {"C", "CA", "N", "O", "OXT", "OT2"};
        bbAts.addAll(Arrays.asList(names));
        backboneNames = Collections.unmodifiableSet(bbAts);
    }

    /**
     * <p>
     * Constructor for PDBFilter.</p>
     *
     * @param files             a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public PDBFilter(List<File> files, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(files, molecularAssembly, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
        readFile = files.get(0);
    }

    /**
     * Parse the PDB File from a URL.
     *
     * @param file              a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public PDBFilter(File file, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssembly, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
        readFile = file;
    }

    /**
     * Parse the PDB File from a URL.
     *
     * @param file                a {@link java.io.File} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forceField          a {@link ffx.potential.parameters.ForceField} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public PDBFilter(File file, List<MolecularAssembly> molecularAssemblies,
                     ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssemblies, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
        readFile = file;
    }

    /**
     * Mutate a residue at the PDB file is being parsed.
     *
     * @param chainID the Chain ID of the residue to mutate.
     * @param resID   the Residue ID of the residue to mutate.
     * @param name    the 3-letter code of the amino acid to mutate to.
     */
    public void mutate(Character chainID, int resID, String name) {
        logger.info(String.format(" Mutating chain %c residue %d to %s.", chainID, resID, name));
        mutate = true;
        if (mutations == null) {
            mutations = new ArrayList<>();
        }
        mutations.add(new Mutation(resID, chainID, name));
    }

    /**
     * <p>mutate.</p>
     *
     * @param mutation a {@link ffx.potential.parsers.PDBFilter.Mutation} object.
     */
    public void mutate(Mutation mutation) {
        mutations.add(mutation);
    }

    /**
     * <p>mutate.</p>
     *
     * @param mutations a {@link java.util.List} object.
     */
    public void mutate(List<Mutation> mutations) {
        this.mutations.addAll(mutations);
    }

    /**
     * <p>
     * clearSegIDs</p>
     */
    public void clearSegIDs() {
        segIDs.clear();
    }

    /**
     * Specify the alternate location.
     *
     * @param molecularAssembly The MolecularAssembly to populate.
     * @param altLoc            The alternate location to use.
     */
    public void setAltID(MolecularAssembly molecularAssembly, Character altLoc) {
        setMolecularSystem(molecularAssembly);
        currentAltLoc = altLoc;
    }

    /**
     * <p>setSymOp.</p>
     *
     * @param symOp a int.
     */
    public void setSymOp(int symOp) {
        this.nSymOp = symOp;
    }

    /**
     * Get the list of alternate locations encountered.
     *
     * @return the alternate location list.
     */
    public List<Character> getAltLocs() {
        return altLocs;
    }

    /**
     * <p>Setter for the field <code>listMode</code>.</p>
     *
     * @param set a boolean.
     */
    public void setListMode(boolean set) {
        listMode = set;
        listOutput = new ArrayList<>();
    }

    /**
     * <p>setModelNumbering.</p>
     *
     * @param set a boolean.
     */
    public void setModelNumbering(boolean set, int modelsWritten) {
        this.modelsWritten = modelsWritten;
    }

    /**
     * <p>Setter for the field <code>noVersioning</code>.</p>
     *
     * @param set a boolean.
     */
    public void setNoVersioning(boolean set) {
        noVersioning = set;
    }

    /**
     * <p>Getter for the field <code>listOutput</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<String> getListOutput() {
        return listOutput;
    }

    /**
     * <p>clearListOutput.</p>
     */
    public void clearListOutput() {
        listOutput.clear();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Parse the PDB File
     */
    @Override
    public boolean readFile() {
        // First atom is #1, to match xyz file format
        int xyzIndex = 1;
        setFileRead(false);
        systems.add(activeMolecularAssembly);

        List<String> conects = new ArrayList<>();
        List<String> links = new ArrayList<>();
        List<String> ssbonds = new ArrayList<>();
        List<String> structs = new ArrayList<>();
        BufferedReader br = null;
        try {
            for (File file : files) {
                currentFile = file;
                if (mutate) {
                    List<Character> chainIDs = new ArrayList<>();
                    try (BufferedReader cr = new BufferedReader(new FileReader(file))) {
                        String line = cr.readLine();
                        while (line != null) {
                            // Replace all tabs w/ 4x spaces
                            line = line.replaceAll("\t", "    ");
                            String identity = line;
                            if (line.length() > 6) {
                                identity = line.substring(0, 6);
                            }
                            identity = identity.trim().toUpperCase();
                            Record record;
                            try {
                                record = Record.valueOf(identity);
                            } catch (Exception e) {
                                /**
                                 * Continue until the record is recognized.
                                 */
                                line = cr.readLine();
                                continue;
                            }
                            switch (record) {
                                case ANISOU:
                                case HETATM:
                                case ATOM:
                                    char c22 = line.charAt(21);
                                    boolean idFound = false;
                                    for (Character chainID : chainIDs) {
                                        if (c22 == chainID) {
                                            idFound = true;
                                            break;
                                        }
                                    }
                                    if (!idFound) {
                                        chainIDs.add(c22);
                                    }
                                    break;
                            }
                            line = cr.readLine();
                        }
                        for (Mutation mtn : mutations) {
                            if (!chainIDs.contains(mtn.chainChar)) {
                                if (chainIDs.size() == 1) {
                                    logger.warning(String.format(" Chain ID %c for "
                                            + "mutation not found: only one chain %c "
                                            + "found.", mtn.chainChar, chainIDs.get(0)));
                                    mtn = new Mutation(mtn.resID, chainIDs.get(0), mtn.resName);
                                } else {
                                    logger.warning(String.format(" Chain ID %c for "
                                            + "mutation not found: mutation will not "
                                            + "proceed.", mtn.chainChar));
                                }
                            }
                        }
                    } catch (IOException ex) {
                        logger.finest(String.format(" Exception %s in parsing file to find chain IDs", ex.toString()));
                    }
                }
                /**
                 * Check that the current file exists and that we can read it.
                 */
                if (currentFile == null || !currentFile.exists() || !currentFile.canRead()) {
                    return false;
                }
                /**
                 * Open the current file for parsing.
                 */
                FileReader fr = new FileReader(currentFile);
                br = new BufferedReader(fr);
                /**
                 * Echo the alternate location being parsed.
                 */
                if (currentAltLoc == 'A') {
                    logger.info(format(" Reading %s", currentFile.getName()));
                } else {
                    logger.info(format(" Reading %s alternate location %s",
                            currentFile.getName(), currentAltLoc));
                }
                /**
                 * Reset the current chain and segID.
                 */
                currentChainID = null;
                currentSegID = null;
                boolean containsInsCode = false;
                /**
                 * Read the first line of the file.
                 */
                String line = br.readLine();
                /**
                 * Parse until END or ENDMDL is found, or to the end of the
                 * file.
                 */
                while (line != null) {
                    // Replace all tabs w/ 4x spaces
                    line = line.replaceAll("\t", "    ");
                    String identity = line;
                    if (line.length() > 6) {
                        identity = line.substring(0, 6);
                    }
                    identity = identity.trim().toUpperCase();
                    Record record;
                    try {
                        record = Record.valueOf(identity);
                    } catch (Exception e) {
                        /**
                         * Continue until the record is recognized.
                         */
                        line = br.readLine();
                        continue;
                    }
                    /**
                     * Switch on the known record.
                     */
                    switch (record) {
                        case ENDMDL:
                        case END:
                            /**
                             * Setting "line" to null will exit the loop.
                             */
                            line = null;
                            continue;
                        case DBREF:
// =============================================================================
//  1 -  6       Record name   "DBREF "
//  8 - 11       IDcode        idCode             ID code of this entry.
// 13            Character     chainID            Chain  identifier.
// 15 - 18       Integer       seqBegin           Initial sequence number of the
//                                                PDB sequence segment.
// 19            AChar         insertBegin        Initial  insertion code of the
//                                                PDB  sequence segment.
// 21 - 24       Integer       seqEnd             Ending sequence number of the
//                                                PDB  sequence segment.
// 25            AChar         insertEnd          Ending insertion code of the
//                                                PDB  sequence segment.
// 27 - 32       LString       database           Sequence database name.
// 34 - 41       LString       dbAccession        Sequence database accession code.
// 43 - 54       LString       dbIdCode           Sequence  database identification code.
// 56 - 60       Integer       dbseqBegin         Initial sequence number of the
//                                                database seqment.
// 61            AChar         idbnsBeg           Insertion code of initial residue of the
//                                                segment, if PDB is the reference.
// 63 - 67       Integer       dbseqEnd           Ending sequence number of the
//                                                database segment.
// 68            AChar         dbinsEnd           Insertion code of the ending residue of
//                                                the segment, if PDB is the reference.
// =============================================================================
                            Character chainID = line.substring(12, 13).toUpperCase().charAt(0);
                            int seqBegin = Integer.parseInt(line.substring(14, 18).trim());
                            int seqEnd = Integer.parseInt(line.substring(20, 24).trim());
                            int[] seqRange = dbref.get(chainID);
                            if (seqRange == null) {
                                seqRange = new int[2];
                                dbref.put(chainID, seqRange);
                            }
                            seqRange[0] = seqBegin;
                            seqRange[1] = seqEnd;
                            break;
                        case SEQRES:
// =============================================================================
//  1 -  6        Record name    "SEQRES"
//  8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
//                                            current  chain. Starts at 1 and increments
//                                            by one  each line. Reset to 1 for each chain.
// 12             Character      chainID      Chain identifier. This may be any single
//                                            legal  character, including a blank which is
//                                            is used if there is only one chain.
// 14 - 17        Integer        numRes       Number of residues in the chain.
//                                            This  value is repeated on every record.
// 20 - 22        Residue name   resName      Residue name.
// 24 - 26        Residue name   resName      Residue name.
// 28 - 30        Residue name   resName      Residue name.
// 32 - 34        Residue name   resName      Residue name.
// 36 - 38        Residue name   resName      Residue name.
// 40 - 42        Residue name   resName      Residue name.
// 44 - 46        Residue name   resName      Residue name.
// 48 - 50        Residue name   resName      Residue name.
// 52 - 54        Residue name   resName      Residue name.
// 56 - 58        Residue name   resName      Residue name.
// 60 - 62        Residue name   resName      Residue name.
// 64 - 66        Residue name   resName      Residue name.
// 68 - 70        Residue name   resName      Residue name.
// =============================================================================
                            activeMolecularAssembly.addHeaderLine(line);
                            chainID = line.substring(11, 12).toUpperCase().charAt(0);
                            int serNum = Integer.parseInt(line.substring(7, 10).trim());
                            String[] chain = seqres.get(chainID);
                            int numRes = Integer.parseInt(line.substring(13, 17).trim());
                            if (chain == null) {
                                chain = new String[numRes];
                                seqres.put(chainID, chain);
                            }
                            int resID = (serNum - 1) * 13;
                            int end = line.length();
                            for (int start = 19; start + 3 <= end; start += 4) {
                                String res = line.substring(start, start + 3).trim();
                                if (res == null || res.length() < 1) {
                                    break;
                                }
                                chain[resID++] = res;
                            }
                            break;
                        case MODRES:
                            String modResName = line.substring(12, 15).trim();
                            String stdName = line.substring(24, 27).trim();
                            modres.put(modResName.toUpperCase(), stdName.toUpperCase());
                            activeMolecularAssembly.addHeaderLine(line);
// =============================================================================
//  1 -  6        Record name     "MODRES"
//  8 - 11        IDcode          idCode         ID code of this entry.
// 13 - 15        Residue name    resName        Residue name used in this entry.
// 17             Character       chainID        Chain identifier.
// 19 - 22        Integer         seqNum         Sequence number.
// 23             AChar           iCode          Insertion code.
// 25 - 27        Residue name    stdRes         Standard residue name.
// 30 - 70        String          comment        Description of the residue modification.
// =============================================================================
                            break;
                        case ANISOU:
// =============================================================================
//  1 - 6        Record name   "ANISOU"
//  7 - 11       Integer       serial         Atom serial number.
// 13 - 16       Atom          name           Atom name.
// 17            Character     altLoc         Alternate location indicator
// 18 - 20       Residue name  resName        Residue name.
// 22            Character     chainID        Chain identifier.
// 23 - 26       Integer       resSeq         Residue sequence number.
// 27            AChar         iCode          Insertion code.
// 29 - 35       Integer       u[0][0]        U(1,1)
// 36 - 42       Integer       u[1][1]        U(2,2)
// 43 - 49       Integer       u[2][2]        U(3,3)
// 50 - 56       Integer       u[0][1]        U(1,2)
// 57 - 63       Integer       u[0][2]        U(1,3)
// 64 - 70       Integer       u[1][2]        U(2,3)
// 77 - 78       LString(2)    element        Element symbol, right-justified.
// 79 - 80       LString(2)    charge         Charge on the atom.
// =============================================================================
                            boolean deleteAnisou = properties.getBoolean("delete-anisou", false);
                            if (deleteAnisou) {
                                break;
                            }
                            Integer serial = Hybrid36.decode(5, line.substring(6, 11));
                            Character altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                            if (!altLocs.contains(altLoc)) {
                                altLocs.add(altLoc);
                            }
                            if (!altLoc.equals(' ') && !altLoc.equals('A')
                                    && !altLoc.equals(currentAltLoc)) {
                                break;
                            }
                            double adp[] = new double[6];
                            adp[0] = Integer.valueOf(line.substring(28, 35).trim()) * 1.0e-4;
                            adp[1] = Integer.valueOf(line.substring(35, 42).trim()) * 1.0e-4;
                            adp[2] = Integer.valueOf(line.substring(42, 49).trim()) * 1.0e-4;
                            adp[3] = Integer.valueOf(line.substring(49, 56).trim()) * 1.0e-4;
                            adp[4] = Integer.valueOf(line.substring(56, 63).trim()) * 1.0e-4;
                            adp[5] = Integer.valueOf(line.substring(63, 70).trim()) * 1.0e-4;
                            if (atoms.containsKey(serial)) {
                                Atom a = atoms.get(serial);
                                a.setAltLoc(altLoc);
                                a.setAnisou(adp);
                            } else {
                                logger.info(format(" No ATOM record for ANISOU serial number %d has been found.", serial));
                                logger.info(format(" This ANISOU record will be ignored:\n %s", line));
                            }
                            break;
                        case ATOM:
// =============================================================================
//  1 -  6        Record name   "ATOM  "
//  7 - 11        Integer       serial       Atom serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
// =============================================================================
                            String name;
                            String resName;
                            String segID;
                            int resSeq;
                            boolean printAtom;
                            double d[];
                            double occupancy;
                            double tempFactor;
                            Atom newAtom;
                            Atom returnedAtom;
                            // If it's a misnamed water, it will fall through to HETATM.
                            if (!line.substring(17, 20).trim().equals("HOH")) {
                                serial = Hybrid36.decode(5, line.substring(6, 11));
                                name = line.substring(12, 16).trim();
                                if (name.toUpperCase().contains("1H") || name.toUpperCase().contains("2H")
                                        || name.toUpperCase().contains("3H")) {
                                    // VERSION3_2 is presently just a placeholder for "anything non-standard".
                                    fileStandard = VERSION3_2;
                                }
                                altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                                if (!altLocs.contains(altLoc)) {
                                    altLocs.add(altLoc);
                                }
                                if (!altLoc.equals(' ') && !altLoc.equals('A')
                                        && !altLoc.equals(currentAltLoc)) {
                                    break;
                                }
                                resName = line.substring(17, 20).trim();
                                chainID = line.substring(21, 22).charAt(0);
                                segID = getSegID(chainID);
                                resSeq = Hybrid36.decode(4, line.substring(22, 26));

                                char insertionCode = line.charAt(26);
                                if (insertionCode != ' ' && !containsInsCode) {
                                    containsInsCode = true;
                                    logger.warning(" FFX support for files with "
                                            + "insertion codes is experimental. "
                                            + "Residues will be renumbered to "
                                            + "eliminate insertion codes (52A "
                                            + "becomes 53, 53 becomes 54, etc)");
                                }

                                int offset = inscodeCount.getOrDefault(chainID, 0);
                                String pdbResNum = String.format("%c%d%c", chainID, resSeq, insertionCode);
                                if (!pdbToNewResMap.containsKey(pdbResNum)) {
                                    if (insertionCode != ' ') {
                                        ++offset;
                                        inscodeCount.put(chainID, offset);
                                    }
                                    resSeq += offset;
                                    if (offset != 0) {
                                        logger.info(String.format(" Chain %c "
                                                        + "residue %s-%s renumbered to %c %s-%d",
                                                chainID, pdbResNum.substring(1).trim(),
                                                resName, chainID, resName, resSeq));
                                    }
                                    String newNum = String.format("%c%d", chainID, resSeq);
                                    pdbToNewResMap.put(pdbResNum, newNum);
                                } else {
                                    resSeq += offset;
                                }

                                printAtom = false;
                                if (mutate) {
                                    boolean doBreak = false;
                                    for (Mutation mtn : mutations) {
                                        if (chainID == mtn.chainChar
                                                && resSeq == mtn.resID) {
                                            String atomName = name.toUpperCase();
                                            /*if (atomName.equals("N") || atomName.equals("C")
                                                    || atomName.equals("O") || atomName.equals("CA")) {*/
                                            if (backboneNames.contains(atomName)) {
                                                printAtom = true;
                                                resName = mtn.resName;
                                            } else {
                                                logger.info(String.format(" Deleting atom %s of %s %d",
                                                        atomName, resName, resSeq));
                                                doBreak = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (doBreak) {
                                        break;
                                    }
                                }
                                d = new double[3];
                                d[0] = Double.valueOf(line.substring(30, 38).trim());
                                d[1] = Double.valueOf(line.substring(38, 46).trim());
                                d[2] = Double.valueOf(line.substring(46, 54).trim());
                                occupancy = 1.0;
                                tempFactor = 1.0;
                                try {
                                    occupancy = Double.valueOf(line.substring(54, 60).trim());
                                    tempFactor = Double.valueOf(line.substring(60, 66).trim());
                                } catch (NumberFormatException | StringIndexOutOfBoundsException e) {
                                    // Use default values.
                                    if (printMissingFields) {
                                        logger.info(" Missing occupancy and b-factors set to 1.0.");
                                        printMissingFields = false;
                                    } else if (logger.isLoggable(Level.FINE)) {
                                        logger.fine(" Missing occupancy and b-factors set to 1.0.");
                                    }
                                }
                                newAtom = new Atom(0, name, altLoc, d, resName, resSeq,
                                        chainID, occupancy, tempFactor, segID);
                                // Check if this is a modified residue.
                                if (modres.containsKey(resName.toUpperCase())) {
                                    newAtom.setModRes(true);
                                }

                                returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
                                if (returnedAtom != newAtom) {
                                    // A previously added atom has been retained.
                                    atoms.put(serial, returnedAtom);
                                    if (logger.isLoggable(Level.FINE)) {
                                        logger.fine(returnedAtom + " has been retained over\n" + newAtom);
                                    }
                                } else {
                                    // The new atom has been added.
                                    atoms.put(serial, newAtom);
                                    // Check if the newAtom took the xyzIndex of a previous alternate conformer.
                                    if (newAtom.getIndex() == 0) {
                                        newAtom.setXyzIndex(xyzIndex++);
                                    }
                                    if (printAtom) {
                                        logger.info(newAtom.toString());
                                    }
                                }
                                break;
                            }
                        case HETATM:
// =============================================================================
//  1 - 6        Record name    "HETATM"
//  7 - 11       Integer        serial        Atom serial number.
// 13 - 16       Atom           name          Atom name.
// 17            Character      altLoc        Alternate location indicator.
// 18 - 20       Residue name   resName       Residue name.
// 22            Character      chainID       Chain identifier.
// 23 - 26       Integer        resSeq        Residue sequence number.
// 27            AChar          iCode         Code for insertion of residues.
// 31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
// 39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
// 47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
// 55 - 60       Real(6.2)      occupancy     Occupancy.
// 61 - 66       Real(6.2)      tempFactor    Temperature factor.
// 77 - 78       LString(2)     element       Element symbol; right-justified.
// 79 - 80       LString(2)     charge        Charge on the atom.
// =============================================================================
                            serial = Hybrid36.decode(5, line.substring(6, 11));
                            name = line.substring(12, 16).trim();
                            altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                            if (!altLocs.contains(altLoc)) {
                                altLocs.add(altLoc);
                            }
                            if (!altLoc.equals(' ') && !altLoc.equals('A')
                                    && !altLoc.equals(currentAltLoc)) {
                                break;
                            }
                            resName = line.substring(17, 20).trim();
                            chainID = line.substring(21, 22).charAt(0);
                            segID = getSegID(chainID);
                            resSeq = Hybrid36.decode(4, line.substring(22, 26));

                            char insertionCode = line.charAt(26);
                            if (insertionCode != ' ' && !containsInsCode) {
                                containsInsCode = true;
                                logger.warning(" FFX support for files with "
                                        + "insertion codes is experimental. "
                                        + "Residues will be renumbered to "
                                        + "eliminate insertion codes (52A "
                                        + "becomes 53, 53 becomes 54, etc)");
                            }

                            int offset = inscodeCount.getOrDefault(chainID, 0);
                            String pdbResNum = String.format("%c%d%c", chainID, resSeq, insertionCode);
                            if (!pdbToNewResMap.containsKey(pdbResNum)) {
                                if (insertionCode != ' ') {
                                    ++offset;
                                    inscodeCount.put(chainID, offset);
                                }
                                resSeq += offset;
                                if (offset != 0) {
                                    logger.info(String.format(" Chain %c "
                                                    + "molecule %s-%s renumbered to %c %s-%d",
                                            chainID, pdbResNum.substring(1).trim(),
                                            resName, chainID, resName, resSeq));
                                }
                                String newNum = String.format("%c%d", chainID, resSeq);
                                pdbToNewResMap.put(pdbResNum, newNum);
                            } else {
                                resSeq += offset;
                            }

                            d = new double[3];
                            d[0] = Double.valueOf(line.substring(30, 38).trim());
                            d[1] = Double.valueOf(line.substring(38, 46).trim());
                            d[2] = Double.valueOf(line.substring(46, 54).trim());
                            occupancy = 1.0;
                            tempFactor = 1.0;
                            try {
                                occupancy = Double.valueOf(line.substring(54, 60).trim());
                                tempFactor = Double.valueOf(line.substring(60, 66).trim());
                            } catch (NumberFormatException | StringIndexOutOfBoundsException e) {
                                // Use default values.
                                if (printMissingFields) {
                                    logger.info(" Missing occupancy and b-factors set to 1.0.");
                                    printMissingFields = false;
                                } else if (logger.isLoggable(Level.FINE)) {
                                    logger.fine(" Missing occupancy and b-factors set to 1.0.");
                                }
                            }
                            newAtom = new Atom(0, name, altLoc, d, resName, resSeq, chainID,
                                    occupancy, tempFactor, segID);
                            newAtom.setHetero(true);
                            // Check if this is a modified residue.
                            if (modres.containsKey(resName.toUpperCase())) {
                                newAtom.setModRes(true);
                            }
                            returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
                            if (returnedAtom != newAtom) {
                                // A previously added atom has been retained.
                                atoms.put(serial, returnedAtom);
                                if (logger.isLoggable(Level.FINE)) {
                                    logger.fine(returnedAtom + " has been retained over\n" + newAtom);
                                }
                            } else {
                                // The new atom has been added.
                                atoms.put(serial, newAtom);
                                newAtom.setXyzIndex(xyzIndex++);
                            }
                            break;
                        case CRYST1:
// =============================================================================
// The CRYST1 record presents the unit cell parameters, space group, and Z
// value. If the structure was not determined by crystallographic means, CRYST1
// simply provides the unitary values, with an appropriate REMARK.
//
//  7 - 15       Real(9.3)     a              a (Angstroms).
// 16 - 24       Real(9.3)     b              b (Angstroms).
// 25 - 33       Real(9.3)     c              c (Angstroms).
// 34 - 40       Real(7.2)     alpha          alpha (degrees).
// 41 - 47       Real(7.2)     beta           beta (degrees).
// 48 - 54       Real(7.2)     gamma          gamma (degrees).
// 56 - 66       LString       sGroup         Space  group.
// 67 - 70       Integer       z              Z value.
// =============================================================================
                            double aaxis = Double.valueOf(line.substring(6, 15).trim());
                            double baxis = Double.valueOf(line.substring(15, 24).trim());
                            double caxis = Double.valueOf(line.substring(24, 33).trim());
                            double alpha = Double.valueOf(line.substring(33, 40).trim());
                            double beta = Double.valueOf(line.substring(40, 47).trim());
                            double gamma = Double.valueOf(line.substring(47, 54).trim());
                            int limit = 66;
                            if (line.length() < 66) {
                                limit = line.length();
                            }
                            String sg = line.substring(55, limit).trim();
                            properties.addProperty("a-axis", aaxis);
                            properties.addProperty("b-axis", baxis);
                            properties.addProperty("c-axis", caxis);
                            properties.addProperty("alpha", alpha);
                            properties.addProperty("beta", beta);
                            properties.addProperty("gamma", gamma);
                            properties.addProperty("spacegroup", SpaceGroup.pdb2ShortName(sg));
                            break;
                        case CONECT:
// =============================================================================
//  7 - 11        Integer        serial       Atom  serial number
// 12 - 16        Integer        serial       Serial number of bonded atom
// 17 - 21        Integer        serial       Serial number of bonded atom
// 22 - 26        Integer        serial       Serial number of bonded atom
// 27 - 31        Integer        serial       Serial number of bonded atom
//
// CONECT records involving atoms for which the coordinates are not present
// in the entry (e.g., symmetry-generated) are not given.
// CONECT records involving atoms for which the coordinates are missing due
// to disorder, are also not provided.
// =============================================================================
                            conects.add(line);
                            break;
                        case LINK:
// =============================================================================
// The LINK records specify connectivity between residues that is not implied by
// the primary structure. Connectivity is expressed in terms of the atom names.
// They also include the distance associated with the each linkage following the
// symmetry operations at the end of each record.
// 13 - 16         Atom           name1           Atom name.
// 17              Character      altLoc1         Alternate location indicator.
// 18 - 20         Residue name   resName1        Residue  name.
// 22              Character      chainID1        Chain identifier.
// 23 - 26         Integer        resSeq1         Residue sequence number.
// 27              AChar          iCode1          Insertion code.
// 43 - 46         Atom           name2           Atom name.
// 47              Character      altLoc2         Alternate location indicator.
// 48 - 50         Residue name   resName2        Residue name.
// 52              Character      chainID2        Chain identifier.
// 53 - 56         Integer        resSeq2         Residue sequence number.
// 57              AChar          iCode2          Insertion code.
// 60 - 65         SymOP          sym1            Symmetry operator atom 1.
// 67 - 72         SymOP          sym2            Symmetry operator atom 2.
// 74 – 78         Real(5.2)      Length          Link distance
// =============================================================================
                            Character a1 = line.charAt(16);
                            Character a2 = line.charAt(46);
                            if (a1 != a2) {
                                // logger.info(format(" Ignoring LINK record as alternate locations do not match\n %s.", line));
                                break;
                            }
                            if (currentAltLoc == 'A') {
                                if ((a1 == ' ' || a1 == 'A')
                                        && (a2 == ' ' || a2 == 'A')) {
                                    links.add(line);
                                }
                            } else if (a1 == currentAltLoc && a2 == currentAltLoc) {
                                links.add(line);
                            }
                            break;
                        case SSBOND:
// =============================================================================
// The SSBOND record identifies each disulfide bond in protein and polypeptide
// structures by identifying the two residues involved in the bond.
// The disulfide bond distance is included after the symmetry operations at
// the end of the SSBOND record.
//
//  8 - 10        Integer         serNum       Serial number.
// 12 - 14        LString(3)      "CYS"        Residue name.
// 16             Character       chainID1     Chain identifier.
// 18 - 21        Integer         seqNum1      Residue sequence number.
// 22             AChar           icode1       Insertion code.
// 26 - 28        LString(3)      "CYS"        Residue name.
// 30             Character       chainID2     Chain identifier.
// 32 - 35        Integer         seqNum2      Residue sequence number.
// 36             AChar           icode2       Insertion code.
// 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
// 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
// 74 – 78        Real(5.2)      Length        Disulfide bond distance
//
// If SG of cysteine is disordered then there are possible alternate linkages.
// wwPDB practice is to put together all possible SSBOND records. This is
// problematic because the alternate location identifier is not specified in
// the SSBOND record.
//
// Notes:
// SSBOND records may be invalid if chain IDs are reused.
// SSBOND records are applied by FFX to the A conformer (not alternate conformers).
// =============================================================================
                            if (currentAltLoc == 'A') {
                                ssbonds.add(line);
                            }
                            break;
                        case HELIX:
// =============================================================================
// HELIX records are used to identify the position of helices in the molecule.
// Helices are named, numbered, and classified by type. The residues where the
// helix begins and ends are noted, as well as the total length.
//
//  8 - 10        Integer        serNum        Serial number of the helix. This starts
//                                             at 1  and increases incrementally.
// 12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
//                                             number, each helix is given an
//                                             alphanumeric character helix identifier.
// 16 - 18        Residue name   initResName   Name of the initial residue.
// 20             Character      initChainID   Chain identifier for the chain containing
//                                             this  helix.
// 22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
// 26             AChar          initICode     Insertion code of the initial residue.
// 28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
// 32             Character      endChainID    Chain identifier for the chain containing
//                                             this  helix.
// 34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
// 38             AChar          endICode      Insertion code of the terminal residue.
// 39 - 40        Integer        helixClass    Helix class (see below).
// 41 - 70        String         comment       Comment about this helix.
// 72 - 76        Integer        length        Length of this helix.
//
//                                      CLASS NUMBER
// TYPE OF  HELIX                     (COLUMNS 39 - 40)
// --------------------------------------------------------------
// Right-handed alpha (default)                1
// Right-handed omega                          2
// Right-handed pi                             3
// Right-handed gamma                          4
// Right-handed 3 - 10                         5
// Left-handed alpha                           6
// Left-handed omega                           7
// Left-handed gamma                           8
// 2 - 7 ribbon/helix                          9
// Polyproline                                10
// =============================================================================
                        case SHEET:
// =============================================================================
// SHEET records are used to identify the position of sheets in the molecule.
// Sheets are both named and numbered. The residues where the sheet begins and
// ends are noted.
//
//  8 - 10        Integer       strand         Strand  number which starts at 1 for each
//                                             strand within a sheet and increases by one.
// 12 - 14        LString(3)    sheetID        Sheet  identifier.
// 15 - 16        Integer       numStrands     Number  of strands in sheet.
// 18 - 20        Residue name  initResName    Residue  name of initial residue.
// 22             Character     initChainID    Chain identifier of initial residue
//                                             in strand.
// 23 - 26        Integer       initSeqNum     Sequence number of initial residue
//                                             in strand.
// 27             AChar         initICode      Insertion code of initial residue
//                                             in  strand.
// 29 - 31        Residue name  endResName     Residue name of terminal residue.
// 33             Character     endChainID     Chain identifier of terminal residue.
// 34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
// 38             AChar         endICode       Insertion code of terminal residue.
// 39 - 40        Integer       sense          Sense of strand with respect to previous
//                                             strand in the sheet. 0 if first strand,
//                                             1 if  parallel,and -1 if anti-parallel.
// 42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
// 46 - 48        Residue name  curResName     Registration.  Residue name in current strand
// 50             Character     curChainId     Registration. Chain identifier in
//                                             current strand.
// 51 - 54        Integer       curResSeq      Registration.  Residue sequence number
//                                             in current strand.
// 55             AChar         curICode       Registration. Insertion code in
//                                             current strand.
// 57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
// 61 - 63        Residue name  prevResName    Registration.  Residue name in
//                                             previous strand.
// 65             Character     prevChainId    Registration.  Chain identifier in
//                                             previous  strand.
// 66 - 69        Integer       prevResSeq     Registration. Residue sequence number
//                                             in previous strand.
// 70             AChar         prevICode      Registration.  Insertion code in
//                                             previous strand.
// =============================================================================
                            structs.add(line);
                            break;
                        case MODEL: // Currently, no handling in initial read.
                            break;
                        case MTRIX1:
// ================================================================================
// MTRIXn (n = 1, 2, or 3) records present transformations expressing non-crystallographic symmetry.
// MTRIXn will appear only when such transformations are required to generate an entire asymmetric unit,
// such as a large viral structure.
//
//  8 - 10        Integer       serial         Serial number.
// 11 - 20        Real(10.6)    m[n][1]        Mn1
// 21 - 30        Real(10.6)    m[n][2]        Mn2
// 31 - 40        Real(10.6)    m[n][3]        Mn3
// 21 - 30        Real(10.6)    v[n]           Vn
// 60             Integer       iGiven         1 if coordinates for the representations which are
//                                              approximately related by the transformations of the
//                                              molecule are contained in the entry. Otherwise, blank.
// =================================================================================
                            StringBuilder MTRX1 = new StringBuilder(line.substring(11, 55));
                            properties.addProperty("MTRIX1", MTRX1);
                            break;
                        case MTRIX2:
                            StringBuilder MTRX2 = new StringBuilder(line.substring(11, 55));
                            properties.addProperty("MTRIX2", MTRX2);
                            break;
                        case MTRIX3:
                            StringBuilder MTRX3 = new StringBuilder(line.substring(11, 55));
                            properties.addProperty("MTRIX3", MTRX3);
                            break;
                        case REMARK:
// =================================================================================
// REMARK 350: presents all transformations, both crystallographic and non-crystallographic, needed to
// generate the biomolecule. These transformations operate on the coordinates in the entry. Both author
// and computational descriptions of assemblies are provided, if applicable. For strict ncs case where
// more than one assembly presents in asymmetric unit, only one chain with unit matrix will reported in
// REMARK 350, the other chain will be generated by rotation and translation.
//
// 20 - 23        Integer       serial         Serial number.
// 24 - 33        Real(10.6)    m[n][1]        Mn1
// 34 - 43        Real(10.6)    m[n][2]        Mn2
// 44 - 53        Real(10.6)    m[n][3]        Mn3
// 59 - 68        Real(10.6)    v[n]           Vn
// =================================================================================
                            if (line.length() >= 68) {
                                String remarkType = line.substring(7, 10).trim();
                                if (remarkType.matches("\\d+") && Integer.parseInt(remarkType) == 350 && line.substring(13, 18).equalsIgnoreCase("BIOMT")) {
                                    properties.addProperty("BIOMTn", new StringBuilder(line.substring(24, 68)));
                                }
                            }
                            break;
                        default:
                            break;
                    }
                    line = br.readLine();
                }
                br.close();
            }
            xyzIndex--;
            setFileRead(true);
        } catch (IOException e) {
            logger.exiting(PDBFilter.class.getName(), "readFile", e);
            return false;
        }

        /**
         * Locate disulfide bonds; bond parameters are assigned below.
         */
        List<Bond> ssBondList = locateDisulfideBonds(ssbonds);

        /**
         * Record the number of atoms read in from the PDB file before applying
         * algorithms that may build new atoms.
         */
        int pdbAtoms = activeMolecularAssembly.getAtomArray().length;

        /**
         * Build missing backbone atoms in loops.
         */
        buildMissingResidues(xyzIndex);

        /**
         * Assign atom types. Missing side-chains atoms and missing hydrogens
         * will be built in.
         */
        assignAtomTypes();

        /**
         * Assign disulfide bonds parameters and log their creation.
         */
        buildDisulfideBonds(ssBondList);

        /**
         * Finally, re-number the atoms if missing atoms were created.
         */
        if (pdbAtoms != activeMolecularAssembly.getAtomArray().length) {
            numberAtoms();
        }

        return true;
    }

    /**
     * <p>setIgnoreInactiveAtoms.</p>
     *
     * @param ignoreInactiveAtoms a boolean.
     */
    public void setIgnoreInactiveAtoms(boolean ignoreInactiveAtoms) {
        this.ignoreUnusedAtoms = ignoreInactiveAtoms;
    }

    /**
     * Sets whether this PDBFilter should log each time it saves to a file.
     *
     * @param logWrites a boolean.
     */
    public void setLogWrites(boolean logWrites) {
        this.logWrites = logWrites;
    }

    /**
     * <p>writeFileWithHeader.</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header   a {@link java.lang.String} object.
     * @param append   a boolean.
     * @return a boolean.
     */
    public boolean writeFileWithHeader(File saveFile, String header, boolean append) {
        FileWriter fw;
        BufferedWriter bw;
        if (standardizeAtomNames) {
            renameAtomsToPDBStandard(activeMolecularAssembly);
        }

        try {
            File newFile = saveFile;
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (listMode) {
                listOutput.add(header);
            } else {
                fw = new FileWriter(newFile, append);
                bw = new BufferedWriter(fw);
                bw.write(header);
                bw.newLine();
                bw.close();
            }
        } catch (Exception e) {
            String message = "Exception writing to file: " + saveFile.toString();
            logger.log(Level.WARNING, message, e);
            return false;
        }
        if (writeFile(saveFile, true)) {
            logger.log(Level.INFO, " Wrote PDB to file {0}", saveFile.getPath());
            return true;
        } else {
            logger.log(Level.INFO, " Error writing to file {0}", saveFile.getPath());
            return false;
        }
    }

    /**
     * <p>writeFileWithHeader.</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header   a {@link java.lang.String} object.
     * @return a boolean.
     */
    public boolean writeFileWithHeader(File saveFile, String header) {
        return writeFileWithHeader(saveFile, header, true);
    }

    /**
     * <p>writeFileWithHeader.</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header   a {@link java.lang.StringBuilder} object.
     * @return a boolean.
     */
    public boolean writeFileWithHeader(File saveFile, StringBuilder header) {
        return writeFileWithHeader(saveFile, header.toString());
    }

    /**
     * <p>
     * writeFile</p>
     *
     * @param saveFile    a {@link java.io.File} object.
     * @param append      Whether to append to saveFile (vs over-write).
     * @param printLinear Whether to print atoms linearly or by element.
     * @param writeEnd    True if this is the final model.
     * @return Success of writing.
     */
    public boolean writeFile(File saveFile, boolean append, boolean printLinear, boolean writeEnd) {
        return writeFile(saveFile, append, printLinear, Collections.emptySet(), writeEnd);
    }

    /**
     * <p>
     * writeFile</p>
     *
     * @param saveFile    a {@link java.io.File} object to save to.
     * @param append      Whether to append to saveFile (vs over-write).
     * @param printLinear Whether to print atoms linearly or by element.
     * @param toExclude   A {@link java.util.Set} of {@link ffx.potential.bonded.Atom}s to exclude from writing.
     * @param writeEnd    True if this is the final model.
     * @return Success of writing.
     */
    public boolean writeFile(File saveFile, boolean append, boolean printLinear, Set<Atom> toExclude, boolean writeEnd) {
        if (standardizeAtomNames) {
            renameAtomsToPDBStandard(activeMolecularAssembly);
        }

        final Set<Atom> atomExclusions = toExclude == null ? Collections.emptySet() : toExclude;

        if (saveFile == null) {
            return false;
        }

        if (vdwH) {
            logger.info(" Printing hydrogens to van der Waals centers instead of nuclear locations.");
        }

        if (nSymOp != 0) {
            logger.info(String.format(" Printing atoms with symmetry operator %s\n",
                    activeMolecularAssembly.getCrystal().spaceGroup.getSymOp(nSymOp).toString()));
        }

        /**
         * Create StringBuilders for ATOM, ANISOU and TER records that can be
         * reused.
         */
        StringBuilder sb = new StringBuilder("ATOM  ");
        StringBuilder anisouSB = new StringBuilder("ANISOU");
        StringBuilder terSB = new StringBuilder("TER   ");
        StringBuilder model = null;
        for (int i = 6; i < 80; i++) {
            sb.append(' ');
            anisouSB.append(' ');
            terSB.append(' ');
        }
        FileWriter fw;
        BufferedWriter bw;
        try {
            File newFile = saveFile;
            if (!append) {
                if (!noVersioning) {
                    newFile = version(saveFile);
                }
            } else if (modelsWritten >= 0) {
                model = new StringBuilder(String.format("MODEL     %-4d", ++modelsWritten));
                for (int i = 15; i < 80; i++) {
                    model.append(' ');
                }
            }
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (logWrites) {
                logger.log(Level.INFO, " Saving {0}", activeMolecularAssembly.getName());
            }
            fw = new FileWriter(newFile, append);
            bw = new BufferedWriter(fw);
            /**
             * Will come before CRYST1 and ATOM records, but after anything
             * written by writeFileWithHeader (particularly X-ray refinement
             * statistics).
             */
            String[] headerLines = activeMolecularAssembly.getHeaderLines();
            for (String line : headerLines) {
                bw.write(String.format("%s\n", line));
            }
            if (model != null) {
                if (!listMode) {
                    bw.write(model.toString());
                    bw.newLine();
                } else {
                    listOutput.add(model.toString());
                }
            }
// =============================================================================
// The CRYST1 record presents the unit cell parameters, space group, and Z
// value. If the structure was not determined by crystallographic means, CRYST1
// simply provides the unitary values, with an appropriate REMARK.
//
//  7 - 15       Real(9.3)     a              a (Angstroms).
// 16 - 24       Real(9.3)     b              b (Angstroms).
// 25 - 33       Real(9.3)     c              c (Angstroms).
// 34 - 40       Real(7.2)     alpha          alpha (degrees).
// 41 - 47       Real(7.2)     beta           beta (degrees).
// 48 - 54       Real(7.2)     gamma          gamma (degrees).
// 56 - 66       LString       sGroup         Space  group.
// 67 - 70       Integer       z              Z value.
// =============================================================================
            Crystal crystal = activeMolecularAssembly.getCrystal();
            if (crystal != null && !crystal.aperiodic()) {
                Crystal c = crystal.getUnitCell();
                if (!listMode) {
                    bw.write(format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s\n", c.a, c.b, c.c, c.alpha, c.beta,
                            c.gamma, padRight(c.spaceGroup.pdbName, 10)));
                } else {
                    listOutput.add(format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s", c.a, c.b, c.c, c.alpha, c.beta,
                            c.gamma, padRight(c.spaceGroup.pdbName, 10)));
                }
            }
// =============================================================================
// The SSBOND record identifies each disulfide bond in protein and polypeptide
// structures by identifying the two residues involved in the bond.
// The disulfide bond distance is included after the symmetry operations at
// the end of the SSBOND record.
//
//  8 - 10        Integer         serNum       Serial number.
// 12 - 14        LString(3)      "CYS"        Residue name.
// 16             Character       chainID1     Chain identifier.
// 18 - 21        Integer         seqNum1      Residue sequence number.
// 22             AChar           icode1       Insertion code.
// 26 - 28        LString(3)      "CYS"        Residue name.
// 30             Character       chainID2     Chain identifier.
// 32 - 35        Integer         seqNum2      Residue sequence number.
// 36             AChar           icode2       Insertion code.
// 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
// 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
// 74 – 78        Real(5.2)      Length        Disulfide bond distance
//
// If SG of cysteine is disordered then there are possible alternate linkages.
// wwPDB practice is to put together all possible SSBOND records. This is
// problematic because the alternate location identifier is not specified in
// the SSBOND record.
// =============================================================================
            int serNum = 1;
            Polymer[] polymers = activeMolecularAssembly.getChains();
            if (polymers != null) {
                for (Polymer polymer : polymers) {
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        if (residue.getName().equalsIgnoreCase("CYS")) {
                            List<Atom> cysAtoms = residue.getAtomList().stream().
                                    filter(a -> !atomExclusions.contains(a)).
                                    collect(Collectors.toList());
                            Atom SG1 = null;
                            for (Atom atom : cysAtoms) {
                                String atName = atom.getName().toUpperCase();
                                if (atName.equals("SG") || atName.equals("SH") || atom.getAtomType().atomicNumber == 16) {
                                    SG1 = atom;
                                    break;
                                }
                            }
                            List<Bond> bonds = SG1.getBonds();
                            for (Bond bond : bonds) {
                                Atom SG2 = bond.get1_2(SG1);
                                if (SG2.getAtomType().atomicNumber == 16 && !atomExclusions.contains(SG2)) {
                                    if (SG1.getIndex() < SG2.getIndex()) {
                                        bond.energy(false);
                                        if (!listMode) {
                                            bw.write(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n",
                                                    serNum++,
                                                    SG1.getChainID().toString(), Hybrid36.encode(4, SG1.getResidueNumber()),
                                                    SG2.getChainID().toString(), Hybrid36.encode(4, SG2.getResidueNumber()),
                                                    "", bond.getValue()));
                                        } else {
                                            listOutput.add(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n",
                                                    serNum++,
                                                    SG1.getChainID().toString(), Hybrid36.encode(4, SG1.getResidueNumber()),
                                                    SG2.getChainID().toString(), Hybrid36.encode(4, SG2.getResidueNumber()),
                                                    "", bond.getValue()));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
// =============================================================================
//
//  7 - 11        Integer       serial       Atom serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
// =============================================================================
//         1         2         3         4         5         6         7
//123456789012345678901234567890123456789012345678901234567890123456789012345678
//ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
//ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
            MolecularAssembly molecularAssemblies[] = this.getMolecularAssemblys();
            int serial = 1;
            // Loop over biomolecular chains
            if (polymers != null) {
                for (Polymer polymer : polymers) {
                    currentSegID = polymer.getName();
                    currentChainID = polymer.getChainID();
                    sb.setCharAt(21, currentChainID);
                    // Loop over residues
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        String resName = residue.getName();
                        if (resName.length() > 3) {
                            resName = resName.substring(0, 3);
                        }
                        int resID = residue.getResidueNumber();
                        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                        sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                        // Loop over atoms
                        /*ArrayList<Atom> residueAtoms = residue.getAtomList();
                        ArrayList<Atom> backboneAtoms = residue.getBackboneAtoms();*/
                        List<Atom> residueAtoms = residue.getAtomList().stream().
                                filter(a -> !atomExclusions.contains(a)).
                                collect(Collectors.toList());
                        List<Atom> backboneAtoms = residue.getBackboneAtoms().stream().
                                filter(a -> !atomExclusions.contains(a)).
                                collect(Collectors.toList());
                        boolean altLocFound = false;
                        for (Atom atom : backboneAtoms) {
                            writeAtom(atom, serial++, sb, anisouSB, bw);
                            Character altLoc = atom.getAltLoc();
                            if (altLoc != null && !altLoc.equals(' ')) {
                                altLocFound = true;
                            }
                            residueAtoms.remove(atom);
                        }
                        for (Atom atom : residueAtoms) {
                            writeAtom(atom, serial++, sb, anisouSB, bw);
                            Character altLoc = atom.getAltLoc();
                            if (altLoc != null && !altLoc.equals(' ')) {
                                altLocFound = true;
                            }
                        }
                        // Write out alternate conformers
                        if (altLocFound) {
                            for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                                MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                                Polymer altPolymer = altMolecularAssembly.getPolymer(currentChainID,
                                        currentSegID, false);
                                Residue altResidue = altPolymer.getResidue(resName, resID, false);
                                backboneAtoms = altResidue.getBackboneAtoms();
                                residueAtoms = altResidue.getAtomList();
                                for (Atom atom : backboneAtoms) {
                                    if (atom.getAltLoc() != null
                                            && !atom.getAltLoc().equals(' ')
                                            && !atom.getAltLoc().equals('A')) {
                                        writeAtom(atom, serial++, sb, anisouSB, bw);
                                    }
                                    residueAtoms.remove(atom);
                                }
                                for (Atom atom : residueAtoms) {
                                    if (atom.getAltLoc() != null
                                            && !atom.getAltLoc().equals(' ')
                                            && !atom.getAltLoc().equals('A')) {
                                        writeAtom(atom, serial++, sb, anisouSB, bw);
                                    }
                                }
                            }
                        }
                    }
                    terSB.replace(6, 11, String.format("%5s", Hybrid36.encode(5, serial++)));
                    terSB.replace(12, 16, "    ");
                    terSB.replace(16, 26, sb.substring(16, 26));
                    if (!listMode) {
                        bw.write(terSB.toString());
                        bw.newLine();
                    } else {
                        listOutput.add(terSB.toString());
                    }
                }
            }
            sb.replace(0, 6, "HETATM");
            sb.setCharAt(21, 'A');
            int resID = 1;
            Polymer polymer = activeMolecularAssembly.getPolymer('A', "A", false);
            if (polymer != null) {
                ArrayList<Residue> residues = polymer.getResidues();
                for (Residue residue : residues) {
                    int resID2 = residue.getResidueNumber();
                    if (resID2 >= resID) {
                        resID = resID2 + 1;
                    }
                }
            }

            /**
             * Loop over molecules, ions and then water.
             */
            ArrayList<MSNode> molecules = activeMolecularAssembly.getMolecules();
            for (int i = 0; i < molecules.size(); i++) {
                Molecule molecule = (Molecule) molecules.get(i);
                Character chainID = molecule.getChainID();
                sb.setCharAt(21, chainID);
                String resName = molecule.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                //ArrayList<Atom> moleculeAtoms = molecule.getAtomList();
                List<Atom> moleculeAtoms = molecule.getAtomList().stream().
                        filter(a -> !atomExclusions.contains(a)).
                        collect(Collectors.toList());
                boolean altLocFound = false;
                for (Atom atom : moleculeAtoms) {
                    writeAtom(atom, serial++, sb, anisouSB, bw);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altmolecule = altMolecularAssembly.getMolecules().get(i);
                        moleculeAtoms = altmolecule.getAtomList();
                        for (Atom atom : moleculeAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeAtom(atom, serial++, sb, anisouSB, bw);
                            }
                        }
                    }
                }
                resID++;
            }

            ArrayList<MSNode> ions = activeMolecularAssembly.getIons();
            for (int i = 0; i < ions.size(); i++) {
                Molecule ion = (Molecule) ions.get(i);
                Character chainID = ion.getChainID();
                sb.setCharAt(21, chainID);
                String resName = ion.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                //ArrayList<Atom> ionAtoms = ion.getAtomList();
                List<Atom> ionAtoms = ion.getAtomList().stream().
                        filter(a -> !atomExclusions.contains(a)).
                        collect(Collectors.toList());
                boolean altLocFound = false;
                for (Atom atom : ionAtoms) {
                    writeAtom(atom, serial++, sb, anisouSB, bw);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altion = altMolecularAssembly.getIons().get(i);
                        ionAtoms = altion.getAtomList();
                        for (Atom atom : ionAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeAtom(atom, serial++, sb, anisouSB, bw);
                            }
                        }
                    }
                }
                resID++;
            }

            ArrayList<MSNode> waters = activeMolecularAssembly.getWaters();
            for (int i = 0; i < waters.size(); i++) {
                Molecule water = (Molecule) waters.get(i);
                Character chainID = water.getChainID();
                sb.setCharAt(21, chainID);
                String resName = water.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                //ArrayList<Atom> waterAtoms = water.getAtomList();
                List<Atom> waterAtoms = water.getAtomList().stream().
                        filter(a -> !atomExclusions.contains(a)).
                        collect(Collectors.toList());
                boolean altLocFound = false;
                for (Atom atom : waterAtoms) {
                    writeAtom(atom, serial++, sb, anisouSB, bw);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altwater = altMolecularAssembly.getWaters().get(i);
                        waterAtoms = altwater.getAtomList();
                        for (Atom atom : waterAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeAtom(atom, serial++, sb, anisouSB, bw);
                            }
                        }
                    }
                }
                resID++;
            }

            if (writeEnd) {
                String end = model != null ? "ENDMDL" : "END";
                if (!listMode) {
                    bw.write(end);
                    bw.newLine();
                } else {
                    listOutput.add(end);
                }
            }
            bw.close();
        } catch (Exception e) {
            String message = "Exception writing to file: " + saveFile.toString();
            logger.log(Level.WARNING, message, e);
            return false;
        }
        return true;
    }

    /**
     * <p>writeSIFTFile.</p>
     *
     * @param saveFile    a {@link java.io.File} object.
     * @param append      a boolean.
     * @param resAndScore an array of {@link java.lang.String} objects.
     * @return a boolean.
     */
    public boolean writeSIFTFile(File saveFile, boolean append, String[] resAndScore) {
        if (saveFile == null) {
            return false;
        }

        if (vdwH) {
            logger.info(" Printing hydrogens to van der Waals centers instead of nuclear locations.");
        }

        if (nSymOp != 0) {
            logger.info(String.format(" Printing atoms with symmetry operator %s",
                    activeMolecularAssembly.getCrystal().spaceGroup.getSymOp(nSymOp).toString()));
        }

        /**
         * Create StringBuilders for ATOM, ANISOU and TER records that can be
         * reused.
         */
        StringBuilder sb = new StringBuilder("ATOM  ");
        StringBuilder anisouSB = new StringBuilder("ANISOU");
        StringBuilder terSB = new StringBuilder("TER   ");
        for (int i = 6; i < 80; i++) {
            sb.append(' ');
            anisouSB.append(' ');
            terSB.append(' ');
        }
        FileWriter fw;
        BufferedWriter bw;
        try {
            File newFile = saveFile;
            if (!append && !noVersioning) {
                newFile = version(saveFile);
            }
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (logWrites) {
                logger.log(Level.INFO, " Saving {0}", newFile.getName());
            }
            fw = new FileWriter(newFile, append);
            bw = new BufferedWriter(fw);
// =============================================================================
// The CRYST1 record presents the unit cell parameters, space group, and Z
// value. If the structure was not determined by crystallographic means, CRYST1
// simply provides the unitary values, with an appropriate REMARK.
//
//  7 - 15       Real(9.3)     a              a (Angstroms).
// 16 - 24       Real(9.3)     b              b (Angstroms).
// 25 - 33       Real(9.3)     c              c (Angstroms).
// 34 - 40       Real(7.2)     alpha          alpha (degrees).
// 41 - 47       Real(7.2)     beta           beta (degrees).
// 48 - 54       Real(7.2)     gamma          gamma (degrees).
// 56 - 66       LString       sGroup         Space  group.
// 67 - 70       Integer       z              Z value.
// =============================================================================
            Crystal crystal = activeMolecularAssembly.getCrystal();
            if (crystal != null && !crystal.aperiodic()) {
                Crystal c = crystal.getUnitCell();
                if (!listMode) {
                    bw.write(format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s\n", c.a, c.b, c.c, c.alpha, c.beta,
                            c.gamma, padRight(c.spaceGroup.pdbName, 10)));
                } else {
                    listOutput.add(format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s", c.a, c.b, c.c, c.alpha, c.beta,
                            c.gamma, padRight(c.spaceGroup.pdbName, 10)));
                }
            }
// =============================================================================
// The SSBOND record identifies each disulfide bond in protein and polypeptide
// structures by identifying the two residues involved in the bond.
// The disulfide bond distance is included after the symmetry operations at
// the end of the SSBOND record.
//
//  8 - 10        Integer         serNum       Serial number.
// 12 - 14        LString(3)      "CYS"        Residue name.
// 16             Character       chainID1     Chain identifier.
// 18 - 21        Integer         seqNum1      Residue sequence number.
// 22             AChar           icode1       Insertion code.
// 26 - 28        LString(3)      "CYS"        Residue name.
// 30             Character       chainID2     Chain identifier.
// 32 - 35        Integer         seqNum2      Residue sequence number.
// 36             AChar           icode2       Insertion code.
// 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
// 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
// 74 – 78        Real(5.2)      Length        Disulfide bond distance
//
// If SG of cysteine is disordered then there are possible alternate linkages.
// wwPDB practice is to put together all possible SSBOND records. This is
// problematic because the alternate location identifier is not specified in
// the SSBOND record.
// =============================================================================
            int serNum = 1;
            Polymer polymers[] = activeMolecularAssembly.getChains();
            if (polymers != null) {
                for (Polymer polymer : polymers) {
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        if (residue.getName().equalsIgnoreCase("CYS")) {
                            List<Atom> cysAtoms = residue.getAtomList();
                            Atom SG1 = null;
                            for (Atom atom : cysAtoms) {
                                if (atom.getName().equalsIgnoreCase("SG")) {
                                    SG1 = atom;
                                    break;
                                }
                            }
                            List<Bond> bonds = SG1.getBonds();
                            for (Bond bond : bonds) {
                                Atom SG2 = bond.get1_2(SG1);
                                if (SG2.getName().equalsIgnoreCase("SG")) {
                                    if (SG1.getIndex() < SG2.getIndex()) {
                                        bond.energy(false);
                                        if (!listMode) {
                                            bw.write(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n",
                                                    serNum++,
                                                    SG1.getChainID().toString(), Hybrid36.encode(4, SG1.getResidueNumber()),
                                                    SG2.getChainID().toString(), Hybrid36.encode(4, SG2.getResidueNumber()),
                                                    "", bond.getValue()));
                                        } else {
                                            listOutput.add(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n",
                                                    serNum++,
                                                    SG1.getChainID().toString(), Hybrid36.encode(4, SG1.getResidueNumber()),
                                                    SG2.getChainID().toString(), Hybrid36.encode(4, SG2.getResidueNumber()),
                                                    "", bond.getValue()));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
// =============================================================================
//
//  7 - 11        Integer       serial       Atom serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
// =============================================================================
//         1         2         3         4         5         6         7
//123456789012345678901234567890123456789012345678901234567890123456789012345678
//ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
//ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
            MolecularAssembly molecularAssemblies[] = this.getMolecularAssemblys();
            int serial = 1;
            // Loop over biomolecular chains
            if (polymers != null) {
                for (Polymer polymer : polymers) {
                    currentSegID = polymer.getName();
                    currentChainID = polymer.getChainID();
                    sb.setCharAt(21, currentChainID);
                    // Loop over residues
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        String resName = residue.getName();
                        if (resName.length() > 3) {
                            resName = resName.substring(0, 3);
                        }
                        int resID = residue.getResidueNumber();
                        int i = 0;
                        String[] entries = null;
                        for (; i < resAndScore.length; i++) {
                            entries = resAndScore[i].split("\\t");
                            if (!entries[0].equals(entries[0].replaceAll("\\D+", ""))) {
                                String[] subEntries = entries[0].split("[^0-9]");
                                entries[0] = subEntries[0];
                            }
                            if (entries[0].equals(String.valueOf(resID))
                                    && !".".equals(entries[1])) {
                                break;
                            }
                        }
                        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                        sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                        // Loop over atoms
                        ArrayList<Atom> residueAtoms = residue.getAtomList();
                        boolean altLocFound = false;
                        for (Atom atom : residueAtoms) {
                            if (i != resAndScore.length) {
                                writeSIFTAtom(atom, serial++, sb, anisouSB, bw, entries[1]);
                            } else {
                                writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                            }
                            Character altLoc = atom.getAltLoc();
                            if (altLoc != null && !altLoc.equals(' ')) {
                                altLocFound = true;
                            }
                        }
                        // Write out alternate conformers
                        if (altLocFound) {
                            for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                                MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                                Polymer altPolymer = altMolecularAssembly.getPolymer(currentChainID,
                                        currentSegID, false);
                                Residue altResidue = altPolymer.getResidue(resName, resID, false);
                                residueAtoms = altResidue.getAtomList();
                                for (Atom atom : residueAtoms) {
                                    if (atom.getAltLoc() != null
                                            && !atom.getAltLoc().equals(' ')
                                            && !atom.getAltLoc().equals('A')) {
                                        if (i != resAndScore.length) {
                                            writeSIFTAtom(atom, serial++, sb, anisouSB, bw, entries[1]);
                                        } else {
                                            writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    terSB.replace(6, 11, String.format("%5s", Hybrid36.encode(5, serial++)));
                    terSB.replace(12, 16, "    ");
                    terSB.replace(16, 26, sb.substring(16, 26));
                    if (!listMode) {
                        bw.write(terSB.toString());
                        bw.newLine();
                    } else {
                        listOutput.add(terSB.toString());
                    }
                }
            }
            sb.replace(0, 6, "HETATM");
            sb.setCharAt(21, 'A');
            int resID = 1;
            Polymer polymer = activeMolecularAssembly.getPolymer('A', "A", false);
            if (polymer != null) {
                ArrayList<Residue> residues = polymer.getResidues();
                for (Residue residue : residues) {
                    int resID2 = residue.getResidueNumber();
                    if (resID2 >= resID) {
                        resID = resID2 + 1;
                    }
                }
            }

            /**
             * Loop over molecules, ions and then water.
             */
            ArrayList<MSNode> molecules = activeMolecularAssembly.getMolecules();
            for (int i = 0; i < molecules.size(); i++) {
                Molecule molecule = (Molecule) molecules.get(i);
                Character chainID = molecule.getChainID();
                sb.setCharAt(21, chainID);
                String resName = molecule.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                ArrayList<Atom> moleculeAtoms = molecule.getAtomList();
                boolean altLocFound = false;
                for (Atom atom : moleculeAtoms) {
                    writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altmolecule = altMolecularAssembly.getMolecules().get(i);
                        moleculeAtoms = altmolecule.getAtomList();
                        for (Atom atom : moleculeAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                            }
                        }
                    }
                }
                resID++;
            }

            ArrayList<MSNode> ions = activeMolecularAssembly.getIons();
            for (int i = 0; i < ions.size(); i++) {
                Molecule ion = (Molecule) ions.get(i);
                Character chainID = ion.getChainID();
                sb.setCharAt(21, chainID);
                String resName = ion.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                ArrayList<Atom> ionAtoms = ion.getAtomList();
                boolean altLocFound = false;
                for (Atom atom : ionAtoms) {
                    writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altion = altMolecularAssembly.getIons().get(i);
                        ionAtoms = altion.getAtomList();
                        for (Atom atom : ionAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                            }
                        }
                    }
                }
                resID++;
            }

            ArrayList<MSNode> waters = activeMolecularAssembly.getWaters();
            for (int i = 0; i < waters.size(); i++) {
                Molecule water = (Molecule) waters.get(i);
                Character chainID = water.getChainID();
                sb.setCharAt(21, chainID);
                String resName = water.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                ArrayList<Atom> waterAtoms = water.getAtomList();
                boolean altLocFound = false;
                for (Atom atom : waterAtoms) {
                    writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                    Character altLoc = atom.getAltLoc();
                    if (altLoc != null && !altLoc.equals(' ')) {
                        altLocFound = true;
                    }
                }
                // Write out alternate conformers
                if (altLocFound) {
                    for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                        MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                        MSNode altwater = altMolecularAssembly.getWaters().get(i);
                        waterAtoms = altwater.getAtomList();
                        for (Atom atom : waterAtoms) {
                            if (atom.getAltLoc() != null
                                    && !atom.getAltLoc().equals(' ')
                                    && !atom.getAltLoc().equals('A')) {
                                writeSIFTAtom(atom, serial++, sb, anisouSB, bw, null);
                            }
                        }
                    }
                }
                resID++;
            }

            if (!listMode) {
                bw.write("END");
                bw.newLine();
            } else {
                listOutput.add("END");
            }
            bw.close();
        } catch (Exception e) {
            String message = "Exception writing to file: " + saveFile.toString();
            logger.log(Level.WARNING, message, e);
            return false;
        }
        return true;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Write out the Atomic information in PDB format.
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        return writeFile(saveFile, append, false, true);
    }

    /**
     * Locate disulfide bonds based on SSBOND records.
     *
     * @param ssbonds
     * @return An ArrayList of Bond instances for SS-Bonds.
     */
    private List<Bond> locateDisulfideBonds(List<String> ssbonds) {
        List<Bond> ssBondList = new ArrayList<>();
        for (String ssbond : ssbonds) {
            // =============================================================================
            // The SSBOND record identifies each disulfide bond in protein and polypeptide
            // structures by identifying the two residues involved in the bond.
            // The disulfide bond distance is included after the symmetry operations at
            // the end of the SSBOND record.
            //
            //  8 - 10        Integer         serNum       Serial number.
            // 12 - 14        LString(3)      "CYS"        Residue name.
            // 16             Character       chainID1     Chain identifier.
            // 18 - 21        Integer         seqNum1      Residue sequence number.
            // 22             AChar           icode1       Insertion code.
            // 26 - 28        LString(3)      "CYS"        Residue name.
            // 30             Character       chainID2     Chain identifier.
            // 32 - 35        Integer         seqNum2      Residue sequence number.
            // 36             AChar           icode2       Insertion code.
            // 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
            // 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
            // 74 – 78        Real(5.2)      Length        Disulfide bond distance
            //
            // If SG of cysteine is disordered then there are possible alternate linkages.
            // wwPDB practice is to put together all possible SSBOND records. This is
            // problematic because the alternate location identifier is not specified in
            // the SSBOND record.
            // =============================================================================
            try {
                char c1ch = ssbond.charAt(15);
                char c2ch = ssbond.charAt(29);
                Polymer c1 = activeMolecularAssembly.getChain(String.format("%c", c1ch));
                Polymer c2 = activeMolecularAssembly.getChain(String.format("%c", c2ch));

                Polymer[] chains = activeMolecularAssembly.getChains();
                if (c1 == null) {
                    c1 = chains[0];
                }
                if (c2 == null) {
                    c2 = chains[0];
                }

                String origResNum1 = ssbond.substring(17, 21).trim();
                char insChar1 = ssbond.charAt(21);
                String origResNum2 = ssbond.substring(31, 35).trim();
                char insChar2 = ssbond.charAt(35);

                String pdbResNum1 = String.format("%c%s%c", c1ch, origResNum1, insChar1);
                String pdbResNum2 = String.format("%c%s%c", c2ch, origResNum2, insChar2);
                String resnum1 = pdbToNewResMap.get(pdbResNum1);
                String resnum2 = pdbToNewResMap.get(pdbResNum2);
                if (resnum1 == null) {
                    logger.warning(String.format(" Could not find residue %s for SS-bond %s", pdbResNum1, ssbond));
                    continue;
                }
                if (resnum2 == null) {
                    logger.warning(String.format(" Could not find residue %s for SS-bond %s", pdbResNum2, ssbond));
                    continue;
                }

                Residue r1 = c1.getResidue(Integer.parseInt(resnum1.substring(1)));
                Residue r2 = c2.getResidue(Integer.parseInt(resnum2.substring(1)));
                /*Residue r1 = c1.getResidue(Hybrid36.decode(4, ssbond.substring(17, 21)));
                Residue r2 = c2.getResidue(Hybrid36.decode(4, ssbond.substring(31, 35)));*/
                List<Atom> atoms1 = r1.getAtomList();
                List<Atom> atoms2 = r2.getAtomList();
                Atom SG1 = null;
                Atom SG2 = null;
                for (Atom atom : atoms1) {
                    if (atom.getName().equalsIgnoreCase("SG")) {
                        SG1 = atom;
                        break;
                    }
                }
                for (Atom atom : atoms2) {
                    if (atom.getName().equalsIgnoreCase("SG")) {
                        SG2 = atom;
                        break;
                    }
                }
                if (SG1 == null) {
                    logger.warning(String.format(" SG atom 1 of SS-bond %s is null", ssbond));
                }
                if (SG2 == null) {
                    logger.warning(String.format(" SG atom 2 of SS-bond %s is null", ssbond));
                }
                if (SG1 == null || SG2 == null) {
                    continue;
                }
                double d = VectorMath.dist(SG1.getXYZ(null), SG2.getXYZ(null));
                if (d < 5.0) {
                    r1.setName("CYX");
                    r2.setName("CYX");
                    for (Atom atom : atoms1) {
                        atom.setResName("CYX");
                    }
                    for (Atom atom : atoms2) {
                        atom.setResName("CYX");
                    }
                    Bond bond = new Bond(SG1, SG2);
                    ssBondList.add(bond);
                } else {
                    String message = format("Ignoring [%s]\n due to distance %8.3f A.", ssbond, d);
                    logger.log(Level.WARNING, message);
                }
            } catch (Exception e) {
                String message = format("Ignoring [%s]", ssbond);
                logger.log(Level.WARNING, message, e);
            }
        }
        return ssBondList;
    }

    /**
     * Assign parameters to disulfide bonds.
     *
     * @param ssBondList
     */
    private void buildDisulfideBonds(List<Bond> ssBondList) {
        StringBuilder sb = new StringBuilder(" Disulfide Bonds:");
        for (Bond bond : ssBondList) {
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            int c[] = new int[2];
            c[0] = a1.getAtomType().atomClass;
            c[1] = a2.getAtomType().atomClass;
            String key = BondType.sortKey(c);
            BondType bondType = forceField.getBondType(key);
            if (bondType == null) {
                logger.severe(format("No BondType for key: %s\n %s\n %s", key, a1, a2));
            } else {
                bond.setBondType(bondType);
            }
            double d = VectorMath.dist(a1.getXYZ(null), a2.getXYZ(null));
            Polymer c1 = activeMolecularAssembly.getChain(a1.getSegID());
            Polymer c2 = activeMolecularAssembly.getChain(a2.getSegID());
            Residue r1 = c1.getResidue(a1.getResidueNumber());
            Residue r2 = c2.getResidue(a2.getResidueNumber());
            sb.append(format("\n S-S distance of %6.2f for %s and %s.", d, r1.toString(), r2.toString()));
            bondList.add(bond);
        }
        if (ssBondList.size() > 0) {
            logger.info(sb.toString());
        }
    }

    /**
     * Currently builds missing internal loops based on information in DBREF and
     * SEQRES records.
     * <p>
     * Known limitations include: 1) No building n- and c-terminal loops. 2) No
     * support for DBREF1 or DBREF2 records. 3) Incomplete optimization scheme
     * to position the loops.
     *
     * @param xyzIndex
     * @return xyzIndex updated based on built atoms.
     */
    private int buildMissingResidues(int xyzIndex) {

        /**
         * Only build loops if the buildLoops flag is true.
         */
        if (!properties.getBoolean("buildLoops", false)) {
            return xyzIndex;
        }
        Polymer polymers[] = activeMolecularAssembly.getChains();
        for (Polymer polymer : polymers) {

            Character chainID = polymer.getChainID();
            String resNames[] = seqres.get(chainID);
            int seqRange[] = dbref.get(chainID);
            if (resNames == null || seqRange == null) {
                continue;
            }
            int seqBegin = seqRange[0];
            int seqEnd = seqRange[1];
            logger.info(format("\n Checking for missing residues in chain %s between residues %d and %d.",
                    polymer.toString(), seqBegin, seqEnd));

            int firstResID = polymer.getFirstResidue().getResidueNumber();
            for (int i = 0; i < resNames.length; i++) {
                int currentID = seqBegin + i;

                Residue currentResidue = polymer.getResidue(currentID);
                if (currentResidue != null) {
                    continue;
                }

                if (currentID <= firstResID) {
                    logger.info(format(" Residue %d is missing, but is at the beginning of the chain.",
                            currentID));
                    continue;
                }

                Residue previousResidue = polymer.getResidue(currentID - 1);
                if (previousResidue == null) {
                    logger.info(format(" Residue %d is missing, but could not be build (previous residue missing).",
                            currentID));
                    continue;
                }

                Residue nextResidue = null;
                for (int j = currentID + 1; j <= seqEnd; j++) {
                    nextResidue = polymer.getResidue(j);
                    if (nextResidue != null) {
                        break;
                    }
                }
                /**
                 * Residues at the end of the chain are not currently built.
                 */
                if (nextResidue == null) {
                    logger.info(format(" Residue %d is missing, but is at the end of the chain.",
                            currentID));
                    break;
                }
                /**
                 * Find the previous carbonyl carbon and next nitrogen.
                 */
                Atom C = (Atom) previousResidue.getAtomNode("C");
                Atom N = (Atom) nextResidue.getAtomNode("N");
                if (C == null || N == null) {
                    logger.info(format(" Residue %d is missing, but bonding atoms are missing (C or N).",
                            currentID));
                    continue;
                }

                /**
                 * Build the missing residue.
                 */
                currentResidue = polymer.getResidue(resNames[i], currentID, true);

                double vector[] = new double[3];
                int count = 3 * (nextResidue.getResidueNumber() - previousResidue.getResidueNumber());
                VectorMath.diff(N.getXYZ(null), C.getXYZ(null), vector);
                VectorMath.scalar(vector, 1.0 / count, vector);

                double nXYZ[] = new double[3];
                VectorMath.sum(C.getXYZ(null), vector, nXYZ);
                nXYZ[0] += Math.random() - 0.5;
                nXYZ[1] += Math.random() - 0.5;
                nXYZ[2] += Math.random() - 0.5;
                Atom newN = new Atom(xyzIndex++, "N", C.getAltLoc(), nXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newN);

                double caXYZ[] = new double[3];
                VectorMath.scalar(vector, 2.0, vector);
                VectorMath.sum(C.getXYZ(null), vector, caXYZ);
                caXYZ[0] += Math.random() - 0.5;
                caXYZ[1] += Math.random() - 0.5;
                caXYZ[2] += Math.random() - 0.5;
                Atom newCA = new Atom(xyzIndex++, "CA", C.getAltLoc(), caXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newCA);

                double cXYZ[] = new double[3];
                VectorMath.scalar(vector, 1.5, vector);
                VectorMath.sum(C.getXYZ(null), vector, cXYZ);
                cXYZ[0] += Math.random() - 0.5;
                cXYZ[1] += Math.random() - 0.5;
                cXYZ[2] += Math.random() - 0.5;
                Atom newC = new Atom(xyzIndex++, "C", C.getAltLoc(), cXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newC);

                double oXYZ[] = new double[3];
                vector[0] = Math.random() - 0.5;
                vector[1] = Math.random() - 0.5;
                vector[2] = Math.random() - 0.5;
                VectorMath.sum(cXYZ, vector, oXYZ);
                Atom newO = new Atom(xyzIndex++, "O", C.getAltLoc(), oXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newO);
                logger.info(String.format(" Building residue %8s.", currentResidue.toString()));
            }
        }
        return xyzIndex;
    }

    /**
     * <p>
     * numberAtoms</p>
     */
    private void numberAtoms() {
        int index = 1;
        for (Atom a : activeMolecularAssembly.getAtomArray()) {
            a.setXyzIndex(index++);
        }
        index--;
        if (logger.isLoggable(Level.INFO)) {
            logger.info(String.format(" Total number of atoms: %d\n", index));
        }

        Polymer[] polymers = activeMolecularAssembly.getChains();
        if (polymers != null) {
            for (Polymer p : polymers) {
                List<Residue> residues = p.getResidues();
                for (Residue r : residues) {
                    r.reOrderAtoms();
                }
            }
        }
        List<MSNode> molecules = activeMolecularAssembly.getMolecules();
        for (MSNode n : molecules) {
            MSGroup m = (MSGroup) n;
            m.reOrderAtoms();
        }
        List<MSNode> waters = activeMolecularAssembly.getWaters();
        for (MSNode n : waters) {
            MSGroup m = (MSGroup) n;
            m.reOrderAtoms();
        }
        List<MSNode> ions = activeMolecularAssembly.getIons();
        for (MSNode n : ions) {
            MSGroup m = (MSGroup) n;
            m.reOrderAtoms();
        }

    }

    /**
     * Assign force field atoms types to common chemistries using "biotype"
     * records.
     */
    private void assignAtomTypes() {
        /**
         * Create a list to store bonds defined by PDB atom names.
         */
        bondList = new ArrayList<>();

        /**
         * To Do: Look for cyclic peptides and disulfides.
         */
        Polymer[] polymers = activeMolecularAssembly.getChains();

        /**
         * Loop over chains.
         */
        if (polymers != null) {
            logger.info(format("\n Assigning atom types for %d chains.", polymers.length));
            for (Polymer polymer : polymers) {
                List<Residue> residues = polymer.getResidues();
                int numberOfResidues = residues.size();
                /**
                 * Check if all residues are known amino acids.
                 */
                boolean isProtein = true;
                for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
                    Residue residue = residues.get(residueNumber);
                    String name = residue.getName().toUpperCase();
                    boolean aa = false;
                    for (AminoAcid3 amino : aminoAcidList) {
                        if (amino.toString().equalsIgnoreCase(name)) {
                            aa = true;
                            checkHydrogenAtomNames(residue, fileStandard);
                            break;
                        }
                    }
                    // Check for a patch.
                    if (!aa) {
                        logger.info(" Checking for non-standard amino acid patch " + name);
                        HashMap<String, AtomType> types = forceField.getAtomTypes(name);
                        if (types.isEmpty()) {
                            isProtein = false;
                            break;
                        } else {
                            logger.info(" Patch found for non-standard amino acid " + name);
                        }
                    }
                }

                /**
                 * If all the residues in this chain have known amino acids
                 * names, then attempt to assign atom types.
                 */
                if (isProtein) {
                    try {
                        logger.info(format(" Amino acid chain %s", polymer.getName()));
                        double dist = properties.getDouble("chainbreak", 3.0);
                        // Detect main chain breaks!
                        List<List<Residue>> subChains = findChainBreaks(residues, dist);
                        for (List<Residue> subChain : subChains) {
                            assignAminoAcidAtomTypes(subChain);
                        }
                    } catch (MissingHeavyAtomException missingHeavyAtomException) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(missingHeavyAtomException));
                        logger.severe(missingHeavyAtomException.toString());
                    } catch (MissingAtomTypeException missingAtomTypeException) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(missingAtomTypeException));
                        logger.severe(missingAtomTypeException.toString());
                    }
                    continue;
                }

                /**
                 * Check if all residues have known nucleic acids names.
                 */
                boolean isNucleicAcid = true;
                for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
                    Residue residue = residues.get(residueNumber);
                    String name = residue.getName().toUpperCase();
                    /**
                     * Convert 1 and 2-character nucleic acid names to
                     * 3-character names.
                     */
                    if (name.length() == 1) {
                        if (name.equals("A")) {
                            name = NucleicAcid3.ADE.toString();
                        } else if (name.equals("C")) {
                            name = NucleicAcid3.CYT.toString();
                        } else if (name.equals("G")) {
                            name = NucleicAcid3.GUA.toString();
                        } else if (name.equals("T")) {
                            name = NucleicAcid3.THY.toString();
                        } else if (name.equals("U")) {
                            name = NucleicAcid3.URI.toString();
                        }
                    } else if (name.length() == 2) {
                        if (name.equals("YG")) {
                            name = NucleicAcid3.YYG.toString();
                        }
                    }
                    residue.setName(name);
                    NucleicAcid3 nucleicAcid = null;
                    for (NucleicAcid3 nucleic : nucleicAcidList) {
                        String nuc3 = nucleic.toString();
                        nuc3 = nuc3.substring(nuc3.length() - 3);
                        if (nuc3.equalsIgnoreCase(name)) {
                            nucleicAcid = nucleic;
                            break;
                        }
                    }
                    if (nucleicAcid == null) {
                        logger.info(format("Nucleic acid was not recognized %s.", name));
                        isNucleicAcid = false;
                        break;
                    }
                }

                /**
                 * If all the residues in this chain have known nucleic acids
                 * names, then attempt to assign atom types.
                 */
                if (isNucleicAcid) {
                    try {
                        logger.info(format(" Nucleic acid chain %s", polymer.getName()));
                        logger.info(String.format(" EXPERIMENTAL: Finding chain breaks for possible nucleic acid chain %s", polymer.getName()));
                        double dist = properties.getDouble("chainbreak", 4.0);
                        // Detect main chain breaks!
                        List<List<Residue>> subChains = findChainBreaks(residues, dist);

                        for (List<Residue> subChain : subChains) {
                            assignNucleicAcidAtomTypes(subChain, forceField, bondList);
                        }
                    } catch (MissingHeavyAtomException | MissingAtomTypeException e) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(e));
                        logger.severe(e.toString());
                    }
                }
            }
        }

        // Assign ion atom types.
        ArrayList<MSNode> ions = activeMolecularAssembly.getIons();
        if (ions != null && ions.size() > 0) {
            logger.info(format(" Assigning atom types for %d ions.", ions.size()));
            for (MSNode m : ions) {
                Molecule ion = (Molecule) m;
                String name = ion.getResidueName().toUpperCase();
                HetAtoms hetatm = HetAtoms.valueOf(name);
                Atom atom = ion.getAtomList().get(0);
                if (ion.getAtomList().size() != 1) {
                    logger.severe(format(" Check residue %s of chain %s.", ion.toString(), ion.getChainID()));
                }
                try {
                    switch (hetatm) {
                        case NA:
                            atom.setAtomType(findAtomType(2003));
                            break;
                        case K:
                            atom.setAtomType(findAtomType(2004));
                            break;
                        case MG:
                        case MG2:
                            atom.setAtomType(findAtomType(2005));
                            break;
                        case CA:
                        case CA2:
                            atom.setAtomType(findAtomType(2006));
                            break;
                        case CL:
                            atom.setAtomType(findAtomType(2007));
                            break;
                        case ZN:
                        case ZN2:
                            atom.setAtomType(findAtomType(2008));
                            break;
                        case BR:
                            atom.setAtomType(findAtomType(2009));
                            break;
                        default:
                            logger.severe(format(" Check residue %s of chain %s.", ion.toString(), ion.getChainID()));
                    }
                } catch (Exception e) {
                    logger.log(Level.INFO, Utilities.stackTraceToString(e));
                    String message = "Error assigning atom types.";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }
        // Assign water atom types.
        ArrayList<MSNode> water = activeMolecularAssembly.getWaters();
        if (water != null && water.size() > 0) {
            logger.info(format(" Assigning atom types for %d waters.", water.size()));
            for (MSNode m : water) {
                Molecule wat = (Molecule) m;
                try {
                    Atom O = buildHeavy(wat, "O", null, 2001);
                    Atom H1 = buildHydrogen(wat, "H1", O, 0.96e0, null, 109.5e0, null, 120.0e0, 0, 2002);
                    H1.setHetero(true);
                    Atom H2 = buildHydrogen(wat, "H2", O, 0.96e0, H1, 109.5e0, null, 120.0e0, 0, 2002);
                    H2.setHetero(true);
                } catch (Exception e) {
                    logger.log(Level.INFO, Utilities.stackTraceToString(e));
                    String message = "Error assigning atom types to a water.";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }

        // Assign small molecule atom types.
        ArrayList<MSNode> molecules = activeMolecularAssembly.getMolecules();
        for (MSNode m : molecules) {
            Molecule molecule = (Molecule) m;
            String moleculeName = molecule.getResidueName();
            logger.info(" Attempting to patch " + moleculeName);
            ArrayList<Atom> moleculeAtoms = molecule.getAtomList();
            boolean patched = true;
            HashMap<String, AtomType> types = forceField.getAtomTypes(moleculeName);
            /**
             * Assign atom types for all known atoms.
             */
            for (Atom atom : moleculeAtoms) {
                String atomName = atom.getName().toUpperCase();
                AtomType atomType = types.get(atomName);
                if (atomType == null) {
                    logger.info(" No atom type was found for " + atomName + " of " + moleculeName + ".");
                    patched = false;
                    break;
                } else {
                    logger.fine(" " + atom.toString() + " -> " + atomType.toString());
                    atom.setAtomType(atomType);
                    types.remove(atomName);
                }
            }
            /**
             * Create missing hydrogen atoms. Check for missing heavy atoms.
             */
            if (patched && !types.isEmpty()) {
                for (AtomType type : types.values()) {
                    if (type.atomicNumber != 1) {
                        logger.info(" Missing heavy atom " + type.name);
                        patched = false;
                        break;
                    }
                }
            }
            // Create bonds between known atoms.
            if (patched) {
                for (Atom atom : moleculeAtoms) {
                    String atomName = atom.getName();
                    String bonds[] = forceField.getBonds(moleculeName, atomName);
                    if (bonds != null) {
                        for (String name : bonds) {
                            Atom atom2 = molecule.getAtom(name);
                            if (atom2 != null && !atom.isBonded(atom2)) {
                                buildBond(atom, atom2);
                            }
                        }
                    }
                }
            }
            // Create missing hydrogen atoms.
            if (patched && !types.isEmpty()) {
                // Create a hashmap of the molecule's atoms
                HashMap<String, Atom> atomMap = new HashMap<String, Atom>();
                for (Atom atom : moleculeAtoms) {
                    atomMap.put(atom.getName().toUpperCase(), atom);
                }
                for (String atomName : types.keySet()) {
                    AtomType type = types.get(atomName);
                    String bonds[] = forceField.getBonds(moleculeName, atomName.toUpperCase());
                    if (bonds == null || bonds.length != 1) {
                        patched = false;
                        logger.info(" Check biotype for hydrogen " + type.name + ".");
                        break;
                    }
                    // Get the heavy atom the hydrogen is bonded to.
                    Atom ia = atomMap.get(bonds[0].toUpperCase());
                    Atom hydrogen = new Atom(0, atomName, ia.getAltLoc(), new double[3],
                            ia.getResidueName(), ia.getResidueNumber(), ia.getChainID(),
                            ia.getOccupancy(), ia.getTempFactor(), ia.getSegID());
                    logger.fine(" Created hydrogen " + atomName + ".");
                    hydrogen.setAtomType(type);
                    hydrogen.setHetero(true);
                    molecule.addMSNode(hydrogen);
                    int valence = ia.getAtomType().valence;
                    List<Bond> aBonds = ia.getBonds();
                    int numBonds = aBonds.size();
                    /**
                     * Try to find the following configuration: ib-ia-ic
                     */
                    Atom ib = null;
                    Atom ic = null;
                    Atom id = null;
                    if (numBonds > 0) {
                        Bond bond = aBonds.get(0);
                        ib = bond.get1_2(ia);
                    }
                    if (numBonds > 1) {
                        Bond bond = aBonds.get(1);
                        ic = bond.get1_2(ia);
                    }
                    if (numBonds > 2) {
                        Bond bond = aBonds.get(2);
                        id = bond.get1_2(ia);
                    }

                    /**
                     * Building the hydrogens depends on hybridization and the
                     * locations of other bonded atoms.
                     */
                    logger.fine(" Bonding " + atomName + " to " + ia.getName()
                            + " (" + numBonds + " of " + valence + ").");
                    switch (valence) {
                        case 4:
                            switch (numBonds) {
                                case 3:
                                    // Find the average coordinates of atoms ib, ic and id.
                                    double b[] = ib.getXYZ(null);
                                    double c[] = ib.getXYZ(null);
                                    double d[] = ib.getXYZ(null);
                                    double a[] = new double[3];
                                    a[0] = (b[0] + c[0] + d[0]) / 3.0;
                                    a[1] = (b[1] + c[1] + d[1]) / 3.0;
                                    a[2] = (b[2] + c[2] + d[2]) / 3.0;
                                    // Place the hydrogen at chiral position #1.
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 1);
                                    double e1[] = new double[3];
                                    hydrogen.getXYZ(e1);
                                    double ret[] = new double[3];
                                    diff(a, e1, ret);
                                    double l1 = r(ret);
                                    // Place the hydrogen at chiral position #2.
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, -1);
                                    double e2[] = new double[3];
                                    hydrogen.getXYZ(e2);
                                    diff(a, e2, ret);
                                    double l2 = r(ret);
                                    // Revert to #1 if it is farther from the average.
                                    if (l1 > l2) {
                                        hydrogen.setXYZ(e1);
                                    }
                                    break;
                                case 2:
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 0);
                                    break;
                                case 1:
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, null, 0.0, 0);
                                    break;
                                case 0:
                                    intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                    break;
                                default:
                                    logger.info(" Check biotype for hydrogen " + atomName + ".");
                                    patched = false;
                            }
                            break;
                        case 3:
                            switch (numBonds) {
                                case 2:
                                    intxyz(hydrogen, ia, 1.0, ib, 120.0, ic, 0.0, 0);
                                    break;
                                case 1:
                                    intxyz(hydrogen, ia, 1.0, ib, 120.0, null, 0.0, 0);
                                    break;
                                case 0:
                                    intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                    break;
                                default:
                                    logger.info(" Check biotype for hydrogen " + atomName + ".");
                                    patched = false;
                            }
                            break;
                        case 2:
                            switch (numBonds) {
                                case 1:
                                    intxyz(hydrogen, ia, 1.0, ib, 120.0, null, 0.0, 0);
                                    break;
                                case 0:
                                    intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                    break;
                                default:
                                    logger.info(" Check biotype for hydrogen " + atomName + ".");
                                    patched = false;
                            }
                            break;
                        case 1:
                            switch (numBonds) {
                                case 0:
                                    intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                    break;
                                default:
                                    logger.info(" Check biotype for hydrogen " + atomName + ".");
                                    patched = false;
                            }
                            break;
                        default:
                            logger.info(" Check biotype for hydrogen " + atomName + ".");
                            patched = false;
                    }
                    if (!patched) {
                        break;
                    } else {
                        buildBond(ia, hydrogen);
                    }
                }
            }
            if (!patched) {
                logger.log(Level.WARNING, format(" Deleting unrecognized molecule %s.", m.toString()));
                activeMolecularAssembly.deleteMolecule((Molecule) m);
            } else {
                logger.info(" Patch for " + moleculeName + " succeeded.");
            }
        }
        resolvePolymerLinks(molecules);
    }

    /**
     * Resolves links between polymeric hetero groups; presently only functional
     * for cyclic molecules.
     *
     * @param molecules List of Molecules in the molecular assembly.
     */
    private void resolvePolymerLinks(List<MSNode> molecules) {
        CompositeConfiguration ffProps = forceField.getProperties();
        for (String polyLink : ffProps.getStringArray("polymerlink")) {
            logger.info(" Experimental: linking a cyclic hetero group: " + polyLink);

            // Format: polymerlink resname atom1 atom2 [cyclize]
            String[] toks = polyLink.split("\\s+");
            String resName = toks[0];
            String name1 = toks[1];
            String name2 = toks[2];
            int cyclicLen = 0;
            if (toks.length > 3) {
                cyclicLen = Integer.parseInt(toks[3]);
            }

            ArrayList<Molecule> matches = new ArrayList<>();
            for (MSNode node : molecules) {
                Molecule m = (Molecule) node;
                if (m.getResidueName().equalsIgnoreCase(resName)) {
                    matches.add(m);
                }
            }

//            Molecule[] matches = molecules.stream().
//                    filter((Molecule m) -> m.getResidueName().equalsIgnoreCase(resName)).
//                    toArray(Molecule[]::new);

            for (int i = 0; i < matches.size(); i++) {
                Molecule mi = matches.get(i);
                int ii = i + 1;
                if (cyclicLen < 1) {
                    logger.severe(" No current support for polymeric, non-cyclic hetero groups");
                    // Would probably split by chain.
                } else {
                    Molecule next;
                    if (ii % cyclicLen == 0) {
                        next = matches.get(ii - cyclicLen);
                        logger.info(String.format(" Cyclizing molecule %s to %s", mi, next));
                    } else {
                        next = matches.get(ii);
                        logger.info(String.format(" Extending chain from %s to %s.", mi, next));
                    }
                    Atom from = mi.getAtomByName(name1, true);
                    Atom to = next.getAtomByName(name2, true);
                    buildBond(from, to);
                }
            }
        }
    }

    private List<List<Residue>> findChainBreaks(List<Residue> residues, double cutoff) {
        List<List<Residue>> subChains = new ArrayList<List<Residue>>();

        /**
         * Chain-start atom: N (amino)/O5' (nucleic)
         * Chain-end atom:   C (amino)/O3' (nucleic)
         */
        Residue.ResidueType rType = residues.get(0).getResidueType();
        String startAtName;
        String endAtName;
        switch (rType) {
            case AA:
                startAtName = "N";
                endAtName = "C";
                break;
            case NA:
                startAtName = "O5\'";
                endAtName = "O3\'";
                break;
            case UNK:
            default:
                logger.fine(" Not attempting to find chain breaks for chain with residue " + residues.get(0).toString());
                List<List<Residue>> retList = new ArrayList<>();
                retList.add(residues);
                return retList;
        }

        List<Residue> subChain = null;
        Residue previousResidue = null;
        Atom priorEndAtom = null;
        StringBuilder sb = new StringBuilder(" Chain Breaks:");

        for (Residue residue : residues) {
            List<Atom> resAtoms = residue.getAtomList();
            if (priorEndAtom == null) {
                /**
                 * Initialization.
                 */
                subChain = new ArrayList<Residue>();
                subChain.add(residue);
                subChains.add(subChain);
            } else {
                /**
                 * Find the start atom of the current residue.
                 */
                Atom startAtom = null;
                for (Atom a : resAtoms) {
                    if (a.getName().equalsIgnoreCase(startAtName)) {
                        startAtom = a;
                        break;
                    }
                }
                if (startAtom == null) {
                    subChain.add(residue);
                    continue;
                }
                /**
                 * Compute the distance between the previous carbonyl carbon and
                 * the current nitrogen.
                 */
                double r = VectorMath.dist(priorEndAtom.getXYZ(null), startAtom.getXYZ(null));
                if (r > cutoff) {
                    /**
                     * Start a new chain.
                     */
                    subChain = new ArrayList<Residue>();
                    subChain.add(residue);
                    subChains.add(subChain);
                    char ch1 = previousResidue.getChainID();
                    char ch2 = residue.getChainID();
                    sb.append(format("\n C-N distance of %6.2f A for %c-%s and %c-%s.",
                            r, ch1, previousResidue.toString(), ch2, residue.toString()));
                } else {
                    /**
                     * Continue the current chain.
                     */
                    subChain.add(residue);
                }
            }

            /**
             * Save the carbonyl carbon.
             */
            for (Atom a : resAtoms) {
                if (a.getName().equalsIgnoreCase(endAtName)) {
                    priorEndAtom = a;
                    break;
                }
            }
            previousResidue = residue;
        }

        if (subChains.size() > 1) {
            logger.info(sb.toString());
        }

        return subChains;
    }

    /**
     * Assign atom types to an amino acid polymer.
     *
     * @param residues The residues to assign atom types to.
     * @throws MissingHeavyAtomException
     * @since 1.0
     */
    private void assignAminoAcidAtomTypes(List<Residue> residues)
            throws MissingHeavyAtomException, MissingAtomTypeException {
        /**
         * Loop over amino acid residues.
         */
        int numberOfResidues = residues.size();
        for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
            Residue residue = residues.get(residueNumber);
            Residue previousResidue = null;
            Residue nextResidue = null;
            if (residueNumber > 0) {
                previousResidue = residues.get(residueNumber - 1);
            }
            if (residueNumber < numberOfResidues - 1) {
                nextResidue = residues.get(residueNumber + 1);
            }
            AminoAcidUtils.assignAminoAcidAtomTypes(residue, previousResidue, nextResidue, forceField, bondList);
        }
    }

    /**
     * <p>
     * writeAtom</p>
     *
     * @param atom     a {@link ffx.potential.bonded.Atom} object.
     * @param serial   a int.
     * @param sb       a {@link java.lang.StringBuilder} object.
     * @param anisouSB a {@link java.lang.StringBuilder} object.
     * @param bw       a {@link java.io.BufferedWriter} object.
     * @throws java.io.IOException if any.
     */
    private void writeAtom(Atom atom, int serial, StringBuilder sb,
                           StringBuilder anisouSB, BufferedWriter bw)
            throws IOException {
        if (ignoreUnusedAtoms && !atom.getUse()) {
            return;
        }
        String name = atom.getName();
        if (name.length() > 4) {
            name = name.substring(0, 4);
        } else if (name.length() == 1) {
            name = name + "  ";
        } else if (name.length() == 2) {
            if (atom.getAtomType().valence == 0) {
                name = name + "  ";
            } else {
                name = name + " ";
            }
        }
        double xyz[] = vdwH ? atom.getRedXYZ() : atom.getXYZ(null);
        if (nSymOp != 0) {
            Crystal crystal = activeMolecularAssembly.getCrystal();
            SymOp symOp = crystal.spaceGroup.getSymOp(nSymOp);
            double[] newXYZ = new double[xyz.length];
            crystal.applySymOp(xyz, newXYZ, symOp);
            xyz = newXYZ;
        }
        sb.replace(6, 16, String.format("%5s " + padLeft(name.toUpperCase(), 4), Hybrid36.encode(5, serial)));
        Character altLoc = atom.getAltLoc();
        if (altLoc != null) {
            sb.setCharAt(16, altLoc);
        } else {
            sb.setCharAt(16, ' ');
        }

        /*sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                xyz[0], xyz[1], xyz[2], atom.getOccupancy(), atom.getTempFactor()));*/
        /**
         * On the following code: #1: StringBuilder.replace will allow for
         * longer strings, expanding the StringBuilder's length if necessary.
         * #2: sb was never re-initialized, so if there was overflow, sb would
         * continue to be > 80 characters long, resulting in broken PDB files
         * #3: It may be wiser to have XYZ coordinates result in shutdown, not
         * truncation of coordinates. #4: Excessive B-factors aren't much of an
         * issue; if the B-factor is past 999.99, that's the difference between
         * "density extends to Venus" and "density extends to Pluto".
         */
        StringBuilder decimals = new StringBuilder();
        for (int i = 0; i < 3; i++) {
            try {
                decimals.append(StringUtils.fwFpDec(xyz[i], 8, 3));
            } catch (IllegalArgumentException ex) {
                String newValue = StringUtils.fwFpTrunc(xyz[i], 8, 3);
                logger.info(String.format(" XYZ %d coordinate %8.3f for atom %s "
                                + "overflowed bounds of 8.3f string specified by PDB "
                                + "format; truncating value to %s", i, xyz[i], atom.toString(),
                        newValue));
                decimals.append(newValue);
            }
        }
        try {
            decimals.append(StringUtils.fwFpDec(atom.getOccupancy(), 6, 2));
        } catch (IllegalArgumentException ex) {
            logger.severe(String.format(" Occupancy %f for atom %s is impossible; "
                    + "value must be between 0 and 1", atom.getOccupancy(), atom.toString()));
        }
        try {
            decimals.append(StringUtils.fwFpDec(atom.getTempFactor(), 6, 2));
        } catch (IllegalArgumentException ex) {
            String newValue = StringUtils.fwFpTrunc(atom.getTempFactor(), 6, 2);
            logger.info(String.format(" Atom temp factor %6.2f for atom %s overflowed "
                    + "bounds of 6.2f string specified by PDB format; truncating "
                    + "value to %s", atom.getTempFactor(), atom.toString(), newValue));
            decimals.append(newValue);
        }
        sb.replace(30, 66, decimals.toString());

        name = Atom.ElementSymbol.values()[atom.getAtomicNumber() - 1].toString();
        name = name.toUpperCase();
        if (atom.isDeuterium()) {
            name = "D";
        }
        sb.replace(76, 78, padLeft(name, 2));
        sb.replace(78, 80, String.format("%2d", 0));
        if (!listMode) {
            bw.write(sb.toString());
            bw.newLine();
        } else {
            listOutput.add(sb.toString());
        }
// =============================================================================
//  1 - 6        Record name   "ANISOU"
//  7 - 11       Integer       serial         Atom serial number.
// 13 - 16       Atom          name           Atom name.
// 17            Character     altLoc         Alternate location indicator
// 18 - 20       Residue name  resName        Residue name.
// 22            Character     chainID        Chain identifier.
// 23 - 26       Integer       resSeq         Residue sequence number.
// 27            AChar         iCode          Insertion code.
// 29 - 35       Integer       u[0][0]        U(1,1)
// 36 - 42       Integer       u[1][1]        U(2,2)
// 43 - 49       Integer       u[2][2]        U(3,3)
// 50 - 56       Integer       u[0][1]        U(1,2)
// 57 - 63       Integer       u[0][2]        U(1,3)
// 64 - 70       Integer       u[1][2]        U(2,3)
// 77 - 78       LString(2)    element        Element symbol, right-justified.
// 79 - 80       LString(2)    charge         Charge on the atom.
// =============================================================================
        double[] anisou = atom.getAnisou(null);
        if (anisou != null) {
            anisouSB.replace(6, 80, sb.substring(6, 80));
            anisouSB.replace(28, 70, String.format("%7d%7d%7d%7d%7d%7d",
                    (int) (anisou[0] * 1e4), (int) (anisou[1] * 1e4),
                    (int) (anisou[2] * 1e4), (int) (anisou[3] * 1e4),
                    (int) (anisou[4] * 1e4), (int) (anisou[5] * 1e4)));
            if (!listMode) {
                bw.write(anisouSB.toString());
                bw.newLine();
            } else {
                listOutput.add(anisouSB.toString());
            }
        }
    }

    private void writeSIFTAtom(Atom atom, int serial, StringBuilder sb,
                               StringBuilder anisouSB, BufferedWriter bw, String siftScore)
            throws IOException {
        String name = atom.getName();
        if (name.length() > 4) {
            name = name.substring(0, 4);
        } else if (name.length() == 1) {
            name = name + "  ";
        } else if (name.length() == 2) {
            if (atom.getAtomType().valence == 0) {
                name = name + "  ";
            } else {
                name = name + " ";
            }
        }
        double xyz[] = vdwH ? atom.getRedXYZ() : atom.getXYZ(null);
        if (nSymOp != 0) {
            Crystal crystal = activeMolecularAssembly.getCrystal();
            SymOp symOp = crystal.spaceGroup.getSymOp(nSymOp);
            double[] newXYZ = new double[xyz.length];
            crystal.applySymOp(xyz, newXYZ, symOp);
            xyz = newXYZ;
        }
        sb.replace(6, 16, String.format("%5s " + padLeft(name.toUpperCase(), 4), Hybrid36.encode(5, serial)));
        Character altLoc = atom.getAltLoc();
        if (altLoc != null) {
            sb.setCharAt(16, altLoc);
        } else {
            sb.setCharAt(16, ' ');
        }
        if (siftScore == null) {
            sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                    xyz[0], xyz[1], xyz[2], atom.getOccupancy(), 110.0));
        } else {
            sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                    xyz[0], xyz[1], xyz[2], atom.getOccupancy(), (1 + (-1 * Float.parseFloat(siftScore))) * 100));
        }
        name = Atom.ElementSymbol.values()[atom.getAtomicNumber() - 1].toString();
        name = name.toUpperCase();
        if (atom.isDeuterium()) {
            name = "D";
        }
        sb.replace(76, 78, padLeft(name, 2));
        sb.replace(78, 80, String.format("%2d", 0));
        if (!listMode) {
            bw.write(sb.toString());
            bw.newLine();
        } else {
            listOutput.add(sb.toString());
        }
// =============================================================================
//  1 - 6        Record name   "ANISOU"
//  7 - 11       Integer       serial         Atom serial number.
// 13 - 16       Atom          name           Atom name.
// 17            Character     altLoc         Alternate location indicator
// 18 - 20       Residue name  resName        Residue name.
// 22            Character     chainID        Chain identifier.
// 23 - 26       Integer       resSeq         Residue sequence number.
// 27            AChar         iCode          Insertion code.
// 29 - 35       Integer       u[0][0]        U(1,1)
// 36 - 42       Integer       u[1][1]        U(2,2)
// 43 - 49       Integer       u[2][2]        U(3,3)
// 50 - 56       Integer       u[0][1]        U(1,2)
// 57 - 63       Integer       u[0][2]        U(1,3)
// 64 - 70       Integer       u[1][2]        U(2,3)
// 77 - 78       LString(2)    element        Element symbol, right-justified.
// 79 - 80       LString(2)    charge         Charge on the atom.
// =============================================================================
        double[] anisou = atom.getAnisou(null);
        if (anisou != null) {
            anisouSB.replace(6, 80, sb.substring(6, 80));
            anisouSB.replace(28, 70, String.format("%7d%7d%7d%7d%7d%7d",
                    (int) (anisou[0] * 1e4), (int) (anisou[1] * 1e4),
                    (int) (anisou[2] * 1e4), (int) (anisou[3] * 1e4),
                    (int) (anisou[4] * 1e4), (int) (anisou[5] * 1e4)));
            if (!listMode) {
                bw.write(anisouSB.toString());
                bw.newLine();
            } else {
                listOutput.add(anisouSB.toString());
            }
        }
    }

    @Override
    public int getSnapshot() {
        return modelsRead;
    }

    @Override
    public int countNumModels(){
        Pattern model = Pattern.compile("MODEL");
        int numModels = 0;
        boolean eof = true;
        for (MolecularAssembly system : systems) {
            try {
                BufferedReader currentReader;
                if (readers.containsKey(system)) {
                    currentReader = readers.get(system);
                    if (!currentReader.ready()) {
                        currentReader = new BufferedReader(new FileReader(readFile));
                        readers.put(system, currentReader);
                    }
                } else {
                    currentReader = new BufferedReader(new FileReader(readFile));
                    readers.put(system, currentReader);
                }
                // Skip to appropriate model.
                String line = currentReader.readLine();
                while (line != null) {
                    line = line.trim();
                    Matcher match = model.matcher(line);
                    if (match.find()) {
                        numModels++;
                        eof = false;
                    }
                    line = currentReader.readLine();
                }
                if (eof) {
                    logger.log(Level.INFO, String.format(" End of file reached for %s", readFile));
                    currentReader.close();
                    return numModels;
                }
            } catch (IOException ex) {
                logger.info(String.format(" Exception in parsing frame %d of %s:"
                        + " %s", modelsRead, system.toString(), ex.toString()));
            }
        }
        return numModels;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readNext() {
        return readNext(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readNext(boolean resetPosition) {
        return readNext(resetPosition, false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readNext(boolean resetPosition, boolean print) {
        // ^ is beginning of line, \\s+ means "one or more whitespace", (\\d+) means match and capture one or more digits.
        Pattern modelPatt = Pattern.compile("^MODEL\\s+(\\d+)");
        modelsRead = resetPosition ? 1 : modelsRead + 1;
        boolean eof = true;
        for (MolecularAssembly system : systems) {
            try {
                BufferedReader currentReader;
                if (readers.containsKey(system)) {
                    currentReader = readers.get(system);
                    if (!currentReader.ready()) {
                        currentReader = new BufferedReader(new FileReader(readFile));
                        readers.put(system, currentReader);
                    }
                } else {
                    currentReader = new BufferedReader(new FileReader(readFile));
                    readers.put(system, currentReader);
                }
                // Skip to appropriate model.
                String line = currentReader.readLine();
                while (line != null) {
                    line = line.trim();
                    Matcher m = modelPatt.matcher(line);
                    if (m.find()) {
                        int modelNum = Integer.parseInt(m.group(1));
                        if (modelNum == modelsRead) {
                            if (print) {
                                logger.log(Level.INFO, String.format(" Reading model %d for %s", modelNum, currentFile));
                            }
                            eof = false;
                            break;
                        }
                    }
                    line = currentReader.readLine();
                }
                if (eof) {
                    logger.log(Level.INFO, String.format(" End of file reached for %s", readFile));
                    currentReader.close();
                    return false;
                }

                // Begin parsing the model.
                boolean modelDone = false;
                line = currentReader.readLine();
                while (line != null) {
                    line = line.trim();
                    String recID = line.substring(0, Math.min(6, line.length())).trim();
                    try {
                        Record record = Record.valueOf(recID);
                        boolean hetatm = true;
                        switch (record) {
// =============================================================================
//
//  7 - 11        Integer       serial       Atom serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
// =============================================================================
//         1         2         3         4         5         6         7
//123456789012345678901234567890123456789012345678901234567890123456789012345678
//ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
//ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
                            case ATOM:
                                hetatm = false;
                            case HETATM:
                                if (!line.substring(17, 20).trim().equals("HOH")) {
                                    //int serial = Hybrid36.decode(5, line.substring(6, 11));
                                    String name = line.substring(12, 16).trim();
                                    if (name.toUpperCase().contains("1H") || name.toUpperCase().contains("2H")
                                            || name.toUpperCase().contains("3H")) {
                                        // VERSION3_2 is presently just a placeholder for "anything non-standard".
                                        fileStandard = VERSION3_2;
                                    }
                                    Character altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                                    if (!altLoc.equals(' ') && !altLoc.equals('A')
                                            && !altLoc.equals(currentAltLoc)) {
                                        break;
                                    }
                                    String resName = line.substring(17, 20).trim();
                                    Character chainID = line.substring(21, 22).charAt(0);

                                    List<String> segIDList = segidMap.get(chainID);
                                    if (segIDList == null) {
                                        logger.log(Level.WARNING, String.format(" No "
                                                + "known segment ID corresponds to "
                                                + "chain ID %s", chainID.toString()));
                                        break;
                                    }

                                    String segID = segIDList.get(0);
                                    if (segIDList.size() > 1) {
                                        logger.log(Level.WARNING, String.format(" "
                                                        + "Multiple segment IDs correspond to"
                                                        + "chain ID %s; assuming %s",
                                                chainID.toString(), segID));
                                    }

                                    int resSeq = Hybrid36.decode(4, line.substring(22, 26));

                                    double[] d = new double[3];
                                    d[0] = Double.valueOf(line.substring(30, 38).trim());
                                    d[1] = Double.valueOf(line.substring(38, 46).trim());
                                    d[2] = Double.valueOf(line.substring(46, 54).trim());
                                    double occupancy = 1.0;
                                    double tempFactor = 1.0;
                                    Atom newAtom = new Atom(0, name, altLoc, d, resName, resSeq,
                                            chainID, occupancy, tempFactor, segID);
                                    newAtom.setHetero(hetatm);
                                    // Check if this is a modified residue.
                                    if (modres.containsKey(resName.toUpperCase())) {
                                        newAtom.setModRes(true);
                                    }

                                    Atom returnedAtom = activeMolecularAssembly.findAtom(newAtom);
                                    if (returnedAtom != null) {
                                        returnedAtom.setXYZ(d);
                                        double[] retXYZ = new double[3];
                                        returnedAtom.getXYZ(retXYZ);
                                    } else {
                                        String message = String.format(" "
                                                        + "Could not find atom %s in assembly",
                                                newAtom.toString());
                                        if (dieOnMissingAtom) {
                                            logger.severe(message);
                                        } else {
                                            logger.warning(message);
                                        }
                                    }
                                    break;
                                }
                            case ENDMDL:
                            case END: // Technically speaking, END should be at the end of the file, not end of the model.
                                logger.log(Level.FINE, String.format(" Model %d successfully read", modelsRead));
                                modelDone = true;
                            default:
                                break;
                        }
                    } catch (Exception ex) {
                        // Do nothing; it's not an ATOM/HETATM line.
                    }
                    if (modelDone) {
                        break;
                    }
                    line = currentReader.readLine();

                }
                return true;
            } catch (IOException ex) {
                logger.info(String.format(" Exception in parsing frame %d of %s:"
                        + " %s", modelsRead, system.toString(), ex.toString()));
            }
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void closeReader() {
        // Java 8 stuff that Netbeans suggested. Faster than for loop?
        systems.stream().forEach((system) -> {
            BufferedReader br = readers.get(system);
            try {
                br.close();
            } catch (IOException ex) {
                logger.warning(String.format(" Exception in closing system %s: %s", system.toString(), ex.toString()));
            }
        });
    }

    /**
     * Ensures proper naming of hydrogens according to latest PDB format.
     * Presently mostly guesses at which hydrogens to re-assign, which may cause
     * chirality errors for prochiral hydrogens. If necessary, we will implement
     * more specific mapping.
     *
     * @param residue      Residue to examine.
     * @param fileStandard PDB File Standard to use.
     */
    public static void checkHydrogenAtomNames(Residue residue, PDBFileStandard fileStandard) {
        switch (fileStandard) {
            case VERSION3_3:
                return;
            case VERSION3_2:
            default:
                break;
        }
        // May have to get position.
        String residueType = residue.getName().toUpperCase();
        ArrayList<Atom> resAtoms = residue.getAtomList();
        for (Atom atom : resAtoms) {
            if (atom == null) {
                continue;
            }
            String atomName = atom.getName().toUpperCase();
            // Handles situations such as 1H where it should be NA_H1, etc.
            if (atomName.contains("H")) {
                try {
                    String firstChar = atomName.substring(0, 1);
                    Integer.parseInt(firstChar);
                    atomName = atomName.substring(1);
                    atomName = atomName.concat(firstChar);
                    atom.setName(atomName);
                } catch (NumberFormatException e) {
                    // Do nothing.
                }
            }
        }
        // Ensures proper hydrogen assignment; for example, Gln should have HB2,
        // HB3 instead of HB1, HB2.
        ArrayList<Atom> betas;
        ArrayList<Atom> gammas;
        ArrayList<Atom> deltas;
        ArrayList<Atom> epsilons;
        ArrayList<Atom> zetas;
        String atomName;
        Atom OH;
        Atom HH;
        Atom HG;
        Atom HD2;
        switch (getAminoAcid(residueType)) {
            case GLY:
                ArrayList<Atom> alphas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HA")) {
                        alphas.add(atom);
                    }
                }
                renameGlycineAlphaHydrogens(residue, alphas);
                break;
            case ALA:
                // No known errors with alanine
                break;
            case VAL:
                // No known errors with valine
                break;
            case LEU:
                betas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HB")) {
                        betas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case ILE:
                ArrayList<Atom> ileAtoms = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HG1")) {
                        ileAtoms.add(atom);
                    }
                }
                renameIsoleucineHydrogens(residue, ileAtoms);
                break;
            case SER:
                betas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HB")) {
                        betas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case THR:
                Atom HG1 = (Atom) residue.getAtomNode("HG1");
                if (HG1 == null) {
                    for (Atom atom : resAtoms) {
                        atomName = atom.getName().toUpperCase();
                        // Gets first HG-containing name of length < 4
                        // Length < 4 avoids bringing in HG21, HG22, or HG23.
                        if (atomName.length() < 4 && atomName.contains("HG")) {
                            atom.setName("HG1");
                            break;
                        }
                    }
                }
                break;
            case CYS:
                betas = new ArrayList<>();
                HG = (Atom) residue.getAtomNode("HG");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (HG == null && atomName.contains("HG")) {
                        HG = atom;
                        HG.setName("HG");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case CYX:
                // I pray this is never important, because I don't have an example CYX to work from.
                break;
            case CYD:
                betas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HB")) {
                        betas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case PRO:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                deltas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                renameDeltaHydrogens(residue, deltas, 23);
                break;
            case PHE:
                betas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                Atom HZ = (Atom) residue.getAtomNode("HZ");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (HZ == null && atomName.contains("HZ")) {
                        HZ = atom;
                        HZ.setName("HZ");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameDeltaHydrogens(residue, deltas, 12);
                renameEpsilonHydrogens(residue, epsilons, 12);
                break;
            case TYR:
                betas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                HH = (Atom) residue.getAtomNode("HH");
                OH = (Atom) residue.getAtomNode("OH");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (HH == null && atomName.contains("HH")) {
                        HH = atom;
                        HH.setName("HH");
                    } else if (OH == null && atomName.contains("O") && atomName.contains("H")) {
                        OH = atom;
                        OH.setName("OH");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameDeltaHydrogens(residue, deltas, 12);
                renameEpsilonHydrogens(residue, epsilons, 12);
                break;
            case TYD:
                betas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                OH = (Atom) residue.getAtomNode("OH");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (OH == null && atomName.contains("O") && atomName.contains("H")) {
                        OH = atom;
                        OH.setName("OH");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameDeltaHydrogens(residue, deltas, 12);
                renameEpsilonHydrogens(residue, epsilons, 12);
                break;
            case TRP:
                betas = new ArrayList<>();
                epsilons = new ArrayList<>();
                zetas = new ArrayList<>();
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HH2 = (Atom) residue.getAtomNode("HH2");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (atomName.contains("HZ")) {
                        zetas.add(atom);
                    } else if (HD1 == null && atomName.contains("HD")) {
                        HD1 = atom;
                        HD1.setName("HD1");
                    } else if (HH2 == null && atomName.contains("HH")) {
                        HH2 = atom;
                        HH2.setName("HH2");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameEpsilonHydrogens(residue, epsilons, 13);
                renameZetaHydrogens(residue, zetas, 23);
                break;
            case HIS:
                betas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameDeltaHydrogens(residue, deltas, 12);
                renameEpsilonHydrogens(residue, epsilons, 12);
                break;
            case HID:
                betas = new ArrayList<>();
                deltas = new ArrayList<>();
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (HE1 == null && atomName.contains("HE")) {
                        HE1 = atom;
                        HE1.setName("HE1");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameDeltaHydrogens(residue, deltas, 12);
                break;
            case HIE:
                betas = new ArrayList<>();
                epsilons = new ArrayList<>();
                HD2 = (Atom) residue.getAtomNode("HD2");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (HD2 == null && atomName.contains("HD")) {
                        HD2 = atom;
                        HD2.setName("HD2");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameEpsilonHydrogens(residue, epsilons, 12);
                break;
            case ASP:
                betas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    if (atom.getName().toUpperCase().contains("HB")) {
                        betas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case ASH:
                betas = new ArrayList<>();
                HD2 = (Atom) residue.getAtomNode("HD2");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (HD2 == null && atomName.contains("HD")) {
                        HD2 = atom;
                        HD2.setName("HD2");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                break;
            case ASN:
                betas = new ArrayList<>();
                ArrayList<Atom> HD2s = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HD")) {
                        HD2s.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameAsparagineHydrogens(residue, HD2s);
                break;
            case GLU:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                break;
            case GLH:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (HE2 == null && atomName.contains("HE")) {
                        HE2 = atom;
                        HE2.setName("HE2");
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                break;
            case GLN:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                epsilons = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                renameGlutamineHydrogens(residue, epsilons);
                break;
            case MET:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                // Epsilons should not break, as they are 1-3.
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                break;
            case LYS:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                // Zetas are 1-3, should not break.
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                renameDeltaHydrogens(residue, deltas, 23);
                renameEpsilonHydrogens(residue, epsilons, 23);
                break;
            case LYD:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                deltas = new ArrayList<>();
                epsilons = new ArrayList<>();
                zetas = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (atomName.contains("HE")) {
                        epsilons.add(atom);
                    } else if (atomName.contains("HZ")) {
                        zetas.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                renameDeltaHydrogens(residue, deltas, 23);
                renameEpsilonHydrogens(residue, epsilons, 23);
                renameZetaHydrogens(residue, zetas, 12);
                break;
            case ARG:
                betas = new ArrayList<>();
                gammas = new ArrayList<>();
                deltas = new ArrayList<>();
                Atom HE = (Atom) residue.getAtomNode("HE");
                ArrayList<Atom> HHn = new ArrayList<>();
                for (Atom atom : resAtoms) {
                    atomName = atom.getName().toUpperCase();
                    if (atomName.contains("HB")) {
                        betas.add(atom);
                    } else if (atomName.contains("HG")) {
                        gammas.add(atom);
                    } else if (atomName.contains("HD")) {
                        deltas.add(atom);
                    } else if (HE == null && atomName.contains("HE")) {
                        HE = atom;
                        HE.setName("HE");
                    } else if (atomName.contains("HH")) {
                        HHn.add(atom);
                    }
                }
                renameBetaHydrogens(residue, betas, 23);
                renameGammaHydrogens(residue, gammas, 23);
                renameDeltaHydrogens(residue, deltas, 23);
                renameArginineHydrogens(residue, HHn);
                break;
            case ORN:
            case AIB:
            case PCA:
            case UNK:
            default:
                // I am currently unaware of how these amino acids are typically
                // labeled under older PDB standards.
                break;
        }
    }

    private Atom buildHeavy(MSGroup residue, String atomName, Atom bondedTo, int key)
            throws MissingHeavyAtomException {
        try {
            return BondedUtils.buildHeavy(residue, atomName, bondedTo, key, forceField, bondList);
        } catch (MissingHeavyAtomException ex) {
            bondedTo.removeFromParent();
            residue.destroy();
            return null;
        }
    }

    private Atom buildHydrogen(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
                               Atom ic, double angle2, int chiral, int lookUp) {
        return BondedUtils.buildHydrogen(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, lookUp, forceField, bondList);
    }

    private Bond buildBond(Atom a1, Atom a2) {
        return BondedUtils.buildBond(a1, a2, forceField, bondList);
    }

    /**
     * Determine the atom type based on a biotype key.
     *
     * @param key The biotype key.
     * @return The atom type.
     * @since 1.0
     */
    private AtomType findAtomType(int key) {
        return BondedUtils.findAtomType(key, forceField);
    }

    /**
     * Convert possibly duplicate chain IDs into unique segIDs.
     *
     * @param c chain ID just read.
     * @return a unique segID.
     */
    private String getSegID(Character c) {
        if (c.equals(' ')) {
            c = 'A';
        }

        // If the chain ID has not changed, return the existing segID.
        if (c.equals(currentChainID)) {
            return currentSegID;
        }

        // Loop through existing segIDs to find the first one that is unused.
        int n = segIDs.size();
        int count = 0;
        for (int i = 0; i < n; i++) {
            String segID = segIDs.get(i);
            if (segID.endsWith(c.toString())) {
                count++;
            }
        }

        // If the count is greater than 0, then append it.
        String newSegID;
        if (count == 0) {
            newSegID = c.toString();
        } else {
            newSegID = Integer.toString(count) + c.toString();
        }

        segIDs.add(newSegID);
        currentChainID = c;
        currentSegID = newSegID;

        if (segidMap.containsKey(c)) {
            segidMap.get(c).add(newSegID);
        } else {
            List<String> newChainList = new ArrayList<>();
            newChainList.add(newSegID);
            segidMap.put(c, newChainList);
        }

        return newSegID;
    }

    /**
     * Finds Atoms bonded to a given Atom that match a certain atomic number.
     *
     * @param atom    Atom to search from.
     * @param element Atomic number to search for.
     * @return Bonded atoms of an element.
     */
    private static List<Atom> findBondedAtoms(Atom atom, int element) {
        return findBondedAtoms(atom, null, element);
    }

    /**
     * Finds Atoms bonded to a given Atom that match a certain atomic number that do not match an excluded atom.
     *
     * @param atom      Atom to search from.
     * @param toExclude Atom to exclude from search.
     * @param element   Atomic number to search for.
     * @return Bonded atoms of an element.
     */
    private static List<Atom> findBondedAtoms(Atom atom, Atom toExclude, int element) {
        return atom.getBonds().stream().
                map((Bond b) -> b.get1_2(atom)).
                filter((Atom a) -> a != toExclude).
                filter((Atom a) -> a.getAtomType().atomicNumber == element).
                collect(Collectors.toList());
    }

    /**
     * Checks if there is an Atom of a given atomic number bonded to the provided Atom.
     *
     * @param atom    Atom to search from.
     * @param element Atomic number to search for.
     * @return If bonded atoms of given element exist.
     */
    private static boolean hasAttachedAtom(Atom atom, int element) {
        return atom.getBonds().stream().
                map((Bond b) -> b.get1_2(atom)).
                anyMatch((Atom a) -> a.getAtomType().atomicNumber == element);
    }

    /**
     * Checks if atom a1 is bonded to atom a2.
     *
     * @param a1 An Atom.
     * @param a2 Another Atom.
     * @return If a1 is bonded to a2.
     */
    private static boolean atomAttachedToAtom(Atom a1, Atom a2) {
        assert a1 != a2;
        return a1.getBonds().stream().
                anyMatch((Bond b) -> b.get1_2(a1) == a2);
    }

    /**
     * Finds all Atoms belonging to a Residue of a given atomic number.
     *
     * @param residue Residue to search in.
     * @param element Atomic number to search for.
     * @return A list of matching Atoms.
     */
    private static List<Atom> findAtomsOfElement(Residue residue, int element) {
        return residue.getAtomList().stream().
                filter((Atom a) -> a.getAtomType().atomicNumber == element).
                collect(Collectors.toList());
    }

    /**
     * Renames Atoms to PDB standard using bonding patterns, atomic numbers, and residue types.
     * <p>
     * Will not work if a force field definition botches its atomic numbers.
     * <p>
     * Only implemented for amino acids and nucleic acids at this time.
     *
     * @param assembly MolecularAssembly to fix.
     */
    public static void renameAtomsToPDBStandard(MolecularAssembly assembly) {
        Polymer[] polys = assembly.getChains();
        if (polys != null && polys.length > 0) {
            for (Polymer polymer : polys) {
                for (Residue residue : polymer.getResidues()) {
                    switch (residue.getResidueType()) {
                        case AA:
                            renameAminoAcidToPDBStandard(residue);
                            break;
                        case NA:
                            renameNucleicAcidToPDBStandard(residue);
                            break;
                        case UNK:
                        default:
                            break;
                    }
                }
            }
        }
    }

    /**
     * Finds the alpha carbon of a residue, and handles any C-terminal ACE caps while at it.
     *
     * @param residue Find the alpha carbon of.
     * @param N       The residue's backbone nitrogen.
     * @return The alpha carbon.
     */
    private static Atom getAlphaCarbon(Residue residue, Atom N) {
        List<Atom> resAtoms = residue.getAtomList();
        List<Atom> caCandidates = findBondedAtoms(N, 6).stream().
                filter((Atom carbon) -> resAtoms.contains(carbon)).
                collect(Collectors.toList());


        switch (residue.getAminoAcid3()) {
            case PRO: {
                Atom CA = null;
                Atom CD = null;
                Atom aceC = null;
                for (Atom caCand : caCandidates) {
                    if (hasAttachedAtom(caCand, 8)) {
                        aceC = caCand;
                    } else {
                        List<Atom> attachedH = findBondedAtoms(caCand, 1);
                        if (attachedH.size() == 1) {
                            CA = caCand;
                        } else if (attachedH.size() == 2) {
                            CD = caCand;
                        } else {
                            throw new IllegalArgumentException(String.format(" Error in parsing proline %s", residue));
                        }
                    }
                }
                assert CA != null && CD != null;
                if (aceC != null) {
                    nameAcetylCap(residue, aceC);
                }
                return CA;
            }
            default: {
                if (caCandidates.size() == 1) {
                    return caCandidates.get(0);
                } else {
                    Atom CA = null;
                    Atom aceC = null;
                    for (Atom caCand : caCandidates) {
                        if (hasAttachedAtom(caCand, 8)) {
                            aceC = caCand;
                        } else {
                            CA = caCand;
                        }
                    }
                    nameAcetylCap(residue, aceC);
                    return CA;
                }
            }
        }
    }

    /**
     * Names the atoms in an N-terminal acetyl ACE capping group.
     *
     * @param residue Residue containing an acetyl cap.
     * @param aceC    The acetyl group's C atom.
     */
    private static void nameAcetylCap(Residue residue, Atom aceC) {
        logger.warning(String.format(" Probable ACE cap attached to residue %s; duplicate atom names may result!", residue));
        aceC.setName("C");
        findBondedAtoms(aceC, 8).get(0).setName("O");
        Atom CH3 = findBondedAtoms(aceC, 6).get(0);
        CH3.setName("CH3");
        List<Atom> ntermHs = findBondedAtoms(CH3, 1);
        for (int i = 0; i < 3; i++) {
            ntermHs.get(i).setName(String.format("H%d", (i + 1)));
        }
    }

    /**
     * Renames the Atoms in an amino acid to PDB standard.
     *
     * @param residue Residue to fix atom names of.
     */
    private static void renameAminoAcidToPDBStandard(Residue residue) {
        if (residue.getChainID() == null) {
            residue.setChainID('Z');
        }
        final Atom N = findNitrogenAtom(residue);
        AminoAcid3 aa3 = residue.getAminoAcid3();
        if (N != null) {
            N.setName("N");

            Atom CA = getAlphaCarbon(residue, N);
            CA.setName("CA");

            List<Atom> has = findBondedAtoms(CA, 1);
            switch (aa3) {
                case NME:
                    // Do all renaming here then return out of the method.
                    findBondedAtoms(N, 1).get(0).setName("H");
                    CA.setName("CH3");
                    for (int i = 1; i <= 3; i++) {
                        has.get(i - 1).setName(String.format("H%d", i));
                    }
                    return;
                case GLY:
                    has.get(0).setName("HA2");
                    has.get(1).setName("HA3");
                    break;
                default:
                    has.get(0).setName("HA");
                    break;
            }

            Atom C = null;
            Atom CB = null;
            for (Atom carb : findBondedAtoms(CA, 6)) {
                // Second check is because of serine/threonine OG bonded straight to CB.
                if (hasAttachedAtom(carb, 8) && !hasAttachedAtom(carb, 1)) {
                    C = carb;
                    C.setName("C");
                } else {
                    CB = carb;
                    CB.setName("CB");
                }
            }
            if (C == null) {
                throw new IllegalArgumentException(String.format(" Could not find carbonyl carbon for residue %s!", residue));
            }
            if (CB == null && aa3 != AminoAcid3.GLY) {
                throw new IllegalArgumentException(String.format(" Could not find beta carbon for residue %s!", residue));
            }

            List<Atom> ctermOxygens = findBondedAtoms(C, 8);
            switch (ctermOxygens.size()) {
                case 1:
                    ctermOxygens.get(0).setName("O");
                    break;
                case 2:
                    Atom O = null;
                    for (Atom oxy : ctermOxygens) {
                        if (oxy.getBonds().size() == 2) {
                            O = oxy;
                            O.setName("O");
                            findBondedAtoms(O, 1).get(0).setName("HO");
                        }
                    }
                    if (O == null) {
                        ctermOxygens.get(0).setName("O");
                        ctermOxygens.get(1).setName("OXT");
                    }
            }

            List<Atom> amideProtons = findBondedAtoms(N, 1);
            switch (amideProtons.size()) {
                case 1:
                    amideProtons.get(0).setName("H");
                    break;
                default:
                    // Should catch both N-termini and proline.
                    for (int i = 1; i <= amideProtons.size(); i++) {
                        amideProtons.get(i - 1).setName(String.format("H%d", i));
                    }
                    break;
            }

            /**
             * All common atoms are now named: N, H[1-3], CA, HA[2-3], CB, C, O[XT], [HO]
             */
            renameCommonAminoAcids(residue, aa3, CA, CB);
        } else {
            switch (aa3) {
                case ACE: {
                    Atom O = findAtomsOfElement(residue, 8).get(0);
                    O.setName("O");
                    Atom C = findBondedAtoms(O, 6).get(0);
                    C.setName("C");
                    Atom CH3 = findBondedAtoms(C, 6).get(0);
                    CH3.setName("CH3");
                    List<Atom> hydrogens = findBondedAtoms(CH3, 1);
                    for (int i = 1; i <= 3; i++) {
                        hydrogens.get(i - 1).setName(String.format("H%d", i));
                    }
                }
                break;
                default:
                    throw new IllegalArgumentException(String.format(" Could not find nitrogen atom for residue %s!", residue));
            }
        }
    }

    /**
     * Renames atoms in common amino acids to PDB standard.
     *
     * @param residue Residue to perform renaming for.
     * @param aa3     Its AA3 code.
     * @param CA      Its alpha carbon.
     * @param CB      Its beta carbon.
     */
    private static void renameCommonAminoAcids(Residue residue, AminoAcid3 aa3, Atom CA, Atom CB) {
        switch (aa3) {
            case ALA: {
                renameAlkyl(CB, CA, 1, 'B');
            }
            break;
            case CYS:
            case CYD: {
                Atom SG = renameAlkyl(CB, CA, 2, 'B').get();
                SG.setName("SG");
                if (hasAttachedAtom(SG, 1)) {
                    assert aa3 == AminoAcid3.CYS;
                    findBondedAtoms(SG, 1).get(0).setName("HG");
                } else if (hasAttachedAtom(SG, 16)) {
                    logger.finer(String.format(" SG atom %s likely part of a disulfide bond.", SG));
                } else {
                    residue.setName("CYD");
                }
            }
            break;
            case ASP:
            case ASH: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                List<Atom> ODs = findBondedAtoms(CG, 8);

                int protonatedOD = -1; // -1: Deprotonated ASP. 0/1: Index of protonated oxygen (ASH).
                for (int i = 0; i < 2; i++) {
                    if (hasAttachedAtom(ODs.get(i), 1)) {
                        protonatedOD = i;
                        break;
                    }
                }

                switch (protonatedOD) {
                    case -1:
                        ODs.get(0).setName("OD1");
                        ODs.get(1).setName("OD2");
                        break;
                    case 0:
                        if (aa3 != AminoAcid3.ASH) {
                            residue.setName("ASH");
                        }
                        ODs.get(0).setName("OD2");
                        findBondedAtoms(ODs.get(0), 1).get(0).setName("HD2");
                        ODs.get(1).setName("OD1");
                        break;
                    case 1:
                        if (aa3 != AminoAcid3.ASH) {
                            residue.setName("ASH");
                        }
                        ODs.get(1).setName("OD2");
                        findBondedAtoms(ODs.get(1), 1).get(0).setName("HD2");
                        ODs.get(0).setName("OD1");
                        break;
                }
            }
            break;
            case GLU:
            case GLH: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
                CD.setName("CD");
                List<Atom> OEs = findBondedAtoms(CD, 8);

                int protonatedOE = -1; // If it remains -1, deprotonated ASP, else ASH.
                for (int i = 0; i < 2; i++) {
                    if (hasAttachedAtom(OEs.get(i), 1)) {
                        protonatedOE = i;
                        break;
                    }
                }

                switch (protonatedOE) {
                    case -1:
                        OEs.get(0).setName("OE1");
                        OEs.get(1).setName("OE2");
                        break;
                    case 0:
                        if (aa3 != AminoAcid3.GLH) {
                            residue.setName("GLH");
                        }
                        OEs.get(0).setName("OE2");
                        findBondedAtoms(OEs.get(0), 1).get(0).setName("HE2");
                        OEs.get(1).setName("OE1");
                        break;
                    case 1:
                        if (aa3 != AminoAcid3.GLH) {
                            residue.setName("GLH");
                        }
                        OEs.get(1).setName("OE2");
                        findBondedAtoms(OEs.get(1), 1).get(0).setName("HE2");
                        OEs.get(0).setName("OE1");
                        break;
                }
            }
            break;
            case PHE: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                List<Atom> CDs = findBondedAtoms(CG, CB, 6);

                Atom CZ = null;
                for (int i = 1; i <= 2; i++) {
                    Atom CD = CDs.get(i - 1);
                    Atom CE = renameBranchedAlkyl(CD, CG, 0, i, 'D').get();
                    CZ = renameBranchedAlkyl(CE, CD, 0, i, 'E').get();
                }
                CZ.setName("CZ");
                findBondedAtoms(CZ, 1).get(0).setName("HZ");
            }
            break;
            case GLY:
                break;
            case HIS:
            case HIE:
            case HID: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");

                Atom CD2 = findBondedAtoms(CG, 6).stream().
                        filter((Atom a) -> a != CB).
                        findAny().get();
                CD2.setName("CD2");
                findBondedAtoms(CD2, 1).get(0).setName("HD2");

                Atom NE2 = findBondedAtoms(CD2, 7).get(0);
                NE2.setName("NE2");
                List<Atom> HE2 = findBondedAtoms(NE2, 1);
                boolean epsProtonated = (HE2 != null && !HE2.isEmpty());
                if (epsProtonated) {
                    HE2.get(0).setName("HE2");
                }

                Atom CE1 = findBondedAtoms(NE2, CD2, 6).get(0);
                CE1.setName("CE1");
                findBondedAtoms(CE1, 1).get(0).setName("HE1");

                Atom ND1 = findBondedAtoms(CG, 7).get(0);
                ND1.setName("ND1");
                List<Atom> HD1 = findBondedAtoms(ND1, 1);
                boolean deltaProtonated = (HD1 != null && !HD1.isEmpty());
                if (deltaProtonated) {
                    HD1.get(0).setName("HD1");
                }

                // All constant atoms found: now check protonation state.
                if (epsProtonated && deltaProtonated) {
                    assert aa3 == AminoAcid3.HIS;
                } else if (epsProtonated) {
                    residue.setName("HIE");
                } else if (deltaProtonated) {
                    residue.setName("HID");
                } else {
                    throw new IllegalArgumentException(String.format(" Histidine residue %s is doubly deprotonated!", residue));
                }
            }
            break;
            case ILE: {
                findBondedAtoms(CB, 1).get(0).setName("HB");
                List<Atom> CGs = findBondedAtoms(CB, CA, 6);

                for (Atom CG : CGs) {
                    List<Atom> HGs = findBondedAtoms(CG, 1);
                    int numHGs = HGs.size();
                    if (numHGs == 3) {
                        renameBranchedAlkyl(CG, CB, 1, 2, 'G');
                    } else if (numHGs == 2) {
                        Atom CD1 = renameBranchedAlkyl(CG, CB, 2, 1, 'G').get();
                        renameBranchedAlkyl(CD1, CG, 1, 1, 'D');
                    } else {
                        throw new IllegalArgumentException(String.format(" Isoleucine residue %s had %d gamma hydrogens, expecting 2-3!", residue, numHGs));
                    }
                }
            }
            break;
            case LYS:
            case LYD: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
                Atom CE = renameAlkyl(CD, CG, 2, 'D').get();
                Atom NZ = renameAlkyl(CE, CD, 2, 'E').get();
                // For a very brief period, NZ will be named CZ.
                renameAlkyl(NZ, CE, 1, 'Z');
                NZ.setName("NZ");
                int numH = findBondedAtoms(NZ, 1).size();
                switch (numH) {
                    case 2:
                        residue.setName("LYD");
                        break;
                    case 3:
                        assert aa3 == AminoAcid3.LYS;
                        break;
                    default:
                        throw new IllegalArgumentException(String.format(" Lysine residue %s had %d amine protons, expecting 2-3!", residue, numH));
                }
            }
            break;
            case LEU: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                findBondedAtoms(CG, 1).get(0).setName("HG");
                List<Atom> CDs = findBondedAtoms(CG, CB, 6);

                for (int i = 0; i < 2; i++) {
                    renameBranchedAlkyl(CDs.get(i), CG, 1, (i + 1), 'D');
                }
            }
            break;
            case MET: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom SD = renameAlkyl(CG, CB, 2, 'G').get();
                Atom CE = renameAlkyl(SD, CG, 0, 'D').get();
                // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
                SD.setName("SD");
                renameAlkyl(CE, SD, 1, 'E');
            }
            break;
            case ASN: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                findBondedAtoms(CG, 8).get(0).setName("OD1");
                Atom ND2 = findBondedAtoms(CG, 7).get(0);
                renameBranchedAlkyl(ND2, CG, 1, 2, 'D');
                // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
                ND2.setName("ND2");
            }
            break;
            case PRO: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
                Atom N = renameAlkyl(CD, CG, 2, 'D').get();
                assert N.getName().equals("N");
            }
            break;
            case GLN: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
                CD.setName("CD");

                findBondedAtoms(CD, 8).get(0).setName("OE1");
                Atom NE2 = findBondedAtoms(CD, 7).get(0);
                renameBranchedAlkyl(NE2, CD, 1, 2, 'E');
                // Once again, briefly misnamed atom because I'm kludging it through renameAlkyl.
                NE2.setName("NE2");
            }
            break;
            case ARG: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                Atom CD = renameAlkyl(CG, CB, 2, 'G').get();
                Atom NE = renameAlkyl(CD, CG, 2, 'D').get();
                Atom CZ = renameAlkyl(NE, CD, 0, 'E').get();
                NE.setName("NE");
                CZ.setName("CZ");

                List<Atom> NHs = findBondedAtoms(CZ, NE, 7);
                assert NHs.size() == 2;
                for (int i = 0; i < 2; i++) {
                    Atom NHx = NHs.get(i);
                    renameBranchedAlkyl(NHx, CZ, 1, (i + 1), 'H');
                    NHx.setName(String.format("NH%d", (i + 1)));
                }
            }
            break;
            case SER: {
                Atom OG = renameAlkyl(CB, CA, 2, 'B').get();
                renameAlkyl(OG, CB, 0, 'G');
                OG.setName("OG");
            }
            break;
            case THR: {
                CB.setName("CB"); // Should be unnecessary.
                findBondedAtoms(CB, 1).get(0).setName("HB");

                Atom OG1 = findBondedAtoms(CB, 8).get(0);
                OG1.setName("OG1");
                findBondedAtoms(OG1, 1).get(0).setName("HG1");

                Atom CG2 = findBondedAtoms(CB, CA, 6).get(0);
                renameBranchedAlkyl(CG2, CB, 1, 2, 'G');
            }
            break;
            case VAL: {
                CB.setName("CB"); // Should be unnecessary.
                findBondedAtoms(CB, 1).get(0).setName("HB");

                List<Atom> CGs = findBondedAtoms(CB, CA, 6);

                assert CGs.size() == 2;
                for (int i = 0; i < 2; i++) {
                    Atom CGx = CGs.get(i);
                    renameBranchedAlkyl(CGx, CB, 1, (i + 1), 'G');
                }
            }
            break;
            case TRP: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                List<Atom> CDs = findBondedAtoms(CG, CB, 6);
                Atom CD1 = null;
                Atom CD2 = null;

                for (Atom CDx : CDs) {
                    if (hasAttachedAtom(CDx, 1)) {
                        CD1 = CDx;
                    } else {
                        CD2 = CDx;
                        CD2.setName("CD2");
                    }
                }
                Atom NE1 = renameBranchedAlkyl(CD1, CG, 0, 1, 'D').get();
                Atom CE2 = renameBranchedAlkyl(NE1, CD1, 0, 1, 'E').get();
                NE1.setName("NE1");
                CE2.setName("CE2");

                Atom CZ2 = findBondedAtoms(CE2, CD2, 6).get(0);
                Atom CH2 = renameBranchedAlkyl(CZ2, CE2, 0, 2, 'Z').get();
                Atom CZ3 = renameBranchedAlkyl(CH2, CZ2, 0, 2, 'H').get();
                Atom CE3 = renameBranchedAlkyl(CZ3, CH2, 0, 3, 'Z').get();
                if (CD2 != renameBranchedAlkyl(CE3, CZ3, 0, 3, 'E').get()) {
                    throw new IllegalArgumentException(String.format(" Error in cyclizing tryptophan %s!", residue));
                }
            }
            break;
            case TYR:
            case TYD: {
                Atom CG = renameAlkyl(CB, CA, 2, 'B').get();
                CG.setName("CG");
                List<Atom> CDs = findBondedAtoms(CG, CB, 6);
                Atom CZ = null;

                assert CDs.size() == 2;
                for (int i = 1; i <= 2; i++) {
                    Atom CDx = CDs.get(i - 1);
                    Atom CEx = renameBranchedAlkyl(CDx, CG, 0, i, 'D').get();
                    CZ = renameBranchedAlkyl(CEx, CDx, 0, i, 'E').get();
                }

                CZ.setName("CZ");
                Atom OH = findBondedAtoms(CZ, 8).get(0);
                OH.setName("OH");
                if (hasAttachedAtom(OH, 1)) {
                    assert aa3 == AminoAcid3.TYR;
                    findBondedAtoms(OH, 1).get(0).setName("HH");
                } else {
                    residue.setName("TYD");
                }
            }
            break;
            default:
                throw new IllegalArgumentException((String.format(" Amino acid %s not recognized!", residue)));
        }
    }

    /**
     * Renames an atom, its bonded hydrogens, and returns the next atom in the chain.
     * <p>
     * If applied to an atom that is not a carbon, it will be misnamed as a carbon, so fix that afterwards.
     *
     * @param carbon       Alkyl carbon to rename.
     * @param priorAtom    Prior atom in the chain.
     * @param protonOffset Number of the first hydrogen (such as 2 for HB2-3).
     * @param posName      Name of the position (such as B for CB).
     * @return Next atom in the chain if present.
     */
    private static Optional<Atom> renameAlkyl(Atom carbon, Atom priorAtom, int protonOffset, char posName) {
        carbon.setName(String.format("C%c", posName));
        List<Atom> hydrogens = findBondedAtoms(carbon, 1);
        int numH = hydrogens.size();
        if (numH == 1) {
            hydrogens.get(0).setName(String.format("H%c", posName));
        } else {
            for (int i = 0; i < numH; i++) {
                hydrogens.get(i).setName(String.format("H%c%d", posName, i + protonOffset));
            }
        }

        return carbon.getBonds().stream().
                map((Bond b) -> b.get1_2(carbon)).
                filter((Atom a) -> a != priorAtom).
                filter((Atom a) -> !hydrogens.contains(a)).
                findAny();
    }

    /**
     * Renames a numbered carbon, its bonded hydrogens, and returns the next atom in the chain.
     * <p>
     * If applied to an atom that is not a carbon, it will be misnamed as a carbon, so fix that afterwards.
     * <p>
     * This is for carbons like PHE CD1 and CD2.
     *
     * @param carbon       Alkyl carbon to rename.
     * @param priorAtom    Prior atom in the chain.
     * @param protonOffset Number of the first hydrogen (such as 2 for HB2-3).
     * @param branchNum    Index of the branch.
     * @param posName      Name of the position (such as B for CB).
     * @return Next atom in the chain if present.
     */
    private static Optional<Atom> renameBranchedAlkyl(Atom carbon, Atom priorAtom, int protonOffset, int branchNum, char posName) {
        carbon.setName(String.format("C%c%d", posName, branchNum));
        List<Atom> hydrogens = findBondedAtoms(carbon, 1);
        int numH = hydrogens.size();
        if (numH == 1) {
            hydrogens.get(0).setName(String.format("H%c%d", posName, branchNum));
        } else {
            for (int i = 0; i < numH; i++) {
                hydrogens.get(i).setName(String.format("H%c%d%d", posName, branchNum, i + protonOffset));
            }
        }

        return carbon.getBonds().stream().
                map((Bond b) -> b.get1_2(carbon)).
                filter((Atom a) -> a != priorAtom).
                filter((Atom a) -> !hydrogens.contains(a)).
                findAny();
    }

    /**
     * Finds the backbone nitrogen of a residue.
     *
     * @param residue Amino acid residue to search for.
     * @return backbone nitrogen.
     */
    private static Atom findNitrogenAtom(Residue residue) {
        assert residue.getResidueType() == Residue.ResidueType.AA;

        // Will filter out amide N from NME caps at the end of the method.
        List<Atom> nitrogenCandidates = new ArrayList<>(2);

        switch (residue.getAminoAcid3()) {
            case LYS:
            case LYD: {
                /**
                 * For lysine: find the nitrogen bonded to a carbon that does not have two protons.
                 */
                List<Atom> nitrogens = findAtomsOfElement(residue, 7);
                for (Atom nitrogen : nitrogens) {
                    List<Atom> carbons = findBondedAtoms(nitrogen, 6);
                    if (carbons.size() == 2) {
                        nitrogenCandidates.add(nitrogen);
                    } else if (findBondedAtoms(carbons.get(0), 1).size() < 2) {
                        nitrogenCandidates.add(nitrogen);
                    }
                }
                if (nitrogenCandidates.isEmpty()) {
                    throw new IllegalArgumentException(String.format(" Could not identify N atom of residue %s!", residue));
                }
            }
            break;
            // Arginine and histidine can be handled very similarly.
            case ARG:
            case HIS:
            case HIE:
            case HID: {
                /**
                 * Easiest to the carbon bonded to all the sidechain nitrogens,
                 * then find the nitrogen not thus bonded.
                 */
                List<Atom> nitrogens = findAtomsOfElement(residue, 7);
                Atom commonC = findAtomsOfElement(residue, 6).stream().
                        filter((Atom carbon) -> findBondedAtoms(carbon, 7).size() >= 2).
                        findAny().get();
                nitrogenCandidates = nitrogens.stream().
                        filter((Atom nitr) -> !atomAttachedToAtom(nitr, commonC)).
                        collect(Collectors.toList());
            }
            break;
            case ASN:
            case GLN: {
                /**
                 * Find a bonded carbon that is not bonded to an oxygen.
                 * Both N and ND/NE have an attached carbonyl carbon.
                 * Only N will have CA attached.
                 */
                List<Atom> nitrogens = findAtomsOfElement(residue, 7);
                for (Atom nitrogen : nitrogens) {
                    List<Atom> bondedCarbs = findBondedAtoms(nitrogen, 6);
                    for (Atom carbon : bondedCarbs) {
                        if (!hasAttachedAtom(carbon, 8)) {
                            nitrogenCandidates.add(nitrogen);
                        }
                    }
                }
                if (nitrogenCandidates.isEmpty()) {
                    throw new IllegalArgumentException(String.format(" Could not identify N atom of residue %s!", residue));
                }
            }
            break;
            case TRP: {
                /**
                 * For tryptophan:
                 * If at an N-terminus, there will be only one bonded carbon.
                 * Else, one carbon will be a carbonyl carbon.
                 */
                List<Atom> nitrogens = findAtomsOfElement(residue, 7);
                for (Atom nitrogen : nitrogens) {
                    List<Atom> bondedCarbs = findBondedAtoms(nitrogen, 6);
                    if (bondedCarbs.size() == 1) {
                        nitrogenCandidates.add(nitrogen);
                    }
                    for (Atom carbon : bondedCarbs) {
                        if (hasAttachedAtom(carbon, 8)) {
                            nitrogenCandidates.add(nitrogen);
                        }
                    }
                }
                if (nitrogenCandidates.isEmpty()) {
                    throw new IllegalArgumentException(String.format(" Could not identify N atom of residue %s!", residue));
                }
            }
            break;
            case ACE:
                return null;
            default:
                /**
                 * All others should only have one nitrogen atom.
                 */
                nitrogenCandidates = findAtomsOfElement(residue, 7);
                break;
        }

        switch (nitrogenCandidates.size()) {
            case 0:
                logger.warning(" Did not find any atoms that might be the amide nitrogen for residue " + residue.toString());
                return null;
            case 1:
                return nitrogenCandidates.get(0);
            case 2:
                logger.warning(String.format(" Probable NME C-terminal cap attached to residue %s, some atom names may be duplicated!", residue));
                Atom N = null;
                for (Atom nitro : nitrogenCandidates) {
                    nitro.setName("N");
                    Optional<Atom> capMethyl = findBondedAtoms(nitro, 6).stream().
                            filter((Atom carb) -> findBondedAtoms(carb, 1).size() == 3).
                            findAny();
                    if (capMethyl.isPresent()) {
                        findBondedAtoms(nitro, 1).get(0).setName("H");
                        Atom theCap = capMethyl.get();
                        theCap.setName("CH3");
                        List<Atom> capHydrogens = findBondedAtoms(theCap, 1);
                        for (int i = 0; i < 3; i++) {
                            capHydrogens.get(i).setName(String.format("H%d", i + 1));
                        }
                    } else {
                        N = nitro;
                    }
                }
                return N;
            default:
                throw new IllegalArgumentException(String.format(" Could not definitely identify amide nitrogen for residue %s", residue));
        }
    }

    /**
     * Sorts toCompare by distance to the reference Atom, returning a sorted array.
     *
     * @param reference Atom to compare distances to.
     * @param toCompare Atoms to sort by distance (not modified).
     * @return Sorted array of atoms in toCompare.
     */
    private static Atom[] sortAtomsByDistance(Atom reference, List<Atom> toCompare) {
        Atom[] theAtoms = toCompare.toArray(new Atom[toCompare.size()]);
        sortAtomsByDistance(reference, theAtoms);
        return theAtoms;
    }

    /**
     * In-place sorts toCompare by distance to the reference Atom. Modifies toCompare.
     *
     * @param reference Atom to compare distances to.
     * @param toCompare Atoms to sort (in-place) by distance.
     */
    private static void sortAtomsByDistance(Atom reference, Atom[] toCompare) {
        final double[] refXYZ = reference.getXYZ(new double[3]);
        Arrays.sort(toCompare, Comparator.comparingDouble(a -> {
            double[] atXYZ = a.getXYZ(new double[3]);
            return VectorMath.dist2(refXYZ, atXYZ);
        }));
    }

    private static void renameNucleicAcidToPDBStandard(Residue residue) {
        if (residue.getChainID() == null) {
            residue.setChainID('Z');
        }
        assert residue.getResidueType() == Residue.ResidueType.NA;
        NucleicAcid3 na3 = residue.getNucleicAcid3(true);
        residue.setName(na3.toString());
        // logger.info(residue.toString());
        switch (na3) {
            case ADE:
            case DAD:
            case CYT:
            case DCY:
            case GUA:
            case DGU:
            case THY:
            case DTY:
            case URI:
                renameCommonNucleicAcid(residue, na3);
                break;
            default:
                logger.info(" Could not rename atoms for nonstandard nucleic acid " + na3.toString());
                break;
        }
    }

    private static void renameCommonNucleicAcid(Residue residue, NucleicAcid3 na3) {
        Optional<Atom> optO4s = findNucleotideO4s(residue);
        if (optO4s.isPresent()) {
            // Name O4', which is the unique ether oxygen.
            Atom O4s = optO4s.get();
            O4s.setName("O4\'");

            // C1' is bonded to a nitrogen (at least for non-abasic sites), C4' isn't
            List<Atom> bondedC = findBondedAtoms(O4s, 6);
            Atom C4s = null;
            Atom C1s = null;
            // Will need the first base nitrogen (N1/N9), and H1' later.
            Atom N19 = null;
            Atom H1s = null;
            for (Atom c : bondedC) {
                if (hasAttachedAtom(c, 7)) {
                    C1s = c;
                    C1s.setName("C1\'");
                    H1s = findBondedAtoms(C1s, 1).get(0);
                    H1s.setName("H1\'");
                    N19 = findBondedAtoms(C1s, 7).get(0);
                } else {
                    C4s = c;
                    C4s.setName("C4\'");
                    findBondedAtoms(C4s, 1).get(0).setName("H4\'");
                }
            }
            assert C4s != null && C1s != null;

            Atom C2s = findBondedAtoms(C1s, 6).get(0);
            C2s.setName("C2\'");

            bondedC = findBondedAtoms(C4s, 6);
            Atom C5s = null;
            Atom C3s = null;
            Atom O3s = null;
            for (Atom c : bondedC) {
                if (c.getBonds().stream().anyMatch(b -> b.get1_2(c) == C2s)) {
                    C3s = c;
                    C3s.setName("C3\'");
                    O3s = findBondedAtoms(C3s, 8).get(0);
                    O3s.setName("O3\'");
                    findBondedAtoms(C3s, 1).get(0).setName("H3\'");
                    if (hasAttachedAtom(O3s, 1)) {
                        findBondedAtoms(O3s, 1).get(0).setName("HO3\'");
                    } // Else, handle the possibility of 3'-P cap later.
                } else {
                    C5s = c;
                    C5s.setName("C5\'");
                    List<Atom> allH5List = findBondedAtoms(C5s, 1);
                    Atom[] allH5s = allH5List.toArray(new Atom[allH5List.size()]);
                    sortAtomsByDistance(O4s, allH5s);
                    allH5s[0].setName("H5\'");
                    allH5s[1].setName("H5\'\'");
                }
            }

            if (hasAttachedAtom(C2s, 8)) {
                Atom O2s = findBondedAtoms(C2s, 8).get(0);
                O2s.setName("O2\'");
                findBondedAtoms(O2s, 1).get(0).setName("HO2\'");
                findBondedAtoms(C2s, 1).get(0).setName("H2\'");
            } else {
                List<Atom> bothH2List = findBondedAtoms(C2s, 1);
                Atom[] bothH2 = bothH2List.toArray(new Atom[bothH2List.size()]);
                sortAtomsByDistance(H1s, bothH2);
                // Best-guess assignment, but is sometimes the other way around.
                bothH2[0].setName("H2\'\'");
                bothH2[1].setName("H2\'");
            }

            // logger.info(String.format(" C5\' null: %b", C5s == null));
            Atom O5s = findBondedAtoms(C5s, 8).get(0);
            O5s.setName("O5\'");

            if (hasAttachedAtom(O5s, 1)) {
                findBondedAtoms(O5s, 1).get(0).setName("HO5\'");
            } else if (hasAttachedAtom(O5s, 15)) {
                Atom P = findBondedAtoms(O5s, 15).get(0);
                P.setName("P");
                List<Atom> bondedO = findBondedAtoms(P, O5s, 8);
                List<Atom> thisResO = bondedO.stream().
                        filter(o -> residue.getAtomList().contains(o)).
                        collect(Collectors.toList());
                int nBonded = bondedO.size();
                int nRes = thisResO.size();

                if (nBonded == 0) {
                    // Do nothing.
                } else if (nBonded == nRes) {
                    Atom OP1 = bondedO.get(0);
                    OP1.setName("OP1");

                    // OP2 is approximately +120 degrees from OP1, OP3 is -120 degrees.
                    final double[] xyzC5s = C5s.getXYZ(new double[3]);
                    final double[] xyzO5s = O5s.getXYZ(new double[3]);
                    final double[] xyzP = P.getXYZ(new double[3]);
                    final double[] xyzOP1 = OP1.getXYZ(new double[3]);
                    double dihedral = VectorMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzOP1);
                    double twoPiOver3 = 2.0 * Math.PI / 3.0;
                    double target = Crystal.modToRange(dihedral + twoPiOver3, -Math.PI, Math.PI);
                    List<Atom> otherO = bondedO.stream().
                            filter(o -> o != OP1).
                            sorted(Comparator.comparingDouble((Atom o) -> {
                                double[] xyzO = o.getXYZ(new double[3]);
                                double dihedO = VectorMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzO);
                                double diff = dihedO - target;
                                double twoPi = 2 * Math.PI;
                                diff = Crystal.modToRange(diff, 0, twoPi);
                                diff = diff < Math.PI ? diff : twoPi - diff;
                                return diff;
                            })).collect(Collectors.toList());
                    for (int i = 0; i < otherO.size(); i++) {
                        otherO.get(i).setName(String.format("OP%d", i + 2));
                    }
                } else {
                    Atom nextO3s = bondedO.stream().
                            filter(o -> !residue.getAtomList().contains(o)).
                            findAny().get();

                    // OP1 is approximately +120 degrees from next O3', OP2 is -120 degrees.
                    final double[] xyzC5s = C5s.getXYZ(new double[3]);
                    final double[] xyzO5s = O5s.getXYZ(new double[3]);
                    final double[] xyzP = P.getXYZ(new double[3]);
                    final double[] xyzNextO3s = nextO3s.getXYZ(new double[3]);
                    double dihedral = VectorMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzNextO3s);
                    double twoPiOver3 = 2.0 * Math.PI / 3.0;
                    double target = Crystal.modToRange(dihedral + twoPiOver3, -Math.PI, Math.PI);
                    List<Atom> otherO = bondedO.stream().
                            filter(o -> o != nextO3s).
                            sorted(Comparator.comparingDouble((Atom o) -> {
                                double[] xyzO = o.getXYZ(new double[3]);
                                double dihedO = VectorMath.dihedralAngle(xyzC5s, xyzO5s, xyzP, xyzO);
                                double diff = dihedO - target;
                                double twoPi = 2 * Math.PI;
                                diff = Crystal.modToRange(diff, 0, twoPi);
                                diff = diff < Math.PI ? diff : twoPi - diff;
                                return diff;
                            })).collect(Collectors.toList());
                    for (int i = 0; i < otherO.size(); i++) {
                        otherO.get(i).setName(String.format("OP%d", i + 1));
                    }
                }

                for (Atom op : bondedO) {
                    if (hasAttachedAtom(op, 1)) {
                        findBondedAtoms(op, 1).get(0).setName("H" + op.getName());
                    }
                }
            }
            renameCommonNucleobase(residue, N19, C1s, na3);
        } else {
            logger.warning(" Could not find O4\' for residue " + residue);
        }
    }

    /**
     * Renames atoms common to all standard pyrimidines (C, T, U)
     *
     * @param N1  The N1 atom.
     * @param C1s The C1' atom.
     * @return A Map containing Atoms important to finding and naming base-unique atoms.
     */
    private static Map<String, Atom> renameCommonPyrimidine(Atom N1, Atom C1s) {
        Map<String, Atom> keyAtoms = new HashMap<>();
        N1.setName("N1");
        for (Atom c : findBondedAtoms(N1, C1s, 6)) {
            if (hasAttachedAtom(c, 8)) {
                Atom C2 = c;
                C2.setName("C2");
                findBondedAtoms(C2, 8).get(0).setName("O2");
                Atom N3 = findBondedAtoms(C2, N1, 7).get(0);
                N3.setName("N3");
                keyAtoms.put("N3", N3);
                Atom C4 = findBondedAtoms(N3, C2, 6).get(0);
                C4.setName("C4");
                keyAtoms.put("C4", C4);
                Atom C5 = findBondedAtoms(C4, 6).get(0);
                C5.setName("C5");
                keyAtoms.put("C5", C5);
                if (hasAttachedAtom(C5, 1)) {
                    findBondedAtoms(C5, 1).get(0).setName("H5");
                }
            } else {
                Atom C6 = c;
                C6.setName("C6");
                findBondedAtoms(C6, 1).get(0).setName("H6");
            }
        }
        /**
         * Common atoms:
         * N1, C2, O2, N3, C4, C5, C6, H6
         */
        return keyAtoms;
    }

    /**
     * Renames atoms common to all standard purines (A, G)
     *
     * @param N9  The N9 atom.
     * @param C1s The C1' atom.
     * @return A Map containing Atoms important to finding and naming base-unique atoms.
     */
    private static Map<String, Atom> renameCommonPurine(Atom N9, Atom C1s) {
        Map<String, Atom> keyAtoms = new HashMap<>(10);
        N9.setName("N9");
        for (Atom c : findBondedAtoms(N9, C1s, 6)) {
            if (hasAttachedAtom(c, 1)) {
                Atom C8 = c;
                C8.setName("C8");
                findBondedAtoms(C8, 1).get(0).setName("H8");
                Atom N7 = findBondedAtoms(C8, N9, 7).get(0);
                N7.setName("N7");
                Atom C5 = findBondedAtoms(N7, C8, 6).get(0);
                C5.setName("C5");
            } else {
                Atom C4 = c;
                C4.setName("C4");
                Atom N3 = findBondedAtoms(C4, N9, 7).get(0);
                N3.setName("N3");
                Atom C2 = findBondedAtoms(N3, C4, 6).get(0);
                C2.setName("C2");
                keyAtoms.put("C2", C2);
                Atom N1 = findBondedAtoms(C2, N3, 7).get(0);
                N1.setName("N1"); // And not, say, "largest non-nuclear explosion ever".
                keyAtoms.put("N1", N1);
                Atom C6 = findBondedAtoms(N1, C2, 6).get(0);
                C6.setName("C6");
                keyAtoms.put("C6", C6);
            }
        }
        /**
         * Common atoms:
         * N1, C2, N3, C4, C5, C6, N7, C8, H8, N9.
         */
        return keyAtoms;
    }

    /**
     * Renames the atoms of the common nucleobases (A, C, G, T, U, and deoxy variants).
     *
     * @param residue Nucleic acid to fix atom names of.
     * @param N19     N1 of pyrimidines, N9 of purines.
     * @param C1s     C1' of the ribose sugar.
     * @param na3     Identity of the nucleic acid.
     */
    private static void renameCommonNucleobase(Residue residue, Atom N19, Atom C1s, NucleicAcid3 na3) {
        switch (na3) {
            case ADE:
            case DAD: {
                Map<String, Atom> purineBase = renameCommonPurine(N19, C1s);
                // Unique to A: H2, N6, H6[12]
                findBondedAtoms(purineBase.get("C2"), 1).get(0).setName("H2");
                Atom C6 = purineBase.get("C6");
                Atom N1 = purineBase.get("N1");
                Atom N6 = findBondedAtoms(C6, N1, 7).get(0);
                N6.setName("N6");
                List<Atom> allH6List = findBondedAtoms(N6, 1);
                Atom[] allH6 = sortAtomsByDistance(N1, allH6List);
                allH6[0].setName("H61");
                allH6[1].setName("H62");
            }
            break;
            case CYT:
            case DCY: {
                Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
                // Unique to C: N4, H4[12]
                Atom C4 = pyrimidineBase.get("C4");
                Atom N3 = pyrimidineBase.get("N3");
                Atom N4 = findBondedAtoms(C4, N3, 7).get(0);
                N4.setName("N4");
                Atom[] allH4 = sortAtomsByDistance(N3, findBondedAtoms(N4, 1));
                allH4[0].setName("H41");
                allH4[1].setName("H42");
            }
            break;
            case GUA:
            case DGU: {
                Map<String, Atom> purineBase = renameCommonPurine(N19, C1s);
                // Unique to G: H1, N2, H2[12], O6
                Atom N1 = purineBase.get("N1");
                Atom C2 = purineBase.get("C2");
                Atom C6 = purineBase.get("C6");
                findBondedAtoms(N1, 1).get(0).setName("H1");
                Atom N2 = findBondedAtoms(C2, N1, 7).
                        stream().
                        filter(n -> hasAttachedAtom(n, 1)).
                        findAny().get();
                N2.setName("N2");
                Atom[] allH2 = sortAtomsByDistance(N1, findBondedAtoms(N2, 1));
                allH2[0].setName("H21");
                allH2[1].setName("H22");
                findBondedAtoms(C6, 8).get(0).setName("O6");
            }
            break;
            case URI: {
                Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
                // Unique to U: H3, O4
                findBondedAtoms(pyrimidineBase.get("N3"), 1).get(0).setName("H3");
                findBondedAtoms(pyrimidineBase.get("C4"), 8).get(0).setName("O4");
            }
            break;
            case THY:
            case DTY: {
                Map<String, Atom> pyrimidineBase = renameCommonPyrimidine(N19, C1s);
                // Unique to T: H3, O4, C7
                findBondedAtoms(pyrimidineBase.get("N3"), 1).get(0).setName("H3");
                findBondedAtoms(pyrimidineBase.get("C4"), 8).get(0).setName("O4");
                Atom C5 = pyrimidineBase.get("C5");
                for (Atom c : findBondedAtoms(C5, 6)) {
                    List<Atom> bondedH = findBondedAtoms(c, 1);
                    if (bondedH != null && bondedH.size() == 3) {
                        c.setName("C7");
                        for (int i = 0; i < 3; i++) {
                            bondedH.get(i).setName(String.format("H7%d", i + 1));
                        }
                        break;
                    }
                }
            }
            break;
        }
    }

    /**
     * Find the O4' of a nucleic acid Residue. This is fairly unique in standard nucleotides,
     * as O4' is the only ether oxygen (bonded to two carbons).
     *
     * @param residue Residue to find O4' of.
     * @return O4'.
     */
    private static Optional<Atom> findNucleotideO4s(Residue residue) {
        assert residue.getResidueType() == Residue.ResidueType.NA;
        return findAtomsOfElement(residue, 8).
                stream().
                filter(o -> findBondedAtoms(o, 6).size() == 2).
                findAny();
    }

    public static class Mutation {

        /**
         * Residue ID of the residue to mutate.
         */
        final int resID;
        /**
         * Residue name after mutation.
         */
        final String resName;
        /**
         * Character for the chain ID of the residue that will be mutated.
         */
        final char chainChar;

        public Mutation(int resID, char chainChar, String newResName) {
            newResName = newResName.toUpperCase();
            if (newResName == null || newResName.length() != 3) {
                logger.log(Level.WARNING, String.format("Invalid mutation target: %s.", newResName));
            }
            try {
                AminoAcid3.valueOf(newResName);
            } catch (IllegalArgumentException ex) {
                logger.log(Level.WARNING, String.format("Invalid mutation target: %s.", newResName));
            }
            this.resID = resID;
            this.chainChar = chainChar;
            this.resName = newResName;
        }

        public Mutation(char chain, int res, String newName) {
            this(res, chain, newName);
        }
    }

    /**
     * PDB records that are recognized.
     */
    private enum Record {

        ANISOU,
        ATOM,
        CONECT,
        CRYST1,
        DBREF,
        END,
        MODEL,
        ENDMDL,
        HELIX,
        HETATM,
        LINK,
        MTRIX1,
        MTRIX2,
        MTRIX3,
        MODRES,
        SEQRES,
        SHEET,
        SSBOND,
        REMARK
    }

    /**
     * Presently, VERSION3_3 is default, and VERSION3_2 is anything
     * non-standard.
     */
    public enum PDBFileStandard {

        VERSION2_3,
        VERSION3_0,
        VERSION3_1,
        VERSION3_2,
        VERSION3_3
    }

    /**
     * Common HETATOM labels for water and ions.
     */
    public enum HetAtoms {

        BR,
        CA,
        CA2,
        CL,
        K,
        MG,
        MG2,
        NA,
        HOH,
        H2O,
        WAT,
        ZN,
        ZN2
    }

}
