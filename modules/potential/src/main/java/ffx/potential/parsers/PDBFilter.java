/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.numerics.VectorMath;
import ffx.potential.MolecularAssembly;
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

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
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
 * @see <a href="http://www.wwpdb.org/documentation/format32/v3.2.html"> PDB
 * format 3.2</a>
 *
 * @author Michael J. Schnieders
 *
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
     *
     * The expectation is for chain naming from A-Z, then from 0-9. For large
     * systems, chain names are sometimes reused due to limitations in the PDB
     * format.
     *
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
    /**
     * Residue ID of the residue to mutate.
     */
    private int mutateResID = 0;
    /**
     * Residue name after mutation.
     */
    private String mutateToResname = null;
    /**
     * Character for the chain ID of the residue that will be mutated.
     */
    private Character mutateChainID = null;
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
     * <p>
     * Constructor for PDBFilter.</p>
     *
     * @param files a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public PDBFilter(List<File> files, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        super(files, molecularAssembly, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
    }

    /**
     * Parse the PDB File from a URL.
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public PDBFilter(File file, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssembly, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
    }

    /**
     * Parse the PDB File from a URL.
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public PDBFilter(File file, List<MolecularAssembly> molecularAssemblies,
            ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssemblies, forceField, properties);
        bondList = new ArrayList<>();
        this.fileType = FileType.PDB;
    }

    /**
     * Mutate a residue at the PDB file is being parsed.
     *
     * @param chainID the Chain ID of the residue to mutate.
     * @param resID the Residue ID of the residue to mutate.
     * @param name the 3-letter code of the amino acid to mutate to.
     */
    public void mutate(Character chainID, int resID, String name) {
        if (name != null && name.length() == 3) {
            logger.info(String.format(" Mutating chain %c residue %d to %s.", chainID, resID, name));
            mutate = true;
            mutateResID = resID;
            mutateChainID = chainID;
            mutateToResname = name;
        }
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
     * @param altLoc The alternate location to use.
     */
    public void setAltID(MolecularAssembly molecularAssembly, Character altLoc) {
        setMolecularSystem(molecularAssembly);
        currentAltLoc = altLoc;
    }

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

    public void setListMode(boolean set) {
        listMode = set;
        listOutput = new ArrayList<>();
    }

    public ArrayList<String> getListOutput() {
        return listOutput;
    }

    public void clearListOutput() {
        listOutput.clear();
    }
    
    /**
     * {@inheritDoc}
     *
     * Parse the PDB File
     *
     * @return true if the file is read successfully.
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
                        if (!chainIDs.contains(mutateChainID)) {
                            if (chainIDs.size() == 1) {
                                logger.warning(String.format(" Chain ID %c for "
                                        + "mutation not found: only one chain %c "
                                        + "found.", mutateChainID, chainIDs.get(0)));
                                mutateChainID = chainIDs.get(0);
                            } else {
                                logger.warning(String.format(" Chain ID %c for "
                                        + "mutation not found: mutation will not "
                                        + "proceed.", mutateChainID));
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
                 * Parse until END or ENDMDL is found, or to the end of the file.
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
                            adp[0] = new Integer(line.substring(28, 35).trim()) * 1.0e-4;
                            adp[1] = new Integer(line.substring(35, 42).trim()) * 1.0e-4;
                            adp[2] = new Integer(line.substring(42, 49).trim()) * 1.0e-4;
                            adp[3] = new Integer(line.substring(49, 56).trim()) * 1.0e-4;
                            adp[4] = new Integer(line.substring(56, 63).trim()) * 1.0e-4;
                            adp[5] = new Integer(line.substring(63, 70).trim()) * 1.0e-4;
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
                                if (mutate && chainID.equals(mutateChainID) && mutateResID == resSeq) {
                                    String atomName = name.toUpperCase();
                                    if (atomName.equals("N") || atomName.equals("C")
                                            || atomName.equals("O") || atomName.equals("CA")) {
                                        printAtom = true;
                                        resName = mutateToResname;
                                    } else {
                                        logger.info(String.format(" Deleting atom %s of %s %d",
                                                atomName, resName, resSeq));
                                        break;
                                    }
                                }
                                d = new double[3];
                                d[0] = new Double(line.substring(30, 38).trim());
                                d[1] = new Double(line.substring(38, 46).trim());
                                d[2] = new Double(line.substring(46, 54).trim());
                                occupancy = 1.0;
                                tempFactor = 1.0;
                                try {
                                    occupancy = new Double(line.substring(54, 60).trim());
                                    tempFactor = new Double(line.substring(60, 66).trim());
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
                                    if (newAtom.xyzIndex == 0) {
                                        newAtom.setXYZIndex(xyzIndex++);
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
                            d[0] = new Double(line.substring(30, 38).trim());
                            d[1] = new Double(line.substring(38, 46).trim());
                            d[2] = new Double(line.substring(46, 54).trim());
                            occupancy = 1.0;
                            tempFactor = 1.0;
                            try {
                                occupancy = new Double(line.substring(54, 60).trim());
                                tempFactor = new Double(line.substring(60, 66).trim());
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
                                newAtom.setXYZIndex(xyzIndex++);
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
                            double aaxis = new Double(line.substring(6, 15).trim());
                            double baxis = new Double(line.substring(15, 24).trim());
                            double caxis = new Double(line.substring(24, 33).trim());
                            double alpha = new Double(line.substring(33, 40).trim());
                            double beta = new Double(line.substring(40, 47).trim());
                            double gamma = new Double(line.substring(47, 54).trim());
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
                                logger.info(format(" Ignoring LINK record as alternate locations do not match\n %s.", line));
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

    public void setIgnoreInactiveAtoms(boolean ignoreInactiveAtoms) {
        this.ignoreUnusedAtoms = ignoreInactiveAtoms;
    }
    
    /**
     * Sets whether this PDBFilter should log each time it saves to a file.
     * @param logWrites 
     */
    public void setLogWrites(boolean logWrites) {
        this.logWrites = logWrites;
    }

    /**
     * <p>
     * writeFileWithHeader</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header a {@link java.lang.StringBuilder} object.
     * @return a boolean.
     */
    public boolean writeFileWithHeader(File saveFile, StringBuilder header) {
        return writeFileWithHeader(saveFile, header, true);
    }

    /**
     * <p>
     * writeFileWithHeader</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header a {@link java.lang.StringBuilder} object.
     * @param append a boolean.
     * @return a boolean.
     */
    public boolean writeFileWithHeader(File saveFile, StringBuilder header, boolean append) {
        FileWriter fw;
        BufferedWriter bw;
        if (header.charAt(header.length() - 1) != '\n') {
            header.append("\n");
        }
        try {
            File newFile = saveFile;
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (!listMode) {
                fw = new FileWriter(newFile, append);
                bw = new BufferedWriter(fw);
                bw.write(header.toString());
                bw.close();
            } else {
                listOutput.add(header.toString());
            }
        } catch (Exception e) {
            String message = "Exception writing to file: " + saveFile.toString();
            logger.log(Level.WARNING, message, e);
            return false;
        }
        return writeFile(saveFile, true);
    }

    /**
     * <p>
     * writeFile</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param append a {@link java.lang.StringBuilder} object.
     * @param printLinear Whether to print atoms linearly or by element
     * @return Success of writing.
     */
    public boolean writeFile(File saveFile, boolean append, boolean printLinear) {
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
                newFile = version(saveFile);
            }
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (logWrites) {
                logger.log(Level.INFO, " Saving {0}", newFile.getName());
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
                                    if (SG1.xyzIndex < SG2.xyzIndex) {
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
                        ArrayList<Atom> residueAtoms = residue.getAtomList();
                        ArrayList<Atom> backboneAtoms = residue.getBackboneAtoms();
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
            ArrayList<Molecule> molecules = activeMolecularAssembly.getMolecules();
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
                ArrayList<Atom> ionAtoms = ion.getAtomList();
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
                ArrayList<Atom> waterAtoms = water.getAtomList();
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
            if (!append) {
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
                                    if (SG1.xyzIndex < SG2.xyzIndex) {
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
            ArrayList<Molecule> molecules = activeMolecularAssembly.getMolecules();
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
     *
     * Write out the Atomic information in PDB format.
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        return writeFile(saveFile, append, false);
    }

    /**
     * Locate disulfide bonds based on SSBOND records.
     *
     * @param ssbonds
     *
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
                if (d < 3.0) {
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
     *
     * Known limitations include: 1) No building n- and c-terminal loops. 2) No
     * support for DBREF1 or DBREF2 records. 3) Incomplete optimization scheme
     * to position the loops.
     *
     * @param xyzIndex
     *
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
            a.xyzIndex = index++;
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
        List<Molecule> molecules = activeMolecularAssembly.getMolecules();
        for (Molecule n : molecules) {
            n.reOrderAtoms();
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
                            checkHydrogenAtomNames(residue);
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
                        logger.severe(missingHeavyAtomException.toString());
                    } catch (MissingAtomTypeException missingAtomTypeException) {
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
                        assignNucleicAcidAtomTypes(residues, forceField, bondList);
                    } catch (MissingHeavyAtomException | MissingAtomTypeException e) {
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
                    String message = "Error assigning atom types to a water.";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }

        // Assign small molecule atom types.
        ArrayList<Molecule> molecules = activeMolecularAssembly.getMolecules();
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
    }

    private List<List<Residue>> findChainBreaks(List<Residue> residues, double cutoff) {
        List<List<Residue>> subChains = new ArrayList<List<Residue>>();

        List<Residue> subChain = null;
        Residue previousResidue = null;
        Atom pC = null;
        StringBuilder sb = new StringBuilder(" Chain Breaks:");

        for (Residue residue : residues) {
            List<Atom> resAtoms = residue.getAtomList();
            if (pC == null) {
                /**
                 * Initialization.
                 */
                subChain = new ArrayList<Residue>();
                subChain.add(residue);
                subChains.add(subChain);
            } else {
                /**
                 * Find the Nitrogen of the current residue.
                 */
                Atom N = null;
                for (Atom a : resAtoms) {
                    if (a.getName().equalsIgnoreCase("N")) {
                        N = a;
                        break;
                    }
                }
                if (N == null) {
                    subChain.add(residue);
                    continue;
                }
                /**
                 * Compute the distance between the previous carbonyl carbon and
                 * the current nitrogen.
                 */
                double r = VectorMath.dist(pC.getXYZ(null), N.getXYZ(null));
                if (r > cutoff) {
                    /**
                     * Start a new chain.
                     */
                    subChain = new ArrayList<Residue>();
                    subChain.add(residue);
                    subChains.add(subChain);
                    sb.append(format("\n C-N distance of %6.2f A for %s and %s.",
                            r, previousResidue.toString(), residue.toString()));
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
                if (a.getName().equalsIgnoreCase("C")) {
                    pC = a;
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
     *
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     *
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
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param serial a int.
     * @param sb a {@link java.lang.StringBuilder} object.
     * @param anisouSB a {@link java.lang.StringBuilder} object.
     * @param bw a {@link java.io.BufferedWriter} object.
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
        sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                xyz[0], xyz[1], xyz[2], atom.getOccupancy(), atom.getTempFactor()));
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
    public boolean readNext() {
        return readNext(false);
    }
    
    @Override
    public boolean readNext(boolean resetPosition) {
        // ^ is beginning of line, \\s+ means "one or more whitespace", (\\d+) means match and capture one or more digits.
        Pattern modelPatt = Pattern.compile("^MODEL\\s+(\\d+)");
        modelsRead = resetPosition ? 1 : modelsRead + 1;
        boolean eof = true;
        
        for (MolecularAssembly system : systems) {
            File file = system.getFile();
            currentFile = file;
            try {
                BufferedReader currentReader;
                if (readers.containsKey(system)) {
                    currentReader = readers.get(system);
                    if (!currentReader.ready()) {
                        currentReader = new BufferedReader(new FileReader(currentFile));
                        readers.put(system, currentReader);
                    }
                } else {
                    currentReader = new BufferedReader(new FileReader(currentFile));
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
                            logger.log(Level.INFO, String.format(" Reading model %d for %s", modelNum, currentFile));
                            eof = false;
                            break;
                        }
                    }
                    line = currentReader.readLine();
                }

                if (eof) {
                    logger.log(Level.INFO, String.format(" End of file reached for %s", file));
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
                                    d[0] = new Double(line.substring(30, 38).trim());
                                    d[1] = new Double(line.substring(38, 46).trim());
                                    d[2] = new Double(line.substring(46, 54).trim());
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
     * @param residue
     */
    private void checkHydrogenAtomNames(Residue residue) {
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
        return BondedUtils.buildHeavy(residue, atomName, bondedTo, key, forceField, bondList);
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
