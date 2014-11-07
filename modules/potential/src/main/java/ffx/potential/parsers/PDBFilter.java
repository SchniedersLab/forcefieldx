/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.numerics.VectorMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
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
import ffx.potential.bonded.Residue.ResiduePosition;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.utilities.Hybrid36;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
import static ffx.potential.bonded.AminoAcidUtils.buildAIB;
import static ffx.potential.bonded.AminoAcidUtils.buildAlanine;
import static ffx.potential.bonded.AminoAcidUtils.buildArginine;
import static ffx.potential.bonded.AminoAcidUtils.buildAsparagine;
import static ffx.potential.bonded.AminoAcidUtils.buildAspartate;
import static ffx.potential.bonded.AminoAcidUtils.buildCysteine;
import static ffx.potential.bonded.AminoAcidUtils.buildCystine;
import static ffx.potential.bonded.AminoAcidUtils.buildDeprotonatedCysteine;
import static ffx.potential.bonded.AminoAcidUtils.buildDeprotonatedLysine;
import static ffx.potential.bonded.AminoAcidUtils.buildDeprotonatedTyrosine;
import static ffx.potential.bonded.AminoAcidUtils.buildGlutamate;
import static ffx.potential.bonded.AminoAcidUtils.buildGlutamine;
import static ffx.potential.bonded.AminoAcidUtils.buildHistidine;
import static ffx.potential.bonded.AminoAcidUtils.buildIsoleucine;
import static ffx.potential.bonded.AminoAcidUtils.buildLeucine;
import static ffx.potential.bonded.AminoAcidUtils.buildLysine;
import static ffx.potential.bonded.AminoAcidUtils.buildMethionine;
import static ffx.potential.bonded.AminoAcidUtils.buildNeutralAsparticAcid;
import static ffx.potential.bonded.AminoAcidUtils.buildNeutralGlutamicAcid;
import static ffx.potential.bonded.AminoAcidUtils.buildNeutralHistidineD;
import static ffx.potential.bonded.AminoAcidUtils.buildNeutralHistidineE;
import static ffx.potential.bonded.AminoAcidUtils.buildOrnithine;
import static ffx.potential.bonded.AminoAcidUtils.buildPCA;
import static ffx.potential.bonded.AminoAcidUtils.buildPhenylalanine;
import static ffx.potential.bonded.AminoAcidUtils.buildProline;
import static ffx.potential.bonded.AminoAcidUtils.buildSerine;
import static ffx.potential.bonded.AminoAcidUtils.buildThreonine;
import static ffx.potential.bonded.AminoAcidUtils.buildTryptophan;
import static ffx.potential.bonded.AminoAcidUtils.buildTyrosine;
import static ffx.potential.bonded.AminoAcidUtils.buildValine;
import static ffx.potential.bonded.AminoAcidUtils.cType;
import static ffx.potential.bonded.AminoAcidUtils.caType;
import static ffx.potential.bonded.AminoAcidUtils.cbType;
import static ffx.potential.bonded.AminoAcidUtils.haType;
import static ffx.potential.bonded.AminoAcidUtils.hnType;
import static ffx.potential.bonded.AminoAcidUtils.nType;
import static ffx.potential.bonded.AminoAcidUtils.oType;
import static ffx.potential.bonded.AminoAcidUtils.removeH1_H2_H3;
import static ffx.potential.bonded.AminoAcidUtils.removeOXT_OT2;
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
import static ffx.potential.bonded.NucleicAcidUtils.c1Typ;
import static ffx.potential.bonded.NucleicAcidUtils.c2Typ;
import static ffx.potential.bonded.NucleicAcidUtils.c3Typ;
import static ffx.potential.bonded.NucleicAcidUtils.c4Typ;
import static ffx.potential.bonded.NucleicAcidUtils.c5Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h1Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h21Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h22Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h3Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h3tTyp;
import static ffx.potential.bonded.NucleicAcidUtils.h4Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h51Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h52Typ;
import static ffx.potential.bonded.NucleicAcidUtils.h5tTyp;
import static ffx.potential.bonded.NucleicAcidUtils.o2Typ;
import static ffx.potential.bonded.NucleicAcidUtils.o3Typ;
import static ffx.potential.bonded.NucleicAcidUtils.o4Typ;
import static ffx.potential.bonded.NucleicAcidUtils.o5Typ;
import static ffx.potential.bonded.NucleicAcidUtils.opTyp;
import static ffx.potential.bonded.NucleicAcidUtils.pTyp;
import static ffx.potential.bonded.Residue.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidHeavyAtoms;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidList;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcid;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcidNumber;
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
 * @since 1.0
 *
 */
public final class PDBFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(PDBFilter.class.getName());

    /**
     * PDB records that are recognized.
     */
    private enum Record {

        ANISOU, ATOM, CONECT, CRYST1, END, HELIX, HETATM, LINK, SHEET,
        SSBOND, REMARK
    };
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
     * systems, chain names are sometimes reused due to limitations in the PBD
     * format.
     *
     * However, we define segIDs to always be unique. For the first A-Z,0-9
     * series chainID == segID. Then, for second A-Z,0-9 series, the segID =
     * 1A-1Z,10-19, and for the third series segID = 2A-2Z,20-29, and so on.
     */
    private final List<String> segIDs = new ArrayList<>();
    private Character currentChainID = null;
    private String currentSegID = null;
    private boolean mutate = false;
    private int mutateResID = 0;
    private String mutateToResname = null;
    private Character mutateChainID = null;
    private boolean print = true;
    private PDBFileStandard fileStandard = VERSION3_3; // Assume current standard.
    /**
     * If true, output is directed into arrayOutput instead of the file.
     */
    private boolean listMode = false;
    private ArrayList<String> listOutput = new ArrayList<>();

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

        return newSegID;
    }
    /**
     * Keep track of ATOM record serial numbers to match them with ANISOU
     * records.
     */
    private final HashMap<Integer, Atom> atoms = new HashMap<>();

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
     * Specify the alternate location.
     *
     * @param molecularAssembly The MolecularAssembly to populate.
     * @param altLoc The alternate location to use.
     */
    public void setAltID(MolecularAssembly molecularAssembly, Character altLoc) {
        setMolecularSystem(molecularAssembly);
        currentAltLoc = altLoc;
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
                /**
                 * Read the first line of the file.
                 */
                String line = br.readLine();
                /**
                 * Parse until END is found or to the end of the file.
                 */
                while (line != null) {
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
                        case END:
                            /**
                             * Setting "line" to null will exit the loop.
                             */
                            line = null;
                            continue;
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
                            Character chainID;
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
                                } catch (NumberFormatException e) {
                                    // Use default values.
                                    if (print) {
                                        logger.warning(" No values for occupancy or b-factors; defaulting to 1.00 (further warnings suppressed).");
                                        print = false;
                                    } else if (logger.isLoggable(Level.FINE)) {
                                        logger.fine(" No values for occupancy or b-factors; defaulting to 1.00.");
                                    }
                                }
                                newAtom = new Atom(0, name, altLoc, d, resName, resSeq,
                                        chainID, occupancy, tempFactor, segID);

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
                            d = new double[3];
                            d[0] = new Double(line.substring(30, 38).trim());
                            d[1] = new Double(line.substring(38, 46).trim());
                            d[2] = new Double(line.substring(46, 54).trim());
                            occupancy = 1.0;
                            tempFactor = 1.0;
                            try {
                                occupancy = new Double(line.substring(54, 60).trim());
                                tempFactor = new Double(line.substring(60, 66).trim());
                            } catch (NumberFormatException e) {
                                // Use default values.
                                if (print) {
                                    logger.warning(" No values for occupancy or b-factors; defaulting to 1.00 (further warnings suppressed).");
                                    print = false;
                                } else if (logger.isLoggable(Level.FINE)) {
                                    logger.fine(" No values for occupancy or b-factors; defaulting to 1.00.");
                                }
                            }
                            newAtom = new Atom(0, name, altLoc, d, resName, resSeq, chainID,
                                    occupancy, tempFactor, segID);
                            newAtom.setHetero(true);
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
                            } else {
                                if (a1 == currentAltLoc && a2 == currentAltLoc) {
                                    links.add(line);
                                }
                            }
                            break;
                        case SSBOND:
                            /**
                             * SSBOND records may be invalid if chain IDs are
                             * reused.
                             *
                             * They are applied to the A conformer and not
                             * alternate conformers.
                             */
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
        List<Bond> ssBondList = new ArrayList<>();
        for (String ssbond : ssbonds) {
            try {
                Polymer c1 = activeMolecularAssembly.getChain(ssbond.substring(15, 16));
                Polymer c2 = activeMolecularAssembly.getChain(ssbond.substring(29, 30));
                Residue r1 = c1.getResidue(Hybrid36.decode(4, ssbond.substring(17, 21)));
                Residue r2 = c2.getResidue(Hybrid36.decode(4, ssbond.substring(31, 35)));
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
                    SG1.getName();
                }
                if (SG2 == null) {
                    SG2.getName();
                }
                double d = VectorMath.dist(SG1.getXYZ(), SG2.getXYZ());
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

        int pdbAtoms = activeMolecularAssembly.getAtomArray().length;
        assignAtomTypes();

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
            double d = VectorMath.dist(a1.getXYZ(), a2.getXYZ());
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

        /**
         * Finally, re-number the atoms if missing atoms were created.
         */
        if (pdbAtoms != activeMolecularAssembly.getAtomArray().length) {
            numberAtoms();
        }

        return true;
    }

    /**
     * <p>
     * numberAtoms</p>
     */
    public void numberAtoms() {
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
    public void assignAtomTypes() {
        /**
         * Create a new List to store bonds determined based on PDB atom names.
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
                if (!residues.isEmpty()) {
                    //renameNTerminusHydrogens(residues.get(0)); Not safe to use until it distinguishes between true N-termini and N-terminal residues in general.
                }
                for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
                    Residue residue = residues.get(residueNumber);
                    String name = residue.getName().toUpperCase();
                    boolean aa = false;
                    for (AminoAcid3 amino : aminoAcidList) {
                        if (amino.toString().equalsIgnoreCase(name)) {
                            aa = true;
                            renameNonstandardHydrogens(residue);
                            break;
                        }
                    }
                    // Check for a patch.
                    if (!aa) {
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
                        assignNucleicAcidAtomTypes(residues);
                    } catch (MissingHeavyAtomException missingHeavyAtomException) {
                        logger.severe(missingHeavyAtomException.toString());
                    } catch (MissingAtomTypeException missingAtomTypeException) {
                        logger.severe(missingAtomTypeException.toString());
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
                                    double b[] = ib.getXYZ();
                                    double c[] = ib.getXYZ();
                                    double d[] = ib.getXYZ();
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

                /**
                 * Compute the distance between the previous carbonyl carbon and
                 * the current nitrogen.
                 */
                double r = VectorMath.dist(pC.getXYZ(), N.getXYZ());
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
     * Assign atom types for a nucleic acid polymer.
     *
     * @param residues
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     */
    private void assignNucleicAcidAtomTypes(List<Residue> residues)
            throws MissingHeavyAtomException, MissingAtomTypeException {
        /**
         * A reference to the O3* atom of the previous base.
         */
        Atom pO3s = null;
        /**
         * Loop over residues.
         */
        int numberOfResidues = residues.size();
        for (int residueNumber = 0; residueNumber
                < numberOfResidues; residueNumber++) {
            /**
             * Match the residue name to a known nucleic acid residue.
             */
            Residue residue = residues.get(residueNumber);
            String residueName = residue.getName().toUpperCase();
            NucleicAcid3 nucleicAcid = null;
            int naNumber = -1;
            for (NucleicAcid3 nucleic : nucleicAcidList) {
                naNumber++;
                String nuc3 = nucleic.toString();
                nuc3 = nuc3.substring(nuc3.length() - 3);
                if (nuc3.equalsIgnoreCase(residueName)) {
                    nucleicAcid = nucleic;
                    break;
                }
            }
            /**
             * Do atom name conversions.
             */
            List<Atom> resAtoms = residue.getAtomList();
            int natoms = resAtoms.size();
            for (int i = 0; i < natoms; i++) {
                Atom atom = resAtoms.get(i);
                String name = atom.getName();
                name = name.replace('*', '\'');
                //name = name.replace('D', 'H');
                atom.setName(name);
            }

            /**
             * Check if this is a 3' phosphate being listed as its own residue.
             */
            /*if (residue.getAtomList().size() == 1) {
             Atom P3s = (Atom) residue.getAtomNode("P");
             if (P3s != null) {
             Residue prevResidue = residue.getPreviousResidue();
             if (prevResidue != null) {
             Atom O2sPrev = (Atom) prevResidue.getAtomNode("O2\'");
             if (O2sPrev == null) {
             P3s = buildHeavy(prevResidue, "P3s", null, 1247);
             } else {
             P3s = buildHeavy(prevResidue, "P3s", null, 1235);
             }
             } else {
             return;
             }
             } else {
             return;
             }
             }*/
            /**
             * Check if the sugar is deoxyribose and change the residue name if
             * necessary.
             */
            boolean isDNA = false;
            Atom O2s = (Atom) residue.getAtomNode("O2\'");
            if (O2s == null) {
                /**
                 * Assume deoxyribose (DNA) since there is an O2* atom.
                 */
                isDNA = true;
                if (!residueName.startsWith("D")) {
                    switch (nucleicAcid) {
                        case ADE:
                            nucleicAcid = NucleicAcid3.DAD;
                            residueName = "DAD";
                            residue.setName(residueName);
                            break;
                        case CYT:
                            nucleicAcid = NucleicAcid3.DCY;
                            residueName = "DCY";
                            residue.setName(residueName);
                            break;
                        case GUA:
                            nucleicAcid = NucleicAcid3.DGU;
                            residueName = "DGU";
                            residue.setName(residueName);
                            break;
                        case THY:
                            nucleicAcid = NucleicAcid3.DTY;
                            residueName = "DTY";
                            residue.setName(residueName);
                            break;
                        default:
                    }
                }
            } else {
                /**
                 * Assume ribose (RNA) since there is an O2* atom.
                 */
                if (residueName.startsWith("D")) {
                    switch (nucleicAcid) {
                        case DAD:
                            nucleicAcid = NucleicAcid3.ADE;
                            residueName = "ADE";
                            residue.setName(residueName);
                            break;
                        case DCY:
                            nucleicAcid = NucleicAcid3.CYT;
                            residueName = "CYT";
                            residue.setName(residueName);
                            break;
                        case DGU:
                            nucleicAcid = NucleicAcid3.GUA;
                            residueName = "GUA";
                            residue.setName(residueName);
                            break;
                        case DTY:
                            nucleicAcid = NucleicAcid3.THY;
                            residueName = "THY";
                            residue.setName(residueName);
                            break;
                        default:
                    }
                }
            }

            /**
             * Set a position flag.
             */
            ResiduePosition position = MIDDLE_RESIDUE;
            if (residueNumber == 0) {
                position = FIRST_RESIDUE;
            } else if (residueNumber == numberOfResidues - 1) {
                position = LAST_RESIDUE;
            }
            /**
             * Build the phosphate atoms of the current residue.
             */
            Atom P = null;
            Atom O5s = null;
            if (position == FIRST_RESIDUE) {
                /**
                 * The 5' O5' oxygen of the nucleic acid is generally terminated
                 * by 1.) A phosphate group PO3 (-3). 2.) A hydrogen.
                 *
                 * If the base has phosphate atom we will assume a PO3 group.
                 */
                P = (Atom) residue.getAtomNode("P");
                if (P != null) {
                    if (isDNA) {
                        P = buildHeavy(residue, "P", null, 1247);
                        buildHeavy(residue, "OP1", P, 1248);
                        buildHeavy(residue, "OP2", P, 1248);
                        buildHeavy(residue, "OP3", P, 1248);
                        O5s = buildHeavy(residue, "O5\'", P, 1246);
                    } else {
                        P = buildHeavy(residue, "P", null, 1235);
                        buildHeavy(residue, "OP1", P, 1236);
                        buildHeavy(residue, "OP2", P, 1236);
                        buildHeavy(residue, "OP3", P, 1236);
                        O5s = buildHeavy(residue, "O5\'", P, 1234);
                    }
                } else {
                    if (isDNA) {
                        O5s = buildHeavy(residue, "O5\'", P, 1244);
                    } else {
                        O5s = buildHeavy(residue, "O5\'", P, 1232);
                    }
                }
            } else {
                P = buildHeavy(residue, "P", pO3s, pTyp[naNumber]);
                buildHeavy(residue, "OP1", P, opTyp[naNumber]);
                buildHeavy(residue, "OP2", P, opTyp[naNumber]);
                O5s = buildHeavy(residue, "O5\'", P, o5Typ[naNumber]);
            }
            /**
             * Build the ribose sugar atoms of the current base.
             */
            Atom C5s = buildHeavy(residue, "C5\'", O5s, c5Typ[naNumber]);
            Atom C4s = buildHeavy(residue, "C4\'", C5s, c4Typ[naNumber]);
            Atom O4s = buildHeavy(residue, "O4\'", C4s, o4Typ[naNumber]);
            Atom C1s = buildHeavy(residue, "C1\'", O4s, c1Typ[naNumber]);
            Atom C3s = buildHeavy(residue, "C3\'", C4s, c3Typ[naNumber]);
            Atom C2s = buildHeavy(residue, "C2\'", C3s, c2Typ[naNumber]);
            buildBond(C2s, C1s);
            Atom O3s = null;
            if (position == LAST_RESIDUE || numberOfResidues == 1) {
                if (isDNA) {
                    O3s = buildHeavy(residue, "O3\'", C3s, 1249);
                } else {
                    O3s = buildHeavy(residue, "O3\'", C3s, 1237);
                }
            } else {
                O3s = buildHeavy(residue, "O3\'", C3s, o3Typ[naNumber]);
            }
            if (!isDNA) {
                O2s = buildHeavy(residue, "O2\'", C2s, o2Typ[naNumber]);
            }
            /**
             * Build the backbone hydrogen atoms.
             */
            if (position == FIRST_RESIDUE && P == null) {
                buildHydrogen(residue, "H5T", O5s, 1.00e0, C5s, 109.5e0, C4s, 180.0e0, 0, h5tTyp[naNumber]);
            }
            buildHydrogen(residue, "H5\'1", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, 1, h51Typ[naNumber]);
            buildHydrogen(residue, "H5\'2", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, -1, h52Typ[naNumber]);
            buildHydrogen(residue, "H4\'", C4s, 1.09e0, C5s, 109.5e0, C3s, 109.5e0, -1, h4Typ[naNumber]);
            buildHydrogen(residue, "H3\'", C3s, 1.09e0, C4s, 109.5e0, C2s, 109.5e0, -1, h3Typ[naNumber]);
            if (isDNA) {
                buildHydrogen(residue, "H2\'1", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber]);
                buildHydrogen(residue, "H2\'2", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, 1, h22Typ[naNumber]);
            } else {
                buildHydrogen(residue, "H2\'", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber]);
                // Add the O2' Methyl for OMC and OMG
                if (nucleicAcid == NucleicAcid3.OMC || nucleicAcid == NucleicAcid3.OMG) {
                    Atom CM2 = buildHeavy(residue, "CM2", O2s, 1427);
                    Atom HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, O2s, 109.5e0, C2s, 0.0e0, 0, 1428);
                    buildHydrogen(residue, "HM22", CM2, 1.08e0, O2s, 109.5e0, HM21, 109.5e0, 1, 1429);
                    buildHydrogen(residue, "HM23", CM2, 1.08e0, O2s, 109.5e0, HM21, 109.5e0, -1, 1430);
                } else {
                    buildHydrogen(residue, "HO\'", O2s, 1.00e0, C2s, 109.5e0, C3s, 180.0e0, 0, h22Typ[naNumber]);
                }
            }
            buildHydrogen(residue, "H1\'", C1s, 1.09e0, O4s, 109.5e0, C2s, 109.5e0, -1, h1Typ[naNumber]);
            if (position == LAST_RESIDUE || numberOfResidues == 1) {
                buildHydrogen(residue, "H3T", O3s, 1.00e0, C3s, 109.5e0, C4s, 180.0e0, 0, h3tTyp[naNumber]);
                // Else, if it is terminated by a 3' phosphate cap:
                // Will need to see how PDB would label a 3' phosphate cap.
            }
            /**
             * Build the nucleic acid base.
             */
            try {
                assignNucleicAcidBaseAtomTypes(nucleicAcid, residue, C1s, O4s, C2s);
            } catch (MissingHeavyAtomException missingHeavyAtomException) {
                logger.throwing(PDBFilter.class.getName(), "assignNucleicAcidAtomTypes", missingHeavyAtomException);
                throw missingHeavyAtomException;
            }

            /**
             * Do some checks on the current base to make sure all atoms have
             * been assigned an atom type.
             */
            resAtoms = residue.getAtomList();
            for (Atom atom : resAtoms) {
                AtomType atomType = atom.getAtomType();
                if (atomType == null) {
                    MissingAtomTypeException missingAtomTypeException = new MissingAtomTypeException(residue, atom);
                    logger.throwing(PDBFilter.class.getName(), "assignNucleicAcidAtomTypes", missingAtomTypeException);
                    throw missingAtomTypeException;
                }
                int numberOfBonds = atom.getNumBonds();
                if (numberOfBonds != atomType.valence) {
                    if (atom == O3s && numberOfBonds == atomType.valence - 1
                            && position != LAST_RESIDUE && numberOfResidues != 1) {
                        continue;
                    }
                    logger.log(Level.WARNING, format(" An atom for residue %s has the wrong number of bonds:\n %s",
                            residueName, atom.toString()));
                    logger.log(Level.WARNING, format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
                }
            }

            /**
             * Save a reference to the current O3* oxygen.
             */
            pO3s = O3s;
        }
    }

    /**
     * Assign atom types to the nucleic acid base.
     *
     * @param nucleicAcid The nucleic acid base to use.
     * @param residue The residue node.
     * @param C1s The CS* attachement atom.
     *
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     *
     * @since 1.0
     */
    private void assignNucleicAcidBaseAtomTypes(NucleicAcid3 nucleicAcid, Residue residue, Atom C1s,
            Atom O4s, Atom C2s)
            throws MissingHeavyAtomException {
        double glyco = 0;
        switch (nucleicAcid) {
            case ADE:
                Atom N9,
                 C8,
                 N7,
                 C5,
                 C6,
                 N6,
                 N1,
                 C2,
                 N3,
                 C4;
                N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1017);
                C8 = buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180, 0, 1021);
                N7 = buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1020);
                C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1019);
                C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1025);
                N6 = buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1027);
                N1 = buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1024);
                C2 = buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1023);
                N3 = buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1022);
                C4 = buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1018);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1030);
                buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1028);
                buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1029);
                buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1026);
                break;
            case M1MA:
                Atom CM1;
                Atom HM11;
                N9 = buildHeavy(residue, "N9", C1s, 1605);
                C8 = buildHeavy(residue, "C8", N9, 1609);
                N7 = buildHeavy(residue, "N7", C8, 1608);
                C5 = buildHeavy(residue, "C5", N7, 1607);
                C6 = buildHeavy(residue, "C6", C5, 1613);
                N6 = buildHeavy(residue, "N6", C6, 1615);
                N1 = buildHeavy(residue, "N1", C6, 1612);
                C2 = buildHeavy(residue, "C2", N1, 1611);
                N3 = buildHeavy(residue, "N3", C2, 1610);
                C4 = buildHeavy(residue, "C4", N3, 1606);
                CM1 = buildHeavy(residue, "CM1", N1, 1619);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1614);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 109.5e0, C4, 180.0e0, 0, 1623);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1618);
                buildHydrogen(residue, "HN61", N6, 1.00e0, C6, 109.5e0, C5, 0.0e0, 0, 1616);
                buildHydrogen(residue, "HN62", N6, 1.00e0, C6, 109.5e0, C5, 109.5e0, 0, 1617);
                HM11 = buildHydrogen(residue, "HM11", CM1, 1.08e0, N1, 109.5e0, C2, 0.0e0, 0, 1620);
                buildHydrogen(residue, "HM12", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, 1, 1621);
                buildHydrogen(residue, "HM13", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, -1, 1622);
                break;
            case CYT:
            case OMC:
                Atom O2;
                Atom N4;
                N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1078);
                C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco + 180, 0, 1079);
                O2 = buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1084);
                N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180.0, 0, 1080);
                C4 = buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1081);
                N4 = buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, O2, 180.0, 0, 1085);
                C5 = buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1082);
                C6 = buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1083);
                buildBond(C6, N1);
                buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1086);
                buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1087);
                buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1088);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1089);
                break;
            case M5MC:
                Atom CM5;
                Atom HM51;
                N1 = buildHeavy(residue, "N1", C1s, 1508);
                C2 = buildHeavy(residue, "C2", N1, 1509);
                O2 = buildHeavy(residue, "O2", C2, 1514);
                N3 = buildHeavy(residue, "N3", C2, 1510);
                C4 = buildHeavy(residue, "C4", N3, 1511);
                N4 = buildHeavy(residue, "N4", C4, 1515);
                C5 = buildHeavy(residue, "C5", C4, 1512);
                C6 = buildHeavy(residue, "C6", C5, 1513);
                CM5 = buildHeavy(residue, "CM5", C5, 1519);
                buildBond(C6, N1);
                buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1516);
                buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, C5, 0.0e0, 0, 1517);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1518);
                HM51 = buildHydrogen(residue, "HM51", CM5, 1.08e0, C5, 109.5e0, C4, 0.0e0, 0, 1520);
                buildHydrogen(residue, "HM52", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, 1, 1521);
                buildHydrogen(residue, "HM53", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, -1, 1522);
                break;
            case GUA:
            case OMG:
                Atom O6;
                Atom N2;
                N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1047);
                C8 = buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1051);
                N7 = buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1050);
                C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1049);
                C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1055);
                O6 = buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1060);
                N1 = buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1054);
                C2 = buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1053);
                N2 = buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1057);
                N3 = buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1052);
                C4 = buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1048);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1061);
                buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1056);
                buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1058);
                buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1059);
                break;
            case YYG:
                Atom C10,
                 C11,
                 C12,
                 C3,
                 C13,
                 C14,
                 C15,
                 C16,
                 O17,
                 O18,
                 C19,
                 N20,
                 C21,
                 O22,
                 O23,
                 C24,
                 H31,
                 H101,
                 H191,
                 H241;
                N9 = buildHeavy(residue, "N9", C1s, 1640);
                C8 = buildHeavy(residue, "C8", N9, 1644);
                N7 = buildHeavy(residue, "N7", C8, 1643);
                C5 = buildHeavy(residue, "C5", N7, 1642);
                C6 = buildHeavy(residue, "C6", C5, 1648);
                O6 = buildHeavy(residue, "O6", C6, 1650);
                N1 = buildHeavy(residue, "N1", C6, 1647);
                C2 = buildHeavy(residue, "C2", N1, 1646);
                N2 = buildHeavy(residue, "N2", C2, 1649);
                N3 = buildHeavy(residue, "N3", C2, 1645);
                C3 = buildHeavy(residue, "C3", N3, 1652);
                C4 = buildHeavy(residue, "C4", N3, 1641);
                C11 = buildHeavy(residue, "C11", N2, 1657);
                C10 = buildHeavy(residue, "C10", C11, 1658);
                C12 = buildHeavy(residue, "C12", C11, 1656);
                C13 = buildHeavy(residue, "C13", C12, 1662);
                C14 = buildHeavy(residue, "C14", C13, 1665);
                C15 = buildHeavy(residue, "C15", C14, 1668);
                C16 = buildHeavy(residue, "C16", C15, 1675);
                O17 = buildHeavy(residue, "O17", C16, 1676);
                O18 = buildHeavy(residue, "O18", C16, 1674);
                C19 = buildHeavy(residue, "C19", O18, 1670);
                N20 = buildHeavy(residue, "N20", C15, 1677);
                C21 = buildHeavy(residue, "C21", N20, 1679);
                O22 = buildHeavy(residue, "O22", C21, 1680);
                O23 = buildHeavy(residue, "O23", C21, 1681);
                C24 = buildHeavy(residue, "C24", O23, 1682);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildBond(N1, C12);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1651);
                H31 = buildHydrogen(residue, "H31", C3, 1.08e0, N3, 109.5e0, C4, 0.0e0, 0, 1653);
                buildHydrogen(residue, "H32", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, 1, 1654);
                buildHydrogen(residue, "H33", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, -1, 1655);
                H101 = buildHydrogen(residue, "H101", C10, 1.08e0, C11, 109.5e0, N2, 0.0e0, 0, 1659);
                buildHydrogen(residue, "H102", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, 1, 1660);
                buildHydrogen(residue, "H103", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, -1, 1661);
                buildHydrogen(residue, "H131", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, 1, 1663);
                buildHydrogen(residue, "H132", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, -1, 1664);
                buildHydrogen(residue, "H141", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, 1, 1666);
                buildHydrogen(residue, "H142", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, -1, 1667);
                buildHydrogen(residue, "H15", C15, 1.08e0, C14, 109.5e0, O18, 180.e0, 0, 1669);
                H191 = buildHydrogen(residue, "H191", C19, 1.08e0, O18, 109.5e0, C16, 0.0e0, 0, 1671);
                buildHydrogen(residue, "H192", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, 1, 1672);
                buildHydrogen(residue, "H193", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, -1, 1673);
                buildHydrogen(residue, "HN2", N20, 1.00e0, C15, 109.5e0, O22, 180.0e0, 0, 1678);
                H241 = buildHydrogen(residue, "H241", C24, 1.08e0, O23, 109.5e0, C21, 0.0e0, 0, 1683);
                buildHydrogen(residue, "H242", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, 1, 1684);
                buildHydrogen(residue, "H243", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, -1, 1685);
                break;
            case M2MG:
                Atom CM2;
                Atom HM21;
                N9 = buildHeavy(residue, "N9", C1s, 1316);
                C8 = buildHeavy(residue, "C8", N9, 1320);
                N7 = buildHeavy(residue, "N7", C8, 1319);
                C5 = buildHeavy(residue, "C5", N7, 1318);
                C6 = buildHeavy(residue, "C6", C5, 1324);
                O6 = buildHeavy(residue, "O6", C6, 1328);
                N1 = buildHeavy(residue, "N1", C6, 1323);
                C2 = buildHeavy(residue, "C2", N1, 1322);
                N2 = buildHeavy(residue, "N2", C2, 1326);
                N3 = buildHeavy(residue, "N3", C2, 1321);
                C4 = buildHeavy(residue, "C4", N3, 1317);
                CM2 = buildHeavy(residue, "CM2", N2, 1330);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1329);
                buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1325);
                buildHydrogen(residue, "H2", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1327);
                HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1331);
                buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1332);
                buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1333);
                break;
            case M2G:
                N9 = buildHeavy(residue, "N9", C1s, 1379);
                C8 = buildHeavy(residue, "C8", N9, 1383);
                N7 = buildHeavy(residue, "N7", C8, 1382);
                C5 = buildHeavy(residue, "C5", N7, 1381);
                C6 = buildHeavy(residue, "C6", C5, 1387);
                O6 = buildHeavy(residue, "O6", C6, 1390);
                N1 = buildHeavy(residue, "N1", C6, 1386);
                C2 = buildHeavy(residue, "C2", N1, 1385);
                N2 = buildHeavy(residue, "N2", C2, 1389);
                N3 = buildHeavy(residue, "N3", C2, 1384);
                C4 = buildHeavy(residue, "C4", N3, 1380);
                CM1 = buildHeavy(residue, "CM1", N2, 1392);
                CM2 = buildHeavy(residue, "CM2", N2, 1396);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1391);
                buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1388);
                HM11 = buildHydrogen(residue, "HM11", CM1, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1393);
                buildHydrogen(residue, "HM12", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, 1, 1394);
                buildHydrogen(residue, "HM13", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, -1, 1395);
                HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1397);
                buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1398);
                buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1399);
                break;
            case M7MG:
                Atom CM7;
                Atom HM71;
                N9 = buildHeavy(residue, "N9", C1s, 1539);
                C8 = buildHeavy(residue, "C8", N9, 1543);
                N7 = buildHeavy(residue, "N7", C8, 1542);
                C5 = buildHeavy(residue, "C5", N7, 1541);
                C6 = buildHeavy(residue, "C6", C5, 1547);
                O6 = buildHeavy(residue, "O6", C6, 1552);
                N1 = buildHeavy(residue, "N1", C6, 1546);
                C2 = buildHeavy(residue, "C2", N1, 1545);
                N2 = buildHeavy(residue, "N2", C2, 1549);
                N3 = buildHeavy(residue, "N3", C2, 1544);
                C4 = buildHeavy(residue, "C4", N3, 1540);
                CM7 = buildHeavy(residue, "CM7", N7, 1555);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H81", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, 1, 1553);
                buildHydrogen(residue, "H82", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, -1, 1554);
                buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1548);
                buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1550);
                buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N3, 0.0e0, 0, 1551);
                HM71 = buildHydrogen(residue, "HM71", CM7, 1.08e0, N7, 109.5e0, C8, 0.0e0, 0, 1556);
                buildHydrogen(residue, "HM72", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, 1, 1557);
                buildHydrogen(residue, "HM73", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, -1, 1558);
                break;
            case URI:
                Atom O4;
                N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1106);
                C2 = buildHeavy(residue, "C2", N1, 1.38, C1s, 117.1, O4s, glyco, 0, 1107);
                O2 = buildHeavy(residue, "O2", C2, 1.22, N1, 123.2, C1s, 0.0, 0, 1112);
                N3 = buildHeavy(residue, "N3", C2, 1.37, N1, 114.8, C1s, 180.0, 0, 1108);
                C4 = buildHeavy(residue, "C4", N3, 1.38, C2, 127.0, N1, 0.0, 0, 1109);
                O4 = buildHeavy(residue, "O4", C4, 1.23, N3, 119.8, C2, 180.0, 0, 1114);
                C5 = buildHeavy(residue, "C5", C4, 1.44, N3, 114.7, C2, 0.0, 0, 1110);
                C6 = buildHeavy(residue, "C6", C5, 1.34, O4, 119.2, C4, 0.0, 0, 1111);
                buildBond(C6, N1);
                buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1113);
                buildHydrogen(residue, "H5", C5, 1.08e0, C4, 120.4e0, N3, 180.0e0, 0, 1115);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1116);
                break;
            case PSU:
                // C1s bonds to C5 in PsuedoUridine
                C5 = buildHeavy(residue, "C5", C1s, 1485);
                C6 = buildHeavy(residue, "C6", C5, 1486);
                N1 = buildHeavy(residue, "N1", C6, 1481);
                C2 = buildHeavy(residue, "C2", N1, 1482);
                O2 = buildHeavy(residue, "O2", C2, 1487);
                N3 = buildHeavy(residue, "N3", C2, 1483);
                C4 = buildHeavy(residue, "C4", N3, 1484);
                O4 = buildHeavy(residue, "O4", C4, 1489);
                buildBond(C4, C5);
                buildHydrogen(residue, "H1", N1, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1491);
                buildHydrogen(residue, "H3", N3, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1488);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 120.0e0, C1s, 0.0e0, 0, 1490);
                break;
            case H2U:
                N1 = buildHeavy(residue, "N1", C1s, 1350);
                C2 = buildHeavy(residue, "C2", N1, 1351);
                O2 = buildHeavy(residue, "O2", C2, 1356);
                N3 = buildHeavy(residue, "N3", C2, 1352);
                C4 = buildHeavy(residue, "C4", N3, 1353);
                O4 = buildHeavy(residue, "O4", C4, 1358);
                C5 = buildHeavy(residue, "C5", C4, 1354);
                C6 = buildHeavy(residue, "C6", C5, 1355);
                buildBond(C6, N1);
                buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1357);
                buildHydrogen(residue, "H51", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, 1, 1359);
                buildHydrogen(residue, "H52", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, -1, 1360);
                buildHydrogen(residue, "H61", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, 1, 1361);
                buildHydrogen(residue, "H62", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, -1, 1362);
                break;
            case M5MU:
                Atom C5M;
                Atom H5M1;
                N1 = buildHeavy(residue, "N1", C1s, 1575);
                C2 = buildHeavy(residue, "C2", N1, 1576);
                O2 = buildHeavy(residue, "O2", C2, 1581);
                N3 = buildHeavy(residue, "N3", C2, 1577);
                C4 = buildHeavy(residue, "C4", N3, 1578);
                O4 = buildHeavy(residue, "O4", C4, 1583);
                C5 = buildHeavy(residue, "C5", C4, 1579);
                C6 = buildHeavy(residue, "C6", C5, 1580);
                C5M = buildHeavy(residue, "C5M", C5, 1585);
                buildBond(C6, N1);
                buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1582);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1584);
                H5M1 = buildHydrogen(residue, "H5M1", C5M, 1.08e0, C5, 109.5e0, C6, 0.0e0, 0, 1586);
                buildHydrogen(residue, "H5M2", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, 1, 1587);
                buildHydrogen(residue, "H5M3", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, -1, 1588);
                break;
            case DAD:
                N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1132);
                C8 = buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180, 0, 1136);
                N7 = buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1135);
                C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1134);
                C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1140);
                N6 = buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1142);
                N1 = buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1139);
                C2 = buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1138);
                N3 = buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1137);
                C4 = buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1133);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1145);
                buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1143);
                buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1144);
                buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1141);
                break;
            case DCY:
                N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1191);
                C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco, 0, 1192);
                O2 = buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1197);
                N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180, 0, 1193);
                C4 = buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1194);
                N4 = buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, C2, 180.0, 0, 1198);
                C5 = buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1195);
                C6 = buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1196);
                buildBond(C6, N1);
                buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1199);
                buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1200);
                buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1201);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1202);
                break;
            case DGU:
                N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1161);
                C8 = buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1165);
                N7 = buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1164);
                C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1163);
                C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1169);
                O6 = buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1174);
                N1 = buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1168);
                C2 = buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1167);
                N2 = buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1171);
                N3 = buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1166);
                C4 = buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1162);
                buildBond(C4, C5);
                buildBond(C4, N9);
                buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1175);
                buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1170);
                buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1172);
                buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1173);
                break;
            case DTY:
                Atom C7;
                Atom H;
                N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1218);
                C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.1, O4s, glyco, 0, 1219);
                O2 = buildHeavy(residue, "O2", C2, 1.22, N1, 122.9, C1s, 0.0, 0, 1224);
                N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 115.4, C1s, 180.0, 0, 1220);
                C4 = buildHeavy(residue, "C4", N3, 1.38, C2, 126.4, N1, 0.0, 0, 1221);
                O4 = buildHeavy(residue, "O4", C4, 1.23, N3, 120.5, C2, 180.0, 0, 1226);
                C5 = buildHeavy(residue, "C5", C4, 1.44, N3, 114.1, C2, 0.0, 0, 1222);
                C7 = buildHeavy(residue, "C7", C5, 1.50, C4, 117.5, N3, 180.0, 0, 1227);
                C6 = buildHeavy(residue, "C6", C5, 1.34, C4, 120.8, N3, 0.0, 0, 1223);
                buildBond(C6, N1);
                buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.8e0, N1, 180.0e0, 0, 1225);
                H = buildHydrogen(residue, "H71", C7, 1.09e0, C5, 109.5e0, C4, 0.0e0, 0, 1228);
                buildHydrogen(residue, "H72", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, 1, 1228);
                buildHydrogen(residue, "H73", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, -1, 1228);
                buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1229);
                break;
        }
    }

    /**
     * Check for missing heavy atoms. This check ignores special terminating
     * groups like FOR, NH2, etc.
     *
     * @param aminoAcidNumber
     * @param aminoAcid
     * @param position
     * @param residue
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     */
    private void checkForMissingHeavyAtoms(int aminoAcidNumber, AminoAcid3 aminoAcid,
            ResiduePosition position, Residue residue) throws MissingHeavyAtomException {
        int expected = aminoAcidHeavyAtoms[aminoAcidNumber];
        if (aminoAcid != AminoAcid3.GLY && expected >= 4) {
            int actual = 0;
            List<Atom> resAtoms = residue.getAtomList();
            for (Atom atom : resAtoms) {
                String label = atom.getName().toUpperCase();
                if (!(label.equalsIgnoreCase("OXT") || label.equalsIgnoreCase("OT2"))) {
                    if (!label.startsWith("H") && !label.startsWith("D")) {
                        actual++;
                    }
                }
            }
            if (actual != expected) {
                Atom N = (Atom) residue.getAtomNode("N");
                if (N == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("N", null, null);
                    throw e;
                }
                Atom CA = (Atom) residue.getAtomNode("CA");
                if (CA == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("CA", null, null);
                    throw e;
                }
                Atom C = (Atom) residue.getAtomNode("C");
                if (C == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("C", null, null);
                    throw e;
                }
                Atom O = (Atom) residue.getAtomNode("O");
                if (O == null && position == LAST_RESIDUE) {
                    O = (Atom) residue.getAtomNode("OT1");
                }
                if (O == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("O", null, null);
                    throw e;
                }
            }
        }
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
            assignAminoAcidAtomTypes(residue, previousResidue, nextResidue);
        }
    }

    private void assignAminoAcidAtomTypes(Residue residue, Residue previousResidue, Residue nextResidue)
            throws MissingHeavyAtomException, MissingAtomTypeException {

        String residueName = residue.getName().toUpperCase();

        int j = 1;
        ResiduePosition position = MIDDLE_RESIDUE;
        if (previousResidue == null) {
            j = 0;
            position = FIRST_RESIDUE;
        } else if (nextResidue == null) {
            j = 2;
            position = LAST_RESIDUE;
            /**
             * If the last residue only contains a nitrogen turn it into an NH2
             * group.
             */
            Atom N = (Atom) residue.getAtomNode("N");
            if (residue.getAtomNodeList().size() == 1 && N != null) {
                residueName = "NH2".intern();
                residue.setName(residueName);
            }
        }

        AminoAcid3 aminoAcid = getAminoAcid(residueName);
        int aminoAcidNumber = getAminoAcidNumber(residueName);
        /**
         * Non-standard Amino Acid; use ALA backbone types.
         */
        boolean nonStandard = false;
        if (aminoAcid == AminoAcid3.UNK) {
            aminoAcidNumber = getAminoAcidNumber("ALA");
            nonStandard = true;
        }

        /**
         * Only the last residue in a chain should have an OXT/OT2 atom.
         */
        if (nextResidue != null) {
            removeOXT_OT2(residue);
        }

        /**
         * Only the first nitrogen should have H1, H2 and H3 atoms, unless it's
         * an NME cap.
         */
        if (previousResidue != null) {
            removeH1_H2_H3(aminoAcid, residue);
        }

        /**
         * Check for missing heavy atoms. This check ignores special terminating
         * groups like FOR, NH2, etc.
         */
        if (!nonStandard) {
            try {
                checkForMissingHeavyAtoms(aminoAcidNumber, aminoAcid, position, residue);
            } catch (MissingHeavyAtomException e) {
                logger.log(Level.INFO, " {0} could not be parsed.", residue.toString());
                throw e;
            }
        }

        Atom pC = null;
        Atom pCA = null;
        if (previousResidue != null) {
            pC = (Atom) previousResidue.getAtomNode("C");
            pCA = (Atom) previousResidue.getAtomNode("CA");
        }

        /**
         * Backbone heavy atoms.
         */
        Atom N = (Atom) residue.getAtomNode("N");
        if (N != null) {
            N.setAtomType(findAtomType(nType[j][aminoAcidNumber]));
            if (position != FIRST_RESIDUE) {
                buildBond(pC, N);
            }
        }

        Atom CA = null;
        Atom C = null;
        Atom O = null;
        if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NH2)) {
            if (aminoAcid == AminoAcid3.ACE || aminoAcid == AminoAcid3.NME) {
                CA = buildHeavy(residue, "CH3", N, caType[j][aminoAcidNumber]);
            } else {
                CA = buildHeavy(residue, "CA", N, caType[j][aminoAcidNumber]);
            }
            if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NME)) {
                C = buildHeavy(residue, "C", CA, cType[j][aminoAcidNumber]);
                O = (Atom) residue.getAtomNode("O");
                if (O == null) {
                    O = (Atom) residue.getAtomNode("OT1");
                }
                AtomType atomType = findAtomType(oType[j][aminoAcidNumber]);
                if (O == null) {
                    MissingHeavyAtomException missingHeavyAtom = new MissingHeavyAtomException("O", atomType, C);
                    throw missingHeavyAtom;
                }
                O.setAtomType(atomType);
                buildBond(C, O);
            }
        }
        /**
         * Nitrogen hydrogen atoms.
         */
        AtomType atomType = findAtomType(hnType[j][aminoAcidNumber]);
        switch (position) {
            case FIRST_RESIDUE:
                switch (aminoAcid) {
                    case PRO:
                        buildHydrogenAtom(residue, "H2", N, 1.02, CA, 109.5, C, 0.0, 0, atomType);
                        buildHydrogenAtom(residue, "H3", N, 1.02, CA, 109.5, C, -120.0, 0, atomType);
                        break;
                    case PCA:
                        buildHydrogenAtom(residue, "H", N, 1.02, CA, 109.5, C, -60.0, 0, atomType);
                        break;
                    case ACE:
                        break;
                    default:
                        buildHydrogenAtom(residue, "H1", N, 1.02, CA, 109.5, C, 180.0, 0, atomType);
                        buildHydrogenAtom(residue, "H2", N, 1.02, CA, 109.5, C, 60.0, 0, atomType);
                        buildHydrogenAtom(residue, "H3", N, 1.02, CA, 109.5, C, -60.0, 0, atomType);
                }
                break;
            case LAST_RESIDUE:
                switch (aminoAcid) {
                    case NH2:
                        buildHydrogenAtom(residue, "H1", N, 1.02, pC, 119.0, pCA, 0.0, 0, atomType);
                        buildHydrogenAtom(residue, "H2", N, 1.02, pC, 119.0, pCA, 180.0, 0, atomType);
                        break;
                    case NME:
                        buildHydrogenAtom(residue, "H", N, 1.02, pC, 118.0, CA, 121.0, 1, atomType);
                        break;
                    default:
                        buildHydrogenAtom(residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1, atomType);
                }
                break;
            default:
                // Mid-chain nitrogen hydrogen.
                buildHydrogenAtom(residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1, atomType);
        }
        /**
         * C-alpha hydrogen atoms.
         */
        String haName = "HA";
        if (aminoAcid == AminoAcid3.GLY) {
            haName = "HA2";
        }
        atomType = findAtomType(haType[j][aminoAcidNumber]);
        switch (position) {
            case FIRST_RESIDUE:
                switch (aminoAcid) {
                    case FOR:
                        buildHydrogenAtom(residue, "H", C, 1.12, O, 0.0, null, 0.0, 0, atomType);
                        break;
                    case ACE:
                        buildHydrogenAtom(residue, "H1", CA, 1.10, C, 109.5, O, 180.0, 0, atomType);
                        buildHydrogenAtom(residue, "H2", CA, 1.10, C, 109.5, O, 60.0, 0, atomType);
                        buildHydrogenAtom(residue, "H3", CA, 1.10, C, 109.5, O, -60.0, 0, atomType);
                        break;
                    default:
                        buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1, atomType);
                        break;
                }
                break;
            case LAST_RESIDUE:
                switch (aminoAcid) {
                    case NME:
                        buildHydrogenAtom(residue, "H1", CA, 1.10, N, 109.5, pC, 180.0, 0, atomType);
                        buildHydrogenAtom(residue, "H2", CA, 1.10, N, 109.5, pC, 60.0, 0, atomType);
                        buildHydrogenAtom(residue, "H3", CA, 1.10, N, 109.5, pC, -60.0, 0, atomType);
                        break;
                    default:
                        buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1, atomType);
                }
                break;
            default:
                buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.0, -1, atomType);
        }
        /**
         * Build the amino acid side chain.
         */
        assignAminoAcidSideChain(position, aminoAcid, residue, CA, N, C);

        /**
         * Build the terminal oxygen if the residue is not NH2 or NME.
         */
        if (position == LAST_RESIDUE && !(aminoAcid == AminoAcid3.NH2 || aminoAcid == AminoAcid3.NME)) {
            atomType = findAtomType(oType[2][aminoAcidNumber]);
            Atom OXT = (Atom) residue.getAtomNode("OXT");
            if (OXT == null) {
                OXT = (Atom) residue.getAtomNode("OT2");
                if (OXT != null) {
                    OXT.setName("OXT");
                }
            }
            if (OXT == null) {
                String resName = C.getResidueName();
                int resSeq = C.getResidueNumber();
                Character chainID = C.getChainID();
                Character altLoc = C.getAltLoc();
                String segID = C.getSegID();
                double occupancy = C.getOccupancy();
                double tempFactor = C.getTempFactor();
                OXT = new Atom(0, "OXT", altLoc, new double[3], resName, resSeq, chainID,
                        occupancy, tempFactor, segID);
                OXT.setAtomType(atomType);
                residue.addMSNode(OXT);
                intxyz(OXT, C, 1.25, CA, 117.0, O, 126.0, 1);
            } else {
                OXT.setAtomType(atomType);
            }
            buildBond(C, OXT);
        }
        /**
         * Do some checks on the current residue to make sure all atoms have
         * been assigned an atom type.
         */
        List<Atom> resAtoms = residue.getAtomList();
        for (Atom atom : resAtoms) {
            atomType = atom.getAtomType();
            if (atomType == null) {
                MissingAtomTypeException missingAtomTypeException = new MissingAtomTypeException(residue, atom);
                throw missingAtomTypeException;
            }
            int numberOfBonds = atom.getNumBonds();
            if (numberOfBonds != atomType.valence) {
                if (atom == C && numberOfBonds == atomType.valence - 1 && position != LAST_RESIDUE) {
                    continue;
                }
                logger.warning(format(" An atom for residue %s has the wrong number of bonds:\n %s",
                        residueName, atom.toString()));
                logger.warning(format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
            }
        }
    }

    /**
     * Assign atom types to a single amino acid side chain.
     *
     * @param position The position of this amino acid in the chain.
     * @param aminoAcid The amino acid to use.
     * @param residue The residue node.
     * @param CA The C-alpha carbon of this residue.
     * @param N The peptide nitrogen of this residue.
     * @param C The peptide carbonyl carbon.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException this
     * exception is thrown if when heavy is atom is missing that cannot be
     * built.
     */
    public void assignAminoAcidSideChain(ResiduePosition position, AminoAcid3 aminoAcid, Residue residue,
            Atom CA, Atom N, Atom C) throws MissingHeavyAtomException {
        int k = cbType[aminoAcid.ordinal()];
        switch (aminoAcid) {
            case GLY:
                switch (position) {
                    case FIRST_RESIDUE:
                        k = haType[0][k];
                        break;
                    case LAST_RESIDUE:
                        k = haType[2][k];
                        break;
                    default:
                        k = haType[1][k];

                }
                buildHydrogen(residue, "HA3", CA, 1.10, N, 109.5, C, 109.5, 1, k);
                break;
            case ALA:
                buildAlanine(residue, CA, N, C, k, forceField, bondList);
                break;
            case VAL:
                buildValine(residue, CA, N, C, k, forceField, bondList);
                break;
            case LEU:
                buildLeucine(residue, CA, N, C, k, forceField, bondList);
                break;
            case ILE:
                buildIsoleucine(residue, CA, N, C, k, forceField, bondList);
                break;
            case SER:
                buildSerine(residue, CA, N, C, k, forceField, bondList);
                break;
            case THR:
                buildThreonine(residue, CA, N, C, k, forceField, bondList);
                break;
            case CYS:
                buildCysteine(residue, CA, N, C, k, forceField, bondList);
                break;
            case CYX:
                buildCystine(residue, CA, N, C, k, forceField, bondList);
                break;
            case CYD:
                buildDeprotonatedCysteine(residue, CA, N, C, k, forceField, bondList);
                break;
            case PRO:
                buildProline(residue, CA, N, C, k, position, forceField, bondList);
                break;
            case PHE:
                buildPhenylalanine(residue, CA, N, C, k, forceField, bondList);
                break;
            case TYR:
                buildTyrosine(residue, CA, N, C, k, forceField, bondList);
                break;
            case TYD:
                buildDeprotonatedTyrosine(residue, CA, N, C, k, forceField, bondList);
                break;
            case TRP:
                buildTryptophan(residue, CA, N, C, k, forceField, bondList);
                break;
            case HIS:
                buildHistidine(residue, CA, N, C, k, forceField, bondList);
                break;
            case HID:
                buildNeutralHistidineD(residue, CA, N, C, k, forceField, bondList);
                break;
            case HIE:
                buildNeutralHistidineE(residue, CA, N, C, k, forceField, bondList);
                break;
            case ASP:
                buildAspartate(residue, CA, N, C, k, forceField, bondList);
                break;
            case ASH:
                buildNeutralAsparticAcid(residue, CA, N, C, k, forceField, bondList);
                break;
            case ASN:
                buildAsparagine(residue, CA, N, C, k, forceField, bondList);
                break;
            case GLU:
                buildGlutamate(residue, CA, N, C, k, forceField, bondList);
                break;
            case GLH:
                buildNeutralGlutamicAcid(residue, CA, N, C, k, forceField, bondList);
                break;
            case GLN:
                buildGlutamine(residue, CA, N, C, k, forceField, bondList);
                break;
            case MET:
                buildMethionine(residue, CA, N, C, k, forceField, bondList);
                break;
            case LYS:
                buildLysine(residue, CA, N, C, k, forceField, bondList);
                break;
            case LYD:
                buildDeprotonatedLysine(residue, CA, N, C, k, forceField, bondList);
                break;
            case ARG:
                buildArginine(residue, CA, N, C, k, forceField, bondList);
                break;
            case ORN:
                buildOrnithine(residue, CA, N, C, k, forceField, bondList);
                break;
            case AIB:
                buildAIB(residue, CA, N, C, k, forceField, bondList);
                break;
            case PCA:
                buildPCA(residue, CA, N, C, k, forceField, bondList);
                break;
            case UNK:
                String residueName = residue.getName();
                logger.log(Level.INFO, " Patching side-chain {0}", residueName);
                HashMap<String, AtomType> types = forceField.getAtomTypes(residueName);
                HashMap<AtomType, AtomType> typeMap = new HashMap<>();
                if (!types.isEmpty()) {
                    boolean patched = true;
                    ArrayList<Atom> residueAtoms = residue.getAtomList();
                    // Assign atom types for side-chain atoms.
                    for (Atom atom : residueAtoms) {
                        String atomName = atom.getName().toUpperCase();
                        AtomType internalType = atom.getAtomType();
                        // Map the internal type to the new type.
                        if (internalType != null) {
                            AtomType newType = types.get(atomName);
                            if (newType != null) {
                                typeMap.put(newType, internalType);
                                types.remove(atomName);
                            }
                        } else {
                            AtomType atomType = types.get(atomName);
                            if (atomType == null) {
                                logger.log(Level.INFO, " No atom type was found for {0} of {1}.",
                                        new Object[]{atomName, residueName});
                                patched = false;
                                break;
                            } else {
                                atom.setAtomType(atomType);
                                types.remove(atomName);
                            }
                        }
                    }

                    forceField.patchClassesAndTypes(typeMap);

                    // Create a new multipole type for HA with the correct frame.
                    Atom HA = (Atom) residue.getAtomNode("HA");
                    Atom CAlpha = (Atom) residue.getAtomNode("CA");
                    Atom CBeta = (Atom) residue.getAtomNode("CB");
                    int frame[] = new int[3];
                    frame[0] = HA.getAtomType().type;
                    frame[1] = CAlpha.getAtomType().type;
                    frame[2] = forceField.getAtomType("Alanine", "CB").type;
                    MultipoleType multipoleType = forceField.getMultipoleType(frame[0] + " " + frame[1] + " " + frame[2]);

                    frame[2] = CBeta.getAtomType().type;
                    multipoleType = new MultipoleType(multipoleType.charge, multipoleType.dipole,
                            multipoleType.quadrupole, frame, multipoleType.frameDefinition);
                    forceField.addForceFieldType(multipoleType);

                    // Check for missing heavy atoms.
                    if (patched && !types.isEmpty()) {
                        for (AtomType type : types.values()) {
                            if (type.atomicNumber != 1) {
                                logger.log(Level.INFO, " Missing heavy atom {0}", type.name);
                                patched = false;
                                break;
                            }
                        }
                    }
                    // Create bonds between known atoms.
                    if (patched) {
                        for (Atom atom : residueAtoms) {
                            String atomName = atom.getName();
                            String bonds[] = forceField.getBonds(residueName, atomName);
                            if (bonds != null) {
                                for (String name : bonds) {
                                    Atom atom2 = (Atom) residue.getAtomNode(name);
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
                        HashMap<String, Atom> atomMap = new HashMap<>();
                        for (Atom atom : residueAtoms) {
                            atomMap.put(atom.getName().toUpperCase(), atom);
                        }
                        for (String atomName : types.keySet()) {
                            AtomType type = types.get(atomName);
                            String bonds[] = forceField.getBonds(residueName, atomName.toUpperCase());
                            if (bonds == null || bonds.length != 1) {
                                patched = false;
                                logger.log(Level.INFO, " Check biotype for hydrogen {0}.", type.name);
                                break;
                            }
                            // Get the heavy atom the hydrogen is bonded to.
                            Atom ia = atomMap.get(bonds[0].toUpperCase());
                            Atom hydrogen = new Atom(0, atomName, ia.getAltLoc(), new double[3],
                                    ia.getResidueName(), ia.getResidueNumber(), ia.getChainID(),
                                    ia.getOccupancy(), ia.getTempFactor(), ia.getSegID());
                            logger.log(Level.FINE, " Created hydrogen {0}.", atomName);
                            hydrogen.setAtomType(type);
                            hydrogen.setHetero(true);
                            residue.addMSNode(hydrogen);
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
                             * Building the hydrogens depends on hybridization
                             * and the locations of other bonded atoms.
                             */
                            logger.log(Level.FINE, " Bonding {0} to {1} ({2} of {3}).",
                                    new Object[]{atomName, ia.getName(), numBonds, valence});
                            switch (valence) {
                                case 4:
                                    switch (numBonds) {
                                        case 3:
                                            // Find the average coordinates of atoms ib, ic and id.
                                            double b[] = ib.getXYZ();
                                            double c[] = ib.getXYZ();
                                            double d[] = ib.getXYZ();
                                            double a[] = new double[3];
                                            a[0] = (b[0] + c[0] + d[0]) / 3.0;
                                            a[1] = (b[1] + c[1] + d[1]) / 3.0;
                                            a[2] = (b[2] + c[2] + d[2]) / 3.0;

                                            // Place the hydrogen at chiral position #1.
                                            intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 0);
                                            double e1[] = hydrogen.getXYZ();
                                            double ret[] = new double[3];
                                            diff(a, e1, ret);
                                            double l1 = r(ret);

                                            // Place the hydrogen at chiral position #2.
                                            intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 1);
                                            double e2[] = hydrogen.getXYZ();
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
                                            logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                                            patched = false;
                                    }
                                    break;
                                case 3:
                                    switch (numBonds) {
                                        case 2:
                                            intxyz(hydrogen, ia, 1.0, ib, 120.0, ic, 180.0, 0);
                                            break;
                                        case 1:
                                            intxyz(hydrogen, ia, 1.0, ib, 120.0, null, 0.0, 0);
                                            break;
                                        case 0:
                                            intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                            break;
                                        default:
                                            logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
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
                                            logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                                            patched = false;
                                    }
                                    break;
                                case 1:
                                    switch (numBonds) {
                                        case 0:
                                            intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                                            break;
                                        default:
                                            logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
                                            patched = false;
                                    }
                                    break;
                                default:
                                    logger.log(Level.INFO, " Check biotype for hydrogen {0}.", atomName);
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
                        logger.log(Level.SEVERE, format(" Could not patch %s.", residueName));
                    } else {
                        logger.log(Level.INFO, " Patch for {0} succeeded.", residueName);
                    }
                } else {

                    switch (position) {
                        case FIRST_RESIDUE:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 355);
                            break;
                        case LAST_RESIDUE:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 506);
                            break;
                        default:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 6);
                    }
                }
                break;
        }
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
        FileWriter fw;
        BufferedWriter bw;
        try {
            File newFile = saveFile;
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            if (!listMode) {
                fw = new FileWriter(newFile, false);
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
     * {@inheritDoc}
     *
     * Write out the Atomic information in PDB format.
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        if (saveFile == null) {
            return false;
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
            logger.log(Level.INFO, " Saving {0}", newFile.getName());
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
                        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                        sb.replace(22, 26, String.format("%4s", Hybrid36.encode(4, resID)));
                        // Loop over atoms
                        ArrayList<Atom> residueAtoms = residue.getAtomList();
                        boolean altLocFound = false;
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
                                residueAtoms = altResidue.getAtomList();
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
            logger.log(Level.INFO, " Saving {0}", newFile.getName());
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
                            if (!entries[0].equals(entries[0].replaceAll("\\D+",""))) {
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
    public void writeAtom(Atom atom, int serial, StringBuilder sb,
            StringBuilder anisouSB, BufferedWriter bw)
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
        double xyz[] = atom.getXYZ();
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
        double[] anisou = atom.getAnisou();
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

    public void writeSIFTAtom(Atom atom, int serial, StringBuilder sb,
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
        double xyz[] = atom.getXYZ();
        sb.replace(6, 16, String.format("%5s " + padLeft(name.toUpperCase(), 4), Hybrid36.encode(5, serial)));
        Character altLoc = atom.getAltLoc();
        if (altLoc != null) {
            sb.setCharAt(16, altLoc);
        } else {
            sb.setCharAt(16, ' ');
        }
        if (siftScore == null) {
            sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                xyz[0], xyz[1], xyz[2], atom.getOccupancy(), 2.00));
        } else {
            sb.replace(30, 66, String.format("%8.3f%8.3f%8.3f%6.2f%6.2f",
                xyz[0], xyz[1], xyz[2], atom.getOccupancy(), Float.parseFloat(siftScore)));
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
        double[] anisou = atom.getAnisou();
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

    private enum HetAtoms {

        HOH, H2O, WAT, NA, K, MG, MG2, CA, CA2, CL, BR, ZN, ZN2
    };

    /**
     * Ensures proper naming of hydrogens according to latest PDB format.
     * Presently mostly guesses at which hydrogens to re-assign, which may cause
     * chirality errors for prochiral hydrogens. If necessary, we will implement
     * more specific mapping.
     *
     * @param residue
     */
    private void renameNonstandardHydrogens(Residue residue) {
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
            // Handles situations such as 1H where it should be H1, etc.
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
                // Pray, for I have no examples to work from.
                break;
        }
    }

    private Atom buildHeavy(MSGroup residue, String atomName, Atom bondedTo, int key)
            throws MissingHeavyAtomException {
        return BondedUtils.buildHeavy(residue, atomName, bondedTo, key, forceField, bondList);
    }

    private Atom buildHeavy(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, int lookUp) {
        return BondedUtils.buildHeavy(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, lookUp,
                forceField, bondList);
    }

    private Atom buildHydrogen(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, int lookUp) {
        return BondedUtils.buildHydrogen(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, lookUp, forceField, bondList);
    }

    private Atom buildHydrogenAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, AtomType atomType) {
        return BondedUtils.buildHydrogenAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType,
                forceField, bondList);
    }

    public Bond buildBond(Atom a1, Atom a2) {
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

    // Presently, VERSION3_3 is default, and VERSION3_2 is anything non-standard.
    public enum PDBFileStandard {

        VERSION3_3, VERSION3_2, VERSION3_1, VERSION3_0, VERSION2_3;
    }
}
