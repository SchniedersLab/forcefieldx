/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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

import java.io.*;
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
import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.ResidueEnumerations.NucleicAcid3;
import ffx.potential.bonded.*;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.parameters.*;
import ffx.potential.parsers.PDBFilter.ResiduePosition;
import ffx.utilities.Hybrid36;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
import static ffx.potential.ResidueEnumerations.aminoAcidHeavyAtoms;
import static ffx.potential.ResidueEnumerations.aminoAcidList;
import static ffx.potential.ResidueEnumerations.nucleicAcidList;
import static ffx.potential.parsers.INTFilter.intxyz;
import static ffx.potential.parsers.PDBFilter.ResiduePosition.*;

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
    private List<Character> altLocs = new ArrayList<Character>();
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
    private List<String> segIDs = new ArrayList<String>();
    private Character currentChainID = null;
    private String currentSegID = null;
    private boolean mutate = false;
    private int mutateResID = 0;
    private String mutateToResname = null;
    private Character mutateChainID = null;

    /**
     * Mutate a residue at the PDB file is being parsed.
     */
    public void mutate(Character chainID, int resID, String name) {
        if (name != null && name.length() == 3) {
            mutate = true;
            mutateResID = resID;
            mutateChainID = chainID;
            mutateToResname = name;
        }
    }

    /**
     * <p>clearSegIDs</p>
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
        String newSegID = null;
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
    private HashMap<Integer, Atom> atoms = new HashMap<Integer, Atom>();

    /**
     * <p>Constructor for PDBFilter.</p>
     *
     * @param files a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.bonded.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public PDBFilter(List<File> files, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        super(files, molecularAssembly, forceField, properties);
        this.fileType = FileType.PDB;
    }

    /**
     * Parse the PDB File from a URL.
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.bonded.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public PDBFilter(File file, MolecularAssembly molecularAssembly,
            ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssembly, forceField, properties);
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
        this.fileType = FileType.PDB;
    }

    /**
     * <p>pdbForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String pdbForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".pdb";
    }

    /**
     * <p>cifForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String cifForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".cif";
    }

    /**
     * Specify the alternate location.
     *
     * @param molecularAssembly The MolecularAssembly to populate.
     * @param altLoc The alternate location to use.
     */
    public void setAltID(MolecularAssembly molecularAssembly,
            Character altLoc) {
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
     */
    @Override
    public boolean readFile() {
        // First atom is #1, to match xyz file format
        int xyzIndex = 1;
        setFileRead(false);
        systems.add(activeMolecularAssembly);

        List<String> conects = new ArrayList<String>();
        List<String> links = new ArrayList<String>();
        List<String> ssbonds = new ArrayList<String>();
        List<String> structs = new ArrayList<String>();
        BufferedReader br = null;
        try {
            for (int i = 0; i < files.size(); i++) {
                currentFile = files.get(i);
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
                    Record record = null;
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
                            Integer serial = new Integer(Hybrid36.decode(5, line.substring(6, 11)));
                            Character altLoc = new Character(line.substring(16, 17).toUpperCase().charAt(0));
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
                            serial = new Integer(Hybrid36.decode(5, line.substring(6, 11)));
                            String name = line.substring(12, 16).trim();
                            altLoc = new Character(line.substring(16, 17).toUpperCase().charAt(0));
                            if (!altLocs.contains(altLoc)) {
                                altLocs.add(altLoc);
                            }
                            if (!altLoc.equals(' ') && !altLoc.equals('A')
                                    && !altLoc.equals(currentAltLoc)) {
                                break;
                            }
                            String resName = line.substring(17, 20).trim();
                            Character chainID = line.substring(21, 22).charAt(0);
                            String segID = getSegID(chainID);
                            int resSeq = new Integer(Hybrid36.decode(4, line.substring(22, 26)));
                            if (mutate && chainID.equals(mutateChainID) && mutateResID == resSeq) {
                                logger.info(" Mutate code reached");
                                String atomName = name.toUpperCase();
                                if (atomName.equals("N") || atomName.equals("C")
                                        || atomName.equals("O") || atomName.equals("CA")) {
                                    resName = mutateToResname;
                                } else {
                                    break;
                                }
                            }
                            double d[] = new double[3];
                            d[0] = new Double(line.substring(30, 38).trim());
                            d[1] = new Double(line.substring(38, 46).trim());
                            d[2] = new Double(line.substring(46, 54).trim());
                            double occupancy = new Double(line.substring(54, 60).trim());
                            double tempFactor = new Double(line.substring(60, 66).trim());
                            Atom newAtom = new Atom(0, name, altLoc, d, resName, resSeq,
                                    chainID, occupancy, tempFactor, segID);
                            Atom returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
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
                            }
                            break;
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
                            serial = new Integer(Hybrid36.decode(5, line.substring(6, 11)));
                            name = line.substring(12, 16).trim();
                            altLoc = new Character(line.substring(16, 17).toUpperCase().charAt(0));
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
                            resSeq = new Integer(Hybrid36.decode(4, line.substring(22, 26)));
                            d = new double[3];
                            d[0] = new Double(line.substring(30, 38).trim());
                            d[1] = new Double(line.substring(38, 46).trim());
                            d[2] = new Double(line.substring(46, 54).trim());
                            occupancy = new Double(line.substring(54, 60).trim());
                            tempFactor = new Double(line.substring(60, 66).trim());
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
        List<Bond> ssBondList = new ArrayList<Bond>();
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
     * <p>numberAtoms</p>
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
    public void assignAtomTypes() {
        /**
         * Create a new List to store bonds determined based on PDB atom names.
         */
        bondList = new ArrayList<Bond>();

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
                    continue;
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
                    logger.log(Level.WARNING, format(" An atom for residue %s has the wrong number of bonds:\n %s", residueName, atom.toString()));
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
                C2 = buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 167);
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
        Atom pC = null;
        Atom pCA = null;
        /**
         * Loop over amino acid residues.
         */
        int numberOfResidues = residues.size();
        for (int residueNumber = 0; residueNumber
                < numberOfResidues; residueNumber++) {
            Residue residue = residues.get(residueNumber);
            String residueName = residue.getName().toUpperCase();
            int j = 1;
            ResiduePosition position = MIDDLE_RESIDUE;
            if (residueNumber == 0) {
                j = 0;
                position = FIRST_RESIDUE;
            } else if (residueNumber == numberOfResidues - 1) {
                j = 2;
                position = LAST_RESIDUE;
                /**
                 * If the last residue only contains a nitrogen turn it into an
                 * NH2 group.
                 */
                Atom N = (Atom) residue.getAtomNode("N");
                if (residue.getAtomNodeList().size() == 1 && N != null) {
                    residueName = "NH2".intern();
                    residue.setName(residueName);
                }
            }
            /**
             * Only the last residue in a chain should have an OXT/OT2 atom.
             */
            if (position != LAST_RESIDUE && numberOfResidues > 1) {
                Atom OXT = (Atom) residue.getAtomNode("OXT");
                if (OXT != null) {
                    residue.deleteAtom(OXT);
                }
                Atom OT2 = (Atom) residue.getAtomNode("OT2");
                if (OT2 != null) {
                    residue.deleteAtom(OT2);
                }
            }

            AminoAcid3 aminoAcid = AminoAcid3.UNK;
            int aminoAcidNumber = -1;
            for (AminoAcid3 amino : aminoAcidList) {
                aminoAcidNumber++;
                if (amino.toString().equalsIgnoreCase(residueName)) {
                    aminoAcid = amino;
                    break;
                }
            }

            /**
             * Only the first nitrogen should have H1, H2 and H3 atoms, unless
             * it's an NME cap.
             */
            if (position != FIRST_RESIDUE && numberOfResidues > 1 && aminoAcid != AminoAcid3.NME) {
                Atom H1 = (Atom) residue.getAtomNode("H1");
                if (H1 != null) {
                    residue.deleteAtom(H1);
                }
                Atom H2 = (Atom) residue.getAtomNode("H2");
                if (H2 != null) {
                    residue.deleteAtom(H2);
                }
                Atom H3 = (Atom) residue.getAtomNode("H3");
                if (H3 != null) {
                    residue.deleteAtom(H3);
                }
            }

            /**
             * Non-standard Amino Acid; use ALA backbone types.
             */
            boolean nonStandard = false;
            if (aminoAcid == AminoAcid3.UNK) {
                residueName = "ALA";
                aminoAcidNumber = -1;
                nonStandard = true;
                for (AminoAcid3 amino : aminoAcidList) {
                    aminoAcidNumber++;
                    if (amino.toString().equalsIgnoreCase(residueName)) {
                        break;
                    }
                }
                residueName = residue.getName().toUpperCase();
            }

            /**
             * Check for missing heavy atoms. This check ignores special
             * terminating groups like FOR, NH2, etc.
             */
            int expected = aminoAcidHeavyAtoms[aminoAcidNumber];
            if (aminoAcid != AminoAcid3.GLY && expected >= 4 && !nonStandard) {
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
                    /**
                     * if (aminoAcid == AminoAcid3.ALA && actual == 4) {
                     * residueName = "GLY".intern();
                     * residue.setName(residueName); } else if (actual == 5) {
                     * residueName = "ALA".intern();
                     * residue.setName(residueName); }
                     */
                }
            }

            aminoAcid = AminoAcid3.UNK;
            aminoAcidNumber = -1;
            for (AminoAcid3 amino : aminoAcidList) {
                aminoAcidNumber++;
                if (amino.toString().equalsIgnoreCase(residueName)) {
                    aminoAcid = amino;
                    break;
                }
            }

            /**
             * Non-standard Amino Acid; use ALA backbone types.
             */
            if (aminoAcid == AminoAcid3.UNK) {
                residueName = "ALA";
                aminoAcidNumber = -1;
                for (AminoAcid3 amino : aminoAcidList) {
                    aminoAcidNumber++;
                    if (amino.toString().equalsIgnoreCase(residueName)) {
                        break;
                    }
                }
                residueName = residue.getName().toUpperCase();
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
                    logger.warning(format(" An atom for residue %s has the wrong number of bonds:\n %s", residueName, atom.toString()));
                    logger.warning(format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
                }
            }
            /**
             * Remember the current C-alpha and carbonyl C atoms for use with
             * the next residue.
             */
            pCA = CA;
            pC = C;
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
     *
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     */
    private void assignAminoAcidSideChain(ResiduePosition position, AminoAcid3 aminoAcid, Residue residue,
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
                buildAlanine(residue, CA, N, C, k);
                break;
            case VAL:
                buildValine(residue, CA, N, C, k);
                break;
            case LEU:
                buildLeucine(residue, CA, N, C, k);
                break;
            case ILE:
                buildIsoleucine(residue, CA, N, C, k);
                break;
            case SER:
                buildSerine(residue, CA, N, C, k);
                break;
            case THR:
                buildThreonine(residue, CA, N, C, k);
                break;
            case CYS:
                buildCysteine(residue, CA, N, C, k);
                break;
            case CYX:
                buildCystine(residue, CA, N, C, k);
                break;
            case CYD:
                buildDeprotonatedCysteine(residue, CA, N, C, k);
                break;
            case PRO:
                buildProline(residue, CA, N, C, k, position);
                break;
            case PHE:
                buildPhenylalanine(residue, CA, N, C, k);
                break;
            case TYR:
                buildTyrosine(residue, CA, N, C, k);
                break;
            case TYD:
                buildDeprotonatedTyrosine(residue, CA, N, C, k);
                break;
            case TRP:
                buildTryptophan(residue, CA, N, C, k);
                break;
            case HIS:
                buildHistidine(residue, CA, N, C, k);
                break;
            case HID:
                buildNeutralHistidineD(residue, CA, N, C, k);
                break;
            case HIE:
                buildNeutralHistidineE(residue, CA, N, C, k);
                break;
            case ASP:
                buildAspartate(residue, CA, N, C, k);
                break;
            case ASH:
                buildNeutralAsparticAcid(residue, CA, N, C, k);
                break;
            case ASN:
                buildAsparagine(residue, CA, N, C, k);
                break;
            case GLU:
                buildGlutamate(residue, CA, N, C, k);
                break;
            case GLH:
                buildNeutralGlutamicAcid(residue, CA, N, C, k);
                break;
            case GLN:
                buildGlutamine(residue, CA, N, C, k);
                break;
            case MET:
                buildMethionine(residue, CA, N, C, k);
                break;
            case LYS:
                buildLysine(residue, CA, N, C, k);
                break;
            case LYD:
                buildDeprotonatedLysine(residue, CA, N, C, k);
                break;
            case ARG:
                buildArginine(residue, CA, N, C, k);
                break;
            case ORN:
                buildOrnithine(residue, CA, N, C, k);
                break;
            case AIB:
                buildAIB(residue, CA, N, C, k);
                break;
            case PCA:
                buildPCA(residue, CA, N, C, k);
                break;
            case UNK:
                String residueName = residue.getName();
                logger.info(" Patching side-chain " + residueName);
                HashMap<String, AtomType> types = forceField.getAtomTypes(residueName);
                HashMap<AtomType, AtomType> typeMap = new HashMap<AtomType, AtomType>();
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
                                logger.info(" No atom type was found for " + atomName + " of " + residueName + ".");
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
                                logger.info(" Missing heavy atom " + type.name);
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
                        HashMap<String, Atom> atomMap = new HashMap<String, Atom>();
                        for (Atom atom : residueAtoms) {
                            atomMap.put(atom.getName().toUpperCase(), atom);
                        }
                        for (String atomName : types.keySet()) {
                            AtomType type = types.get(atomName);
                            String bonds[] = forceField.getBonds(residueName, atomName.toUpperCase());
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
                                            logger.info(" Check biotype for hydrogen " + atomName + ".");
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
                        logger.log(Level.SEVERE, format(" Could not patch %s.", residueName));
                    } else {
                        logger.info(" Patch for " + residueName + " succeeded.");
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
     * This exception is thrown when an atom type could not be assigned.
     */
    private class MissingAtomTypeException extends Exception {

        public final Residue residue;
        public final Atom atom;

        public MissingAtomTypeException(Residue residue, Atom atom) {
            super();
            this.residue = residue;
            this.atom = atom;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(super.toString());
            sb.append(format(" Atom %s", atom.toString()));
            if (residue != null) {
                sb.append(format("\n of residue %s", residue.toString()));
            }
            sb.append("\n could not be assigned an atom type.\n");
            return sb.toString();
        }
    }

    /**
     * This exception is thrown when a heavy atom is not found.
     */
    private class MissingHeavyAtomException extends Exception {

        public final String atomName;
        public final AtomType atomType;
        public final Atom bondedTo;

        public MissingHeavyAtomException(String atomName, AtomType atomType, Atom bondedTo) {
            super();
            this.atomName = atomName;
            this.atomType = atomType;
            this.bondedTo = bondedTo;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(super.toString());
            if (atomType != null) {
                sb.append(format("\n An atom of type\n %s\n", atomType.toString()));
            } else {
                sb.append(format("\n Atom %s", atomName));
            }
            sb.append(" was not found");
            if (bondedTo != null) {
                sb.append(format(" bonded to atom %s ", bondedTo));
            }
            sb.append(".\n");
            return sb.toString();
        }
    }

    /**
     * Determine the atom type based on a biotype key.
     *
     * @param key The biotype key.
     * @return The atom type.
     * @since 1.0
     */
    private AtomType findAtomType(int key) {
        BioType bioType = forceField.getBioType(Integer.toString(key));
        if (bioType != null) {
            AtomType atomType = forceField.getAtomType(Integer.toString(bioType.atomType));
            if (atomType != null) {
                return atomType;
            } else {
                logger.severe(format("The atom type %s was not found for biotype %s.", bioType.atomType,
                        bioType.toString()));
            }
        }
        /*
         * else { logger.severe(format("The biotype %s was not found.",
         * bioType.toString())); }
         */
        return null;
    }

    /**
     * <p>padRight</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padRight(String s, int n) {
        return String.format("%-" + n + "s", s);
    }

    /**
     * <p>padLeft</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padLeft(String s, int n) {
        return String.format("%" + n + "s", s);
    }

    /**
     * <p>writeFileWithHeader</p>
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
            fw = new FileWriter(newFile, false);
            bw = new BufferedWriter(fw);
            bw.write(header.toString());
            bw.close();
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
                bw.write(format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s\n", c.a, c.b, c.c, c.alpha, c.beta,
                        c.gamma, padRight(c.spaceGroup.pdbName, 10)));
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
                                        bw.write(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n",
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
                    bw.write(terSB.toString());
                    bw.newLine();
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

            bw.write("END");
            bw.newLine();
            bw.close();
        } catch (Exception e) {
            String message = "Exception writing to file: " + saveFile.toString();
            logger.log(Level.WARNING, message, e);
            return false;
        }
        return true;
    }

    /**
     * <p>writeAtom</p>
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
        bw.write(sb.toString());
        bw.newLine();
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
            bw.write(anisouSB.toString());
            bw.newLine();
        }
    }

    /**
     * The location of a residue within a chain.
     */
    public enum ResiduePosition {

        FIRST_RESIDUE, MIDDLE_RESIDUE, LAST_RESIDUE
    };

    private enum HetAtoms {

        HOH, H2O, WAT, NA, K, MG, MG2, CA, CA2, CL, BR, ZN, ZN2
    };
    /**
     * Biotype keys for nucleic acid backbone atom types. These are consistent
     * with parameter files from TINKER v. 6.1 (June 2012).
     */
    private static final int o5Typ[] = {
        1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0,
        1300, 1334, 1363, 1400, 1431, 1465, 1492, 1523, 1559, 1589, 1624
    };
    private static final int c5Typ[] = {
        1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0,
        1301, 1335, 1364, 1401, 1432, 1466, 1493, 1524, 1560, 1590, 1625
    };
    private static final int h51Typ[] = {
        1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0,
        1302, 1336, 1365, 1402, 1433, 1467, 1494, 1525, 1561, 1591, 1626
    };
    private static final int h52Typ[] = {
        1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0,
        1303, 1337, 1366, 1403, 1434, 1468, 1495, 1526, 1562, 1592, 1627
    };
    private static final int c4Typ[] = {
        1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0,
        1304, 1338, 1367, 1404, 1435, 1469, 1496, 1527, 1563, 1593, 1628
    };
    private static final int h4Typ[] = {
        1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0,
        1305, 1339, 1368, 1405, 1436, 1470, 1497, 1528, 1564, 1594, 1629
    };
    private static final int o4Typ[] = {
        1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0,
        1306, 1340, 1369, 1406, 1437, 1471, 1498, 1529, 1565, 1595, 1630
    };
    private static final int c1Typ[] = {
        1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0,
        1307, 1341, 1370, 1407, 1438, 1472, 1499, 1530, 1566, 1596, 1631
    };
    private static final int h1Typ[] = {
        1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0,
        1308, 1342, 1371, 1408, 1439, 1473, 1500, 1531, 1567, 1597, 1632
    };
    private static final int c3Typ[] = {
        1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0,
        1309, 1343, 1372, 1409, 1440, 1474, 1501, 1532, 1568, 1598, 1633
    };
    private static final int h3Typ[] = {
        1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0,
        1310, 1344, 1373, 1410, 1441, 1475, 1502, 1533, 1569, 1599, 1634
    };
    private static final int c2Typ[] = {
        1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0,
        1311, 1345, 1374, 1411, 1442, 1476, 1503, 1534, 1570, 1600, 1635
    };
    private static final int h21Typ[] = {
        1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0,
        1312, 1346, 1375, 1412, 1443, 1477, 1504, 1535, 1571, 1601, 1636
    };
    private static final int o2Typ[] = {
        1014, 1044, 1075, 1103, 0, 0, 0, 0, 0, 0, 0, 0,
        1313, 1347, 1376, 1413, 1444, 1478, 1505, 1536, 1572, 1602, 1637
    };
    private static final int h22Typ[] = {
        1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0,
        1314, 1348, 1377, 0, 0, 1479, 1506, 1537, 1573, 1603, 1638
    };
    private static final int o3Typ[] = {
        1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0,
        1315, 1349, 1378, 1414, 1445, 1480, 1507, 1538, 1574, 1604, 1639
    };
    private static final int pTyp[] = {
        1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0,
        1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230
    };
    private static final int opTyp[] = {
        1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0,
        1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231
    };
    private static final int h5tTyp[] = {
        1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0,
        1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233
    };
    private static final int h3tTyp[] = {
        1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0,
        1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238
    };
    /**
     * Biotype keys for amino acid backbone atom types. These are consistent
     * with parameter files from TINKER v. 6.1 (June 2012).
     * <br>
     * xType[0][..] are for N-terminal residues.
     * <br>
     * xType[1][..] are mid-chain residues.
     * <br>
     * xType[2][..] are for C-terminal residues.
     */
    private static final int nType[][] = {
        {
            403, 409, 415, 421, 427, 433, 439, 445,
            451, 457, 463, 471, 477, 483, 489, 495,
            501, 507, 513, 519, 525, 531, 537, 543,
            549, 555, 561, 567, 573, 579, 391, 762,
            0, 0, 0, 0, 0, 403
        }, {
            1, 7, 15, 27, 41, 55, 65, 77,
            87, 96, 105, 116, 131, 147, 162, 185,
            202, 218, 234, 244, 256, 268, 280, 294,
            308, 321, 337, 353, 370, 384, 391, 0,
            0, 0, 0, 0, 0, 1
        }, {
            584, 590, 596, 602, 608, 614, 620, 626,
            632, 638, 644, 649, 655, 661, 667, 673,
            679, 685, 691, 697, 703, 709, 715, 721,
            727, 733, 739, 745, 751, 757, 0, 0,
            0, 0, 773, 775, 777, 584
        }
    };
    private static final int caType[][] = {
        {
            404, 410, 416, 422, 428, 434, 440, 446,
            452, 458, 464, 472, 478, 484, 490, 496,
            502, 508, 514, 520, 526, 532, 538, 544,
            550, 556, 562, 568, 574, 580, 392, 0,
            0, 767, 0, 0, 0, 404
        }, {
            2, 8, 16, 28, 42, 56, 66, 78,
            88, 97, 106, 117, 132, 148, 163, 186,
            203, 219, 235, 245, 257, 269, 281, 295,
            309, 322, 338, 354, 371, 385, 392, 0,
            0, 0, 0, 0, 0, 2
        }, {
            585, 591, 597, 603, 609, 615, 621, 627,
            633, 639, 645, 650, 656, 662, 668, 674,
            680, 686, 692, 698, 704, 710, 716, 722,
            728, 734, 740, 746, 752, 758, 0, 0,
            0, 0, 0, 0, 779, 585
        }
    };
    private static final int cType[][] = {
        {
            405, 411, 417, 423, 429, 435, 441, 447,
            453, 459, 465, 473, 479, 485, 491, 497,
            503, 509, 515, 521, 527, 533, 539, 545,
            551, 557, 563, 569, 575, 581, 393, 0,
            764, 769, 0, 0, 0, 405
        }, {
            3, 9, 17, 29, 43, 57, 67, 79,
            89, 98, 107, 118, 133, 149, 164, 187,
            204, 220, 236, 246, 258, 270, 282, 296,
            310, 323, 339, 355, 372, 386, 393, 0,
            0, 0, 0, 0, 0, 3
        }, {
            586, 592, 598, 604, 610, 616, 622, 628,
            634, 640, 646, 651, 657, 663, 669, 675,
            681, 687, 693, 699, 705, 711, 717, 723,
            729, 735, 741, 747, 753, 759, 0, 0,
            0, 0, 771, 0, 0, 586
        }
    };
    private static final int hnType[][] = {
        {
            406, 412, 418, 424, 430, 436, 442, 448,
            454, 460, 466, 474, 480, 486, 492, 498,
            504, 510, 516, 522, 528, 534, 540, 546,
            552, 558, 564, 570, 576, 582, 394, 763,
            0, 0, 0, 0, 0, 406
        }, {
            4, 10, 18, 30, 44, 58, 68, 80,
            90, 99, 0, 119, 134, 150, 165, 188,
            205, 221, 237, 247, 259, 271, 283, 297,
            311, 324, 340, 356, 373, 387, 394, 0,
            0, 0, 0, 0, 0, 4
        }, {
            587, 593, 599, 605, 611, 617, 623, 629,
            635, 641, 0, 652, 658, 664, 670, 676,
            682, 688, 694, 700, 706, 712, 718, 724,
            730, 736, 742, 748, 754, 760, 0, 0,
            0, 0, 774, 776, 778, 587}
    };
    private static final int oType[][] = {
        {
            407, 413, 419, 425, 431, 437, 443, 449,
            455, 461, 467, 475, 481, 487, 493, 499,
            505, 511, 517, 523, 529, 535, 541, 547,
            553, 559, 565, 571, 577, 583, 395, 0,
            766, 770, 0, 0, 0, 407
        }, {
            5, 11, 19, 31, 45, 59, 69, 81,
            91, 100, 108, 120, 135, 151, 166, 189,
            206, 222, 238, 248, 260, 272, 284, 298,
            312, 325, 341, 357, 374, 388, 395, 0,
            0, 0, 0, 0, 0, 5
        }, {
            588, 594, 600, 606, 612, 618, 624, 630,
            636, 642, 647, 653, 659, 665, 671, 677,
            683, 689, 695, 701, 707, 713, 719, 725,
            731, 737, 743, 749, 755, 761, 0, 0,
            0, 0, 772, 0, 0, 588
        }
    };
    private static final int haType[][] = {
        {
            408, 414, 420, 426, 432, 438, 444, 450,
            456, 462, 468, 476, 482, 488, 494, 500,
            506, 512, 518, 524, 530, 536, 542, 548,
            554, 560, 566, 572, 578, 0, 396, 0,
            765, 768, 0, 0, 0, 408},
        {
            6, 12, 20, 32, 46, 60, 70, 82,
            92, 101, 109, 121, 136, 152, 167, 190,
            207, 223, 239, 249, 261, 273, 285, 299,
            313, 326, 342, 358, 375, 0, 396, 0,
            0, 0, 0, 0, 0, 6},
        {
            589, 595, 601, 607, 613, 619, 625, 631,
            637, 643, 648, 654, 660, 666, 672, 678,
            684, 690, 696, 702, 708, 714, 720, 726,
            732, 738, 744, 750, 756, 0, 0, 0,
            0, 0, 0, 0, 780, 589
        }
    };
    private static final int cbType[] = {
        0, 13, 21, 33, 47, 61, 71, 83,
        93, 102, 110, 122, 137, 153, 168, 191,
        208, 224, 240, 250, 262, 274, 286, 300,
        314, 327, 343, 359, 376, 389, 397, 0,
        0, 0, 0, 0, 0, 0
    };

    private void buildAlanine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom HB1 = buildHydrogen(residue, "HB1", CB, 1.11, CA, 109.4, N, 180.0, 0, k + 1);
        buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, HB1, 109.4, 1, k + 1);
        buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, HB1, 109.4, -1, k + 1);
    }

    private void buildValine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG1 = buildHeavy(residue, "CG1", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, CG1, 109.5, -1, k + 4);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, CG1, 109.4, 1, k + 1);
        Atom HG11 = buildHydrogen(residue, "HG11", CG1, 1.11, CB, 109.4, CA, 180.0, 0, k + 3);
        Atom HG12 = buildHydrogen(residue, "HG12", CG1, 1.11, CB, 109.4, HG11, 109.4, 1, k + 3);
        Atom HG13 = buildHydrogen(residue, "HG13", CG1, 1.11, CB, 109.4, HG11, 109.4, -1, k + 3);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 109.4, CA, 180.0, 0, k + 5);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 109.4, HG21, 109.4, 1, k + 5);
        Atom HG23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 109.4, HG21, 109.4, -1, k + 5);
    }

    private void buildLeucine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.54, CB, 109.5, CD1, 109.5, -1, k + 6);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG = buildHydrogen(residue, "HG", CG, 1.11, CB, 109.4, CD1, 109.4, 1, k + 3);
        Atom HD11 = buildHydrogen(residue, "HD11", CD1, 1.11, CG, 109.4, CB, 180.0, 0, k + 5);
        Atom HD12 = buildHydrogen(residue, "HD12", CD1, 1.11, CG, 109.4, HD11, 109.4, 1, k + 5);
        Atom HD13 = buildHydrogen(residue, "HD13", CD1, 1.11, CG, 109.4, HD11, 109.4, -1, k + 5);
        Atom HD21 = buildHydrogen(residue, "HD21", CD2, 1.11, CG, 109.4, CB, 180.0, 0, k + 7);
        Atom HD22 = buildHydrogen(residue, "HD22", CD2, 1.11, CG, 109.4, HD21, 109.4, 1, k + 7);
        Atom HD23 = buildHydrogen(residue, "HD23", CD2, 1.11, CG, 109.4, HD21, 109.4, -1, k + 7);
    }

    private void buildIsoleucine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom CG1 = buildHeavy(residue, "CG1", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, CG1, 109.5, -1, k + 4);
        Atom CD1 = buildHeavy(residue, "CD1", CG1, 1.54, CB, 109.5, CA, 180, 0, k + 6);
        //  CD1 = setHeavy(residue, "CD", CG1, 1.54, CB, 109.5, CA, 180, 0, k + 6);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, CG1, 109.4, 1, k + 1);
        Atom HG12 = buildHydrogen(residue, "HG12", CG1, 1.11, CB, 109.4, CD1, 109.4, 1, k + 3);
        Atom HG13 = buildHydrogen(residue, "HG13", CG1, 1.11, CB, 109.4, CD1, 109.4, -1, k + 3);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 110.0, CG1, 180.0, 0, k + 5);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 110.0, HG21, 109.0, 1, k + 5);
        Atom HG23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 110.0, HG21, 109.0, -1, k + 5);
        Atom HD11 = buildHydrogen(residue, "HD11", CD1, 1.11, CG1, 110.0, CB, 180.0, 0, k + 7);
        Atom HD12 = buildHydrogen(residue, "HD12", CD1, 1.11, CG1, 110.0, HD11, 109.0, 1, k + 7);
        Atom HD13 = buildHydrogen(residue, "HD13", CD1, 1.11, CG1, 110.0, HD11, 109.0, -1, k + 7);
    }

    private Residue buildSerine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom OG = buildHeavy(residue, "OG", CB, 1.41, CA, 107.5, N, 180, 0, k + 2);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, OG, 106.7, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, OG, 106.7, -1, k + 1);
        Atom HG = buildHydrogen(residue, "HG", OG, 0.94, CB, 106.9, CA, 180.0, 0, k + 3);
        return residue;
    }

    private Residue buildThreonine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom OG1 = buildHeavy(residue, "OG1", CB, 1.41, CA, 107.5, N, 180, 0, k + 2);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, OG1, 107.7, 1, k + 4);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, OG1, 106.7, -1, k + 1);
        Atom HG1 = buildHydrogen(residue, "HG1", OG1, 0.94, CB, 106.9, CA, 180.0, 0, k + 3);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 110.0, CA, 180.0, 0, k + 5);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 110.0, HG21, 109.0, 1, k + 5);
        Atom HD23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 110.0, HG21, 109.0, -1, k + 5);
        return residue;
    }

    private Residue buildCysteine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1);
        Atom HG = buildHydrogen(residue, "HG", SG, 1.34, CB, 96.0, CA, 180.0, 0, k + 3);
        return residue;
    }

    private Residue buildDeprotonatedCysteine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1);
        return residue;
    }

    private Residue buildCystine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1);
        List<Atom> resAtoms = residue.getAtomList();
        for (Atom atom : resAtoms) {
            atom.setResName("CYS");
        }
        residue.setName("CYS");
        return residue;
    }

    private Residue buildProline(Residue residue, Atom CA, Atom N, Atom C, int k, ResiduePosition position) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 107.0, C, 109.5, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 107.0, N, 180, 0, k + 2);
        Atom CD;
        if (position == FIRST_RESIDUE) {
            CD = buildHeavy(residue, "CD", CG, 1.54, CB, 107.0, CA, 180, 0, 469);
        } else {
            CD = buildHeavy(residue, "CD", CG, 1.54, CB, 107.0, CA, 180, 0, k + 4);
        }
        buildBond(CD, N);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        if (position == FIRST_RESIDUE) {
            buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, N, 109.4, 1, 470);
            buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, N, 109.4, -1, 470);
        } else {
            buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, N, 109.4, 1, k + 5);
            buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, N, 109.4, -1, k + 5);
        }
        return residue;
    }

    private Residue buildPhenylalanine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7);
        buildBond(CE2, CZ);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.11, CG, 120.0, CE1, 120.0, 1, k + 4);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.11, CG, 120.0, CE2, 120.0, 1, k + 4);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.11, CD1, 120.0, CZ, 120.0, 1, k + 6);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.11, CD2, 120.0, CZ, 120.0, 1, k + 6);
        Atom HZ = buildHydrogen(residue, "HZ", CZ, 1.11, CE1, 120.0, CE2, 120.0, 1, k + 8);
        return residue;
    }

    private Residue buildTyrosine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7);
        buildBond(CE2, CZ);
        Atom OH = buildHeavy(residue, "OH", CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, k + 8);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 120.0, CE1, 120.0, 1, k + 4);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 120.0, CE2, 120.0, 1, k + 4);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, k + 6);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, k + 6);
        Atom HH = buildHydrogen(residue, "HH", OH, 0.97, CZ, 108.0, CE2, 0.0, 0, k + 9);
        return residue;
    }

    private Residue buildDeprotonatedTyrosine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7);
        buildBond(CE2, CZ);
        Atom OH = buildHeavy(residue, "OH", CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, k + 8);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 120.0, CE1, 120.0, 1, k + 4);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 120.0, CE2, 120.0, 1, k + 4);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, k + 6);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, k + 6);
        return residue;
    }

    private Residue buildTryptophan(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 0.0, 0, k + 2);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, CD1, 108.0, 1, k + 5);
        Atom NE1 = buildHeavy(residue, "NE1", CD1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 6);
        Atom CE2 = buildHeavy(residue, "CE2", NE1, 1.35, CD1, 108.0, CG, 0.0, 0, k + 8);
        buildBond(CE2, CD2);
        Atom CE3 = buildHeavy(residue, "CE3", CD2, 1.35, CE2, 120.0, NE1, 180.0, 0, k + 9);
        Atom CZ2 = buildHeavy(residue, "CZ2", CE2, 1.35, CD2, 120.0, CE3, 0.0, 0, k + 11);
        Atom CZ3 = buildHeavy(residue, "CZ3", CE3, 1.35, CD2, 120.0, CE2, 0.0, 0, k + 13);
        Atom CH2 = buildHeavy(residue, "CH2", CZ2, 1.35, CE2, 120.0, CD2, 0.0, 0, k + 15);
        buildBond(CH2, CZ3);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 126.0, NE1, 126.0, 1, k + 4);
        Atom HE1 = buildHydrogen(residue, "HE1", NE1, 1.05, CD1, 126.0, CE2, 126.0, 1, k + 7);
        Atom HE3 = buildHydrogen(residue, "HE3", CE3, 1.10, CD1, 120.0, CZ3, 120.0, 1, k + 10);
        Atom HZ2 = buildHydrogen(residue, "HZ2", CZ2, 1.10, CE2, 120.0, CH2, 120.0, 1, k + 12);
        Atom HZ3 = buildHydrogen(residue, "HZ3", CZ3, 1.10, CE3, 120.0, CH2, 120.0, 1, k + 14);
        Atom HH2 = buildHydrogen(residue, "HH2", CH2, 1.10, CZ2, 120.0, CZ3, 120.0, 1, k + 16);
        return residue;
    }

    private Residue buildHistidine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 5);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 7);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 9);
        buildBond(NE2, CE1);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 6);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 8);
        Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 10);
        return residue;
    }

    private Residue buildNeutralHistidineD(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 5);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 7);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 9);
        buildBond(NE2, CE1);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 6);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 8);
        return residue;
    }

    private Residue buildNeutralHistidineE(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 4);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 6);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 8);
        buildBond(NE2, CE1);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 5);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 7);
        Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 9);
        return residue;
    }

    private Residue buildAspartate(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.25, CB, 117.0, CA, 0.0, 0, k + 3);
        Atom OD2 = buildHeavy(residue, "OD2", CG, 1.25, CB, 117.0, OD1, 126.0, 1, k + 3);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1);
        return residue;
    }

    private Residue buildNeutralAsparticAcid(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.25, CB, 117.0, CA, 0.0, 0, k + 3);
        Atom OD2 = buildHeavy(residue, "OD2", CG, 1.25, CB, 117.0, OD1, 126.0, 1, k + 4);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1);
        Atom HD2 = buildHydrogen(residue, "HD2", OD2, 0.98, CG, 108.7, OD1, 0.0, 0, k + 5);
        return residue;
    }

    private Residue buildAsparagine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.22, CB, 122.5, CA, 180, 0, k + 3);
        Atom ND2 = buildHeavy(residue, "ND2", CG, 1.34, CB, 112.7, OD1, 124.0, 1, k + 4);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1);
        Atom HD21 = buildHydrogen(residue, "HD21", ND2, 1.02, CG, 119.0, CB, 0.0, 0, k + 5);
        Atom HD22 = buildHydrogen(residue, "HD22", ND2, 1.02, CG, 119.0, HD21, 120.0, 1, k + 5);
        return residue;
    }

    private Residue buildGlutamate(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.25, CG, 117.0, CB, 180, 0, k + 5);
        Atom OE2 = buildHeavy(residue, "OE2", CD, 1.25, CG, 117.0, OE1, 126.0, 1, k + 5);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3);
        return residue;
    }

    private Residue buildNeutralGlutamicAcid(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.25, CG, 117.0, CB, 180, 0, k + 5);
        Atom OE2 = buildHeavy(residue, "OE2", CD, 1.25, CG, 117.0, OE1, 126.0, 1, k + 6);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3);
        Atom HE2 = buildHydrogen(residue, "HE2", OE2, 0.98, CD, 108.7, OE1, 0.0, 0, k + 7);
        return residue;
    }

    public Residue buildGlutamine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.22, CG, 122.5, CB, 180, 0, k + 5);
        Atom NE2 = buildHeavy(residue, "NE2", CD, 1.34, CG, 112.7, OE1, 124.0, 1, k + 6);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3);
        Atom HE21 = buildHydrogen(residue, "HE21", NE2, 1.02, CD, 119.0, CG, 0.0, 0, k + 7);
        Atom HE22 = buildHydrogen(residue, "HE22", NE2, 1.02, CD, 119.0, HE21, 120.0, 1, k + 7);
        return residue;
    }

    public Residue buildMethionine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom SD = buildHeavy(residue, "SD", CG, 1.82, CB, 109.0, CA, 180, 0, k + 4);
        Atom CE = buildHeavy(residue, "CE", SD, 1.82, CG, 96.3, CB, 180, 0, k + 5);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, SD, 112.0, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, SD, 112.0, -1, k + 3);
        Atom HE1 = buildHydrogen(residue, "HE1", CE, 1.11, SD, 112.0, CG, 180.0, 0, k + 6);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, SD, 112.0, HE1, 109.4, 1, k + 6);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, SD, 112.0, HE1, 109.4, -1, k + 6);
        return residue;
    }

    public Residue buildLysine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom CE = buildHeavy(residue, "CE", CD, 1.54, CG, 109.5, CB, 180, 0, k + 6);
        Atom NZ = buildHeavy(residue, "NZ", CE, 1.50, CD, 109.5, CG, 180, 0, k + 8);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, CE, 109.4, 1, k + 5);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, CE, 109.4, -1, k + 5);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, CD, 109.4, NZ, 108.8, 1, k + 7);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, CD, 109.4, NZ, 108.8, -1, k + 7);
        Atom HZ1 = buildHydrogen(residue, "HZ1", NZ, 1.02, CE, 109.5, CD, 180.0, 0, k + 9);
        Atom HZ2 = buildHydrogen(residue, "HZ2", NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, k + 9);
        Atom HZ3 = buildHydrogen(residue, "HZ3", NZ, 1.02, CE, 109.5, HZ1, 109.5, -1, k + 9);
        return residue;
    }

    public Residue buildDeprotonatedLysine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom CE = buildHeavy(residue, "CE", CD, 1.54, CG, 109.5, CB, 180, 0, k + 6);
        Atom NZ = buildHeavy(residue, "NZ", CE, 1.50, CD, 109.5, CG, 180, 0, k + 8);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, CE, 109.4, 1, k + 5);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, CE, 109.4, -1, k + 5);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, CD, 109.4, NZ, 108.8, 1, k + 7);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, CD, 109.4, NZ, 108.8, -1, k + 7);
        Atom HZ1 = buildHydrogen(residue, "HZ1", NZ, 1.02, CE, 109.5, CD, 180.0, 0, k + 9);
        Atom HZ2 = buildHydrogen(residue, "HZ2", NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, k + 9);
        return residue;
    }

    public Residue buildArginine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom NE = buildHeavy(residue, "NE", CD, 1.45, CG, 109.5, CB, 180, 0, k + 6);
        Atom CZ = buildHeavy(residue, "CZ", NE, 1.35, CD, 120.0, CG, 180, 0, k + 8);
        Atom NH1 = buildHeavy(residue, "NH1", CZ, 1.35, NE, 120.0, CD, 180, 0, k + 9);
        Atom NH2 = buildHeavy(residue, "NH2", CZ, 1.35, NE, 120.0, NH1, 120.0, 1, k + 9);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, NE, 109.4, 1, k + 5);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, NE, 109.4, -1, k + 5);
        Atom HE = buildHydrogen(residue, "HE", NE, 1.02, CD, 120.0, CZ, 120.0, 1, k + 7);
        Atom HH11 = buildHydrogen(residue, "HH11", NH1, 1.02, CZ, 120.0, NE, 180.0, 0, k + 10);
        Atom HH12 = buildHydrogen(residue, "HH12", NH1, 1.02, CZ, 120.0, HH11, 120.0, 1, k + 10);
        Atom HH21 = buildHydrogen(residue, "HH21", NH2, 1.02, CZ, 120.0, NE, 180.0, 0, k + 10);
        Atom HH22 = buildHydrogen(residue, "HH22", NH2, 1.02, CZ, 120.0, HH21, 120.0, 1, k + 10);
        return residue;
    }

    public Residue buildOrnithine(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom NE = buildHeavy(residue, "NE", CD, 1.50, CG, 109.5, CB, 180, 0, k + 6);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, NE, 109.4, 1, k + 5);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, NE, 109.4, -1, k + 5);
        Atom HE1 = buildHydrogen(residue, "HE1", NE, 1.02, CD, 109.5, CG, 180.0, 0, k + 7);
        Atom HE2 = buildHydrogen(residue, "HE2", NE, 1.02, CD, 109.5, HE1, 109.5, 1, k + 7);
        Atom HE3 = buildHydrogen(residue, "HE3", NE, 1.02, CD, 109.5, HE1, 109.5, -1, k + 7);
        return residue;
    }

    public Residue buildAIB(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB1 = buildHeavy(residue, "CB1", CA, 1.54, N, 109.5, C, 107.8, -1, k);
        Atom CB2 = buildHeavy(residue, "CB1", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom HB11 = buildHydrogen(residue, "HB11", CB1, 1.11, CA, 109.4, N, 180.0, 0, k + 1);
        Atom HB12 = buildHydrogen(residue, "HB12", CB1, 1.11, CA, 109.4, HB11, 109.4, 1, k + 1);
        Atom HB13 = buildHydrogen(residue, "HB13", CB1, 1.11, CA, 109.4, HB11, 109.4, -1, k + 1);
        Atom HG21 = buildHydrogen(residue, "HG21", CB2, 1.11, CA, 109.4, N, 180.0, 0, k + 1);
        Atom HG22 = buildHydrogen(residue, "HG22", CB2, 1.11, CA, 109.4, HG21, 109.4, 1, k + 1);
        Atom HG23 = buildHydrogen(residue, "HG23", CB2, 1.11, CA, 109.4, HG21, 109.4, -1, k + 1);
        return residue;
    }

    public Residue buildPCA(Residue residue, Atom CA, Atom N, Atom C, int k) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4);
        Atom OE = buildHeavy(residue, "OE", CD, 1.22, N, 126.0, CG, 126.0, 1, k + 5);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3);
        return residue;
    }

    private Atom buildHeavy(MSGroup residue, String atomName, Atom bondedTo, int key)
            throws MissingHeavyAtomException {
        Atom atom = (Atom) residue.getAtomNode(atomName);
        AtomType atomType = findAtomType(key);
        if (atom == null) {
            MissingHeavyAtomException missingHeavyAtom = new MissingHeavyAtomException(atomName, atomType, bondedTo);
            throw missingHeavyAtom;
        }
        atom.setAtomType(atomType);
        if (bondedTo != null) {
            buildBond(atom, bondedTo);
        }
        return atom;
    }

    private Atom buildHeavy(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, int lookUp) {
        AtomType atomType = findAtomType(lookUp);
        return buildHeavyAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType);
    }

    private Atom buildHeavyAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, AtomType atomType) {
        Atom atom = (Atom) residue.getAtomNode(atomName);
        if (atomType == null) {
            return null;
        }
        if (atom == null) {
            String resName = ia.getResidueName();
            int resSeq = ia.getResidueNumber();
            Character chainID = ia.getChainID();
            Character altLoc = ia.getAltLoc();
            String segID = ia.getSegID();
            double occupancy = ia.getOccupancy();
            double tempFactor = ia.getTempFactor();
            atom = new Atom(0, atomName, altLoc, new double[3], resName, resSeq, chainID,
                    occupancy, tempFactor, segID);
            residue.addMSNode(atom);
            intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
        }
        atom.setAtomType(atomType);
        buildBond(ia, atom);
        return atom;
    }

    private Atom buildHydrogen(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, int lookUp) {
        AtomType atomType = findAtomType(lookUp);
        return buildHydrogenAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType);
    }

    private Atom buildHydrogenAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, AtomType atomType) {
        if (atomType == null) {
            return null;
        }
        Atom atom = (Atom) residue.getAtomNode(atomName);
        // It may be a Deuterium
        if (atom == null) {
            String dAtomName = atomName.replaceFirst("H", "D");
            atom = (Atom) residue.getAtomNode(dAtomName);
        }
        if (atom == null) {
            String resName = ia.getResidueName();
            int resSeq = ia.getResidueNumber();
            Character chainID = ia.getChainID();
            Character altLoc = ia.getAltLoc();
            String segID = ia.getSegID();
            double occupancy = ia.getOccupancy();
            double tempFactor = ia.getTempFactor();
            atom = new Atom(0, atomName, altLoc, new double[3], resName, resSeq, chainID,
                    occupancy, tempFactor, segID);
            residue.addMSNode(atom);
            intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
        }
        atom.setAtomType(atomType);
        buildBond(ia, atom);
        return atom;
    }

    private Bond buildBond(Atom a1, Atom a2) {
        Bond bond = new Bond(a1, a2);
        int c[] = new int[2];
        c[0] = a1.getAtomType().atomClass;
        c[1] = a2.getAtomType().atomClass;
        String key = BondType.sortKey(c);
        BondType bondType = forceField.getBondType(key);
        if (bondType == null) {
            logger.severe(format("No BondType for key: %s\n %s\n %s\n %s\n %s", key,
                    a1.toString(), a1.getAtomType().toString(),
                    a2.toString(), a2.getAtomType().toString()));
        } else {
            bond.setBondType(bondType);
        }
        bondList.add(bond);
        return bond;
    }
}
