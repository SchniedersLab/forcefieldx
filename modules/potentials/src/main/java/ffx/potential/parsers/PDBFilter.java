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

import static ffx.potential.parsers.INTFilter.intxyz;
import static ffx.potential.parsers.PDBFilter.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.parsers.PDBFilter.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.parsers.PDBFilter.ResiduePosition.LAST_RESIDUE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import ffx.crystal.SpaceGroup;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSGroup;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.BioType;
import ffx.potential.parsers.PDBFilter.ResiduePosition;

/**
 * The PDBFilter class parses data from a Protein DataBank (*.PDB) file. The
 * following records are recognized: ANISOU, ATOM, CONECT, CRYST1, HELIX,
 * HETATM, LINK, SHEET, SSBOND, TURN, REMARK. The rest are ignored.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class PDBFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(PDBFilter.class.getName());
    private String pdbURL = null;

    public static String pdbForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".pdb.gz";
    }

    public static String cifForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".cif.gz";
    }
    /**
     * Keep track of altLoc Characters.
     */
    private Vector<Character> altLocs = new Vector<Character>();
    /**
     * Keep track of ATOM record serial numbers to match them with ANISOU
     * records.
     */
    private HashMap<Integer, Atom> atoms = new HashMap<Integer, Atom>();

    /**
     * Default Constructor
     */
    public PDBFilter() {
        super();
    }

    /**
     * Parse the PDB File from a local disk
     */
    public PDBFilter(MolecularAssembly f) {
        super(f);
    }

    /**
     * Parse the PDB File from a URL.
     */
    public PDBFilter(MolecularAssembly f, String pdb, String vrml) {
        super(f);
        pdbURL = pdb;
    }

    private enum Card {

        ANISOU, ATOM, CONECT, CRYST1, HELIX, HETATM, LINK, SHEET, SSBOND, TURN, REMARK
    };

    /**
     * Parse the PDB File
     */
    @Override
    public boolean readFile() {
        FileWriter fw = null;
        BufferedReader br = null;
        BufferedWriter bw = null;
        molecularAssembly.setFileType(FileType.PDB);
        try {
            setFileRead(false);
            if (pdbURL == null) {
                // Open a data stream to the PDB file
                File pdbFile = molecularAssembly.getFile();
                if (pdbFile == null || !pdbFile.exists() || !pdbFile.canRead()) {
                    return false;
                }
                FileReader fr = new FileReader(pdbFile);
                br = new BufferedReader(fr);
                if (logger.isLoggable(Level.INFO)) {
                    logger.info(" Opening " + pdbFile.getName());
                }
            } else {
                try {
                    URL url = new URL(pdbURL);
                    GZIPInputStream is = new GZIPInputStream(url.openStream());
                    br = new BufferedReader(new InputStreamReader(is));
                    int retry = 0;
                    while (!br.ready() && retry < 10) {
                        synchronized (this) {
                            if (logger.isLoggable(Level.INFO)) {
                                logger.info("Waiting on Network");
                            }
                            wait(50);
                            retry++;
                        }
                    }
                } catch (Exception e) {
                    logger.exiting(PDBFilter.class.getName(), "readFile", e);
                    return false;
                }
                // The downloaded PBD file will be echoed to the local file
                // system (if the file does not already exist)
                File pdbFile = molecularAssembly.getFile();
                if (pdbFile != null && !pdbFile.exists()) {
                    fw = new FileWriter(pdbFile);
                    bw = new BufferedWriter(fw);
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(" Saving to: " + pdbFile.getAbsolutePath());
                    }
                }
            }
            // First atom is #1, to match xyz file format
            int xyzIndex = 1;
            boolean useSegID = false;
            String[] connect;
            ArrayList<String[]> links = new ArrayList<String[]>();
            Vector<String[]> structs = new Vector<String[]>();
            // While the END parameter is not read in, load atoms
            String pdbLine = br.readLine();
            while ((pdbLine != null) && (!pdbLine.startsWith("END"))) {
                int len = pdbLine.length();
                String identity = pdbLine;
                if (len > 6) {
                    identity = pdbLine.substring(0, 6).trim().toUpperCase().intern();
                }
                Card card = null;
                try {
                    card = Card.valueOf(identity);
                } catch (Exception e) {
                    card = Card.REMARK;
                }
                switch (card) {
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
                        Integer serial = new Integer(pdbLine.substring(6, 11).trim());
                        String name = pdbLine.substring(12, 16).trim().intern();
                        Character altLoc = new Character(pdbLine.substring(16, 17).toUpperCase().charAt(0));
                        String resName = pdbLine.substring(17, 20).trim().intern();
                        String chainID = null;
                        if (!useSegID) {
                            chainID = pdbLine.substring(21, 22).intern();
                        } else {
                            chainID = pdbLine.substring(72, 76).intern();
                        }
                        if (chainID.equalsIgnoreCase(" ")) {
                            chainID = "Blank".intern();
                        }
                        int resSeq = Integer.decode(pdbLine.substring(22, 26).trim()).intValue();
                        double adp[] = new double[6];
                        adp[0] = new Integer(pdbLine.substring(28, 35).trim()) * 1.0e-4;
                        adp[1] = new Integer(pdbLine.substring(35, 42).trim()) * 1.0e-4;
                        adp[2] = new Integer(pdbLine.substring(42, 49).trim()) * 1.0e-4;
                        adp[3] = new Integer(pdbLine.substring(49, 56).trim()) * 1.0e-4;
                        adp[4] = new Integer(pdbLine.substring(56, 63).trim()) * 1.0e-4;
                        adp[5] = new Integer(pdbLine.substring(63, 70).trim()) * 1.0e-4;
                        if (atoms.containsKey(serial)) {
                            Atom a = atoms.get(serial);
                            a.setAltLoc(altLoc);
                            a.setAnisou(adp);
                        } else {
                            logger.info("No ATOM record for ANISOU serial number " + serial + ".");
                            logger.info("The following ANISOU record will be ignored:\n" + pdbLine);
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
                        serial = new Integer(pdbLine.substring(6, 11).trim());
                        name = pdbLine.substring(12, 16).trim().intern();
                        altLoc = new Character(pdbLine.substring(16, 17).toUpperCase().charAt(0));
                        if (!altLocs.contains(altLoc)) {
                            altLocs.add(altLoc);
                        }
                        resName = pdbLine.substring(17, 20).trim().intern();
                        chainID = null;
                        if (!useSegID) {
                            chainID = pdbLine.substring(21, 22).intern();
                        } else {
                            chainID = pdbLine.substring(72, 76).intern();
                        }
                        if (chainID.equalsIgnoreCase(" ")) {
                            chainID = "Blank".intern();
                        }
                        resSeq = new Integer(pdbLine.substring(22, 26).trim());
                        double d[] = new double[3];
                        d[0] = new Double(pdbLine.substring(30, 38).trim());
                        d[1] = new Double(pdbLine.substring(38, 46).trim());
                        d[2] = new Double(pdbLine.substring(46, 54).trim());
                        double occupancy = new Double(pdbLine.substring(54, 60).trim());
                        double tempFactor = new Double(pdbLine.substring(60, 66).trim());
                        Atom a = new Atom(0, name, altLoc, d, resName, resSeq, chainID, occupancy, tempFactor);
                        Atom prev = (Atom) molecularAssembly.contains(a);
                        if (prev != null) {
                            atoms.put(serial, prev);
                            prev.addAltLoc(altLoc, d, occupancy, tempFactor);
                        } else {
                            a.setXYZIndex(xyzIndex++);
                            atoms.put(serial, a);
                            molecularAssembly.addMSNode(a);
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
                        serial = new Integer(pdbLine.substring(6, 11).trim());
                        name = pdbLine.substring(12, 16).trim().intern();
                        altLoc = new Character(pdbLine.substring(16, 17).toUpperCase().charAt(0));
                        if (!altLocs.contains(altLoc)) {
                            altLocs.add(altLoc);
                        }
                        resName = pdbLine.substring(17, 20).trim().intern();
                        chainID = null;
                        if (!useSegID) {
                            chainID = pdbLine.substring(21, 22).intern();
                        } else {
                            chainID = pdbLine.substring(72, 76).intern();
                        }
                        if (chainID.equalsIgnoreCase(" ")) {
                            chainID = "Blank".intern();
                        }
                        resSeq = new Integer(pdbLine.substring(22, 26).trim());
                        d = new double[3];
                        d[0] = new Double(pdbLine.substring(30, 38).trim());
                        d[1] = new Double(pdbLine.substring(38, 46).trim());
                        d[2] = new Double(pdbLine.substring(46, 54).trim());
                        occupancy = new Double(pdbLine.substring(54, 60).trim());
                        tempFactor = new Double(pdbLine.substring(60, 66).trim());
                        a = new Atom(0, name, altLoc, d, resName, resSeq, chainID, occupancy, tempFactor);
                        prev = (Atom) molecularAssembly.contains(a);
                        if (prev != null) {
                            atoms.put(serial, prev);
                            prev.addAltLoc(altLoc, d, occupancy, tempFactor);
                        } else {
                            a.setXYZIndex(xyzIndex++);
                            a.setHetero(true);
                            atoms.put(serial, a);
                            molecularAssembly.addMSNode(a);
                        }
                        break;
                    case CONECT:
                        connect = new String[2];
                        connect[0] = pdbLine.substring(7, 11).trim();
                        connect[1] = pdbLine.substring(12, 16).trim();
                        links.add(connect);
                        break;
                    case CRYST1:
// =============================================================================
//  7 - 15       Real(9.3)     a              a (Angstroms).
// 16 - 24       Real(9.3)     b              b (Angstroms).
// 25 - 33       Real(9.3)     c              c (Angstroms).
// 34 - 40       Real(7.2)     alpha          alpha (degrees).
// 41 - 47       Real(7.2)     beta           beta (degrees).
// 48 - 54       Real(7.2)     gamma          gamma (degrees).
// 56 - 66       LString       sGroup         Space  group.
// 67 - 70       Integer       z              Z value.
// =============================================================================
                        double aaxis = new Double(pdbLine.substring(6, 15).trim());
                        double baxis = new Double(pdbLine.substring(15, 24).trim());
                        double caxis = new Double(pdbLine.substring(24, 33).trim());
                        double alpha = new Double(pdbLine.substring(33, 40).trim());
                        double beta = new Double(pdbLine.substring(40, 47).trim());
                        double gamma = new Double(pdbLine.substring(47, 54).trim());
                        int limit = 66;
                        if (len < 66) {
                            limit = len;
                        }
                        String sg = pdbLine.substring(55, limit).trim();
                        properties.addProperty("a-axis", aaxis);
                        properties.addProperty("b-axis", baxis);
                        properties.addProperty("c-axis", caxis);
                        properties.addProperty("alpha", alpha);
                        properties.addProperty("beta", beta);
                        properties.addProperty("gamma", gamma);
                        properties.addProperty("spacegroup", SpaceGroup.pdb2ShortName(sg));
                    case SSBOND:
                        connect = new String[6];
                        connect[0] = pdbLine.substring(15, 16);
                        // Polymers
                        connect[1] = pdbLine.substring(29, 30);
                        connect[2] = pdbLine.substring(17, 21).trim();
                        // Residues
                        connect[3] = pdbLine.substring(31, 35).trim();
                        connect[4] = new String("SG");
                        // Atoms
                        connect[5] = new String("SG");
                        links.add(connect);
                        break;
                    case LINK:
                        connect = new String[6];
                        connect[0] = pdbLine.substring(21, 22);
                        // Polymers
                        connect[1] = pdbLine.substring(51, 52);
                        connect[2] = pdbLine.substring(22, 26).trim();
                        // Residues
                        connect[3] = pdbLine.substring(52, 56).trim();
                        connect[4] = pdbLine.substring(12, 16).trim();
                        // Atoms
                        connect[5] = pdbLine.substring(42, 46).trim();
                        links.add(connect);
                        break;
                    case HELIX:
                        String[] struct = new String[6];
                        struct[0] = pdbLine.substring(0, 6).trim(); // HELIX
                        struct[1] = pdbLine.substring(19, 20); // Polymers
                        struct[2] = pdbLine.substring(31, 32);
                        struct[3] = pdbLine.substring(21, 25).trim(); // Residue
                        struct[4] = pdbLine.substring(33, 37).trim();
                        struct[5] = pdbLine.substring(38, 40).trim();
                        structs.add(struct);
                        break;
                    case SHEET:
                        struct = new String[6];
                        struct[0] = pdbLine.substring(0, 6).trim(); // SHEET
                        struct[1] = pdbLine.substring(21, 22); // Polymers
                        struct[2] = pdbLine.substring(32, 33);
                        struct[3] = pdbLine.substring(22, 26).trim(); // Residue
                        struct[4] = pdbLine.substring(33, 37).trim();
                        struct[5] = pdbLine.substring(38, 40).trim(); // Strand
                        structs.add(struct);
                        break;
                    case TURN:
                        struct = new String[6];
                        struct[0] = pdbLine.substring(0, 6).trim(); // TURN
                        struct[1] = pdbLine.substring(19, 20); // Polymers
                        struct[2] = pdbLine.substring(30, 31);
                        struct[3] = pdbLine.substring(20, 24).trim(); // Residue
                        struct[4] = pdbLine.substring(31, 35).trim();
                        structs.add(struct);
                        break;
                    default:
                        /**
                         * Do nothing for the other cards.
                         */
                        break;
                }
                if (bw != null) {
                    bw.write(pdbLine);
                    bw.newLine();
                }
                pdbLine = br.readLine();
            }
            if (bw != null) {
                bw.flush();
                bw.close();
            }
            br.close();
            xyzIndex--;
            if (logger.isLoggable(Level.INFO)) {
                logger.info(" Read " + xyzIndex + " atoms");
                StringBuffer altLocString = new StringBuffer(" Alternate locations [ ");
                for (Character c : altLocs) {
                    altLocString.append("(" + c + ") ");
                }
                altLocString.append("]");
                logger.info(altLocString.toString());
            }
            // Assign Secondary Structure Based on PDB Info
            for (String[] s : structs) {
                if (s[1].equalsIgnoreCase(" ")) {
                    s[1] = "Blank".intern();
                }
                Polymer p = molecularAssembly.getPolymer(s[1], false);
                if (p != null) {
                    for (int i = Integer.parseInt(s[3]); i <= Integer.parseInt(s[4]); i++) {
                        Residue r = p.getResidue(i);
                        if (r != null) {
                            r.setSSType(Residue.SSType.valueOf(s[0]));
                        }
                    }
                }
            }
            setFileRead(true);
        } catch (IOException e) {
            logger.exiting(PDBFilter.class.getName(), "readFile", e);
            return false;
        }
        // Now read in VRML data
        /*
        String name = molecularAssembly.getName() + ".wrl";
        File vrmlFile =
        new File(molecularAssembly.getFile().getParent() + File.separator + name);
        if (vrmlFile == null || !vrmlFile.exists()) {
        vrmlFile = new File(molecularAssembly.getFile().getParent() + File.separator + name);
        int retry = 0;
        boolean done = false;
        logger.info("Downloading VRML");
        try {
        URL copyURL = new URL(vrmlURL);
        logger.info(copyURL.toString());
        InputStreamReader ir = new InputStreamReader(copyURL.openStream());
        BufferedReader vbr = new BufferedReader(ir);
        BufferedWriter vbw = new BufferedWriter(new FileWriter(vrmlFile));
        while (!done && retry < 10) {
        while (!vbr.ready() && retry < 10) {
        synchronized (this) {
        logger.info("Waiting on Network");
        wait(1000);
        retry++;
        }
        }
        if (retry < 10) {
        String data = vbr.readLine();
        vbw.write(data);
        vbw.newLine();
        if (data.toUpperCase().indexOf("END OF OUTPUT") >= 0) {
        done = true;
        }
        }
        }
        vbr.close();
        vbw.flush();
        vbw.close();
        } catch (Exception e) {
        logger.severe(e.toString());
        done = false;
        } finally {
        if (retry == 100 || !done) {
        logger.warning("VRML download failed.");
        vrmlFile.delete();
        }
        }
        }
        try {
        logger.info("Loading VRML");
        if (vrmlFile.exists()) {
        molecularAssembly.setVRML(vrmlFile);
        } else if (vrmlURL != null) {
        URL url = new URL(vrmlURL);
        molecularAssembly.setVRML(url);
        }
        } catch (Exception e) {
        logger.warning("VRML Failed to Load\n" + e);
        return false;
        } */

        assignAtomTypes();
        renumberAtoms();

        return true;
    }

    public void renumberAtoms() {
        int index = 1;
        for (Atom a : molecularAssembly.getAtomList()) {
            a.xyzIndex = index++;
        }
        index--;
        if (logger.isLoggable(Level.INFO)) {
            logger.info(String.format(" Total number of atoms: %d", index));
        }
    }

    /**
     * Assign force field atoms types to common chemistries using
     * "biotype" records.
     */
    public void assignAtomTypes() {
        /**
         * Create a new List to store bonds determined based on PDB atom names.
         */
        bondList = new ArrayList<Bond>();

        /**
         * To Do: Look for cyclic peptides and disulfides.
         */
        String[] chainNames = molecularAssembly.getChainNames();

        /**
         * Loop over chains.
         */
        for (String chain : chainNames) {
            Polymer polymer = molecularAssembly.getPolymer(chain, false);
            ArrayList<Residue> residues = polymer.getResidues();
            int numberOfResidues = residues.size();

            /**
             * Check if all residues are known amino acids.
             */
            boolean isProtein = true;
            for (int residueNumber = 0; residueNumber < numberOfResidues; residueNumber++) {
                Residue residue = residues.get(residueNumber);
                String name = residue.getName().toUpperCase();
                AminoAcid3 aminoAcid = AminoAcid3.UNK;
                for (int a = 0; a < numberOfKnownAminoAcids; a++) {
                    AminoAcid3 amino = knownAminoAcids[a];
                    if (amino.toString().equalsIgnoreCase(name)) {
                        aminoAcid = amino;
                        break;
                    }
                }
                if (aminoAcid == AminoAcid3.UNK) {
                    isProtein = false;
                    break;
                }
            }

            /**
             * If all the residues in this chain have known amino acids names,
             * then attempt to assign atom types.
             */
            if (isProtein) {
                try {
                    assignAminoAcidAtomTypes(residues);
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(" Atom type assignment completed for amino acid chain " + chain + ".");
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
                NucleicAcid3 nucleicAcid = NucleicAcid3.UNK;
                for (int a = 0; a < numberOfKnownNucleicAcids; a++) {
                    NucleicAcid3 nucleic = knownNucleicAcids[a];
                    if (nucleic.toString().equalsIgnoreCase(name)) {
                        nucleicAcid = nucleic;
                        break;
                    }
                }
                if (nucleicAcid == NucleicAcid3.UNK) {
                    isNucleicAcid = false;
                    break;
                }
            }

            /**
             * If all the residues in this chain have known nucleic acids names,
             * then attempt to assign atom types.
             */
            if (isNucleicAcid) {
                try {
                    assignNucleicAcidAtomTypes(residues);
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info("Atom type assignment completed for nucleic acid chain " + chain + ".");
                    }
                } catch (MissingHeavyAtomException missingHeavyAtomException) {
                    logger.severe(missingHeavyAtomException.toString());
                } catch (MissingAtomTypeException missingAtomTypeException) {
                    logger.severe(missingAtomTypeException.toString());
                }
                continue;
            }
        }

        // Assign ion atom types.
        ArrayList<MSNode> ions = molecularAssembly.getIons();
        for (MSNode m : ions) {
            Molecule ion = (Molecule) m;
            String name = ion.getResidueName().toUpperCase();
            HETATOMS hetatm = HETATOMS.valueOf(name);
            Atom atom = ion.getAtomList().get(0);
            if (ion.getAtomList().size() != 1) {
                logger.severe("Check residue " + ion.toString() + " of chain " + ion.getPolymerName() + ".");
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
                        atom.setAtomType(findAtomType(2003));
                        break;
                    default:
                        logger.severe("Check residue " + ion.getResidueName() + " of chain " + ion.getPolymerName() + ".");
                }
            } catch (Exception e) {
                String message = "Error assigning atom types.";
                logger.log(Level.SEVERE, message, e);
            }
        }
        // Assign water atom types.
        ArrayList<MSNode> water = molecularAssembly.getWaters();
        for (MSNode m : water) {
            Molecule wat = (Molecule) m;
            try {
                Atom O = setHeavyAtom(wat, "O", null, 2001);
                Atom H1 = setHydrogenAtom(wat, "H1", O, 0.96e0, null, 109.5e0, null, 120.0e0, 0, 2002);
                Atom H2 = setHydrogenAtom(wat, "H2", O, 0.96e0, H1, 109.5e0, null, 120.0e0, 0, 2002);
            } catch (Exception e) {
                String message = "Error assigning atom types to a water.";
                logger.log(Level.SEVERE, message, e);
            }
        }
        // Assign small molecule atom types.
        ArrayList<MSNode> molecules = molecularAssembly.getMolecules();
        for (MSNode m : molecules) {
            molecularAssembly.deleteMolecule((Molecule) m);
        }
    }

    /**
     * Assign atom types for a nucleic acid polymer.
     *
     * @param residues
     * @throws ffx.potential.parsers.PDBFilter.MissingHeavyAtomException
     */
    private void assignNucleicAcidAtomTypes(ArrayList<Residue> residues)
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
            NucleicAcid3 nucleicAcid = NucleicAcid3.UNK;
            int naNumber = -1;
            for (int n = 0; n
                    < numberOfKnownNucleicAcids; n++) {
                NucleicAcid3 amino = knownNucleicAcids[n];
                if (amino.toString().equalsIgnoreCase(residueName)) {
                    nucleicAcid = amino;
                    naNumber = n;
                    break;
                }
            }
            /**
             * Check if the sugar is deoxyribose and change the residue
             * name if necessary.
             */
            boolean isDNA = false;
            Atom O2s = (Atom) residue.getAtomNode("O2*");
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
                 * The 5' O5* oxygen of the nucleic acid is generally
                 * terminated by
                 * 1.) A phosphate group PO3 (-3).
                 * 2.) A hydrogen.
                 *
                 * If the base has phosphate atom we will assume a PO3 group.
                 */
                P = (Atom) residue.getAtomNode("P");
                if (P != null) {
                    if (isDNA) {
                        P = setHeavyAtom(residue, "P", null, 1247);
                        setHeavyAtom(residue, "OP1", P, 1248);
                        setHeavyAtom(residue, "OP2", P, 1248);
                        setHeavyAtom(residue, "OP3", P, 1248);
                        O5s = setHeavyAtom(residue, "O5*", P, 1246);
                    } else {
                        P = setHeavyAtom(residue, "P", null, 1235);
                        setHeavyAtom(residue, "OP1", P, 1236);
                        setHeavyAtom(residue, "OP2", P, 1236);
                        setHeavyAtom(residue, "OP3", P, 1236);
                        O5s = setHeavyAtom(residue, "O5*", P, 1234);
                    }
                } else {
                    if (isDNA) {
                        O5s = setHeavyAtom(residue, "O5*", P, 1244);
                    } else {
                        O5s = setHeavyAtom(residue, "O5*", P, 1232);
                    }
                }
            } else {
                P = setHeavyAtom(residue, "P", pO3s, pTyp[naNumber]);
                setHeavyAtom(residue, "OP1", P, opTyp[naNumber]);
                setHeavyAtom(residue, "OP2", P, opTyp[naNumber]);
                O5s = setHeavyAtom(residue, "O5*", P, o5Typ[naNumber]);
            }
            /**
             * Build the ribose sugar atoms of the current base.
             */
            Atom C5s = setHeavyAtom(residue, "C5*", O5s, c5Typ[naNumber]);
            Atom C4s = setHeavyAtom(residue, "C4*", C5s, c4Typ[naNumber]);
            Atom O4s = setHeavyAtom(residue, "O4*", C4s, o4Typ[naNumber]);
            Atom C1s = setHeavyAtom(residue, "C1*", O4s, c1Typ[naNumber]);
            Atom C3s = setHeavyAtom(residue, "C3*", C4s, c3Typ[naNumber]);
            Atom C2s = setHeavyAtom(residue, "C2*", C3s, c2Typ[naNumber]);
            bond(C2s, C1s);
            Atom O3s = null;
            if (position == LAST_RESIDUE) {
                if (isDNA) {
                    O3s = setHeavyAtom(residue, "O3*", C3s, 1249);
                } else {
                    O3s = setHeavyAtom(residue, "O3*", C3s, 1237);
                }
            } else {
                O3s = setHeavyAtom(residue, "O3*", C3s, o3Typ[naNumber]);
            }
            if (!isDNA) {
                O2s = setHeavyAtom(residue, "O2*", C2s, o2Typ[naNumber]);
            }
            /**
             * Build the backbone hydrogen atoms.
             */
            if (position == FIRST_RESIDUE && P == null) {
                setHydrogenAtom(residue, "H5T", O5s, 1.00e0, C5s, 109.5e0, C4s, 180.0e0, 0, h5tTyp[naNumber]);
            }
            setHydrogenAtom(residue, "H5*1", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, 1, h51Typ[naNumber]);
            setHydrogenAtom(residue, "H5*2", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, -1, h52Typ[naNumber]);
            setHydrogenAtom(residue, "H4*", C4s, 1.09e0, C5s, 109.5e0, C3s, 109.5e0, -1, h4Typ[naNumber]);
            setHydrogenAtom(residue, "H3*", C3s, 1.09e0, C4s, 109.5e0, C2s, 109.5e0, -1, h3Typ[naNumber]);
            if (isDNA) {
                setHydrogenAtom(residue, "H2*1", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber]);
                setHydrogenAtom(residue, "H2*2", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, 1, h22Typ[naNumber]);
            } else {
                setHydrogenAtom(residue, "H2*", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber]);
                setHydrogenAtom(residue, "HO*", O2s, 1.00e0, C2s, 109.5e0, C3s, 180.0e0, 0, h22Typ[naNumber]);
            }
            setHydrogenAtom(residue, "H1*", C1s, 1.09e0, O4s, 109.5e0, C2s, 109.5e0, -1, h1Typ[naNumber]);
            if (position == LAST_RESIDUE) {
                setHydrogenAtom(residue, "H3T", O3s, 1.00e0, C3s, 109.5e0, C4s, 180.0e0, 0, h3tTyp[naNumber]);
            }
            /**
             * Build the nucleic acid base.
             */
            try {
                assignNucleicAcidBaseAtomTypes(nucleicAcid, residue, C1s);
            } catch (MissingHeavyAtomException missingHeavyAtomException) {
                logger.throwing(PDBFilter.class.getName(), "assignNucleicAcidAtomTypes", missingHeavyAtomException);
                throw missingHeavyAtomException;
            }

            /**
             * Do some checks on the current base to make sure all atoms
             * have been assigned an atom type.
             */
            ArrayList<Atom> atoms = residue.getAtomList();
            for (Atom atom : atoms) {
                AtomType atomType = atom.getAtomType();
                if (atomType == null) {
                    MissingAtomTypeException missingAtomTypeException = new MissingAtomTypeException(residue, atom);
                    logger.throwing(PDBFilter.class.getName(), "assignNucleicAcidAtomTypes", missingAtomTypeException);
                    throw missingAtomTypeException;
                }
                int numberOfBonds = atom.getNumBonds();
                if (numberOfBonds != atomType.valence) {
                    if (atom == O3s && numberOfBonds == atomType.valence - 1 && position != LAST_RESIDUE) {
                        continue;
                    }
                    System.out.println("An atom for residue " + residueName
                            + " has the wrong number of bonds.\n" + atom.toString());
                    System.out.println("Expected: " + atomType.valence + " Actual: " + numberOfBonds);
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
    private void assignNucleicAcidBaseAtomTypes(NucleicAcid3 nucleicAcid, Residue residue, Atom C1s)
            throws MissingHeavyAtomException {
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
                N9 = setHeavyAtom(residue, "N9", C1s, 1017);
                C8 = setHeavyAtom(residue, "C8", N9, 1021);
                N7 = setHeavyAtom(residue, "N7", C8, 1020);
                C5 = setHeavyAtom(residue, "C5", N7, 1019);
                C6 = setHeavyAtom(residue, "C6", C5, 1025);
                N6 = setHeavyAtom(residue, "N6", C6, 1027);
                N1 = setHeavyAtom(residue, "N1", C6, 1024);
                C2 = setHeavyAtom(residue, "C2", N1, 1023);
                N3 = setHeavyAtom(residue, "N3", C2, 1022);
                C4 = setHeavyAtom(residue, "C4", N3, 1018);
                bond(C4, C5);
                bond(C4, N9);
                setHydrogenAtom(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1030);
                setHydrogenAtom(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1028);
                setHydrogenAtom(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1029);
                setHydrogenAtom(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1026);
                break;
            case CYT:
                Atom O2;
                Atom N4;
                N1 = setHeavyAtom(residue, "N1", C1s, 1078);
                C2 = setHeavyAtom(residue, "C2", N1, 1079);
                O2 = setHeavyAtom(residue, "O2", C2, 1084);
                N3 = setHeavyAtom(residue, "N3", C2, 1080);
                C4 = setHeavyAtom(residue, "C4", N3, 1081);
                N4 = setHeavyAtom(residue, "N4", C4, 1085);
                C5 = setHeavyAtom(residue, "C5", C4, 1082);
                C6 = setHeavyAtom(residue, "C6", C5, 1083);
                bond(C6, N1);
                setHydrogenAtom(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1086);
                setHydrogenAtom(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1087);
                setHydrogenAtom(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1088);
                setHydrogenAtom(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1089);
                break;
            case GUA:
                Atom O6;
                Atom N2;
                N9 = setHeavyAtom(residue, "N9", C1s, 1047);
                C8 = setHeavyAtom(residue, "C8", N9, 1051);
                N7 = setHeavyAtom(residue, "N7", C8, 1050);
                C5 = setHeavyAtom(residue, "C5", N7, 1049);
                C6 = setHeavyAtom(residue, "C6", C5, 1055);
                O6 = setHeavyAtom(residue, "O6", C6, 1060);
                N1 = setHeavyAtom(residue, "N1", C6, 1054);
                C2 = setHeavyAtom(residue, "C2", N1, 1053);
                N2 = setHeavyAtom(residue, "N2", C2, 1057);
                N3 = setHeavyAtom(residue, "N3", C2, 1052);
                C4 = setHeavyAtom(residue, "C4", N3, 1048);
                bond(C4, C5);
                bond(C4, N9);
                setHydrogenAtom(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1061);
                setHydrogenAtom(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1056);
                setHydrogenAtom(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1058);
                setHydrogenAtom(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1059);
                break;
            case URI:
                Atom O4;
                N1 = setHeavyAtom(residue, "N1", C1s, 1106);
                C2 = setHeavyAtom(residue, "C2", N1, 1107);
                O2 = setHeavyAtom(residue, "O2", C2, 1112);
                N3 = setHeavyAtom(residue, "N3", C2, 1108);
                C4 = setHeavyAtom(residue, "C4", N3, 1109);
                O4 = setHeavyAtom(residue, "O4", C4, 1114);
                C5 = setHeavyAtom(residue, "C5", C4, 1110);
                C6 = setHeavyAtom(residue, "C6", C5, 1111);
                bond(C6, N1);
                setHydrogenAtom(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1113);
                setHydrogenAtom(residue, "H5", C5, 1.08e0, C4, 120.4e0, N3, 180.0e0, 0, 1115);
                setHydrogenAtom(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1116);
                break;
            case DAD:
                N9 = setHeavyAtom(residue, "N9", C1s, 1132);
                C8 = setHeavyAtom(residue, "C8", N9, 1136);
                N7 = setHeavyAtom(residue, "N7", C8, 1135);
                C5 = setHeavyAtom(residue, "C5", N7, 1134);
                C6 = setHeavyAtom(residue, "C6", C5, 1140);
                N6 = setHeavyAtom(residue, "N6", C6, 1142);
                N1 = setHeavyAtom(residue, "N1", C6, 1139);
                C2 = setHeavyAtom(residue, "C2", N1, 1138);
                N3 = setHeavyAtom(residue, "N3", C2, 1137);
                C4 = setHeavyAtom(residue, "C4", N3, 1133);
                bond(C4, C5);
                bond(C4, N9);
                setHydrogenAtom(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1145);
                setHydrogenAtom(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1143);
                setHydrogenAtom(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1144);
                setHydrogenAtom(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1141);
                break;
            case DCY:
                N1 = setHeavyAtom(residue, "N1", C1s, 1191);
                C2 = setHeavyAtom(residue, "C2", N1, 1192);
                O2 = setHeavyAtom(residue, "O2", C2, 1197);
                N3 = setHeavyAtom(residue, "N3", C2, 1193);
                C4 = setHeavyAtom(residue, "C4", N3, 1194);
                N4 = setHeavyAtom(residue, "N4", C4, 1198);
                C5 = setHeavyAtom(residue, "C5", C4, 1195);
                C6 = setHeavyAtom(residue, "C6", C5, 1196);
                bond(C6, N1);
                setHydrogenAtom(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1199);
                setHydrogenAtom(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1200);
                setHydrogenAtom(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1201);
                setHydrogenAtom(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1202);
                break;
            case DGU:
                N9 = setHeavyAtom(residue, "N9", C1s, 1161);
                C8 = setHeavyAtom(residue, "C8", N9, 1165);
                N7 = setHeavyAtom(residue, "N7", C8, 1164);
                C5 = setHeavyAtom(residue, "C5", N7, 1163);
                C6 = setHeavyAtom(residue, "C6", C5, 1169);
                O6 = setHeavyAtom(residue, "O6", C6, 1174);
                N1 = setHeavyAtom(residue, "N1", C6, 1168);
                C2 = setHeavyAtom(residue, "C2", N1, 1167);
                N2 = setHeavyAtom(residue, "N2", C2, 1171);
                N3 = setHeavyAtom(residue, "N3", C2, 1166);
                C4 = setHeavyAtom(residue, "C4", N3, 1162);
                bond(C4, C5);
                bond(C4, N9);
                setHydrogenAtom(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1175);
                setHydrogenAtom(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1170);
                setHydrogenAtom(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1172);
                setHydrogenAtom(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1173);
                break;
            case DTY:
                Atom C7;
                Atom H;
                N1 = setHeavyAtom(residue, "N1", C1s, 1218);
                C2 = setHeavyAtom(residue, "C2", N1, 1219);
                O2 = setHeavyAtom(residue, "O2", C2, 1224);
                N3 = setHeavyAtom(residue, "N3", O2, 1220);
                C4 = setHeavyAtom(residue, "C4", C2, 1221);
                O4 = setHeavyAtom(residue, "O4", C4, 1226);
                C5 = setHeavyAtom(residue, "C5", C4, 1222);
                C7 = setHeavyAtom(residue, "C7", C5, 1227);
                C6 = setHeavyAtom(residue, "C6", C5, 1223);
                bond(C6, N1);
                setHydrogenAtom(residue, "H3", N3, 1.00e0, C2, 116.8e0, N1, 180.0e0, 0, 1225);
                H = setHydrogenAtom(residue, "H71", C7, 1.09e0, C5, 109.5e0, C4, 0.0e0, 0, 1228);
                setHydrogenAtom(residue, "H72", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, 1, 1228);
                setHydrogenAtom(residue, "H73", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, -1, 1228);
                setHydrogenAtom(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1229);
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
    private void assignAminoAcidAtomTypes(ArrayList<Residue> residues)
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
                 * If the lest residue only contains a nitrogen turn it into
                 * an NH2 group.
                 */
                Atom N = (Atom) residue.getAtomNode("N");
                if (residue.getAtomNodeList().size() == 1 && N != null) {
                    residueName = "NH2".intern();
                    residue.setName(residueName);
                }
            }
            AminoAcid3 aminoAcid = AminoAcid3.UNK;
            int aminoAcidNumber = -1;
            for (int a = 0; a
                    < numberOfKnownAminoAcids; a++) {
                AminoAcid3 amino = knownAminoAcids[a];
                if (amino.toString().equalsIgnoreCase(residueName)) {
                    aminoAcid = amino;
                    aminoAcidNumber = a;
                    break;
                }
            }
            /**
             * Check for missing heavy atoms.
             *
             * This check ignores special terminating groups like
             * FOR, NH2, etc.
             */
            int expected = aminoAcidHeavyAtoms[aminoAcidNumber];
            if (aminoAcid != AminoAcid3.GLY && expected >= 4) {
                int actual = 0;
                ArrayList<Atom> atoms = residue.getAtomList();
                for (Atom atom : atoms) {
                    String label = atom.getName().toUpperCase();
                    if (!(label.equalsIgnoreCase("OXT") || label.equalsIgnoreCase("OT2"))) {
                        if (!label.startsWith("H")) {
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
                    if (aminoAcid == AminoAcid3.ALA && actual == 4) {
                        residueName = "GLY".intern();
                        residue.setName(residueName);
                    } else if (actual == 5) {
                        residueName = "ALA".intern();
                        residue.setName(residueName);
                    }
                }
            }
            aminoAcid = AminoAcid3.UNK;
            aminoAcidNumber = -1;
            for (int a = 0; a
                    < numberOfKnownAminoAcids; a++) {
                AminoAcid3 amino = knownAminoAcids[a];
                if (amino.toString().equalsIgnoreCase(residueName)) {
                    aminoAcid = amino;
                    aminoAcidNumber = a;
                    break;
                }
            }

            /**
             * Backbone heavy atoms.
             */
            Atom N = (Atom) residue.getAtomNode("N");
            N.setAtomType(findAtomType(nType[j][aminoAcidNumber]));
            if (position != FIRST_RESIDUE) {
                bond(pC, N);
            }
            Atom CA = null;
            Atom C = null;
            Atom O = null;
            if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NH2)) {
                CA = setHeavyAtom(residue, "CA", N, caType[j][aminoAcidNumber]);
                if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NME)) {
                    C = setHeavyAtom(residue, "C", CA, cType[j][aminoAcidNumber]);
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
                    bond(C, O);
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
                            setHydrogenAtom(residue, "H2", N, 1.01e0, CA, 109.5e0, C, 0.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H3", N, 1.01e0, CA, 109.5e0, C, -120.0e0, 0, atomType);
                            break;
                        case PCA:
                            setHydrogenAtom(residue, "H", N, 1.01e0, CA, 109.5e0, C, -60.0e0, 0, atomType);
                            break;
                        default:
                            setHydrogenAtom(residue, "H1", N, 1.01e0, CA, 109.5e0, C, 180.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H2", N, 1.01e0, CA, 109.5e0, C, 60.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H3", N, 1.01e0, CA, 109.5e0, C, -60.0e0, 0, atomType);
                    }
                    break;
                case LAST_RESIDUE:
                    switch (aminoAcid) {
                        case NH2:
                            setHydrogenAtom(residue, "H1", N, 1.01e0, pC, 120.9e0, pCA, 0.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H2", N, 1.01e0, pC, 120.3e0, pCA, 180.0e0, 0, atomType);
                            break;
                        case NME:
                            setHydrogenAtom(residue, "H", N, 1.01e0, C, 119.0e0, CA, 119.0e0, 1, atomType);
                            break;
                        default:
                            setHydrogenAtom(residue, "H", N, 1.01e0, pC, 119.0e0, CA, 119.0e0, 1, atomType);
                    }
                    break;
                default:
                    // Mid-chain nitrogen hydrogen.
                    setHydrogenAtom(residue, "H", N, 1.01e0, pC, 119.0e0, CA, 119.0e0, 1, atomType);
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
                            setHydrogenAtom(residue, "H", C, 1.12e0, O, 0.0e0, null, 0.0e0, 0, atomType);
                            break;
                        case ACE:
                            setHydrogenAtom(residue, "H1", CA, 1.10e0, C, 109.5e0, O, 180.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H2", CA, 1.10e0, C, 109.5e0, O, 60.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H3", CA, 1.10e0, C, 109.5e0, O, -60.0e0, 0, atomType);
                            break;
                        default:
                            setHydrogenAtom(residue, haName, CA, 1.10e0, N, 109.5e0, C, 109.5e0, -1, atomType);
                            break;
                    }
                    break;
                case LAST_RESIDUE:
                    switch (aminoAcid) {
                        case NME:
                            setHydrogenAtom(residue, "H1", CA, 1.10e0, N, 109.5e0, pC, 180.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H2", CA, 1.10e0, N, 109.5e0, pC, 60.0e0, 0, atomType);
                            setHydrogenAtom(residue, "H3", CA, 1.10e0, N, 109.5e0, pC, -60.0e0, 0, atomType);
                            break;
                        default:
                            setHydrogenAtom(residue, haName, CA, 1.10e0, N, 109.5e0, C, 109.5e0, -1, atomType);
                    }
                    break;
                default:
                    setHydrogenAtom(residue, haName, CA, 1.10e0, N, 109.5e0, C, 109.0e0, -1, atomType);
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
                    OXT = new Atom("OXT", atomType, new double[3]);
                    residue.addMSNode(OXT);
                    OXT.setOccupancy(C.getOccupancy());
                    OXT.setTempFactor(C.getTempFactor());
                    intxyz(OXT, C, 1.25e0, CA, 117.0e0, O, 126.0, 1);
                } else {
                    OXT.setAtomType(atomType);
                }
                bond(C, OXT);
            }
            /**
             * Do some checks on the current residue to make sure all atoms
             * have been assigned an atom type.
             */
            ArrayList<Atom> atoms = residue.getAtomList();
            for (Atom atom : atoms) {
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
                    System.out.println("An atom for residue " + residueName
                            + " has the wrong number of bonds.\n" + atom.toString());
                    System.out.println("Expected: " + atomType.valence + " Actual: " + numberOfBonds);
                }
            }
            /**
             * Remember the current C-alpha and carboxyl C atoms for use
             * with the next residue.
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
        switch (aminoAcid) {
            case GLY:
                switch (position) {
                    case FIRST_RESIDUE:
                        setHydrogenAtom(residue, "HA3", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 355);
                        break;
                    case LAST_RESIDUE:
                        setHydrogenAtom(residue, "HA3", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 506);
                        break;
                    default:
                        setHydrogenAtom(residue, "HA3", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 6);
                }
                break;
            case ALA:
                Atom CB = setHeavyAtom(residue, "CB", CA, 13);
                setHydrogenAtom(residue, "HB1", CB, 1.10e0, CA, 110.2e0, N, 180.0e0, 0, 14);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 110.2e0, N, 60.0e0, 0, 14);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 110.2e0, N, -60.0e0, 0, 14);
                break;
            case VAL:
                CB = setHeavyAtom(residue, "CB", CA, 21);
                Atom CG1 = setHeavyAtom(residue, "CG1", CB, 23);
                Atom CG2 = setHeavyAtom(residue, "CG2", CB, 25);
                setHydrogenAtom(residue, "HB", CB, 1.10e0, CA, 107.0e0, CG1, 108.2e0, 1, 22);
                setHydrogenAtom(residue, "HG11", CG1, 1.10e0, CB, 111.6e0, CA, 180.0e0, 0, 24);
                setHydrogenAtom(residue, "HG12", CG1, 1.10e0, CB, 111.6e0, CA, 60.0e0, 0, 24);
                setHydrogenAtom(residue, "HG13", CG1, 1.10e0, CB, 111.6e0, CA, -60.0e0, 0, 24);
                setHydrogenAtom(residue, "HG21", CG2, 1.10e0, CB, 111.6e0, CA, 180.0e0, 0, 26);
                setHydrogenAtom(residue, "HG22", CG2, 1.10e0, CB, 111.6e0, CA, 60.0e0, 0, 26);
                setHydrogenAtom(residue, "HG23", CG2, 1.10e0, CB, 111.6e0, CA, -60.0e0, 0, 26);
                break;
            case LEU:
                CB = setHeavyAtom(residue, "CB", CA, 33);
                Atom CG = setHeavyAtom(residue, "CG", CB, 35);
                Atom CD1 = setHeavyAtom(residue, "CD1", CG, 37);
                Atom CD2 = setHeavyAtom(residue, "CD2", CG, 39);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 34);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 34);
                setHydrogenAtom(residue, "HG", CG, 1.10e0, CB, 107.0e0, CD1, 108.2e0, 1, 36);
                setHydrogenAtom(residue, "HD11", CD1, 1.10e0, CG, 111.6e0, CB, 180.0e0, 0, 38);
                setHydrogenAtom(residue, "HD12", CD1, 1.10e0, CG, 111.6e0, CB, 60.0e0, 0, 38);
                setHydrogenAtom(residue, "HD13", CD1, 1.10e0, CG, 111.6e0, CB, -60.0e0, 0, 38);
                setHydrogenAtom(residue, "HD21", CD2, 1.10e0, CG, 111.6e0, CB, 180.0e0, 0, 40);
                setHydrogenAtom(residue, "HD22", CD2, 1.10e0, CG, 111.6e0, CB, 60.0e0, 0, 40);
                setHydrogenAtom(residue, "HD23", CD2, 1.10e0, CG, 111.6e0, CB, -60.0e0, 0, 40);
                break;
            case ILE:
                CB = setHeavyAtom(residue, "CB", CA, 47);
                CG1 = setHeavyAtom(residue, "CG1", CB, 49);
                CG2 = setHeavyAtom(residue, "CG2", CB, 51);
                try {
                    CD1 = setHeavyAtom(residue, "CD1", CG1, 53);
                } catch (MissingHeavyAtomException missingHeavyAtomException) {
                    CD1 = setHeavyAtom(residue, "CD", CG1, 53);
                }
                setHydrogenAtom(residue, "HB", CB, 1.10e0, CA, 107.0e0, CG1, 108.2e0, -1, 48);
                setHydrogenAtom(residue, "HG12", CG1, 1.10e0, CB, 109.5e0, CD1, 109.5e0, 1, 50);
                setHydrogenAtom(residue, "HG13", CG1, 1.10e0, CB, 109.5e0, CD1, 109.5e0, -1, 50);
                setHydrogenAtom(residue, "HG21", CG2, 1.10e0, CB, 111.6e0, CA, 180.0e0, 0, 52);
                setHydrogenAtom(residue, "HG22", CG2, 1.10e0, CB, 111.6e0, CA, 60.0e0, 0, 52);
                setHydrogenAtom(residue, "HG23", CG2, 1.10e0, CB, 111.6e0, CA, -60.0e0, 0, 52);
                setHydrogenAtom(residue, "HD11", CD1, 1.10e0, CG1, 111.6e0, CB, 180.0e0, 0, 54);
                setHydrogenAtom(residue, "HD12", CD1, 1.10e0, CG1, 111.6e0, CB, 60.0e0, 0, 54);
                setHydrogenAtom(residue, "HD13", CD1, 1.10e0, CG1, 111.6e0, CB, -60.0e0, 0, 54);
                break;
            case SER:
                CB = setHeavyAtom(residue, "CB", CA, 61);
                Atom OG = setHeavyAtom(residue, "OG", CB, 63);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 109.2e0, OG, 109.5e0, 1, 62);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 109.2e0, OG, 109.5e0, -1, 62);
                setHydrogenAtom(residue, "HG", OG, 0.94e0, CB, 106.9e0, CA, 180.0e0, 0, 64);
                break;
            case THR:
                CB = setHeavyAtom(residue, "CB", CA, 71);
                Atom OG1 = setHeavyAtom(residue, "OG1", CB, 73);
                CG2 = setHeavyAtom(residue, "CG2", CB, 75);
                setHydrogenAtom(residue, "HB", CB, 1.10e0, CA, 107.0e0, OG1, 108.2e0, -1, 72);
                setHydrogenAtom(residue, "HG1", OG1, 0.94e0, CB, 106.9e0, CA, 180.0e0, 0, 74);
                setHydrogenAtom(residue, "HG21", CG2, 1.10e0, CB, 111.6e0, CA, 180.0e0, 0, 76);
                setHydrogenAtom(residue, "HG22", CG2, 1.10e0, CB, 111.6e0, CA, 60.0e0, 0, 76);
                setHydrogenAtom(residue, "HG23", CG2, 1.10e0, CB, 111.6e0, CA, -60.0e0, 0, 76);
                break;
            case CYS:
                CB = setHeavyAtom(residue, "CB", CA, 83);
                Atom SG = setHeavyAtom(residue, "SG", CB, 85);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 109.5e0, SG, 107.5e0, 1, 84);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 109.5e0, SG, 107.5e0, -1, 84);
                setHydrogenAtom(residue, "HG", SG, 1.34e0, CB, 96.0e0, CA, 180.0e0, 0, 86);
                break;
            case CYX:
                CB = setHeavyAtom(residue, "CB", CA, 93);
                SG = setHeavyAtom(residue, "SG", CB, 95);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 109.5e0, SG, 107.5e0, 1, 94);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 109.5e0, SG, 107.5e0, -1, 94);
                break;
            case PRO:
                CB = setHeavyAtom(residue, "CB", CA, 101);
                CG = setHeavyAtom(residue, "CG", CB, 103);
                Atom CD = null;
                if (position == FIRST_RESIDUE) {
                    CD = setHeavyAtom(residue, "CD", CG, 410);
                } else {
                    CD = setHeavyAtom(residue, "CD", CG, 105);
                }
                bond(CD, N);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 111.2e0, CG, 111.2e0, 1, 102);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 111.2e0, CG, 111.2e0, -1, 102);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 111.2e0, CD, 111.2e0, 1, 104);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 111.2e0, CD, 111.2e0, -1, 104);
                if (position == FIRST_RESIDUE) {
                    setHydrogenAtom(residue, "HD2", CD, 1.10e0, CG, 111.2e0, N, 111.2e0, 1, 411);
                    setHydrogenAtom(residue, "HD3", CD, 1.10e0, CG, 111.2e0, N, 111.2e0, -1, 411);
                } else {
                    setHydrogenAtom(residue, "HD2", CD, 1.10e0, CG, 111.2e0, N, 111.2e0, 1, 106);
                    setHydrogenAtom(residue, "HD3", CD, 1.10e0, CG, 111.2e0, N, 111.2e0, -1, 106);
                }
                break;
            case PHE:
                CB = setHeavyAtom(residue, "CB", CA, 113);
                CG = setHeavyAtom(residue, "CG", CB, 115);
                CD1 = setHeavyAtom(residue, "CD1", CG, 116);
                CD2 = setHeavyAtom(residue, "CD2", CG, 116);
                Atom CE1 = setHeavyAtom(residue, "CE1", CD1, 118);
                Atom CE2 = setHeavyAtom(residue, "CE2", CD2, 118);
                Atom CZ = setHeavyAtom(residue, "CZ", CE1, 120);
                bond(CE2, CZ);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 114);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 114);
                setHydrogenAtom(residue, "HD1", CD1, 1.10e0, CG, 120.0e0, CB, 0.0e0, 0, 117);
                setHydrogenAtom(residue, "HD2", CD2, 1.10e0, CG, 120.0e0, CB, 0.0e0, 0, 117);
                setHydrogenAtom(residue, "HE1", CE1, 1.10e0, CD1, 120.0e0, CG, 180.0e0, 0, 119);
                setHydrogenAtom(residue, "HE2", CE2, 1.10e0, CD2, 120.0e0, CG, 180.0e0, 0, 119);
                setHydrogenAtom(residue, "HZ", CZ, 1.10e0, CE2, 120.0e0, CD2, 180.0e0, 0, 121);
                break;
            case TYR:
                CB = setHeavyAtom(residue, "CB", CA, 128);
                CG = setHeavyAtom(residue, "CG", CB, 130);
                CD1 = setHeavyAtom(residue, "CD1", CG, 131);
                CD2 = setHeavyAtom(residue, "CD2", CG, 131);
                CE1 = setHeavyAtom(residue, "CE1", CD1, 133);
                CE2 = setHeavyAtom(residue, "CE2", CD2, 133);
                CZ = setHeavyAtom(residue, "CZ", CE1, 135);
                bond(CE2, CZ);
                Atom OH = setHeavyAtom(residue, "OH", CZ, 136);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 129);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 129);
                setHydrogenAtom(residue, "HD1", CD1, 1.10e0, CG, 120.0e0, CB, 0.0e0, 0, 132);
                setHydrogenAtom(residue, "HD2", CD2, 1.10e0, CG, 120.0e0, CB, 0.0e0, 0, 132);
                setHydrogenAtom(residue, "HE1", CE1, 1.10e0, CD1, 120.0e0, CG, 180.0e0, 0, 134);
                setHydrogenAtom(residue, "HE2", CE2, 1.10e0, CD2, 120.0e0, CG, 180.0e0, 0, 134);
                setHydrogenAtom(residue, "HH", OH, 0.97e0, CZ, 108.0e0, CE2, 0.0e0, 0, 137);
                break;
            case TRP:
                CB = setHeavyAtom(residue, "CB", CA, 144);
                CG = setHeavyAtom(residue, "CG", CB, 146);
                CD1 = setHeavyAtom(residue, "CD1", CG, 147);
                CD2 = setHeavyAtom(residue, "CD2", CG, 149);
                Atom NE1 = setHeavyAtom(residue, "NE1", CD1, 150);
                CE2 = setHeavyAtom(residue, "CE2", NE1, 152);
                bond(CE2, CD2);
                Atom CE3 = setHeavyAtom(residue, "CE3", CD2, 153);
                Atom CZ2 = setHeavyAtom(residue, "CZ2", CE2, 155);
                Atom CZ3 = setHeavyAtom(residue, "CZ3", CE3, 157);
                Atom CH2 = setHeavyAtom(residue, "CH2", CZ3, 159);
                bond(CH2, CZ2);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 145);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 145);
                setHydrogenAtom(residue, "HD1", CD1, 1.09e0, CG, 126.0e0, CB, 0.0e0, 0, 148);
                setHydrogenAtom(residue, "HE1", NE1, 1.01e0, CD1, 126.3e0, CG, 180.0e0, 0, 151);
                setHydrogenAtom(residue, "HE3", CE3, 1.09e0, CZ3, 120.0e0, CH2, 180.0e0, 0, 154);
                setHydrogenAtom(residue, "HZ2", CZ2, 1.09e0, CH2, 120.0e0, CZ3, 180.0e0, 0, 156);
                setHydrogenAtom(residue, "HZ3", CZ3, 1.09e0, CH2, 120.0e0, CZ2, 180.0e0, 0, 158);
                setHydrogenAtom(residue, "HH2", CH2, 1.09e0, CZ3, 120.0e0, CE3, 180.0e0, 0, 160);
                break;
            case HIS:
                CB = setHeavyAtom(residue, "CB", CA, 167);
                CG = setHeavyAtom(residue, "CG", CB, 169);
                Atom ND1 = setHeavyAtom(residue, "ND1", CG, 170);
                CD2 = setHeavyAtom(residue, "CD2", CG, 172);
                CE1 = setHeavyAtom(residue, "CE1", ND1, 174);
                Atom NE2 = setHeavyAtom(residue, "NE2", CE1, 176);
                bond(NE2, CD2);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 168);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 168);
                setHydrogenAtom(residue, "HD1", ND1, 1.02e0, CE1, 126.0e0, NE2, 180.0e0, 0, 171);
                setHydrogenAtom(residue, "HD2", CD2, 1.09e0, NE2, 126.0e0, CE1, 180.0e0, 0, 173);
                setHydrogenAtom(residue, "HE1", CE1, 1.09e0, NE2, 126.0e0, CD2, 180.0e0, 0, 175);
                setHydrogenAtom(residue, "HE2", NE2, 1.02e0, CE1, 126.0e0, ND1, 180.0e0, 0, 177);
                break;
            case HID:
                CB = setHeavyAtom(residue, "CB", CA, 184);
                CG = setHeavyAtom(residue, "CG", CB, 186);
                ND1 = setHeavyAtom(residue, "ND1", CG, 187);
                CD2 = setHeavyAtom(residue, "CD2", CG, 189);
                CE1 = setHeavyAtom(residue, "CE1", ND1, 191);
                NE2 = setHeavyAtom(residue, "NE2", CE1, 193);
                bond(NE2, CD2);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 185);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 185);
                setHydrogenAtom(residue, "HD1", ND1, 1.02e0, CE1, 126.0e0, NE2, 180.0e0, 0, 188);
                setHydrogenAtom(residue, "HD2", CD2, 1.09e0, NE2, 126.0e0, CE1, 180.0e0, 0, 190);
                setHydrogenAtom(residue, "HE1", CE1, 1.09e0, NE2, 126.0e0, CD2, 180.0e0, 0, 192);
                break;
            case HIE:
                CB = setHeavyAtom(residue, "CB", CA, 200);
                CG = setHeavyAtom(residue, "CG", CB, 202);
                ND1 = setHeavyAtom(residue, "ND1", CG, 203);
                CD2 = setHeavyAtom(residue, "CD2", CG, 204);
                CE1 = setHeavyAtom(residue, "CE1", ND1, 206);
                NE2 = setHeavyAtom(residue, "NE2", CE1, 208);
                bond(NE2, CD2);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 201);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 201);
                setHydrogenAtom(residue, "HD2", CD2, 1.09e0, NE2, 126.0e0, CE1, 180.0e0, 0, 205);
                setHydrogenAtom(residue, "HE1", CE1, 1.09e0, NE2, 126.0e0, CD2, 180.0e0, 0, 207);
                setHydrogenAtom(residue, "HE2", NE2, 1.02e0, CE1, 126.0e0, ND1, 180.0e0, 0, 209);
                break;
            case ASP:
                CB = setHeavyAtom(residue, "CB", CA, 216);
                CG = setHeavyAtom(residue, "CG", CB, 218);
                Atom OD1 = setHeavyAtom(residue, "OD1", CG, 219);
                Atom OD2 = setHeavyAtom(residue, "OD2", CG, 219);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 217);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 217);
                break;
            case ASN:
                CB = setHeavyAtom(residue, "CB", CA, 226);
                CG = setHeavyAtom(residue, "CG", CB, 228);
                OD1 = setHeavyAtom(residue, "OD1", CG, 229);
                Atom ND2 = setHeavyAtom(residue, "ND2", CG, 230);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 227);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 227);
                setHydrogenAtom(residue, "HD21", ND2, 1.01e0, CG, 120.9e0, CB, 0.0e0, 0, 231);
                setHydrogenAtom(residue, "HD22", ND2, 1.01e0, CG, 120.3e0, CB, 180.0e0, 0, 231);
                break;
            case GLU:
                CB = setHeavyAtom(residue, "CB", CA, 238);
                CG = setHeavyAtom(residue, "CG", CB, 240);
                CD = setHeavyAtom(residue, "CD", CG, 242);
                Atom OE1 = setHeavyAtom(residue, "OE1", CD, 243);
                Atom OE2 = setHeavyAtom(residue, "OE2", CD, 243);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 239);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 239);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, 1, 241);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, -1, 241);
                break;
            case GLN:
                CB = setHeavyAtom(residue, "CB", CA, 250);
                CG = setHeavyAtom(residue, "CG", CB, 252);
                CD = setHeavyAtom(residue, "CD", CG, 254);
                OE1 = setHeavyAtom(residue, "OE1", CD, 255);
                NE2 = setHeavyAtom(residue, "NE2", CD, 256);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 251);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 251);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, 1, 253);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, -1, 253);
                setHydrogenAtom(residue, "HE21", NE2, 1.01e0, CD, 120.9e0, CG, 0.0e0, 0, 257);
                setHydrogenAtom(residue, "HE22", NE2, 1.01e0, CD, 120.3e0, CG, 180.0e0, 0, 257);
                break;
            case MET:
                CB = setHeavyAtom(residue, "CB", CA, 264);
                CG = setHeavyAtom(residue, "CG", CB, 266);
                Atom SD = setHeavyAtom(residue, "SD", CG, 268);
                Atom CE = setHeavyAtom(residue, "CE", SD, 269);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 265);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 265);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, SD, 109.5e0, 1, 267);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, SD, 109.5e0, -1, 267);
                setHydrogenAtom(residue, "HE1", CE, 1.10e0, SD, 110.2e0, CG, 180.0e0, 0, 270);
                setHydrogenAtom(residue, "HE2", CE, 1.10e0, SD, 110.2e0, CG, 60.0e0, 0, 270);
                setHydrogenAtom(residue, "HE3", CE, 1.10e0, SD, 110.2e0, CG, -60.0e0, 0, 270);
                break;
            case LYS:
                CB = setHeavyAtom(residue, "CB", CA, 277);
                CG = setHeavyAtom(residue, "CG", CB, 279);
                CD = setHeavyAtom(residue, "CD", CG, 281);
                CE = setHeavyAtom(residue, "CE", CD, 283);
                Atom NZ = setHeavyAtom(residue, "NZ", CE, 285);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 278);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 278);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, 1, 280);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, -1, 280);
                setHydrogenAtom(residue, "HD2", CD, 1.10e0, CG, 109.5e0, CE, 109.5e0, 1, 282);
                setHydrogenAtom(residue, "HD3", CD, 1.10e0, CG, 109.5e0, CE, 109.5e0, -1, 282);
                setHydrogenAtom(residue, "HE2", CE, 1.10e0, CD, 110.9e0, NZ, 107.3e0, 1, 284);
                setHydrogenAtom(residue, "HE3", CE, 1.10e0, CD, 110.9e0, NZ, 107.3e0, -1, 284);
                setHydrogenAtom(residue, "HZ1", NZ, 1.04e0, CE, 110.5e0, CD, 180.0e0, 0, 286);
                setHydrogenAtom(residue, "HZ2", NZ, 1.04e0, CE, 110.5e0, CD, 60.0e0, 0, 286);
                setHydrogenAtom(residue, "HZ3", NZ, 1.04e0, CE, 110.5e0, CD, -60.0e0, 0, 286);
                break;
            case ARG:
                CB = setHeavyAtom(residue, "CB", CA, 293);
                CG = setHeavyAtom(residue, "CG", CB, 295);
                CD = setHeavyAtom(residue, "CD", CG, 297);
                Atom NE = setHeavyAtom(residue, "NE", CD, 299);
                CZ = setHeavyAtom(residue, "CZ", NE, 301);
                Atom NH1 = setHeavyAtom(residue, "NH1", CZ, 302);
                Atom NH2 = setHeavyAtom(residue, "NH2", CZ, 302);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 294);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 294);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, 1, 296);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, -1, 296);
                setHydrogenAtom(residue, "HD2", CD, 1.10e0, CG, 109.5e0, NE, 109.5e0, 1, 298);
                setHydrogenAtom(residue, "HD3", CD, 1.10e0, CG, 109.5e0, NE, 109.5e0, -1, 298);
                setHydrogenAtom(residue, "HE", NE, 1.01e0, CD, 118.5e0, CZ, 120.0e0, 1, 300);
                setHydrogenAtom(residue, "HH11", NH1, 1.01e0, CZ, 122.5e0, NE, 0.0e0, 0, 303);
                setHydrogenAtom(residue, "HH12", NH1, 1.01e0, CZ, 118.5e0, NE, 180.0e0, 0, 303);
                setHydrogenAtom(residue, "HH21", NH2, 1.01e0, CZ, 122.5e0, NE, 0.0e0, 0, 303);
                setHydrogenAtom(residue, "HH22", NH2, 1.01e0, CZ, 118.5e0, NE, 180.0e0, 0, 303);
                break;
            case ORN:
                CB = setHeavyAtom(residue, "CB", CA, 310);
                CG = setHeavyAtom(residue, "CG", CB, 312);
                CD = setHeavyAtom(residue, "CD", CG, 314);
                NE = setHeavyAtom(residue, "NE", CD, 316);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, 1, 311);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 107.9e0, CG, 110.0e0, -1, 311);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, 1, 313);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 109.5e0, CD, 109.5e0, -1, 313);
                setHydrogenAtom(residue, "HD2", CD, 1.10e0, CG, 109.5e0, NE, 109.5e0, 1, 315);
                setHydrogenAtom(residue, "HD3", CD, 1.10e0, CG, 109.5e0, NE, 109.5e0, -1, 315);
                setHydrogenAtom(residue, "HE1", NE, 1.04e0, CD, 110.5e0, CG, 180.0e0, 0, 317);
                setHydrogenAtom(residue, "HE2", NE, 1.04e0, CD, 110.5e0, CG, 60.0e0, 0, 317);
                setHydrogenAtom(residue, "HE3", NE, 1.04e0, CD, 110.5e0, CG, -60.0e0, 0, 317);
                break;
            case AIB:
                Atom CB1 = setHeavyAtom(residue, "CB1", CA, 323);
                Atom CB2 = setHeavyAtom(residue, "CB1", CA, 323);
                setHydrogenAtom(residue, "HB11", CB1, 1.10e0, CA, 110.2e0, N, 180.0e0, 0, 324);
                setHydrogenAtom(residue, "HB12", CB1, 1.10e0, CA, 110.2e0, N, 60.0e0, 0, 324);
                setHydrogenAtom(residue, "HB13", CB1, 1.10e0, CA, 110.2e0, N, -60.0e0, 0, 324);
                setHydrogenAtom(residue, "HB21", CB2, 1.10e0, CA, 110.2e0, N, 180.0e0, 0, 324);
                setHydrogenAtom(residue, "HG22", CB2, 1.10e0, CA, 110.2e0, N, 60.0e0, 0, 324);
                setHydrogenAtom(residue, "HG23", CB2, 1.10e0, CA, 110.2e0, N, -60.0e0, 0, 324);
                break;
            case PCA:
                CB = setHeavyAtom(residue, "CB", CA, 331);
                CG = setHeavyAtom(residue, "CG", CB, 333);
                CD = setHeavyAtom(residue, "CD", CG, 335);
                Atom OE = setHeavyAtom(residue, "OE", CD, 336);
                setHydrogenAtom(residue, "HB2", CB, 1.10e0, CA, 111.2e0, CG, 111.2e0, 1, 332);
                setHydrogenAtom(residue, "HB3", CB, 1.10e0, CA, 111.2e0, CG, 111.2e0, -1, 332);
                setHydrogenAtom(residue, "HG2", CG, 1.10e0, CB, 111.2e0, CD, 111.2e0, 1, 334);
                setHydrogenAtom(residue, "HG3", CG, 1.10e0, CB, 111.2e0, CD, 111.2e0, -1, 334);
                break;
            case UNK:
                switch (position) {
                    case FIRST_RESIDUE:
                        setHydrogenAtom(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 355);
                        break;
                    case LAST_RESIDUE:
                        setHydrogenAtom(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 506);
                        break;
                    default:
                        setHydrogenAtom(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 6);
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
            StringBuffer sb = new StringBuffer(super.toString());
            sb.append("Atom " + atom.toString());
            if (residue != null) {
                sb.append("\nof residue " + residue.toString());
            }
            sb.append("\ncould not be assigned an atom type.\n");
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
            StringBuffer sb = new StringBuffer(super.toString());
            if (atomType != null) {
                sb.append("\nAn atom of type:\n" + atomType.toString() + "\n");
            } else {
                sb.append("\nAtom " + atomName + " ");
            }
            sb.append("was not found");
            if (bondedTo != null) {
                sb.append(" bonded to atom " + bondedTo);
            }
            sb.append(".\n");
            return sb.toString();
        }
    }

    private Atom setHeavyAtom(MSGroup residue, String atomName, Atom bondedTo, int key)
            throws MissingHeavyAtomException {
        Atom atom = (Atom) residue.getAtomNode(atomName);
        AtomType atomType = findAtomType(key);
        if (atom == null) {
            MissingHeavyAtomException missingHeavyAtom = new MissingHeavyAtomException(atomName, atomType, bondedTo);
            throw missingHeavyAtom;
        }
        atom.setAtomType(atomType);
        if (bondedTo != null) {
            bond(atom, bondedTo);
        }
        return atom;
    }

    private Atom setHydrogenAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, int lookUp) {
        AtomType atomType = findAtomType(lookUp);
        return setHydrogenAtom(residue, atomName, ia, bond, ib, angle1, ic, angle2, chiral, atomType);
    }

    private Atom setHydrogenAtom(MSGroup residue, String atomName, Atom ia, double bond, Atom ib, double angle1,
            Atom ic, double angle2, int chiral, AtomType atomType) {
        if (atomType == null) {
            return null;
        }
        Atom atom = (Atom) residue.getAtomNode(atomName);
        if (atom == null) {
            atom = new Atom(atomName, atomType, new double[3]);
            residue.addMSNode(atom);
            atom.setOccupancy(ia.getOccupancy());
            atom.setTempFactor(ia.getTempFactor());
            intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
        } else {
            atom.setAtomType(atomType);
        }
        bond(ia, atom);
        return atom;
    }

    private void bond(Atom a1, Atom a2) {
        Bond bond = new Bond(a1, a2);
        int c[] = new int[2];
        c[0] = a1.getAtomType().atomClass;
        c[1] = a2.getAtomType().atomClass;
        String key = BondType.sortKey(c);
        BondType bondType = forceField.getBondType(key);
        if (bondType == null) {
            System.out.println("No BondType for key: " + key + "\n" + a1.toString() + "\n" + a2.toString());
        } else {
            bond.setBondType(bondType);
        }
        bondList.add(bond);
    }

    /**
     * Determine the atom type based on a biotype key.
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
                System.out.println("The atom type " + bioType.atomType + " was not found.");
            }
        } else {
            //System.out.println("The biotype look-up " + lookUp + " was not found.");
        }
        return null;
    }

    public static String padRight(String s, int n) {
        return String.format("%1$-" + n + "s", s);
    }

    public static String padLeft(String s, int n) {
        return String.format("%1$#" + n + "s", s);
    }

    /**
     * Write out the Atomic information in PDB format.
     *
     * @return <code>true</code> if the read was successful.
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        if (saveFile == null) {
            return false;
        }
        // Create StringBuffers for ATOM, ANISOU and TER records.
        StringBuffer sb = new StringBuffer("ATOM  ");
        StringBuffer anisouSB = new StringBuffer("ANISOU");
        StringBuffer terSB = new StringBuffer("TER   ");
        for (int i = 6; i < 80; i++) {
            sb.append(' ');
            anisouSB.append(' ');
            terSB.append(' ');
        }
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {
            File newFile = saveFile;
            if (!append) {
                newFile = version(saveFile);
            }
            molecularAssembly.setFile(newFile);
            molecularAssembly.setName(newFile.getName());
            fw = new FileWriter(newFile, append);
            bw = new BufferedWriter(fw);
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
//         1         2         3         4         5         6         7
//123456789012345678901234567890123456789012345678901234567890123456789012345678
//ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
//ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
            int serial = 1;
            // Loop over biomolecular chains
            String chains[] = molecularAssembly.getChainNames();
            for (String chain : chains) {
                if (chain.equalsIgnoreCase("Blank")) {
                    sb.setCharAt(21, ' ');
                } else {
                    sb.setCharAt(21, chain.toUpperCase().charAt(0));
                }
                Polymer polymer = molecularAssembly.getPolymer(chain, false);
                // Loop over residues
                ArrayList<Residue> residues = polymer.getResidues();
                for (Residue residue : residues) {
                    String resName = residue.getName();
                    if (resName.length() > 3) {
                        resName = resName.substring(0, 3);
                    }
                    int resID = residue.getResidueNumber();
                    sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                    sb.replace(22, 26, String.format("%4d", resID));
                    // Loop over atoms
                    ArrayList<Atom> residueAtoms = residue.getAtomList();
                    for (Atom atom : residueAtoms) {
                        String name = atom.getID();
                        if (name.length() > 4) {
                            name = name.substring(0, 4);
                        } else if (name.length() == 1) {
                            name = name + "  ";
                        } else if (name.length() == 2) {
                            name = name + " ";
                        }
                        double xyz[] = atom.getXYZ();
                        sb.replace(6, 16, String.format("%5d " + padLeft(name.toUpperCase(), 4), serial++));
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
                }
                terSB.replace(6, 11, Integer.toString(serial++));
                terSB.replace(12, 16, "    ");
                terSB.replace(16, 26, sb.substring(16, 26));
                bw.write(terSB.toString());
                bw.newLine();
            }
            sb.replace(0, 6, "HETATM");
            // Loop over molecules, ions and then water.
            ArrayList<MSNode> molecules = new ArrayList<MSNode>();
            molecules.addAll(molecularAssembly.getMolecules());
            molecules.addAll(molecularAssembly.getIons());
            molecules.addAll(molecularAssembly.getWaters());
            for (MSNode node : molecules) {
                Molecule molecule = (Molecule) node;
                String chain = molecule.getPolymerName();
                if (chain.equalsIgnoreCase("Blank")) {
                    sb.setCharAt(21, ' ');
                } else {
                    sb.setCharAt(21, chain.toUpperCase().charAt(0));
                }
                String resName = molecule.getResidueName();
                if (resName.length() > 3) {
                    resName = resName.substring(0, 3);
                }
                int resID = molecule.getResidueNumber();
                sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
                sb.replace(22, 26, String.format("%4d", resID));
                // Loop over atoms
                ArrayList<Atom> residueAtoms = molecule.getAtomList();
                for (Atom atom : residueAtoms) {
                    String name = atom.getID();
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
                    sb.replace(6, 16, String.format("%5d " + padLeft(name.toUpperCase(), 4), serial++));
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
     * The location of a residue within a chain.
     */
    public enum ResiduePosition {

        FIRST_RESIDUE, MIDDLE_RESIDUE, LAST_RESIDUE
    };

    public enum AminoAcid1 {

        G, A, V, L, I, S, T, C, c, P, F, Y, W, H, U, Z,
        D, N, E, Q, M, K, R, O, B, J, f, a, n, m, X
    };

    public enum AminoAcid3 {

        GLY, ALA, VAL, LEU, ILE, SER, THR, CYS, CYX, PRO, PHE, TYR, TRP, HIS, HID, HIE,
        ASP, ASN, GLU, GLN, MET, LYS, ARG, ORN, AIB, PCA, FOR, ACE, NH2, NME, UNK
    };
    public final int aminoAcidHeavyAtoms[] = {4, 5, 7, 8, 8, 6, 7, 6, 6, 7, 11, 12, 14, 10, 10, 10,
        8, 8, 9, 9, 8, 9, 11, 8, 6, 8, 0, 0, 0, 0, 0};
    static AminoAcid3 knownAminoAcids[] = AminoAcid3.values();
    static int numberOfKnownAminoAcids = knownAminoAcids.length;

    public enum NucleicAcid1 {

        A, G, C, U, D, B, I, T, O, W, H, X
    };

    public enum NucleicAcid3 {

        ADE, GUA, CYT, URI, DAD, DGU, DCY, DTY, MP1, DP2, TP3, UNK
    };
    static NucleicAcid3 knownNucleicAcids[] = NucleicAcid3.values();
    static int numberOfKnownNucleicAcids = knownNucleicAcids.length;

    private enum HETATOMS {

        HOH, H2O, WAT, NA, K, MG, MG2, CA, CA2, CL
    };
    /**
     * Biotype keys for nucleic acid backbone atom types.
     */
    private final int o5Typ[] = {1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0};
    private final int c5Typ[] = {1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0};
    private final int h51Typ[] = {1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0};
    private final int h52Typ[] = {1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0};
    private final int c4Typ[] = {1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0};
    private final int h4Typ[] = {1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0};
    private final int o4Typ[] = {1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0};
    private final int c1Typ[] = {1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0};
    private final int h1Typ[] = {1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0};
    private final int c3Typ[] = {1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0};
    private final int h3Typ[] = {1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0};
    private final int c2Typ[] = {1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0};
    private final int h21Typ[] = {1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0};
    private final int h22Typ[] = {1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0};
    private final int o3Typ[] = {1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0};
    private final int o2Typ[] = {1014, 1044, 1075, 1103, 0, 0, 0, 0, 0, 0, 0, 0};
    private final int pTyp[] = {1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0};
    private final int opTyp[] = {1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0};
    private final int h5tTyp[] = {1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0};
    private final int h3tTyp[] = {1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0};
    /**
     * Biotype keys for amino acid backbone atom types.
     *
     * ntyp[0][..] are for N-terminal residues.
     * ntyp[1][..] are mid-chain residues.
     * ntyp[2][..] are for C-terminal residues.
     */
    private static final int nType[][] = {
        {350, 356, 362, 368, 374, 380, 386, 392, 398, 404, 412, 418, 424, 430, 436, 442, 448,
            454, 460, 466, 472, 478, 484, 490, 496, 325, 0, 0, 0, 0, 350},
        {1, 7, 15, 27, 41, 55, 65, 77, 87, 96, 107, 122, 138, 161, 178, 194, 210, 220,
            232, 244, 258, 271, 287, 304, 318, 325, 0, 0, 0, 0, 1},
        {501, 507, 513, 519, 525, 531, 537, 543, 549, 555, 560, 566, 572, 578, 584, 590, 596,
            602, 608, 614, 620, 626, 632, 638, 644, 0, 0, 0, 344, 346, 501}};
    private static final int caType[][] = {
        {351, 357, 363, 369, 375, 381, 387, 393, 399, 405, 413, 419, 425, 431, 437, 443, 449,
            455, 461, 467, 473, 479, 485, 491, 497, 326, 0, 340, 0, 0, 351},
        {2, 8, 16, 28, 42, 56, 66, 78, 88, 97, 108, 123, 139, 162, 179, 195, 211, 221,
            233, 245, 259, 272, 288, 305, 319, 326, 0, 0, 0, 0, 2},
        {502, 508, 514, 520, 526, 532, 538, 544, 550, 556, 561, 567, 573, 579, 585, 591, 597,
            603, 609, 615, 621, 627, 633, 639, 645, 0, 0, 0, 0, 348, 502}};
    private static final int cType[][] = {{
            352, 358, 364, 370, 376, 382, 388, 394, 400, 406, 414, 420, 426, 432, 438, 444, 450,
            456, 462, 468, 474, 480, 486, 492, 498, 327, 337, 342, 0, 0, 352},
        {3, 9, 17, 29, 43, 57, 67, 79, 89, 98, 109, 124, 140, 163, 180, 196, 212, 222,
            234, 246, 260, 273, 289, 306, 320, 327, 0, 0, 0, 0, 3},
        {503, 509, 515, 521, 527, 533, 539, 545, 551, 557, 562, 568, 574, 580, 586, 592, 598,
            604, 610, 616, 622, 628, 634, 640, 646, 0, 0, 0, 0, 0, 503}};
    private static final int hnType[][] = {{
            353, 359, 365, 371, 377, 383, 389, 395, 401, 407, 415, 421, 427, 433, 439, 445, 451,
            457, 463, 469, 475, 481, 487, 493, 499, 328, 0, 0, 0, 0, 353},
        {4, 10, 18, 30, 44, 58, 68, 80, 90, 0, 110, 125, 141, 164, 181, 197, 213, 223,
            235, 247, 261, 274, 290, 307, 321, 328, 0, 0, 0, 0, 4},
        {504, 510, 516, 522, 528, 534, 540, 546, 552, 0, 563, 569, 575, 581, 587, 593, 599,
            605, 611, 617, 623, 629, 635, 641, 647, 0, 0, 0, 345, 347, 504}};
    private static final int oType[][] = {{
            354, 360, 366, 372, 378, 384, 390, 396, 402, 408, 416, 422, 428, 434, 440, 446, 452,
            458, 464, 470, 476, 482, 488, 494, 500, 329, 339, 343, 0, 0, 354},
        {5, 11, 19, 31, 45, 59, 69, 81, 91, 99, 111, 126, 142, 165, 182, 198, 214, 224,
            236, 248, 262, 275, 291, 308, 322, 329, 0, 0, 0, 0, 5},
        {505, 511, 517, 523, 529, 535, 541, 547, 553, 558, 564, 570, 576, 582, 588, 594, 600,
            606, 612, 618, 624, 630, 636, 642, 648, 0, 0, 0, 0, 0, 505}};
    private static final int haType[][] = {{
            355, 361, 367, 373, 379, 385, 391, 397, 403, 409, 417, 423, 429, 435, 441, 447, 453,
            459, 465, 471, 477, 483, 489, 495, 0, 330, 338, 341, 0, 0, 355},
        {6, 12, 20, 32, 46, 60, 70, 82, 92, 100, 112, 127, 143, 166, 183, 199, 215, 225,
            237, 249, 263, 276, 292, 309, 0, 330, 0, 0, 0, 0, 6},
        {506, 512, 518, 524, 530, 536, 542, 548, 554, 559, 565, 571, 577, 583, 589, 595, 601,
            607, 613, 619, 625, 631, 637, 643, 0, 0, 0, 0, 0, 349, 506}};
}
