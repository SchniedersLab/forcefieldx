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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;

import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Hybrid36;
import static ffx.potential.bonded.BondedUtils.numberAtoms;
import static ffx.potential.bonded.PolymerUtils.assignAtomTypes;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_2;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_3;
import static ffx.utilities.StringUtils.padLeft;
import static ffx.utilities.StringUtils.padRight;

/**
 * The BioJavaFilter class parses data from a Biojava 5 Structure object.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class BioJavaFilter extends ConversionFilter {

    private static final Logger logger = Logger.getLogger(BioJavaFilter.class.getName());
    private final Structure structure;
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
     * systems, chain names are sometimes reused due to limitations in the PBD
     * format.
     * <p>
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
    private PDBFilter.PDBFileStandard fileStandard = VERSION3_3; // Assume current standard.
    /**
     * If true, output is directed into arrayOutput instead of the file.
     */
    private boolean listMode = false;
    private ArrayList<String> listOutput = new ArrayList<>();

    /**
     * <p>Constructor for BiojavaFilter.</p>
     *
     * @param structure         a {@link org.biojava.nbio.structure.Structure} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forcefield        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public BioJavaFilter(Structure structure, MolecularAssembly molecularAssembly, ForceField forcefield, CompositeConfiguration properties) {
        super(structure, molecularAssembly, forcefield, properties);
        this.structure = structure;
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
    }

    /**
     * <p>Constructor for BiojavaFilter.</p>
     *
     * @param structures        a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param forcefield        a {@link ffx.potential.parameters.ForceField} object.
     * @param properties        a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public BioJavaFilter(List<Structure> structures, MolecularAssembly molecularAssembly, ForceField forcefield, CompositeConfiguration properties) {
        super(structures, molecularAssembly, forcefield, properties);
        if (structures != null && !structures.isEmpty()) {
            this.structure = structures.get(0);
        } else {
            structure = null;
        }
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
    }

    /**
     * <p>Constructor for BiojavaFilter.</p>
     *
     * @param structure           a {@link org.biojava.nbio.structure.Structure} object.
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param forcefield          a {@link ffx.potential.parameters.ForceField} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public BioJavaFilter(Structure structure, List<MolecularAssembly> molecularAssemblies, ForceField forcefield, CompositeConfiguration properties) {
        super(structure, molecularAssemblies, forcefield, properties);
        this.structure = structure;
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
    }

    /**
     * Mutate a residue at the PDB file is being parsed.
     *
     * @param chainID the Chain ID of the residue to mutate.
     * @param resID   the Residue ID of the residue to mutate.
     * @param name    the 3-letter code of the amino acid to mutate to.
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
     * Get the list of alternate locations encountered.
     *
     * @return the alternate location list.
     */
    public List<Character> getAltLocs() {
        return altLocs;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Parse the Biojava Structure
     */
    @Override
    public boolean convert() {
        // First atom is #1, to match xyz file format
        int xyzIndex = 1;
        setConverted(false);
        systems.add(activeMolecularAssembly);

        if (mutate) {
            List<Character> chainIDs = new ArrayList<>();
            for (Chain chain : structure.getChains()) {
                chainIDs.add(chain.getChainID().charAt(0));
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
        }

        // Echo the alternate location being parsed.
        if (currentAltLoc == 'A') {
            logger.info(String.format(" Reading %s", structure.getName()));
        } else {
            logger.info(String.format(" Reading %s alternate location %s",
                    structure.getName(), currentAltLoc));
        }

        org.biojava.nbio.structure.Atom[] bjAtoms = StructureTools.getAllAtomArray(structure);
        // Reset the current chain and segID.
        currentChainID = null;
        currentSegID = null;
        PDBCrystallographicInfo cInfo = structure.getCrystallographicInfo();

        if (cInfo != null) {
            // I do not think we need to check if it already has these properties, but it can be done.
            properties.addProperty("a-axis", cInfo.getA());
            properties.addProperty("b-axis", cInfo.getB());
            properties.addProperty("c-axis", cInfo.getC());
            properties.addProperty("alpha", cInfo.getAlpha());
            properties.addProperty("beta", cInfo.getBeta());
            properties.addProperty("gamma", cInfo.getGamma());
            properties.addProperty("spacegroup", cInfo.getSpaceGroup().getShortSymbol());
        }

        for (org.biojava.nbio.structure.Atom atom : bjAtoms) {
            String name = atom.getName().toUpperCase().trim();
            double[] xyz = new double[3];
            xyz[0] = atom.getX();
            xyz[1] = atom.getY();
            xyz[2] = atom.getZ();
            char altLoc = atom.getAltLoc();
            if (!altLocs.contains(altLoc)) {
                altLocs.add(altLoc);
            }
            if (altLoc != ' ' && altLoc != 'A' && altLoc != currentAltLoc) {
                break;
            }

            if (name.contains("1H") || name.toUpperCase().contains("2H")
                    || name.toUpperCase().contains("3H")) {
                // VERSION3_2 is presently just a placeholder for "anything non-standard".
                fileStandard = VERSION3_2;
            }

            Group group = atom.getGroup();
            ResidueNumber resnum = group.getResidueNumber();
            int resSeq = resnum.getSeqNum();
            String resName = group.getPDBName().trim().toUpperCase();

            Chain chain = group.getChain();
            char chainID = chain.getChainID().charAt(0);
            String segID = getSegID(chainID);

            boolean printAtom = false;
            if (mutate && chainID == mutateChainID && mutateResID == resSeq) {
                if (name.equals("N") || name.equals("C") || name.equals("O")
                        || name.equals("CA")) {
                    printAtom = true;
                    name = mutateToResname;
                } else {
                    logger.info(String.format(" Deleting atom %s of %s %d", name, resName, resSeq));
                    break;
                }
            }

            // Biojava does not maintain ANISOU records.
            Atom newAtom = new Atom(0,
                    name, altLoc, xyz, resName, resSeq, chainID, atom.getOccupancy(),
                    atom.getTempFactor(), segID);

            /* Biojava sets at least some capping groups, and possibly nonstandard amino acids to be heteroatoms. */
            boolean hetatm = true;
            for (AminoAcid3 aa3Name : AminoAcid3.values()) {
                if (aa3Name.name().equals(resName)) {
                    hetatm = false;
                    break;
                }
            }
            newAtom.setHetero(hetatm);

            Atom returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
            if (returnedAtom != newAtom) {
                atoms.put(atom.getPDBserial(), returnedAtom);
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(String.format("%s has been retained over\n%s",
                            returnedAtom.toString(), newAtom.toString()));
                }
            } else {
                atoms.put(atom.getPDBserial(), newAtom);
                if (newAtom.getIndex() == 0) {
                    newAtom.setXyzIndex(xyzIndex++);
                }
                if (printAtom) {
                    logger.info(newAtom.toString());
                }
            }
        }

        // Locate disulfide bonds; bond parameters are assigned below.
        // List<Bond> ssBondList = locateDisulfideBonds(ssbonds, activeMolecularAssembly, pdbToNewResMap);

        int pdbAtoms = activeMolecularAssembly.getAtomArray().length;

        // Build missing backbone atoms in loops.
        // buildMissingResidues(xyzIndex, activeMolecularAssembly, seqRes, dbRef);

        // Assign atom types. Missing side-chains atoms and missing hydrogens will be built in.
        bondList = assignAtomTypes(activeMolecularAssembly, fileStandard);

        // Assign disulfide bonds parameters and log their creation.
        // buildDisulfideBonds(ssBondList, activeMolecularAssembly, bondList);

        // Finally, re-number the atoms if missing atoms were created.
        if (pdbAtoms != activeMolecularAssembly.getAtomArray().length) {
            numberAtoms(activeMolecularAssembly);
        }

        return true;
    }

    /**
     * <p>
     * writeFileWithHeader</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param header   a {@link java.lang.StringBuilder} object.
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
     * <p>
     * writeFile</p>
     *
     * @param saveFile    a {@link java.io.File} object.
     * @param append      a {@link java.lang.StringBuilder} object.
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
        return writeFile(saveFile, append, false);
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
    public void writeAtom(Atom atom, int serial, StringBuilder sb, StringBuilder anisouSB, BufferedWriter bw)
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
}
