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

import org.apache.commons.configuration.CompositeConfiguration;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SSBond;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.numerics.VectorMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
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
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter.HetAtoms;
import ffx.utilities.Hybrid36;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
import ffx.potential.bonded.AminoAcidUtils;
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
import ffx.potential.bonded.NucleicAcidUtils;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidList;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcid;
import static ffx.potential.bonded.ResidueEnumerations.nucleicAcidList;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_2;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_3;
import static ffx.utilities.StringUtils.padLeft;
import static ffx.utilities.StringUtils.padRight;

/**
 * The BiojavaFilter class parses data from a Biojava Structure object.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 *
 */
public class BiojavaFilter extends ConversionFilter {

    private static final Logger logger = Logger.getLogger(BiojavaFilter.class.getName());
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
    private PDBFilter.PDBFileStandard fileStandard = VERSION3_3; // Assume current standard.
    /**
     * If true, output is directed into arrayOutput instead of the file.
     */
    private boolean listMode = false;
    private ArrayList<String> listOutput = new ArrayList<>();
    /**
     * Don't output atoms which fail Atom.getUse().
     */
    private boolean ignoreInactiveAtoms = false;

    public BiojavaFilter(Structure structure, MolecularAssembly molecularAssembly, ForceField forcefield, CompositeConfiguration properties) {
        super(structure, molecularAssembly, forcefield, properties);
        this.structure = structure;
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
    }

    public BiojavaFilter(List<Structure> structures, MolecularAssembly molecularAssembly, ForceField forcefield, CompositeConfiguration properties) {
        super(structures, molecularAssembly, forcefield, properties);
        if (structures != null && !structures.isEmpty()) {
            this.structure = structures.get(0);
        } else {
            structure = null;
        }
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
    }

    public BiojavaFilter(Structure structure, List<MolecularAssembly> molecularAssemblies, ForceField forcefield, CompositeConfiguration properties) {
        super(structure, molecularAssemblies, forcefield, properties);
        this.structure = structure;
        this.dataType = Utilities.DataType.BIOJAVA;
        this.fileType = Utilities.FileType.PDB;
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
     * Parse the Biojava Structure
     *
     * @return true if the structure is successfully converted
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
        /**
         * Echo the alternate location being parsed.
         */
        if (currentAltLoc == 'A') {
            logger.info(String.format(" Reading %s", structure.getName()));
        } else {
            logger.info(String.format(" Reading %s alternate location %s",
                    structure.getName(), currentAltLoc));
        }

        org.biojava.nbio.structure.Atom[] bjAtoms = StructureTools.getAllAtomArray(structure);
        int nAtoms = bjAtoms.length;
        Atom[] ffxAtoms = new Atom[nAtoms];
        /**
         * Reset the current chain and segID.
         */
        currentChainID = null;
        currentSegID = null;
        

        if (structure.isCrystallographic()) {
            // I do not think we need to check if it already has these properties,
            // but it can be done.
            PDBCrystallographicInfo cInfo = structure.getCrystallographicInfo();
            properties.addProperty("a-axis", cInfo.getA());
            properties.addProperty("b-axis", cInfo.getB());
            properties.addProperty("c-axis", cInfo.getC());
            properties.addProperty("alpha", cInfo.getAlpha());
            properties.addProperty("beta", cInfo.getBeta());
            properties.addProperty("gamma", cInfo.getGamma());
            properties.addProperty("spacegroup", SpaceGroup.pdb2ShortName(cInfo.getSpaceGroup().getShortSymbol()));
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
            Character insCode = atom.getGroup().getResidueNumber().getInsCode();
            int pdbSer = atom.getPDBserial();

            boolean printAtom = false;
            if (mutate && chainID == mutateChainID && mutateResID == resSeq) {
                if (name.equals("N") || name.equals("C") || name.equals("O")
                        || name.equals("CA")) {
                    printAtom = true;
                    name = mutateToResname;
                } else {
                    logger.info(String.format(" Deleting atom %s of %s %d",
                            name, resName, resSeq));
                    break;
                }
            }

            Atom newAtom = new Atom(0,
                    name, altLoc, xyz, resName, resSeq, chainID, atom.getOccupancy(),
                    atom.getTempFactor(), segID, true);
            if (insCode != null) {
                newAtom.setInsCode(insCode);
            }
            newAtom.setPDBserial(pdbSer);

            /* Biojava sets at least some capping groups, and possibly nonstandard
             amino acids to be heteroatoms. May wish to revisit at some point. */
            boolean hetatm = true;
            for (AminoAcid3 aa3Name : AminoAcid3.values()) {
                if (aa3Name.name().equals(resName)) {
                    hetatm = false;
                    break;
                }
            }
            if (hetatm) {
                for (NucleicAcid3 na3Name : NucleicAcid3.values()) {
                    if (na3Name.name().equals(resName)) {
                        hetatm = false;
                        break;
                    }
                }
            }
            newAtom.setHetero(hetatm);

            // Biojava essentially ignores ANISOU records.
            Atom returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
            if (returnedAtom != newAtom) {
                atoms.put(atom.getPDBserial(), returnedAtom);
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(String.format("%s has been retained over\n%s",
                            returnedAtom.toString(), newAtom.toString()));
                }
            } else {
                atoms.put(atom.getPDBserial(), newAtom);
                if (newAtom.xyzIndex == 0) {
                    newAtom.setXYZIndex(xyzIndex++);
                }
                if (printAtom) {
                    logger.info(newAtom.toString());
                }
            }
        }

        List<Bond> ssBondList = new ArrayList<>();
        for (SSBond ssbond : structure.getSSBonds()) {
            Polymer c1 = activeMolecularAssembly.getChain(ssbond.getChainID1());
            Polymer c2 = activeMolecularAssembly.getChain(ssbond.getChainID2());
            int rn1;
            int rn2;
            try {
                rn1 = Integer.parseInt(ssbond.getResnum1());
                rn2 = Integer.parseInt(ssbond.getResnum1());
            } catch (NumberFormatException ex) {
                logger.warning(String.format(" Could not parse SSbond %d", ssbond.getSerNum()));
                continue;
            }
            Residue r1 = c1.getResidue(rn1);
            Residue r2 = c2.getResidue(rn2);
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
                String message = String.format("Ignoring [%s]\n due to distance %8.3f A.", ssbond.toString(), d);
                logger.log(Level.WARNING, message);
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
            String key = ffx.potential.parameters.BondType.sortKey(c);
            ffx.potential.parameters.BondType bondType = forceField.getBondType(key);
            if (bondType == null) {
                logger.severe(String.format("No BondType for key: %s\n %s\n %s", key, a1, a2));
            } else {
                bond.setBondType(bondType);
            }
            double d = VectorMath.dist(a1.getXYZ(null), a2.getXYZ(null));
            Polymer c1 = activeMolecularAssembly.getChain(a1.getSegID());
            Polymer c2 = activeMolecularAssembly.getChain(a2.getSegID());
            Residue r1 = c1.getResidue(a1.getResidueNumber());
            Residue r2 = c2.getResidue(a2.getResidueNumber());
            sb.append(String.format("\n S-S distance of %6.2f for %s and %s.", d, r1.toString(), r2.toString()));
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
    
    public static Atom readAtom(org.biojava.nbio.structure.Atom atom, String segID) {
        String name = atom.getName().toUpperCase().trim();
        double[] xyz = new double[3];
        xyz[0] = atom.getX();
        xyz[1] = atom.getY();
        xyz[2] = atom.getZ();
        char altLoc = atom.getAltLoc();

        Group group = atom.getGroup();
        ResidueNumber resnum = group.getResidueNumber();
        int resSeq = resnum.getSeqNum();
        String resName = group.getPDBName().trim().toUpperCase();

        Chain chain = group.getChain();
        char chainID = chain.getChainID().charAt(0);

        boolean printAtom = false;

        Atom newAtom = new Atom(0,
                name, altLoc, xyz, resName, resSeq, chainID, atom.getOccupancy(),
                atom.getTempFactor(), segID, true);

        /* Biojava sets at least some capping groups, and possibly nonstandard
         amino acids to be heteroatoms. */
        boolean hetatm = true;
        for (AminoAcid3 aa3Name : AminoAcid3.values()) {
            if (aa3Name.name().equals(resName)) {
                hetatm = false;
                break;
            }
        }
        newAtom.setHetero(hetatm);
        return newAtom;
    }
    
    public static void addAtomToAssembly(org.biojava.nbio.structure.Atom atom, String segID, MolecularAssembly assembly) {
        Atom newAtom = readAtom(atom, segID);
        Atom returnedAtom = (Atom) assembly.addMSNode(newAtom);
        if (returnedAtom == newAtom) {
            returnedAtom.setXYZIndex(assembly.getMaxXYZIndex() + 1);
        }
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

        Polymer[] polymers = activeMolecularAssembly.getPolymers();
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
        Polymer[] polymers = activeMolecularAssembly.getPolymers();

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
                            AminoAcidUtils.renameToProtonationState(residue);
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
                        NucleicAcidUtils.assignNucleicAcidAtomTypes(residues, forceField, bondList);
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
                            Atom atom2 = molecule.getMoleculeAtom(name);
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
                            ia.getOccupancy(), ia.getTempFactor(), ia.getSegID(), true);
                    logger.fine(" Created hydrogen " + atomName + ".");
                    hydrogen.setAtomType(type);
                    hydrogen.setHetero(true);
                    molecule.addMSNode(hydrogen);
                    int valence = ia.getAtomType().valence;
                    List<Bond> aBonds = ia.getFFXBonds();
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
            Polymer polymers[] = activeMolecularAssembly.getPolymers();
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
                            List<Bond> bonds = SG1.getFFXBonds();
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
                    currentChainID = polymer.getChainIDChar();
                    sb.setCharAt(21, currentChainID);
                    // Loop over residues
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        String resName = residue.getName();
                        if (resName.length() > 3) {
                            resName = resName.substring(0, 3);
                        }
                        int resID = residue.getResidueIndex();
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
                    int resID2 = residue.getResidueIndex();
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

    /*public Structure writeToStructure(String header) {
     return writeToStructure(activeMolecularAssembly, header);
     }

     public static Structure writeToStructure(MolecularAssembly assembly, String header) {
     Structure structure = new StructureImpl();
     for (Polymer polymer : assembly.getChains()) {
     Chain chain = new ChainImpl();
     for (Residue residue : polymer.getResidues()) {
     Group group;
     switch (residue.getResidueType()) {
     case AA:
     group = new AminoAcidImpl();
     break;
     case NA:
     group = new NucleotideImpl();
     break;
     default:
     group = new HetatomImpl();
     break;
     }
     for (Atom atom : residue.getAtomList()) {
     org.biojava.nbio.structure.Atom bjAtom = new AtomImpl();
     group.addAtom(bjAtom);
     }
     chain.addGroup(group);
     }
     structure.addChain(chain);
     }

     for (Molecule molecule : assembly.getMolecules()) {
     for (Atom atom : molecule.getAtomList()) {

     }
     }
     }*/
    public boolean writeSIFTFile(File saveFile, boolean append, String[] resAndScore) {
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
            Polymer polymers[] = activeMolecularAssembly.getPolymers();
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
                            List<Bond> bonds = SG1.getFFXBonds();
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
                    currentChainID = polymer.getChainIDChar();
                    sb.setCharAt(21, currentChainID);
                    // Loop over residues
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        String resName = residue.getName();
                        if (resName.length() > 3) {
                            resName = resName.substring(0, 3);
                        }
                        int resID = residue.getResidueIndex();
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
                    int resID2 = residue.getResidueIndex();
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
    
    /*public static PDBCrystallographicInfo generateCrystalInfo(MolecularAssembly assembly) {
        Crystal crystal = assembly.getCrystal();
        if (crystal.aperiodic()) {
            return null;
        } else {
            PDBCrystallographicInfo cInfo = new PDBCrystallographicInfo();
            CrystalCell cell = new CrystalCell(crystal.a, crystal.b, crystal.c, 
                    crystal.alpha, crystal.beta, crystal.gamma);
            cInfo.setCrystalCell(cell);
            SpaceGroup sgroup = crystal.spaceGroup;
            int id = sgroup.number;
            int multiplicity = sgroup.getNumberOfSymOps();
        }
        if (structure.isCrystallographic()) {
            // I do not think we need to check if it already has these properties,
            // but it can be done.
            PDBCrystallographicInfo cInfo = structure.getCrystallographicInfo();
            properties.addProperty("a-axis", cInfo.getA());
            properties.addProperty("b-axis", cInfo.getB());
            properties.addProperty("c-axis", cInfo.getC());
            properties.addProperty("alpha", cInfo.getAlpha());
            properties.addProperty("beta", cInfo.getBeta());
            properties.addProperty("gamma", cInfo.getGamma());
            properties.addProperty("spacegroup", SpaceGroup.pdb2ShortName(cInfo.getSpaceGroup().getShortSymbol()));
        }
    }
    
    public static void setCrystalProperties(MolecularAssembly assembly, PDBCrystallographicInfo cInfo) {
        
    }*/

    public void setIgnoreInactiveAtoms(boolean ignoreInactiveAtoms) {
        this.ignoreInactiveAtoms = ignoreInactiveAtoms;
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
        if (ignoreInactiveAtoms && !atom.isActive()) {
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
        //sb.replace(78, 80, String.format("%2d", 0));
        sb.replace(78, 80, String.format("%2d", atom.getFormalCharge()));
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
        double xyz[] = vdwH ? atom.getRedXYZ() : atom.getXYZ(null);
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
            String minusNumbers = atomName.replaceAll("[1-9]", "");
            if (minusNumbers.startsWith("D") || minusNumbers.startsWith("H")) {
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
}
