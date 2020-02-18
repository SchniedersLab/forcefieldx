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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter.PDBFileStandard;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dist;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;
import static ffx.potential.bonded.AminoAcidUtils.assignAminoAcidAtomTypes;
import static ffx.potential.bonded.Bond.logNoBondType;
import static ffx.potential.bonded.BondedUtils.buildBond;
import static ffx.potential.bonded.BondedUtils.buildHeavy;
import static ffx.potential.bonded.BondedUtils.buildHydrogen;
import static ffx.potential.bonded.BondedUtils.findAtomType;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static ffx.potential.bonded.NamingUtils.checkHydrogenAtomNames;
import static ffx.potential.bonded.NucleicAcidUtils.assignNucleicAcidAtomTypes;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidList;
import static ffx.potential.bonded.ResidueEnumerations.nucleicAcidList;

/**
 * Utilities for creating polymers.
 *
 * @author Michael Schnieders
 * @since 1.0
 */
public class PolymerUtils {

    private static final Logger logger = Logger.getLogger(PolymerUtils.class.getName());

    private static final double DEFAULT_AA_CHAINBREAK = 3.0;
    private static final double DEFAULT_NA_CHAINBREAK_MULT = 2.0;

    /**
     * Assign force field atoms types to common chemistries using "biotype" records.
     *
     * @param molecularAssembly MolecularAssembly to operate on.
     * @param fileStandard      PDB file standard to follow.
     * @return An ArrayList of created bonds.
     */
    public static ArrayList<Bond> assignAtomTypes(MolecularAssembly molecularAssembly,
                                                  PDBFileStandard fileStandard) {
        // Create a list to store bonds defined by PDB atom names.
        ArrayList<Bond> bondList = new ArrayList<>();

        ForceField forceField = molecularAssembly.getForceField();
        CompositeConfiguration properties = molecularAssembly.getProperties();

        // To Do: Look for cyclic peptides and disulfides.
        Polymer[] polymers = molecularAssembly.getChains();

        // Loop over chains.
        if (polymers != null) {
            logger.info(format("\n Assigning atom types for %d chains.", polymers.length));
            for (Polymer polymer : polymers) {
                List<Residue> residues = polymer.getResidues();
                boolean isProtein = true;
                // Check if all residues are known amino acids.
                for (Residue residue : residues) {
                    String name = residue.getName().toUpperCase();
                    boolean aa = false;
                    for (ResidueEnumerations.AminoAcid3 amino : aminoAcidList) {
                        if (amino.toString().equalsIgnoreCase(name)) {
                            aa = true;
                            checkHydrogenAtomNames(residue, fileStandard);
                            break;
                        }
                    }
                    // Check for a patch.
                    if (!aa) {
                        if (logger.isLoggable(Level.FINE)) {
                            logger.fine(" Checking for non-standard amino acid patch " + name);
                        }
                        HashMap<String, AtomType> types = forceField.getAtomTypes(name);
                        if (types.isEmpty()) {
                            isProtein = false;
                            break;
                        } else {
                            if (logger.isLoggable(Level.FINE)) {
                                logger.fine(" Patch found for non-standard amino acid " + name);
                            }
                        }
                    }
                }

                // If all the residues in this chain have known amino acids names, then attempt to assign atom types.
                if (isProtein) {
                    try {
                        logger.info(format(" Amino acid chain %s", polymer.getName()));
                        double dist = properties.getDouble("chainbreak", DEFAULT_AA_CHAINBREAK);
                        // Detect main chain breaks!
                        List<List<Residue>> subChains = findChainBreaks(residues, dist);
                        for (List<Residue> subChain : subChains) {
                            assignAminoAcidAtomTypes(subChain, forceField, bondList);
                        }
                    } catch (BondedUtils.MissingHeavyAtomException missingHeavyAtomException) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(missingHeavyAtomException));
                        logger.severe(missingHeavyAtomException.toString());
                    } catch (BondedUtils.MissingAtomTypeException missingAtomTypeException) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(missingAtomTypeException));
                        logger.severe(missingAtomTypeException.toString());
                    }
                    continue;
                }

                // Check if all residues have known nucleic acids names.
                boolean isNucleicAcid = true;
                for (Residue residue : residues) {
                    String name = residue.getName().toUpperCase();
                    // Convert 1 and 2-character nucleic acid names to 3-character names.
                    switch (name) {
                        case "A":
                            name = ResidueEnumerations.NucleicAcid3.ADE.toString();
                            break;
                        case "C":
                            name = ResidueEnumerations.NucleicAcid3.CYT.toString();
                            break;
                        case "G":
                            name = ResidueEnumerations.NucleicAcid3.GUA.toString();
                            break;
                        case "T":
                            name = ResidueEnumerations.NucleicAcid3.THY.toString();
                            break;
                        case "U":
                            name = ResidueEnumerations.NucleicAcid3.URI.toString();
                            break;
                        case "YG":
                            name = ResidueEnumerations.NucleicAcid3.YYG.toString();
                            break;
                        case "DA":
                            name = ResidueEnumerations.NucleicAcid3.DAD.toString();
                            break;
                        case "DC":
                            name = ResidueEnumerations.NucleicAcid3.DCY.toString();
                            break;
                        case "DG":
                            name = ResidueEnumerations.NucleicAcid3.DGU.toString();
                            break;
                        case "DT":
                            name = ResidueEnumerations.NucleicAcid3.DTY.toString();
                            break;
                    }
                    residue.setName(name);
                    ResidueEnumerations.NucleicAcid3 nucleicAcid = null;
                    for (ResidueEnumerations.NucleicAcid3 nucleic : nucleicAcidList) {
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

                /*
                  If all the residues in this chain have known nucleic acids
                  names, then attempt to assign atom types.
                 */
                if (isNucleicAcid) {
                    try {
                        logger.info(format(" Nucleic acid chain %s", polymer.getName()));

                        if (logger.isLoggable(Level.FINE)) {
                            logger.fine(format(" EXPERIMENTAL: Finding chain breaks for possible nucleic acid chain %s", polymer.getName()));
                        }
                        double dist = properties.getDouble("chainbreak", DEFAULT_AA_CHAINBREAK);
                        dist *= DEFAULT_NA_CHAINBREAK_MULT;
                        // Detect main chain breaks!
                        List<List<Residue>> subChains = findChainBreaks(residues, dist);

                        for (List<Residue> subChain : subChains) {
                            assignNucleicAcidAtomTypes(subChain, forceField, bondList);
                        }
                    } catch (BondedUtils.MissingHeavyAtomException | BondedUtils.MissingAtomTypeException e) {
                        logger.log(Level.INFO, Utilities.stackTraceToString(e));
                        logger.severe(e.toString());
                    }
                }
            }
        }

        // Assign ion atom types.
        ArrayList<MSNode> ions = molecularAssembly.getIons();
        if (ions != null && ions.size() > 0) {
            logger.info(format(" Assigning atom types for %d ions.", ions.size()));
            for (MSNode m : ions) {
                Molecule ion = (Molecule) m;
                String name = ion.getResidueName().toUpperCase();
                NamingUtils.HetAtoms hetatm = NamingUtils.HetAtoms.parse(name);
                Atom atom = ion.getAtomList().get(0);
                if (ion.getAtomList().size() != 1) {
                    logger.severe(format(" Check residue %s of chain %s.", ion.toString(), ion.getChainID()));
                }
                try {
                    switch (hetatm) {
                        case NA:
                            atom.setAtomType(findAtomType(2004, forceField));
                            break;
                        case K:
                            atom.setAtomType(findAtomType(2005, forceField));
                            break;
                        case MG:
                        case MG2:
                            atom.setAtomType(findAtomType(2008, forceField));
                            break;
                        case CA:
                        case CA2:
                            atom.setAtomType(findAtomType(2009, forceField));
                            break;
                        case ZN:
                        case ZN2:
                            atom.setAtomType(findAtomType(2010, forceField));
                            break;
                        case CL:
                            atom.setAtomType(findAtomType(2013, forceField));
                            break;
                        case BR:
                            atom.setAtomType(findAtomType(2012, forceField));
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
        ArrayList<MSNode> water = molecularAssembly.getWaters();
        if (water != null && water.size() > 0) {
            logger.info(format(" Assigning atom types for %d waters.", water.size()));
            for (MSNode m : water) {
                Molecule wat = (Molecule) m;
                try {
                    Atom O = buildHeavy(wat, "O", null, 2001, forceField, bondList);
                    Atom H1 = buildHydrogen(wat, "H1", O, 0.96e0, null, 109.5e0, null, 120.0e0, 0, 2002, forceField, bondList);
                    Atom H2 = buildHydrogen(wat, "H2", O, 0.96e0, H1, 109.5e0, null, 120.0e0, 0, 2002, forceField, bondList);
                    O.setHetero(true);
                    H1.setHetero(true);
                    H2.setHetero(true);
                } catch (Exception e) {
                    logger.log(Level.INFO, Utilities.stackTraceToString(e));
                    String message = " Error assigning atom types to a water.";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }

        // Assign small molecule atom types.
        ArrayList<MSNode> molecules = molecularAssembly.getMolecules();
        for (MSNode m : molecules) {
            Molecule molecule = (Molecule) m;
            String moleculeName = molecule.getResidueName();
            logger.info(" Attempting to patch " + moleculeName);
            ArrayList<Atom> moleculeAtoms = molecule.getAtomList();
            boolean patched = true;
            HashMap<String, AtomType> types = forceField.getAtomTypes(moleculeName);
            // Assign atom types for all known atoms.
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
            // Create missing hydrogen atoms. Check for missing heavy atoms.
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
                    String[] bonds = forceField.getBonds(moleculeName, atomName);
                    if (bonds != null) {
                        for (String name : bonds) {
                            Atom atom2 = molecule.getAtom(name);
                            if (atom2 != null && !atom.isBonded(atom2)) {
                                buildBond(atom, atom2, forceField, bondList);
                            }
                        }
                    }
                }
            }
            // Create missing hydrogen atoms.
            if (patched && !types.isEmpty()) {
                // Create a map of the molecule's atoms
                Map<String, Atom> atomMap = new HashMap<>();
                for (Atom atom : moleculeAtoms) {
                    atomMap.put(atom.getName().toUpperCase(), atom);
                }
                for (String atomName : types.keySet()) {
                    AtomType type = types.get(atomName);
                    String[] bonds = forceField.getBonds(moleculeName, atomName.toUpperCase());
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
                    // Try to find the following configuration: ib-ia-ic
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

                    // Building the hydrogens depends on hybridization and the locations of other bonded atoms.
                    logger.fine(" Bonding " + atomName + " to " + ia.getName()
                            + " (" + numBonds + " of " + valence + ").");
                    switch (valence) {
                        case 4:
                            switch (numBonds) {
                                case 3:
                                    // Find the average coordinates of atoms ib, ic and id.
                                    double[] b = ib.getXYZ(null);
                                    double[] c = ic.getXYZ(null);
                                    double[] d = id.getXYZ(null);
                                    double[] a = new double[3];
                                    a[0] = (b[0] + c[0] + d[0]) / 3.0;
                                    a[1] = (b[1] + c[1] + d[1]) / 3.0;
                                    a[2] = (b[2] + c[2] + d[2]) / 3.0;
                                    // Place the hydrogen at chiral position #1.
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 1);
                                    double[] e1 = new double[3];
                                    hydrogen.getXYZ(e1);
                                    double[] ret = new double[3];
                                    diff(a, e1, ret);
                                    double l1 = r(ret);
                                    // Place the hydrogen at chiral position #2.
                                    intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, -1);
                                    double[] e2 = new double[3];
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
                            if (numBonds == 0) {
                                intxyz(hydrogen, ia, 1.0, null, 0.0, null, 0.0, 0);
                            } else {
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
                        buildBond(ia, hydrogen, forceField, bondList);
                    }
                }
            }
            if (!patched) {
                logger.log(Level.WARNING, format(" Deleting unrecognized molecule %s.", m.toString()));
                molecularAssembly.deleteMolecule((Molecule) m);
            } else {
                logger.info(" Patch for " + moleculeName + " succeeded.");
            }
        }
        resolvePolymerLinks(molecules, molecularAssembly, bondList);
        return bondList;
    }

    /**
     * Locate disulfide bonds based on SSBOND records.
     *
     * @param ssbonds           List of SSBOND records.
     * @param molecularAssembly The MolecularAssembly to operate on.
     * @param pdbToNewResMap    Maps chainIDResNumInsCode to renumbered chainIDResNum. For example,
     *                          residue 52A in chain C might be renumbered to residue 53, and mapped as
     *                          "C52A" to "C53".
     * @return An ArrayList of Bond instances for SS-Bonds.
     */
    public static List<Bond> locateDisulfideBonds(List<String> ssbonds, MolecularAssembly molecularAssembly, Map<String, String> pdbToNewResMap) {
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
                Polymer c1 = molecularAssembly.getChain(format("%c", c1ch));
                Polymer c2 = molecularAssembly.getChain(format("%c", c2ch));

                Polymer[] chains = molecularAssembly.getChains();
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

                String pdbResNum1 = format("%c%s%c", c1ch, origResNum1, insChar1);
                String pdbResNum2 = format("%c%s%c", c2ch, origResNum2, insChar2);
                String resnum1 = pdbToNewResMap.get(pdbResNum1);
                String resnum2 = pdbToNewResMap.get(pdbResNum2);
                if (resnum1 == null) {
                    logger.warning(format(" Could not find residue %s for SS-bond %s", pdbResNum1, ssbond));
                    continue;
                }
                if (resnum2 == null) {
                    logger.warning(format(" Could not find residue %s for SS-bond %s", pdbResNum2, ssbond));
                    continue;
                }

                Residue r1 = c1.getResidue(parseInt(resnum1.substring(1)));
                Residue r2 = c2.getResidue(parseInt(resnum2.substring(1)));
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
                    logger.warning(format(" SG atom 1 of SS-bond %s is null", ssbond));
                }
                if (SG2 == null) {
                    logger.warning(format(" SG atom 2 of SS-bond %s is null", ssbond));
                }
                if (SG1 == null || SG2 == null) {
                    continue;
                }
                double d = dist(SG1.getXYZ(null), SG2.getXYZ(null));
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
     * @param ssBondList        List of SSBOND records.
     * @param molecularAssembly MolecularAssembly to operate on.
     * @param bondList          Add new SS-Bonds to this list.
     */
    public static void buildDisulfideBonds(List<Bond> ssBondList, MolecularAssembly molecularAssembly, ArrayList<Bond> bondList) {
        StringBuilder sb = new StringBuilder(" Disulfide Bonds:");
        ForceField forceField = molecularAssembly.getForceField();
        for (Bond bond : ssBondList) {
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            int[] c = new int[2];
            c[0] = a1.getAtomType().atomClass;
            c[1] = a2.getAtomType().atomClass;
            String key = BondType.sortKey(c);
            BondType bondType = forceField.getBondType(key);
            if (bondType == null) {
                logNoBondType(a1, a2, key);
            } else {
                bond.setBondType(bondType);
            }
            double d = dist(a1.getXYZ(null), a2.getXYZ(null));
            Polymer c1 = molecularAssembly.getChain(a1.getSegID());
            Polymer c2 = molecularAssembly.getChain(a2.getSegID());
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
     * Known limitations include:
     * 1) No building n- and c-terminal loops.
     * 2) No support for DBREF1 or DBREF2 records.
     * 3) Incomplete optimization scheme to position the loops.
     *
     * @param xyzIndex          XYZ index to begin from.
     * @param molecularAssembly MolecularAssembly to operate on.
     * @param seqres            Map of SEQRES entries.
     * @param dbref             Map of DBREF entries.
     * @return xyzIndex updated based on built atoms.
     */
    public static int buildMissingResidues(int xyzIndex, MolecularAssembly molecularAssembly,
                                           Map<Character, String[]> seqres,
                                           Map<Character, int[]> dbref) {

        // Only build loops if the buildLoops flag is true.
        CompositeConfiguration properties = molecularAssembly.getProperties();
        if (!properties.getBoolean("buildLoops", false)) {
            return xyzIndex;
        }

        Polymer[] polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            Character chainID = polymer.getChainID();
            String[] resNames = seqres.get(chainID);
            int[] seqRange = dbref.get(chainID);
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
                // Residues at the end of the chain are not currently built.
                if (nextResidue == null) {
                    logger.info(format(" Residue %d is missing, but is at the end of the chain.",
                            currentID));
                    break;
                }
                // Find the previous carbonyl carbon and next nitrogen.
                Atom C = (Atom) previousResidue.getAtomNode("C");
                Atom N = (Atom) nextResidue.getAtomNode("N");
                if (C == null || N == null) {
                    logger.info(format(" Residue %d is missing, but bonding atoms are missing (C or N).",
                            currentID));
                    continue;
                }

                // Build the missing residue.
                currentResidue = polymer.getResidue(resNames[i], currentID, true);

                double[] vector = new double[3];
                int count = 3 * (nextResidue.getResidueNumber() - previousResidue.getResidueNumber());
                diff(N.getXYZ(null), C.getXYZ(null), vector);
                scalar(vector, 1.0 / count, vector);

                double[] nXYZ = new double[3];
                sum(C.getXYZ(null), vector, nXYZ);
                nXYZ[0] += random() - 0.5;
                nXYZ[1] += random() - 0.5;
                nXYZ[2] += random() - 0.5;
                Atom newN = new Atom(xyzIndex++, "N", C.getAltLoc(), nXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newN);

                double[] caXYZ = new double[3];
                scalar(vector, 2.0, vector);
                sum(C.getXYZ(null), vector, caXYZ);
                caXYZ[0] += Math.random() - 0.5;
                caXYZ[1] += Math.random() - 0.5;
                caXYZ[2] += Math.random() - 0.5;
                Atom newCA = new Atom(xyzIndex++, "CA", C.getAltLoc(), caXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newCA);

                double[] cXYZ = new double[3];
                scalar(vector, 1.5, vector);
                sum(C.getXYZ(null), vector, cXYZ);
                cXYZ[0] += Math.random() - 0.5;
                cXYZ[1] += Math.random() - 0.5;
                cXYZ[2] += Math.random() - 0.5;
                Atom newC = new Atom(xyzIndex++, "C", C.getAltLoc(), cXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newC);

                double[] oXYZ = new double[3];
                vector[0] = Math.random() - 0.5;
                vector[1] = Math.random() - 0.5;
                vector[2] = Math.random() - 0.5;
                sum(cXYZ, vector, oXYZ);
                Atom newO = new Atom(xyzIndex++, "O", C.getAltLoc(), oXYZ, resNames[i], currentID,
                        chainID, 1.0, C.getTempFactor(), C.getSegID(), true);
                currentResidue.addMSNode(newO);
                logger.info(format(" Building residue %8s.", currentResidue.toString()));
            }
        }
        return xyzIndex;
    }

    /**
     * Resolves links between polymeric hetero groups; presently only functional
     * for cyclic molecules.
     *
     * @param molecules         List of Molecules in the molecular assembly.
     * @param molecularAssembly MolecularAssembly to operate on.
     * @param bondList          Add created bonds to this list.
     */
    public static void resolvePolymerLinks(List<MSNode> molecules, MolecularAssembly molecularAssembly,
                                           ArrayList<Bond> bondList) {

        ForceField forceField = molecularAssembly.getForceField();
        CompositeConfiguration properties = molecularAssembly.getProperties();

        for (String polyLink : properties.getStringArray("polymerlink")) {
            logger.info(" Experimental: linking a cyclic hetero group: " + polyLink);

            // Format: polymerlink resname atom1 atom2 [cyclize]
            String[] toks = polyLink.split("\\s+");
            String resName = toks[0];
            String name1 = toks[1];
            String name2 = toks[2];
            int cyclicLen = 0;
            if (toks.length > 3) {
                cyclicLen = parseInt(toks[3]);
            }

            ArrayList<Molecule> matches = new ArrayList<>();
            for (MSNode node : molecules) {
                Molecule m = (Molecule) node;
                if (m.getResidueName().equalsIgnoreCase(resName)) {
                    matches.add(m);
                }
            }

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
                        logger.info(format(" Cyclizing molecule %s to %s", mi, next));
                    } else {
                        next = matches.get(ii);
                        logger.info(format(" Extending chain from %s to %s.", mi, next));
                    }
                    Atom from = mi.getAtomByName(name1, true);
                    Atom to = next.getAtomByName(name2, true);
                    buildBond(from, to, forceField, bondList);
                }
            }
        }
    }

    public static List<List<Residue>> findChainBreaks(List<Residue> residues, double cutoff) {
        List<List<Residue>> subChains = new ArrayList<>();

        // Chain-start atom: N (amino)/O5* (nucleic)
        // Chain-end atom:   C (amino)/O3* (nucleic)
        Residue.ResidueType rType = residues.get(0).getResidueType();
        String startAtName;
        String endAtName;
        switch (rType) {
            case AA:
                startAtName = "N";
                endAtName = "C";
                break;
            case NA:
                boolean namedStar = residues.stream().
                        flatMap((Residue r) -> r.getAtomList().stream()).
                        anyMatch((Atom a) -> a.getName().equals("O5*"));
                if (namedStar) {
                    startAtName = "O5*";
                    endAtName = "O3*";
                } else {
                    startAtName = "O5\'";
                    endAtName = "O3\'";
                }
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
                // Initialization.
                subChain = new ArrayList<>();
                subChain.add(residue);
                subChains.add(subChain);
            } else {
                // Find the start atom of the current residue.
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
                // Compute the distance between the previous carbonyl carbon and the current nitrogen.
                double r = dist(priorEndAtom.getXYZ(null), startAtom.getXYZ(null));
                if (r > cutoff) {
                    // Start a new chain.
                    subChain = new ArrayList<>();
                    subChain.add(residue);
                    subChains.add(subChain);
                    char ch1 = previousResidue.getChainID();
                    char ch2 = residue.getChainID();
                    sb.append(format("\n C-N distance of %6.2f A for %c-%s and %c-%s.",
                            r, ch1, previousResidue.toString(), ch2, residue.toString()));
                } else {
                    // Continue the current chain.
                    subChain.add(residue);
                }
            }

            // Save the carbonyl carbon.
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

}
