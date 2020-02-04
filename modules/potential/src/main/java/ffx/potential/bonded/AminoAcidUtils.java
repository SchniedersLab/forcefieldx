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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.bonded.Residue.AA3;
import ffx.potential.bonded.Residue.ResiduePosition;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dihedralAngle;
import static ffx.numerics.math.VectorMath.r;
import static ffx.potential.bonded.BondedUtils.buildBond;
import static ffx.potential.bonded.BondedUtils.buildHeavy;
import static ffx.potential.bonded.BondedUtils.buildHydrogen;
import static ffx.potential.bonded.BondedUtils.buildHydrogenAtom;
import static ffx.potential.bonded.BondedUtils.findAtomType;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static ffx.potential.bonded.Residue.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.bonded.ResidueEnumerations.aminoAcidHeavyAtoms;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcid;
import static ffx.potential.bonded.ResidueEnumerations.getAminoAcidNumber;

/**
 * Utilities for creating Amino Acid residues.
 *
 * @author Michael Schnieders
 * @since 1.0
 */
public class AminoAcidUtils {

    /**
     * Private constructor.
     */
    private AminoAcidUtils() {
    }

    private static final Logger logger = Logger.getLogger(AminoAcidUtils.class.getName());

    /**
     * Assign atom types to an amino acid polymer.
     *
     * @param residues The residues to assign atom types to.
     * @throws MissingHeavyAtomException A needed heavy atom was not found.
     * @throws MissingAtomTypeException  An atom type could not be found.
     * @since 1.0
     */
    public static void assignAminoAcidAtomTypes(List<Residue> residues, ForceField forceField, ArrayList<Bond> bondList)
            throws MissingHeavyAtomException, MissingAtomTypeException {
        // Loop over amino acid residues.
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
     * <p>assignAminoAcidAtomTypes.</p>
     *
     * @param residue         a {@link ffx.potential.bonded.Residue} object.
     * @param previousResidue a {@link ffx.potential.bonded.Residue} object.
     * @param nextResidue     a {@link ffx.potential.bonded.Residue} object.
     * @param forceField      a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList        a {@link java.util.ArrayList} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     * @throws ffx.potential.bonded.BondedUtils.MissingAtomTypeException  if any.
     */
    public static void assignAminoAcidAtomTypes(Residue residue, Residue previousResidue,
                                                Residue nextResidue, ForceField forceField, ArrayList<Bond> bondList)
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
            // If the last residue only contains a nitrogen turn it into an NH2 group.
            Atom N = (Atom) residue.getAtomNode("N");
            if (residue.getAtomNodeList().size() == 1 && N != null) {
                residueName = "NH2".intern();
                residue.setName(residueName);
            }
        }

        AminoAcid3 aminoAcid = getAminoAcid(residueName);
        int aminoAcidNumber = getAminoAcidNumber(residueName);
        // Non-standard Amino Acid; use ALA backbone types.
        boolean nonStandard = false;
        if (aminoAcid == AminoAcid3.UNK) {
            aminoAcidNumber = getAminoAcidNumber("ALA");
            nonStandard = true;
        }

        // Only the last residue in a chain should have an OXT/OT2 atom.
        if (nextResidue != null) {
            removeOXT_OT2(residue);
        }

        // Only the first nitrogen should have H1, H2 and H3 atoms, unless it's an NME cap.
        if (previousResidue != null) {
            removeH1_H2_H3(aminoAcid, residue);
        }

        // Check for missing heavy atoms. This check ignores special terminating groups like FOR, NH2, etc.
        if (!nonStandard) {
            try {
                checkForMissingHeavyAtoms(aminoAcidNumber, aminoAcid, position, residue);
            } catch (BondedUtils.MissingHeavyAtomException e) {
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

        // Backbone heavy atoms.
        Atom N = (Atom) residue.getAtomNode("N");
        if (N != null) {
            N.setAtomType(findAtomType(AA_N[j][aminoAcidNumber], forceField));
            if (position != FIRST_RESIDUE) {
                buildBond(pC, N, forceField, bondList);
            }
        }

        Atom CA = null;
        Atom C = null;
        Atom O = null;
        if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NH2)) {
            if (aminoAcid == AminoAcid3.ACE || aminoAcid == AminoAcid3.NME) {
                CA = buildHeavy(residue, "CH3", N, AA_CA[j][aminoAcidNumber], forceField, bondList);
            } else {
                CA = buildHeavy(residue, "CA", N, AA_CA[j][aminoAcidNumber], forceField, bondList);
            }
            if (!(position == LAST_RESIDUE && aminoAcid == AminoAcid3.NME)) {
                C = buildHeavy(residue, "C", CA, AA_C[j][aminoAcidNumber], forceField, bondList);
                O = (Atom) residue.getAtomNode("O");
                if (O == null) {
                    O = (Atom) residue.getAtomNode("OT1");
                }
                if (O == null) {
                    // Build the carbonyl oxygen; assume the final residue has beta-sheet psi
                    double psi = 135.0;
                    if (nextResidue != null && N != null) {
                        Atom nextN = (Atom) nextResidue.getAtomNode("N");
                        psi = toDegrees(dihedralAngle(N.getXYZ(null), CA.getXYZ(null), C.getXYZ(null), nextN.getXYZ(null)));
                    }
                    O = buildHeavy(residue, "O", C, 1.25, CA, 117.0, N, psi - 180.0, 0,
                            AA_O[j][aminoAcidNumber], forceField, bondList);
                } else {
                    AtomType atomType = findAtomType(AA_O[j][aminoAcidNumber], forceField);
                    O.setAtomType(atomType);
                    buildBond(C, O, forceField, bondList);
                }
            }
        }

        // Nitrogen hydrogen atoms.
        AtomType atomType = findAtomType(AA_HN[j][aminoAcidNumber], forceField);
        switch (position) {
            case FIRST_RESIDUE:
                switch (aminoAcid) {
                    case PRO:
                        // Check if, for some reason, the structure has an H1 but not an H2 or H3.
                        List<Atom> resAtoms = residue.getAtomList();
                        Optional<Atom> H1 = resAtoms.stream().filter((Atom a) -> a.getName().equals("H1")).findAny();
                        Optional<Atom> H2 = resAtoms.stream().filter((Atom a) -> a.getName().equals("H2")).findAny();
                        Optional<Atom> H3 = resAtoms.stream().filter((Atom a) -> a.getName().equals("H3")).findAny();
                        if (H1.isPresent() && H2.isPresent() && H3.isPresent()) {
                            logger.severe(String.format(" Proline residue %s somehow has three N-terminal hydrogens!", residue));
                        }
                        if (H1.isPresent() && H2.isPresent()) {
                            H1.get().setName("H3");
                        } else if (H1.isPresent() && H3.isPresent()) {
                            H1.get().setName("H2");
                        } // Else, all is hunky-dory!
                        buildHydrogenAtom(residue, "H2", N, 1.02, CA, 109.5, C, 0.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H3", N, 1.02, CA, 109.5, C, -120.0, 0,
                                atomType, forceField, bondList);
                        break;
                    case PCA:
                        buildHydrogenAtom(residue, "H", N, 1.02, CA, 109.5, C, -60.0, 0,
                                atomType, forceField, bondList);
                        break;
                    case ACE:
                        break;
                    default:
                        buildHydrogenAtom(residue, "H1", N, 1.02, CA, 109.5, C, 180.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H2", N, 1.02, CA, 109.5, C, 60.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H3", N, 1.02, CA, 109.5, C, -60.0, 0,
                                atomType, forceField, bondList);
                }
                break;
            case LAST_RESIDUE:
                switch (aminoAcid) {
                    case NH2:
                        buildHydrogenAtom(residue, "H1", N, 1.02, pC, 119.0, pCA, 0.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H2", N, 1.02, pC, 119.0, pCA, 180.0, 0,
                                atomType, forceField, bondList);
                        break;
                    case NME:
                        buildHydrogenAtom(residue, "H", N, 1.02, pC, 118.0, CA, 121.0, 1,
                                atomType, forceField, bondList);
                        break;
                    default:
                        buildHydrogenAtom(residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1,
                                atomType, forceField, bondList);
                }
                break;
            default:
                // Mid-chain nitrogen hydrogen.
                buildHydrogenAtom(residue, "H", N, 1.02, pC, 119.0, CA, 119.0, 1,
                        atomType, forceField, bondList);
        }

        // C-alpha hydrogen atoms.
        String haName = "HA";
        if (aminoAcid == AminoAcid3.GLY) {
            haName = "HA2";
        }
        atomType = findAtomType(AA_HA[j][aminoAcidNumber], forceField);
        switch (position) {
            case FIRST_RESIDUE:
                switch (aminoAcid) {
                    case FOR:
                        buildHydrogenAtom(residue, "H", C, 1.12, O, 0.0, null, 0.0, 0,
                                atomType, forceField, bondList);
                        break;
                    case ACE:
                        buildHydrogenAtom(residue, "H1", CA, 1.10, C, 109.5, O, 180.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H2", CA, 1.10, C, 109.5, O, 60.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H3", CA, 1.10, C, 109.5, O, -60.0, 0,
                                atomType, forceField, bondList);
                        break;
                    default:
                        buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1,
                                atomType, forceField, bondList);
                        break;
                }
                break;
            case LAST_RESIDUE:
                switch (aminoAcid) {
                    case NME:
                        buildHydrogenAtom(residue, "H1", CA, 1.10, N, 109.5, pC, 180.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H2", CA, 1.10, N, 109.5, pC, 60.0, 0,
                                atomType, forceField, bondList);
                        buildHydrogenAtom(residue, "H3", CA, 1.10, N, 109.5, pC, -60.0, 0,
                                atomType, forceField, bondList);
                        break;
                    default:
                        buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.5, -1,
                                atomType, forceField, bondList);
                }
                break;
            default:
                buildHydrogenAtom(residue, haName, CA, 1.10, N, 109.5, C, 109.0, -1,
                        atomType, forceField, bondList);
        }

        // Build the amino acid side chain.
        assignAminoAcidSideChain(position, aminoAcid, residue, CA, N, C, forceField, bondList);

        // Build the terminal oxygen if the residue is not NH2 or NME.
        if (position == LAST_RESIDUE && !(aminoAcid == AminoAcid3.NH2 || aminoAcid == AminoAcid3.NME)) {
            atomType = findAtomType(AA_O[2][aminoAcidNumber], forceField);
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
            buildBond(C, OXT, forceField, bondList);
        }

        // Do some checks on the current residue to make sure all atoms have been assigned an atom type.
        List<Atom> resAtoms = residue.getAtomList();
        for (Atom atom : resAtoms) {
            atomType = atom.getAtomType();
            if (atomType == null) {
                /**
                 * Sometimes, with deuterons, a proton has been constructed in
                 * its place, so we have a "dummy" deuteron still hanging
                 * around.
                 */
                String protonEq = atom.getName().replaceFirst("D", "H");
                Atom correspH = (Atom) residue.getAtomNode(protonEq);
                if (correspH == null || correspH.getAtomType() == null) {
                    MissingAtomTypeException missingAtomTypeException
                            = new MissingAtomTypeException(residue, atom);
                    logger.warning("MissingAtomTypeException incoming from 393.");
                    throw missingAtomTypeException;
                } else {
                    correspH.setName(atom.getName());
                    atom.removeFromParent();
                    atom = correspH;
                    atomType = atom.getAtomType();
                }

            }
            int numberOfBonds = atom.getNumBonds();
            if (numberOfBonds != atomType.valence) {
                if (atom == C && numberOfBonds == atomType.valence - 1 && position != LAST_RESIDUE) {
                    continue;
                }
                logger.warning(format(" An atom for residue %s has the wrong number of bonds:\n %s",
                        residueName, atom.toString()));
                logger.info(format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
                for (Bond bond : atom.getBonds()) {
                    logger.info(" " + bond.toString());
                }
            }
        }
    }

    /**
     * Assign atom types to a single amino acid side chain.
     *
     * @param position   The position of this amino acid in the chain.
     * @param aminoAcid  The amino acid to use.
     * @param residue    The residue node.
     * @param CA         The C-alpha carbon of this residue.
     * @param N          The peptide nitrogen of this residue.
     * @param C          The peptide carbonyl carbon.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     */
    private static void assignAminoAcidSideChain(ResiduePosition position, AminoAcid3 aminoAcid, Residue residue,
                                                 Atom CA, Atom N, Atom C, ForceField forceField, ArrayList<Bond> bondList) {
        int k = AA_CB[aminoAcid.ordinal()];
        switch (aminoAcid) {
            case GLY:
                switch (position) {
                    case FIRST_RESIDUE:
                        k = AA_HA[0][k];
                        break;
                    case LAST_RESIDUE:
                        k = AA_HA[2][k];
                        break;
                    default:
                        k = AA_HA[1][k];

                }
                buildHydrogen(residue, "HA3", CA, 1.10, N, 109.5, C, 109.5, 1, k, forceField, bondList);
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
                if (!types.isEmpty()) {
                    boolean patched = true;
                    ArrayList<Atom> residueAtoms = residue.getAtomList();
                    // Assign atom types for side-chain atoms.
                    for (Atom atom : residueAtoms) {
                        String atomName = atom.getName().toUpperCase();
                        AtomType type = atom.getAtomType();
                        if (type == null) {
                            type = types.get(atomName);
                            atom.setAtomType(type);
                            types.remove(atomName);
                        }
                    }
                    // Create bonds between known atoms.
                    for (Atom atom : residueAtoms) {
                        String atomName = atom.getName();
                        String[] bonds = forceField.getBonds(residueName, atomName);
                        if (bonds != null) {
                            for (String name : bonds) {
                                Atom atom2 = (Atom) residue.getAtomNode(name);
                                if (atom2 != null && !atom.isBonded(atom2)) {
                                    buildBond(atom, atom2, forceField, bondList);
                                }
                            }
                        }
                    }

                    // Create missing hydrogen atoms.
                    if (!types.isEmpty()) {
                        // Create a hashmap of the molecule's atoms
                        HashMap<String, Atom> atomMap = new HashMap<>();
                        for (Atom atom : residueAtoms) {
                            atomMap.put(atom.getName().toUpperCase(), atom);
                        }
                        for (String atomName : types.keySet()) {
                            AtomType type = types.get(atomName);
                            String[] bonds = forceField.getBonds(residueName, atomName.toUpperCase());
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
                            logger.log(Level.FINE, " Bonding {0} to {1} ({2} of {3}).",
                                    new Object[]{atomName, ia.getName(), numBonds, valence});
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
                                            intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 0);
                                            double e1[] = hydrogen.getXYZ(null);
                                            double ret[] = new double[3];
                                            diff(a, e1, ret);
                                            double l1 = r(ret);

                                            // Place the hydrogen at chiral position #2.
                                            intxyz(hydrogen, ia, 1.0, ib, 109.5, ic, 109.5, 1);
                                            double e2[] = hydrogen.getXYZ(null);
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
                                buildBond(ia, hydrogen, forceField, bondList);
                            }
                        }
                    }
                    if (!patched) {
                        logger.log(Level.SEVERE, format(" Could not patch %s.", residueName));
                    } else {
                        logger.log(Level.INFO, " Patch for {0} succeeded.", residueName);
                        residueAtoms = residue.getAtomList();
                        // Assign atom types for side-chain atoms.
                        double charge = 0.0;
                        for (Atom atom : residueAtoms) {
                            logger.info(atom.toString() + " -> " + atom.getAtomType().toString());
                        }
                    }
                } else {
                    switch (position) {
                        case FIRST_RESIDUE:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 355,
                                    forceField, bondList);
                            break;
                        case LAST_RESIDUE:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 506,
                                    forceField, bondList);
                            break;
                        default:
                            buildHydrogen(residue, "HA2", CA, 1.10e0, N, 109.5e0, C, 109.5e0, 1, 6,
                                    forceField, bondList);
                    }
                }
                break;
        }
    }

    /**
     * <p>copyResidue.</p>
     *
     * @param fromResidue a {@link ffx.potential.bonded.Residue} object.
     * @param toResidue   a {@link ffx.potential.bonded.Residue} object.
     */
    public static void copyResidue(Residue fromResidue, Residue toResidue) {
        String resName = fromResidue.getName();
        AA3 res = AA3.valueOf(resName);
        ArrayList<String> atomNames = new ArrayList<>();
        switch (res) {
            case ALA:
                atomNames.addAll(Arrays.asList(alanineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case GLY:
                atomNames.addAll(Arrays.asList(glycineAtoms));
                atomNames.addAll(Arrays.asList(glycineBackboneAtoms));
                break;
            case PRO:
                atomNames.addAll(Arrays.asList(prolineAtoms));
                atomNames.addAll(Arrays.asList(prolineBackboneAtoms));
                break;
            case VAL:
                atomNames.addAll(Arrays.asList(valineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case LEU:
                atomNames.addAll(Arrays.asList(leucineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case ILE:
                atomNames.addAll(Arrays.asList(isoleucineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case SER:
                atomNames.addAll(Arrays.asList(serineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case THR:
                atomNames.addAll(Arrays.asList(threonineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case CYS:
                atomNames.addAll(Arrays.asList(cysteineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case PHE:
                atomNames.addAll(Arrays.asList(phenylalanineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case TYR:
                atomNames.addAll(Arrays.asList(tyrosineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case TRP:
                atomNames.addAll(Arrays.asList(tryptophanAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case HIS:
                atomNames.addAll(Arrays.asList(histidineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case ASP:
                atomNames.addAll(Arrays.asList(aspartateAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case ASN:
                atomNames.addAll(Arrays.asList(asparagineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case GLU:
                atomNames.addAll(Arrays.asList(glutamateAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case GLN:
                atomNames.addAll(Arrays.asList(glutamineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case MET:
                atomNames.addAll(Arrays.asList(methionineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case LYS:
                atomNames.addAll(Arrays.asList(lysineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            case ARG:
                atomNames.addAll(Arrays.asList(arginineAtoms));
                atomNames.addAll(Arrays.asList(backboneAtoms));
                break;
            default:
                atomNames = null;
        }
        for (String atomName : atomNames) {
            copyCoordinates(fromResidue, toResidue, atomName);
        }
    }

    /**
     * Only the first nitrogen should have H1, H2 and H3 atoms, unless it's an
     * NME cap.
     *
     * @param aminoAcid 3-letter amino acid name.
     * @param residue   the amino acid Residue.
     */
    public static void removeH1_H2_H3(AminoAcid3 aminoAcid, Residue residue) {
        if (aminoAcid != AminoAcid3.NME) {
            if (aminoAcid != AminoAcid3.NH2) {
                Atom H1 = (Atom) residue.getAtomNode("H1");
                if (H1 != null) {
                    residue.deleteAtom(H1);
                }
                Atom H2 = (Atom) residue.getAtomNode("H2");
                if (H2 != null) {
                    residue.deleteAtom(H2);
                }
            }
            Atom H3 = (Atom) residue.getAtomNode("H3");
            if (H3 != null) {
                residue.deleteAtom(H3);
            }
        }
    }

    /**
     * Only the last residue in a chain should have an OXT/OT2 atom.
     *
     * @param residue the amino acid residue.
     */
    public static void removeOXT_OT2(Residue residue) {
        Atom OXT = (Atom) residue.getAtomNode("OXT");
        if (OXT != null) {
            residue.deleteAtom(OXT);
        }
        Atom OT2 = (Atom) residue.getAtomNode("OT2");
        if (OT2 != null) {
            residue.deleteAtom(OT2);
        }
    }

    /**
     * <p>buildAlanine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildAlanine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                       ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom HB1 = buildHydrogen(residue, "HB1", CB, 1.11, CA, 109.4, N, 180.0, 0, k + 1, forceField, bondList);
        buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, HB1, 109.4, 1, k + 1, forceField, bondList);
        buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, HB1, 109.4, -1, k + 1, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildValine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildValine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                      ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG1 = buildHeavy(residue, "CG1", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, CG1, 109.5, -1, k + 4, forceField, bondList);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, CG1, 109.4, 1, k + 1, forceField, bondList);
        Atom HG11 = buildHydrogen(residue, "HG11", CG1, 1.11, CB, 109.4, CA, 180.0, 0, k + 3, forceField, bondList);
        Atom HG12 = buildHydrogen(residue, "HG12", CG1, 1.11, CB, 109.4, HG11, 109.4, 1, k + 3, forceField, bondList);
        Atom HG13 = buildHydrogen(residue, "HG13", CG1, 1.11, CB, 109.4, HG11, 109.4, -1, k + 3, forceField, bondList);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 109.4, CA, 180.0, 0, k + 5, forceField, bondList);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 109.4, HG21, 109.4, 1, k + 5, forceField, bondList);
        Atom HG23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 109.4, HG21, 109.4, -1, k + 5, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildLeucine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildLeucine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                       ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.54, CB, 109.5, CD1, 109.5, -1, k + 6, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG = buildHydrogen(residue, "HG", CG, 1.11, CB, 109.4, CD1, 109.4, 1, k + 3, forceField, bondList);
        Atom HD11 = buildHydrogen(residue, "HD11", CD1, 1.11, CG, 109.4, CB, 180.0, 0, k + 5, forceField, bondList);
        Atom HD12 = buildHydrogen(residue, "HD12", CD1, 1.11, CG, 109.4, HD11, 109.4, 1, k + 5, forceField, bondList);
        Atom HD13 = buildHydrogen(residue, "HD13", CD1, 1.11, CG, 109.4, HD11, 109.4, -1, k + 5, forceField, bondList);
        Atom HD21 = buildHydrogen(residue, "HD21", CD2, 1.11, CG, 109.4, CB, 180.0, 0, k + 7, forceField, bondList);
        Atom HD22 = buildHydrogen(residue, "HD22", CD2, 1.11, CG, 109.4, HD21, 109.4, 1, k + 7, forceField, bondList);
        Atom HD23 = buildHydrogen(residue, "HD23", CD2, 1.11, CG, 109.4, HD21, 109.4, -1, k + 7, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildIsoleucine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildIsoleucine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                          ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom CG1 = buildHeavy(residue, "CG1", CB, 1.54, CA, 109.5, N, 0, 0, k + 2, forceField, bondList);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, CG1, 109.5, 1, k + 4, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG1, 1.54, CB, 109.5, CA, 180, 0, k + 6, forceField, bondList);
        //  CD1 = setHeavy(residue, "CD", CG1, 1.54, CB, 109.5, CA, 180, 0, k + 6, forceField, bondList);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, CG2, 109.4, 1, k + 1, forceField, bondList);
        Atom HG12 = buildHydrogen(residue, "HG12", CG1, 1.11, CB, 109.4, CD1, 109.4, 1, k + 3, forceField, bondList);
        Atom HG13 = buildHydrogen(residue, "HG13", CG1, 1.11, CB, 109.4, CD1, 109.4, -1, k + 3, forceField, bondList);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 110.0, CG1, 180.0, 0, k + 5, forceField, bondList);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 110.0, HG21, 109.0, 1, k + 5, forceField, bondList);
        Atom HG23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 110.0, HG21, 109.0, -1, k + 5, forceField, bondList);
        Atom HD11 = buildHydrogen(residue, "HD11", CD1, 1.11, CG1, 110.0, CB, 180.0, 0, k + 7, forceField, bondList);
        Atom HD12 = buildHydrogen(residue, "HD12", CD1, 1.11, CG1, 110.0, HD11, 109.0, 1, k + 7, forceField, bondList);
        Atom HD13 = buildHydrogen(residue, "HD13", CD1, 1.11, CG1, 110.0, HD11, 109.0, -1, k + 7, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildSerine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildSerine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                      ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom OG = buildHeavy(residue, "OG", CB, 1.41, CA, 107.5, N, 180, 0, k + 2, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, OG, 106.7, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, OG, 106.7, -1, k + 1, forceField, bondList);
        Atom HG = buildHydrogen(residue, "HG", OG, 0.94, CB, 106.9, CA, 180.0, 0, k + 3, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildThreonine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildThreonine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom OG1 = buildHeavy(residue, "OG1", CB, 1.41, CA, 107.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CG2 = buildHeavy(residue, "CG2", CB, 1.54, CA, 109.5, OG1, 107.7, 1, k + 4, forceField, bondList);
        Atom HB = buildHydrogen(residue, "HB", CB, 1.11, CA, 109.4, OG1, 106.7, -1, k + 1, forceField, bondList);
        Atom HG1 = buildHydrogen(residue, "HG1", OG1, 0.94, CB, 106.9, CA, 180.0, 0, k + 3, forceField, bondList);
        Atom HG21 = buildHydrogen(residue, "HG21", CG2, 1.11, CB, 110.0, CA, 180.0, 0, k + 5, forceField, bondList);
        Atom HG22 = buildHydrogen(residue, "HG22", CG2, 1.11, CB, 110.0, HG21, 109.0, 1, k + 5, forceField, bondList);
        Atom HD23 = buildHydrogen(residue, "HG23", CG2, 1.11, CB, 110.0, HG21, 109.0, -1, k + 5, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildCysteine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildCysteine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                        ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1, forceField, bondList);
        Atom HG = buildHydrogen(residue, "HG", SG, 1.34, CB, 96.0, CA, 180.0, 0, k + 3, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDeprotonatedCysteine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildDeprotonatedCysteine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildCystine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildCystine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                       ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom SG = buildHeavy(residue, "SG", CB, 1.82, CA, 109.0, N, 180, 0, k + 2, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, SG, 112.0, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, SG, 112.0, -1, k + 1, forceField, bondList);
        List<Atom> resAtoms = residue.getAtomList();
        for (Atom atom : resAtoms) {
            atom.setResName("CYS");
        }
        residue.setName("CYS");
        return residue;
    }

    /**
     * <p>buildProline.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param position   a {@link ffx.potential.bonded.Residue.ResiduePosition} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildProline(Residue residue, Atom CA, Atom N, Atom C, int k,
                                       ResiduePosition position, ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.5247, N, 104.0, C, 109.5, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.5247, CA, 104.0, N, 30.0, 0, k + 2, forceField, bondList);
        int cdKey = position == FIRST_RESIDUE ? 469 : k + 4;
        //Atom CD = buildHeavy(residue, "CD", CG, 1.5247, CB, 104.0, CA, 30.0, 0, cdKey, forceField, bondList);
        // Initial fix attempt
        Atom CD = buildHeavy(residue, "CD", N, 1.5247, CA, 104.0, CB, 0, 0, cdKey, forceField, bondList);
        /* Old code
         Atom CG = buildHeavy(residue, "CG", CB, 1.54, N, 107.0, CA, 180, 0, k + 2, forceField, bondList);
         Atom CD;
         if (position == FIRST_RESIDUE, ForceField forceField, ArrayList<Bond> bondList) {
         CD = buildHeavy(residue, "CD", CG, 1.54, CA, 107.0, CB, 180, 0, 469, forceField, bondList);
         } else {
         CD = buildHeavy(residue, "CD", CG, 1.54, CA, 107.0, CB, 180, 0, k + 4, forceField, bondList);
         }
         buildBond(CD, N, forceField, bondList);
         */
        buildBond(CD, CG, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        if (position == FIRST_RESIDUE) {
            buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, N, 109.4, 1, 470, forceField, bondList);
            buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, N, 109.4, -1, 470, forceField, bondList);
        } else {
            buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, N, 109.4, 1, k + 5, forceField, bondList);
            buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, N, 109.4, -1, k + 5, forceField, bondList);
        }
        return residue;
    }

    /**
     * <p>buildPhenylalanine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildPhenylalanine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                             ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 180, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7, forceField, bondList);
        buildBond(CE2, CZ, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.11, CG, 120.0, CE1, 120.0, 1, k + 4, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.11, CG, 120.0, CE2, 120.0, 1, k + 4, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.11, CD1, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.11, CD2, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        Atom HZ = buildHydrogen(residue, "HZ", CZ, 1.11, CE1, 120.0, CE2, 120.0, 1, k + 8, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildTyrosine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildTyrosine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                        ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 62, 0, k + 2, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 90, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7, forceField, bondList);
        buildBond(CE2, CZ, forceField, bondList);
        Atom OH = buildHeavy(residue, "OH", CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, k + 8, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 120.0, CE1, 120.0, 1, k + 4, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 120.0, CE2, 120.0, 1, k + 4, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        Atom HH = buildHydrogen(residue, "HH", OH, 0.97, CZ, 108.0, CE2, 0.0, 0, k + 9, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDeprotonatedTyrosine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildDeprotonatedTyrosine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 62, 0, k + 2, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.39, CB, 120.0, CA, 90, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.39, CB, 120.0, CD1, 120.0, 1, k + 3, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", CD1, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CE2 = buildHeavy(residue, "CE2", CD2, 1.39, CG, 120.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom CZ = buildHeavy(residue, "CZ", CE1, 1.39, CD1, 120.0, CG, 0.0, 0, k + 7, forceField, bondList);
        buildBond(CE2, CZ, forceField, bondList);
        Atom OH = buildHeavy(residue, "OH", CZ, 1.36, CE2, 120.0, CE1, 120.0, 1, k + 8, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 120.0, CE1, 120.0, 1, k + 4, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 120.0, CE2, 120.0, 1, k + 4, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, CD1, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE2, 1.10, CD2, 120.0, CZ, 120.0, 1, k + 6, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildTryptophan.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildTryptophan(Residue residue, Atom CA, Atom N, Atom C, int k,
                                          ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 62, 0, k + 2, forceField, bondList);
        Atom CD1 = buildHeavy(residue, "CD1", CG, 1.35, CB, 126.0, CA, -90, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, CD1, 108.0, 1, k + 5, forceField, bondList);
        Atom NE1 = buildHeavy(residue, "NE1", CD1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 6, forceField, bondList);
        Atom CE2 = buildHeavy(residue, "CE2", NE1, 1.35, CD1, 108.0, CG, 0.0, 0, k + 8, forceField, bondList);
        buildBond(CE2, CD2, forceField, bondList);
        Atom CE3 = buildHeavy(residue, "CE3", CD2, 1.35, CE2, 120.0, NE1, 180.0, 0, k + 9, forceField, bondList);
        Atom CZ2 = buildHeavy(residue, "CZ2", CE2, 1.35, CD2, 120.0, CE3, 0.0, 0, k + 11, forceField, bondList);
        Atom CZ3 = buildHeavy(residue, "CZ3", CE3, 1.35, CD2, 120.0, CE2, 0.0, 0, k + 13, forceField, bondList);
        Atom CH2 = buildHeavy(residue, "CH2", CZ2, 1.35, CE2, 120.0, CD2, 0.0, 0, k + 15, forceField, bondList);
        buildBond(CH2, CZ3, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", CD1, 1.10, CG, 126.0, NE1, 126.0, 1, k + 4, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", NE1, 1.05, CD1, 126.0, CE2, 126.0, 1, k + 7, forceField, bondList);
        Atom HE3 = buildHydrogen(residue, "HE3", CE3, 1.10, CD1, 120.0, CZ3, 120.0, 1, k + 10, forceField, bondList);
        Atom HZ2 = buildHydrogen(residue, "HZ2", CZ2, 1.10, CE2, 120.0, CH2, 120.0, 1, k + 12, forceField, bondList);
        Atom HZ3 = buildHydrogen(residue, "HZ3", CZ3, 1.10, CE3, 120.0, CH2, 120.0, 1, k + 14, forceField, bondList);
        Atom HH2 = buildHydrogen(residue, "HH2", CH2, 1.10, CZ2, 120.0, CZ3, 120.0, 1, k + 16, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildHistidine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildHistidine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 5, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 7, forceField, bondList);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 9, forceField, bondList);
        buildBond(NE2, CE1, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 6, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 8, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 10, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildNeutralHistidineD.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildNeutralHistidineD(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                 ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 5, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 7, forceField, bondList);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 9, forceField, bondList);
        buildBond(NE2, CE1, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD1 = buildHydrogen(residue, "HD1", ND1, 1.02, CG, 126.0, CB, 0.0, 0, k + 4, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 6, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 8, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildNeutralHistidineE.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildNeutralHistidineE(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                 ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 109.5, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.50, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom ND1 = buildHeavy(residue, "ND1", CG, 1.35, CB, 126.0, CA, 180, 0, k + 3, forceField, bondList);
        Atom CD2 = buildHeavy(residue, "CD2", CG, 1.35, CB, 126.0, ND1, 108.0, 1, k + 4, forceField, bondList);
        Atom CE1 = buildHeavy(residue, "CE1", ND1, 1.35, CG, 108.0, CD2, 0.0, 0, k + 6, forceField, bondList);
        Atom NE2 = buildHeavy(residue, "NE2", CD2, 1.35, CG, 108.0, ND1, 0.0, 0, k + 8, forceField, bondList);
        buildBond(NE2, CE1, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD2, 1.10, CG, 126.0, NE2, 126.0, 1, k + 5, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE1, 1.10, ND1, 126.0, NE2, 126.0, 1, k + 7, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", NE2, 1.02, CD2, 126.0, CE1, 126.0, 1, k + 9, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildAspartate.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildAspartate(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2, forceField, bondList);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.25, CB, 117.0, CA, 0.0, 0, k + 3, forceField, bondList);
        Atom OD2 = buildHeavy(residue, "OD2", CG, 1.25, CB, 117.0, OD1, 126.0, 1, k + 3, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildNeutralAsparticAcid.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildNeutralAsparticAcid(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                   ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2, forceField, bondList);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.25, CB, 117.0, CA, 0.0, 0, k + 3, forceField, bondList);
        Atom OD2 = buildHeavy(residue, "OD2", CG, 1.25, CB, 117.0, OD1, 126.0, 1, k + 4, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", OD2, 0.98, CG, 108.7, OD1, 0.0, 0, k + 5, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildAsparagine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildAsparagine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                          ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.51, CA, 107.8, N, 180, 0, k + 2, forceField, bondList);
        Atom OD1 = buildHeavy(residue, "OD1", CG, 1.22, CB, 122.5, CA, 180, 0, k + 3, forceField, bondList);
        Atom ND2 = buildHeavy(residue, "ND2", CG, 1.34, CB, 112.7, OD1, 124.0, 1, k + 4, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 107.9, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 107.9, -1, k + 1, forceField, bondList);
        Atom HD21 = buildHydrogen(residue, "HD21", ND2, 1.02, CG, 119.0, CB, 0.0, 0, k + 5, forceField, bondList);
        Atom HD22 = buildHydrogen(residue, "HD22", ND2, 1.02, CG, 119.0, HD21, 120.0, 1, k + 5, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildGlutamate.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildGlutamate(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4, forceField, bondList);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.25, CG, 117.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom OE2 = buildHeavy(residue, "OE2", CD, 1.25, CG, 117.0, OE1, 126.0, 1, k + 5, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildNeutralGlutamicAcid.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildNeutralGlutamicAcid(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                   ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4, forceField, bondList);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.25, CG, 117.0, CB, 180, 0, k + 5, forceField, bondList);
        Atom OE2 = buildHeavy(residue, "OE2", CD, 1.25, CG, 117.0, OE1, 126.0, 1, k + 6, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", OE2, 0.98, CD, 108.7, OE1, 0.0, 0, k + 7, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildGlutamine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildGlutamine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.51, CB, 107.8, CA, 180, 0, k + 4, forceField, bondList);
        Atom OE1 = buildHeavy(residue, "OE1", CD, 1.22, CG, 122.5, CB, 180, 0, k + 5, forceField, bondList);
        Atom NE2 = buildHeavy(residue, "NE2", CD, 1.34, CG, 112.7, OE1, 124.0, 1, k + 6, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 107.9, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 107.9, -1, k + 3, forceField, bondList);
        Atom HE21 = buildHydrogen(residue, "HE21", NE2, 1.02, CD, 119.0, CG, 0.0, 0, k + 7, forceField, bondList);
        Atom HE22 = buildHydrogen(residue, "HE22", NE2, 1.02, CD, 119.0, HE21, 120.0, 1, k + 7, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildMethionine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildMethionine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                          ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom SD = buildHeavy(residue, "SD", CG, 1.82, CB, 109.0, CA, 180, 0, k + 4, forceField, bondList);
        Atom CE = buildHeavy(residue, "CE", SD, 1.82, CG, 96.3, CB, 180, 0, k + 5, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, SD, 112.0, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, SD, 112.0, -1, k + 3, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", CE, 1.11, SD, 112.0, CG, 180.0, 0, k + 6, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, SD, 112.0, HE1, 109.4, 1, k + 6, forceField, bondList);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, SD, 112.0, HE1, 109.4, -1, k + 6, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildLysine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildLysine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                      ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom CE = buildHeavy(residue, "CE", CD, 1.54, CG, 109.5, CB, 180, 0, k + 6, forceField, bondList);
        Atom NZ = buildHeavy(residue, "NZ", CE, 1.50, CD, 109.5, CG, 180, 0, k + 8, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, CE, 109.4, 1, k + 5, forceField, bondList);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, CE, 109.4, -1, k + 5, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, CD, 109.4, NZ, 108.8, 1, k + 7, forceField, bondList);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, CD, 109.4, NZ, 108.8, -1, k + 7, forceField, bondList);
        Atom HZ1 = buildHydrogen(residue, "HZ1", NZ, 1.02, CE, 109.5, CD, 180.0, 0, k + 9, forceField, bondList);
        Atom HZ2 = buildHydrogen(residue, "HZ2", NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, k + 9, forceField, bondList);
        Atom HZ3 = buildHydrogen(residue, "HZ3", NZ, 1.02, CE, 109.5, HZ1, 109.5, -1, k + 9, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDeprotonatedLysine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildDeprotonatedLysine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                                  ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom CE = buildHeavy(residue, "CE", CD, 1.54, CG, 109.5, CB, 180, 0, k + 6, forceField, bondList);
        Atom NZ = buildHeavy(residue, "NZ", CE, 1.50, CD, 109.5, CG, 180, 0, k + 8, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, CE, 109.4, 1, k + 5, forceField, bondList);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, CE, 109.4, -1, k + 5, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", CE, 1.11, CD, 109.4, NZ, 108.8, 1, k + 7, forceField, bondList);
        Atom HE3 = buildHydrogen(residue, "HE3", CE, 1.11, CD, 109.4, NZ, 108.8, -1, k + 7, forceField, bondList);
        Atom HZ1 = buildHydrogen(residue, "HZ1", NZ, 1.02, CE, 109.5, CD, 180.0, 0, k + 9, forceField, bondList);
        Atom HZ2 = buildHydrogen(residue, "HZ2", NZ, 1.02, CE, 109.5, HZ1, 109.5, 1, k + 9, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildArginine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildArginine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                        ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom NE = buildHeavy(residue, "NE", CD, 1.45, CG, 109.5, CB, 180, 0, k + 6, forceField, bondList);
        Atom CZ = buildHeavy(residue, "CZ", NE, 1.35, CD, 120.0, CG, 180, 0, k + 8, forceField, bondList);
        Atom NH1 = buildHeavy(residue, "NH1", CZ, 1.35, NE, 120.0, CD, 180, 0, k + 9, forceField, bondList);
        Atom NH2 = buildHeavy(residue, "NH2", CZ, 1.35, NE, 120.0, NH1, 120.0, 1, k + 9, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, NE, 109.4, 1, k + 5, forceField, bondList);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, NE, 109.4, -1, k + 5, forceField, bondList);
        Atom HE = buildHydrogen(residue, "HE", NE, 1.02, CD, 120.0, CZ, 120.0, 1, k + 7, forceField, bondList);
        Atom HH11 = buildHydrogen(residue, "HH11", NH1, 1.02, CZ, 120.0, NE, 180.0, 0, k + 10, forceField, bondList);
        Atom HH12 = buildHydrogen(residue, "HH12", NH1, 1.02, CZ, 120.0, HH11, 120.0, 1, k + 10, forceField, bondList);
        Atom HH21 = buildHydrogen(residue, "HH21", NH2, 1.02, CZ, 120.0, NE, 180.0, 0, k + 10, forceField, bondList);
        Atom HH22 = buildHydrogen(residue, "HH22", NH2, 1.02, CZ, 120.0, HH21, 120.0, 1, k + 10, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildOrnithine.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildOrnithine(Residue residue, Atom CA, Atom N, Atom C, int k,
                                         ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom NE = buildHeavy(residue, "NE", CD, 1.50, CG, 109.5, CB, 180, 0, k + 6, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        Atom HD2 = buildHydrogen(residue, "HD2", CD, 1.11, CG, 109.4, NE, 109.4, 1, k + 5, forceField, bondList);
        Atom HD3 = buildHydrogen(residue, "HD3", CD, 1.11, CG, 109.4, NE, 109.4, -1, k + 5, forceField, bondList);
        Atom HE1 = buildHydrogen(residue, "HE1", NE, 1.02, CD, 109.5, CG, 180.0, 0, k + 7, forceField, bondList);
        Atom HE2 = buildHydrogen(residue, "HE2", NE, 1.02, CD, 109.5, HE1, 109.5, 1, k + 7, forceField, bondList);
        Atom HE3 = buildHydrogen(residue, "HE3", NE, 1.02, CD, 109.5, HE1, 109.5, -1, k + 7, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildAIB.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildAIB(Residue residue, Atom CA, Atom N, Atom C, int k,
                                   ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB1 = buildHeavy(residue, "CB1", CA, 1.54, N, 109.5, C, 107.8, -1, k, forceField, bondList);
        Atom CB2 = buildHeavy(residue, "CB1", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom HB11 = buildHydrogen(residue, "HB11", CB1, 1.11, CA, 109.4, N, 180.0, 0, k + 1, forceField, bondList);
        Atom HB12 = buildHydrogen(residue, "HB12", CB1, 1.11, CA, 109.4, HB11, 109.4, 1, k + 1, forceField, bondList);
        Atom HB13 = buildHydrogen(residue, "HB13", CB1, 1.11, CA, 109.4, HB11, 109.4, -1, k + 1, forceField, bondList);
        Atom HG21 = buildHydrogen(residue, "HG21", CB2, 1.11, CA, 109.4, N, 180.0, 0, k + 1, forceField, bondList);
        Atom HG22 = buildHydrogen(residue, "HG22", CB2, 1.11, CA, 109.4, HG21, 109.4, 1, k + 1, forceField, bondList);
        Atom HG23 = buildHydrogen(residue, "HG23", CB2, 1.11, CA, 109.4, HG21, 109.4, -1, k + 1, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildPCA.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param CA         a {@link ffx.potential.bonded.Atom} object.
     * @param N          a {@link ffx.potential.bonded.Atom} object.
     * @param C          a {@link ffx.potential.bonded.Atom} object.
     * @param k          a int.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public static Residue buildPCA(Residue residue, Atom CA, Atom N, Atom C, int k,
                                   ForceField forceField, ArrayList<Bond> bondList) {
        Atom CB = buildHeavy(residue, "CB", CA, 1.54, N, 109.5, C, 107.8, 1, k, forceField, bondList);
        Atom CG = buildHeavy(residue, "CG", CB, 1.54, CA, 109.5, N, 180, 0, k + 2, forceField, bondList);
        Atom CD = buildHeavy(residue, "CD", CG, 1.54, CB, 109.5, CA, 180, 0, k + 4, forceField, bondList);
        Atom OE = buildHeavy(residue, "OE", CD, 1.22, N, 126.0, CG, 126.0, 1, k + 5, forceField, bondList);
        Atom HB2 = buildHydrogen(residue, "HB2", CB, 1.11, CA, 109.4, CG, 109.4, 1, k + 1, forceField, bondList);
        Atom HB3 = buildHydrogen(residue, "HB3", CB, 1.11, CA, 109.4, CG, 109.4, -1, k + 1, forceField, bondList);
        Atom HG2 = buildHydrogen(residue, "HG2", CG, 1.11, CB, 109.4, CD, 109.4, 1, k + 3, forceField, bondList);
        Atom HG3 = buildHydrogen(residue, "HG3", CG, 1.11, CB, 109.4, CD, 109.4, -1, k + 3, forceField, bondList);
        return residue;
    }

    /**
     * Check for missing heavy atoms. This check ignores special terminating
     * groups like FOR, NH2, etc.
     *
     * @param aminoAcidNumber a int.
     * @param aminoAcid       a {@link ffx.potential.bonded.ResidueEnumerations.AminoAcid3} object.
     * @param position        a {@link ffx.potential.bonded.Residue.ResiduePosition} object.
     * @param residue         a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static void checkForMissingHeavyAtoms(int aminoAcidNumber, AminoAcid3 aminoAcid,
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
                    logger.warning(format(" Residue %c-%s is missing its N-terminal amide nitrogen", residue.getChainID(), residue));
                    throw e;
                }
                Atom CA = (Atom) residue.getAtomNode("CA");
                if (CA == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("CA", null, null);
                    logger.warning(format(" Residue %c-%s is missing its alpha carbon", residue.getChainID(), residue));
                    throw e;
                }
                Atom C = (Atom) residue.getAtomNode("C");
                if (C == null) {
                    MissingHeavyAtomException e = new MissingHeavyAtomException("C", null, null);
                    logger.warning(format(" Residue %c-%s is missing its C-terminal carboxyl carbon", residue.getChainID(), residue));
                    throw e;
                }
            }
        }
    }

    /**
     * Copy coordinates of an atom from one residue to another.
     *
     * @param fromResidue Use the coordinates from this residue.
     * @param toResidue   Send the coordinates to this residue.
     * @param atomName    The name of the atom whose coordinates will be updated.
     */
    private static void copyCoordinates(Residue fromResidue, Residue toResidue, String atomName) {
        Atom fromAtom;
        if (fromResidue.getAtomNode(atomName) != null) {
            fromAtom = (Atom) fromResidue.getAtomNode(atomName);
        } else {
            fromAtom = (Atom) fromResidue.getAtomNode("H1");
        }

        Atom toAtom = (Atom) toResidue.getAtomNode(atomName);
        toAtom.setXYZ(fromAtom.getXYZ(null));
    }

    /**
     * Constant <code>nCapBackboneAtoms</code>
     */
    public static final String[] nCapBackboneAtoms = {"N", "H1", "H2", "H3", "CA", "HA", "C", "O"};
    /**
     * Constant <code>backboneAtoms</code>
     */
    public static final String[] backboneAtoms = {"N", "H", "CA", "HA", "C", "O"};
    /**
     * Constant <code>glycineBackboneAtoms</code>
     */
    private static final String[] glycineBackboneAtoms = {"N", "H", "CA", "HA2", "HA3", "C", "O"};
    /**
     * Constant <code>prolineBackboneAtoms</code>
     */
    private static final String[] prolineBackboneAtoms = {"N", "CA", "HA", "C", "O"};
    /**
     * Constant <code>alanineAtoms</code>
     */
    private static final String[] alanineAtoms = {"CB", "HB1", "HB2", "HB3"};
    /**
     * Constant <code>glycineAtoms</code>
     */
    private static final String[] glycineAtoms = {"HA2"};
    /**
     * Constant <code>valineAtoms</code>
     */
    private static final String[] valineAtoms = {"CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23"};
    /**
     * Constant <code>leucineAtoms</code>
     */
    private static final String[] leucineAtoms = {"CB", "HB2", "HB3", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23"};
    /**
     * Constant <code>isoleucineAtoms</code>
     */
    private static final String[] isoleucineAtoms = {"CB", "HB", "CG1", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23", "CD1", "HD11", "HD12", "HD13"};
    /**
     * Constant <code>serineAtoms</code>
     */
    private static final String[] serineAtoms = {"CB", "HB2", "HB3", "OG", "HG"};
    /**
     * Constant <code>threonineAtoms</code>
     */
    private static final String[] threonineAtoms = {"CB", "HB", "OG1", "HG1", "CG2", "HG21", "HG22", "HG23"};
    /**
     * Constant <code>cysteineAtoms</code>
     */
    private static final String[] cysteineAtoms = {"CB", "HB2", "HB3", "SG", "HG"};
    /**
     * Constant <code>prolineAtoms</code>
     */
    private static final String[] prolineAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3"};
    /**
     * Constant <code>phenylalanineAtoms</code>
     */
    private static final String[] phenylalanineAtoms = {"CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "HE2", "CZ", "HZ"};
    /**
     * Constant <code>tyrosineAtoms</code>
     */
    private static final String[] tyrosineAtoms = {"CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "HE2", "CZ", "OH", "HH"};
    /**
     * Constant <code>tryptophanAtoms</code>
     */
    private static final String[] tryptophanAtoms = {"CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "NE1", "HE1", "CE2", "CE3", "HE3", "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2"};
    /**
     * Constant <code>histidineAtoms</code>
     */
    private static final String[] histidineAtoms = {"CB", "HB2", "HB3", "CG", "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"};
    /**
     * Constant <code>aspartateAtoms</code>
     */
    private static final String[] aspartateAtoms = {"CB", "HB2", "HB3", "CG", "OD1", "OD2"};
    /**
     * Constant <code>asparagineAtoms</code>
     */
    private static final String[] asparagineAtoms = {"CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"};
    /**
     * Constant <code>glutamateAtoms</code>
     */
    private static final String[] glutamateAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"};
    /**
     * Constant <code>glutamineAtoms</code>
     */
    private static final String[] glutamineAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"};
    /**
     * Constant <code>methionineAtoms</code>
     */
    private static final String[] methionineAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "SD", "CE", "HE1", "HE2", "HE3"};
    /**
     * Constant <code>lysineAtoms</code>
     */
    private static final String[] lysineAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3", "NZ", "HZ1", "HZ2", "HZ3"};
    /**
     * Constant <code>arginineAtoms</code>
     */
    private static final String[] arginineAtoms = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"};

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
    public static final int[][] AA_N = {
            {
                    403, 409, 415, 421, 427, 433, 439, 445,
                    451, 457, 463, 471, 477, 483, 489, 495,
                    501, 507, 513, 519, 525, 531, 537, 543,
                    549, 555, 561, 567, 573, 579, 391, 762,
                    0, 0, 0, 0, 0, 403
            },
            {
                    1, 7, 15, 27, 41, 55, 65, 77,
                    87, 96, 105, 116, 131, 147, 162, 185,
                    202, 218, 234, 244, 256, 268, 280, 294,
                    308, 321, 337, 353, 370, 384, 391, 0,
                    0, 0, 0, 0, 0, 1
            },
            {
                    584, 590, 596, 602, 608, 614, 620, 626,
                    632, 638, 644, 649, 655, 661, 667, 673,
                    679, 685, 691, 697, 703, 709, 715, 721,
                    727, 733, 739, 745, 751, 757, 0, 0,
                    0, 0, 773, 775, 777, 584
            }};

    /**
     * Constant <code>AA_CA</code>
     */
    public static final int[][] AA_CA = {
            {
                    404, 410, 416, 422, 428, 434, 440, 446,
                    452, 458, 464, 472, 478, 484, 490, 496,
                    502, 508, 514, 520, 526, 532, 538, 544,
                    550, 556, 562, 568, 574, 580, 392, 0,
                    0, 767, 0, 0, 0, 404
            },
            {
                    2, 8, 16, 28, 42, 56, 66, 78,
                    88, 97, 106, 117, 132, 148, 163, 186,
                    203, 219, 235, 245, 257, 269, 281, 295,
                    309, 322, 338, 354, 371, 385, 392, 0,
                    0, 0, 0, 0, 0, 2
            },
            {
                    585, 591, 597, 603, 609, 615, 621, 627,
                    633, 639, 645, 650, 656, 662, 668, 674,
                    680, 686, 692, 698, 704, 710, 716, 722,
                    728, 734, 740, 746, 752, 758, 0, 0,
                    0, 0, 0, 0, 779, 585
            }};

    /**
     * Constant <code>AA_C</code>
     */
    public static final int[][] AA_C = {
            {
                    405, 411, 417, 423, 429, 435, 441, 447,
                    453, 459, 465, 473, 479, 485, 491, 497,
                    503, 509, 515, 521, 527, 533, 539, 545,
                    551, 557, 563, 569, 575, 581, 393, 0,
                    764, 769, 0, 0, 0, 405
            },
            {
                    3, 9, 17, 29, 43, 57, 67, 79,
                    89, 98, 107, 118, 133, 149, 164, 187,
                    204, 220, 236, 246, 258, 270, 282, 296,
                    310, 323, 339, 355, 372, 386, 393, 0,
                    0, 0, 0, 0, 0, 3
            },
            {
                    586, 592, 598, 604, 610, 616, 622, 628,
                    634, 640, 646, 651, 657, 663, 669, 675,
                    681, 687, 693, 699, 705, 711, 717, 723,
                    729, 735, 741, 747, 753, 759, 0, 0,
                    0, 0, 771, 0, 0, 586
            }};

    /**
     * Constant <code>AA_HN</code>
     */
    public static final int[][] AA_HN = {
            {
                    406, 412, 418, 424, 430, 436, 442, 448,
                    454, 460, 466, 474, 480, 486, 492, 498,
                    504, 510, 516, 522, 528, 534, 540, 546,
                    552, 558, 564, 570, 576, 582, 394, 763,
                    0, 0, 0, 0, 0, 406
            },
            {
                    4, 10, 18, 30, 44, 58, 68, 80,
                    90, 99, 0, 119, 134, 150, 165, 188,
                    205, 221, 237, 247, 259, 271, 283, 297,
                    311, 324, 340, 356, 373, 387, 394, 0,
                    0, 0, 0, 0, 0, 4
            },
            {
                    587, 593, 599, 605, 611, 617, 623, 629,
                    635, 641, 0, 652, 658, 664, 670, 676,
                    682, 688, 694, 700, 706, 712, 718, 724,
                    730, 736, 742, 748, 754, 760, 0, 0,
                    0, 0, 774, 776, 778, 587
            }};

    /**
     * Constant <code>AA_O</code>
     */
    public static final int[][] AA_O = {
            {
                    407, 413, 419, 425, 431, 437, 443, 449,
                    455, 461, 467, 475, 481, 487, 493, 499,
                    505, 511, 517, 523, 529, 535, 541, 547,
                    553, 559, 565, 571, 577, 583, 395, 0,
                    766, 770, 0, 0, 0, 407
            },
            {
                    5, 11, 19, 31, 45, 59, 69, 81,
                    91, 100, 108, 120, 135, 151, 166, 189,
                    206, 222, 238, 248, 260, 272, 284, 298,
                    312, 325, 341, 357, 374, 388, 395, 0,
                    0, 0, 0, 0, 0, 5
            },
            {
                    588, 594, 600, 606, 612, 618, 624, 630,
                    636, 642, 647, 653, 659, 665, 671, 677,
                    683, 689, 695, 701, 707, 713, 719, 725,
                    731, 737, 743, 749, 755, 761, 0, 0,
                    0, 0, 772, 0, 0, 588
            }};

    /**
     * Constant <code>AA_HA</code>
     */
    public static final int[][] AA_HA = {
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
            }};

    /**
     * Constant <code>AA_CB</code>
     */
    public static final int[] AA_CB = {
            0, 13, 21, 33, 47, 61, 71, 83,
            93, 102, 110, 122, 137, 153, 168, 191,
            208, 224, 240, 250, 262, 274, 286, 300,
            314, 327, 343, 359, 376, 389, 397, 0,
            0, 0, 0, 0, 0, 0
    };

}
