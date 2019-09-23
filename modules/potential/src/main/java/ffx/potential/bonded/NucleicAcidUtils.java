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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.math3.util.FastMath;

import ffx.numerics.math.VectorMath;
import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import static ffx.potential.bonded.BondedUtils.buildBond;
import static ffx.potential.bonded.BondedUtils.buildHeavy;
import static ffx.potential.bonded.BondedUtils.buildHydrogen;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static ffx.potential.bonded.Residue.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.bonded.ResidueEnumerations.nucleicAcidList;

/**
 * Utilities for creating Nucleic Acid residues.
 *
 * @author Michael Schnieders
 * @since 1.0
 */
public class NucleicAcidUtils {

    private static final Logger logger = Logger.getLogger(NucleicAcidUtils.class.getName());

    /**
     * Substitutes from old PDB nucleic acid atom names to new ones.
     */
    private static final Map<String, String> newNucleicAcidAtomNames;

    static {
        Map<String, String> naAtomSubs = new HashMap<>();
        naAtomSubs.put("HO'", "HO2\'");
        naAtomSubs.put("H3T", "HO3\'");
        naAtomSubs.put("H5\'1", "H5\'");
        naAtomSubs.put("H5\'2", "H5\'\'");
        naAtomSubs.put("H5T", "HO5\'");
        naAtomSubs.put("H2\'1", "H2\'");
        naAtomSubs.put("H2\'2", "H2\'\'");
        newNucleicAcidAtomNames = Collections.unmodifiableMap(naAtomSubs);
    }

    /**
     * Assign atom types for a nucleic acid polymer.
     *
     * @param residues   A list of residues that form the nucleic acid polymer.
     * @param forceField The ForceField in use.
     * @param bondList   A list of created bonds.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     * @throws ffx.potential.bonded.BondedUtils.MissingAtomTypeException  if any.
     */
    public static void assignNucleicAcidAtomTypes(List<Residue> residues,
                                                  ForceField forceField, ArrayList<Bond> bondList)
            throws MissingHeavyAtomException, MissingAtomTypeException {

        // A reference to the O3* atom of the previous base.
        Atom pSugarO3 = null;

        // Loop over residues.
        int numberOfResidues = residues.size();
        for (int residueNumber = 0; residueNumber
                < numberOfResidues; residueNumber++) {

            // Match the residue name to a known nucleic acid residue.
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

            // Do atom name conversions.
            List<Atom> resAtoms = residue.getAtomList();
            int natoms = resAtoms.size();
            for (Atom atom : resAtoms) {
                String name = atom.getName();
                name = name.replace('*', '\'');
                name = newNucleicAcidAtomNames.getOrDefault(name, name);
                atom.setName(name);
            }

            // Set a position flag.
            Residue.ResiduePosition position = MIDDLE_RESIDUE;
            boolean lastRes = false;
            if (residueNumber == 0) {
                position = FIRST_RESIDUE;
            }
            if (residueNumber == numberOfResidues - 1) {
                lastRes = true;
                if (position != FIRST_RESIDUE) {
                    position = LAST_RESIDUE;
                }
            }

            // Check if this is a lone 3'-terminal phosphate cap; if so, will skip a lot of logic.
            boolean threePrimeCap = false;
            if (natoms < 7) {
                boolean unrecognized = false;
                for (Atom atom : resAtoms) {
                    String aname = atom.getName();
                    // Recognized: P, OP[1-3], O5' (will be renamed to OP3), HOP[1-3], and DOP[1-3] (deuterated)
                    if (!(aname.equals("P") || aname.matches("[HD]?OP[1-3]") || aname.matches("O5\'"))) {
                        unrecognized = true;
                        break;
                    }
                }
                threePrimeCap = !unrecognized;
                logger.info(" three prime cap " + threePrimeCap);
            }

            // Check if the sugar is deoxyribose and change the residue name if necessary.
            boolean isDNA = false;
            Residue resToCheck = threePrimeCap ? residues.get(residueNumber - 1) : residue;
            Atom sugarO2 = (Atom) resToCheck.getAtomNode("O2\'");

            if (sugarO2 == null) {
                // Assume deoxyribose (DNA) since there is an O2* atom.
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
                // Assume ribose (RNA) since there is no O2* atom.
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

            // Build the phosphate atoms of the current residue.
            Atom phosphate = null;
            Atom sugarO5 = null;
            Atom sugarO3 = null;

            if (threePrimeCap) {
                logger.info(String.format(" EXPERIMENTAL: Adding 3\'-terminal phosphate cap %s", residue));
                Residue priorResidue = residues.get(residueNumber - 1);
                int phosType = isDNA ? 1252 : 1240;
                int oxygenType = isDNA ? 1253 : 1241;

                phosphate = buildHeavy(residue, "P", pSugarO3, phosType, forceField, bondList);
                buildHeavy(residue, "OP1", phosphate, oxygenType, forceField, bondList);
                buildHeavy(residue, "OP2", phosphate, oxygenType, forceField, bondList);
                switch (natoms) {
                    case 3: {
                        buildOP3(residue, phosphate, oxygenType, forceField, bondList, priorResidue, false);
                    }
                    break;
                    case 4: {
                        MSNode bogusO5s = residue.getAtomNode("O5\'");
                        if (bogusO5s != null) {
                            bogusO5s.setName("OP3");
                        }
                        buildHeavy(residue, "OP3", phosphate, oxygenType, forceField, bondList);
                    }
                    break;
                    case 5:
                    case 6: {
                        logger.severe(" Currently, FFX does not support partially-protonated " +
                                "3\'-terminal phosphate caps from PDB files!");
                    }
                }
            } else {
                if (position == FIRST_RESIDUE) {
                    /**
                     * The 5' O5' oxygen of the nucleic acid is generally terminated
                     * by 1.) A phosphate group PO3 (-3). 2.) A hydrogen.
                     *
                     * If the base has phosphate atom we will assume a PO3 group.
                     */
                    phosphate = (Atom) residue.getAtomNode("P");
                    if (phosphate != null) {
                        Residue nextRes = lastRes ? null : residues.get(residueNumber + 1);
                        if (isDNA) {
                            phosphate = buildHeavy(residue, "P", null, 1247, forceField, bondList);
                            buildHeavy(residue, "OP1", phosphate, 1248, forceField, bondList);
                            buildHeavy(residue, "OP2", phosphate, 1248, forceField, bondList);
                            buildOP3(residue, phosphate, 1248, forceField, bondList, nextRes, true);
                            sugarO5 = buildHeavy(residue, "O5\'", phosphate, 1246, forceField, bondList);
                        } else {
                            phosphate = buildHeavy(residue, "P", null, 1235, forceField, bondList);
                            buildHeavy(residue, "OP1", phosphate, 1236, forceField, bondList);
                            buildHeavy(residue, "OP2", phosphate, 1236, forceField, bondList);
                            buildOP3(residue, phosphate, 1236, forceField, bondList, nextRes, true);
                            sugarO5 = buildHeavy(residue, "O5\'", phosphate, 1234, forceField, bondList);
                        }
                    } else if (isDNA) {
                        sugarO5 = buildHeavy(residue, "O5\'", phosphate, 1244, forceField, bondList);
                    } else {
                        sugarO5 = buildHeavy(residue, "O5\'", phosphate, 1232, forceField, bondList);
                    }
                } else {
                    phosphate = buildHeavy(residue, "P", pSugarO3, NA_P[naNumber], forceField, bondList);
                    buildHeavy(residue, "OP1", phosphate, NA_OP[naNumber], forceField, bondList);
                    buildHeavy(residue, "OP2", phosphate, NA_OP[naNumber], forceField, bondList);
                    sugarO5 = buildHeavy(residue, "O5\'", phosphate, NA_O5[naNumber], forceField, bondList);
                }

                // Build the ribose sugar atoms of the current base.
                Atom sugarC5 = buildHeavy(residue, "C5\'", sugarO5, NA_C5[naNumber], forceField, bondList);
                Atom sugarC4 = buildHeavy(residue, "C4\'", sugarC5, NA_C4[naNumber], forceField, bondList);
                Atom sugarO4 = buildHeavy(residue, "O4\'", sugarC4, NA_O4[naNumber], forceField, bondList);
                Atom sugarC1 = buildHeavy(residue, "C1\'", sugarO4, NA_C1[naNumber], forceField, bondList);
                Atom sugarC3 = buildHeavy(residue, "C3\'", sugarC4, NA_C3[naNumber], forceField, bondList);
                Atom sugarC2 = buildHeavy(residue, "C2\'", sugarC3, NA_C2[naNumber], forceField, bondList);
                buildBond(sugarC2, sugarC1, forceField, bondList);
                if (position == LAST_RESIDUE || numberOfResidues == 1) {
                    if (isDNA) {
                        sugarO3 = buildHeavy(residue, "O3\'", sugarC3, 1249, forceField, bondList);
                    } else {
                        sugarO3 = buildHeavy(residue, "O3\'", sugarC3, 1237, forceField, bondList);
                    }
                } else {
                    boolean nextResIsCap = (residues.get(residueNumber + 1).getAtomList().size() < 7);
                    int o3Type = NA_O3[naNumber];
                    if (nextResIsCap) {
                        logger.fine(" Applying a 3\'-terminal-phos-cap O3\' atom type to residue " + residue.toString());
                        o3Type = isDNA ? 1251 : 1239;
                    }
                    sugarO3 = buildHeavy(residue, "O3\'", sugarC3, o3Type, forceField, bondList);
                }
                if (!isDNA) {
                    sugarO2 = buildHeavy(residue, "O2\'", sugarC2, NA_O2[naNumber], forceField, bondList);
                }

                // Build the backbone hydrogen atoms.
                if (position == FIRST_RESIDUE && phosphate == null) {
                    buildHydrogen(residue, "HO5\'", sugarO5, 1.00e0, sugarC5, 109.5e0,
                            sugarC4, 180.0e0, 0, NA_HO5T[naNumber], forceField, bondList);
                }
                buildHydrogen(residue, "H5\'", sugarC5, 1.09e0, sugarO5, 109.5e0,
                        sugarC4, 109.5e0, 1, NA_H51[naNumber], forceField, bondList);
                buildHydrogen(residue, "H5\'\'", sugarC5, 1.09e0, sugarO5, 109.5e0,
                        sugarC4, 109.5e0, -1, NA_H52[naNumber], forceField, bondList);
                buildHydrogen(residue, "H4\'", sugarC4, 1.09e0, sugarC5, 109.5e0,
                        sugarC3, 109.5e0, -1, NA_H4[naNumber], forceField, bondList);
                buildHydrogen(residue, "H3\'", sugarC3, 1.09e0, sugarC4, 109.5e0,
                        sugarC2, 109.5e0, -1, NA_H3[naNumber], forceField, bondList);
                if (isDNA) {
                    buildHydrogen(residue, "H2\'", sugarC2, 1.09e0, sugarC3, 109.5e0,
                            sugarC1, 109.5e0, -1, NA_H21[naNumber], forceField, bondList);
                    buildHydrogen(residue, "H2\'\'", sugarC2, 1.09e0, sugarC3, 109.5e0,
                            sugarC1, 109.5e0, 1, NA_H22[naNumber], forceField, bondList);
                } else {
                    buildHydrogen(residue, "H2\'", sugarC2, 1.09e0, sugarC3, 109.5e0,
                            sugarC1, 109.5e0, -1, NA_H21[naNumber], forceField, bondList);
                    // Add the NA_O2' Methyl for OMC and OMG
                    // Note: may be completely out-of-date with the current Chemical Component Dictionary.
                    if (nucleicAcid == NucleicAcid3.OMC || nucleicAcid == NucleicAcid3.OMG) {
                        Atom CM2 = buildHeavy(residue, "CM2", sugarO2, 1427, forceField, bondList);
                        Atom HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, sugarO2, 109.5e0,
                                sugarC2, 0.0e0, 0, 1428, forceField, bondList);
                        buildHydrogen(residue, "HM22", CM2, 1.08e0, sugarO2, 109.5e0,
                                HM21, 109.5e0, 1, 1429, forceField, bondList);
                        buildHydrogen(residue, "HM23", CM2, 1.08e0, sugarO2, 109.5e0,
                                HM21, 109.5e0, -1, 1430, forceField, bondList);
                    } else {
                        buildHydrogen(residue, "HO2\'", sugarO2, 1.00e0, sugarC2, 109.5e0,
                                sugarC3, 180.0e0, 0, NA_H22[naNumber], forceField, bondList);
                    }
                }
                buildHydrogen(residue, "H1\'", sugarC1, 1.09e0, sugarO4, 109.5e0,
                        sugarC2, 109.5e0, -1, NA_H1[naNumber], forceField, bondList);
                if (position == LAST_RESIDUE || numberOfResidues == 1) {
                    buildHydrogen(residue, "HO3\'", sugarO3, 1.00e0, sugarC3, 109.5e0,
                            sugarC4, 180.0e0, 0, NA_HO3T[naNumber], forceField, bondList);
                    // Else, if it is terminated by a 3' phosphate cap:
                    // Will need to see how PDB would label a 3' phosphate cap.
                }

                // Build the nucleic acid base.
                assignNucleicAcidBaseAtomTypes(nucleicAcid, residue, sugarC1, sugarO4, sugarC2, forceField, bondList);
            }

            // Do some checks on the current base to make sure all atoms have been assigned an atom type.
            resAtoms = residue.getAtomList();
            for (Atom atom : resAtoms) {
                AtomType atomType = atom.getAtomType();
                if (atomType == null) {
                    throw new MissingAtomTypeException(residue, atom);
                }
                int numberOfBonds = atom.getNumBonds();
                if (numberOfBonds != atomType.valence) {
                    if (atom == sugarO3 && numberOfBonds == atomType.valence - 1
                            && position != LAST_RESIDUE && numberOfResidues != 1) {
                        continue;
                    }
                    logger.log(Level.WARNING, format(" An atom for residue %s has the wrong number of bonds:\n %s",
                            residueName, atom.toString()));
                    logger.log(Level.WARNING, format(" Expected: %d Actual: %d.", atomType.valence, numberOfBonds));
                }
            }

            // Save a reference to the current O3* oxygen.
            pSugarO3 = sugarO3;
        }
    }

    /**
     * Assign atom types to the nucleic acid base.
     *
     * @param nucleicAcid The nucleic acid base to use.
     * @param residue     The residue node.
     * @param C1s         The C1* attachment atom.
     * @param O4s         The O4* attachment atom.
     * @param C2s         The C2* attachment atom.
     * @param forceField  The ForceField in use.
     * @param bondList    List of created bonds.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     * @since 1.0
     */
    private static void assignNucleicAcidBaseAtomTypes(
            NucleicAcid3 nucleicAcid, Residue residue, Atom C1s,
            Atom O4s, Atom C2s, ForceField forceField, ArrayList<Bond> bondList)
            throws MissingHeavyAtomException {
        double glyco = 0;
        switch (nucleicAcid) {
            case ADE:
                buildADE(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case M1MA:
                buildM1MA(residue, C1s, forceField, bondList);
                break;
            case CYT:
                buildCYT(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case OMC:
                buildOMC(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case M5MC:
                buildM5MC(residue, C1s, forceField, bondList);
                break;
            case GUA:
                buildGUA(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case OMG:
                buildOMG(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case YYG:
                buildYYG(residue, C1s, forceField, bondList);
                break;
            case M2MG:
                buildM2MG(residue, C1s, forceField, bondList);
                break;
            case M2G:
                buildM2G(residue, C1s, forceField, bondList);
                break;
            case M7MG:
                buildM7MG(residue, C1s, forceField, bondList);
                break;
            case URI:
                buildURI(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case PSU:
                buildPSU(residue, C1s, forceField, bondList);
                break;
            case H2U:
                buildH2U(residue, C1s, forceField, bondList);
                break;
            case M5MU:
                buildM5MU(residue, C1s, forceField, bondList);
                break;
            case DAD:
                buildDAD(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case DCY:
                buildDCY(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case DGU:
                buildDGU(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
            case DTY:
                buildDTY(residue, C1s, O4s, C2s, glyco, forceField, bondList);
                break;
        }
    }

    /**
     * Builds a missing OP3 atom. Can be applied either to 5'-terminal or 3'-terminal phosphates.
     *
     * @param residue         Residue to build an OP3 for.
     * @param phosphate       The phosphate atom.
     * @param aType           The atom type to use.
     * @param forceField      Force field in use.
     * @param bondList        List of bonds.
     * @param adjacentResidue A Residue either 5' or 3' of residue.
     * @param at5prime        If this residue is at the 5' terminus (vs. 3').
     * @throws MissingHeavyAtomException Thrown if a needed heavy atom is missing.
     */
    private static void buildOP3(Residue residue, Atom phosphate, int aType,
                                 ForceField forceField, ArrayList<Bond> bondList, Residue adjacentResidue,
                                 boolean at5prime) throws MissingHeavyAtomException {
        List<Atom> resAts = residue.getAtomList();

        Atom P = null;
        Atom OP1 = null;
        Atom OP2 = null;
        Atom riboOxygen = null;
        Atom riboCarbon = null;

        boolean foundOP3 = false;

        // Label for a break statement: OP3 causes a break from a switch to the surrounding for loop.
        ATLOOP:
        for (Atom at : resAts) {
            String atName = at.getName().toUpperCase();
            switch (atName) {
                case "OP3": {
                    buildHeavy(residue, "OP3", phosphate, aType, forceField, bondList);
                    foundOP3 = true;
                    break ATLOOP;
                }
                case "OP1": {
                    OP1 = at;
                    break;
                }
                case "OP2": {
                    OP2 = at;
                    break;
                }
                case "P": {
                    P = at;
                    break;
                }
                case "O5\'": {
                    riboOxygen = at;
                    break;
                }
                case "C5\'": {
                    riboCarbon = at;
                    break;
                }
            }
        }

        if (!at5prime) {
            riboOxygen = (Atom) adjacentResidue.getAtomNode("O3\'");
            riboCarbon = (Atom) adjacentResidue.getAtomNode("C3\'");
            if (riboCarbon == null || riboOxygen == null) {
                logger.severe(format(" Could not find either O3\' " +
                        "or C3\' in residue %s prior to presumed 3\' " +
                        "phosphate cap %s", adjacentResidue, residue));
            }
        }

        if (!foundOP3) {
            logger.info(format(" EXPERIMENTAL: OP3 of residue %s not found, being rebuilt based on ideal geometry", residue));
            if (P == null || OP1 == null || OP2 == null) {
                throw new IllegalArgumentException(format(" Attempted to build OP3 for residue %s, but one of P, OP1, OP2, O5\', or C5\' were null", residue));
            }

            if (at5prime) {
                if (riboOxygen == null) {
                    logger.fine(" Attempting to find O5\' of subsequent residue");
                    riboOxygen = getNextResAtom(residue, adjacentResidue, "O5\'");
                }
                if (riboCarbon == null) {
                    logger.fine(" Attempting to find C5\' of subsequent residue");
                    riboCarbon = getNextResAtom(residue, adjacentResidue, "C5\'");
                }
            }

            // Borrow the P-OP1 distance.
            double[] xyzOP1 = new double[3];
            double[] xyzP = new double[3];
            xyzOP1 = OP1.getXYZ(xyzOP1);
            xyzP = P.getXYZ(xyzP);
            double distance = VectorMath.dist(xyzP, xyzOP1);

            // Borrow the O5'-P-OP1 angle.
            double[] xyzRiboO = new double[3];
            xyzRiboO = riboOxygen.getXYZ(xyzRiboO);
            double angle = FastMath.toDegrees(VectorMath.bondAngle(xyzRiboO, xyzP, xyzOP1));

            // Borrow the C5'-O5'-P-OP1 dihedral (which will have +120 or +240 degrees added on).
            double[] xyzRiboC = new double[3];
            xyzRiboC = riboCarbon.getXYZ(xyzRiboC);
            double dihedralOP1 = FastMath.toDegrees(VectorMath.dihedralAngle(xyzRiboC, xyzRiboO, xyzP, xyzOP1));

            Atom OP3 = buildHydrogen(residue, "OP3", P, distance, riboOxygen, angle, riboCarbon, dihedralOP1 + 120, 0, aType, forceField, bondList);

            // Measure OP3-OP2 distance for test dihedral + 120 degrees.
            double[] xyzOP2 = new double[3];
            xyzOP2 = OP2.getXYZ(xyzOP2);
            double[] xyzChiral1 = new double[3];
            xyzChiral1 = OP3.getXYZ(xyzChiral1);
            double distChiral1 = VectorMath.dist(xyzChiral1, xyzOP2);

            // Measure OP3-OP2 distance for test dihedral + 240 degrees.
            //intxyz(OP3, P, pOP1, OP1, op1Pop2, O5s, 109.4, 1);
            intxyz(OP3, P, distance, riboOxygen, angle, riboCarbon, dihedralOP1 + 240, 0);
            double[] xyzChiral2 = new double[3];
            xyzChiral2 = OP3.getXYZ(xyzChiral2);
            double distChiral2 = VectorMath.dist(xyzChiral2, xyzOP2);

            logger.fine(String.format(" Bond: %10.5f Angle: %10.5f Dihedral: %10.5f", distance, angle, dihedralOP1));
            logger.fine(String.format(" OP2 position: %s", Arrays.toString(xyzOP2)));
            logger.fine(String.format(" Position 1: %10.6g %10.6g %10.6g with distance %10.6g", xyzChiral1[0], xyzChiral1[1], xyzChiral1[2], distChiral1));
            logger.fine(String.format(" Position 2: %10.6g %10.6g %10.6g with distance %10.6g", xyzChiral2[0], xyzChiral2[1], xyzChiral2[2], distChiral2));
            if (distChiral1 > distChiral2) {
                logger.fine(" Picked dihedral +120");
                OP3.setXYZ(xyzChiral1);
            } else {
                logger.fine(" Picked dihedral +240");
                OP3.setXYZ(xyzChiral2);
            }
        }
    }

    /**
     * Intended for 5'-terminal phosphate caps listed as their own residue: find a given atom name in the next residue.
     *
     * @param residue     Current residue.
     * @param nextResidue Next residue to search in.
     * @param atomName    Atom name to look for.
     * @return The requested atom.
     */
    private static Atom getNextResAtom(Residue residue, Residue nextResidue, String atomName) {
        if (nextResidue == null) {
            throw new IllegalArgumentException(String.format(" Residue %s: No subsequent residue to find atom %s in!", residue, atomName));
        }
        List<Atom> nrAtoms = nextResidue.getAtomList();
        for (Atom atom : nrAtoms) {
            if (atom.getName().equalsIgnoreCase(atomName)) {
                return atom;
            }
        }
        throw new IllegalArgumentException(String.format(" Residue %s: Subsequent residue %s did not contain an atom %s", residue, nextResidue, atomName));
    }

    /**
     * <p>buildADE.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildADE(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1017, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180.0, 0, 1021, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1020, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1019, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1025, forceField, bondList);
        Atom N6 = buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1027, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1024, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1023, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1022, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1018, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1030, forceField, bondList);
        buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1028, forceField, bondList);
        buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1029, forceField, bondList);
        buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1026, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM1MA.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM1MA(Residue residue, Atom C1s, ForceField forceField,
                                     ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1605, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1609, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1608, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1607, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1613, forceField, bondList);
        Atom N6 = buildHeavy(residue, "N6", C6, 1615, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1612, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1611, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1610, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1606, forceField, bondList);
        Atom CM1 = buildHeavy(residue, "CM1", N1, 1619, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1614, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 109.5e0, C4, 180.0e0, 0, 1623, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1618, forceField, bondList);
        buildHydrogen(residue, "HN61", N6, 1.00e0, C6, 109.5e0, C5, 0.0e0, 0, 1616, forceField, bondList);
        buildHydrogen(residue, "HN62", N6, 1.00e0, C6, 109.5e0, C5, 109.5e0, 0, 1617, forceField, bondList);
        Atom HM11 = buildHydrogen(residue, "HM11", CM1, 1.08e0, N1, 109.5e0, C2, 0.0e0, 0, 1620, forceField, bondList);
        buildHydrogen(residue, "HM12", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, 1, 1621, forceField, bondList);
        buildHydrogen(residue, "HM13", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, -1, 1622, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildOMC.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildOMC(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        return buildCYT(residue, C1s, O4s, C2s, glyco, forceField, bondList);
    }

    /**
     * <p>buildCYT.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildCYT(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1078, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco + 180, 0, 1079, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1084, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180.0, 0, 1080, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1081, forceField, bondList);
        Atom N4 = buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, O2, 180.0, 0, 1085, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1082, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1083, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1086, forceField, bondList);
        buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1087, forceField, bondList);
        buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1088, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1089, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM5MC.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM5MC(Residue residue, Atom C1s, ForceField forceField,
                                     ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N1 = buildHeavy(residue, "N1", C1s, 1508, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1509, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1514, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1510, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1511, forceField, bondList);
        Atom N4 = buildHeavy(residue, "N4", C4, 1515, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1512, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1513, forceField, bondList);
        Atom CM5 = buildHeavy(residue, "CM5", C5, 1519, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1516, forceField, bondList);
        buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, C5, 0.0e0, 0, 1517, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1518, forceField, bondList);
        Atom HM51 = buildHydrogen(residue, "HM51", CM5, 1.08e0, C5, 109.5e0, C4, 0.0e0, 0, 1520, forceField, bondList);
        buildHydrogen(residue, "HM52", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, 1, 1521, forceField, bondList);
        buildHydrogen(residue, "HM53", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, -1, 1522, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildOMG.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildOMG(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        return buildGUA(residue, C1s, O4s, C2s, glyco, forceField, bondList);
    }

    /**
     * <p>buildGUA.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildGUA(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1047, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1051, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1050, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1049, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1055, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1060, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1054, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1053, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1057, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1052, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1048, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1061, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1056, forceField, bondList);
        buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1058, forceField, bondList);
        buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1059, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildYYG.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildYYG(Residue residue, Atom C1s, ForceField forceField,
                                    ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1640, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1644, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1643, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1642, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1648, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1650, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1647, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1646, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1649, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1645, forceField, bondList);
        Atom C3 = buildHeavy(residue, "C3", N3, 1652, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1641, forceField, bondList);
        Atom C11 = buildHeavy(residue, "C11", N2, 1657, forceField, bondList);
        Atom C10 = buildHeavy(residue, "C10", C11, 1658, forceField, bondList);
        Atom C12 = buildHeavy(residue, "C12", C11, 1656, forceField, bondList);
        Atom C13 = buildHeavy(residue, "C13", C12, 1662, forceField, bondList);
        Atom C14 = buildHeavy(residue, "C14", C13, 1665, forceField, bondList);
        Atom C15 = buildHeavy(residue, "C15", C14, 1668, forceField, bondList);
        Atom C16 = buildHeavy(residue, "C16", C15, 1675, forceField, bondList);
        Atom O17 = buildHeavy(residue, "O17", C16, 1676, forceField, bondList);
        Atom O18 = buildHeavy(residue, "O18", C16, 1674, forceField, bondList);
        Atom C19 = buildHeavy(residue, "C19", O18, 1670, forceField, bondList);
        Atom N20 = buildHeavy(residue, "N20", C15, 1677, forceField, bondList);
        Atom C21 = buildHeavy(residue, "C21", N20, 1679, forceField, bondList);
        Atom O22 = buildHeavy(residue, "O22", C21, 1680, forceField, bondList);
        Atom O23 = buildHeavy(residue, "O23", C21, 1681, forceField, bondList);
        Atom C24 = buildHeavy(residue, "C24", O23, 1682, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildBond(N1, C12, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1651, forceField, bondList);
        Atom H31 = buildHydrogen(residue, "H31", C3, 1.08e0, N3, 109.5e0, C4, 0.0e0, 0, 1653, forceField, bondList);
        buildHydrogen(residue, "H32", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, 1, 1654, forceField, bondList);
        buildHydrogen(residue, "H33", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, -1, 1655, forceField, bondList);
        Atom H101 = buildHydrogen(residue, "H101", C10, 1.08e0, C11, 109.5e0, N2, 0.0e0, 0, 1659, forceField, bondList);
        buildHydrogen(residue, "H102", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, 1, 1660, forceField, bondList);
        buildHydrogen(residue, "H103", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, -1, 1661, forceField, bondList);
        buildHydrogen(residue, "H131", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, 1, 1663, forceField, bondList);
        buildHydrogen(residue, "H132", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, -1, 1664, forceField, bondList);
        buildHydrogen(residue, "H141", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, 1, 1666, forceField, bondList);
        buildHydrogen(residue, "H142", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, -1, 1667, forceField, bondList);
        buildHydrogen(residue, "H15", C15, 1.08e0, C14, 109.5e0, O18, 180.e0, 0, 1669, forceField, bondList);
        Atom H191 = buildHydrogen(residue, "H191", C19, 1.08e0, O18, 109.5e0, C16, 0.0e0, 0, 1671, forceField, bondList);
        buildHydrogen(residue, "H192", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, 1, 1672, forceField, bondList);
        buildHydrogen(residue, "H193", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, -1, 1673, forceField, bondList);
        buildHydrogen(residue, "HN2", N20, 1.00e0, C15, 109.5e0, O22, 180.0e0, 0, 1678, forceField, bondList);
        Atom H241 = buildHydrogen(residue, "H241", C24, 1.08e0, O23, 109.5e0, C21, 0.0e0, 0, 1683, forceField, bondList);
        buildHydrogen(residue, "H242", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, 1, 1684, forceField, bondList);
        buildHydrogen(residue, "H243", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, -1, 1685, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM2MG.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM2MG(Residue residue, Atom C1s, ForceField forceField,
                                     ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1316, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1320, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1319, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1318, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1324, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1328, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1323, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1322, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1326, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1321, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1317, forceField, bondList);
        Atom CM2 = buildHeavy(residue, "CM2", N2, 1330, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1329, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1325, forceField, bondList);
        buildHydrogen(residue, "H2", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1327, forceField, bondList);
        Atom HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1331, forceField, bondList);
        buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1332, forceField, bondList);
        buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1333, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM2G.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM2G(Residue residue, Atom C1s, ForceField forceField,
                                    ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1379, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1383, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1382, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1381, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1387, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1390, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1386, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1385, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1389, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1384, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1380, forceField, bondList);
        Atom CM1 = buildHeavy(residue, "CM1", N2, 1392, forceField, bondList);
        Atom CM2 = buildHeavy(residue, "CM2", N2, 1396, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1391, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1388, forceField, bondList);
        Atom HM11 = buildHydrogen(residue, "HM11", CM1, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1393, forceField, bondList);
        buildHydrogen(residue, "HM12", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, 1, 1394, forceField, bondList);
        buildHydrogen(residue, "HM13", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, -1, 1395, forceField, bondList);
        Atom HM21 = buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1397, forceField, bondList);
        buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1398, forceField, bondList);
        buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1399, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM7MG.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM7MG(Residue residue, Atom C1s, ForceField forceField,
                                     ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        Atom N9 = buildHeavy(residue, "N9", C1s, 1539, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1543, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1542, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1541, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1547, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1552, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1546, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1545, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1549, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1544, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1540, forceField, bondList);
        Atom CM7 = buildHeavy(residue, "CM7", N7, 1555, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H81", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, 1, 1553, forceField, bondList);
        buildHydrogen(residue, "H82", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, -1, 1554, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1548, forceField, bondList);
        buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1550, forceField, bondList);
        buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N3, 0.0e0, 0, 1551, forceField, bondList);
        Atom HM71 = buildHydrogen(residue, "HM71", CM7, 1.08e0, N7, 109.5e0, C8, 0.0e0, 0, 1556, forceField, bondList);
        buildHydrogen(residue, "HM72", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, 1, 1557, forceField, bondList);
        buildHydrogen(residue, "HM73", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, -1, 1558, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildURI.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildURI(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1106, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.38, C1s, 117.1, O4s, glyco, 0, 1107, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1.22, N1, 123.2, C1s, 0.0, 0, 1112, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.37, N1, 114.8, C1s, 180.0, 0, 1108, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.38, C2, 127.0, N1, 0.0, 0, 1109, forceField, bondList);
        Atom O4 = buildHeavy(residue, "O4", C4, 1.23, N3, 119.8, C2, 180.0, 0, 1114, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1.44, N3, 114.7, C2, 0.0, 0, 1110, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.34, O4, 119.2, C4, 0.0, 0, 1111, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1113, forceField, bondList);
        buildHydrogen(residue, "H5", C5, 1.08e0, C4, 120.4e0, N3, 180.0e0, 0, 1115, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1116, forceField, bondList);

        return residue;
    }

    /**
     * <p>buildPSU.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildPSU(Residue residue, Atom C1s, ForceField forceField,
                                    ArrayList<Bond> bondList) throws MissingHeavyAtomException {
        // C1s bonds to C5 in PsuedoUridine
        Atom C5 = buildHeavy(residue, "C5", C1s, 1485, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1486, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1481, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1482, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1487, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1483, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1484, forceField, bondList);
        Atom O4 = buildHeavy(residue, "O4", C4, 1489, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1491, forceField, bondList);
        buildHydrogen(residue, "H3", N3, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1488, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 120.0e0, C1s, 0.0e0, 0, 1490, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildH2U.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildH2U(Residue residue, Atom C1s, ForceField forceField,
                                    ArrayList<Bond> bondList) throws MissingHeavyAtomException {

        Atom N1 = buildHeavy(residue, "N1", C1s, 1350, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1351, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1356, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1352, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1353, forceField, bondList);
        Atom O4 = buildHeavy(residue, "O4", C4, 1358, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1354, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1355, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1357, forceField, bondList);
        buildHydrogen(residue, "H51", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, 1, 1359, forceField, bondList);
        buildHydrogen(residue, "H52", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, -1, 1360, forceField, bondList);
        buildHydrogen(residue, "H61", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, 1, 1361, forceField, bondList);
        buildHydrogen(residue, "H62", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, -1, 1362, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildM5MU.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException if any.
     */
    private static Residue buildM5MU(Residue residue, Atom C1s, ForceField forceField,
                                     ArrayList<Bond> bondList) throws MissingHeavyAtomException {

        Atom N1 = buildHeavy(residue, "N1", C1s, 1575, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1576, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1581, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1577, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1578, forceField, bondList);
        Atom O4 = buildHeavy(residue, "O4", C4, 1583, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1579, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1580, forceField, bondList);
        Atom C5M = buildHeavy(residue, "C5M", C5, 1585, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1582, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1584, forceField, bondList);
        Atom H5M1 = buildHydrogen(residue, "H5M1", C5M, 1.08e0, C5, 109.5e0, C6, 0.0e0, 0, 1586, forceField, bondList);
        buildHydrogen(residue, "H5M2", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, 1, 1587, forceField, bondList);
        buildHydrogen(residue, "H5M3", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, -1, 1588, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDAD.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildDAD(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {

        Atom N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1132, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180, 0, 1136, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1135, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1134, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1140, forceField, bondList);
        Atom N6 = buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1142, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1139, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1138, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1137, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1133, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1145, forceField, bondList);
        buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1143, forceField, bondList);
        buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1144, forceField, bondList);
        buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1141, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDCY.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildDCY(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {

        Atom N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1191, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco, 0, 1192, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1197, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180, 0, 1193, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1194, forceField, bondList);
        Atom N4 = buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, C2, 180.0, 0, 1198, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1195, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1196, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1199, forceField, bondList);
        buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1200, forceField, bondList);
        buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1201, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1202, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDGU.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildDGU(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {

        Atom N9 = buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1161, forceField, bondList);
        Atom C8 = buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1165, forceField, bondList);
        Atom N7 = buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1164, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1163, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1169, forceField, bondList);
        Atom O6 = buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1174, forceField, bondList);
        Atom N1 = buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1168, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1167, forceField, bondList);
        Atom N2 = buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1171, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1166, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1162, forceField, bondList);
        buildBond(C4, C5, forceField, bondList);
        buildBond(C4, N9, forceField, bondList);
        buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1175, forceField, bondList);
        buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1170, forceField, bondList);
        buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1172, forceField, bondList);
        buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1173, forceField, bondList);
        return residue;
    }

    /**
     * <p>buildDTY.</p>
     *
     * @param residue    a {@link ffx.potential.bonded.Residue} object.
     * @param C1s        a {@link ffx.potential.bonded.Atom} object.
     * @param O4s        a {@link ffx.potential.bonded.Atom} object.
     * @param C2s        a {@link ffx.potential.bonded.Atom} object.
     * @param glyco      a double.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param bondList   a {@link java.util.ArrayList} object.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    private static Residue buildDTY(Residue residue, Atom C1s, Atom O4s, Atom C2s, double glyco,
                                    ForceField forceField, ArrayList<Bond> bondList) {
        Atom N1 = buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1218, forceField, bondList);
        Atom C2 = buildHeavy(residue, "C2", N1, 1.37, C1s, 117.1, O4s, glyco, 0, 1219, forceField, bondList);
        Atom O2 = buildHeavy(residue, "O2", C2, 1.22, N1, 122.9, C1s, 0.0, 0, 1224, forceField, bondList);
        Atom N3 = buildHeavy(residue, "N3", C2, 1.38, N1, 115.4, C1s, 180.0, 0, 1220, forceField, bondList);
        Atom C4 = buildHeavy(residue, "C4", N3, 1.38, C2, 126.4, N1, 0.0, 0, 1221, forceField, bondList);
        Atom O4 = buildHeavy(residue, "O4", C4, 1.23, N3, 120.5, C2, 180.0, 0, 1226, forceField, bondList);
        Atom C5 = buildHeavy(residue, "C5", C4, 1.44, N3, 114.1, C2, 0.0, 0, 1222, forceField, bondList);
        Atom C7 = buildHeavy(residue, "C7", C5, 1.50, C4, 117.5, N3, 180.0, 0, 1227, forceField, bondList);
        Atom C6 = buildHeavy(residue, "C6", C5, 1.34, C4, 120.8, N3, 0.0, 0, 1223, forceField, bondList);
        buildBond(C6, N1, forceField, bondList);
        buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.8e0, N1, 180.0e0, 0, 1225, forceField, bondList);
        Atom H = buildHydrogen(residue, "H71", C7, 1.09e0, C5, 109.5e0, C4, 0.0e0, 0, 1228, forceField, bondList);
        buildHydrogen(residue, "H72", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, 1, 1228, forceField, bondList);
        buildHydrogen(residue, "H73", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, -1, 1228, forceField, bondList);
        buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1229, forceField, bondList);
        return residue;
    }

    /**
     * Biotype keys for nucleic acid backbone atom types. These are consistent
     * with parameter files from TINKER v. 6.1 (June 2012).
     */
    public static final int[] NA_O5 = {
            1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0,
            1300, 1334, 1363, 1400, 1431, 1465, 1492, 1523, 1559, 1589, 1624
    };
    /**
     * Constant <code>NA_C5</code>
     */
    public static final int[] NA_C5 = {
            1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0,
            1301, 1335, 1364, 1401, 1432, 1466, 1493, 1524, 1560, 1590, 1625
    };
    /**
     * Constant <code>NA_H51</code>
     */
    public static final int[] NA_H51 = {
            1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0,
            1302, 1336, 1365, 1402, 1433, 1467, 1494, 1525, 1561, 1591, 1626
    };
    /**
     * Constant <code>NA_H52</code>
     */
    public static final int[] NA_H52 = {
            1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0,
            1303, 1337, 1366, 1403, 1434, 1468, 1495, 1526, 1562, 1592, 1627
    };
    /**
     * Constant <code>NA_C4</code>
     */
    public static final int[] NA_C4 = {
            1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0,
            1304, 1338, 1367, 1404, 1435, 1469, 1496, 1527, 1563, 1593, 1628
    };
    /**
     * Constant <code>NA_H4</code>
     */
    public static final int[] NA_H4 = {
            1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0,
            1305, 1339, 1368, 1405, 1436, 1470, 1497, 1528, 1564, 1594, 1629
    };
    /**
     * Constant <code>NA_O4</code>
     */
    public static final int[] NA_O4 = {
            1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0,
            1306, 1340, 1369, 1406, 1437, 1471, 1498, 1529, 1565, 1595, 1630
    };
    /**
     * Constant <code>NA_C1</code>
     */
    public static final int[] NA_C1 = {
            1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0,
            1307, 1341, 1370, 1407, 1438, 1472, 1499, 1530, 1566, 1596, 1631
    };
    /**
     * Constant <code>NA_H1</code>
     */
    public static final int[] NA_H1 = {
            1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0,
            1308, 1342, 1371, 1408, 1439, 1473, 1500, 1531, 1567, 1597, 1632
    };
    /**
     * Constant <code>NA_C3</code>
     */
    public static final int[] NA_C3 = {
            1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0,
            1309, 1343, 1372, 1409, 1440, 1474, 1501, 1532, 1568, 1598, 1633
    };
    /**
     * Constant <code>NA_H3</code>
     */
    public static final int[] NA_H3 = {
            1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0,
            1310, 1344, 1373, 1410, 1441, 1475, 1502, 1533, 1569, 1599, 1634
    };
    /**
     * Constant <code>NA_C2</code>
     */
    public static final int[] NA_C2 = {
            1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0,
            1311, 1345, 1374, 1411, 1442, 1476, 1503, 1534, 1570, 1600, 1635
    };
    /**
     * Constant <code>NA_H21</code>
     */
    public static final int[] NA_H21 = {
            1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0,
            1312, 1346, 1375, 1412, 1443, 1477, 1504, 1535, 1571, 1601, 1636
    };
    /**
     * Constant <code>NA_O2</code>
     */
    public static final int[] NA_O2 = {
            1014, 1044, 1075, 1103, 0, 0, 0, 0, 0, 0, 0, 0,
            1313, 1347, 1376, 1413, 1444, 1478, 1505, 1536, 1572, 1602, 1637
    };
    /**
     * Constant <code>NA_H22</code>
     */
    public static final int[] NA_H22 = {
            1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0,
            1314, 1348, 1377, 0, 0, 1479, 1506, 1537, 1573, 1603, 1638
    };
    /**
     * Constant <code>NA_O3</code>
     */
    public static final int[] NA_O3 = {
            1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0,
            1315, 1349, 1378, 1414, 1445, 1480, 1507, 1538, 1574, 1604, 1639
    };
    /**
     * Constant <code>NA_P</code>
     */
    public static final int[] NA_P = {
            1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0,
            1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230
    };
    /**
     * Constant <code>NA_OP</code>
     */
    public static final int[] NA_OP = {
            1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0,
            1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231
    };
    /**
     * Should be NA_HO5' (' replaced by T in the name).
     *
     * Constant <code>NA_HO5T</code>
     */
    public static final int[] NA_HO5T = {
            1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0,
            1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233
    };
    /**
     * Should be NA_HO3' (' replaced by T in the name).
     *
     * Constant <code>NA_HO3T</code>
     */
    public static final int[] NA_HO3T = {
            1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0,
            1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238
    };

}
