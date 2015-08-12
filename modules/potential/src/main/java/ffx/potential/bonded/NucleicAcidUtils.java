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
package ffx.potential.bonded;

import static ffx.potential.bonded.Residue.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.bonded.ResidueEnumerations.nucleicAcidList;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import static java.lang.String.format;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utilities for creating Nucleic Acid residues.
 *
 * @author Michael Schnieders
 */
public class NucleicAcidUtils {
    private static final Logger logger = Logger.getLogger(NucleicAcidUtils.class.getName());

    /**
     * Biotype keys for nucleic acid backbone atom types. These are consistent
     * with parameter files from TINKER v. 6.1 (June 2012).
     */
    public static final int o5Typ[] = {
        1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0,
        1300, 1334, 1363, 1400, 1431, 1465, 1492, 1523, 1559, 1589, 1624
    };
    public static final int c5Typ[] = {
        1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0,
        1301, 1335, 1364, 1401, 1432, 1466, 1493, 1524, 1560, 1590, 1625
    };
    public static final int h51Typ[] = {
        1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0,
        1302, 1336, 1365, 1402, 1433, 1467, 1494, 1525, 1561, 1591, 1626
    };
    public static final int h52Typ[] = {
        1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0,
        1303, 1337, 1366, 1403, 1434, 1468, 1495, 1526, 1562, 1592, 1627
    };
    public static final int c4Typ[] = {
        1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0,
        1304, 1338, 1367, 1404, 1435, 1469, 1496, 1527, 1563, 1593, 1628
    };
    public static final int h4Typ[] = {
        1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0,
        1305, 1339, 1368, 1405, 1436, 1470, 1497, 1528, 1564, 1594, 1629
    };
    public static final int o4Typ[] = {
        1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0,
        1306, 1340, 1369, 1406, 1437, 1471, 1498, 1529, 1565, 1595, 1630
    };
    public static final int c1Typ[] = {
        1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0,
        1307, 1341, 1370, 1407, 1438, 1472, 1499, 1530, 1566, 1596, 1631
    };
    public static final int h1Typ[] = {
        1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0,
        1308, 1342, 1371, 1408, 1439, 1473, 1500, 1531, 1567, 1597, 1632
    };
    public static final int c3Typ[] = {
        1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0,
        1309, 1343, 1372, 1409, 1440, 1474, 1501, 1532, 1568, 1598, 1633
    };
    public static final int h3Typ[] = {
        1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0,
        1310, 1344, 1373, 1410, 1441, 1475, 1502, 1533, 1569, 1599, 1634
    };
    public static final int c2Typ[] = {
        1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0,
        1311, 1345, 1374, 1411, 1442, 1476, 1503, 1534, 1570, 1600, 1635
    };
    public static final int h21Typ[] = {
        1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0,
        1312, 1346, 1375, 1412, 1443, 1477, 1504, 1535, 1571, 1601, 1636
    };
    public static final int o2Typ[] = {
        1014, 1044, 1075, 1103, 0, 0, 0, 0, 0, 0, 0, 0,
        1313, 1347, 1376, 1413, 1444, 1478, 1505, 1536, 1572, 1602, 1637
    };
    public static final int h22Typ[] = {
        1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0,
        1314, 1348, 1377, 0, 0, 1479, 1506, 1537, 1573, 1603, 1638
    };
    public static final int o3Typ[] = {
        1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0,
        1315, 1349, 1378, 1414, 1445, 1480, 1507, 1538, 1574, 1604, 1639
    };
    public static final int pTyp[] = {
        1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0,
        1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230
    };
    public static final int opTyp[] = {
        1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0,
        1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231
    };
    public static final int h5tTyp[] = {
        1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0,
        1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233
    };
    public static final int h3tTyp[] = {
        1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0,
        1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238
    };

    /**
     * Assign atom types for a nucleic acid polymer.
     *
     * @param residues
     * @param forceField
     * @param bondList
     * @throws ffx.potential.bonded.BondedUtils.MissingHeavyAtomException
     * @throws ffx.potential.bonded.BondedUtils.MissingAtomTypeException
     */
    public static void assignNucleicAcidAtomTypes(List<Residue> residues, ForceField forceField, ArrayList<Bond> bondList)
            throws BondedUtils.MissingHeavyAtomException, BondedUtils.MissingAtomTypeException {
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
            ResidueEnumerations.NucleicAcid3 nucleicAcid = null;
            int naNumber = -1;
            for (ResidueEnumerations.NucleicAcid3 nucleic : nucleicAcidList) {
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
             P3s = BondedUtils.buildHeavy(prevResidue, "P3s", null, 1247, forceField, bondList);
             } else {
             P3s = BondedUtils.buildHeavy(prevResidue, "P3s", null, 1235, forceField, bondList);
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
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.DAD;
                            residueName = "DAD";
                            residue.setName(residueName);
                            break;
                        case CYT:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.DCY;
                            residueName = "DCY";
                            residue.setName(residueName);
                            break;
                        case GUA:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.DGU;
                            residueName = "DGU";
                            residue.setName(residueName);
                            break;
                        case THY:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.DTY;
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
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.ADE;
                            residueName = "ADE";
                            residue.setName(residueName);
                            break;
                        case DCY:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.CYT;
                            residueName = "CYT";
                            residue.setName(residueName);
                            break;
                        case DGU:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.GUA;
                            residueName = "GUA";
                            residue.setName(residueName);
                            break;
                        case DTY:
                            nucleicAcid = ResidueEnumerations.NucleicAcid3.THY;
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
            Residue.ResiduePosition position = MIDDLE_RESIDUE;
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
                        P = BondedUtils.buildHeavy(residue, "P", null, 1247, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP1", P, 1248, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP2", P, 1248, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP3", P, 1248, forceField, bondList);
                        O5s = BondedUtils.buildHeavy(residue, "O5\'", P, 1246, forceField, bondList);
                    } else {
                        P = BondedUtils.buildHeavy(residue, "P", null, 1235, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP1", P, 1236, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP2", P, 1236, forceField, bondList);
                        BondedUtils.buildHeavy(residue, "OP3", P, 1236, forceField, bondList);
                        O5s = BondedUtils.buildHeavy(residue, "O5\'", P, 1234, forceField, bondList);
                    }
                } else {
                    if (isDNA) {
                        O5s = BondedUtils.buildHeavy(residue, "O5\'", P, 1244, forceField, bondList);
                    } else {
                        O5s = BondedUtils.buildHeavy(residue, "O5\'", P, 1232, forceField, bondList);
                    }
                }
            } else {
                P = BondedUtils.buildHeavy(residue, "P", pO3s, pTyp[naNumber], forceField, bondList);
                BondedUtils.buildHeavy(residue, "OP1", P, opTyp[naNumber], forceField, bondList);
                BondedUtils.buildHeavy(residue, "OP2", P, opTyp[naNumber], forceField, bondList);
                O5s = BondedUtils.buildHeavy(residue, "O5\'", P, o5Typ[naNumber], forceField, bondList);
            }
            /**
             * Build the ribose sugar atoms of the current base.
             */
            Atom C5s = BondedUtils.buildHeavy(residue, "C5\'", O5s, c5Typ[naNumber], forceField, bondList);
            Atom C4s = BondedUtils.buildHeavy(residue, "C4\'", C5s, c4Typ[naNumber], forceField, bondList);
            Atom O4s = BondedUtils.buildHeavy(residue, "O4\'", C4s, o4Typ[naNumber], forceField, bondList);
            Atom C1s = BondedUtils.buildHeavy(residue, "C1\'", O4s, c1Typ[naNumber], forceField, bondList);
            Atom C3s = BondedUtils.buildHeavy(residue, "C3\'", C4s, c3Typ[naNumber], forceField, bondList);
            Atom C2s = BondedUtils.buildHeavy(residue, "C2\'", C3s, c2Typ[naNumber], forceField, bondList);
            BondedUtils.buildBond(C2s, C1s, forceField, bondList);
            Atom O3s = null;
            if (position == LAST_RESIDUE || numberOfResidues == 1) {
                if (isDNA) {
                    O3s = BondedUtils.buildHeavy(residue, "O3\'", C3s, 1249, forceField, bondList);
                } else {
                    O3s = BondedUtils.buildHeavy(residue, "O3\'", C3s, 1237, forceField, bondList);
                }
            } else {
                O3s = BondedUtils.buildHeavy(residue, "O3\'", C3s, o3Typ[naNumber], forceField, bondList);
            }
            if (!isDNA) {
                O2s = BondedUtils.buildHeavy(residue, "O2\'", C2s, o2Typ[naNumber], forceField, bondList);
            }
            /**
             * Build the backbone hydrogen atoms.
             */
            if (position == FIRST_RESIDUE && P == null) {
                BondedUtils.buildHydrogen(residue, "H5T", O5s, 1.00e0, C5s, 109.5e0, C4s, 180.0e0, 0, h5tTyp[naNumber], forceField, bondList);
            }
            BondedUtils.buildHydrogen(residue, "H5\'1", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, 1, h51Typ[naNumber], forceField, bondList);
            BondedUtils.buildHydrogen(residue, "H5\'2", C5s, 1.09e0, O5s, 109.5e0, C4s, 109.5e0, -1, h52Typ[naNumber], forceField, bondList);
            BondedUtils.buildHydrogen(residue, "H4\'", C4s, 1.09e0, C5s, 109.5e0, C3s, 109.5e0, -1, h4Typ[naNumber], forceField, bondList);
            BondedUtils.buildHydrogen(residue, "H3\'", C3s, 1.09e0, C4s, 109.5e0, C2s, 109.5e0, -1, h3Typ[naNumber], forceField, bondList);
            if (isDNA) {
                BondedUtils.buildHydrogen(residue, "H2\'1", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber], forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H2\'2", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, 1, h22Typ[naNumber], forceField, bondList);
            } else {
                BondedUtils.buildHydrogen(residue, "H2\'", C2s, 1.09e0, C3s, 109.5e0, C1s, 109.5e0, -1, h21Typ[naNumber], forceField, bondList);
                // Add the O2' Methyl for OMC and OMG
                if (nucleicAcid == ResidueEnumerations.NucleicAcid3.OMC || nucleicAcid == ResidueEnumerations.NucleicAcid3.OMG) {
                    Atom CM2 = BondedUtils.buildHeavy(residue, "CM2", O2s, 1427, forceField, bondList);
                    Atom HM21 = BondedUtils.buildHydrogen(residue, "HM21", CM2, 1.08e0, O2s, 109.5e0, C2s, 0.0e0, 0, 1428, forceField, bondList);
                    BondedUtils.buildHydrogen(residue, "HM22", CM2, 1.08e0, O2s, 109.5e0, HM21, 109.5e0, 1, 1429, forceField, bondList);
                    BondedUtils.buildHydrogen(residue, "HM23", CM2, 1.08e0, O2s, 109.5e0, HM21, 109.5e0, -1, 1430, forceField, bondList);
                } else {
                    BondedUtils.buildHydrogen(residue, "HO\'", O2s, 1.00e0, C2s, 109.5e0, C3s, 180.0e0, 0, h22Typ[naNumber], forceField, bondList);
                }
            }
            BondedUtils.buildHydrogen(residue, "H1\'", C1s, 1.09e0, O4s, 109.5e0, C2s, 109.5e0, -1, h1Typ[naNumber], forceField, bondList);
            if (position == LAST_RESIDUE || numberOfResidues == 1) {
                BondedUtils.buildHydrogen(residue, "H3T", O3s, 1.00e0, C3s, 109.5e0, C4s, 180.0e0, 0, h3tTyp[naNumber], forceField, bondList);
                // Else, if it is terminated by a 3' phosphate cap:
                // Will need to see how PDB would label a 3' phosphate cap.
            }
            /**
             * Build the nucleic acid base.
             */
            try {
                assignNucleicAcidBaseAtomTypes(nucleicAcid, residue, C1s, O4s, C2s, forceField, bondList);
            } catch (BondedUtils.MissingHeavyAtomException missingHeavyAtomException) {
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
                    String protonEq = atom.getName().replaceFirst("D", "H");
                    Atom correspH = (Atom) residue.getAtomNode(protonEq, true);
                    if (correspH == null || correspH.getAtomType() == null) {
                        MissingAtomTypeException missingAtomTypeException
                                = new MissingAtomTypeException(residue, atom);
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
    private static void assignNucleicAcidBaseAtomTypes(ResidueEnumerations.NucleicAcid3 nucleicAcid, Residue residue, Atom C1s,
            Atom O4s, Atom C2s, ForceField forceField, ArrayList<Bond> bondList)
            throws BondedUtils.MissingHeavyAtomException {
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
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1017, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180, 0, 1021, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1020, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1019, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1025, forceField, bondList);
                N6 = BondedUtils.buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1027, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1024, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1023, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1022, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1018, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1030, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1028, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1029, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1026, forceField, bondList);
                break;
            case M1MA:
                Atom CM1;
                Atom HM11;
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1605, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1609, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1608, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1607, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1613, forceField, bondList);
                N6 = BondedUtils.buildHeavy(residue, "N6", C6, 1615, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1612, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1611, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1610, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1606, forceField, bondList);
                CM1 = BondedUtils.buildHeavy(residue, "CM1", N1, 1619, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1614, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 109.5e0, C4, 180.0e0, 0, 1623, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1618, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HN61", N6, 1.00e0, C6, 109.5e0, C5, 0.0e0, 0, 1616, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HN62", N6, 1.00e0, C6, 109.5e0, C5, 109.5e0, 0, 1617, forceField, bondList);
                HM11 = BondedUtils.buildHydrogen(residue, "HM11", CM1, 1.08e0, N1, 109.5e0, C2, 0.0e0, 0, 1620, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM12", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, 1, 1621, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM13", CM1, 1.08e0, N1, 109.5e0, HM11, 109.5e0, -1, 1622, forceField, bondList);
                break;
            case CYT:
            case OMC:
                Atom O2;
                Atom N4;
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1078, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco + 180, 0, 1079, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1084, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180.0, 0, 1080, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1081, forceField, bondList);
                N4 = BondedUtils.buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, O2, 180.0, 0, 1085, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1082, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1083, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1086, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1087, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1088, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1089, forceField, bondList);
                break;
            case M5MC:
                Atom CM5;
                Atom HM51;
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1508, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1509, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1514, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1510, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1511, forceField, bondList);
                N4 = BondedUtils.buildHeavy(residue, "N4", C4, 1515, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1512, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1513, forceField, bondList);
                CM5 = BondedUtils.buildHeavy(residue, "CM5", C5, 1519, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1516, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, C5, 0.0e0, 0, 1517, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1518, forceField, bondList);
                HM51 = BondedUtils.buildHydrogen(residue, "HM51", CM5, 1.08e0, C5, 109.5e0, C4, 0.0e0, 0, 1520, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM52", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, 1, 1521, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM53", CM5, 1.08e0, C5, 109.5e0, HM51, 109.5e0, -1, 1522, forceField, bondList);
                break;
            case GUA:
            case OMG:
                Atom O6;
                Atom N2;
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1047, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1051, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1050, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1049, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1055, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1060, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1054, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1053, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1057, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1052, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1048, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1061, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1056, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1058, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1059, forceField, bondList);
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
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1640, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1644, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1643, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1642, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1648, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1650, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1647, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1646, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1649, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1645, forceField, bondList);
                C3 = BondedUtils.buildHeavy(residue, "C3", N3, 1652, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1641, forceField, bondList);
                C11 = BondedUtils.buildHeavy(residue, "C11", N2, 1657, forceField, bondList);
                C10 = BondedUtils.buildHeavy(residue, "C10", C11, 1658, forceField, bondList);
                C12 = BondedUtils.buildHeavy(residue, "C12", C11, 1656, forceField, bondList);
                C13 = BondedUtils.buildHeavy(residue, "C13", C12, 1662, forceField, bondList);
                C14 = BondedUtils.buildHeavy(residue, "C14", C13, 1665, forceField, bondList);
                C15 = BondedUtils.buildHeavy(residue, "C15", C14, 1668, forceField, bondList);
                C16 = BondedUtils.buildHeavy(residue, "C16", C15, 1675, forceField, bondList);
                O17 = BondedUtils.buildHeavy(residue, "O17", C16, 1676, forceField, bondList);
                O18 = BondedUtils.buildHeavy(residue, "O18", C16, 1674, forceField, bondList);
                C19 = BondedUtils.buildHeavy(residue, "C19", O18, 1670, forceField, bondList);
                N20 = BondedUtils.buildHeavy(residue, "N20", C15, 1677, forceField, bondList);
                C21 = BondedUtils.buildHeavy(residue, "C21", N20, 1679, forceField, bondList);
                O22 = BondedUtils.buildHeavy(residue, "O22", C21, 1680, forceField, bondList);
                O23 = BondedUtils.buildHeavy(residue, "O23", C21, 1681, forceField, bondList);
                C24 = BondedUtils.buildHeavy(residue, "C24", O23, 1682, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildBond(N1, C12, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1651, forceField, bondList);
                H31 = BondedUtils.buildHydrogen(residue, "H31", C3, 1.08e0, N3, 109.5e0, C4, 0.0e0, 0, 1653, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H32", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, 1, 1654, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H33", C3, 1.08e0, N3, 109.5e0, H31, 109.5e0, -1, 1655, forceField, bondList);
                H101 = BondedUtils.buildHydrogen(residue, "H101", C10, 1.08e0, C11, 109.5e0, N2, 0.0e0, 0, 1659, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H102", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, 1, 1660, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H103", C10, 1.08e0, C11, 109.5e0, H101, 109.5e0, -1, 1661, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H131", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, 1, 1663, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H132", C13, 1.08e0, C12, 109.5e0, C14, 109.5e0, -1, 1664, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H141", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, 1, 1666, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H142", C14, 1.08e0, C13, 109.5e0, C15, 109.5e0, -1, 1667, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H15", C15, 1.08e0, C14, 109.5e0, O18, 180.e0, 0, 1669, forceField, bondList);
                H191 = BondedUtils.buildHydrogen(residue, "H191", C19, 1.08e0, O18, 109.5e0, C16, 0.0e0, 0, 1671, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H192", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, 1, 1672, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H193", C19, 1.08e0, O18, 109.5e0, H191, 109.5e0, -1, 1673, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HN2", N20, 1.00e0, C15, 109.5e0, O22, 180.0e0, 0, 1678, forceField, bondList);
                H241 = BondedUtils.buildHydrogen(residue, "H241", C24, 1.08e0, O23, 109.5e0, C21, 0.0e0, 0, 1683, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H242", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, 1, 1684, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H243", C24, 1.08e0, O23, 109.5e0, H241, 109.5e0, -1, 1685, forceField, bondList);
                break;
            case M2MG:
                Atom CM2;
                Atom HM21;
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1316, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1320, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1319, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1318, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1324, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1328, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1323, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1322, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1326, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1321, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1317, forceField, bondList);
                CM2 = BondedUtils.buildHeavy(residue, "CM2", N2, 1330, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1329, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1325, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H2", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1327, forceField, bondList);
                HM21 = BondedUtils.buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1331, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1332, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1333, forceField, bondList);
                break;
            case M2G:
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1379, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1383, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1382, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1381, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1387, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1390, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1386, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1385, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1389, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1384, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1380, forceField, bondList);
                CM1 = BondedUtils.buildHeavy(residue, "CM1", N2, 1392, forceField, bondList);
                CM2 = BondedUtils.buildHeavy(residue, "CM2", N2, 1396, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1391, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1388, forceField, bondList);
                HM11 = BondedUtils.buildHydrogen(residue, "HM11", CM1, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1393, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM12", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, 1, 1394, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM13", CM1, 1.08e0, N2, 109.5e0, HM11, 109.5e0, -1, 1395, forceField, bondList);
                HM21 = BondedUtils.buildHydrogen(residue, "HM21", CM2, 1.08e0, N2, 109.5e0, C2, 0.0e0, 0, 1397, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM22", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, 1, 1398, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM23", CM2, 1.08e0, N2, 109.5e0, HM21, 109.5e0, -1, 1399, forceField, bondList);
                break;
            case M7MG:
                Atom CM7;
                Atom HM71;
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1539, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1543, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1542, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1541, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1547, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1552, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1546, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1545, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1549, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1544, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1540, forceField, bondList);
                CM7 = BondedUtils.buildHeavy(residue, "CM7", N7, 1555, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H81", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, 1, 1553, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H82", C8, 1.08e0, N7, 109.5e0, N9, 109.5e0, -1, 1554, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1548, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1550, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N3, 0.0e0, 0, 1551, forceField, bondList);
                HM71 = BondedUtils.buildHydrogen(residue, "HM71", CM7, 1.08e0, N7, 109.5e0, C8, 0.0e0, 0, 1556, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM72", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, 1, 1557, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "HM73", CM7, 1.08e0, N7, 109.5e0, HM71, 109.5e0, -1, 1558, forceField, bondList);
                break;
            case URI:
                Atom O4;
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1106, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.38, C1s, 117.1, O4s, glyco, 0, 1107, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1.22, N1, 123.2, C1s, 0.0, 0, 1112, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.37, N1, 114.8, C1s, 180.0, 0, 1108, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.38, C2, 127.0, N1, 0.0, 0, 1109, forceField, bondList);
                O4 = BondedUtils.buildHeavy(residue, "O4", C4, 1.23, N3, 119.8, C2, 180.0, 0, 1114, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1.44, N3, 114.7, C2, 0.0, 0, 1110, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.34, O4, 119.2, C4, 0.0, 0, 1111, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1113, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H5", C5, 1.08e0, C4, 120.4e0, N3, 180.0e0, 0, 1115, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1116, forceField, bondList);
                break;
            case PSU:
                // C1s bonds to C5 in PsuedoUridine
                C5 = BondedUtils.buildHeavy(residue, "C5", C1s, 1485, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1486, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1481, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1482, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1487, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1483, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1484, forceField, bondList);
                O4 = BondedUtils.buildHeavy(residue, "O4", C4, 1489, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1491, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H3", N3, 1.00e0, C2, 120.0e0, O2, 0.0e0, 0, 1488, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 120.0e0, C1s, 0.0e0, 0, 1490, forceField, bondList);
                break;
            case H2U:
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1350, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1351, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1356, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1352, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1353, forceField, bondList);
                O4 = BondedUtils.buildHeavy(residue, "O4", C4, 1358, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1354, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1355, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1357, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H51", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, 1, 1359, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H52", C5, 1.08e0, C4, 109.5e0, C6, 109.5e0, -1, 1360, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H61", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, 1, 1361, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H62", C6, 1.08e0, C5, 109.5e0, N1, 109.5e0, -1, 1362, forceField, bondList);
                break;
            case M5MU:
                Atom C5M;
                Atom H5M1;
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1575, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1576, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1581, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1577, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1578, forceField, bondList);
                O4 = BondedUtils.buildHeavy(residue, "O4", C4, 1583, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1579, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1580, forceField, bondList);
                C5M = BondedUtils.buildHeavy(residue, "C5M", C5, 1585, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.5e0, N1, 180.0e0, 0, 1582, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 118.6e0, C4, 180.0e0, 0, 1584, forceField, bondList);
                H5M1 = BondedUtils.buildHydrogen(residue, "H5M1", C5M, 1.08e0, C5, 109.5e0, C6, 0.0e0, 0, 1586, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H5M2", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, 1, 1587, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H5M3", C5M, 1.08e0, C5, 109.5e0, H5M1, 109.5e0, -1, 1588, forceField, bondList);
                break;
            case DAD:
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1132, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1.37, C1s, 128.4, O4s, glyco + 180, 0, 1136, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1.30, N9, 113.8, C1s, 180.0, 0, 1135, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1.39, C8, 104.0, N9, 0.0, 0, 1134, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.40, N7, 132.4, C8, 180.0, 0, 1140, forceField, bondList);
                N6 = BondedUtils.buildHeavy(residue, "N6", C6, 1.34, C5, 123.5, N7, 0.0, 0, 1142, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1.35, C5, 117.4, N7, 180.0, 0, 1139, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.33, C6, 118.8, C5, 0.0, 0, 1138, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.32, N1, 129.2, C6, 0.0, 0, 1137, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.35, C2, 110.9, N1, 0.0, 0, 1133, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.1e0, C5, 180.0e0, 0, 1145, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H61", N6, 1.00e0, C6, 120.0e0, N7, 180.0e0, 0, 1143, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H62", N6, 1.00e0, C6, 120.0e0, N7, 0.0e0, 0, 1144, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H2", C2, 1.08e0, N3, 115.4e0, C4, 180.0e0, 0, 1141, forceField, bondList);
                break;
            case DCY:
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1191, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.37, C1s, 117.8, O4s, glyco, 0, 1192, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1.24, N1, 118.9, C1s, 0.0, 0, 1197, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.38, N1, 118.7, C1s, 180, 0, 1193, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.34, C2, 120.6, N1, 0.0, 0, 1194, forceField, bondList);
                N4 = BondedUtils.buildHeavy(residue, "N4", C4, 1.32, N3, 118.3, C2, 180.0, 0, 1198, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1.43, N3, 121.6, C2, 0.0, 0, 1195, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.36, C4, 116.9, N3, 0.0, 0, 1196, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H41", N4, 1.00e0, C4, 120.0e0, N3, 0.0e0, 0, 1199, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H42", N4, 1.00e0, C4, 120.0e0, N3, 180.0e0, 0, 1200, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H5", C5, 1.08e0, C4, 121.6e0, N3, 180.0e0, 0, 1201, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1202, forceField, bondList);
                break;
            case DGU:
                N9 = BondedUtils.buildHeavy(residue, "N9", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1161, forceField, bondList);
                C8 = BondedUtils.buildHeavy(residue, "C8", N9, 1.38, C1s, 128.4, O4s, glyco + 180, 0, 1165, forceField, bondList);
                N7 = BondedUtils.buildHeavy(residue, "N7", C8, 1.31, N9, 114.0, C1s, 180.0, 0, 1164, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", N7, 1.39, C8, 103.8, N9, 0.0, 0, 1163, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.40, N7, 130.1, C8, 180.0, 0, 1169, forceField, bondList);
                O6 = BondedUtils.buildHeavy(residue, "O6", C6, 1.23, C5, 128.8, N7, 0.0, 0, 1174, forceField, bondList);
                N1 = BondedUtils.buildHeavy(residue, "N1", C6, 1.40, C5, 111.4, N7, 180.0, 0, 1168, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.38, C6, 125.2, C5, 0.0, 0, 1167, forceField, bondList);
                N2 = BondedUtils.buildHeavy(residue, "N2", C2, 1.34, N1, 116.1, C6, 180.0, 0, 1171, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.33, N1, 123.3, O6, 0.0, 0, 1166, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.36, C2, 112.3, N1, 0.0, 0, 1162, forceField, bondList);
                BondedUtils.buildBond(C4, C5, forceField, bondList);
                BondedUtils.buildBond(C4, N9, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H8", C8, 1.08e0, N7, 123.0e0, C5, 180.0e0, 0, 1175, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H1", N1, 1.00e0, C6, 117.4e0, C5, 180.0e0, 0, 1170, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H21", N2, 1.00e0, C2, 120.0e0, N1, 0.0e0, 0, 1172, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H22", N2, 1.00e0, C2, 120.0e0, N1, 180.0e0, 0, 1173, forceField, bondList);
                break;
            case DTY:
                Atom C7;
                Atom H;
                N1 = BondedUtils.buildHeavy(residue, "N1", C1s, 1.48, O4s, 108.1, C2s, 113.7, 1, 1218, forceField, bondList);
                C2 = BondedUtils.buildHeavy(residue, "C2", N1, 1.37, C1s, 117.1, O4s, glyco, 0, 1219, forceField, bondList);
                O2 = BondedUtils.buildHeavy(residue, "O2", C2, 1.22, N1, 122.9, C1s, 0.0, 0, 1224, forceField, bondList);
                N3 = BondedUtils.buildHeavy(residue, "N3", C2, 1.38, N1, 115.4, C1s, 180.0, 0, 1220, forceField, bondList);
                C4 = BondedUtils.buildHeavy(residue, "C4", N3, 1.38, C2, 126.4, N1, 0.0, 0, 1221, forceField, bondList);
                O4 = BondedUtils.buildHeavy(residue, "O4", C4, 1.23, N3, 120.5, C2, 180.0, 0, 1226, forceField, bondList);
                C5 = BondedUtils.buildHeavy(residue, "C5", C4, 1.44, N3, 114.1, C2, 0.0, 0, 1222, forceField, bondList);
                C7 = BondedUtils.buildHeavy(residue, "C7", C5, 1.50, C4, 117.5, N3, 180.0, 0, 1227, forceField, bondList);
                C6 = BondedUtils.buildHeavy(residue, "C6", C5, 1.34, C4, 120.8, N3, 0.0, 0, 1223, forceField, bondList);
                BondedUtils.buildBond(C6, N1, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H3", N3, 1.00e0, C2, 116.8e0, N1, 180.0e0, 0, 1225, forceField, bondList);
                H = BondedUtils.buildHydrogen(residue, "H71", C7, 1.09e0, C5, 109.5e0, C4, 0.0e0, 0, 1228, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H72", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, 1, 1228, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H73", C7, 1.09e0, C5, 109.5e0, H, 109.5e0, -1, 1228, forceField, bondList);
                BondedUtils.buildHydrogen(residue, "H6", C6, 1.08e0, C5, 119.4e0, C4, 180.0e0, 0, 1229, forceField, bondList);
                break;
        }
    }
}
