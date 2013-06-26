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
package ffx.potential;

import java.util.logging.Logger;

import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Residue;

import static ffx.potential.parsers.INTFilter.intxyz;

/**
 * The Rotamer Library Class manages a library of side-chain Rotamers.
 *
 * @author Ava M. Lynn
 * @author Shibo Gao
 */
public class RotamerLibrary {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    private static final int numberOfAminoAcids = AminoAcid3.values().length;
    private static final Rotamer[][] rotamerCache = new Rotamer[numberOfAminoAcids][];
    private static LibraryName libraryName = LibraryName.PonderAndRichards;
    private static int evaluatedPermutations = 0;

    public static void setLibrary(LibraryName name) {
        libraryName = name;
        for (int i = 0; i < numberOfAminoAcids; i++) {
            rotamerCache[i] = null;
        }
    }

    public static LibraryName getLibrary() {
        return libraryName;
    }

    public static Rotamer[] getRotamers(Residue residue) {
        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
        return getRotamers(name);
    }

    /**
     * Return an array of Rotamers for the given amino acid.
     *
     * @param name The name of the amino acid.
     * @return An array of Rotamers.
     */
    public static Rotamer[] getRotamers(AminoAcid3 name) {
        Rotamer[] rotamers = null;
        switch (libraryName) {
            case PonderAndRichards:
                rotamers = getPonderAndRichardsRotamers(name);
                break;
            case Richardson:
                rotamers = getRichardsonRotamers(name);
                break;
        }
        return rotamers;
    }

    private static Rotamer[] getPonderAndRichardsRotamers(AminoAcid3 name) {
        int n = name.ordinal();
        if (rotamerCache[n] != null) {
            return rotamerCache[n];
        }
        switch (name) {
            case VAL:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 173.5, 9.0);
                rotamerCache[n][1] = new Rotamer(name, -63.4, 8.1);
                rotamerCache[n][2] = new Rotamer(name, 69.3, 9.6);
                break;
            case LEU:
                rotamerCache[n] = new Rotamer[4];
                rotamerCache[n][0] = new Rotamer(name, -64.9, 8.2, 176.0, 9.9);
                rotamerCache[n][1] = new Rotamer(name, -176.4, 10.2, 63.1, 8.2);
                rotamerCache[n][2] = new Rotamer(name, -165.3, 10.0, 168.2, 34.2);
                rotamerCache[n][3] = new Rotamer(name, 44.3, 20.0, 60.4, 18.8);
                break;
            case ILE:
                rotamerCache[n] = new Rotamer[5];
                rotamerCache[n][0] = new Rotamer(name, -60.9, 7.5, 168.7, 11.6);
                rotamerCache[n][1] = new Rotamer(name, -59.6, 9.6, -64.1, 14.3);
                rotamerCache[n][2] = new Rotamer(name, 61.7, 5.0, 163.8, 16.4);
                rotamerCache[n][3] = new Rotamer(name, -166.6, 10.1, 166.0, 8.9);
                rotamerCache[n][4] = new Rotamer(name, -174.8, 24.9, 72.1, 10.5);
                break;
            case SER:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 64.7, 16.1);
                rotamerCache[n][1] = new Rotamer(name, -69.7, 14.6);
                rotamerCache[n][2] = new Rotamer(name, -176.1, 20.2);
                break;
            case THR:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 62.7, 8.5);
                rotamerCache[n][1] = new Rotamer(name, -59.7, 9.4);
                rotamerCache[n][2] = new Rotamer(name, -169.5, 6.6);
                break;
            case CYS:
            case CYD:
                rotamerCache[n] = null;
//                rotamerCache[n] = new Rotamer[3];
//                rotamerCache[n][0] = new Rotamer(name, -65.2, 10.1);
//                rotamerCache[n][1] = new Rotamer(name, -179.6, 9.5);
//                rotamerCache[n][2] = new Rotamer(name, 63.5, 9.6);
                break;
            case PRO:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 24.0, 8.0);
                rotamerCache[n][1] = new Rotamer(name, 0.0, 8.0);
                rotamerCache[n][2] = new Rotamer(name, -24.0, 8.0);
                break;
            case PHE:
                rotamerCache[n] = new Rotamer[4];
                rotamerCache[n][0] = new Rotamer(name, -66.3, 10.2, 94.3, 19.5);
                rotamerCache[n][1] = new Rotamer(name, -179.2, 9.3, 7.8, 8.9);
                rotamerCache[n][2] = new Rotamer(name, 66.0, 12.0, 90.7, 9.4);
                rotamerCache[n][3] = new Rotamer(name, -71.9, 16.3, -0.4, 26.1);
                break;
            case TYR:
            case TYD:
                rotamerCache[n] = new Rotamer[4];
                rotamerCache[n][0] = new Rotamer(name, -66.5, 11.4, 96.6, 21.8);
                rotamerCache[n][1] = new Rotamer(name, -179.7, 12.6, 71.9, 13.4);
                rotamerCache[n][2] = new Rotamer(name, 63.3, 9.4, 89.1, 13.0);
                rotamerCache[n][3] = new Rotamer(name, -67.2, 13.2, -1.0, 20.1);
                break;
            case TRP:
                rotamerCache[n] = new Rotamer[6];
                rotamerCache[n][0] = new Rotamer(name, -70.4, 7.0, 100.5, 18.2);
                rotamerCache[n][1] = new Rotamer(name, 64.8, 13.0, -88.9, 5.3);
                rotamerCache[n][2] = new Rotamer(name, -177.3, 7.9, -95.1, 7.6);
                rotamerCache[n][3] = new Rotamer(name, -179.5, 3.4, 87.5, 3.8);
                rotamerCache[n][4] = new Rotamer(name, -73.3, 6.5, -87.7, 8.1);
                rotamerCache[n][5] = new Rotamer(name, 62.2, 10.0, 112.5, 15.0);
                break;
            case HIS:
            case HIE:
            case HID:
                rotamerCache[n] = new Rotamer[6];
                rotamerCache[n][0] = new Rotamer(name, -62.8, 10.0, -74.3, 17.2);
                rotamerCache[n][1] = new Rotamer(name, -175.2, 15.4, -88.7, 43.5);
                rotamerCache[n][2] = new Rotamer(name, -69.8, 5.9, 96.1, 32.2);
                rotamerCache[n][3] = new Rotamer(name, 67.9, 17.4, -80.5, 40.7);
                rotamerCache[n][4] = new Rotamer(name, -177.3, 6.3, 100.5, 14.0);
                rotamerCache[n][5] = new Rotamer(name, 48.8, 10.0, 89.5, 30.0);
                break;
            case ASH:
            case ASP:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, -68.3, 9.2, -25.7, 31.1);
                rotamerCache[n][1] = new Rotamer(name, -169.1, 9.5, 3.9, 38.9);
                rotamerCache[n][2] = new Rotamer(name, 63.7, 9.9, 2.4, 29.4);
                break;
            case ASN:
                rotamerCache[n] = new Rotamer[6];
                rotamerCache[n][0] = new Rotamer(name, -68.3, 12.3, -36.8, 25.2);
                rotamerCache[n][1] = new Rotamer(name, -177.1, 8.8, 1.3, 34.1);
                rotamerCache[n][2] = new Rotamer(name, -67.2, 10.8, 128.8, 24.2);
                rotamerCache[n][3] = new Rotamer(name, 63.9, 3.7, -6.8, 13.5);
                rotamerCache[n][4] = new Rotamer(name, -174.9, 17.9, -156.8, 58.9);
                rotamerCache[n][5] = new Rotamer(name, 63.6, 6.6, 53.8, 17.1);
                break;
            case GLU:
            case GLH:
                rotamerCache[n] = new Rotamer[7];
                rotamerCache[n][0] = new Rotamer(name, -69.6, 19.2, -177.2, 21.7, -11.4, 44.8);
                rotamerCache[n][1] = new Rotamer(name, -176.2, 14.9, 175.4, 10.6, -6.7, 39.0);
                rotamerCache[n][2] = new Rotamer(name, -64.6, 13.5, -69.1, 17.3, -33.4, 27.4);
                rotamerCache[n][3] = new Rotamer(name, -55.6, 10.6, 77.0, 6.8, 25.3, 32.6);
                rotamerCache[n][4] = new Rotamer(name, 69.8, 10.6, -179.0, 23.7, 6.6, 64.2);
                rotamerCache[n][5] = new Rotamer(name, -173.6, 14.6, 70.6, 8.7, 14.0, 37.1);
                rotamerCache[n][6] = new Rotamer(name, 63.0, 4.3, -80.4, 13.9, 16.3, 20.8);
                break;
            case GLN:
                rotamerCache[n] = new Rotamer[10];
                rotamerCache[n][0] = new Rotamer(name, -66.7, 14.1, -178.5, 14.9, -24.0, 38.0);
                rotamerCache[n][1] = new Rotamer(name, -66.7, 14.1, -178.5, 14.9, 156.0, 38.0);
                rotamerCache[n][2] = new Rotamer(name, -174.6, 11.5, -177.7, 17.2, -24.0, 38.0);
                rotamerCache[n][3] = new Rotamer(name, -174.6, 11.5, -177.7, 17.2, 156.0, 38.0);
                rotamerCache[n][4] = new Rotamer(name, -58.7, 11.2, -63.8, 16.1, -46.3, 27.7);
                rotamerCache[n][5] = new Rotamer(name, -51.3, 7.3, -90.4, 22.8, 165.0, 38.2);
                rotamerCache[n][6] = new Rotamer(name, -179.4, 21.5, 67.3, 7.9, 26.8, 38.4);
                rotamerCache[n][7] = new Rotamer(name, 167.5, 14.8, 70.9, 3.7, 174.2, 7.1);
                rotamerCache[n][8] = new Rotamer(name, 70.8, 13.0, -165.6, 9.5, -24.0, 38.0);
                rotamerCache[n][9] = new Rotamer(name, 70.8, 13.0, -165.6, 9.5, 156.0, 38.0);
                break;
            case MET:
                rotamerCache[n] = new Rotamer[13];
                rotamerCache[n][0] = new Rotamer(name, -64.5, 12.7, -68.5, 6.0, -75.6, 14.1);
                rotamerCache[n][1] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, 65.0, 20.0);
                rotamerCache[n][2] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, 180.0, 20.0);
                rotamerCache[n][3] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, -65.0, 20.0);
                rotamerCache[n][4] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, 65.0, 20.0);
                rotamerCache[n][5] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, 65.0, 20.0);
                rotamerCache[n][6] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, -65.0, 20.0);
                rotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 65.0, 20.0);
                rotamerCache[n][8] = new Rotamer(name, -170.0, 24.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][9] = new Rotamer(name, -170.0, 24.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][12] = new Rotamer(name, 61.0, 21.0, 65.0, 20.0, 180.0, 20.0);
                break;
            case LYS:
            case LYD:
                rotamerCache[n] = new Rotamer[12];
                rotamerCache[n][0] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][1] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 65.0, 20.0);
                rotamerCache[n][2] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][3] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][4] = new Rotamer(name, -170.0, 24.0, -65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][5] = new Rotamer(name, -70.0, 21.0, 65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][6] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][8] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, -65.0, 20.0);
                rotamerCache[n][9] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0, -65.0, 20.0);
                break;
            case ARG:
                rotamerCache[n] = new Rotamer[14];
                rotamerCache[n][0] = new Rotamer(name, 61.0, 25.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                rotamerCache[n][1] = new Rotamer(name, 61.0, 25.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][2] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                rotamerCache[n][3] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 90.0, 20.0);
                rotamerCache[n][4] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][5] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, -90.0, 20.0);
                rotamerCache[n][6] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                rotamerCache[n][8] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][9] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                rotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, -90.0, 20.0);
                rotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][12] = new Rotamer(name, -170.0, 21.0, 65.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                rotamerCache[n][13] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                break;
            default:
                // Handles GLY, ALA, CYX, ...
                break;
        }
        return rotamerCache[n];
    }

    private static Rotamer[] getRichardsonRotamers(AminoAcid3 name) {
        int n = name.ordinal();
        if (rotamerCache[n] != null) {
            return rotamerCache[n];
        }
        switch (name) {
            case VAL:
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 64, 0);
                rotamerCache[n][1] = new Rotamer(name, 175, 0);
                rotamerCache[n][2] = new Rotamer(name, -60, 0);
                break;
            case LEU:
                rotamerCache[n] = new Rotamer[5];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 80, 0);
                rotamerCache[n][1] = new Rotamer(name, -177, 0, 65, 0);
                rotamerCache[n][2] = new Rotamer(name, -172, 0, 145, 0);
                rotamerCache[n][3] = new Rotamer(name, -85, 0, 65, 0);
                rotamerCache[n][4] = new Rotamer(name, -65, 0, 175, 0);
                break;
            case ILE:
                rotamerCache[n] = new Rotamer[7];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 100, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 170, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 66, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 165, 0);
                rotamerCache[n][4] = new Rotamer(name, -65, 0, 100, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, 170, 0);
                rotamerCache[n][6] = new Rotamer(name, -57, 0, -60, 0);
                break;
            case SER:
                rotamerCache[n] = new Rotamer[18];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 0, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 60, 0);
                rotamerCache[n][2] = new Rotamer(name, 62, 0, 120, 0);
                rotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0);
                rotamerCache[n][4] = new Rotamer(name, 62, 0, -60, 0);
                rotamerCache[n][5] = new Rotamer(name, 62, 0, -120, 0);
                rotamerCache[n][6] = new Rotamer(name, -177, 0, 0, 0);
                rotamerCache[n][7] = new Rotamer(name, -177, 0, 60, 0);
                rotamerCache[n][8] = new Rotamer(name, -177, 0, 120, 0);
                rotamerCache[n][9] = new Rotamer(name, -177, 0, 180, 0);
                rotamerCache[n][10] = new Rotamer(name, -177, 0, -60, 0);
                rotamerCache[n][11] = new Rotamer(name, -177, 0, -120, 0);
                rotamerCache[n][12] = new Rotamer(name, -65, 0, 0, 0);
                rotamerCache[n][13] = new Rotamer(name, -65, 0, 60, 0);
                rotamerCache[n][14] = new Rotamer(name, -65, 0, 120, 0);
                rotamerCache[n][15] = new Rotamer(name, -65, 0, 180, 0);
                rotamerCache[n][16] = new Rotamer(name, -65, 0, -60, 0);
                rotamerCache[n][17] = new Rotamer(name, -65, 0, -120, 0);
                break;
            case THR:
                rotamerCache[n] = new Rotamer[18];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 0, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 60, 0);
                rotamerCache[n][2] = new Rotamer(name, 62, 0, 120, 0);
                rotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0);
                rotamerCache[n][4] = new Rotamer(name, 62, 0, -60, 0);
                rotamerCache[n][5] = new Rotamer(name, 62, 0, -120, 0);
                rotamerCache[n][6] = new Rotamer(name, -175, 0, 0, 0);
                rotamerCache[n][7] = new Rotamer(name, -175, 0, 60, 0);
                rotamerCache[n][8] = new Rotamer(name, -175, 0, 120, 0);
                rotamerCache[n][9] = new Rotamer(name, -175, 0, 180, 0);
                rotamerCache[n][10] = new Rotamer(name, -175, 0, -60, 0);
                rotamerCache[n][11] = new Rotamer(name, -175, 0, -120, 0);
                rotamerCache[n][12] = new Rotamer(name, -65, 0, 0, 0);
                rotamerCache[n][13] = new Rotamer(name, -65, 0, 60, 0);
                rotamerCache[n][14] = new Rotamer(name, -65, 0, 120, 0);
                rotamerCache[n][15] = new Rotamer(name, -65, 0, 180, 0);
                rotamerCache[n][16] = new Rotamer(name, -65, 0, -60, 0);
                rotamerCache[n][17] = new Rotamer(name, -65, 0, -120, 0);
                break;
            case CYS:
            case CYD:
                rotamerCache[n] = null;
//                rotamerCache[n] = new Rotamer[3];
//                rotamerCache[n][0] = new Rotamer(name, 62, 0);
//                rotamerCache[n][1] = new Rotamer(name, -177, 0);
//                rotamerCache[n][2] = new Rotamer(name, -65, 0);
                break;
            case PHE:
                rotamerCache[n] = new Rotamer[4];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0);
                rotamerCache[n][1] = new Rotamer(name, -177, 0, 80, 0);
                rotamerCache[n][2] = new Rotamer(name, -65, 0, -85, 0);
                rotamerCache[n][3] = new Rotamer(name, -65, 0, -30, 0);
                break;
            case TYR:
                rotamerCache[n] = new Rotamer[8];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0, 0, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 90, 0, 180, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 80, 0, 0, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 80, 0, 180, 0);
                rotamerCache[n][4] = new Rotamer(name, -65, 0, -85, 0, 0, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, -85, 0, 180, 0);
                rotamerCache[n][6] = new Rotamer(name, -65, 0, -30, 0, 0, 0);
                rotamerCache[n][7] = new Rotamer(name, -65, 0, -30, 0, 180, 0);
            case TYD:
                rotamerCache[n] = new Rotamer[4];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0);
                rotamerCache[n][1] = new Rotamer(name, -177, 0, 80, 0);
                rotamerCache[n][2] = new Rotamer(name, -65, 0, -85, 0);
                rotamerCache[n][3] = new Rotamer(name, -65, 0, -30, 0);
                break;
            case TRP:
                rotamerCache[n] = new Rotamer[7];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, -90, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 90, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, -105, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 90, 0);
                rotamerCache[n][4] = new Rotamer(name, -65, 0, -90, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, -5, 0);
                rotamerCache[n][6] = new Rotamer(name, -65, 0, 95, 0);
                break;
            case HIS:
            case HIE:
            case HID:
                rotamerCache[n] = new Rotamer[8];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, -75, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 80, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, -165, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, -80, 0);
                rotamerCache[n][4] = new Rotamer(name, -177, 0, 60, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, -70, 0);
                rotamerCache[n][6] = new Rotamer(name, -65, 0, 165, 0);
                rotamerCache[n][7] = new Rotamer(name, -65, 0, 80, 0);
                break;
            case ASH:
            case ASP:
                rotamerCache[n] = new Rotamer[5];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 10, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 30, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 0, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0);
                rotamerCache[n][4] = new Rotamer(name, -70, 0, -15, 0);
                break;
            case ASN:
                rotamerCache[n] = new Rotamer[7];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, -10, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 30, 0);
                rotamerCache[n][2] = new Rotamer(name, -174, 0, -20, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 30, 0);
                rotamerCache[n][4] = new Rotamer(name, -65, 0, -20, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, -75, 0);
                rotamerCache[n][6] = new Rotamer(name, -65, 0, 120, 0);
                break;
            case GLU:
            case GLH:
                rotamerCache[n] = new Rotamer[8];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, -20, 0);
                rotamerCache[n][1] = new Rotamer(name, 70, 0, -80, 0, 0, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, 10, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 180, 0, 0, 0);
                rotamerCache[n][4] = new Rotamer(name, -177, 0, -80, 0, -25, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, 85, 0, 0, 0);
                rotamerCache[n][6] = new Rotamer(name, -67, 0, -180, 0, -10, 0);
                rotamerCache[n][7] = new Rotamer(name, -65, 0, -65, 0, -40, 0);
                break;
            case GLN:
                rotamerCache[n] = new Rotamer[9];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 20, 0);
                rotamerCache[n][1] = new Rotamer(name, 70, 0, -75, 0, 0, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, -100, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0, 60, 0);
                rotamerCache[n][4] = new Rotamer(name, -177, 0, 180, 0, 0, 0);
                rotamerCache[n][5] = new Rotamer(name, -65, 0, 85, 0, 0, 0);
                rotamerCache[n][6] = new Rotamer(name, -67, 0, 180, 0, -25, 0);
                rotamerCache[n][7] = new Rotamer(name, -65, 0, -65, 0, -40, 0);
                rotamerCache[n][8] = new Rotamer(name, -65, 0, -65, 0, 100, 0);
                break;
            case MET:
                rotamerCache[n] = new Rotamer[13];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 75, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, -75, 0);
                rotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, 75, 0);
                rotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0, 180, 0);
                rotamerCache[n][4] = new Rotamer(name, -177, 0, 180, 0, 75, 0);
                rotamerCache[n][5] = new Rotamer(name, -177, 0, 180, 0, 180, 0);
                rotamerCache[n][6] = new Rotamer(name, -177, 0, 180, 0, -75, 0);
                rotamerCache[n][7] = new Rotamer(name, -67, 0, 180, 0, 75, 0);
                rotamerCache[n][8] = new Rotamer(name, -67, 0, 180, 0, 180, 0);
                rotamerCache[n][9] = new Rotamer(name, -67, 0, 180, 0, -75, 0);
                rotamerCache[n][10] = new Rotamer(name, -65, 0, -65, 0, 103, 0);
                rotamerCache[n][11] = new Rotamer(name, -65, 0, -65, 0, 180, 0);
                rotamerCache[n][12] = new Rotamer(name, -65, 0, -65, 0, -70, 0);
                break;
            case LYS:
            case LYD:
                rotamerCache[n] = new Rotamer[27];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 68, 0, 180, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 65.0, 0);
                rotamerCache[n][2] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0, 180, 0, -65, 0);
                rotamerCache[n][4] = new Rotamer(name, 62, 0, 180, 0, -68, 0, 180, 0);
                rotamerCache[n][5] = new Rotamer(name, -177, 0, 68, 0, 180, 0, 65, 0);
                rotamerCache[n][6] = new Rotamer(name, -177, 0, 68, 0, 180, 0, 180, 0);
                rotamerCache[n][7] = new Rotamer(name, -177, 0, 68, 0, 180, 0, -65, 0);
                rotamerCache[n][8] = new Rotamer(name, -177, 0, 180, 0, 68, 0, 65, 0);
                rotamerCache[n][9] = new Rotamer(name, -177, 0, 180, 0, 68, 0, 180, 0);
                rotamerCache[n][10] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 65, 0);
                rotamerCache[n][11] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][12] = new Rotamer(name, -177, 0, 180, 0, 180, 0, -65, 0);
                rotamerCache[n][13] = new Rotamer(name, -177, 0, 180, 0, -68, 0, 180, 0);
                rotamerCache[n][14] = new Rotamer(name, -177, 0, 180, 0, -68, 0, -65, 0);
                rotamerCache[n][15] = new Rotamer(name, -90, 0, 68, 0, 180, 0, 180);
                rotamerCache[n][16] = new Rotamer(name, -67, 0, 180, 0, 68, 0, -65, 0);
                rotamerCache[n][17] = new Rotamer(name, -67, 0, 180, 0, 68, 0, 180, 0);
                rotamerCache[n][18] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 65, 0);
                rotamerCache[n][19] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][20] = new Rotamer(name, -67, 0, 180, 0, 180, 0, -65, 0);
                rotamerCache[n][21] = new Rotamer(name, -67, 0, 180, 0, -68, 0, 180, 0);
                rotamerCache[n][22] = new Rotamer(name, -67, 0, 180, 0, -68, 0, -65, 0);
                rotamerCache[n][23] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 65, 0);
                rotamerCache[n][24] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 180, 0);
                rotamerCache[n][25] = new Rotamer(name, -62, 0, -68, 0, 180, 0, -65, 0);
                rotamerCache[n][26] = new Rotamer(name, -62, 0, -68, 0, -68, 0, 180, 0);
                break;
            case ARG:
                rotamerCache[n] = new Rotamer[34];
                rotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 65, 0, 85, 0);
                rotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, 65, 0, -175, 0);
                rotamerCache[n][2] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 85, 0);
                rotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][4] = new Rotamer(name, 62, 0, 180, 0, 180, 0, -85, 0);
                rotamerCache[n][5] = new Rotamer(name, 62, 0, 180, 0, -65, 0, 175, 0);
                rotamerCache[n][6] = new Rotamer(name, 62, 0, 180, 0, -65, 0, -85, 0);
                rotamerCache[n][7] = new Rotamer(name, -177, 0, 65, 0, 65, 0, 85, 0);
                rotamerCache[n][8] = new Rotamer(name, -177, 0, 65, 0, 65, 0, -175, 0);
                rotamerCache[n][9] = new Rotamer(name, -177, 0, 65, 0, 180, 0, 85, 0);
                rotamerCache[n][10] = new Rotamer(name, -177, 0, 65, 0, 180, 0, 180, 0);
                rotamerCache[n][11] = new Rotamer(name, -177, 0, 180, 0, 65, 0, 85, 0);
                rotamerCache[n][12] = new Rotamer(name, -177, 0, 180, 0, 65, 0, -175, 0);
                rotamerCache[n][13] = new Rotamer(name, -177, 0, 180, 0, 65, 0, -105, 0);
                rotamerCache[n][14] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 85, 0);
                rotamerCache[n][15] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][16] = new Rotamer(name, -177, 0, 180, 0, 180, 0, -85, 0);
                rotamerCache[n][17] = new Rotamer(name, -177, 0, 180, 0, -65, 0, 105, 0);
                rotamerCache[n][18] = new Rotamer(name, -177, 0, 180, 0, -65, 0, 175, 0);
                rotamerCache[n][19] = new Rotamer(name, -177, 0, 180, 0, -65, 0, -85, 0);
                rotamerCache[n][20] = new Rotamer(name, -67, 0, 180, 0, 65, 0, 85, 0);
                rotamerCache[n][21] = new Rotamer(name, -67, 0, 180, 0, 65, 0, -175, 0);
                rotamerCache[n][22] = new Rotamer(name, -67, 0, 180, 0, 65, 0, -105, 0);
                rotamerCache[n][23] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 85, 0);
                rotamerCache[n][24] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 180, 0);
                rotamerCache[n][25] = new Rotamer(name, -67, 0, 180, 0, 180, 0, -85, 0);
                rotamerCache[n][26] = new Rotamer(name, -67, 0, 180, 0, -65, 0, 105, 0);
                rotamerCache[n][27] = new Rotamer(name, -67, 0, 180, 0, -65, 0, 175, 0);
                rotamerCache[n][28] = new Rotamer(name, -67, 0, -167, 0, -65, 0, -85, 0);
                rotamerCache[n][29] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 85, 0);
                rotamerCache[n][30] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 180, 0);
                rotamerCache[n][31] = new Rotamer(name, -62, 0, -68, 0, 180, 0, -85, 0);
                rotamerCache[n][32] = new Rotamer(name, -62, 0, -68, 0, -65, 0, 175, 0);
                rotamerCache[n][33] = new Rotamer(name, -62, 0, -68, 0, -65, 0, -85, 0);
                break;
            default:
                // Handles GLY, ALA, CYX, PRO, ...
                break;
        }
        return rotamerCache[n];
    }

    public static void applyRotamer(Residue residue, Rotamer rotamer) {
        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
        switch (name) {
            case VAL: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG11 = (Atom) residue.getAtomNode("HG11");
                Atom HG12 = (Atom) residue.getAtomNode("HG12");
                Atom HG13 = (Atom) residue.getAtomNode("HG13");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Bond CG_CB = CB.getBond(CG1);
                Bond HB_CB = CB.getBond(HB);
                Bond HG_CG = HG11.getBond(CG1);
                double dCG_CB = CG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                Angle CG_CB_CA = CG1.getAngle(CB, CA);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG_CG_CB = HG11.getAngle(CG1, CB);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                intxyz(CG1, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG_CB, CA, dCG_CB_CA, CG1, 109.5, -1);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, CG1, 109.4, 1);
                intxyz(HG11, CG1, dHG_CG, CB, dHG_CG_CB, CA, 180.0, 0);
                intxyz(HG12, CG1, dHG_CG, CB, dHG_CG_CB, HG11, 109.4, 1);
                intxyz(HG13, CG1, dHG_CG, CB, dHG_CG_CB, HG11, 109.4, -1);
                intxyz(HG21, CG2, dHG_CG, CB, dHG_CG_CB, CA, 180.0, 0);
                intxyz(HG22, CG2, dHG_CG, CB, dHG_CG_CB, HG21, 109.4, 1);
                intxyz(HG23, CG2, dHG_CG, CB, dHG_CG_CB, HG21, 109.4, -1);
                break;
            }
            case LEU: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                Atom HD11 = (Atom) residue.getAtomNode("HD11");
                Atom HD12 = (Atom) residue.getAtomNode("HD12");
                Atom HD13 = (Atom) residue.getAtomNode("HD13");
                Atom HD21 = (Atom) residue.getAtomNode("HD21");
                Atom HD22 = (Atom) residue.getAtomNode("HD22");
                Atom HD23 = (Atom) residue.getAtomNode("HD23");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG.getBond(CG);
                Bond HD_CD = HD11.getBond(CD1);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD1.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG.getAngle(CG, CB);
                Angle HD_CD_CG = HD11.getAngle(CD1, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CD1, 109.5, -1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG, CG, dHG_CG, CB, dHG_CG_CB, CD1, 109.4, 1);
                intxyz(HD11, CD1, dHD_CD, CG, dHD_CD_CG, CB, 180.0, 0);
                intxyz(HD12, CD1, dHD_CD, CG, dHD_CD_CG, HD11, 109.4, 1);
                intxyz(HD13, CD1, dHD_CD, CG, dHD_CD_CG, HD11, 109.4, -1);
                intxyz(HD21, CD2, dHD_CD, CG, dHD_CD_CG, CB, 180.0, 0);
                intxyz(HD22, CD2, dHD_CD, CG, dHD_CD_CG, HD21, 109.4, 1);
                intxyz(HD23, CD2, dHD_CD, CG, dHD_CD_CG, HD21, 109.4, -1);
                break;
            }
            case ILE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG12 = (Atom) residue.getAtomNode("HG12");
                Atom HG13 = (Atom) residue.getAtomNode("HG13");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Atom HD11 = (Atom) residue.getAtomNode("HD11");
                Atom HD12 = (Atom) residue.getAtomNode("HD12");
                Atom HD13 = (Atom) residue.getAtomNode("HD13");
                Bond CG1_CB = CG1.getBond(CB);
                Bond CG2_CB = CG2.getBond(CB);
                Bond CD1_CG1 = CD1.getBond(CG1);
                Bond HB_CB = HB.getBond(CB);
                Bond HG1_CG = HG12.getBond(CG1);
                Bond HG2_CG = HG22.getBond(CG2);
                Bond HD_CD = HD12.getBond(CD1);
                double dCG1_CB = CG1_CB.bondType.distance;
                double dCG2_CB = CG2_CB.bondType.distance;
                double dCD1_CG1 = CD1_CG1.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG1_CG = HG1_CG.bondType.distance;
                double dHG2_CG = HG2_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG1_CB_CA = CG1.getAngle(CB, CA);
                Angle CG2_CB_CA = CG2.getAngle(CB, CA);
                Angle CD1_CG1_CB = CD1.getAngle(CG1, CB);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG1_CG_CB = HG12.getAngle(CG1, CB);
                Angle HG2_CG_CB = HG21.getAngle(CG2, CB);
                Angle HD_CD1_CG1 = HD11.getAngle(CD1, CG1);
                double dCG1_CB_CA = CG1_CB_CA.angleType.angle[CG1_CB_CA.nh];
                double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
                double dCD1_CG1_CB = CD1_CG1_CB.angleType.angle[CD1_CG1_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG1_CG_CB = HG1_CG_CB.angleType.angle[HG1_CG_CB.nh];
                double dHG2_CG_CB = HG2_CG_CB.angleType.angle[HG2_CG_CB.nh];
                double dHD_CD1_CG1 = HD_CD1_CG1.angleType.angle[HD_CD1_CG1.nh];
                intxyz(CG1, CB, dCG1_CB, CA, dCG1_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG2_CB, CA, dCG2_CB_CA, CG1, 109.5, -1);
                intxyz(CD1, CG1, dCD1_CG1, CB, dCD1_CG1_CB, CA, rotamer.chi2, 0);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, CG1, 109.4, 1);
                intxyz(HG12, CG1, dHG1_CG, CB, dHG1_CG_CB, CD1, 109.4, 1);
                intxyz(HG13, CG1, dHG1_CG, CB, dHG1_CG_CB, CD1, 109.4, -1);
                intxyz(HG21, CG2, dHG2_CG, CB, dHG2_CG_CB, CG1, 180.0, 0);
                intxyz(HG22, CG2, dHG2_CG, CB, dHG2_CG_CB, HG21, 109.0, 1);
                intxyz(HG23, CG2, dHG2_CG, CB, dHG2_CG_CB, HG21, 109.0, -1);
                intxyz(HD11, CD1, dHD_CD, CG1, dHD_CD1_CG1, CB, 180.0, 0);
                intxyz(HD12, CD1, dHD_CD, CG1, dHD_CD1_CG1, HD11, 109.0, 1);
                intxyz(HD13, CD1, dHD_CD, CG1, dHD_CD1_CG1, HD11, 109.0, -1);
                break;
            }
            case SER: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom OG = (Atom) residue.getAtomNode("OG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                Bond OG_CB = OG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_OG = HG.getBond(OG);
                double dOG_CB = OG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_OG = HG_OG.bondType.distance;
                Angle OG_CB_CA = OG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_OG_CB = HG.getAngle(OG, CB);
                double dOG_CB_CA = OG_CB_CA.angleType.angle[OG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_OG_CB = HG_OG_CB.angleType.angle[HG_OG_CB.nh];
                intxyz(OG, CB, dOG_CB, CA, dOG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, OG, 106.7, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, OG, 106.7, -1);
                intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, 180.0, 0);
                if (rotamer.length == 2) {
                    intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, rotamer.chi2, 0);
                } else {
                    intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, 180.0, 0);
                }
                break;
            }
            case THR: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom OG1 = (Atom) residue.getAtomNode("OG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG1 = (Atom) residue.getAtomNode("HG1");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Bond OG1_CB = OG1.getBond(CB);
                Bond CG2_CB = CG2.getBond(CB);
                Bond HB_CB = HB.getBond(CB);
                Bond HG1_OG1 = HG1.getBond(OG1);
                Bond HG2_CG2 = HG21.getBond(CG2);
                double dOG1_CB = OG1_CB.bondType.distance;
                double dCG2_CB = CG2_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG1_OG1 = HG1_OG1.bondType.distance;
                double dHG2_CG2 = HG2_CG2.bondType.distance;
                Angle OG1_CB_CA = OG1.getAngle(CB, CA);
                Angle CG2_CB_CA = CG2.getAngle(CB, CA);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG1_OG1_CB = HG1.getAngle(OG1, CB);
                Angle HG2_CG2_CB = HG21.getAngle(CG2, CB);
                double dOG1_CB_CA = OG1_CB_CA.angleType.angle[OG1_CB_CA.nh];
                double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG1_OG1_CB = HG1_OG1_CB.angleType.angle[HG1_OG1_CB.nh];
                double dHG2_CG2_CB = HG2_CG2_CB.angleType.angle[HG2_CG2_CB.nh];
                intxyz(OG1, CB, dOG1_CB, CA, dOG1_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG2_CB, CA, dCG2_CB_CA, OG1, 107.7, 1);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, OG1, 106.7, -1);
                intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, 180.0, 0);
                if (rotamer.length == 2) {
                    intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, rotamer.chi2, 0);
                } else {
                    intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, 180, 0);
                }
                intxyz(HG21, CG2, dHG2_CG2, CB, dHG2_CG2_CB, CA, 180.0, 0);
                intxyz(HG22, CG2, dHG2_CG2, CB, dHG2_CG2_CB, HG21, 109.0, 1);
                intxyz(HG23, CG2, dHG2_CG2, CB, dHG2_CG2_CB, HG21, 109.0, -1);
                break;
            }
            case CYS:
            case CYX: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom SG = (Atom) residue.getAtomNode("SG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                Bond SG_CB = SG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_SG = HG.getBond(SG);
                double dSG_CB = SG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_SG = HG_SG.bondType.distance;
                Angle SG_CB_CA = SG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_SG_CB = HG.getAngle(SG, CB);
                double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_SG_CB = HG_SG_CB.angleType.angle[HG_SG_CB.nh];
                intxyz(SG, CB, dSG_CB, CA, dSG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, -1);
                intxyz(HG, SG, dHG_SG, CB, dHG_SG_CB, CA, 180.0, 0);
                break;
            }
            case CYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom SG = (Atom) residue.getAtomNode("SG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Bond SG_CB = SG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                double dSG_CB = SG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                Angle SG_CB_CA = SG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                intxyz(SG, CB, dSG_CB, CA, dSG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, -1);
                break;
            }
            case PHE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HZ = (Atom) residue.getAtomNode("HZ");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond CZ_CE1 = CZ.getBond(CE1);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD_CD = HD1.getBond(CD1);
                Bond HE_CE = HE1.getBond(CE1);
                Bond HZ_CZ = HZ.getBond(CZ);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dCZ_CE1 = CZ_CE1.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHZ_CZ = HZ_CZ.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD1.getAngle(CG, CB);
                Angle CE_CD_CG = CE1.getAngle(CD1, CG);
                Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD_CD1_CG = HD1.getAngle(CD1, CG);
                Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
                Angle HZ_CZ_CE1 = HZ.getAngle(CZ, CE1);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD_CD1_CG = HD_CD1_CG.angleType.angle[HD_CD1_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHZ_CZ_CE1 = HZ_CZ_CE1.angleType.angle[HZ_CZ_CE1.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CD1, 120.0, 1);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CZ, CE1, dCZ_CE1, CD1, dCZ_CE1_CD1, CG, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD1_CG, CE1, 120.0, 1);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD1_CG, CE2, 120.0, 1);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, 120.0, 1);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, 120.0, 1);
                intxyz(HZ, CZ, dHZ_CZ, CE1, dHZ_CZ_CE1, CE2, 120.0, 1);
                break;
            }
            case PRO: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, N, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, N, 109.4, -1);
                break;
            }
            case TYR: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom OH = (Atom) residue.getAtomNode("OH");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HH = (Atom) residue.getAtomNode("HH");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond CZ_CE1 = CZ.getBond(CE1);
                Bond OH_CZ = OH.getBond(CZ);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD_CD = HD1.getBond(CD1);
                Bond HE_CE = HE1.getBond(CE1);
                Bond HH_OH = HH.getBond(OH);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dCZ_CE1 = CZ_CE1.bondType.distance;
                double dOH_CZ = OH_CZ.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHH_OH = HH_OH.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD1.getAngle(CG, CB);
                Angle CE_CD_CG = CE1.getAngle(CD1, CG);
                Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
                Angle OH_CZ_CE2 = OH.getAngle(CZ, CE2);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD_CD_CG = HD1.getAngle(CD1, CG);
                Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
                Angle HH_OH_CZ = HH.getAngle(OH, CZ);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
                double dOH_CZ_CE2 = OH_CZ_CE2.angleType.angle[OH_CZ_CE2.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHH_OH_CZ = HH_OH_CZ.angleType.angle[HH_OH_CZ.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CD1, 120.0, 1);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CZ, CE1, dCZ_CE1, CD1, dCZ_CE1_CD1, CG, 0.0, 0);
                intxyz(OH, CZ, dOH_CZ, CE2, dOH_CZ_CE2, CE1, 120.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD_CG, CE1, 120.0, 1);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD_CG, CE2, 120.0, 1);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, 120.0, 1);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, 120.0, 1);
                if (rotamer.length == 3) {
                    intxyz(HH, OH, dHH_OH, CZ, dHH_OH_CZ, CE2, rotamer.chi3, 0);
                } else {
                    intxyz(HH, OH, dHH_OH, CZ, dHH_OH_CZ, CE2, 0.0, 0);
                }
                break;
            }
            case TYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom OH = (Atom) residue.getAtomNode("OH");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond CZ_CE1 = CZ.getBond(CE1);
                Bond OH_CZ = OH.getBond(CZ);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD_CD = HD1.getBond(CD1);
                Bond HE_CE = HE1.getBond(CE1);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dCZ_CE1 = CZ_CE1.bondType.distance;
                double dOH_CZ = OH_CZ.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD1.getAngle(CG, CB);
                Angle CE_CD_CG = CE1.getAngle(CD1, CG);
                Angle CZ_CE1_CD1 = CZ.getAngle(CE1, CD1);
                Angle OH_CZ_CE2 = OH.getAngle(CZ, CE2);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD_CD_CG = HD1.getAngle(CD1, CG);
                Angle HE_CE_CD = HE1.getAngle(CE1, CD1);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dCZ_CE1_CD1 = CZ_CE1_CD1.angleType.angle[CZ_CE1_CD1.nh];
                double dOH_CZ_CE2 = OH_CZ_CE2.angleType.angle[OH_CZ_CE2.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CD1, 120.0, 1);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CZ, CE1, dCZ_CE1, CD1, dCZ_CE1_CD1, CG, 0.0, 0);
                intxyz(OH, CZ, dOH_CZ, CE2, dOH_CZ_CE2, CE1, 120.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD_CG, CE1, 120.0, 1);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD_CG, CE2, 120.0, 1);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, 120.0, 1);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, 120.0, 1);
                break;
            }
            case TRP: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom NE1 = (Atom) residue.getAtomNode("NE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CE3 = (Atom) residue.getAtomNode("CE3");
                Atom CZ2 = (Atom) residue.getAtomNode("CZ2");
                Atom CZ3 = (Atom) residue.getAtomNode("CZ3");
                Atom CH2 = (Atom) residue.getAtomNode("CH2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Atom HZ3 = (Atom) residue.getAtomNode("HZ3");
                Atom HH2 = (Atom) residue.getAtomNode("HH2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD1_CG = CD1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond NE1_CD1 = NE1.getBond(CD1);
                Bond CE2_NE1 = CE2.getBond(NE1);
                Bond CE3_CD2 = CE3.getBond(CD2);
                Bond CZ2_CE2 = CZ2.getBond(CE2);
                Bond CZ3_CE3 = CZ3.getBond(CE3);
                Bond CH2_CZ2 = CH2.getBond(CZ2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_CD1 = HD1.getBond(CD1);
                Bond HE1_NE1 = HE1.getBond(NE1);
                Bond HE3_CE3 = HE3.getBond(CE3);
                Bond HZ2_CZ2 = HZ2.getBond(CZ2);
                Bond HZ3_CZ3 = HZ3.getBond(CZ3);
                Bond HH2_CH2 = HH2.getBond(CH2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD1_CG = CD1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dNE1_CD1 = NE1_CD1.bondType.distance;
                double dCE2_NE1 = CE2_NE1.bondType.distance;
                double dCE3_CD2 = CE3_CD2.bondType.distance;
                double dCZ2_CE2 = CZ2_CE2.bondType.distance;
                double dCZ3_CE3 = CZ3_CE3.bondType.distance;
                double dCH2_CZ2 = CH2_CZ2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_CD1 = HD1_CD1.bondType.distance;
                double dHE1_NE1 = HE1_NE1.bondType.distance;
                double dHE3_CE3 = HE3_CE3.bondType.distance;
                double dHZ2_CZ2 = HZ2_CZ2.bondType.distance;
                double dHZ3_CZ3 = HZ3_CZ3.bondType.distance;
                double dHH2_CH2 = HH2_CH2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD1_CG_CB = CD1.getAngle(CG, CB);
                Angle CD2_CG_CB = CD2.getAngle(CG, CB);
                Angle NE1_CD1_CG = NE1.getAngle(CD1, CG);
                Angle CE2_NE1_CD1 = CE2.getAngle(NE1, CD1);
                Angle CE3_CD2_CE2 = CE3.getAngle(CD2, CE2);
                Angle CZ2_CE2_CD2 = CZ2.getAngle(CE2, CD2);
                Angle CZ3_CE3_CD2 = CZ3.getAngle(CE3, CD2);
                Angle CH2_CZ2_CE2 = CH2.getAngle(CZ2, CE2);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD1_CD1_CG = HD1.getAngle(CD1, CG);
                Angle HE1_NE1_CD1 = HE1.getAngle(NE1, CD1);
                Angle HE3_CE3_CD2 = HE3.getAngle(CE3, CD2);
                Angle HZ2_CZ2_CE2 = HZ2.getAngle(CZ2, CE2);
                Angle HZ3_CZ3_CE3 = HZ3.getAngle(CZ3, CH2);
                Angle HH2_CH2_CZ2 = HH2.getAngle(CH2, CZ3);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD1_CG_CB = CD1_CG_CB.angleType.angle[CD1_CG_CB.nh];
                double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
                double dNE1_CD1_CG = NE1_CD1_CG.angleType.angle[NE1_CD1_CG.nh];
                double dCE2_NE1_CD1 = CE2_NE1_CD1.angleType.angle[CE2_NE1_CD1.nh];
                double dCE3_CD2_CE2 = CE3_CD2_CE2.angleType.angle[CE3_CD2_CE2.nh];
                double dCZ2_CE2_CD2 = CZ2_CE2_CD2.angleType.angle[CZ2_CE2_CD2.nh];
                double dCZ3_CE3_CD2 = CZ3_CE3_CD2.angleType.angle[CZ3_CE3_CD2.nh];
                double dCH2_CZ2_CE2 = CH2_CZ2_CE2.angleType.angle[CH2_CZ2_CE2.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD1_CD1_CG = HD1_CD1_CG.angleType.angle[HD1_CD1_CG.nh];
                double dHE1_NE1_CD1 = HE1_NE1_CD1.angleType.angle[HE1_NE1_CD1.nh];
                double dHE3_CE3_CD2 = HE3_CE3_CD2.angleType.angle[HE3_CE3_CD2.nh];
                double dHZ2_CZ2_CE2 = HZ2_CZ2_CE2.angleType.angle[HZ2_CZ2_CE2.nh];
                double dHZ3_CZ3_CE3 = HZ3_CZ3_CE3.angleType.angle[HZ3_CZ3_CE3.nh];
                double dHH2_CH2_CZ2 = HH2_CH2_CZ2.angleType.angle[HH2_CH2_CZ2.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD1_CG, CB, dCD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, CD1, 108.0, 1);
                intxyz(NE1, CD1, dNE1_CD1, CG, dNE1_CD1_CG, CD2, 0.0, 0);
                intxyz(CE2, NE1, dCE2_NE1, CD1, dCE2_NE1_CD1, CG, 0.0, 0);
                intxyz(CE3, CD2, dCE3_CD2, CE2, dCE3_CD2_CE2, NE1, 180.0, 0);
                intxyz(CZ2, CE2, dCZ2_CE2, CD2, dCZ2_CE2_CD2, CE3, 0.0, 0);
                intxyz(CZ3, CE3, dCZ3_CE3, CD2, dCZ3_CE3_CD2, CE2, 0.0, 0);
                intxyz(CH2, CZ2, dCH2_CZ2, CE2, dCH2_CZ2_CE2, CD2, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, CD1, dHD1_CD1, CG, dHD1_CD1_CG, NE1, 126.0, 1);
                intxyz(HE1, NE1, dHE1_NE1, CD1, dHE1_NE1_CD1, CE2, 126.0, 1);
                intxyz(HE3, CE3, dHE3_CE3, CD2, dHE3_CE3_CD2, CZ3, 120.0, 1);
                intxyz(HZ2, CZ2, dHZ2_CZ2, CE2, dHZ2_CZ2_CE2, CH2, 120.0, 1);
                intxyz(HZ3, CZ3, dHZ3_CZ3, CE3, dHZ3_CZ3_CE3, CH2, 120.0, 1);
                intxyz(HH2, CH2, dHH2_CH2, CZ2, dHH2_CH2_CZ2, CZ3, 120.0, 1);
                break;
            }
            case HIS: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_ND1 = HD1.getBond(ND1);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                Bond HE2_NE2 = HE2.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_ND1 = HD1_ND1.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle ND1_CG_CB = ND1.getAngle(CG, CB);
                Angle CD2_CG_CB = CD2.getAngle(CG, CB);
                Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
                Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD1_ND1_CG = HD1.getAngle(ND1, CG);
                Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
                Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
                Angle HE2_NE2_CD2 = HE2.getAngle(NE2, CD2);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
                double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
                double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
                double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD1_ND1_CG = HD1_ND1_CG.angleType.angle[HD1_ND1_CG.nh];
                double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
                double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
                double dHE2_NE2_CD2 = HE2_NE2_CD2.angleType.angle[HE2_NE2_CD2.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, ND1, 108.0, 1);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, ND1, dHD1_ND1, CG, dHD1_ND1_CG, CB, 0.0, 0);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, 126.0, 1);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, 126.0, 1);
                intxyz(HE2, NE2, dHE2_NE2, CD2, dHE2_NE2_CD2, CE1, 126.0, 1);
                break;
            }
            case HID: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_ND1 = HD1.getBond(ND1);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_ND1 = HD1_ND1.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle ND1_CG_CB = ND1.getAngle(CG, CB);
                Angle CD2_CG_CB = CD2.getAngle(CG, CB);
                Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
                Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD1_ND1_CG = HD1.getAngle(ND1, CG);
                Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
                Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
                double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
                double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
                double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD1_ND1_CG = HD1_ND1_CG.angleType.angle[HD1_ND1_CG.nh];
                double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
                double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, ND1, 108.0, 1);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, ND1, dHD1_ND1, CG, dHD1_ND1_CG, CB, 0.0, 0);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, 126.0, 1);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, 126.0, 1);
                break;
            }
            case HIE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                Bond HE2_NE2 = HE2.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle ND1_CG_CB = ND1.getAngle(CG, CB);
                Angle CD2_CG_CB = CD2.getAngle(CG, CB);
                Angle CE1_ND1_CG = CE1.getAngle(ND1, CG);
                Angle NE2_CD2_CG = NE2.getAngle(CD2, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD2_CD2_CG = HD2.getAngle(CD2, CG);
                Angle HE1_CE1_ND1 = HE1.getAngle(CE1, ND1);
                Angle HE2_NE2_CD2 = HE2.getAngle(NE2, CD2);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dND1_CG_CB = ND1_CG_CB.angleType.angle[ND1_CG_CB.nh];
                double dCD2_CG_CB = CD2_CG_CB.angleType.angle[CD2_CG_CB.nh];
                double dCE1_ND1_CG = CE1_ND1_CG.angleType.angle[CE1_ND1_CG.nh];
                double dNE2_CD2_CG = NE2_CD2_CG.angleType.angle[NE2_CD2_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD2_CD2_CG = HD2_CD2_CG.angleType.angle[HD2_CD2_CG.nh];
                double dHE1_CE1_ND1 = HE1_CE1_ND1.angleType.angle[HE1_CE1_ND1.nh];
                double dHE2_NE2_CD2 = HE2_NE2_CD2.angleType.angle[HE2_NE2_CD2.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, ND1, 108.0, 1);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, 126.0, 1);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, 126.0, 1);
                intxyz(HE2, NE2, dHE2_NE2, CD2, dHE2_NE2_CD2, CE1, 126.0, 1);
                break;
            }
            case ASP: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom OD2 = (Atom) residue.getAtomNode("OD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond OD2_CG = OD2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dOD2_CG = OD2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle OD2_CG_CB = OD2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, 0.0, 0);
                intxyz(OD2, CG, dOD2_CG, CB, dOD2_CG_CB, OD1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                break;
            }
            case ASH: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom OD2 = (Atom) residue.getAtomNode("OD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond OD2_CG = OD2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_OD2 = HD2.getBond(OD2);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dOD2_CG = OD2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_OD2 = HD2_OD2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle OD2_CG_CB = OD2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD2_OD2_CG = HD2.getAngle(OD2, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD2_OD2_CG = HD2_OD2_CG.angleType.angle[HD2_OD2_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OD2, CG, dOD2_CG, CB, dOD2_CG_CB, OD1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                intxyz(HD2, OD2, dHD2_OD2, CG, dHD2_OD2_CG, OD1, 0.0, 0);
                break;
            }
            case ASN: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom ND2 = (Atom) residue.getAtomNode("ND2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD21 = (Atom) residue.getAtomNode("HD21");
                Atom HD22 = (Atom) residue.getAtomNode("HD22");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond ND2_CG = ND2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_ND2 = HD21.getBond(ND2);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dND2_CG = ND2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_ND2 = HD2_ND2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle ND2_CG_CB = ND2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD2_ND2_CG = HD21.getAngle(ND2, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dND2_CG_CB = ND2_CG_CB.angleType.angle[ND2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD2_ND2_CG = HD2_ND2_CG.angleType.angle[HD2_ND2_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(ND2, CG, dND2_CG, CB, dND2_CG_CB, OD1, 124.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                intxyz(HD21, ND2, dHD2_ND2, CG, dHD2_ND2_CG, CB, 0.0, 0);
                intxyz(HD22, ND2, dHD2_ND2, CG, dHD2_ND2_CG, HD21, 120.0, 1);
                break;
            }
            case GLU: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom OE2 = (Atom) residue.getAtomNode("OE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond OE2_CD = OE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dOE2_CD = OE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle OE2_CD_CG = OE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, dOE2_CD, CG, dOE2_CD_CG, OE1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                break;
            }
            case GLH: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom OE2 = (Atom) residue.getAtomNode("OE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond OE2_CD = OE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE2_OE2 = HE2.getBond(OE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dOE2_CD = OE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE2_OE2 = HE2_OE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle OE2_CD_CG = OE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE2_OE2_CD = HE2.getAngle(OE2, CD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE2_OE2_CD = HE2_OE2_CD.angleType.angle[HE2_OE2_CD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, dOE2_CD, CG, dOE2_CD_CG, OE1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                intxyz(HE2, OE2, dHE2_OE2, CD, dHE2_OE2_CD, OE1, 0.0, 0);
                break;
            }
            case GLN: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE21 = (Atom) residue.getAtomNode("HE21");
                Atom HE22 = (Atom) residue.getAtomNode("HE22");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond NE2_CD = NE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE2_NE2 = HE21.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dNE2_CD = NE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle NE2_CD_CG = NE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE2_NE2_CD = HE21.getAngle(NE2, CD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dNE2_CD_CG = NE2_CD_CG.angleType.angle[NE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE2_NE2_CD = HE2_NE2_CD.angleType.angle[HE2_NE2_CD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NE2, CD, dNE2_CD, CG, dNE2_CD_CG, OE1, 124.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                intxyz(HE21, NE2, dHE2_NE2, CD, dHE2_NE2_CD, CG, 0.0, 0);
                intxyz(HE22, NE2, dHE2_NE2, CD, dHE2_NE2_CD, HE21, 120.0, 1);
                break;
            }
            case MET: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom SD = (Atom) residue.getAtomNode("SD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Bond CG_CB = CG.getBond(CB);
                Bond SD_CG = SD.getBond(CG);
                Bond CE_SD = CE.getBond(SD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE_CE = HE1.getBond(CE);
                double dCG_CB = CG_CB.bondType.distance;
                double dSD_CG = SD_CG.bondType.distance;
                double dCE_SD = CE_SD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle SD_CG_CB = SD.getAngle(CG, CB);
                Angle CE_SD_CG = CE.getAngle(SD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE_CE_SD = HE1.getAngle(CE, SD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dSD_CG_CB = SD_CG_CB.angleType.angle[SD_CG_CB.nh];
                double dCE_SD_CG = CE_SD_CG.angleType.angle[CE_SD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE_CE_SD = HE_CE_SD.angleType.angle[HE_CE_SD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(SD, CG, dSD_CG, CB, dSD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, SD, dCE_SD, CG, dCE_SD_CG, CB, rotamer.chi3, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, SD, 112.0, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, SD, 112.0, -1);
                intxyz(HE1, CE, dHE_CE, SD, dHE_CE_SD, CG, 180.0, 0);
                intxyz(HE2, CE, dHE_CE, SD, dHE_CE_SD, HE1, 109.4, 1);
                intxyz(HE3, CE, dHE_CE, SD, dHE_CE_SD, HE1, 109.4, -1);
                break;
            }
            case LYS: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ1 = (Atom) residue.getAtomNode("HZ1");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Atom HZ3 = (Atom) residue.getAtomNode("HZ3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond CE_CD = CE.getBond(CD);
                Bond NZ_CE = NZ.getBond(CE);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_CE = HE2.getBond(CE);
                Bond HZ_NZ = HZ1.getBond(NZ);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dNZ_CE = NZ_CE.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHZ_NZ = HZ_NZ.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle CE_CD_CG = CE.getAngle(CD, CG);
                Angle NZ_CE_CD = NZ.getAngle(CE, CD);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_CE_CD = HE2.getAngle(CE, CD);
                Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, CD, dCE_CD, CG, dCE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, dNZ_CE, CD, dNZ_CE_CD, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, -1);
                intxyz(HE2, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, 1);
                intxyz(HE3, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, -1);
                intxyz(HZ1, NZ, dHZ_NZ, CE, dHZ_NZ_CE, CD, 180.0, 0);
                intxyz(HZ2, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, 1);
                intxyz(HZ3, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, -1);
                break;
            }
            case LYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ1 = (Atom) residue.getAtomNode("HZ1");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond CE_CD = CE.getBond(CD);
                Bond NZ_CE = NZ.getBond(CE);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_CE = HE2.getBond(CE);
                Bond HZ_NZ = HZ1.getBond(NZ);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dNZ_CE = NZ_CE.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHZ_NZ = HZ_NZ.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle CE_CD_CG = CE.getAngle(CD, CG);
                Angle NZ_CE_CD = NZ.getAngle(CE, CD);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_CE_CD = HE2.getAngle(CE, CD);
                Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, CD, dCE_CD, CG, dCE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, dNZ_CE, CD, dNZ_CE_CD, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, -1);
                intxyz(HE2, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, 1);
                intxyz(HE3, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, -1);
                intxyz(HZ1, NZ, dHZ_NZ, CE, dHZ_NZ_CE, CD, 180.0, 0);
                intxyz(HZ2, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, 1);
                break;
            }
            case ARG: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom NE = (Atom) residue.getAtomNode("NE");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom NH1 = (Atom) residue.getAtomNode("NH1");
                Atom NH2 = (Atom) residue.getAtomNode("NH2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE = (Atom) residue.getAtomNode("HE");
                Atom HH11 = (Atom) residue.getAtomNode("HH11");
                Atom HH12 = (Atom) residue.getAtomNode("HH12");
                Atom HH21 = (Atom) residue.getAtomNode("HH21");
                Atom HH22 = (Atom) residue.getAtomNode("HH22");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond NE_CD = NE.getBond(CD);
                Bond CZ_NE = CZ.getBond(NE);
                Bond NH_CZ = NH1.getBond(CZ);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_NE = HE.getBond(NE);
                Bond HH_NH = HH11.getBond(NH1);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dNE_CD = NE_CD.bondType.distance;
                double dCZ_NE = CZ_NE.bondType.distance;
                double dNH_CZ = NH_CZ.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_NE = HE_NE.bondType.distance;
                double dHH_NH = HH_NH.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle NE_CD_CG = NE.getAngle(CD, CG);
                Angle CZ_NE_CD = CZ.getAngle(NE, CD);
                Angle NH_CZ_NE = NH1.getAngle(CZ, NE);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_NE_CD = HE.getAngle(NE, CD);
                Angle HH_NH_CZ = HH11.getAngle(NH1, CZ);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dNE_CD_CG = NE_CD_CG.angleType.angle[NE_CD_CG.nh];
                double dCZ_NE_CD = CZ_NE_CD.angleType.angle[CZ_NE_CD.nh];
                double dNH_CZ_NE = NH_CZ_NE.angleType.angle[NH_CZ_NE.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_NE_CD = HE_NE_CD.angleType.angle[HE_NE_CD.nh];
                double dHH_NH_CZ = HH_NH_CZ.angleType.angle[HH_NH_CZ.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(NE, CD, dNE_CD, CG, dNE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(CZ, NE, dCZ_NE, CD, dCZ_NE_CD, CG, rotamer.chi4, 0);
                intxyz(NH1, CZ, dNH_CZ, NE, dNH_CZ_NE, CD, 180, 0);
                intxyz(NH2, CZ, dNH_CZ, NE, dNH_CZ_NE, NH1, 120.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, NE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, NE, 109.4, -1);
                intxyz(HE, NE, dHE_NE, CD, dHE_NE_CD, CZ, 120.0, 1);
                intxyz(HH11, NH1, dHH_NH, CZ, dHH_NH_CZ, NE, 180.0, 0);
                intxyz(HH12, NH1, dHH_NH, CZ, dHH_NH_CZ, HH11, 120.0, 1);
                intxyz(HH21, NH2, dHH_NH, CZ, dHH_NH_CZ, NE, 180.0, 0);
                intxyz(HH22, NH2, dHH_NH, CZ, dHH_NH_CZ, HH21, 120.0, 1);
                break;
            }
            default:
                break;
        }
    }

    public enum LibraryName {

        PonderAndRichards, Richardson
    }
}

