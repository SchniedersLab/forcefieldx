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

import java.util.ArrayList;
import java.util.logging.Logger;

import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Residue;

import static ffx.potential.parsers.INTFilter.intxyz;

/**
 * The Rotamer Library Class manages a library of side-chain Rotamers.
 *
 * @author Ava M. Lynn
 */
public class RotamerLibrary {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    private static final int numberOfAminoAcids = AminoAcid3.values().length;
    private static final Rotamer[][] rotamerCache = new Rotamer[numberOfAminoAcids][];
    private static LibraryName libraryName = LibraryName.PonderAndRichards;

    public enum LibraryName {

        PonderAndRichards, Richardson
    }

    public static void setLibrary(LibraryName name) {
        libraryName = name;
        for (int i = 0; i < numberOfAminoAcids; i++) {
            rotamerCache[i] = null;
        }
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
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, -65.2, 10.1);
                rotamerCache[n][1] = new Rotamer(name, -179.6, 9.5);
                rotamerCache[n][2] = new Rotamer(name, 63.5, 9.6);
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
                rotamerCache[n] = new Rotamer[3];
                rotamerCache[n][0] = new Rotamer(name, 62, 0);
                rotamerCache[n][1] = new Rotamer(name, -177, 0);
                rotamerCache[n][2] = new Rotamer(name, -65, 0);
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

    /**
     * A brute-force global optimization over side-chain rotamers using a
     * recursive algorithm.
     */
    public static double rotamerOptimization(MolecularAssembly molecularAssembly, ArrayList<Residue> residues,
            double lowEnergy, ArrayList<Integer> optimum) {
        Residue current = residues.remove(0);
        AminoAcid3 name = AminoAcid3.valueOf(current.getName());
        Rotamer[] rotamers = getRotamers(name);
        double currentEnergy = Double.MAX_VALUE;
        if (residues.size() > 0) {
            /**
             * As long as there are more residues, continue the recursion for
             * each rotamer of the current residue.
             */
            if (rotamers == null) {
                /**
                 * Continue to the next residue.
                 */
                currentEnergy = rotamerOptimization(molecularAssembly, residues, lowEnergy, optimum);
                // Add the '-1' flag as a placeholder since this residue has no rotamers.
                if (currentEnergy < lowEnergy) {
                    optimum.add(0, -1);
                }
            } else {
                int minRot = -1;
                for (int i = 0; i < rotamers.length; i++) {
                    applyRotamer(name, current, rotamers[i]);
                    double rotEnergy = rotamerOptimization(molecularAssembly, residues, lowEnergy, optimum);
                    if (rotEnergy < currentEnergy) {
                        currentEnergy = rotEnergy;
                    }
                    if (rotEnergy < lowEnergy) {
                        minRot = i;
                        lowEnergy = rotEnergy;
                    }
                }
                if (minRot > -1) {
                    optimum.add(0, minRot);
                }
            }
        } else {
            /**
             * At the end of the recursion, compute the potential energy for
             * each rotamer of the final residue. If a lower potential energy is
             * discovered, the rotamers of each residue will be collected as the
             * recursion returns up the chain.
             */
            ForceFieldEnergy energy = molecularAssembly.getPotentialEnergy();
            if (rotamers == null) {
                /**
                 * Handle the case where the side-chain has no rotamers.
                 */
                currentEnergy = energy.energy(false, false);
                logger.info(String.format(" Energy: %16.8f", currentEnergy));
                if (currentEnergy < lowEnergy) {
                    optimum.clear();
                    optimum.add(-1);
                }
            } else {
                for (int i = 0; i < rotamers.length; i++) {
                    applyRotamer(name, current, rotamers[i]);
                    double rotEnergy = energy.energy(false, false);
                    logger.info(String.format(" Energy: %16.8f", rotEnergy));
                    if (rotEnergy < currentEnergy) {
                        currentEnergy = rotEnergy;
                    }
                    if (rotEnergy < lowEnergy) {
                        lowEnergy = rotEnergy;
                        optimum.clear();
                        optimum.add(i);
                    }
                }
            }
        }
        residues.add(0, current);
        return currentEnergy;
    }

    public static void applyRotamer(AminoAcid3 name, Residue residue, Rotamer rotamer) {
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
                intxyz(CG1, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CG2, CB, 1.54, CA, 109.5, CG1, 109.5, -1);
                intxyz(HB, CB, 1.11, CA, 109.4, CG1, 109.4, 1);
                intxyz(HG11, CG1, 1.11, CB, 109.4, CA, 180.0, 0);
                intxyz(HG12, CG1, 1.11, CB, 109.4, HG11, 109.4, 1);
                intxyz(HG13, CG1, 1.11, CB, 109.4, HG11, 109.4, -1);
                intxyz(HG21, CG2, 1.11, CB, 109.4, CA, 180.0, 0);
                intxyz(HG22, CG2, 1.11, CB, 109.4, HG21, 109.4, 1);
                intxyz(HG23, CG2, 1.11, CB, 109.4, HG21, 109.4, -1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD1, CG, 1.54, CB, 109.5, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.54, CB, 109.5, CD1, 109.5, -1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG, CG, 1.11, CB, 109.4, CD1, 109.4, 1);
                intxyz(HD11, CD1, 1.11, CG, 109.4, CB, 180.0, 0);
                intxyz(HD12, CD1, 1.11, CG, 109.4, HD11, 109.4, 1);
                intxyz(HD13, CD1, 1.11, CG, 109.4, HD11, 109.4, -1);
                intxyz(HD21, CD2, 1.11, CG, 109.4, CB, 180.0, 0);
                intxyz(HD22, CD2, 1.11, CG, 109.4, HD21, 109.4, 1);
                intxyz(HD23, CD2, 1.11, CG, 109.4, HD21, 109.4, -1);
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
                intxyz(CG1, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CG2, CB, 1.54, CA, 109.5, CG1, 109.5, -1);
                intxyz(CD1, CG1, 1.54, CB, 109.5, CA, rotamer.chi2, 0);
                intxyz(HB, CB, 1.11, CA, 109.4, CG1, 109.4, 1);
                intxyz(HG12, CG1, 1.11, CB, 109.4, CD1, 109.4, 1);
                intxyz(HG13, CG1, 1.11, CB, 109.4, CD1, 109.4, -1);
                intxyz(HG21, CG2, 1.11, CB, 110.0, CG1, 180.0, 0);
                intxyz(HG22, CG2, 1.11, CB, 110.0, HG21, 109.0, 1);
                intxyz(HG23, CG2, 1.11, CB, 110.0, HG21, 109.0, -1);
                intxyz(HD11, CD1, 1.11, CG1, 110.0, CB, 180.0, 0);
                intxyz(HD12, CD1, 1.11, CG1, 110.0, HD11, 109.0, 1);
                intxyz(HD13, CD1, 1.11, CG1, 110.0, HD11, 109.0, -1);
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
                intxyz(OG, CB, 1.41, CA, 107.5, N, rotamer.chi1, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, OG, 106.7, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, OG, 106.7, -1);
                if (rotamer.length == 2) {
                    intxyz(HG, OG, 0.94, CB, 106.9, CA, rotamer.chi2, 0);
                } else {
                    intxyz(HG, OG, 0.94, CB, 106.9, CA, 180.0, 0);
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
                intxyz(OG1, CB, 1.41, CA, 107.5, N, rotamer.chi1, 0);
                intxyz(CG2, CB, 1.54, CA, 109.5, OG1, 107.7, 1);
                intxyz(HB, CB, 1.11, CA, 109.4, OG1, 106.7, -1);
                if (rotamer.length == 2) {
                    intxyz(HG1, OG1, 0.94, CB, 106.9, CA, 180.0, 0);
                } else {
                    intxyz(HG1, OG1, 0.94, CB, 106.9, CA, rotamer.chi2, 0);
                }
                intxyz(HG21, CG2, 1.11, CB, 110.0, CA, 180.0, 0);
                intxyz(HG22, CG2, 1.11, CB, 110.0, HG21, 109.0, 1);
                intxyz(HG23, CG2, 1.11, CB, 110.0, HG21, 109.0, -1);
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
                intxyz(SG, CB, 1.82, CA, 109.0, N, rotamer.chi1, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, SG, 112.0, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, SG, 112.0, -1);
                intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
                break;
            }
            case CYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom SG = (Atom) residue.getAtomNode("SG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                intxyz(SG, CB, 1.82, CA, 109.0, N, rotamer.chi1, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, SG, 112.0, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, SG, 112.0, -1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD1, CG, 1.39, CB, 120.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1);
                intxyz(CE1, CD1, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CE2, CD2, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, CD1, 1.11, CG, 120.0, CE1, 120.0, 1);
                intxyz(HD2, CD2, 1.11, CG, 120.0, CE2, 120.0, 1);
                intxyz(HE1, CE1, 1.11, CD1, 120.0, CZ, 120.0, 1);
                intxyz(HE2, CE2, 1.11, CD2, 120.0, CZ, 120.0, 1);
                intxyz(HZ, CZ, 1.11, CE1, 120.0, CE2, 120.0, 1);
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
                intxyz(CG, CB, 1.54, CA, 107.0, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.54, CB, 107.0, CA, rotamer.chi2, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1);
                intxyz(HD2, CD, 1.11, CG, 109.4, N, 109.4, 1);
                intxyz(HD3, CD, 1.11, CG, 109.4, N, 109.4, -1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD1, CG, 1.39, CB, 120.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1);
                intxyz(CE1, CD1, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CE2, CD2, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0);
                intxyz(OH, CZ, 1.36, CE2, 120.0, CE1, 120.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, CD1, 1.10, CG, 120.0, CE1, 120.0, 1);
                intxyz(HD2, CD2, 1.10, CG, 120.0, CE2, 120.0, 1);
                intxyz(HE1, CE1, 1.10, CD1, 120.0, CZ, 120.0, 1);
                intxyz(HE2, CE2, 1.10, CD2, 120.0, CZ, 120.0, 1);
                if (rotamer.length == 3) {
                    intxyz(HH, OH, 0.97, CZ, 108.0, CE2, rotamer.chi3, 0);
                } else {
                    intxyz(HH, OH, 0.97, CZ, 108.0, CE2, 0.0, 0);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD1, CG, 1.39, CB, 120.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.39, CB, 120.0, CD1, 120.0, 1);
                intxyz(CE1, CD1, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CE2, CD2, 1.39, CG, 120.0, CB, 180, 0);
                intxyz(CZ, CE1, 1.39, CD1, 120.0, CG, 0.0, 0);
                intxyz(OH, CZ, 1.36, CE2, 120.0, CE1, 120.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, CD1, 1.10, CG, 120.0, CE1, 120.0, 1);
                intxyz(HD2, CD2, 1.10, CG, 120.0, CE2, 120.0, 1);
                intxyz(HE1, CE1, 1.10, CD1, 120.0, CZ, 120.0, 1);
                intxyz(HE2, CE2, 1.10, CD2, 120.0, CZ, 120.0, 1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD1, CG, 1.35, CB, 126.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.35, CB, 126.0, CD1, 108.0, 1);
                intxyz(NE1, CD1, 1.35, CG, 108.0, CD2, 0.0, 0);
                intxyz(CE2, NE1, 1.35, CD1, 108.0, CG, 0.0, 0);
                intxyz(CE3, CD2, 1.35, CE2, 120.0, NE1, 180.0, 0);
                intxyz(CZ2, CE2, 1.35, CD2, 120.0, CE3, 0.0, 0);
                intxyz(CZ3, CE3, 1.35, CD2, 120.0, CE2, 0.0, 0);
                intxyz(CH2, CZ2, 1.35, CE2, 120.0, CD2, 0.0, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, CD1, 1.10, CG, 126.0, NE1, 126.0, 1);
                intxyz(HE1, NE1, 1.05, CD1, 126.0, CE2, 126.0, 1);
                intxyz(HE3, CE3, 1.10, CD1, 120.0, CZ3, 120.0, 1);
                intxyz(HZ2, CZ2, 1.10, CE2, 120.0, CH2, 120.0, 1);
                intxyz(HZ3, CZ3, 1.10, CE3, 120.0, CH2, 120.0, 1);
                intxyz(HH2, CH2, 1.10, CZ2, 120.0, CZ3, 120.0, 1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(ND1, CG, 1.35, CB, 126.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1);
                intxyz(CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0);
                intxyz(NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                intxyz(HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1);
                intxyz(HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1);
                intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(ND1, CG, 1.35, CB, 126.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1);
                intxyz(CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0);
                intxyz(NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                intxyz(HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1);
                intxyz(HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1);
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
                intxyz(CG, CB, 1.50, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(ND1, CG, 1.35, CB, 126.0, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, 1.35, CB, 126.0, ND1, 108.0, 1);
                intxyz(CE1, ND1, 1.35, CG, 108.0, CD2, 0.0, 0);
                intxyz(NE2, CD2, 1.35, CG, 108.0, ND1, 0.0, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HD2, CD2, 1.10, CG, 126.0, NE2, 126.0, 1);
                intxyz(HE1, CE1, 1.10, ND1, 126.0, NE2, 126.0, 1);
                intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
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
                intxyz(CG, CB, 1.51, CA, 107.8, N, rotamer.chi1, 0);
                intxyz(OD1, CG, 1.25, CB, 117.0, CA, 0.0, 0);
                intxyz(OD2, CG, 1.25, CB, 117.0, OD1, 126.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1);
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
                intxyz(CG, CB, 1.51, CA, 107.8, N, rotamer.chi1, 0);
                intxyz(OD1, CG, 1.25, CB, 117.0, CA, rotamer.chi2, 0);
                intxyz(OD2, CG, 1.25, CB, 117.0, OD1, 126.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1);
                intxyz(HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0);
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
                intxyz(CG, CB, 1.51, CA, 107.8, N, rotamer.chi1, 0);
                intxyz(OD1, CG, 1.22, CB, 122.5, CA, rotamer.chi2, 0);
                intxyz(ND2, CG, 1.34, CB, 112.7, OD1, 124.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 107.9, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 107.9, -1);
                intxyz(HD21, ND2, 1.02, CG, 119.0, CB, 0.0, 0);
                intxyz(HD22, ND2, 1.02, CG, 119.0, HD21, 120.0, 1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.51, CB, 107.8, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, 1.25, CG, 117.0, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, 1.25, CG, 117.0, OE1, 126.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.51, CB, 107.8, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, 1.25, CG, 117.0, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, 1.25, CG, 117.0, OE1, 126.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1);
                intxyz(HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.51, CB, 107.8, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, 1.22, CG, 122.5, CB, rotamer.chi3, 0);
                intxyz(NE2, CD, 1.34, CG, 112.7, OE1, 124.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 107.9, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 107.9, -1);
                intxyz(HE21, NE2, 1.02, CD, 119.0, CG, 0.0, 0);
                intxyz(HE22, NE2, 1.02, CD, 119.0, HE21, 120.0, 1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(SD, CG, 1.82, CB, 109.0, CA, rotamer.chi2, 0);
                intxyz(CE, SD, 1.82, CG, 96.3, CB, rotamer.chi3, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, SD, 112.0, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, SD, 112.0, -1);
                intxyz(HE1, CE, 1.11, SD, 112.0, CG, 180.0, 0);
                intxyz(HE2, CE, 1.11, SD, 112.0, HE1, 109.4, 1);
                intxyz(HE3, CE, 1.11, SD, 112.0, HE1, 109.4, -1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.54, CB, 109.5, CA, rotamer.chi2, 0);
                intxyz(CE, CD, 1.54, CG, 109.5, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, 1.50, CD, 109.5, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1);
                intxyz(HD2, CD, 1.11, CG, 109.4, CE, 109.4, 1);
                intxyz(HD3, CD, 1.11, CG, 109.4, CE, 109.4, -1);
                intxyz(HE2, CE, 1.11, CD, 109.4, NZ, 108.8, 1);
                intxyz(HE3, CE, 1.11, CD, 109.4, NZ, 108.8, -1);
                intxyz(HZ1, NZ, 1.02, CE, 109.5, CD, 180.0, 0);
                intxyz(HZ2, NZ, 1.02, CE, 109.5, HZ1, 109.5, 1);
                intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.54, CB, 109.5, CA, rotamer.chi2, 0);
                intxyz(CE, CD, 1.54, CG, 109.5, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, 1.50, CD, 109.5, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1);
                intxyz(HD2, CD, 1.11, CG, 109.4, CE, 109.4, 1);
                intxyz(HD3, CD, 1.11, CG, 109.4, CE, 109.4, -1);
                intxyz(HE2, CE, 1.11, CD, 109.4, NZ, 108.8, 1);
                intxyz(HE3, CE, 1.11, CD, 109.4, NZ, 108.8, -1);
                intxyz(HZ1, NZ, 1.02, CE, 109.5, CD, 180.0, 0);
                intxyz(HZ2, NZ, 1.02, CE, 109.5, HZ1, 109.5, 1);
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
                intxyz(CG, CB, 1.54, CA, 109.5, N, rotamer.chi1, 0);
                intxyz(CD, CG, 1.54, CB, 109.5, CA, rotamer.chi2, 0);
                intxyz(NE, CD, 1.45, CG, 109.5, CB, rotamer.chi3, 0);
                intxyz(CZ, NE, 1.35, CD, 120.0, CG, rotamer.chi4, 0);
                intxyz(NH1, CZ, 1.35, NE, 120.0, CD, 180, 0);
                intxyz(NH2, CZ, 1.35, NE, 120.0, NH1, 120.0, 1);
                intxyz(HB2, CB, 1.11, CA, 109.4, CG, 109.4, 1);
                intxyz(HB3, CB, 1.11, CA, 109.4, CG, 109.4, -1);
                intxyz(HG2, CG, 1.11, CB, 109.4, CD, 109.4, 1);
                intxyz(HG3, CG, 1.11, CB, 109.4, CD, 109.4, -1);
                intxyz(HD2, CD, 1.11, CG, 109.4, NE, 109.4, 1);
                intxyz(HD3, CD, 1.11, CG, 109.4, NE, 109.4, -1);
                intxyz(HE, NE, 1.02, CD, 120.0, CZ, 120.0, 1);
                intxyz(HH11, NH1, 1.02, CZ, 120.0, NE, 180.0, 0);
                intxyz(HH12, NH1, 1.02, CZ, 120.0, HH11, 120.0, 1);
                intxyz(HH21, NH2, 1.02, CZ, 120.0, NE, 180.0, 0);
                intxyz(HH22, NH2, 1.02, CZ, 120.0, HH21, 120.0, 1);
                break;
            }
            default:
                break;
        }
    }
}
//    ROTAMER LIBRARY BASED ON STATISTICS FROM PONDER AND RICHARDS
//
//      RESIDUE   CHI1    CHI2    CHI3    CHI4      SIG1  SIG2  SIG3  SIG4
//  1   GLY
//  2   ALA
//  3   VAL      173.5                               9.0
//  4   VAL      -63.4                               8.1
//  5   VAL       69.3                               9.6
//  6   LEU      -64.9   176.0                       8.2   9.9
//  7   LEU     -176.4    63.1                      10.2   8.2
//  8   LEU     -165.3   168.2                      10.0  34.2
//  9   LEU       44.3    60.4                      20.0  18.8
// 10   ILE      -60.9   168.7                       7.5  11.6
// 11   ILE      -59.6   -64.1                       9.6  14.3
// 12   ILE       61.7   163.8                       5.0  16.4
// 13   ILE     -166.6   166.0                      10.1   8.9
// 14   ILE     -174.8    72.1                      24.9  10.5
// 15   SER       64.7                              16.1
// 16   SER      -69.7                              14.6
// 17   SER     -176.1                              20.2
// 18   THR       62.7                               8.5
// 19   THR      -59.7                               9.4
// 20   THR     -169.5                               6.6
// 21   CYS      -65.2                              10.1
// 22   CYS     -179.6                               9.5
// 23   CYS       63.5                               9.6
// 24   PRO       24.0                               8.0
// 25   PRO        0.0                               8.0
// 26   PRO      -24.0                               8.0
// 27   PHE      -66.3    94.3                      10.2  19.5
// 28   PHE     -179.2    78.9                       9.3   8.9
// 29   PHE       66.0    90.7                      12.0   9.4
// 30   PHE      -71.9    -0.4                      16.3  26.1
// 31   TYR      -66.5    96.6                      11.4  21.8
// 32   TYR     -179.7    71.9                      12.6  13.4
// 33   TYR       63.3    89.1                       9.4  13.0
// 34   TYR      -67.2    -1.0                      13.2  20.1
// 35   TRP      -70.4   100.5                       7.0  18.2
// 36   TRP       64.8   -88.9                      13.0   5.3
// 37   TRP     -177.3   -95.1                       7.9   7.6
// 38   TRP     -179.5    87.5                       3.4   3.8
// 39   TRP      -73.3   -87.7                       6.5   8.1
// 40   TRP       62.2   112.5                      10.0  15.0
// 41   HIS      -62.8   -74.3                      10.0  17.2
// 42   HIS     -175.2   -87.7                      15.4  43.5
// 43   HIS      -69.8    96.1                       5.9  32.2
// 44   HIS       67.9   -80.5                      17.4  40.7
// 45   HIS     -177.3   100.5                       6.3  14.0
// 46   HIS       48.8    89.5                      10.0  30.0
// 47   ASP      -68.3   -25.7                       9.2  31.1
// 48   ASP     -169.1     3.9                       9.5  38.9
// 49   ASP       63.7     2.4                       9.9  29.4
// 50   ASN      -68.3   -36.8                      12.3  25.2
// 51   ASN     -177.1     1.3                       8.8  34.1
// 52   ASN      -67.2   128.8                      10.8  24.2
// 53   ASN       63.9    -6.8                       3.7  13.5
// 54   ASN     -174.9  -156.8                      17.9  58.9
// 55   ASN       63.6    53.8                       6.6  17.1
// 56   GLU      -69.6  -177.2   -11.4              19.2  21.7  44.8
// 57   GLU     -176.2   175.4    -6.7              14.9  10.6  39.0
// 58   GLU      -64.6   -69.1   -33.4              13.5  17.3  27.4
// 59   GLU      -55.6    77.0    25.3              10.6   6.8  32.6
// 60   GLU       69.8  -179.0     6.6              10.6  23.7  64.2
// 61   GLU     -173.6    70.6    14.0              14.6   8.7  37.1
// 62   GLU       63.0   -80.4    16.3               4.3  13.9  20.8
// 63   GLN      -66.7  -178.5   -24.0              14.1  14.9  38.0
// 64   GLN      -66.7  -178.5   156.0              14.1  14.9  38.0
// 65   GLN     -174.6  -177.7   -24.0              11.5  17.2  38.0
// 66   GLN     -174.6  -177.7   156.0              11.5  17.2  38.0
// 67   GLN      -58.7   -63.8   -46.3              11.2  16.1  27.7
// 68   GLN      -51.3   -90.4   165.0               7.3  22.8  38.2
// 69   GLN     -179.4    67.3    26.8              21.5   7.9  38.4
// 70   GLN      167.5    70.9   174.2              14.8   3.7   7.1
// 71   GLN       70.8  -165.6   -24.0              13.0   9.5  38.0
// 72   GLN       70.8  -165.6   156.0              13.0   9.5  38.0
// 73   MET      -64.5   -68.5   -75.6              12.7   6.0  14.1
// 74   MET      -78.3  -174.7    65.0               5.4  15.7  20.0
// 75   MET      -78.3  -174.7   180.0               5.4  15.7  20.0
// 76   MET      -78.3  -174.7   -65.0               5.4  15.7  20.0
// 77   MET      178.9   179.0    65.0               8.7  13.4  20.0
// 78   MET      178.9   179.0   180.0               8.7  13.4  20.0
// 79   MET      178.9   179.0   -65.0               8.7  13.4  20.0
// 80   MET      -70.0   -65.0    65.0              21.0  20.0  20.0
// 81   MET     -170.0    65.0   180.0              24.0  20.0  20.0
// 82   MET     -170.0   -65.0   180.0              24.0  20.0  20.0
// 83   MET      -70.0    65.0   180.0              21.0  20.0  20.0
// 84   MET      -70.0   -65.0   180.0              21.0  20.0  20.0
// 85   MET       61.0    65.0   180.0              21.0  20.0  20.0
// 86   LYS     -170.0   180.0    65.0   180.0      24.0  20.0  20.0  20.0
// 87   LYS     -170.0   180.0   180.0    65.0      24.0  20.0  20.0  20.0
// 88   LYS     -170.0   180.0   180.0   180.0      24.0  20.0  20.0  20.0
// 89   LYS     -170.0   180.0   -65.0   180.0      24.0  20.0  20.0  20.0
// 90   LYS     -170.0   -65.0   180.0   180.0      24.0  20.0  20.0  20.0
// 91   LYS      -70.0    65.0   180.0   180.0      21.0  20.0  20.0  20.0
// 92   LYS      -70.0   180.0    65.0   180.0      21.0  20.0  20.0  20.0
// 93   LYS      -70.0   180.0   180.0   180.0      21.0  20.0  20.0  20.0
// 94   LYS      -70.0   180.0   180.0   -65.0      21.0  20.0  20.0  20.0
// 95   LYS      -70.0   180.0   -65.0   180.0      21.0  20.0  20.0  20.0
// 96   LYS      -70.0   -65.0   180.0   180.0      21.0  20.0  20.0  20.0
// 97   LYS      -70.0   -65.0   180.0   -65.0      21.0  20.0  20.0  20.0
// 98   ARG       61.0   180.0    65.0    90.0      25.0  20.0  20.0  20.0
// 99   ARG       61.0   180.0   180.0   180.0      25.0  20.0  20.0  20.0
//100   ARG     -170.0   180.0    65.0    90.0      24.0  20.0  20.0  20.0
//101   ARG     -170.0   180.0   180.0    90.0      24.0  20.0  20.0  20.0
//102   ARG     -170.0   180.0   180.0   180.0      24.0  20.0  20.0  20.0
//103   ARG     -170.0   180.0   180.0   -90.0      24.0  20.0  20.0  20.0
//104   ARG     -170.0   180.0   -65.0   180.0      24.0  20.0  20.0  20.0
//105   ARG      -70.0   180.0    65.0    90.0      21.0  20.0  20.0  20.0
//106   ARG      -70.0   180.0    65.0   180.0      21.0  20.0  20.0  20.0
//107   ARG      -70.0   180.0   180.0   180.0      21.0  20.0  20.0  20.0
//108   ARG      -70.0   180.0   180.0   -90.0      21.0  20.0  20.0  20.0
//109   ARG      -70.0   180.0   -65.0   180.0      21.0  20.0  20.0  20.0
//110   ARG     -170.0    65.0    65.0   180.0      21.0  20.0  20.0  20.0
//111   ARG      -70.0   -65.0   -65.0   180.0      21.0  20.0  20.0  20.0

