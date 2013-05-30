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

import ffx.potential.ResidueEnumerations.AminoAcid3;
import static ffx.potential.ResidueEnumerations.AminoAcid3.ASN;
import static ffx.potential.ResidueEnumerations.AminoAcid3.GLN;
import static ffx.potential.ResidueEnumerations.AminoAcid3.GLU;
import static ffx.potential.ResidueEnumerations.AminoAcid3.LYS;
import static ffx.potential.ResidueEnumerations.AminoAcid3.MET;
import static ffx.potential.ResidueEnumerations.AminoAcid3.TRP;

/**
 * The Rotamer Library Class manages a library of side-chain Rotamers.
 * @author Ava M. Lynn
 */
public class RotamerLibrary {
    private static final int numberOfAminoAcids = 20;
    private static final Rotamer[][] rotamerCache = new Rotamer[numberOfAminoAcids][];
    public static Rotamer[] getRotamer(AminoAcid3 name){
        switch(name){
            case GLY:
            case ALA:
                return null;
            case VAL:
                int n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, 173.5, 9.0);
                    rotamerCache[n][1] = new Rotamer(name, -63.4, 8.1);
                    rotamerCache[n][2] = new Rotamer(name, 69.3, 9.6);
                }
                return rotamerCache[n];
            case LEU:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -64.9, 8.2, 176.0, 9.9);
                    rotamerCache[n][1] = new Rotamer(name, -176.4, 10.2, 63.1, 8.2);
                    rotamerCache[n][2] = new Rotamer(name, -165.3, 10.0, 168.2, 34.2);
                    rotamerCache[n][3] = new Rotamer(name, 44.3, 20.0, 60.4, 18.8);
                }
                return rotamerCache[n];
            case ILE:
                    n = name.ordinal();
                    if (rotamerCache[n] == null) {
                        // Allocate memory and create rotamers
                        rotamerCache[n][0] = new Rotamer(name, -60.9, 7.5, 168.7, 11.6);
                        rotamerCache[n][1] = new Rotamer(name, -59.6, 9.6, -64.1, 14.3);
                        rotamerCache[n][2] = new Rotamer(name, 61.7, 5.0, 163.8, 16.4);
                        rotamerCache[n][3] = new Rotamer(name, -166.6, 10.1, 166.0, 8.9);
                        rotamerCache[n][4] = new Rotamer(name, -174.8, 24.9, 72.1, 10.5);
                }
                return rotamerCache[n];
            case SER:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, 64.7, 16.1);
                    rotamerCache[n][1] = new Rotamer(name, -69.7, 14.6);
                    rotamerCache[n][2] = new Rotamer(name, -176.1, 20.2);
                }
                return rotamerCache[n];
            case THR:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, 62.7, 8.5);
                    rotamerCache[n][1] = new Rotamer(name, -59.7, 9.4);
                    rotamerCache[n][2] = new Rotamer(name, -169.5, 6.6);
                }
                return rotamerCache[n];
            case CYS:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -65.2, 10.1);
                    rotamerCache[n][1] = new Rotamer(name, -179.6, 9.5);
                    rotamerCache[n][2] = new Rotamer(name, 63.5, 9.6);
                }
                return rotamerCache[n];
            case PRO:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, 24.0, 8.0);
                    rotamerCache[n][1] = new Rotamer(name, 0.0, 8.0);
                    rotamerCache[n][2] = new Rotamer(name, -24.0, 8.0);
                }
                return rotamerCache[n];
            case PHE:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -66.3, 10.2, 94.3, 19.5);
                    rotamerCache[n][1] = new Rotamer(name, -179.2, 9.3, 7.8, 8.9);
                    rotamerCache[n][2] = new Rotamer(name, 66.0, 12.0, 90.7, 9.4);
                    rotamerCache[n][3] = new Rotamer(name, -71.9, 16.3, -0.4, 26.1);
                }
                return rotamerCache[n];
            case TYR:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -66.5, 11.4, 96.6, 21.8);
                    rotamerCache[n][1] = new Rotamer(name, -179.7, 12.6, 71.9, 13.4);
                    rotamerCache[n][2] = new Rotamer(name, 63.3, 9.4, 89.1, 13.0);
                    rotamerCache[n][3] = new Rotamer(name, -67.2, 13.2, -1.0, 20.1);
                }
                return rotamerCache[n];
            case TRP:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -70.4, 7.0, 100.5, 18.2);
                    rotamerCache[n][1] = new Rotamer(name, 64.8, 13.0, -88.9, 5.3);
                    rotamerCache[n][2] = new Rotamer(name, -177.3, 7.9, -95.1, 7.6);
                    rotamerCache[n][3] = new Rotamer(name, -179.5, 3.4, 87.5, 3.8);
                    rotamerCache[n][4] = new Rotamer(name, -73.3, 6.5, -87.7, 8.1);
                    rotamerCache[n][5] = new Rotamer(name, 62.2, 10.0, 112.5, 15.0);
                }
                return rotamerCache[n];
            case HIS:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -62.8, 10.0, -74.3, 17.2);
                    rotamerCache[n][1] = new Rotamer(name, -175.2, 15.4, -88.7, 43.5);
                    rotamerCache[n][2] = new Rotamer(name, -69.8, 5.9, 96.1, 32.2);
                    rotamerCache[n][3] = new Rotamer(name, 67.9, 17.4, -80.5, 40.7);
                    rotamerCache[n][4] = new Rotamer(name, -177.3, 6.3, 100.5, 14.0);
                    rotamerCache[n][5] = new Rotamer(name, 48.8, 10.0, 89.5, 30.0);
                }
                return rotamerCache[n];
            case ASH:
            case ASP:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -68.3, 9.2, -25.7, 31.1);
                    rotamerCache[n][1] = new Rotamer(name, -169.1, 9.5, 3.9, 38.9);
                    rotamerCache[n][2] = new Rotamer(name, 63.7, 9.9, 2.4, 29.4);
                }
                return rotamerCache[n];
            case ASN:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -68.3, 12.3, -36.8, 25.2);
                    rotamerCache[n][1] = new Rotamer(name, -177.1, 8.8, 1.3, 34.1);
                    rotamerCache[n][2] = new Rotamer(name, -67.2, 10.8, 128.8, 24.2);
                    rotamerCache[n][3] = new Rotamer(name, 63.9, 3.7, -6.8, 13.5);
                    rotamerCache[n][4] = new Rotamer(name, -174.9, 17.9, -156.8, 58.9);
                    rotamerCache[n][5] = new Rotamer(name, 63.6, 6.6, 53.8, 17.1);
                }
                return rotamerCache[n];
            case GLU:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
                    rotamerCache[n][0] = new Rotamer(name, -69.6, 19.2, -177.2, 21.7, -11.4, 44.8);
                    rotamerCache[n][1] = new Rotamer(name, -176.2, 14.9, 175.4, 10.6, -6.7, 39.0);
                    rotamerCache[n][2] = new Rotamer(name, -64.6, 13.5, -69.1, 17.3, -33.4, 27.4);
                    rotamerCache[n][3] = new Rotamer(name, -55.6, 10.6, 77.0, 6.8, 25.3, 32.6);
                    rotamerCache[n][4] = new Rotamer(name, 69.8, 10.6, -179.0, 23.7, 6.6, 64.2);
                    rotamerCache[n][5] = new Rotamer(name, -173.6, 14.6, 70.6, 8.7, 14.0, 37.1);
                    rotamerCache[n][6] = new Rotamer(name, 63.0, 4.3, -80.4, 13.9, 16.3, 20.8);
                }
                return rotamerCache[n];
            case GLN:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
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
                }
                return rotamerCache[n];
            case MET:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
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
                }
                return rotamerCache[n];
            case LYS:
                n = name.ordinal();
                if (rotamerCache[n] == null) {
                    // Allocate memory and create rotamers
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
                }
                return rotamerCache[n];
        case ARG:
            n = name.ordinal();
            if (rotamerCache[n] == null) {
                // Allocate memory and create rotamers
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
            }
            return rotamerCache[n];
            default:
                // TODO
        }
        return null;
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
    
}
