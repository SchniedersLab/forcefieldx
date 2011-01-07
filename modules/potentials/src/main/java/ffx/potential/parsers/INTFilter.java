/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import static java.lang.Math.*;

import static ffx.numerics.VectorMath.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import java.util.List;

/**
 * The INTFilter class parses TINKER internal coordinate (*.INT) files.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class INTFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(INTFilter.class.getName());
    private static double x[] = new double[3];
    private static double xa[] = new double[3];
    private static double xb[] = new double[3];
    private static double xc[] = new double[3];
    private static double xab[] = new double[3];
    private static double xba[] = new double[3];
    private static double xbc[] = new double[3];
    private static double xac[] = new double[3];
    private static double xt[] = new double[3];
    private static double xu[] = new double[3];
    private static double rab, cosb, sinb, cosg, sing;
    private static double sine, sine2, cosine;
    private static double a, b, c;
    private static double xtmp, ztmp;
    private static double eps = 0.0000001d;

    public INTFilter(List<File> files, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(files, molecularAssembly, forceField, properties);
        fileType = FileType.INT;
    }

    public INTFilter(File file, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssembly, forceField, properties);
        fileType = FileType.INT;
    }

    /**
     * This routine was derived from a similar routine in TINKER.
     */
    public static void intxyz(Atom atom, Atom ia, double bond, Atom ib, double angle1, Atom ic, double angle2, int chiral) {
        angle1 = toRadians(angle1);
        angle2 = toRadians(angle2);
        double zcos0 = cos(angle1);
        double zcos1 = cos(angle2);
        double zsin0 = sin(angle1);
        double zsin1 = sin(angle2);
        // No partners
        if (ia == null) {
            atom.moveTo(0.0d, 0.0d, 0.0d);
        } else if (ib == null) {
            // One partner - place on the z-axis
            ia.getXYZ(xa);
            xa[2] += bond;
            atom.moveTo(xa);
        } else if (ic == null) {
            // Two partners - place in the xz-plane
            ia.getXYZ(xa);
            ib.getXYZ(xb);
            diff(xa, xb, xab);
            rab = r(xab);
            norm(xab, xab);
            cosb = xab[2];
            sinb = sqrt(xab[0] * xab[0] + xab[1] * xab[1]);
            if (sinb == 0.0d) {
                cosg = 1.0d;
                sing = 0.0d;
            } else {
                cosg = xab[1] / sinb;
                sing = xab[0] / sinb;
            }
            xtmp = bond * zsin0;
            ztmp = rab - bond * zcos0;
            x[0] = xb[0] + xtmp * cosg + ztmp * sing * sinb;
            x[1] = xb[1] - xtmp * sing + ztmp * cosg * sinb;
            x[2] = xb[2] + ztmp * cosb;
            atom.moveTo(x);
        } else if (chiral == 0) {
            // General case - with a dihedral
            ia.getXYZ(xa);
            ib.getXYZ(xb);
            ic.getXYZ(xc);
            diff(xa, xb, xab);
            norm(xab, xab);
            diff(xb, xc, xbc);
            norm(xbc, xbc);
            xt[0] = xab[2] * xbc[1] - xab[1] * xbc[2];
            xt[1] = xab[0] * xbc[2] - xab[2] * xbc[0];
            xt[2] = xab[1] * xbc[0] - xab[0] * xbc[1];
            cosine = xab[0] * xbc[0] + xab[1] * xbc[1] + xab[2] * xbc[2];
            sine = sqrt(max(1.0d - cosine * cosine, eps));
            if (abs(cosine) >= 1.0d) {
                logger.warning("Undefined Dihedral");
            }
            scalar(xt, 1.0d / sine, xt);
            xu[0] = xt[1] * xab[2] - xt[2] * xab[1];
            xu[1] = xt[2] * xab[0] - xt[0] * xab[2];
            xu[2] = xt[0] * xab[1] - xt[1] * xab[0];
            x[0] = xa[0] + bond * (xu[0] * zsin0 * zcos1 + xt[0] * zsin0 * zsin1 - xab[0] * zcos0);
            x[1] = xa[1] + bond * (xu[1] * zsin0 * zcos1 + xt[1] * zsin0 * zsin1 - xab[1] * zcos0);
            x[2] = xa[2] + bond * (xu[2] * zsin0 * zcos1 + xt[2] * zsin0 * zsin1 - xab[2] * zcos0);
            atom.moveTo(x);
        } else if (abs(chiral) == 1) {
            ia.getXYZ(xa);
            ib.getXYZ(xb);
            ic.getXYZ(xc);
            diff(xb, xa, xba);
            norm(xba, xba);
            diff(xa, xc, xac);
            norm(xac, xac);
            xt[0] = xba[2] * xac[1] - xba[1] * xac[2];
            xt[1] = xba[0] * xac[2] - xba[2] * xac[0];
            xt[2] = xba[1] * xac[0] - xba[0] * xac[1];
            cosine = xba[0] * xac[0] + xba[1] * xac[1] + xba[2] * xac[2];
            sine2 = max(1.0d - cosine * cosine, eps);
            if (abs(cosine) >= 1.0d) {
                logger.warning("Defining Atom Colinear");
            }
            a = (-zcos1 - cosine * zcos0) / sine2;
            b = (zcos0 + cosine * zcos1) / sine2;
            c = (1.0d + a * zcos1 - b * zcos0) / sine2;
            if (c > eps) {
                c = chiral * sqrt(c);
            } else if (c < -eps) {
                c = sqrt((a * xac[0] + b * xba[0]) * (a * xac[0] + b * xba[0]) + (a * xac[1] + b * xba[1]) * (a * xac[1] + b * xba[1]) + (a * xac[2] + b * xba[2]) * (a * xac[2] + b * xba[2]));
                a /= c;
                b /= c;
                c = 0.0d;
            } else {
                c = 0.0d;
            }
            x[0] = xa[0] + bond * (a * xac[0] + b * xba[0] + c * xt[0]);
            x[1] = xa[1] + bond * (a * xac[1] + b * xba[1] + c * xt[1]);
            x[2] = xa[2] + bond * (a * xac[2] + b * xba[2] + c * xt[2]);
            atom.moveTo(x);
        }
    }

    /**
     * Parse the INT File.
     *
     * @return Returns true on successful read, false otherwise.
     *
     * @since 1.0
     */
    @Override
    public boolean readFile() {
        File intFile = activeMolecularAssembly.getFile();
        if (forceField == null) {
            logger.warning("No force field is associated with " + intFile.toString());
            return false;
        }
        // Open a data stream to the Internal Coordinate file
        try {
            FileReader fr = new FileReader(intFile);
            BufferedReader br = new BufferedReader(fr);
            String data = br.readLine().trim();
            // Read blank lines at the top of the file
            while (data != null && data.length() == 0) {
                data = br.readLine().trim();
            }
            if (data == null) {
                logger.warning("Empty file: " + intFile.toString());
                return false;
            }
            int numberOfAtoms;
            String tokens[] = data.trim().split(" +");
            try {
                numberOfAtoms = Integer.parseInt(tokens[0]);
                if (numberOfAtoms < 1) {
                    logger.warning("Invalid number of atoms: " + numberOfAtoms);
                    return false;
                }
            } catch (Exception e) {
                logger.severe("Error parsing the number of atoms.\n" + e);
                return false;
            }
            if (tokens.length >= 2) {
                tokens = data.trim().split(" +", 2);
                activeMolecularAssembly.setName(tokens[1]);
            }
            logger.info("  Opening " + intFile.getName() + " with " + numberOfAtoms + " atoms");
            double d[] = {0.0d, 0.0d, 0.0d};
            int zi[][] = new int[numberOfAtoms][4];
            double zv[][] = new double[numberOfAtoms][3];
            Vector<int[]> zadd = new Vector<int[]>();
            Vector<int[]> zdel = new Vector<int[]>();
            atomList = new ArrayList<Atom>();
            for (int i = 0; i < numberOfAtoms; i++) {
                // Atom Data
                if (!br.ready()) {
                    return false;
                }
                data = br.readLine();
                if (data == null) {
                    logger.severe("  Check atom " + (i + 1) + " in " + activeMolecularAssembly.getFile().getName());
                    return false;
                }
                tokens = data.trim().split(" +");
                if (tokens == null || tokens.length < 3) {
                    logger.severe("  Check atom " + (i + 1) + " in " + activeMolecularAssembly.getFile().getName());
                    return false;
                }
                // Atom number, name, type
                String name = tokens[1];
                int type = Integer.parseInt(tokens[2]);
                AtomType atomType = forceField.getAtomType(String.valueOf(type));
                if (atomType == null) {
                    logger.severe("  Check atom " + (i + 1) + " in " + activeMolecularAssembly.getFile().getName());
                    return false;
                }
                Atom atom = new Atom(i + 1, name, atomType, d);
                atomList.add(atom);
                // Bond partner and bond value
                if (tokens.length >= 5) {
                    zi[i][0] = Integer.parseInt(tokens[3]);
                    zv[i][0] = Double.parseDouble(tokens[4]);
                } else {
                    zi[i][0] = 0;
                    zv[i][0] = 0.0d;
                }
                // Angle partner and angle value
                if (tokens.length >= 7) {
                    zi[i][1] = Integer.parseInt(tokens[5]);
                    zv[i][1] = Double.parseDouble(tokens[6]);
                } else {
                    zi[i][1] = 0;
                    zv[i][1] = 0.0d;
                }
                // Torsion partner and dihedral value
                if (tokens.length >= 10) {
                    zi[i][2] = Integer.parseInt(tokens[7]);
                    zv[i][2] = Double.parseDouble(tokens[8]);
                    zi[i][3] = Integer.parseInt(tokens[9]);
                } else {
                    zi[i][2] = 0;
                    zv[i][2] = 0.0d;
                    zi[i][3] = 0;
                }
            }
            if (br.ready()) {
                data = br.readLine();
                // Check for a first blank line
                if (data.trim().equalsIgnoreCase("")) {
                    // Parse bond pairs to add until EOF or a blank line is
                    // reached
                    boolean blank = false;
                    while (br.ready() && !blank) {
                        data = br.readLine();
                        if (data.trim().equalsIgnoreCase("")) {
                            blank = true;
                        } else {
                            tokens = data.trim().split(" +");
                            if (tokens.length != 2) {
                                logger.severe("  Check Additional Bond Pair: " + (zadd.size() + 1) + " in " + activeMolecularAssembly.getFile().getName());
                                return false;
                            }
                            int pair[] = new int[2];
                            pair[0] = Integer.parseInt(tokens[0]);
                            pair[1] = Integer.parseInt(tokens[1]);
                            zadd.add(pair);
                        }
                    }
                    // Parse bond pairs to be removed until EOF
                    while (br.ready()) {
                        data = br.readLine();
                        tokens = data.trim().split(" +");
                        if (tokens.length != 2) {
                            logger.severe("  Check Bond Pair to Remove: " + (zadd.size() + 1) + " in " + activeMolecularAssembly.getFile().getName());
                            return false;
                        }
                        int pair[] = new int[2];
                        pair[0] = Integer.parseInt(tokens[0]);
                        pair[1] = Integer.parseInt(tokens[1]);
                        zdel.add(pair);
                    }
                }
            }
            br.close();
            fr.close();
            if (atomList.size() == numberOfAtoms) {
                // Add bonds specified in the Z-matrix
                bondList = new ArrayList<Bond>();
                for (int i = 1; i < numberOfAtoms; i++) {
                    int partner = zi[i][0];
                    boolean del = false;
                    for (int j = 0; j < zdel.size(); j++) {
                        int pair[] = zdel.get(j);
                        if (pair[0] == i + 1 && pair[1] == partner) {
                            del = true;
                        }
                        if (pair[1] == i + 1 && pair[0] == partner) {
                            del = true;
                        }
                    }
                    if (!del) {
                        Atom atom1 = atomList.get(i);
                        Atom atom2 = atomList.get(partner - 1);
                        bondList.add(new Bond(atom1, atom2));
                    }
                }
                // Add additional bonds
                for (int i = 0; i < zadd.size(); i++) {
                    int pair[] = zadd.get(i);
                    Atom atom1 = atomList.get(pair[0] - 1);
                    Atom atom2 = atomList.get(pair[1] - 1);
                    bondList.add(new Bond(atom1, atom2));
                }
                // Determine coordinates from Z-matrix values
                for (int i = 0; i < numberOfAtoms; i++) {
                    Atom atom = atomList.get(i);
                    Atom ia = null;
                    Atom ib = null;
                    Atom ic = null;
                    int[] atoms = zi[i];
                    if (atoms[0] > 0) {
                        ia = atomList.get(atoms[0] - 1);
                    }
                    if (atoms[1] > 0) {
                        ib = atomList.get(atoms[1] - 1);
                    }
                    if (atoms[2] > 0) {
                        ic = atomList.get(atoms[2] - 1);
                    }
                    double bond = zv[i][0];
                    double angle1 = zv[i][1];
                    double angle2 = zv[i][2];
                    int chiral = atoms[3];
                    intxyz(atom, ia, bond, ib, angle1, ic, angle2, chiral);
                }
                return true;
            }
            logger.warning("Reported number of Atoms: " + numberOfAtoms + "\nNumber of Atoms Found: " + atomList.size());
        } catch (IOException e) {
            logger.severe(e.toString());
        }
        return false;
    }

    @Override
    public boolean writeFile(File saveFile, boolean append) {
        /*
         * File xyzfile = getFile(); if (xyzfile == null) { return false; } try
         * { FileWriter fw = new FileWriter(xyzfile); BufferedWriter bw = new
         * BufferedWriter(fw); // XYZ File First Line FSystem M = getFSystem();
         * int numatoms = M.getAtomList().size(); String blanks = new
         * String(" "); int len = (new String("" + numatoms)).length();
         * bw.write(blanks.substring(0, 6 - len) + numatoms + " " + M.toString()
         * + "\n"); Atom a, a2; Bond b; ArrayList bonds; StringBuilder line;
         * StringBuilder lines[] = new StringBuilder[numatoms]; String indexS, id,
         * type, xS, yS, zS; int xi, yi, zi; // XYZ File Atom Lines List atoms =
         * M.getAtomList(); Vector3d offset = M.getOffset(); for (ListIterator
         * li = atoms.listIterator(); li.hasNext(); ) { a = (Atom) li.next();
         * indexS = new String("" + a.getXYZIndex()); line = new
         * StringBuilder(blanks.substring(0, 6 - indexS.length()) + indexS +
         * " "); id = a.getID(); line.append(id + blanks.substring(0, 3 -
         * id.length())); xS = formatCoord.format(a.getX() - offset.x); yS =
         * formatCoord.format(a.getY() - offset.y); zS =
         * formatCoord.format(a.getZ() - offset.z);
         * line.append(blanks.substring(0, 12 - xS.length()) + xS);
         * line.append(blanks.substring(0, 12 - yS.length()) + yS);
         * line.append(blanks.substring(0, 12 - zS.length()) + zS); type = new
         * String("" + a.getAtomType()); line.append(blanks.substring(0, 6 -
         * type.length()) + type); bonds = a.getBonds(); if (bonds != null) {
         * for (ListIterator lj = bonds.listIterator(); lj.hasNext(); ) { b =
         * (Bond) lj.next(); a2 = b.getOtherAtom(a); xS =
         * formatBond.format(a2.xyzindex); line.append(blanks.substring(0, 6 -
         * xS.length()) + xS); } lines[a.getXYZIndex() - 1] = line.append("\n");
         * } } for (int i = 0; i < numatoms; i++) { try {
         * bw.write(lines[i].toString()); } catch (Exception e) {
         * System.out.println("" + i); } } bw.close(); fw.close(); } catch
         * (IOException e) { System.out.println("" + e); return false; } return
         * true;
         */
        return false;
    }
}
