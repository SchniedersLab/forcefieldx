/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import static ffx.potential.bonded.BondedUtils.intxyz;

/**
 * The INTFilter class parses TINKER internal coordinate (*.INT) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class INTFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(INTFilter.class.getName());

    /**
     * <p>
     * Constructor for INTFilter.</p>
     *
     * @param files a {@link java.util.List} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public INTFilter(List<File> files, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(files, molecularAssembly, forceField, properties);
        fileType = FileType.INT;
    }

    /**
     * <p>
     * Constructor for INTFilter.</p>
     *
     * @param file a {@link java.io.File} object.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public INTFilter(File file, MolecularAssembly molecularAssembly,
                     ForceField forceField, CompositeConfiguration properties) {
        super(file, molecularAssembly, forceField, properties);
        fileType = FileType.INT;
    }

    /**
     * {@inheritDoc}
     *
     * Parse the INT File.
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

    /** {@inheritDoc} */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        /*
         * Would need to add something for writing reduced hydrogen coordinates.
         *
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

    /** {@inheritDoc} */
    @Override
    public boolean readNext() {
        return readNext(false);
    }

    /** {@inheritDoc} */
    @Override
    public boolean readNext(boolean resetPosition) {
        return false;
    }

    /** {@inheritDoc} */
    @Override
    public boolean readNext(boolean resetPosition, boolean print) {
        return false;
    }

    /** {@inheritDoc} */
    @Override
    public void closeReader() {
        logger.fine(" Reading trajectories not yet supported for INTFilter");
    }
}
