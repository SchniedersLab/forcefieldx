/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;
import java.util.logging.Logger;

import javax.vecmath.Vector3d;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;

/**
 * The XYZFilter class parses TINKER Cartesian coordinate (*.XYZ) files.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class XYZFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(XYZFilter.class.getName());

    public static boolean readOnto(File newFile, MolecularAssembly oldSystem) {
        if (newFile == null || !newFile.exists() || oldSystem == null) {
            return false;
        }
        try {
            FileReader fr = new FileReader(newFile);
            BufferedReader br = new BufferedReader(fr);
            String data = br.readLine();
            if (data == null) {
                return false;
            }
            String tokens[] = data.trim().split(" +");
            int num_atoms = Integer.parseInt(tokens[0]);
            if (num_atoms != oldSystem.getAtomList().size()) {
                return false;
            }
            double d[][] = new double[num_atoms][3];
            for (int i = 0; i < num_atoms; i++) {
                if (!br.ready()) {
                    return false;
                }
                data = br.readLine();
                if (data == null) {
                    logger.warning("Check atom " + (i + 1));
                    return false;
                }
                tokens = data.trim().split(" +");
                if (tokens == null || tokens.length < 6) {
                    logger.warning("Check atom " + (i + 1));
                    return false;
                }
                d[i][0] = Double.parseDouble(tokens[2]);
                d[i][1] = Double.parseDouble(tokens[3]);
                d[i][2] = Double.parseDouble(tokens[4]);
            }
            ArrayList<Atom> atoms = oldSystem.getAtomList();
            for (Atom a : atoms) {
                int index = a.getXYZIndex() - 1;
                a.setXYZ(d[index]);
            }
            oldSystem.center();
            oldSystem.setFile(newFile);
            br.close();
            fr.close();
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    public XYZFilter() {
        super();
        setType(FileType.XYZ);
    }

    public XYZFilter(MolecularAssembly system) {
        super(system);
        setType(FileType.XYZ);
    }

    public XYZFilter(MolecularAssembly system, ForceField forceField) {
        super(system, forceField);
        setType(FileType.XYZ);
    }

    /**
     * Parse the XYZ File
     */
    @Override
    public boolean readFile() {
        File xyzFile = molecularAssembly.getFile();
        if (forceField == null) {
            logger.warning("No force field is associated with " + xyzFile.toString());
            return false;
        }
        try {
            FileReader fr = new FileReader(xyzFile);
            BufferedReader br = new BufferedReader(fr);
            String data = br.readLine();
            // Read blank lines at the top of the file
            while (data != null && data.trim().equals("")) {
                data = br.readLine();
            }
            if (data == null) {
                return false;
            }
            String tokens[] = data.trim().split(" +", 2);
            int numberOfAtoms = Integer.parseInt(tokens[0]);
            if (numberOfAtoms < 1) {
                return false;
            }
            if (tokens.length == 2) {
                getMolecularSystem().setName(tokens[1]);
            }
            logger.info("\n Opening " + xyzFile.getName() + " with " + numberOfAtoms + " atoms\n");
            // The header line is reasonable - prepare to parse atom lines.
            Hashtable<Integer, Integer> labelHash = new Hashtable<Integer, Integer>();
            int label[] = new int[numberOfAtoms];
            int bonds[][] = new int[numberOfAtoms][8];
            double d[][] = new double[numberOfAtoms][3];
            boolean renumber = false;
            atomList = new ArrayList<Atom>();
            // Loop over the expected number of atoms.
            for (int i = 0; i < numberOfAtoms; i++) {
                if (!br.ready()) {
                    return false;
                }
                data = br.readLine();
                if (data == null) {
                    logger.warning("Check atom " + (i + 1) + " in " + molecularAssembly.getFile().getName());
                    return false;
                }
                tokens = data.trim().split(" +");
                if (tokens == null || tokens.length < 6) {
                    logger.warning("Check atom " + (i + 1) + " in " + molecularAssembly.getFile().getName());
                    return false;
                }
                // Valid number of tokens, so try to parse this line.
                label[i] = Integer.parseInt(tokens[0]);
                // Check for valid atom numbering, or flag for re-numbering.
                if (label[i] != i + 1) {
                    renumber = true;
                }
                String atomName = tokens[1];
                d[i][0] = Double.parseDouble(tokens[2]);
                d[i][1] = Double.parseDouble(tokens[3]);
                d[i][2] = Double.parseDouble(tokens[4]);
                int type = Integer.parseInt(tokens[5]);
                AtomType atomType = forceField.getAtomType(Integer.toString(type));
                if (atomType == null) {
                    logger.warning("Check Atom Type for Atom " + (i + 1) + " in " + molecularAssembly.getFile().getName());
                    return false;
                }
                Atom a = new Atom(i + 1, atomName, atomType, d[i]);
                atomList.add(a);
                // Bond Data
                int numberOfBonds = tokens.length - 6;
                for (int b = 0; b < 8; b++) {
                    if (b < numberOfBonds) {
                        int bond = Integer.parseInt(tokens[6 + b]);
                        bonds[i][b] = bond;
                    } else {
                        bonds[i][b] = 0;
                    }
                }
            }
            // Check if this is an archive.
            if (br.ready()) {
                // Read past blank lines between archive files
                data = br.readLine().trim();
                while (data != null && data.equals("")) {
                    data = br.readLine().trim();
                }
                if (data != null) {
                    tokens = data.split(" +", 2);
                    if (tokens != null && tokens.length > 0) {
                        try {
                            int archiveNumberOfAtoms = Integer.parseInt(tokens[0]);
                            if (archiveNumberOfAtoms == numberOfAtoms) {
                                setType(FileType.ARC);
                            }
                        } catch (Exception e) {
                            tokens = null;
                        }
                    }
                }
            }
            br.close();
            fr.close();
            // Try to renumber
            if (renumber) {
                for (int i = 0; i < numberOfAtoms; i++) {
                    if (labelHash.containsKey(label[i])) {
                        logger.warning("Two atoms have the same index: " + label[i]);
                        return false;
                    }
                    labelHash.put(label[i], i + 1);
                }
                for (int i = 0; i < numberOfAtoms; i++) {
                    int j = -1;
                    while (j < 3 && bonds[i][++j] > 0) {
                        bonds[i][j] = labelHash.get(bonds[i][j]);
                    }
                }
            }
            bondList = new ArrayList<Bond>();
            int c[] = new int[2];
            for (int i = 1; i <= numberOfAtoms; i++) {
                int a1 = i;
                int j = -1;
                while (j < 7 && bonds[i - 1][++j] > 0) {
                    int a2 = bonds[i - 1][j];
                    if (a1 < a2) {
                        if (a1 > numberOfAtoms || a1 < 1 || a2 > numberOfAtoms || a2 < 1) {
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + molecularAssembly.getFile().getName());
                            return false;
                        }
                        // Check for bidirectional connection
                        boolean bidirectional = false;
                        int k = -1;
                        while (k < 7 && bonds[a2 - 1][++k] > 0) {
                            int a3 = bonds[a2 - 1][k];
                            if (a3 == a1) {
                                bidirectional = true;
                                break;
                            }
                        }
                        if (!bidirectional) {
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + molecularAssembly.getFile().getName());
                            return false;
                        }
                        Atom atom1 = atomList.get(a1 - 1);
                        Atom atom2 = atomList.get(a2 - 1);
                        if (atom1 == null || atom2 == null) {
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + molecularAssembly.getFile().getName());
                            return false;
                        }
                        Bond bond = new Bond(atom1, atom2);
                        c[0] = atom1.getAtomType().atomClass;
                        c[1] = atom2.getAtomType().atomClass;
                        String key = BondType.sortKey(c);
                        BondType bondType = forceField.getBondType(key);
                        if (bondType == null) {
                            logger.severe("No BondType for key: " + key);
                        } else {
                            bond.setBondType(bondType);
                        }
                        bondList.add(bond);
                    }
                }
            }
            if (getType() == FileType.ARC) {
                return readtrajectory();
            }
            return true;
        } catch (IOException e) {
            logger.severe(e.toString());
        }
        return false;
    }

    public boolean readtrajectory() {
        // If the first entry was read successfully, reopen the
        // archive file and read the rest of the coordinates
        try {
            logger.info("Trying to parse " + molecularAssembly.getFile() + " as an archive.");
            BufferedReader bin = new BufferedReader(new FileReader(
                    molecularAssembly.getFile()));
            String data = null;
            int numatoms = atomList.size();
            int cycle = 1;
            double coords[][] = new double[numatoms][3];
            // Read past the first N + 1 non-blank lines
            for (int i = 0; i < numatoms + 1; i++) {
                data = bin.readLine();
                while (data != null && data.trim().equals("")) {
                    data = bin.readLine();
                }
            }
            while (bin.ready()) {
                data = bin.readLine();
                // Read past blank lines
                while (data != null && data.trim().equals("")) {
                    data = bin.readLine();
                }
                try {
                    int num = Integer.parseInt(data.trim().split(" +")[0]);
                    if (num != numatoms) {
                        logger.warning(num + " atoms for archive entry " + cycle + " is not equal to " + numatoms + "." + "Only the first " + (cycle - 1) + " entries were read.");
                        return true;
                    }
                } catch (Exception e) {
                    logger.severe(e.toString());
                    return false;
                }
                for (int i = 0; i < numatoms; i++) {
                    data = bin.readLine();
                    // Read past blank lines
                    while (data != null && data.trim().equals("")) {
                        data = bin.readLine();
                    }
                    String[] tokens = data.trim().split(" +");
                    if (tokens == null || tokens.length < 6) {
                        logger.warning("Check atom " + (i + 1) + ", archive entry " + (cycle + 1) + " in " + molecularAssembly.getFile().getName());
                        return false;
                    }
                    coords[i][0] = Double.parseDouble(tokens[2]);
                    coords[i][1] = Double.parseDouble(tokens[3]);
                    coords[i][2] = Double.parseDouble(tokens[4]);
                }
                for (Atom a : atomList) {
                    int i = a.xyzIndex - 1;
                    Vector3d v3d = new Vector3d(coords[i][0], coords[i][1],
                            coords[i][2]);
                    a.addTrajectoryCoords(v3d, cycle);
                }
                cycle++;
            }
            molecularAssembly.setCycles(cycle);
            molecularAssembly.setFileType(FileType.ARC);
            setFileRead(true);
            return true;
        } catch (FileNotFoundException e) {
            logger.warning(e.toString());
        } catch (IOException e) {
            logger.warning(e.toString());
        }
        return false;
    }

    @Override
    public boolean writeFile() {
        File xyzfile = molecularAssembly.getFile();
        if (xyzfile == null) {
            return false;
        }
        try {
            File newFile = version(xyzfile);
            molecularAssembly.setFile(newFile);
            molecularAssembly.setName(newFile.getName());
            FileWriter fw = new FileWriter(newFile);
            BufferedWriter bw = new BufferedWriter(fw);
            // XYZ File First Line
            int numberOfAtoms = molecularAssembly.getAtomList().size();
            String output = String.format("%6d  %s\n", numberOfAtoms,
                    molecularAssembly.toString());
            logger.info("Writing first line of XYZ file:\n" + output);
            bw.write(output);
            Atom a2;
            StringBuffer line;
            StringBuffer lines[] = new StringBuffer[numberOfAtoms];
            // XYZ File Atom Lines
            ArrayList<Atom> atoms = molecularAssembly.getAtomList();
            Vector3d offset = molecularAssembly.getOffset();
            for (Atom a : atoms) {
                line = new StringBuffer(String.format(
                        "%6d  %3s%14.8f%14.8f%14.8f%6d", a.getXYZIndex(), a.getID(), a.getX() - offset.x, a.getY() - offset.y, a.getZ() - offset.z, a.getType()));
                for (Bond b : a.getBonds()) {
                    a2 = b.get1_2(a);
                    line.append(String.format("%6d", a2.xyzIndex));
                }
                lines[a.getXYZIndex() - 1] = line.append("\n");
            }
            try {
                for (int i = 0; i < numberOfAtoms; i++) {
                    bw.write(lines[i].toString());
                }
            } catch (Exception e) {
                logger.severe(
                        "Their was an unexpected error writing to " + getMolecularSystem().toString() + "\n" + e + "\nForce Field Xplor will continue...");
                return false;
            }
            bw.close();
            fw.close();
        } catch (IOException e) {
            logger.severe(
                    "Their was an unexpected error writing to " + getMolecularSystem().toString() + "\n" + e + "\nForce Field Xplor will continue...");
            return false;
        }
        return true;
    }

    public boolean writeP1(Crystal crystal) {
        File xyzfile = molecularAssembly.getFile();
        xyzfile = SystemFilter.version(xyzfile);
        int nSymm = crystal.spaceGroup.symOps.size();
        if (xyzfile == null) {
            return false;
        }
        try {
            FileWriter fw = new FileWriter(xyzfile);
            BufferedWriter bw = new BufferedWriter(fw);
            // XYZ File First Line
            int numberOfAtoms = molecularAssembly.getAtomList().size() * nSymm;
            String output = String.format("%6d  %s\n", numberOfAtoms,
                    molecularAssembly.toString());
            bw.write(output);
            Atom a2;
            StringBuffer line;
            StringBuffer lines[] = new StringBuffer[numberOfAtoms];
            // XYZ File Atom Lines
            ArrayList<Atom> atoms = molecularAssembly.getAtomList();
            Vector3d offset = molecularAssembly.getOffset();
            double xyz[] = new double[3];
            Vector<SymOp> symOps = crystal.spaceGroup.symOps;
            int ii = 0;
            for (SymOp symOp : symOps) {
                int indexOffset = ii * atoms.size();
                ii++;
                for (Atom a : atoms) {
                    int index = a.getXYZIndex() + indexOffset;
                    String id = a.getID();
                    xyz[0] = a.getX() - offset.x;
                    xyz[1] = a.getY() - offset.y;
                    xyz[2] = a.getZ() - offset.z;
                    crystal.applySymOp(xyz, xyz, symOp);
                    int type = a.getType();
                    line = new StringBuffer(String.format(
                            "%6d  %3s%14.8f%14.8f%14.8f%6d", index, id, xyz[0],
                            xyz[1], xyz[2], type));
                    for (Bond b : a.getBonds()) {
                        a2 = b.get1_2(a);
                        line.append(String.format("%7d", a2.xyzIndex + indexOffset));
                    }
                    lines[index - 1] = line.append("\n");
                }
            }
            try {
                for (int i = 0; i < numberOfAtoms; i++) {
                    bw.write(lines[i].toString());
                }
            } catch (Exception e) {
                logger.severe(
                        "Their was an unexpected error writing to " + getMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
                return false;
            }
            bw.close();
            fw.close();
        } catch (IOException e) {
            logger.severe(
                    "Their was an unexpected error writing to " + getMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
            return false;
        }
        return true;
    }
}
