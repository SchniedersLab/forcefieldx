/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.potential.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import javax.vecmath.Vector3d;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;

/**
 * The XYZFilter class parses TINKER Cartesian coordinate (*.XYZ) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class XYZFilter extends SystemFilter {

    private static final Logger logger = Logger.getLogger(XYZFilter.class.getName());
    private BufferedReader bin = null;
    private int snapShot;

    /**
     * <p>
     * readOnto</p>
     *
     * @param newFile a {@link java.io.File} object.
     * @param oldSystem a {@link ffx.potential.MolecularAssembly} object.
     * @return a boolean.
     */
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

    /**
     * <p>
     * Constructor for XYZFilter.</p>
     *
     * @param files a {@link java.util.List} object.
     * @param system a {@link ffx.potential.MolecularAssembly} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public XYZFilter(List<File> files, MolecularAssembly system,
            ForceField forceField, CompositeConfiguration properties) {
        super(files, system, forceField, properties);
        this.fileType = FileType.XYZ;
    }

    /**
     * <p>
     * Constructor for XYZFilter.</p>
     *
     * @param file a {@link java.io.File} object.
     * @param system a {@link ffx.potential.MolecularAssembly} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     */
    public XYZFilter(File file, MolecularAssembly system,
            ForceField forceField, CompositeConfiguration properties) {
        super(file, system, forceField, properties);
        this.fileType = FileType.XYZ;
    }

    /**
     * {@inheritDoc}
     *
     * Parse the XYZ File
     */
    @Override
    public boolean readFile() {
        File xyzFile = activeMolecularAssembly.getFile();

        if (forceField == null) {
            logger.warning(" No force field is associated with " + xyzFile.toString());
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
                getActiveMolecularSystem().setName(tokens[1]);
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
                    logger.warning("Check atom " + (i + 1) + " in " + activeMolecularAssembly.getFile().getName());
                    return false;
                }
                tokens = data.trim().split(" +");
                if (tokens == null || tokens.length < 6) {
                    logger.warning("Check atom " + (i + 1) + " in " + activeMolecularAssembly.getFile().getName());
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
                    StringBuilder message = new StringBuilder("Check atom type ");
                    message.append(type).append(" for Atom ").append(i + 1);
                    message.append(" in ").append(activeMolecularAssembly.getFile().getName());
                    logger.warning(message.toString());
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
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + activeMolecularAssembly.getFile().getName());
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
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + activeMolecularAssembly.getFile().getName());
                            return false;
                        }
                        Atom atom1 = atomList.get(a1 - 1);
                        Atom atom2 = atomList.get(a2 - 1);
                        if (atom1 == null || atom2 == null) {
                            logger.warning("Check the Bond Bewteen " + a1 + " and " + a2 + " in " + activeMolecularAssembly.getFile().getName());
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
            /*
             if (getType() == FileType.ARC) {
             return readtrajectory();
             } */
            return true;
        } catch (IOException e) {
            logger.severe(e.toString());
        }
        return false;
    }

    /**
     * Reads the next snap-shot of an archive into the activeMolecularAssembly.
     * After calling this function, a BufferedReader will remain open until the
     * <code>close</code> method is called.
     *
     * @return true if successful.
     */
    public boolean readNext() {
        try {
            String data = null;
            Atom atoms[] = activeMolecularAssembly.getAtomArray();
            int nSystem = atoms.length;

            if (bin == null || !bin.ready()) {
                bin = new BufferedReader(new FileReader(currentFile));
                // Read past the first N + 1 non-blank lines
                for (int i = 0; i < nSystem + 1; i++) {
                    data = bin.readLine();
                    while (data != null && data.trim().equals("")) {
                        data = bin.readLine();
                    }
                }
                snapShot = 1;
            }

            snapShot++;
            logger.info(String.format(" Reading snapshot %d of %s.",
                    snapShot, activeMolecularAssembly));

            data = bin.readLine();
            // Read past blank lines
            while (data != null && data.trim().equals("")) {
                data = bin.readLine();
            }
            try {
                int nArchive = Integer.parseInt(data.trim().split(" +")[0]);
                if (nArchive != nSystem) {
                    String message = String.format("Number of atoms mismatch (Archive: %d, System: %d).", nArchive, nSystem);
                    logger.warning(message);
                    return false;
                }
            } catch (Exception e) {
                logger.severe(e.toString());
                return false;
            }
            for (int i = 0; i < nSystem; i++) {
                data = bin.readLine();
                // Read past blank lines
                while (data != null && data.trim().equals("")) {
                    data = bin.readLine();
                }
                String[] tokens = data.trim().split(" +");
                if (tokens == null || tokens.length < 6) {
                    String message = String.format("Check atom %d in %s.", (i + 1),
                            currentFile.getName());
                    logger.warning(message);
                    return false;
                }
                double x = Double.parseDouble(tokens[2]);
                double y = Double.parseDouble(tokens[3]);
                double z = Double.parseDouble(tokens[4]);
                int xyzIndex = atoms[i].getXYZIndex();
                if (xyzIndex != i + 1) {
                    String message = String.format("Archive atom index %d being read onto system atom index %d.", i + 1, xyzIndex);
                    logger.warning(message);
                }
                atoms[i].moveTo(x, y, z);
            }
        } catch (FileNotFoundException e) {
            String message = String.format("Exception opening file %s.", currentFile);
            logger.log(Level.WARNING, message, e);
        } catch (IOException e) {
            String message = String.format("Exception reading from file %s.",
                    currentFile);
            logger.log(Level.WARNING, message, e);
        }
        return false;
    }

    /**
     * <p>
     * close</p>
     */
    public void close() {
        if (bin != null) {
            try {
                bin.close();
            } catch (Exception e) {
                String message = String.format("Exception closing file %s.",
                        activeMolecularAssembly.getFile());
                logger.log(Level.WARNING, message, e);
            }
        }
    }

    /**
     * <p>
     * readtrajectory</p>
     *
     * @return a boolean.
     */
    public boolean readtrajectory() {
        // If the first entry was read successfully, reopen the
        // archive file and read the rest of the coordinates
        try {
            logger.info(" Trying to parse " + activeMolecularAssembly.getFile() + " as an archive.");
            bin = new BufferedReader(new FileReader(
                    activeMolecularAssembly.getFile()));
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
                        logger.warning("Check atom " + (i + 1) + ", archive entry " + (cycle + 1) + " in " + activeMolecularAssembly.getFile().getName());
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
            activeMolecularAssembly.setCycles(cycle);
            setFileRead(true);
            return true;
        } catch (FileNotFoundException e) {
            logger.warning(e.toString());
        } catch (IOException e) {
            logger.warning(e.toString());
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        if (saveFile == null) {
            return false;
        }
        try {
            File newFile = saveFile;
            if (!append) {
                newFile = version(saveFile);
            }
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            FileWriter fw = null;
            if (append && !newFile.exists()) {
                fw = new FileWriter(newFile);
            } else {
                fw = new FileWriter(newFile, append);
            }
            BufferedWriter bw = new BufferedWriter(fw);

            // XYZ File First Line
            int numberOfAtoms = activeMolecularAssembly.getAtomList().size();
            String output = format("%7d  %s\n", numberOfAtoms, activeMolecularAssembly.toString());
            bw.write(output);
            Atom a2;
            StringBuilder line;
            StringBuilder lines[] = new StringBuilder[numberOfAtoms];
            // XYZ File Atom Lines
            ArrayList<Atom> atoms = activeMolecularAssembly.getAtomList();
            Vector3d offset = activeMolecularAssembly.getOffset();
            for (Atom a : atoms) {
                line = new StringBuilder(format(
                        "%7d %3s%14.8f%14.8f%14.8f%6d", a.getXYZIndex(), a.getAtomType().name, a.getX() - offset.x, a.getY() - offset.y, a.getZ() - offset.z, a.getType()));
                for (Bond b : a.getBonds()) {
                    a2 = b.get1_2(a);
                    line.append(format("%8d", a2.xyzIndex));
                }
                lines[a.getXYZIndex() - 1] = line.append("\n");
            }
            try {
                for (int i = 0; i < numberOfAtoms; i++) {
                    bw.write(lines[i].toString());
                }
            } catch (Exception e) {
                logger.severe(
                        "Their was an unexpected error writing to " + getActiveMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
                return false;
            }
            bw.close();
            fw.close();
        } catch (IOException e) {
            logger.severe(
                    "Their was an unexpected error writing to " + getActiveMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
            return false;
        }
        return true;
    }

    /**
     * <p>
     * writeFileAsP1</p>
     *
     * @param saveFile a {@link java.io.File} object.
     * @param append a boolean.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @return a boolean.
     */
    public boolean writeFileAsP1(File saveFile, boolean append, Crystal crystal) {
        if (saveFile == null) {
            return false;
        }
        try {
            File newFile = saveFile;
            if (!append) {
                newFile = version(saveFile);
            }
            activeMolecularAssembly.setFile(newFile);
            activeMolecularAssembly.setName(newFile.getName());
            FileWriter fw = new FileWriter(newFile, append);
            BufferedWriter bw = new BufferedWriter(fw);
            int nSymm = crystal.spaceGroup.symOps.size();
            // XYZ File First Line
            int numberOfAtoms = activeMolecularAssembly.getAtomList().size() * nSymm;
            String output = format("%7d %s\n", numberOfAtoms, activeMolecularAssembly.toString());
            bw.write(output);
            Atom a2;
            StringBuilder line;
            StringBuilder lines[] = new StringBuilder[numberOfAtoms];
            // XYZ File Atom Lines
            Atom atoms[] = activeMolecularAssembly.getAtomArray();
            double xyz[] = new double[3];
            for (int iSym = 0; iSym < nSymm; iSym++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSym);
                int indexOffset = iSym * atoms.length;
                for (Atom a : atoms) {
                    int index = a.getXYZIndex() + indexOffset;
                    String id = a.getAtomType().name;
                    xyz[0] = a.getX();
                    xyz[1] = a.getY();
                    xyz[2] = a.getZ();
                    crystal.applySymOp(xyz, xyz, symOp);
                    int type = a.getType();
                    line = new StringBuilder(format(
                            "%7d %3s%14.8f%14.8f%14.8f%6d", index, id, xyz[0],
                            xyz[1], xyz[2], type));
                    for (Bond b : a.getBonds()) {
                        a2 = b.get1_2(a);
                        line.append(format("%8d", a2.xyzIndex + indexOffset));
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
                        "Their was an unexpected error writing to " + getActiveMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
                return false;
            }
            bw.close();
            fw.close();
        } catch (IOException e) {
            logger.severe(
                    "Their was an unexpected error writing to " + getActiveMolecularSystem().toString() + "\n" + e + "\nForce Field X will continue...");
            return false;
        }
        return true;
    }
}
