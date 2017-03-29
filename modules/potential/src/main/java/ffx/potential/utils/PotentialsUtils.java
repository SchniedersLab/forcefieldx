/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.potential.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.Mutation;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;

/**
 * The PotentialsUtils class provides a local implementation, independent of the
 * User Interfaces module, of PotentialsFunctions methods such as file opening.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class PotentialsUtils implements PotentialsFunctions {

    private static final Logger logger = Logger.getLogger(PotentialsUtils.class.getName());
    private final long initTime;
    private long interTime;
    private SystemFilter lastFilter;

    public PotentialsUtils() {
        initTime = System.nanoTime();
        interTime = initTime;
    }

    /**
     * Logs time since this interface was created and the last time this method
     * was called.
     *
     * @return Time since last call (double).
     */
    @Override
    public double time() {
        long currTime = System.nanoTime();
        logger.info(String.format(" Time since interface established: %f", (currTime - initTime) * 1.0E-9));
        double elapsed = (currTime - interTime) * 1.0E-9;
        interTime = currTime;
        logger.info(String.format(" Time since last timer call: %f", elapsed));
        return elapsed;
    }

    /**
     * Returns true (this is the local implementation).
     *
     * @return true
     */
    @Override
    public boolean isLocal() {
        return true;
    }

    /**
     * Opens a file and returns all created MolecularAssembly objects.
     *
     * @param file Filename to open
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] openAll(String file) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }
    
    @Override
    public MolecularAssembly open(String filename) {
        PotentialsFileOpener opener = new PotentialsFileOpener(filename);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", filename);
        }
        return opener.getAssembly();
    }
    
    /**
     * One one file object.
     */
    public MolecularAssembly open(File file) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", file.getName());
        }
        return opener.getAssembly();
    }
    
    /**
     * Mutates file on-the-fly as it is being opened.
     * Used to open files for pHMD in fully-protonated form.
     */
    public MolecularAssembly openWithMutations(File file, List<Mutation> mutations) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.setMutations(mutations);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", file.getName());
        }
        return opener.getAssembly();
    }
    
    /**
     * Open one filename string without printing all the header material.
     */
    public MolecularAssembly openQuietly(String filename) {
        Logger ffxLog = Logger.getLogger("ffx");
        Level prevLev = ffxLog.getLevel();
        ffxLog.setLevel(Level.WARNING);
        MolecularAssembly mola = open(filename);
        ffxLog.setLevel(prevLev);
        return mola;
    }
    
    /**
     * Open one File object without printing all the header material.
     */
    public MolecularAssembly openQuietly(File file) {
        Logger ffxLog = Logger.getLogger("ffx");
        Level prevLev = ffxLog.getLevel();
        ffxLog.setLevel(Level.WARNING);
        MolecularAssembly mola = open(file);
        ffxLog.setLevel(prevLev);
        return mola;
    }

    /**
     * Opens an array of files and returns the created MolecularAssembly
     * objects.
     *
     * @param files Filenames to open.
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] openAll(String[] files) {
        PotentialsFileOpener opener = new PotentialsFileOpener(files);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * Opens a file and returns all created MolecularAssembly objects, setting
     * any underlying Potential to use a certain number of threads.
     *
     * @param file Filename to open
     * @param nThreads Use non-default num threads
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] openAll(String file, int nThreads) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.setNThreads(nThreads);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * Opens an array of files and returns all created MolecularAssembly
     * objects, setting any underlying Potential to use a certain number of
     * threads.
     *
     * @param files Filenames to open.
     * @param nThreads Use non-default num threads
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] open(String[] files, int nThreads) {
        PotentialsFileOpener opener = new PotentialsFileOpener(files);
        opener.setNThreads(nThreads);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * Converts a data structure (such as a Biojava Structure) into one or more
     * MolecularAssembly objects.
     *
     * @param data Structure to convert
     * @return Array of MolecularAssembly
     */
    @Override
    public MolecularAssembly[] convertDataStructure(Object data) {
        try {
            PotentialsDataConverter converter = new PotentialsDataConverter(data);
            converter.run();
            return converter.getAllAssemblies();
        } catch (FileNotFoundException | IllegalArgumentException ex) {
            logger.warning(String.format(" Exception in data structure conversion: %s", ex.toString()));
            return null;
        }
    }

    /**
     * Converts a data structure (such as a Biojava Structure) into one or more
     * MolecularAssembly objects.
     *
     * @param data Structure to convert
     * @param file Source file
     * @return Array of MolecularAssembly
     */
    @Override
    public MolecularAssembly[] convertDataStructure(Object data, File file) {
        try {
            PotentialsDataConverter converter = new PotentialsDataConverter(data, file);
            converter.run();
            return converter.getAllAssemblies();
        } catch (FileNotFoundException | IllegalArgumentException ex) {
            logger.warning(String.format(" Exception in data structure conversion: %s", ex.toString()));
            return null;
        }
    }

    /**
     * Converts a data structure (such as a Biojava Structure) into one or more
     * MolecularAssembly objects.
     *
     * @param data Structure to convert
     * @param filename Source file
     * @return Array of MolecularAssembly
     */
    @Override
    public MolecularAssembly[] convertDataStructure(Object data, String filename) {
        File file = new File(filename);
        if (!file.exists() || file.isDirectory() || !file.canRead()) {
            logger.warning(String.format("%s not a valid file name: file name discarded.", filename));
            return convertDataStructure(data);
        }
        return convertDataStructure(data, file);
    }

    /**
     * Shuts down parallel teams in the force field of the provided
     * MolecularAssembly. Kaminsky's ParallelTeamThreads' run() methods are
     * infinite loops, and because running threads are always GC roots, it is
     * necessary to send them a signal to shut down to enable garbage
     * collection.
     *
     * @param assembly Assembly to close.
     */
    @Override
    public void close(MolecularAssembly assembly) {
        assembly.destroy();
    }

    /**
     * Shuts down parallel teams in the force fields of the provided
     * MolecularAssemblys.
     *
     * @param assemblies Assemblies to close.
     */
    @Override
    public void closeAll(MolecularAssembly[] assemblies) {
        for (MolecularAssembly assembly : assemblies) {
            assembly.destroy();
        }
    }

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
     *
     * @param assembly MolecularAssembly to save
     * @param file Destination .xyz
     */
    @Override
    public void save(MolecularAssembly assembly, File file) {
        saveAsXYZ(assembly, file);
    }

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
     *
     * @param assembly MolecularAssembly to save
     * @param file Destination .xyz
     */
    @Override
    public void saveAsXYZ(MolecularAssembly assembly, File file) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            XYZFilter xyzFilter = new XYZFilter(file, assembly, null, null);
            if (!xyzFilter.writeFile(file, false)) {
                logger.info(String.format(" Save failed for %s", assembly.toString()));
            }
        }
    }

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file as a P1
     * crystal.
     *
     * @param assembly MolecularAssembly to save
     * @param file Destination .xyz
     */
    @Override
    public void saveAsP1(MolecularAssembly assembly, File file) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            XYZFilter filter = new XYZFilter(file, assembly, null, null);
            ForceField forceField = assembly.getForceField();
            final double a = forceField.getDouble(ForceFieldDouble.A_AXIS, 10.0);
            final double b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
            final double c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
            final double alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
            final double beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
            final double gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
            final String spacegroup = forceField.getString(
                    ForceFieldString.SPACEGROUP, "P1");
            Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
            if (!filter.writeFileAsP1(file, false, crystal)) {
                logger.info(String.format(" Save failed for %s", assembly.toString()));
            }
        }
    }

    /**
     * Saves the current state of a MolecularAssembly to a PDB file.
     *
     * @param assembly MolecularAssembly to save
     * @param file Destination .pdb
     */
    @Override
    public void saveAsPDB(MolecularAssembly assembly, File file) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, assembly, null, null);
            if (!pdbFilter.writeFile(file, false)) {
                logger.info(String.format(" Save failed for %s", assembly.toString()));
            }
        }
    }

    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file) {
        savePDBSymMates(assembly, file, "_symMate");
    }

    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file, String suffix) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, assembly, null, null);
            if (!pdbFilter.writeFile(file, false)) {
                logger.info(String.format(" Save failed for %s", assembly.toString()));
            } else {
                Crystal crystal = assembly.getCrystal();
                int nSymOps = crystal.spaceGroup.getNumberOfSymOps();
                String filename = FilenameUtils.removeExtension(file.getName());
                for (int i = 1; i < nSymOps; i++) {
                    pdbFilter.setSymOp(i);
                    String saveFileName = filename + suffix + "_" + i + ".pdb";
                    File saveFile = new File(saveFileName);
                    for (int j = 1; j < 1000; j++) {
                        if (!saveFile.exists()) {
                            break;
                        }
                        saveFile = new File(saveFileName + "_" + j);
                    }
                    StringBuilder symSb = new StringBuilder();
                    String[] symopLines = crystal.spaceGroup.getSymOp(i).toString().split("\\r?\\n");
                    int nLines = symopLines.length;
                    symSb.append("REMARK 350\nREMARK 350 SYMMETRY OPERATORS");
                    for (int j = 0; j < nLines; j++) {
                        symSb.append("\nREMARK 350 ").append(symopLines[j]);
                    }

                    symopLines = crystal.spaceGroup.getSymOp(i).toXYZString().split("\\r?\\n");
                    nLines = symopLines.length;
                    symSb.append("\nREMARK 350\nREMARK 350 SYMMETRY OPERATORS XYZ FORM");
                    for (int j = 0; j < nLines; j++) {
                        symSb.append("\nREMARK 350 ").append(symopLines[j]);
                    }

                    if (saveFile.exists()) {
                        logger.warning(String.format(" Could not successfully version file "
                                + "%s: appending to file %s", saveFileName, saveFile.getName()));
                        if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), true)) {
                            logger.info(String.format(" Save failed for %s", saveFile.getName()));
                        }
                    } else if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), false)) {
                        logger.info(String.format(" Save failed for %s", saveFile.getName()));
                    }
                }
            }
        }
    }

    public void saveAsSIFTPDB(MolecularAssembly assembly, File file, String[] resAndScore) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else if (resAndScore == null) {
            logger.info(" Res and score array was null.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, assembly, null, null);
            if (!pdbFilter.writeSIFTFile(file, false, resAndScore)) {
                logger.info(String.format(" Save failed for %s", assembly.toString()));
            }
        }
    }

    /**
     * Saves the current state of an array of MolecularAssemblys to a PDB file.
     *
     * @param assemblies MolecularAssembly array to save
     * @param file Destination .pdb
     */
    @Override
    public void saveAsPDB(MolecularAssembly[] assemblies, File file) {
        if (assemblies == null) {
            logger.info(" Null array of molecular assemblies to write.");
        } else if (assemblies.length == 0) {
            logger.info(" Zero-length array of molecular assemblies to write.");
        } else if (file == null) {
            logger.info(" No valid file to write to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file,
                    Arrays.asList(assemblies), null, null);
            pdbFilter.writeFile(file, false);
        }
    }

    public void saveAsSIFTPDB(MolecularAssembly[] assemblies, File file, String[] resAndScore) {
        if (assemblies == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else if (resAndScore == null) {
            logger.info(" Res and score array was null.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, Arrays.asList(assemblies), null, null);
            pdbFilter.writeSIFTFile(file, false, resAndScore);
        }
    }

    /**
     * Evaluates the energy of a MolecularAssembly and returns its
     * ForceFieldEnergy object.
     *
     * @param assembly To evaluate
     * @return assembly's ForceFieldEnergy.
     */
    @Override
    public ForceFieldEnergy energy(MolecularAssembly assembly) {
        if (assembly == null) {
            logger.info(" Molecular assembly was null - skipping energy");
            return null;
        } else {
            ForceFieldEnergy energy = assembly.getPotentialEnergy();
            if (energy == null) {
                energy = new ForceFieldEnergy(assembly);
                assembly.setPotential(energy);
            }
            energy.energy(false, true);
            return energy;
        }
    }

    /**
     * Returns the energy of a MolecularAssembly in kcal/mol (as a double) and
     * prints the energy evaluation
     *
     * @param assembly To evaluate energy of
     * @return Potential energy (kcal/mol)
     */
    @Override
    public double returnEnergy(MolecularAssembly assembly) {
        if (assembly == null) {
            logger.info(" Molecular assembly was null - skipping energy");
            return 0.0;
        } else {
            ForceFieldEnergy energy = assembly.getPotentialEnergy();
            if (energy == null) {
                energy = new ForceFieldEnergy(assembly);
                assembly.setPotential(energy);
            }
            return energy.energy(false, true);
        }
    }

    @Override
    public SystemFilter getFilter() {
        return lastFilter;
    }

    public static void analysis(MolecularAssembly molas[]) {
        for (MolecularAssembly mola : molas) {
            analysis(mola);
        }
    }

    public static void analysis(MolecularAssembly mola) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Atom Array: (array_pos, xyz_index, resName, atomName, typeNum, classNum) \n"));
        Atom atoms[] = mola.getAtomArray();
        for (int i = 0; i < atoms.length; i++) {
            String resName = atoms[i].getResidueName();
            String atomName = atoms[i].describe(Atom.Descriptions.TRIM);
            AtomType atomType = atoms[i].getAtomType();
            int typeNum = atomType.type;
            int classNum = atomType.atomClass;
            int xyzIndex = atoms[i].getIndex();
            sb.append(String.format("   %d: %d %s %s %d %d\n", i, xyzIndex, resName, atomName, typeNum, classNum));
        }
        logger.info(sb.toString());
    }

}

/**
 * Returns a FileOpener thread which can be used to create a MolecularAssembly
 * from a file.
 *
 * @param file To be opened
 * @return Opener thread.
 */
/*@Override
 public FileOpener open(String file) {
 return new PotentialsFileOpener(file);
 }*/
/**
 * Returns an array of FileOpener threads which can be used to create
 * MolecularAssembly objects from an array of files.
 *
 * @param filenames To be opened
 * @return Opener threads.
 */
/*@Override
 public FileOpener open(String[] filenames) {
 int numFiles = filenames.length;
 File[] files = new File[numFiles];
 for (int i = 0; i < numFiles; i++) {
 files[i] = new File(filenames[i]);
 }
 return new PotentialsFileOpener(files);
 }

 @Override
 public FileOpener open(File file, String commandDescription) {
 return new PotentialsFileOpener(file);
 }

 @Override
 public FileOpener open(File[] files, String commandDescription) {
 return new PotentialsFileOpener(files);
 }*/
