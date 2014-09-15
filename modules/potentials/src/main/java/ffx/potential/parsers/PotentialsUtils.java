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

// FFX imports
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;

// Java imports
import java.io.File;
import java.util.Arrays;
import java.util.logging.Logger;

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
    
    public PotentialsUtils() {
        initTime = System.nanoTime();
        interTime = initTime;
    }
    
    /**
     * Logs time since this interface was created and the last time this method
     * was called.
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
     * @return true
     */
    @Override
    public boolean isLocal() {
        return true;
    }
    
    /**
     * Opens a file and returns all created MolecularAssembly objects.
     * @param file Filename to open
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] open(String file) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.run();
        return opener.getAllAssemblies();
    }
    
    /**
     * Opens an array of files and returns the created MolecularAssembly objects.
     * @param files Filenames to open.
     * @return Array of MolecularAssembly.
     */
    @Override
    public MolecularAssembly[] open(String[] files) {
        PotentialsFileOpener opener = new PotentialsFileOpener(files);
        opener.run();
        return opener.getAllAssemblies();
    }
    
    /**
     * Does nothing (no higher-level data structure).
     */
    @Override
    public void close() {
        // Is not meaningful in the local implementation.
    }
    
    /**
     * Does nothing (no higher-level data structure).
     */
    @Override
    public void closeAll() {
        // Is not meaningful in the local implementation.
    }
    
    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
     * @param assembly MolecularAssembly to save
     * @param file Destination .xyz
     */
    @Override
    public void save(MolecularAssembly assembly, File file) {
        saveAsXYZ(assembly, file);
    }

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
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

    /**
     * Saves the current state of an array of MolecularAssemblys to a PDB file.
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
    
    /**
     * Evaluates the energy of a MolecularAssembly and returns its ForceFieldEnergy
     * object.
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
}

    
    /**
     * Returns a FileOpener thread which can be used to create a MolecularAssembly
     * from a file.
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