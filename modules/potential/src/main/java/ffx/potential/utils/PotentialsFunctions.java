//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.utils;

import java.io.File;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.parsers.SystemFilter;

/**
 * <p>
 * PotentialsFunctions describes core functionality for many Force Field X algorithms and
 * scripts, such as opening and closing structure files, basic force field energy evaluations,
 * etc.
 *
 * <p>
 * This is implemented in two locations: UIUtils in the User Interfaces package, and in
 * PotentialsUtils in the Potentials package.
 *
 * <p>
 * The UIUtils implementation is the default for Force Field X; on top of the core
 * functionality, it also updates the FFX graphical user interface and tree structure.
 *
 * <p>
 * The PotentialsUtils implementation lacks the extra functionality of the UIUtils
 * implementation, and simply accomplishes the required task. This is used by our tests, and is
 * also potentially useful for third parties who would like to use FFX without its GUI.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface PotentialsFunctions {
    /**
     * Constant <code>logger</code>
     */
    Logger logger = Logger.getLogger(PotentialsFunctions.class.getName());

    /**
     * Performs any necessary shutdown operations on a MolecularAssembly,
     * primarily shutting down Parallel Java threads and closing any other open
     * resources.
     *
     * @param assembly Assembly to close.
     */
    void close(MolecularAssembly assembly);

    /**
     * Performs any necessary shutdown operations on an array of
     * MolecularAssembly, primarily shutting down Parallel Java threads and
     * closing any other open resources.
     *
     * @param assemblies Assemblies to close.
     */
    void closeAll(MolecularAssembly[] assemblies);

    /**
     * Evaluates the energy of a MolecularAssembly and returns its
     * ForceFieldEnergy object.
     *
     * @param assembly To evaluate
     * @return assembly's ForceFieldEnergy.
     */
    ForceFieldEnergy energy(MolecularAssembly assembly);

    /**
     * Returns either the active assembly from the overlying UI, or the "active"
     * molecular assembly from the last used SystemFilter.
     *
     * @return A MolecularAssembly or null
     */
    default MolecularAssembly getActiveAssembly() {
        SystemFilter filt = getFilter();
        if (filt != null) {
            return filt.getActiveMolecularSystem();
        } else {
            return null;
        }
    }

    /**
     * If available, returns CLI arguments; default implementation does not have
     * access to CLI arguments, and throws UnsupportedOperationException.
     *
     * @return CLI arguments
     * @throws java.lang.UnsupportedOperationException If unimplemented
     * @throws java.lang.UnsupportedOperationException if any.
     */
    default List<String> getArguments() throws UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the last SystemFilter created by this (may be null).
     *
     * @return Last SystemFilter
     */
    SystemFilter getFilter();

    /**
     * True if using a local implementation (not in a user interfaces module).
     *
     * @return If a local implementation
     */
    boolean isLocal(); // Return true if the local implementation from Potentials.

    /**
     * <p>open.</p>
     *
     * @param filename a {@link java.lang.String} object.
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    default MolecularAssembly open(String filename) {
        MolecularAssembly[] assemblies = openAll(filename);
        if (assemblies.length > 1) {
            logger.log(Level.WARNING, " Found multiple assemblies in file {0}, opening first.", filename);
        }
        return assemblies[0];
    }

    /**
     * Opens an array of files and returns all created MolecularAssembly
     * objects, setting any underlying Potential to use a certain number of
     * threads. Default implementation simply ignores nThreads.
     *
     * @param files    an array of {@link java.lang.String} objects.
     * @param nThreads Use non-default num threads
     * @return Array of MolecularAssembly.
     */
    default MolecularAssembly[] open(String[] files, int nThreads) {
        return openAll(files);
    }

    /**
     * Opens a file and returns all created MolecularAssembly objects.
     *
     * @param file Filename to open
     * @return Array of MolecularAssembly.
     */
    MolecularAssembly[] openAll(String file);

    /**
     * Opens an array of files and returns the created MolecularAssembly
     * objects.
     *
     * @param files Filenames to open.
     * @return Array of MolecularAssembly.
     */
    MolecularAssembly[] openAll(String[] files);

    /**
     * Opens a file and returns all created MolecularAssembly objects, setting
     * any underlying Potential to use a certain number of threads. Default
     * implementation simply ignores nThreads.
     *
     * @param file     Filename to open
     * @param nThreads Use non-default num threads
     * @return Array of MolecularAssembly.
     */
    default MolecularAssembly[] openAll(String file, int nThreads) {
        return openAll(file);
    }

    /**
     * Returns the energy of a MolecularAssembly in kcal/mol (as a double) and
     * prints the energy evaluation
     *
     * @param assembly To evaluate energy of
     * @return Potential energy (kcal/mol)
     */
    double returnEnergy(MolecularAssembly assembly);

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
     *
     * @param assembly MolecularAssembly to save
     * @param file     Destination .xyz
     */
    void save(MolecularAssembly assembly, File file);

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file as a P1
     * crystal.
     *
     * @param assembly MolecularAssembly to save
     * @param file     Destination .xyz
     */
    void saveAsP1(MolecularAssembly assembly, File file);

    /**
     * Saves the current state of a MolecularAssembly to a PDB file.
     *
     * @param assembly MolecularAssembly to save
     * @param file     Destination .pdb
     */
    void saveAsPDB(MolecularAssembly assembly, File file);

    /**
     * Saves the current state of an array of MolecularAssemblys to a PDB file.
     *
     * @param assemblies MolecularAssembly array to save
     * @param file       Destination .pdb
     */
    void saveAsPDB(MolecularAssembly[] assemblies, File file);

    void saveAsPDB(MolecularAssembly assembly, File file, boolean writeEnd, boolean append);

    /**
     * Saves the current state of a MolecularAssembly to an XYZ file.
     *
     * @param assembly MolecularAssembly to save
     * @param file     Destination .xyz
     */
    void saveAsXYZ(MolecularAssembly assembly, File file);

    /**
     * Saves the symmetry mates of a MolecularAssembly to PDB files.
     *
     * @param assembly To save
     * @param file     Destination file
     */
    void savePDBSymMates(MolecularAssembly assembly, File file);

    /**
     * Saves the symmetry mates of a MolecularAssembly to PDB files.
     *
     * @param assembly To save
     * @param file     Destination file
     * @param suffix   Custom file suffix
     */
    void savePDBSymMates(MolecularAssembly assembly, File file, String suffix);

    /**
     * Logs time elapsed since last call.
     *
     * @return Time.
     */
    double time();

    /**
     * Versions a file, attempting to find an unused filename in the set
     * filename, and filename_1 to filename_999.
     *
     * @param filename To version
     * @return Versioned filename.
     */
    default String versionFile(String filename) {
        if (filename == null) {
            throw new IllegalArgumentException("Filename must not be null!");
        }
        return versionFile(new File(filename)).getName();
    }

    /**
     * Versions a file, attempting to find an unused filename in the set
     * filename, and filename_1 to filename_999.
     *
     * @param file To version
     * @return Versioned file
     */
    default File versionFile(File file) {
        if (file == null) {
            throw new IllegalArgumentException("File must not be null!");
        }
        int counter = 1;
        String filename = file.getName();
        while (file.exists() && counter < 1000) {
            file = new File(String.format("%s_%d", filename, counter++));
        }
        if (file.exists()) {
            throw new IllegalArgumentException(String.format("Could not version file %s", filename));
        }
        return file;
    }
}
