/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.parsers.SystemFilter;

/**
 * The PotentialsFunctions interface specifies utility methods such as opening
 * files into MolecularAssemblys, evaluating energy, and saving assemblies to
 * files. Intended to be analogous to existing Groovy method closures, with both
 * local implementation and a User Interfaces implementation which interacts
 * with our GUI and underlying data structure. This should enable other users to
 * import only Potentials and its dependencies, and slide in their own UI and
 * data structure on top of Potentials.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface PotentialsFunctions {

    public boolean isLocal(); // Return true if the local implementation from Potentials.

    public MolecularAssembly[] open(String file);

    public MolecularAssembly[] open(String[] files);
    
    default public MolecularAssembly[] open(String file, int nThreads) {
        return open(file);
    }
    
    default public MolecularAssembly[] open(String[] file, int nThreads) {
        return open(file);
    }
    
    public MolecularAssembly[] convertDataStructure(Object data);
    
    public MolecularAssembly[] convertDataStructure(Object data, File file);
    
    public MolecularAssembly[] convertDataStructure(Object data, String filename);
    
    /*public MolecularAssembly[] convertDataStructure(Object[] data);
    
    public MolecularAssembly[] convertDataStructure(Object[] data, File file);*/

    public void close(MolecularAssembly assembly);

    public void closeAll(MolecularAssembly[] assemblies);

    public double time();

    public void save(MolecularAssembly assembly, File file);

    public void saveAsXYZ(MolecularAssembly assembly, File file);

    public void saveAsP1(MolecularAssembly assembly, File file);

    public void saveAsPDB(MolecularAssembly assembly, File file);

    public void saveAsPDB(MolecularAssembly[] assemblies, File file);
    
    // public void saveXYZSymMates(MolecularAssembly assembly, File file);
    
    public void savePDBSymMates(MolecularAssembly assembly, File file);
    // Will use default suffix of _symMate
    
    public void savePDBSymMates(MolecularAssembly assembly, File file, String suffix);

    public ForceFieldEnergy energy(MolecularAssembly assembly);

    public double returnEnergy(MolecularAssembly assembly);
    
    /**
     * Returns the last SystemFilter created by this (may be null).
     * @return 
     */
    public SystemFilter getFilter();

    // Subsequent methods were when I was duplicating MainPanel's open() methods,
    // instead of its openWait() methods.
    /*public FileOpener open(String file);
     public FileOpener open(String[] files);
     public FileOpener open(File file, String commandDescription);
     public FileOpener open(File[] files, String commandDescription);*/
    
    /**
     * Versions a file, attempting to find an unused filename in the set filename,
     * and filename_1 to filename_999.
     * 
     * @param filename To version
     * @return Versioned filename.
     */
    default public String versionFile(String filename) {
        if (filename == null) {
            throw new IllegalArgumentException("Filename must not be null!");
        }
        return versionFile(new File(filename)).getName();
    }
    
    /**
     * Versions a file, attempting to find an unused filename in the set filename,
     * and filename_1 to filename_999.
     * 
     * @param file To version
     * @return Versioned file
     */
    default public File versionFile(File file) {
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
