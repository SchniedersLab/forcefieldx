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
package ffx.potential.utils;

import java.io.File;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;

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
    
    public MolecularAssembly[] convertDataStructure(Object data);
    
    public MolecularAssembly[] convertDataStructure(Object data, File file);
    
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

    public ForceFieldEnergy energy(MolecularAssembly assembly);

    public double returnEnergy(MolecularAssembly assembly);

    // Subsequent methods were when I was duplicating MainPanel's open() methods,
    // instead of its openWait() methods.
    /*public FileOpener open(String file);
     public FileOpener open(String[] files);
     public FileOpener open(File file, String commandDescription);
     public FileOpener open(File[] files, String commandDescription);*/
}
