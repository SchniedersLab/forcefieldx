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
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ARCFileFilter;
import ffx.potential.parsers.FileOpener;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.INTFileFilter;
import ffx.potential.parsers.INTFilter;
import ffx.potential.parsers.PDBFileFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFileFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.Keyword;

/**
 * The PotentialsFileOpener class specifies a Runnable object which is
 * constructed with a File and, when run, allows returning any opened
 * MolecularAssembly objects and their associated properties.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class PotentialsFileOpener implements FileOpener {

    private final File file;
    private final Path filepath;
    private final File[] allFiles;
    private final Path[] allPaths;
    private List<MolecularAssembly> assemblies;
    private MolecularAssembly activeAssembly; // Presently, will just be the first element of assemblies.
    private List<CompositeConfiguration> propertyList;
    private CompositeConfiguration activeProperties;

    public PotentialsFileOpener(File file) {
        if (!file.exists() || !file.isFile()) {
            throw new IllegalArgumentException(String.format(" File %s either did not exist or was not a file.", file.getName()));
        }
        this.file = file;
        Path pwdPath;
        Path absPath;
        try {
            pwdPath = Paths.get(new File("").getCanonicalPath());
        } catch (IOException ex) {
            pwdPath = Paths.get(new File("").getAbsolutePath());
        }
        try {
            absPath = Paths.get(file.getCanonicalPath());
        } catch (IOException ex) {
            absPath = Paths.get(file.getAbsolutePath());
        }
        filepath = pwdPath.relativize(absPath);
        allFiles = new File[1];
        allFiles[0] = this.file;
        allPaths = new Path[1];
        allPaths[0] = this.filepath;
        assemblies = new ArrayList<>();
        propertyList = new ArrayList<>();
    }

    public PotentialsFileOpener(String filename) {
        this(new File(filename));
    }

    public PotentialsFileOpener(Path filepath) {
        this(filepath.toString());
    }

    public PotentialsFileOpener(File[] files) {
        if (files == null) {
            throw new IllegalArgumentException(" Array of files to be opened was null.");
        }
        int numFiles = files.length;
        if (numFiles == 0) {
            throw new IllegalArgumentException(" Array of files to be opened was empty.");
        }
        List<File> fileList = new ArrayList<>();
        List<Path> pathList = new ArrayList<>();
        Path pwdPath;
        try {
            pwdPath = Paths.get(new File("").getCanonicalPath());
        } catch (IOException ex) {
            pwdPath = Paths.get(new File("").getAbsolutePath());
        }
        for (File tryFile : files) {
            if (!(tryFile.exists() && tryFile.isFile())) {
                continue;
            }
            Path absPath;
            try {
                absPath = Paths.get(tryFile.getCanonicalPath());
            } catch (IOException ex) {
                absPath = Paths.get(tryFile.getAbsolutePath());
            }
            Path thisPath = pwdPath.relativize(absPath);
            fileList.add(tryFile);
            pathList.add(thisPath);
        }
        int numAccepted = fileList.size();
        if (numAccepted < 1) {
            throw new IllegalArgumentException(" No valid files could be found to open.");
        }
        allFiles = fileList.toArray(new File[numAccepted]);
        allPaths = pathList.toArray(new Path[numAccepted]);
        this.file = allFiles[0];
        this.filepath = allPaths[0];
        assemblies = new ArrayList<>();
        propertyList = new ArrayList<>();
    }

    public PotentialsFileOpener(String[] filenames) {
        if (filenames == null) {
            throw new IllegalArgumentException(" Array of files to be opened was null.");
        }
        int numFiles = filenames.length;
        if (numFiles == 0) {
            throw new IllegalArgumentException(" Array of files to be opened was empty.");
        }

        List<File> fileList = new ArrayList<>();
        List<Path> pathList = new ArrayList<>();
        Path pwdPath;
        try {
            pwdPath = Paths.get(new File("").getCanonicalPath());
        } catch (IOException ex) {
            pwdPath = Paths.get(new File("").getAbsolutePath());
        }
        for (String filename : filenames) {
            try {
                File tryFile = new File(filename);
                if (!(tryFile.exists() && tryFile.isFile())) {
                    continue;
                }
                Path absPath;
                try {
                    absPath = Paths.get(tryFile.getCanonicalPath());
                } catch (IOException ex) {
                    absPath = Paths.get(tryFile.getAbsolutePath());
                }
                Path thisPath = pwdPath.relativize(absPath);
                fileList.add(tryFile);
                pathList.add(thisPath);
            } catch (Exception ex) {
                // Simply continue.
            }
        }
        int numAccepted = fileList.size();
        if (numAccepted < 1) {
            throw new IllegalArgumentException(" No valid files could be found to open.");
        }
        allFiles = fileList.toArray(new File[numAccepted]);
        allPaths = pathList.toArray(new Path[numAccepted]);
        this.file = allFiles[0];
        this.filepath = allPaths[0];
        assemblies = new ArrayList<>();
        propertyList = new ArrayList<>();
    }

    /**
     * At present, parses the PDB, XYZ, INT, or ARC file from the constructor
     * and creates MolecularAssembly and properties objects.
     */
    @Override
    public void run() {
        int numFiles = allFiles.length;
        for (int i = 0; i < numFiles; i++) {
            File fileI = allFiles[i];
            Path pathI = allPaths[i];
            MolecularAssembly assembly = new MolecularAssembly(pathI.toString());
            assembly.setFile(fileI);

            CompositeConfiguration properties = Keyword.loadProperties(fileI);
            ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
            ForceField forceField = forceFieldFilter.parse();
            assembly.setForceField(forceField);
            SystemFilter filter;

            if (new PDBFileFilter().acceptDeep(fileI)) {
                filter = new PDBFilter(fileI, assembly, forceField, properties);
            } else if (new XYZFileFilter().acceptDeep(fileI)) {
                filter = new XYZFilter(fileI, assembly, forceField, properties);
            } else if (new INTFileFilter().acceptDeep(fileI) || new ARCFileFilter().accept(fileI)) {
                filter = new INTFilter(fileI, assembly, forceField, properties);
            } else {
                throw new IllegalArgumentException(String.format(" File %s could not be recognized as a valid PDB, XYZ, INT, or ARC file.", pathI.toString()));
            }
            filter.readFile();
            if (!(filter instanceof PDBFilter)) {
                Utilities.biochemistry(assembly, filter.getAtomList());
            }
            assembly.finalize(true, forceField);
            ForceFieldEnergy energy = new ForceFieldEnergy(assembly);
            assembly.setPotential(energy);
            assemblies.add(assembly);
            propertyList.add(properties);
        }
        activeAssembly = assemblies.get(0);
        activeProperties = propertyList.get(0);
    }

    /**
     * Returns the first MolecularAssembly created by the run() function.
     *
     * @return A MolecularAssembly
     */
    @Override
    public MolecularAssembly getAssembly() {
        return activeAssembly;
    }

    /**
     * Returns all MolecularAssembly objects created by this opener.
     *
     * @return Array of MolecularAssemblys
     */
    @Override
    public MolecularAssembly[] getAllAssemblies() {
        return assemblies.toArray(new MolecularAssembly[assemblies.size()]);
    }

    /**
     * Returns the properties associated with the first MolecularAssembly.
     *
     * @return Active properties
     */
    @Override
    public CompositeConfiguration getProperties() {
        return activeProperties;
    }

    /**
     * Returns the properties of all MolecularAssembly objects created by this
     * opener.
     *
     * @return Array of all properties
     */
    @Override
    public CompositeConfiguration[] getAllProperties() {
        return propertyList.toArray(new CompositeConfiguration[propertyList.size()]);
    }
}
