/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import org.biojava.bio.structure.Structure;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.BiojavaFilter;
import ffx.potential.parsers.FileOpener;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.utilities.Keyword;

/**
 * The PotentialsDataConverter class describes a Runnable object which converts
 * some data structure to Force Field X MolecularAssembly(s).
 *
 * @author Jacob M. Litman
 *
 * @author Michael J. Schnieders
 */
public class PotentialsDataConverter implements FileOpener {

    private static final Logger logger = Logger.getLogger(PotentialsDataConverter.class.getName());

    private final File file; // Used to generate the Properties.
    private final Object dataStructure;
    // Should in the future allow a list of data structures.
    private final Utilities.DataType dataType;
    private List<MolecularAssembly> assemblies;
    private MolecularAssembly activeAssembly; // Presently, will just be the first element of assemblies.
    private List<CompositeConfiguration> propertyList;
    private CompositeConfiguration activeProperties;
    private final static String BIOJAVA_DEFAULT_FILENAME = "UNKNOWN_BIOJAVA_FILE";

    public PotentialsDataConverter(Object data) throws FileNotFoundException {
        this(data, null);
    }

    public PotentialsDataConverter(Object data, File file) throws FileNotFoundException {
        /* If the file is provided, just use that. Else, use a static get<Type>File
         * method to find the file: if that fails, or if the data object is not
         * of a recognized type, throws a relevant exception.
         */
        if (data instanceof Structure) {
            if (file == null) {
                this.file = getBiojavaFile((Structure) data);
            } else {
                this.file = file;
            }
            this.dataStructure = data;
            this.dataType = Utilities.DataType.BIOJAVA;
        } else { // Insert else-ifs for other types here.
            throw new IllegalArgumentException(" Data structure provided was not "
                    + "of a recognized type: must be a Biojava structure");
        }
    }

    /**
     * Switches between various default get-file methods for different data
     * structure types.
     *
     * @param data Data structure to find file file
     * @return Source file
     * @throws FileNotFoundException If no file could be found
     */
    public static File getDefaultFile(Object data) throws FileNotFoundException {
        if (data instanceof Structure) {
            return getBiojavaFile((Structure) data);
        } // Insert else-ifs for other data structures here.
        throw new FileNotFoundException("Could not find a file for data structure.");
    }

    /**
     * Attempt to get the file the Structure was loaded from. Assumes .pdb
     * format, mostly because Biojava can only load from PDB.
     *
     * @param structure A Biojava structure
     * @return pdbcode.pdb, name.pdb, or a default file
     * @throws java.io.FileNotFoundException If no file could be found.
     */
    public static File getBiojavaFile(Structure structure) throws FileNotFoundException {
        String filename = structure.getPDBCode();
        if (filename == null || filename.trim().equals("")) {
            filename = structure.getName();
            if (filename == null || filename.trim().equals("")) {
                filename = BIOJAVA_DEFAULT_FILENAME;
            }
        }
        filename += ".pdb";
        File retFile = new File(filename);
        int counter = 1;
        while (retFile.isDirectory() && counter < 1000) {
            retFile = new File(String.format("%s_%d", filename, counter++));
        }
        if (retFile.isDirectory()) {
            throw new FileNotFoundException(String.format(" Could not find a file "
                    + "for structure %s", structure.toString()));
        }
        return retFile;
    }

    /**
     * Converts the data structure to MolecularAssembly(s).
     */
    @Override
    public void run() {
        if (dataStructure == null || dataType.equals(Utilities.DataType.UNK)) {
            throw new IllegalArgumentException("Object passed was not recognized.");
        }
        assemblies = new ArrayList<>();
        propertyList = new ArrayList<>();
        switch (dataType) {
            case BIOJAVA:
                Structure struct = (Structure) dataStructure;
                String name = struct.getPDBCode();
                CompositeConfiguration properties = Keyword.loadProperties(file);
                MolecularAssembly assembly = new MolecularAssembly(name);
                assembly.setFile(file);
                ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
                ForceField forceField = forceFieldFilter.parse();
                assembly.setForceField(forceField);

                BiojavaFilter filter = new BiojavaFilter(struct, assembly, forceField, properties);
                if (filter.convert()) {
                    filter.applyAtomProperties();
                    assembly.finalize(true, forceField);
                    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(assembly, filter.getCoordRestraints());
                    assembly.setPotential(energy);
                    assemblies.add(assembly);
                    propertyList.add(properties);

                    List<Character> altLocs = filter.getAltLocs();
                    if (altLocs.size() > 1 || altLocs.get(0) != ' ') {
                        StringBuilder altLocString = new StringBuilder("\n Alternate locations found [ ");
                        for (Character c : altLocs) {
                            // Do not report the root conformer.
                            if (c == ' ') {
                                continue;
                            }
                            altLocString.append(format("(%s) ", c));
                        }
                        altLocString.append("]\n");
                        logger.info(altLocString.toString());
                    }

                    /**
                     * Alternate conformers may have different chemistry, so
                     * they each need to be their own MolecularAssembly.
                     */
                    for (Character c : altLocs) {
                        if (c.equals(' ') || c.equals('A')) {
                            continue;
                        }
                        MolecularAssembly newAssembly = new MolecularAssembly(name);
                        newAssembly.setForceField(assembly.getForceField());
                        filter.setAltID(newAssembly, c);
                        filter.clearSegIDs();
                        if (filter.convert()) {
                            String fileName = assembly.getFile().getAbsolutePath();
                            newAssembly.setName(FilenameUtils.getBaseName(fileName) + " " + c);
                            filter.applyAtomProperties();
                            newAssembly.finalize(true, assembly.getForceField());
                            energy = ForceFieldEnergy.energyFactory(newAssembly, filter.getCoordRestraints());
                            newAssembly.setPotential(energy);
                            assemblies.add(newAssembly);
                            properties.addConfiguration(properties);
                        }
                    }
                } else {
                    logger.warning(String.format(" Failed to convert structure %s", dataStructure.toString()));
                }
                activeAssembly = assembly;
                activeProperties = properties;
                break;
            case UNK:
            default:
                throw new IllegalArgumentException("Object passed was not recognized.");
        }
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
     * Returns the i'th MolecularAssembly
     *
     * @param i Index
     * @return The i'th MolecularAssembly
     */
    public MolecularAssembly getAssembly(int i) {
        return assemblies.get(i);
    }

    /**
     * Returns all MolecularAssembly objects created by this converter.
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
     * Returns the properties associated with the i'th MolecularAssembly
     *
     * @param i
     *
     * @return CompositeConfiguration for MolecularAssembly i.
     */
    public CompositeConfiguration getProperties(int i) {
        return propertyList.get(i);
    }

    /**
     * Returns the properties of all MolecularAssembly objects created by this
     * converter.
     *
     * @return Array of all properties
     */
    @Override
    public CompositeConfiguration[] getAllProperties() {
        return propertyList.toArray(new CompositeConfiguration[propertyList.size()]);
    }

}
