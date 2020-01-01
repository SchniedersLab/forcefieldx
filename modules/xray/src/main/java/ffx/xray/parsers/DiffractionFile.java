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
package ffx.xray.parsers;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Logger;

import static org.apache.commons.io.FilenameUtils.isExtension;
import static org.apache.commons.io.FilenameUtils.removeExtension;

import ffx.potential.MolecularAssembly;

/**
 * <p>
 * DiffractionFile class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class DiffractionFile {

    private static final Logger logger = Logger.getLogger(DiffractionFile.class.getName());

    private final String fileName;
    private final double weight;
    private final boolean neutron;
    private final DiffractionFileFilter diffractionFilter;

    /**
     * read in a diffraction file, weight set to 1.0 and neutron value of false
     *
     * @param filename file name to read in
     */
    public DiffractionFile(String filename) {
        this(filename, 1.0, false);
    }

    /**
     * read in a diffraction file, neutron value set to false
     *
     * @param filename file name to read in
     * @param weight   the weight of the data
     */
    public DiffractionFile(String filename, double weight) {
        this(filename, weight, false);
    }

    /**
     * read in a diffraction file
     *
     * @param filename file name to read in
     * @param weight   the weight of the data
     * @param neutron  if true, this is a neutron data set
     */
    public DiffractionFile(String filename, double weight, boolean neutron) {
        File tmp = new File(filename);
        if (!tmp.exists()) {
            logger.severe(" Data file: " + filename + " not found!");
        }

        if (isExtension(filename, "mtz")) {
            diffractionFilter = new MTZFilter();
        } else if (isExtension(filename, new String[]{"cif", "ent", "sf"})) {
            diffractionFilter = new CIFFilter();
        } else if (isExtension(filename, new String[]{"cns", "hkl"})) {
            diffractionFilter = new CNSFilter();
        } else {
            diffractionFilter = null;
        }

        if (diffractionFilter == null) {
            logger.severe("input data file format not recognized!\n Please use an appropriate file extension: [MTZ: mtz] [CIF: cif, ent, sf] [CNS: cns, hkl] to identify your data file type!");
        }

        this.fileName = filename;
        this.weight = weight;
        this.neutron = neutron;
    }

    /**
     * read in a diffraction file based on the molecular assembly fileName,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 fileName will be determined
     */
    public DiffractionFile(MolecularAssembly[] assembly) {
        this(assembly[0], 1.0, false);
    }

    /**
     * <p>
     * Constructor for DiffractionFile.</p>
     *
     * @param assembly an array of {@link ffx.potential.MolecularAssembly}
     *                 objects.
     * @param weight   a double.
     */
    public DiffractionFile(MolecularAssembly[] assembly, double weight) {
        this(assembly[0], weight, false);
    }

    /**
     * <p>
     * Constructor for DiffractionFile.</p>
     *
     * @param assembly an array of {@link ffx.potential.MolecularAssembly}
     *                 objects.
     * @param weight   a double.
     * @param neutron  a boolean.
     */
    public DiffractionFile(MolecularAssembly[] assembly, double weight,
                           boolean neutron) {
        this(assembly[0], weight, neutron);
    }

    /**
     * read in a diffraction file based on the molecular assembly fileName,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 fileName will be determined
     */
    public DiffractionFile(MolecularAssembly assembly) {
        this(assembly, 1.0, false);
    }

    /**
     * read in a diffraction file based on the molecular assembly fileName,
     * using a neutron value of false
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 fileName will be determined
     * @param weight   the weight of the data
     */
    public DiffractionFile(MolecularAssembly assembly, double weight) {
        this(assembly, weight, false);
    }

    /**
     * read in a diffraction file based on the molecular assembly fileName,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 fileName will be determined
     * @param weight   the weight of the data
     * @param neutron  if true, this is a neutron data set
     */
    public DiffractionFile(MolecularAssembly assembly, double weight,
                           boolean neutron) {
        String name = removeExtension(assembly.getFile().getPath());

        File tmp = new File(name + ".mtz");
        if (tmp.exists()) {
            //logger.info("\n Data file: " + tmp.getName());
            diffractionFilter = new MTZFilter();
        } else {
            tmp = new File(name + ".cif");
            if (tmp.exists()) {
                //logger.info("\n Data file: " + tmp.getName());
                diffractionFilter = new CIFFilter();
            } else {
                tmp = new File(name + ".ent");
                if (tmp.exists()) {
                    //logger.info("\n Data file: " + tmp.getName());
                    diffractionFilter = new CIFFilter();
                } else {
                    tmp = new File(name + ".sf");
                    if (tmp.exists()) {
                        //logger.info("\n Data file: " + tmp.getName());
                        diffractionFilter = new CIFFilter();
                    } else {
                        tmp = new File(name + ".cns");
                        if (tmp.exists()) {
                            //logger.info("\n Data file: " + tmp.getName());
                            diffractionFilter = new CNSFilter();
                        } else {
                            tmp = new File(name + ".hkl");
                            if (tmp.exists()) {
                                //logger.info("\n Data file: " + tmp.getName());
                                diffractionFilter = new CNSFilter();
                            } else {
                                logger.severe("no input data found!");
                                diffractionFilter = null;
                            }
                        }
                    }
                }
            }
        }
        String filenameHolder; // Compiler complains if I set this.fileName directly.
        try {
            Path filepath = Paths.get(tmp.getCanonicalPath());
            Path pwdPath = Paths.get(new File("").getCanonicalPath());
            filenameHolder = pwdPath.relativize(filepath).toString();
        } catch (IOException ex) {
            logger.warning(" Relative path to provided data file could not be resolved: using data file name instead.");
            filenameHolder = tmp.getName();
        }
        this.fileName = filenameHolder;
        this.weight = weight;
        this.neutron = neutron;
    }

    /**
     * return the weight of this dataset
     *
     * @return weight wA
     */
    public double getWeight() {
        return weight;
    }

    /**
     * is this a neutron dataset?
     *
     * @return true if neutron
     */
    public boolean isNeutron() {
        return neutron;
    }

    /**
     * <p>getFilename.</p>
     *
     * @return the fileName
     */
    public String getFilename() {
        return fileName;
    }

    /**
     * <p>getDiffractionfilter.</p>
     *
     * @return the diffractionFilter
     */
    public DiffractionFileFilter getDiffractionfilter() {
        return diffractionFilter;
    }

}
