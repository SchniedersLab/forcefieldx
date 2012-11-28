/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
package ffx.xray;

import java.io.File;
import java.util.logging.Logger;

import static org.apache.commons.io.FilenameUtils.*;

import ffx.potential.bonded.MolecularAssembly;

/**
 * <p>DiffractionFile class.</p>
 *
 * @author Tim Fenn
 *
 */
public class DiffractionFile {

    private static final Logger logger = Logger.getLogger(DiffractionFile.class.getName());
    protected final String filename;
    protected final double weight;
    protected final boolean neutron;
    protected final DiffractionFileFilter diffractionfilter;

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
     * @param weight the weight of the data
     */
    public DiffractionFile(String filename, double weight) {
        this(filename, weight, false);
    }

    /**
     * read in a diffraction file
     *
     * @param filename file name to read in
     * @param weight the weight of the data
     * @param neutron if true, this is a neutron data set
     */
    public DiffractionFile(String filename, double weight, boolean neutron) {
        File tmp = new File(filename);
        if (!tmp.exists()) {
            logger.severe("data file: " + filename + " not found!");
        }

        if (isExtension(filename, "mtz")) {
            diffractionfilter = new MTZFilter();
        } else if (isExtension(filename, new String[]{"cif", "ent", "sf"})) {
            diffractionfilter = new CIFFilter();
        } else if (isExtension(filename, new String[]{"cns", "hkl"})) {
            diffractionfilter = new CNSFilter();
        } else {
            diffractionfilter = null;
        }

        if (diffractionfilter == null) {
            logger.severe("input data file format not recognized!\n Please use an appropriate file extension: [MTZ: mtz] [CIF: cif, ent, sf] [CNS: cns, hkl] to identify your data file type!");
        }

        this.filename = filename;
        this.weight = weight;
        this.neutron = neutron;
    }

    /**
     * read in a diffraction file based on the molecular assembly filename,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     */
    public DiffractionFile(MolecularAssembly assembly[]) {
        this(assembly[0], 1.0, false);
    }

    /**
     * <p>Constructor for DiffractionFile.</p>
     *
     * @param assembly an array of
     * {@link ffx.potential.bonded.MolecularAssembly} objects.
     * @param weight a double.
     */
    public DiffractionFile(MolecularAssembly assembly[], double weight) {
        this(assembly[0], weight, false);
    }

    /**
     * <p>Constructor for DiffractionFile.</p>
     *
     * @param assembly an array of
     * {@link ffx.potential.bonded.MolecularAssembly} objects.
     * @param weight a double.
     * @param neutron a boolean.
     */
    public DiffractionFile(MolecularAssembly assembly[], double weight,
            boolean neutron) {
        this(assembly[0], weight, neutron);
    }

    /**
     * read in a diffraction file based on the molecular assembly filename,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     */
    public DiffractionFile(MolecularAssembly assembly) {
        this(assembly, 1.0, false);
    }

    /**
     * read in a diffraction file based on the molecular assembly filename,
     * using a neutron value of false
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     * @param weight the weight of the data
     */
    public DiffractionFile(MolecularAssembly assembly, double weight) {
        this(assembly, weight, false);
    }

    /**
     * read in a diffraction file based on the molecular assembly filename,
     * using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     * @param weight the weight of the data
     * @param neutron if true, this is a neutron data set
     */
    public DiffractionFile(MolecularAssembly assembly, double weight,
            boolean neutron) {
        String name = removeExtension(assembly.getFile().getPath());

        File tmp = new File(name + ".mtz");
        if (tmp.exists()) {
            logger.info("data file: " + tmp.getName());
            diffractionfilter = new MTZFilter();
        } else {
            tmp = new File(name + ".cif");
            if (tmp.exists()) {
                logger.info("data file: " + tmp.getName());
                diffractionfilter = new CIFFilter();
            } else {
                tmp = new File(name + ".ent");
                if (tmp.exists()) {
                    logger.info("data file: " + tmp.getName());
                    diffractionfilter = new CIFFilter();
                } else {
                    tmp = new File(name + ".sf");
                    if (tmp.exists()) {
                        logger.info("data file: " + tmp.getName());
                        diffractionfilter = new CIFFilter();
                    } else {
                        tmp = new File(name + ".cns");
                        if (tmp.exists()) {
                            logger.info("data file: " + tmp.getName());
                            diffractionfilter = new CNSFilter();
                        } else {
                            tmp = new File(name + ".hkl");
                            if (tmp.exists()) {
                                logger.info("data file: " + tmp.getName());
                                diffractionfilter = new CNSFilter();
                            } else {
                                logger.severe("no input data found!");
                                diffractionfilter = null;
                            }
                        }
                    }
                }
            }
        }

        this.filename = tmp.getName();
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
}
