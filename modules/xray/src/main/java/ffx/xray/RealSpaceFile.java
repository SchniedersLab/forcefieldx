/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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

import static org.apache.commons.io.FilenameUtils.isExtension;
import static org.apache.commons.io.FilenameUtils.removeExtension;

import ffx.potential.bonded.MolecularAssembly;

/**
 * <p>RealSpaceFile class.</p>
 *
 * @author Tim Fenn
 *
 */
public class RealSpaceFile {

    private static final Logger logger = Logger.getLogger(RealSpaceFile.class.getName());
    protected final String filename;
    protected final double weight;
    protected final RealSpaceFileFilter realspacefilter;

    /**
     * read in a Real Space density file, weight set to 1.0
     *
     * @param filename file name to read in
     */
    public RealSpaceFile(String filename) {
        this(filename, 1.0);
    }

    /**
     * read in a Real Space density file
     *
     * @param filename file name to read in
     * @param weight the weight of the data
     */
    public RealSpaceFile(String filename, double weight) {
        File tmp = new File(filename);
        if (!tmp.exists()) {
            logger.severe(" Data file: " + filename + " not found!");
        }

        if (isExtension(filename, "map")) {
            realspacefilter = new CCP4MapFilter();
        } else {
            realspacefilter = null;
        }

        this.filename = filename;
        this.weight = weight;
    }

    /**
     * read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0 and neutron value of false
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     */
    public RealSpaceFile(MolecularAssembly assembly[]) {
        this(assembly[0], 1.0);
    }

    /**
     * read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     */
    public RealSpaceFile(MolecularAssembly assembly) {
        this(assembly, 1.0);
    }

    /**
     * read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly} from which
     * a filename will be determined
     * @param weight the weight of the data
     */
    public RealSpaceFile(MolecularAssembly assembly, double weight) {
        String name = removeExtension(assembly.getFile().getPath());

        File tmp = new File(name + ".map");
        if (tmp.exists()) {
            logger.info(" Data file: " + tmp.getName());
            realspacefilter = new CCP4MapFilter();
        } else {
            logger.severe(" No input data was found.");
            realspacefilter = null;
        }

        this.filename = tmp.getName();
        this.weight = weight;
    }

    /**
     * return the weight of this dataset
     *
     * @return weight wA
     */
    public double getWeight() {
        return weight;
    }
}
