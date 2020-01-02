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
package ffx.realspace.parsers;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.io.FilenameUtils.isExtension;
import static org.apache.commons.io.FilenameUtils.removeExtension;

import ffx.potential.MolecularAssembly;

/**
 * <p>
 * RealSpaceFile class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class RealSpaceFile {

    private static final Logger logger = Logger.getLogger(RealSpaceFile.class.getName());

    private final String filename;
    private final double weight;
    private final RealSpaceFileFilter realSpaceFileFilter;

    /**
     * Read in a Real Space density file and set weight set to 1.0.
     *
     * @param filename file name to read in
     */
    public RealSpaceFile(String filename) {
        this(filename, 1.0);
    }

    /**
     * Read in a Real Space density file.
     *
     * @param fileName file name to read in
     * @param weight   the weight of the data
     */
    public RealSpaceFile(String fileName, double weight) {
        File tmp = new File(fileName);
        if (!tmp.exists()) {
            logger.severe(" Data file: " + fileName + " not found!");
        }

        if (isExtension(fileName, "map")) {
            realSpaceFileFilter = new CCP4MapFilter();
        } else {
            realSpaceFileFilter = null;
        }

        this.filename = fileName;
        this.weight = weight;
    }

    /**
     * Read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0 and neutron value of false.
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 filename will be determined.
     */
    public RealSpaceFile(MolecularAssembly[] assembly) {
        this(assembly[0], 1.0);
    }

    /**
     * Read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0.
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 filename will be determined.
     */
    public RealSpaceFile(MolecularAssembly assembly) {
        this(assembly, 1.0);
    }

    /**
     * Read in a Real Space density file based on the molecular assembly
     * filename, using a weight of 1.0.
     *
     * @param assembly {@link ffx.potential.MolecularAssembly} from which a
     *                 filename will be determined
     * @param weight   the weight of the data
     */
    public RealSpaceFile(MolecularAssembly assembly, double weight) {
        String name = removeExtension(assembly.getFile().getPath());

        File tmp = new File(name + ".map");
        if (tmp.exists()) {
            logger.info(" Data file: " + tmp.getName());
            realSpaceFileFilter = new CCP4MapFilter();
        } else {
            logger.severe(" No input data was found.");
            realSpaceFileFilter = null;
        }

        String filenameHolder;
        try {
            Path filePath = Paths.get(tmp.getCanonicalPath());
            Path pwdPath = Paths.get(new File("").getCanonicalPath());
            filenameHolder = pwdPath.relativize(filePath).toString();
        } catch (IOException e) {
            String message
                    = " Relative path to provided data file could not be resolved: using map file name instead.";
            logger.log(Level.WARNING, message, e);
            filenameHolder = tmp.getName();
        }
        this.filename = filenameHolder;
        this.weight = weight;
    }

    /**
     * Return the weight of this dataset.
     *
     * @return weight The weight (wA).
     */
    public double getWeight() {
        return weight;
    }

    /**
     * <p>Getter for the field <code>filename</code>.</p>
     *
     * @return the filename
     */
    public String getFilename() {
        return filename;
    }

    /**
     * <p>Getter for the field <code>realSpaceFileFilter</code>.</p>
     *
     * @return the realSpaceFileFilter
     */
    public RealSpaceFileFilter getRealSpaceFileFilter() {
        return realSpaceFileFilter;
    }
}
