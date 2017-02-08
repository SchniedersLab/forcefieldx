/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.utilities;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.apache.commons.configuration.SystemConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * The Keyword class holds a single Force Field X keyword entry.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class Keyword {

    private static final Logger logger = Logger.getLogger(Keyword.class.getName());
    private String keyword = null;
    private Vector<String> data = null;

    /**
     * <p>
     * Constructor for Keyword.</p>
     *
     * @param k a {@link java.lang.String} object.
     */
    public Keyword(String k) {
        keyword = k;
        data = new Vector<>();
    }

    /**
     * <p>
     * Constructor for Keyword.</p>
     *
     * @param k a {@link java.lang.String} object.
     * @param entry a {@link java.lang.String} object.
     */
    public Keyword(String k, String entry) {
        this(k);
        data.add(entry);
    }

    /**
     * <p>
     * Constructor for Keyword.</p>
     *
     * @param k a {@link java.lang.String} object.
     * @param entry an array of {@link java.lang.String} objects.
     */
    public Keyword(String k, String entry[]) {
        this(k);
        data.addAll(Arrays.asList(entry));
    }

    /**
     * <p>
     * append</p>
     *
     * @param entry a {@link java.lang.String} object.
     */
    public void append(String entry) {
        data.add(entry);
    }

    /**
     * <p>
     * append</p>
     *
     * @param entry an array of {@link java.lang.String} objects.
     */
    public void append(String entry[]) {
        data.addAll(Arrays.asList(entry));
    }

    /**
     * <p>
     * clear</p>
     */
    public void clear() {
        data.clear();
    }

    /**
     * <p>
     * getEntries</p>
     *
     * @return a {@link java.util.Vector} object.
     */
    public Vector<String> getEntries() {
        return data;
    }

    /**
     * <p>
     * getEntry</p>
     *
     * @param i a int.
     * @return a {@link java.lang.String} object.
     */
    public String getEntry(int i) {
        return data.get(i);
    }

    /**
     * <p>
     * Getter for the field <code>keyword</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getKeyword() {
        return keyword;
    }

    /**
     * <p>
     * print</p>
     */
    public void print() {
        logger.info(this.toString());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(keyword + " ");
        for (String s : data) {
            sb.append(s);
        }
        return sb.toString();
    }

    /**
     * This method sets up configuration properties in the following precedence
     * * order:
     *
     * 1.) Structure specific properties (for example pdbname.properties)
     *
     * 2.) Java system properties a.) -Dkey=value from the Java command line b.)
     * System.setProperty("key","value") within Java code.
     *
     * 3.) User specific properties (~/.ffx/ffx.properties)
     *
     * 4.) System wide properties (file defined by environment variable
     * FFX_PROPERTIES)
     *
     * 5.) Internal force field definition.
     *
     * @since 1.0
     * @param file a {@link java.io.File} object.
     * @return a {@link org.apache.commons.configuration.CompositeConfiguration}
     * object.
     */
    public static CompositeConfiguration loadProperties(File file) {
        /**
         * Command line options take precedence.
         */
        CompositeConfiguration properties = new CompositeConfiguration();

        /**
         * Structure specific options are first.
         */
        if (file != null) {
            String structureBasename = FilenameUtils.removeExtension(file.getAbsolutePath());
            String propertyFilename
                    = (new File(structureBasename + ".properties").exists()) ? structureBasename + ".properties"
                    : (new File(structureBasename + ".prop").exists()) ? structureBasename + ".prop"
                    : (new File(structureBasename + ".key").exists()) ? structureBasename + ".key"
                    : null;
            if (propertyFilename != null) {
                File structurePropFile = new File(propertyFilename);
                if (structurePropFile.canRead()) {
                    try {
                        properties.addConfiguration(new PropertiesConfiguration(structurePropFile));
                        properties.addProperty("propertyFile", structurePropFile.getCanonicalPath());
                    } catch (ConfigurationException | IOException e) {
                        logger.log(Level.INFO, " Error loading {0}.", structureBasename);
                    }
                }
            }
        }

        /**
         * Java system properties
         *
         * a.) -Dkey=value from the Java command line
         *
         * b.) System.setProperty("key","value") within Java code.
         */
        properties.addConfiguration(new SystemConfiguration());

        /**
         * User specific options are 3rd.
         */
        String filename = System.getProperty("user.home") + File.separator + ".ffx/ffx.properties";
        File userPropFile = new File(filename);
        if (userPropFile.exists() && userPropFile.canRead()) {
            try {
                properties.addConfiguration(new PropertiesConfiguration(userPropFile));
            } catch (ConfigurationException e) {
                logger.log(Level.INFO, " Error loading {0}.", filename);
            }
        }

        /**
         * System wide options are 2nd to last.
         */
        filename = System.getenv("FFX_PROPERTIES");
        if (filename != null) {
            File systemPropFile = new File(filename);
            if (systemPropFile.exists() && systemPropFile.canRead()) {
                try {
                    properties.addConfiguration(new PropertiesConfiguration(systemPropFile));
                } catch (ConfigurationException e) {
                    logger.log(Level.INFO, " Error loading {0}.", filename);
                }
            }
        }

        /**
         * Echo the interpolated configuration.
         */
        if (logger.isLoggable(Level.FINE)) {
            //Configuration config = properties.interpolatedConfiguration();
            Iterator<String> i = properties.getKeys();
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\n %-30s %s\n", "Property", "Value"));
            while (i.hasNext()) {
                String s = i.next();
                //sb.append(String.format(" %-30s %s\n", s, Arrays.toString(config.getList(s).toArray())));
                sb.append(String.format(" %-30s %s\n", s, Arrays.toString(properties.getList(s).toArray())));
            }
            logger.fine(sb.toString());
        }

        return properties;
    }
}
