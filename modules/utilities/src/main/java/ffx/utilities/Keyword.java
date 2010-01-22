/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.utilities;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.apache.commons.configuration.SystemConfiguration;

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

    public Keyword(String k) {
        keyword = k;
        data = new Vector<String>();
    }

    public Keyword(String k, String entry) {
        this(k);
        data.add(entry);
    }

    public Keyword(String k, String entry[]) {
        this(k);
        for (String s : entry) {
            data.add(s);
        }
    }

    public void append(String entry) {
        data.add(entry);
    }

    public void append(String entry[]) {
        for (String s : entry) {
            data.add(s);
        }
    }

    public void clear() {
        data.clear();
    }

    public Vector<String> getEntries() {
        return data;
    }

    public String getEntry(int i) {
        return data.get(i);
    }

    public String getKeyword() {
        return keyword;
    }

    public void print() {
        logger.info(this.toString());
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer(keyword + " ");
        for (String s : data) {
            sb.append(s);
        }
        return sb.toString();
    }

    /**
     * This method sets up configuration properties in the following precedence
     * order:
     * 1.) Java system properties
     *     a.) -Dkey=value from the Java command line
     *     b.) System.setProperty("key","value") within Java code.
     *
     * 2.) Structure specific properties (for example pdbname.properties)
     *
     * 3.) User specific properties (~/.ffx/ffx.properties)
     *
     * 4.) System wide properties (file defined by environment variable FFX_PROPERTIES)
     *
     * 5.) Internal force field definition.
     */
    public static CompositeConfiguration loadProperties(File file) {
        /**
         * Command line options take precedences.
         */
        CompositeConfiguration properties = new CompositeConfiguration();
        properties.addConfiguration(new SystemConfiguration());

        /**
         * Structure specific options are 2nd.
         */
        if (file != null) {
            String filename = file.getAbsolutePath();
            filename = org.apache.commons.io.FilenameUtils.removeExtension(filename);
            String keyFilename = filename + ".properties";
            File structurePropFile = new File(keyFilename);
            if (structurePropFile.exists() && structurePropFile.canRead()) {
                try {
                    properties.addConfiguration(new PropertiesConfiguration(structurePropFile));
                } catch (Exception e) {
                    logger.info("Error loading " + filename + ".");
                }
            } else {
                keyFilename = filename + ".key";
                structurePropFile = new File(keyFilename);
                if (structurePropFile.exists() && structurePropFile.canRead()) {
                    try {
                        properties.addConfiguration(new PropertiesConfiguration(structurePropFile));
                    } catch (Exception e) {
                        logger.info("Error loading " + filename + ".");
                    }
                }
            }
        }

        /**
         * User specific options are 3rd.
         */
        String filename = System.getProperty("user.home") + File.separator + ".ffx/ffx.properties";
        File userPropFile = new File(filename);
        if (userPropFile.exists() && userPropFile.canRead()) {
            try {
                properties.addConfiguration(new PropertiesConfiguration(userPropFile));
            } catch (Exception e) {
                logger.info("Error loading " + filename + ".");
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
                } catch (Exception e) {
                    logger.info("Error loading " + filename + ".");
                }
            }
        }

        /**
         * Echo the interpolated configuration.
         */
        if (logger.isLoggable(Level.FINE)) {
            Configuration config = properties.interpolatedConfiguration();
            Iterator<String> i = config.getKeys();
            while (i.hasNext()) {
                String s = i.next();
                logger.fine("Key: " + s + ", Value: " + Arrays.toString(config.getList(s).toArray()));
            }
        }

        return properties;
    }
}
