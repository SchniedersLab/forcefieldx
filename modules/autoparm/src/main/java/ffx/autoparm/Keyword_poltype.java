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
package ffx.autoparm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
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
 * @author Gaurav Chattree and Michael J. Schnieders
 * @since 1.0
 *
 */
public class Keyword_poltype {

    private static final Logger logger = Logger.getLogger(Keyword_poltype.class.getName());
    private String keyword = null;
    private Vector<String> data = null;

    /**
     * <p>
     * Constructor for Keyword_poltype.</p>
     *
     * @param k a {@link java.lang.String} object.
     */
    public Keyword_poltype(String k) {
        keyword = k;
        data = new Vector<String>();
    }

    /**
     * <p>
     * Constructor for Keyword_poltype.</p>
     *
     * @param k a {@link java.lang.String} object.
     * @param entry a {@link java.lang.String} object.
     */
    public Keyword_poltype(String k, String entry) {
        this(k);
        data.add(entry);
    }

    /**
     * <p>
     * Constructor for Keyword_poltype.</p>
     *
     * @param k a {@link java.lang.String} object.
     * @param entry an array of {@link java.lang.String} objects.
     */
    public Keyword_poltype(String k, String entry[]) {
        this(k);
        for (String s : entry) {
            data.add(s);
        }
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
        for (String s : entry) {
            data.add(s);
        }
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
     * order: 1.) Java system properties a.) -Dkey=value from the Java command
     * line b.) System.setProperty("key","value") within Java code.
     *
     * 2.) Structure specific properties (for example pdbname.properties)
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
         * Command line options take precedences.
         */
        CompositeConfiguration properties = new CompositeConfiguration();
        properties.addConfiguration(new SystemConfiguration());

        /**
         * Structure specific options are 2nd.
         */
        if (file != null && file.exists() && file.canRead()) {

            try {
                BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
                String prmfilename = br.readLine().split(" +")[1];
                File prmfile = new File(prmfilename);
                if (prmfile.exists() && prmfile.canRead()) {
                    properties.addConfiguration(new PropertiesConfiguration(prmfile));
                    properties.addProperty("propertyFile", prmfile.getCanonicalPath());
                }
            } catch (Exception e1) {
                e1.printStackTrace();
            }
        }

//		/**
//		 * User specific options are 3rd.
//		 */
//		String filename = System.getProperty("user.home") + File.separator
//				+ ".ffx/ffx.properties";
//		File userPropFile = new File(filename);
//		if (userPropFile.exists() && userPropFile.canRead()) {
//			try {
//				properties.addConfiguration(new PropertiesConfiguration(
//						userPropFile));
//			} catch (Exception e) {
//				logger.info("Error loading " + filename + ".");
//			}
//		}
//
//		/**
//		 * System wide options are 2nd to last.
//		 */
//		filename = System.getenv("FFX_PROPERTIES");
//		if (filename != null) {
//			File systemPropFile = new File(filename);
//			if (systemPropFile.exists() && systemPropFile.canRead()) {
//				try {
//					properties.addConfiguration(new PropertiesConfiguration(
//							systemPropFile));
//				} catch (Exception e) {
//					logger.info("Error loading " + filename + ".");
//				}
//			}
//		}
        /**
         * Echo the interpolated configuration.
         */
        if (logger.isLoggable(Level.FINE)) {
            Configuration config = properties.interpolatedConfiguration();
            Iterator<String> i = config.getKeys();
            while (i.hasNext()) {
                String s = i.next();
                logger.fine("Key: " + s + ", Value: "
                        + Arrays.toString(config.getList(s).toArray()));
            }
        }

        return properties;
    }
}
