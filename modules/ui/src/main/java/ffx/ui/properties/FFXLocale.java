/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.ui.properties;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Locale;
import java.util.PropertyResourceBundle;
import java.util.ResourceBundle;
import java.util.logging.Logger;

/**
 * The FFXLocale class will encapsulate internationalization features.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class FFXLocale {

    private static final Logger logger = Logger.getLogger(FFXLocale.class.getName());
    private Locale currentLocale;
    private PropertyResourceBundle ffxLabels;
    private Hashtable<String, String> reverseLookUp = new Hashtable<String, String>();

    /**
     * <p>
     * Constructor for FFXLocale.</p>
     */
    public FFXLocale() {
        currentLocale = Locale.getDefault();
        ffxLabels = (PropertyResourceBundle) ResourceBundle.getBundle(
                "ffx.ui.properties.StringBundle", currentLocale);
        loadHashtable();
    }

    /**
     * <p>
     * Constructor for FFXLocale.</p>
     *
     * @param language a {@link java.lang.String} object.
     * @param country a {@link java.lang.String} object.
     */
    public FFXLocale(String language, String country) {
        setLocale(language, country);
    }

    /**
     * <p>
     * getKey</p>
     *
     * @param string a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getKey(String string) {
        return reverseLookUp.get(string);
    }

    /**
     * <p>
     * getValue</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getValue(String key) {
        return ffxLabels.getString(key).trim();
    }

    /**
     * <p>
     * list</p>
     */
    public void list() {
        for (String value : reverseLookUp.keySet()) {
            String key = reverseLookUp.get(value);
            logger.info("key = " + key + ", " + "value = " + value);
        }
    }

    private void loadHashtable() {
        reverseLookUp.clear();
        Enumeration<String> e = ffxLabels.getKeys();
        while (e.hasMoreElements()) {
            String key = e.nextElement();
            String value = getValue(key);
            reverseLookUp.put(value, key);
        }
    }

    /**
     * <p>
     * setLocale</p>
     *
     * @param language a {@link java.lang.String} object.
     * @param country a {@link java.lang.String} object.
     * @return a boolean.
     */
    public boolean setLocale(String language, String country) {
        Locale locale = new Locale(language, country);
        try {
            ffxLabels = (PropertyResourceBundle) ResourceBundle.getBundle(
                    "ffx.ui.properties.StringBundle", locale);
        } catch (Exception ex) {
            Logger.getLogger("ffx").severe("" + ex);
            return false;
        }
        loadHashtable();
        currentLocale = locale;
        return true;
    }
}
