/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
 */
public class FFXLocale {

    private static final Logger logger = Logger.getLogger(FFXLocale.class.getName());
    private Locale currentLocale;
    private PropertyResourceBundle ffxLabels;
    private Hashtable<String, String> reverseLookUp = new Hashtable<String, String>();

    public FFXLocale() {
        currentLocale = Locale.getDefault();
        ffxLabels = (PropertyResourceBundle) ResourceBundle.getBundle(
                "ffx.ui.properties.StringBundle", currentLocale);
        loadHashtable();
    }

    public FFXLocale(String language, String country) {
        setLocale(language, country);
    }

    public String getKey(String string) {
        return reverseLookUp.get(string);
    }

    public String getValue(String key) {
        return ffxLabels.getString(key).trim();
    }

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
