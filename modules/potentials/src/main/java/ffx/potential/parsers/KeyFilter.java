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
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.logging.Logger;

import ffx.utilities.Keyword;

/**
 * The KeyFilter class parses Force Field X Keyword (*.KEY) and Property
 * (*.PROPERTIES) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class KeyFilter {

    private static final Logger logger = Logger.getLogger(KeyFilter.class.getName());

    /**
     * <p>loadSystemKeywords</p>
     *
     * @return a {@link java.util.Hashtable} object.
     */
    public static Hashtable<String, Keyword> loadSystemKeywords() {
        File f = new File("/etc/ffx.conf");
        Hashtable<String, Keyword> systemKeywords = new Hashtable<String, Keyword>();
        if (f.exists() && f.canRead()) {
            logger.info("Reading /etc/ffx.conf");
            systemKeywords = KeyFilter.open(f, systemKeywords);
        }
        String path = System.getProperty("user.home") + File.separator + ".ffx/ffx.conf";
        f = new File(path);
        if (f.exists() && f.canRead()) {
            logger.info("Reading " + path);
            systemKeywords = KeyFilter.open(f, systemKeywords);
        }
        return systemKeywords;
    }

    /**
     * <p>open</p>
     *
     * @param keyFile a {@link java.io.File} object.
     * @return a {@link java.util.Hashtable} object.
     */
    public static Hashtable<String, Keyword> open(File keyFile) {
        if (keyFile == null || !keyFile.exists() || !keyFile.canRead()) {
            return null;
        }
        Hashtable<String, Keyword> keywordHash = loadSystemKeywords();
        return open(keyFile, keywordHash);
    }

    /**
     * <p>open</p>
     *
     * @param keyFile a {@link java.io.File} object.
     * @param keywordHash a {@link java.util.Hashtable} object.
     * @return a {@link java.util.Hashtable} object.
     */
    public static Hashtable<String, Keyword> open(File keyFile, Hashtable<String, Keyword> keywordHash) {
        if (keyFile == null || !keyFile.exists() || !keyFile.canRead()) {
            return null;
        }
        if (keywordHash == null) {
            keywordHash = new Hashtable<String, Keyword>();
        }
        FileReader fr = null;
        BufferedReader br = null;
        try {
            fr = new FileReader(keyFile);
            br = new BufferedReader(fr);
            Keyword comments = new Keyword("COMMENTS");
            keywordHash.put("COMMENTS", comments);
            while (br.ready()) {
                String s = br.readLine();
                if (s == null) {
                    continue;
                }
                s = s.trim();
                if (s.equals("")) {
                    continue; // Skip blank lines
                }
                // Store comments together
                if (s.startsWith("#") || s.toUpperCase().startsWith("ECHO")) {
                    comments.append(s);
                } else {
                    int firstspace = s.indexOf(" ");
                    String keyword, data;
                    if (firstspace == -1) { // no parameters
                        keyword = s.trim().toUpperCase();
                        // Rattle is special case, because it can be active
                        // without being checked
                        // Valid Key files can have: RATTLE
                        // or RATTLE BONDS
                        // or RATTLE & RATTLE BONDS as seperate lines
                        // Each of these valid cases mean different things...
                        if (keyword.equalsIgnoreCase("rattle")) {
                            data = "RATTLE";
                        } else {
                            data = null;
                        }
                    } else {
                        keyword = s.substring(0, firstspace).toUpperCase();
                        data = s.substring(firstspace).trim();
                    }
                    Keyword kd = keywordHash.get(keyword);
                    if (kd == null) {
                        kd = new Keyword(keyword);
                        keywordHash.put(keyword, kd);
                    }
                    if (data != null) {
                        kd.append(data);
                    }
                    /**
                     * Multipoles and TORTORS are the only keywords that span
                     * multiple lines. Editing these from within Force Field
                     * Xplor seems unlikely, so they are treated as comments.
                     */
                    if (keyword.equalsIgnoreCase("MULTIPOLE")) {
                        int mnum[] = {3, 1, 2, 3};
                        for (int i = 0; i < 4; i++) {
                            if (!br.ready()) {
                                System.out.println("Check for an invalid MULTIPOLE keyword.");
                                return null;
                            }
                            s = br.readLine();
                            if (s == null) {
                                logger.warning("Multipole format error.");
                                return null;
                            }
                            s = s.trim();
                            if (s.split(" +").length != mnum[i]) {
                                logger.warning("Multipole format error.");
                                return null;
                            }
                            kd.append(s);
                        }
                    } else if (keyword.equalsIgnoreCase("TORTORS")) {
                        String res[] = data.split(" +");
                        if (res == null || res.length < 7) {
                            logger.warning("TORTOR format error.");
                            return null;
                        }
                        int xres = Integer.parseInt(res[5]);
                        int yres = Integer.parseInt(res[6]);
                        for (int i = 0; i < xres * yres; i++) {
                            if (!br.ready()) {
                                System.out.println("Check for an invalid TORTOR keyword.");
                                return null;
                            }
                            s = br.readLine();
                            if (s == null) {
                                logger.warning("TORTOR format error.");
                                return null;
                            }
                            s = s.trim();
                            if (s.split(" +").length != 3) {
                                logger.warning("TORTOR format error.");
                                return null;
                            }
                            kd.append(s);
                        }
                    }
                }
            }
            return keywordHash;
        } catch (IOException e) {
            System.err.println("Error reading Key File: " + e);
            return null;
        } finally {
            if (fr != null) {
                try {
                    fr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * <p>Constructor for KeyFilter.</p>
     */
    public KeyFilter() {
    }
}
