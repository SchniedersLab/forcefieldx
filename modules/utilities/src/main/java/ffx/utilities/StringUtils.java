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
package ffx.utilities;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.util.FastMath;

/**
 * <p>StringUtils class.</p>
 *
 * @author Michael Schnieders
 */
public class StringUtils {

    /**
     * Creates a writer for text to a Gzip file.
     *
     * @param file Gzip file to write to.
     * @return A Writer
     * @throws java.io.IOException Thrown if creation of the GZip Writer fails.
     */
    public static Writer createGzipWriter(File file) throws IOException {
        return createGzipWriter(file, Charset.defaultCharset());
    }

    /**
     * Creates a writer for text to a Gzip file.
     *
     * @param file Gzip file to write to.
     * @param cs   Character set to use.
     * @return A Writer
     * @throws java.io.IOException Thrown if creation of the GZip Writer fails.
     */
    public static Writer createGzipWriter(File file, Charset cs) throws IOException {
        /*
         * The BufferedWriter buffers the input.
         * The OutputStreamWriter converts the input to bytes.
         * The GZIPOutputStream compresses the bytes.
         * The FileOutputStream writes bytes to a file.
         */
        return new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file)), cs));
    }

    /**
     * Creates a reader from a Gzip file to text.
     *
     * @param file Gzip file to read from.
     * @return A Reader.
     * @throws java.io.IOException Thrown if creation of the GZip Reader fails.
     */
    public static Reader createGzipReader(File file) throws IOException {
        return createGzipReader(file, Charset.defaultCharset());
    }

    /**
     * Creates a reader from a Gzip file to text.
     *
     * @param file Gzip file to read from.
     * @param cs   Character set to use.
     * @return A Reader.
     * @throws java.io.IOException Thrown if creation of the GZip Reader fails.
     */
    public static Reader createGzipReader(File file, Charset cs) throws IOException {
        /*
         * The BufferedReader buffers the input requests, reading a large chunk at a time and caching it.
         * The InputStreamReader converts the input bytes to characters.
         * The GZIPInputStream decompresses incoming input bytes from GZIP to raw bytes.
         * The FileInputStream reads raw bytes from a (gzipped) file.
         */
        return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)), cs));
    }

    /**
     * Prints a fixed-width decimal, similar to String.format(%width.precf, val),
     * but ensuring the resulting string is never longer than width. If the result
     * would end in a period (such as 14.), will leave off decimal. Will throw
     * an exception if the value cannot be formatted in the specified width.
     *
     * @param val   Value to print
     * @param width Width of field
     * @param prec  Number of decimal places
     * @return Formatted string
     * @throws java.lang.IllegalArgumentException if any.
     */
    public static String fwDec(double val, int width, int prec) throws IllegalArgumentException {
        if (width < 1 || prec < 0) {
            throw new IllegalArgumentException(" Must have width >= 1 and precision >= 0");
        }
        int w1 = width - 1;
        double maxVal = FastMath.pow(10.0, width);
        double minVal = maxVal / -10.0;

        if (val >= maxVal) {
            throw new IllegalArgumentException(String.format(" Value %f exceeded the maximum of %f enforced by width %d", val, maxVal, width));
        } else if (val <= minVal) {
            throw new IllegalArgumentException(String.format(" Value %f is less than the minumum of %f enforced by width %d", val, minVal, width));
        }

        String str = String.format("%" + width + "." + prec + "f", val);
        if (str.charAt(w1) == '.') {
            return " " + str.substring(0, w1);
        } else {
            return str.substring(0, width);
        }
    }

    /**
     * Prints a fixed-width decimal using String.format conventions, throwing an
     * error if the value cannot be formatted within that space.
     *
     * @param val   a double.
     * @param width a int.
     * @param prec  a int.
     * @return a {@link java.lang.String} object.
     * @throws java.lang.IllegalArgumentException If the length of String is greater than the width.
     */
    public static String fwFpDec(double val, int width, int prec) throws IllegalArgumentException {
        String str = String.format("%" + width + "." + prec + "f", val);
        if (str.length() > width) {
            throw new IllegalArgumentException(String.format(" Value %f cannot fit in width %d with precision %d", val, width, prec));
        } else {
            return str;
        }
    }

    /**
     * Prints a fixed-width decimal using String.format conventions, reducing the
     * value if necessary to fit within the width.
     *
     * @param val   a double.
     * @param width a int.
     * @param prec  a int.
     * @return a {@link java.lang.String} object.
     */
    public static String fwFpTrunc(double val, int width, int prec) {
        String str = String.format("%" + width + "." + prec + "f", val);
        if (str.length() > width) {
            StringBuilder sb;
            if (val < 0) {
                sb = new StringBuilder("-");
            } else {
                sb = new StringBuilder("9");
            }
            for (int i = 0; i < (width - prec - 2); i++) {
                sb.append("9");
            }
            sb.append(".");
            for (int i = 0; i < prec; i++) {
                sb.append("9");
            }
            str = sb.toString();
        }
        return str;
    }

    /**
     * <p>
     * padRight</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padRight(String s, int n) {
        return String.format("%-" + n + "s", s);
    }

    /**
     * <p>
     * padLeft</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padLeft(String s, int n) {
        return String.format("%" + n + "s", s);
    }

    /**
     * <p>
     * pdbForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String pdbForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".pdb";
    }

    /**
     * <p>
     * cifForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String cifForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".cif";
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     *
     * @param input  a {@link java.io.InputStream} object.
     * @param name   a {@link java.lang.String} object.
     * @param suffix a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public static String copyInputStreamToTmpFile(final InputStream input,
                                                  String name, final String suffix) throws IOException {
        File tmpFile = null;
        try {
            name = "ffx." + name + ".";
            tmpFile = File.createTempFile(name, suffix);
        } catch (Exception e) {
            System.out.println(" Could not extract a Force Field X library.");
            System.err.println(e.toString());
            System.exit(-1);
        }

        tmpFile.deleteOnExit();
        OutputStream output = null;
        try {
            output = new BufferedOutputStream(new FileOutputStream(tmpFile));
            byte[] buffer = new byte[8192];
            int size;
            while ((size = input.read(buffer)) != -1) {
                output.write(buffer, 0, size);
            }
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
            }
        }

        return tmpFile.toString();
    }

    /**
     * Finds consecutive subranges in an array of ints, and returns their mins
     * and maxes. This can include singletons.
     *
     * Example: [4, 5, 6, 1, 1, 2, 5, 6, 7] would become [4,6],[1,1],[1,2],[5,7]
     *
     * @param set Array of ints to split into consecutive subranges.
     * @return Consecutive subrange mins, maxes
     */
    public static List<int[]> consecutiveInts(int[] set) {
        if (set == null || set.length == 0) {
            return Collections.emptyList();
        }
        List<int[]> allRanges = new ArrayList<>();

        int rangeStart = set[0];
        int rangeEnd = rangeStart;
        for (int i = 1; i < set.length; i++) {
            if (set[i] == rangeEnd + 1) {
                rangeEnd = set[i];
            } else {
                allRanges.add(new int[]{rangeStart, rangeEnd});
                rangeStart = set[i];
                rangeEnd = rangeStart;
            }
        }
        allRanges.add(new int[]{rangeStart, rangeEnd});
        return allRanges;
    }

}
