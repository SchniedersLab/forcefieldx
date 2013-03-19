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
package ffx.utilities;

import java.util.logging.Logger;

/**
 * Java port of the hy36encode() and hy36decode() functions in the hybrid_36.py
 * Python prototype/reference implementation.
 *
 * @author Michael J. Schnieders <br> Derived from code by: Ralf W.
 * Grosse-Kunstleve, Vincent B. Chen, Jeff J. Headd, Sep 2007.
 * @see <a href="http://cci.lbl.gov/hybrid_36/">LBL Hybrid36 Reference</a>
 * @since 1.0
 *
 */
public class Hybrid36 {

    private static final Logger logger = Logger.getLogger(Hybrid36.class.getName());
    private static String digitsBase10 = "0123456789";
    private static String digitsUpper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private static String digitsLower = "0123456789abcdefghijklmnopqrstuvwxyz";
    private static boolean firstCall = true;
    private static String valueOutOfRange = "value out of range.";
    private static String invalidNumberLiteral = "invalid number literal.";
    private static String unsupportedWidth = "unsupported width.";
    private static int[] digitsValuesUpper = new int[128];
    private static int[] digitsValuesLower = new int[128];

    private static String encodePure(String digits, int width, int value) {
        boolean neg = false;
        if (value < 0) {
            neg = true;
            value = -value;
        }
        String buf = "";
        while (true) {
            int rest = value / digits.length();
            buf += digits.charAt(value - rest * digits.length());
            if (rest == 0) {
                break;
            }
            value = rest;
        }
        if (neg) {
            buf += '-';
        }
        String result = "";
        for (int i = buf.length(); i < width; i++) {
            result += " ";
        }
        for (int i = buf.length() - 1; i >= 0; i--) {
            result += buf.charAt(i);
        }
        return result;
    }

    private static int decodePure(int[] digits_values, int digits_size,
            String s) {
        boolean haveMinus = false;
        boolean haveNonBlank = false;
        int value = 0;
        for (int i = 0; i < s.length(); i++) {
            char si = s.charAt(i);
            if (si < 0 || si > 127) {
                throw new Error(invalidNumberLiteral);
            }
            if (si == ' ') {
                if (!haveNonBlank) {
                    continue;
                }
                value *= digits_size;
            } else if (si == '-') {
                if (haveNonBlank) {
                    throw new Error(invalidNumberLiteral);
                }
                haveNonBlank = true;
                haveMinus = true;
                continue;
            } else {
                haveNonBlank = true;
                int dv = digits_values[si];
                if (dv < 0 || dv >= digits_size) {
                    throw new Error(invalidNumberLiteral);
                }
                value *= digits_size;
                value += dv;
            }
        }
        if (haveMinus) {
            value = -value;
        }
        return value;
    }

    /**
     * Hybrid-36 encoder: converts integer value to string result.
     *
     * @param width must be 4 (e.g. for residue sequence numbers) or 5 (e.g. for
     * atom serial numbers).
     * @param value the integer value to be converted.
     * @return a {@link java.lang.String} String of size width.
     */
    public static String encode(int width, int value) {
        int i = value;
        if (width == 4) {
            if (i >= -999) {
                if (i < 10000) {
                    return encodePure(digitsBase10, 4, i);
                }
                i -= 10000;
                if (i < 1213056 /* 26*36**3 */) {
                    i += 466560 /* 10*36**3 */;
                    return encodePure(digitsUpper, 0, i);
                }
                i -= 1213056;
                if (i < 1213056) {
                    i += 466560;
                    return encodePure(digitsLower, 0, i);
                }
            }
        } else if (width == 5) {
            if (i >= -9999) {
                if (i < 100000) {
                    return encodePure(digitsBase10, 5, i);
                }
                i -= 100000;
                if (i < 43670016 /* 26*36**4 */) {
                    i += 16796160 /* 10*36**4 */;
                    return encodePure(digitsUpper, 0, i);
                }
                i -= 43670016;
                if (i < 43670016) {
                    i += 16796160;
                    return encodePure(digitsLower, 0, i);
                }
            }
        } else {
            throw new Error(unsupportedWidth);
        }
        throw new Error(valueOutOfRange);
    }

    /**
     * Hybrid-36 decoder: converts string s to integer result.
     *
     * @param width must be 4 (e.g. for residue sequence numbers) or 5 (e.g. for
     * atom serial numbers)
     * @param s the {@link java.lang.String} to be converted.
     * @return the conversion result.
     */
    public static int decode(int width, String s) {
        String ie_range = "internal error hy36.decode: integer value out of range.";
        if (firstCall) {
            firstCall = false;
            for (int i = 0; i < 128; i++) {
                digitsValuesUpper[i] = -1;
            }
            for (int i = 0; i < 128; i++) {
                digitsValuesLower[i] = -1;
            }
            for (int i = 0; i < 36; i++) {
                int di = (int) digitsUpper.charAt(i);
                if (di < 0 || di > 127) {
                    throw new Error(ie_range);
                }
                digitsValuesUpper[di] = i;
            }
            for (int i = 0; i < 36; i++) {
                int di = (int) digitsLower.charAt(i);
                if (di < 0 || di > 127) {
                    throw new Error(ie_range);
                }
                digitsValuesLower[di] = i;
            }
        }
        if (s.length() == width) {
            int di = (int) s.charAt(0);
            if (di >= 0 && di <= 127) {
                if (digitsValuesUpper[di] >= 10) {
                    int result = decodePure(digitsValuesUpper, 36, s);
                    /* result - 10*36**(width-1) + 10**width */
                    if (width == 4) {
                        result -= 456560;
                    } else if (width == 5) {
                        result -= 16696160;
                    } else {
                        throw new Error(unsupportedWidth);
                    }
                    return result;
                } else if (digitsValuesLower[di] >= 10) {
                    int result = decodePure(digitsValuesLower, 36, s);
                    /* result + 16*36**(width-1) + 10**width */
                    if (width == 4) {
                        result += 756496;
                    } else if (width == 5) {
                        result += 26973856;
                    } else {
                        throw new Error(unsupportedWidth);
                    }
                    return result;
                } else {
                    int result = decodePure(digitsValuesUpper, 10, s);
                    if (!(width == 4 || width == 5)) {
                        throw new Error(unsupportedWidth);
                    }
                    return result;
                }
            }
        }
        throw new Error(invalidNumberLiteral);
    }
}
