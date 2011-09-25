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
package ffx.utilities;

/**
 * Java port of the hy36encode() and hy36decode() functions in the
 * hybrid_36.py Python prototype/reference implementation.
 *
 * @author Michael J. Schnieders
 *         Derived from code by:
 *         Ralf W. Grosse-Kunstleve, Vincent B. Chen, Jeff J. Headd, Sep 2007.
 * @see <a href="http://cci.lbl.gov/hybrid_36/">LBL Hybrid36 Reference</a>
 * @since 1.0
 * @version $Id: $
 */
public class Hybrid36 {

    private static String digitsBase10 = "0123456789";
    private static String digitsUpper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private static String digitsLower = "0123456789abcdefghijklmnopqrstuvwxyz";
    private static boolean firstCall = true;
    private static int[] digitsValuesUpper = new int[128];
    private static int[] digitsValuesLower = new int[128];
    private static String valueOutOfRange = "value out of range.";
    private static String invalidNumberLiteral = "invalid number literal.";
    private static String unsupportedWidth = "unsupported width.";

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
        boolean have_minus = false;
        boolean have_non_blank = false;
        int value = 0;
        for (int i = 0; i < s.length(); i++) {
            char si = s.charAt(i);
            if (si < 0 || si > 127) {
                throw new Error(invalidNumberLiteral);
            }
            if (si == ' ') {
                if (!have_non_blank) {
                    continue;
                }
                value *= digits_size;
            } else if (si == '-') {
                if (have_non_blank) {
                    throw new Error(invalidNumberLiteral);
                }
                have_non_blank = true;
                have_minus = true;
                continue;
            } else {
                have_non_blank = true;
                int dv = digits_values[si];
                if (dv < 0 || dv >= digits_size) {
                    throw new Error(invalidNumberLiteral);
                }
                value *= digits_size;
                value += dv;
            }
        }
        if (have_minus) {
            value = -value;
        }
        return value;
    }

    /**
     * hybrid-36 encoder: converts integer value to string result
     *
     *    width: must be 4 (e.g. for residue sequence numbers)
     *    or 5 (e.g. for atom serial numbers)
     *
     *    value: integer value to be converted
     *
     *    return value: String of size width
     *
     * @param width a int.
     * @param value a int.
     * @return a {@link java.lang.String} object.
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
     * hybrid-36 decoder: converts string s to integer result
     *
     *    width: must be 4 (e.g. for residue sequence numbers)
     *    or 5 (e.g. for atom serial numbers)
     *
     *    s: string to be converted
     *
     *    return value: conversion result
     *
     * @param width a int.
     * @param s a {@link java.lang.String} object.
     * @return a int.
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

    private static void checkString(String result, String expected) {
        if (!result.equals(expected)) {
            System.out.println("ERROR: \"" + result + "\" != \"" + expected + "\"");
        }
    }
    private static void checkInt(int result, int expected) {
        if (result != expected) {
            System.out.println("ERROR: " + result + " != " + expected);
        }
    }
    private static void recycle4(int value, String encoded) {
        String s = encode(4, value);
        checkString(s, encoded);
        int d = decode(4, s);
        checkInt(d, value);
    }
    private static void recycle5(int value, String encoded) {
        String s = encode(5, value);
        checkString(s, encoded);
        int d = decode(5, s);
        checkInt(d, value);
    }
    private static void checkEncodeException(int width, int value, String expected_msg) {
        String msg = "";
        try {
            encode(width, value);
        } catch (Error e) {
            msg = e.toString();
        }
        checkString(msg, "java.lang.Error: " + expected_msg);
    }
    private static void checkDecodeException(int width, String s, String expected_msg) {
        String msg = "";
        try {
            decode(width, s);
        } catch (Error e) {
            msg = e.toString();
        }
        checkString(msg, "java.lang.Error: " + expected_msg);
    }

    private static int random_seed = 13;

    private static int kernighan_and_ritchie_rand() {
        random_seed = random_seed * 1103515245 + 12345;
        int result = (random_seed / 65536) % 32768;
        if (result < 0) {
            result += 32768;
        }
        return result;
    }

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        checkInt(decode(4, "    "), 0);
        checkInt(decode(4, "  -0"), 0);
        recycle4(-999, "-999");
        recycle4(-78, " -78");
        recycle4(-6, "  -6");
        recycle4(0, "   0");
        recycle4(9999, "9999");
        recycle4(10000, "A000");
        recycle4(10001, "A001");
        recycle4(10002, "A002");
        recycle4(10003, "A003");
        recycle4(10004, "A004");
        recycle4(10005, "A005");
        recycle4(10006, "A006");
        recycle4(10007, "A007");
        recycle4(10008, "A008");
        recycle4(10009, "A009");
        recycle4(10010, "A00A");
        recycle4(10011, "A00B");
        recycle4(10012, "A00C");
        recycle4(10013, "A00D");
        recycle4(10014, "A00E");
        recycle4(10015, "A00F");
        recycle4(10016, "A00G");
        recycle4(10017, "A00H");
        recycle4(10018, "A00I");
        recycle4(10019, "A00J");
        recycle4(10020, "A00K");
        recycle4(10021, "A00L");
        recycle4(10022, "A00M");
        recycle4(10023, "A00N");
        recycle4(10024, "A00O");
        recycle4(10025, "A00P");
        recycle4(10026, "A00Q");
        recycle4(10027, "A00R");
        recycle4(10028, "A00S");
        recycle4(10029, "A00T");
        recycle4(10030, "A00U");
        recycle4(10031, "A00V");
        recycle4(10032, "A00W");
        recycle4(10033, "A00X");
        recycle4(10034, "A00Y");
        recycle4(10035, "A00Z");
        recycle4(10036, "A010");
        recycle4(10046, "A01A");
        recycle4(10071, "A01Z");
        recycle4(10072, "A020");
        recycle4(10000 + 36 * 36 - 1, "A0ZZ");
        recycle4(10000 + 36 * 36, "A100");
        recycle4(10000 + 36 * 36 * 36 - 1, "AZZZ");
        recycle4(10000 + 36 * 36 * 36, "B000");
        recycle4(10000 + 26 * 36 * 36 * 36 - 1, "ZZZZ");
        recycle4(10000 + 26 * 36 * 36 * 36, "a000");
        recycle4(10000 + 26 * 36 * 36 * 36 + 35, "a00z");
        recycle4(10000 + 26 * 36 * 36 * 36 + 36, "a010");
        recycle4(10000 + 26 * 36 * 36 * 36 + 36 * 36 - 1, "a0zz");
        recycle4(10000 + 26 * 36 * 36 * 36 + 36 * 36, "a100");
        recycle4(10000 + 26 * 36 * 36 * 36 + 36 * 36 * 36 - 1, "azzz");
        recycle4(10000 + 26 * 36 * 36 * 36 + 36 * 36 * 36, "b000");
        recycle4(10000 + 2 * 26 * 36 * 36 * 36 - 1, "zzzz");
        //
        checkInt(decode(5, "     "), 0);
        checkInt(decode(5, "   -0"), 0);
        recycle5(-9999, "-9999");
        recycle5(-123, " -123");
        recycle5(-45, "  -45");
        recycle5(-6, "   -6");
        recycle5(0, "    0");
        recycle5(12, "   12");
        recycle5(345, "  345");
        recycle5(6789, " 6789");
        recycle5(99999, "99999");
        recycle5(100000, "A0000");
        recycle5(100010, "A000A");
        recycle5(100035, "A000Z");
        recycle5(100036, "A0010");
        recycle5(100046, "A001A");
        recycle5(100071, "A001Z");
        recycle5(100072, "A0020");
        recycle5(100000 + 36 * 36 - 1, "A00ZZ");
        recycle5(100000 + 36 * 36, "A0100");
        recycle5(100000 + 36 * 36 * 36 - 1, "A0ZZZ");
        recycle5(100000 + 36 * 36 * 36, "A1000");
        recycle5(100000 + 36 * 36 * 36 * 36 - 1, "AZZZZ");
        recycle5(100000 + 36 * 36 * 36 * 36, "B0000");
        recycle5(100000 + 2 * 36 * 36 * 36 * 36, "C0000");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 - 1, "ZZZZZ");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36, "a0000");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 - 1, "a000z");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36, "a0010");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 - 1, "a00zz");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36, "a0100");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 - 1, "a0zzz");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36, "a1000");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36 - 1, "azzzz");
        recycle5(100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36, "b0000");
        recycle5(100000 + 2 * 26 * 36 * 36 * 36 * 36 - 1, "zzzzz");
        //
        checkEncodeException(4, -1000, "value out of range.");
        checkEncodeException(4, 2436112, "value out of range.");
        checkEncodeException(5, -10000, "value out of range.");
        checkEncodeException(5, 87440032, "value out of range.");
        //
        checkDecodeException(4, "", "invalid number literal.");
        checkDecodeException(4, "    0", "invalid number literal.");
        checkDecodeException(4, " abc", "invalid number literal.");
        checkDecodeException(4, "abc-", "invalid number literal.");
        checkDecodeException(4, "A=BC", "invalid number literal.");
        checkDecodeException(4, "40a0", "invalid number literal.");
        checkDecodeException(4, "40A0", "invalid number literal.");
        checkDecodeException(5, "", "invalid number literal.");
        checkDecodeException(5, "     0", "invalid number literal.");
        checkDecodeException(5, " abcd", "invalid number literal.");
        checkDecodeException(5, "ABCD-", "invalid number literal.");
        checkDecodeException(5, "a=bcd", "invalid number literal.");
        checkDecodeException(5, "410b0", "invalid number literal.");
        checkDecodeException(5, "410B0", "invalid number literal.");
        //
        checkEncodeException(3, 0, "unsupported width.");
        checkEncodeException(6, 0, "unsupported width.");
        checkDecodeException(3, "AAA", "unsupported width.");
        checkDecodeException(6, "zzzzzz", "unsupported width.");
        //
        int value = -9999;
        while (value < 100000 + 2 * 26 * 36 * 36 * 36 * 36) {
            checkInt(decode(5, encode(5, value)), value);
            value += kernighan_and_ritchie_rand() % 10000;
        }
        System.out.println("OK");
    }
}
