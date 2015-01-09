/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
public class Hybrid36Test {

    /**
     * Test of encode method, of class Hybrid36.
     */
    @Test
    public void testEncode() {
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
        checkEncodeException(3, 0, "unsupported width.");
        checkEncodeException(6, 0, "unsupported width.");
        //
        int value = -9999;
        while (value < 100000 + 2 * 26 * 36 * 36 * 36 * 36) {
            checkInt(Hybrid36.decode(5, Hybrid36.encode(5, value)), value);
            value += kernighan_and_ritchie_rand() % 10000;
        }
    }

    /**
     * Test of decode method, of class Hybrid36.
     */
    @Test
    public void testDecode() {
        checkInt(Hybrid36.decode(4, "    "), 0);
        checkInt(Hybrid36.decode(4, "  -0"), 0);
        checkInt(Hybrid36.decode(5, "     "), 0);
        checkInt(Hybrid36.decode(5, "   -0"), 0);
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
        checkDecodeException(3, "AAA", "unsupported width.");
        checkDecodeException(6, "zzzzzz", "unsupported width.");
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

    private static void checkEncodeException(int width, int value, String expected_msg) {
        String msg = "";
        try {
            Hybrid36.encode(width, value);
        } catch (Error e) {
            msg = e.toString();
        }
        checkString(msg, "java.lang.Error: " + expected_msg);
    }

    private static void checkDecodeException(int width, String s, String expected_msg) {
        String msg = "";
        try {
            Hybrid36.decode(width, s);
        } catch (Error e) {
            msg = e.toString();
        }
        checkString(msg, "java.lang.Error: " + expected_msg);
    }

    private static void checkString(String result, String expected) {
        assertEquals(" checkString", result, expected);
        if (!result.equals(expected)) {
            System.out.println("\"" + result + "\" != \"" + expected + "\"");
        }
    }

    private static void checkInt(int result, int expected) {
        assertEquals(" checkInt", result, expected);
        if (result != expected) {
            System.out.println("" + result + " != " + expected);
        }
    }

    private static void recycle4(int value, String encoded) {
        String s = Hybrid36.encode(4, value);
        checkString(s, encoded);
        int d = Hybrid36.decode(4, s);
        checkInt(d, value);
    }

    private static void recycle5(int value, String encoded) {
        String s = Hybrid36.encode(5, value);
        checkString(s, encoded);
        int d = Hybrid36.decode(5, s);
        checkInt(d, value);
    }
}
