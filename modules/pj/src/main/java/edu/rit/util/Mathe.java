//******************************************************************************
//
// File:    Mathe.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Mathe
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.util;

/**
 * Class Mathe provides useful mathematical operations. The class is named
 * "Mathe" so the compiler won't confuse it with class
 * org.apache.commons.math3.util.FastMath.
 *
 * @author Alan Kaminsky
 * @version 11-Feb-2010
 */
public class Mathe {

// Prevent construction.
    private Mathe() {
    }

// Exported operations.
    /**
     * Compute the integer square root of the integer <TT>x</TT>. The value
     * floor(<TT>x</TT><SUP>1/2</SUP>) is returned. The answer is calculated
     * using an exact integer algorithm taken from:
     * <UL>
     * <LI>
     * J. Crenshaw. Integer square roots.
     * <A HREF="http://www.embedded.com/98/9802fe2.htm"
     * TARGET="_top">http://www.embedded.com/98/9802fe2.htm</A>
     * </UL>
     *
     * @param x Input.
     * @return Floor(<TT>x</TT><SUP>1/2</SUP>).
     * @exception ArithmeticException (unchecked exception) Thrown if <TT>x</TT>
     * &lt; 0.
     */
    public static int sqrt(int x) {
        if (x < 0) {
            throw new ArithmeticException("Mathe.sqrt(): x < 0");
        }
        int rem = 0;
        int root = 0;
        for (int i = 0; i < 16; ++i) {
            root <<= 1;
            rem = (rem << 2) | (x >>> 30);
            x <<= 2;
            ++root;
            if (root <= rem) {
                rem -= root;
                ++root;
            } else {
                --root;
            }
        }
        return root >>> 1;
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		for (int x = 0; x <= 100000000; ++ x)
//			{
//			if ((x % 10000000) == 0) System.out.printf ("%d%n", x);
//			int y = Mathe.sqrt(x);
//			if (y*y > x || x >= (y + 1)*(y + 1))
//				{
//				System.out.printf ("x = %d, sqrt(x) = %d, fail%n", x, y);
//				}
//			}
//		}
}
