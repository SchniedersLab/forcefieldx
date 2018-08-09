//******************************************************************************
//
// File:    Spinner.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.Spinner
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj;

/**
 * Class Spinner provides an object used to implement a spin-wait.
 * <P>
 * The pattern for coding a spin-wait is:
 * <PRE>
 *    if (&lt;condition&gt;)
 *        {
 *        Spinner spinner = new Spinner();
 *        while (&lt;condition&gt;) spinner.spin();
 *        }
 * </PRE> This will wait as long as the <TT>&lt;condition&gt;</TT> is true.
 *
 * @author Alan Kaminsky
 * @version 20-Dec-2007
 */
class Spinner {

// Hidden constants.
    // In a spin-wait, the maximum number of iterations to spin before yielding
    // the CPU.
    static final int MAX_COUNT = 10000;

// Hidden data members.
    // Spin counter.
    volatile int count;

    // 128 bytes of extra padding to avert cache interference.
    private long p0, p1, p2, p3, p4, p5, p6, p7;
    private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.
    /**
     * Construct a new spinner.
     */
    public Spinner() {
    }

// Exported operations.
    /**
     * Spin this spinner. If enough consecutive <TT>spin()</TT> calls occur, the
     * calling thread yields the CPU.
     */
    public void spin() {
        if (count++ > MAX_COUNT) {
            Thread.yield();
            count = 0;
        }
    }

}
