//******************************************************************************
//
// File:    ParallelForLoop.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelForLoop
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj;

/**
 * Class ParallelForLoop is the abstract base class for a parallel for loop that
 * is executed inside a {@linkplain ParallelRegion}. There are four variations
 * of a parallel for loop, depending on the loop index data type and on whether
 * the loop stride is implicit or explicit. There is a subclass for each
 * variation; see the subclasses for further information. The subclasses are:
 * <UL>
 * <LI>
 * Class {@linkplain IntegerForLoop} -- loop index is type <TT>int</TT>, loop
 * stride is implicit (+1).
 * <BR>&nbsp;
 * <LI>
 * Class {@linkplain IntegerStrideForLoop} -- loop index is type <TT>int</TT>,
 * loop stride is explicitly specified.
 * <BR>&nbsp;
 * <LI>
 * Class {@linkplain LongForLoop} -- loop index is type <TT>long</TT>, loop
 * stride is implicit (+1).
 * <BR>&nbsp;
 * <LI>
 * Class {@linkplain LongStrideForLoop} -- loop index is type <TT>long</TT>,
 * loop stride is explicitly specified.
 * <BR>&nbsp;
 * </UL>
 *
 * @author Alan Kaminsky
 * @version 31-May-2007
 */
public abstract class ParallelForLoop
        extends ParallelConstruct {

// Exported constructors.
    /**
     * Construct a new parallel for loop.
     */
    public ParallelForLoop() {
        super();
    }

}
