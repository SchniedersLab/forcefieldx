//******************************************************************************
//
// File:    ObjectOp.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.ObjectOp
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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
package edu.rit.pj.reduction;

/**
 * Class ObjectOp is the abstract base class for a binary operation on object
 * values, used to do reduction in a parallel program.
 * <P>
 * When classes in the Parallel Java Library call the <TT>op(x,y)</TT> method
 * during a reduction operation, the <TT>x</TT> argument is the value of the
 * reduction variable, and the <TT>y</TT> argument is the value to be combined
 * with the reduction variable. The value returned by the <TT>op(x,y)</TT>
 * method is stored back into the reduction variable.
 * <P>
 * The <TT>op(x,y)</TT> method in any subclass of class ObjectOp must obey the
 * following requirements, which are assumed by classes in the Parallel Java
 * Library:
 * <UL>
 * <LI>
 * If as a result of the binary operation the state of the reduction variable
 * will not change, the <TT>op(x,y)</TT> method must return <TT>x</TT>.
 * <P>
 * <LI>
 * If as a result of the binary operation the state of the reduction variable
 * will change, the <TT>op(x,y)</TT> method must return a newly created object
 * containing the desired state.
 * <P>
 * <LI>
 * The <TT>op(x,y)</TT> method must neither change the state of <TT>x</TT> nor
 * change the state of <TT>y</TT>; that is, the <TT>op(x,y)</TT> method must
 * have no side effects.
 * </UL>
 *
 * @param <T> Object data type.
 *
 * @author Alan Kaminsky
 * @version 30-Mar-2008
 */
public abstract class ObjectOp<T>
        extends Op {

// Hidden constructors.
    /**
     * Construct a new object binary operation.
     */
    protected ObjectOp() {
        super();
    }

// Exported operations.
    /**
     * Perform this binary operation.
     *
     * @param x First argument.
     * @param y Second argument.
     *
     * @return (<TT>x</TT> <I>op</I> <TT>y</TT>), where <I>op</I> stands for
     * this binary operation.
     */
    public abstract T op(T x,
            T y);

}
