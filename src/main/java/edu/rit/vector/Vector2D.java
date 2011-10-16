//******************************************************************************
//
// File:    Vector2D.java
// Package: edu.rit.vector
// Unit:    Class edu.rit.vector.Vector2D
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

package edu.rit.vector;

import edu.rit.mp.DoubleBuf;

import edu.rit.util.Range;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Vector2D provides a two-dimensional mathematical vector of type
 * <TT>double</TT>. The vector's components are stored in the fields <TT>x</TT>
 * and <TT>y</TT>.
 * <P>
 * Call one of the following static factory methods to create a buffer for
 * sending and receiving Vector2D items in a message. To reduce the number of
 * bytes in the message, thereby improving performance, each of these methods
 * returns a buffer which treats a vector as a pair of doubles (rather than an
 * object). This buffer does not support reduction. The buffer's
 * <TT>getReductionBuf()</TT> method throws an UnsupportedOperationException.
 * <UL>
 * <LI><TT>doubleBuffer (Vector2D[])</TT>
 * <LI><TT>doubleSliceBuffer (Vector2D[], Range)</TT>
 * <LI><TT>doubleSliceBuffers (Vector2D[], Range[])</TT>
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 04-Feb-2008
 */
public class Vector2D
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = -4617340700183457338L;

// Exported data members.

	/**
	 * X component.
	 */
	public double x;

	/**
	 * Y component.
	 */
	public double y;

// Exported constructors.

	/**
	 * Construct a new vector. The X and Y components are initialized to 0.
	 */
	public Vector2D()
		{
		}

	/**
	 * Construct a new vector with the given X and Y components.
	 *
	 * @param  x  Initial X component.
	 * @param  y  Initial Y component.
	 */
	public Vector2D
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		}

	/**
	 * Construct a new vector that is a copy of the given vector.
	 *
	 * @param  theVector  Vector to copy.
	 */
	public Vector2D
		(Vector2D theVector)
		{
		this.x = theVector.x;
		this.y = theVector.y;
		}

// Exported operations.

	/**
	 * Clear this vector. The X and Y components are set to 0.
	 *
	 * @return  This vector, set to (0,0).
	 */
	public Vector2D clear()
		{
		this.x = 0.0;
		this.y = 0.0;
		return this;
		}

	/**
	 * Assign the given X and Y components to this vector.
	 *
	 * @param  x  X component.
	 * @param  y  Y component.
	 *
	 * @return  This vector, set to (<TT>x</TT>,<TT>y</TT>).
	 */
	public Vector2D assign
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		return this;
		}

	/**
	 * Assign (copy) the given vector to this vector.
	 *
	 * @param  theVector  Vector to copy.
	 *
	 * @return  This vector, set to <TT>theVector</TT>.
	 */
	public Vector2D assign
		(Vector2D theVector)
		{
		this.x = theVector.x;
		this.y = theVector.y;
		return this;
		}

	/**
	 * Add the given vector to this vector.
	 *
	 * @param  theVector  Vector to add.
	 *
	 * @return  This vector, set to <TT>this</TT> + <TT>theVector</TT>.
	 */
	public Vector2D add
		(Vector2D theVector)
		{
		this.x += theVector.x;
		this.y += theVector.y;
		return this;
		}

	/**
	 * Subtract the given vector from this vector.
	 *
	 * @param  theVector  Vector to subtract.
	 *
	 * @return  This vector, set to <TT>this</TT> - <TT>theVector</TT>.
	 */
	public Vector2D sub
		(Vector2D theVector)
		{
		this.x -= theVector.x;
		this.y -= theVector.y;
		return this;
		}

	/**
	 * Multiply this vector by the given scalar.
	 *
	 * @param  a  Scalar.
	 *
	 * @return  This vector, set to <TT>this</TT> * <TT>a</TT>.
	 */
	public Vector2D mul
		(double a)
		{
		this.x *= a;
		this.y *= a;
		return this;
		}

	/**
	 * Divide this vector by the given scalar.
	 *
	 * @param  a  Scalar.
	 *
	 * @return  This vector, set to <TT>this</TT> / <TT>a</TT>.
	 */
	public Vector2D div
		(double a)
		{
		this.x /= a;
		this.y /= a;
		return this;
		}

	/**
	 * Negate this vector.
	 *
	 * @return  This vector, set to -<TT>this</TT>.
	 */
	public Vector2D neg()
		{
		this.x = -this.x;
		this.y = -this.y;
		return this;
		}

	/**
	 * Rotate this vector 90 degrees counterclockwise.
	 *
	 * @return  This vector, rotated.
	 */
	public Vector2D rotate90()
		{
		double tmp = this.x;
		this.x = -this.y;
		this.y = tmp;
		return this;
		}

	/**
	 * Rotate this vector 180 degrees.
	 *
	 * @return  This vector, rotated.
	 */
	public Vector2D rotate180()
		{
		this.x = -this.x;
		this.y = -this.y;
		return this;
		}

	/**
	 * Rotate this vector 270 degrees counterclockwise (90 degrees clockwise).
	 *
	 * @return  This vector, rotated.
	 */
	public Vector2D rotate270()
		{
		double tmp = this.x;
		this.x = this.y;
		this.y = -tmp;
		return this;
		}

	/**
	 * Normalize this vector. This vector is set to a unit vector pointing in
	 * the same direction.
	 *
	 * @return  This vector, set to <TT>this</TT> / <TT>this.mag()</TT>.
	 */
	public Vector2D norm()
		{
		return this.div (this.mag());
		}

	/**
	 * Determine the dot product of this vector and the given vector.
	 *
	 * @param  theVector  Vector.
	 *
	 * @return  Dot product.
	 */
	public double dot
		(Vector2D theVector)
		{
		return this.x * theVector.x + this.y * theVector.y;
		}

	/**
	 * Determine the magnitude of this vector. The magnitude of the vector
	 * (<I>x</I>,<I>y</I>) is
	 * (<I>x</I><SUP>2</SUP>&nbsp;+&nbsp;<I>y</I><SUP>2</SUP>)<SUP>1/2</SUP>.
	 *
	 * @return  Magnitude.
	 */
	public double mag()
		{
		return Math.sqrt (x*x + y*y);
		}

	/**
	 * Determine the squared magnitude of this vector. The squared magnitude of
	 * the vector (<I>x</I>,<I>y</I>) is
	 * (<I>x</I><SUP>2</SUP>&nbsp;+&nbsp;<I>y</I><SUP>2</SUP>).
	 *
	 * @return  Squared magnitude.
	 */
	public double sqrMag()
		{
		return x*x + y*y;
		}

	/**
	 * Determine the argument of this vector. The argument of the vector
	 * (<I>x</I>,<I>y</I>) is tan<SUP>-1</SUP> (<I>y</I>/<I>x</I>).
	 *
	 * @return  Argument.
	 */
	public double arg()
		{
		return Math.atan2 (y, x);
		}

	/**
	 * Determine the distance (magnitude of the difference) between this vector
	 * and the given vector.
	 *
	 * @param  theVector  Vector.
	 *
	 * @return  Distance between this vector and <TT>theVector</TT>.
	 */
	public double dist
		(Vector2D theVector)
		{
		double dx = this.x - theVector.x;
		double dy = this.y - theVector.y;
		return Math.sqrt (dx*dx + dy*dy);
		}

	/**
	 * Determine the squared distance (squared magnitude of the difference)
	 * between this vector and the given vector.
	 *
	 * @param  theVector  Vector.
	 *
	 * @return  Squared distance between this vector and <TT>theVector</TT>.
	 */
	public double sqrDist
		(Vector2D theVector)
		{
		double dx = this.x - theVector.x;
		double dy = this.y - theVector.y;
		return dx*dx + dy*dy;
		}

	/**
	 * Write this vector to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeDouble (x);
		out.writeDouble (y);
		}

	/**
	 * Read this vector from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		x = in.readDouble();
		y = in.readDouble();
		}

	/**
	 * Returns a string version of this vector. The string is of the form
	 * <TT>"(</TT><I>x</I><TT>,</TT><I>y</I><TT>)"</TT>.
	 *
	 * @return  String version.
	 */
	public String toString()
		{
		return "(" + x + "," + y + ")";
		}

	/**
	 * Create a buffer for the entire given vector array. The returned buffer
	 * encompasses all the elements in <TT>theArray</TT>.
 	 * <P>
	 * The returned buffer treats a vector as a pair of doubles (rather than an
	 * object). This buffer does not support reduction. The buffer's
	 * <TT>getReductionBuf()</TT> method throws an
	 * UnsupportedOperationException.
	 *
	 * @param  theArray  Array.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public static DoubleBuf doubleBuffer
		(Vector2D[] theArray)
		{
		return new Vector2DArrayDoubleBuf_1
			(theArray, new Range (0, theArray.length-1));
		}

	/**
	 * Create a buffer for one slice of the given vector array. The returned
	 * buffer encompasses <TT>theRange</TT> of elements in <TT>theArray</TT>.
	 * The range's stride may be 1 or greater than 1.
 	 * <P>
	 * The returned buffer treats a vector as a pair of doubles (rather than an
	 * object). This buffer does not support reduction. The buffer's
	 * <TT>getReductionBuf()</TT> method throws an
	 * UnsupportedOperationException.
	 *
	 * @param  theArray  Array.
	 * @param  theRange  Range of elements to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> does not include
	 *     all the indexes in <TT>theRange</TT>.
	 */
	public static DoubleBuf doubleSliceBuffer
		(Vector2D[] theArray,
		 Range theRange)
		{
		if (theRange.stride() == 1)
			{
			return new Vector2DArrayDoubleBuf_1 (theArray, theRange);
			}
		else
			{
			return new Vector2DArrayDoubleBuf (theArray, theRange);
			}
		}

	/**
	 * Create an array of buffers for multiple slices of the given vector
	 * array. The returned buffer array has the same length as
	 * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
	 * encompasses the elements of <TT>theArray</TT> specified by
	 * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
 	 * <P>
	 * Each buffer in the returned buffer array treats a vector as a pair of
	 * doubles (rather than an object). This buffer does not support reduction.
	 * The buffer's <TT>getReductionBuf()</TT> method throws an
	 * UnsupportedOperationException.
	 *
	 * @param  theArray   Array.
	 * @param  theRanges  Array of ranges of elements to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT>'s allocation does
	 *     not include any element of <TT>theRanges</TT>.
	 */
	public static DoubleBuf[] doubleSliceBuffers
		(Vector2D[] theArray,
		 Range[] theRanges)
		{
		int n = theRanges.length;
		DoubleBuf[] result = new DoubleBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = doubleSliceBuffer (theArray, theRanges[i]);
			}
		return result;
		}

	}
