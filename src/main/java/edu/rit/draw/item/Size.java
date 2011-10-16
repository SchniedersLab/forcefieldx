//******************************************************************************
//
// File:    Size.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Size
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

package edu.rit.draw.item;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Size provides a size (<I>width,height</I>) for a {@linkplain
 * DrawingItem}.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public class Size
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 2935902495697224710L;

	// The width and height.
	double width;
	double height;

// Exported constructors.

	/**
	 * Construct a new size with width = height = 0.
	 */
	public Size()
		{
		}

	/**
	 * Construct a new size with the given width and height.
	 *
	 * @param  width   Width.
	 * @param  height  Height.
	 */
	public Size
		(double width,
		 double height)
		{
		this.width = width;
		this.height = height;
		}

	/**
	 * Construct a new size with the same width and height as the given size.
	 *
	 * @param  theSize  Size to copy.
	 */
	public Size
		(Size theSize)
		{
		this.width = theSize.width;
		this.height = theSize.height;
		}

// Exported operations.

	/**
	 * Returns this size's width.
	 */
	public double width()
		{
		return this.width;
		}

	/**
	 * Returns this size's height.
	 */
	public double height()
		{
		return this.height;
		}

	/**
	 * Returns a new size with its width and height set to the ceiling of this
	 * size's width and height, respectively.
	 */
	public Size ceil()
		{
		return new Size (Math.ceil (this.width), Math.ceil (this.height));
		}

	/**
	 * Returns a new size equal to this size increased by the given amount in
	 * both width and height.
	 *
	 * @param  incr  Increment.
	 */
	public Size add
		(double incr)
		{
		return new Size (this.width + incr, this.height + incr);
		}

	/**
	 * Returns a new size equal to this size increased by the given amounts in
	 * width and height.
	 *
	 * @param  wincr  Width increment.
	 * @param  hincr  Height increment.
	 */
	public Size add
		(double wincr,
		 double hincr)
		{
		return new Size (this.width + wincr, this.height + hincr);
		}

	/**
	 * Returns a new size equal to this size decreased by the given amount in
	 * both width and height.
	 *
	 * @param  decr  Decrement.
	 */
	public Size sub
		(double decr)
		{
		return new Size (this.width - decr, this.height - decr);
		}

	/**
	 * Returns a new size equal to this size decreased by the given amounts in
	 * width and height.
	 *
	 * @param  wdecr  Width decrement.
	 * @param  hdecr  Height decrement.
	 */
	public Size sub
		(double wdecr,
		 double hdecr)
		{
		return new Size (this.width - wdecr, this.height - hdecr);
		}

	/**
	 * Returns a new size equal to this size multiplied by the given scale
	 * factor.
	 *
	 * @param  scale  Scale factor.
	 */
	public Size mul
		(double scale)
		{
		return new Size (this.width * scale, this.height * scale);
		}

	/**
	 * Returns a new size equal to this size multiplied by the given scale
	 * factors.
	 *
	 * @param  wscale  Width scale factor.
	 * @param  hscale  Height scale factor.
	 */
	public Size mul
		(double wscale,
		 double hscale)
		{
		return new Size (this.width * wscale, this.height * hscale);
		}

	/**
	 * Returns a new size equal to this size divided by the given scale factor.
	 *
	 * @param  scale  Scale factor.
	 */
	public Size div
		(double scale)
		{
		return new Size (this.width / scale, this.height / scale);
		}

	/**
	 * Returns a new size equal to this size divided by the given scale factors.
	 *
	 * @param  wscale  Width scale factor.
	 * @param  hscale  Height scale factor.
	 */
	public Size div
		(double wscale,
		 double hscale)
		{
		return new Size (this.width / wscale, this.height / hscale);
		}

	/**
	 * Write this size to the given object output stream.
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
		out.writeDouble (this.width);
		out.writeDouble (this.height);
		}

	/**
	 * Read this size from the given object input stream.
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
		this.width = in.readDouble();
		this.height = in.readDouble();
		}

	/**
	 * Determine if this size is equal to the given object.
	 *
	 * @param  obj  Object to test.
	 *
	 * @return  True if this size is equal to <TT>obj</TT>, false otherwise.
	 */
	public boolean equals
		(Object obj)
		{
		if (! (obj instanceof Size)) return false;
		Size that = (Size) obj;
		return this.width == that.width && this.height == that.height;
		}

	/**
	 * Returns a hash code for this size.
	 */
	public int hashCode()
		{
		long widthbits = Double.doubleToLongBits (this.width);
		long heightbits = Double.doubleToLongBits (this.height);
		return
			((int) (widthbits >>> 32)) +
			((int) (widthbits       )) +
			((int) (heightbits >>> 32)) +
			((int) (heightbits       ));
		}

	/**
	 * Returns a string version of this size.
	 */
	public String toString()
		{
		return "(" + this.width + "," + this.height + ")";
		}

	}
