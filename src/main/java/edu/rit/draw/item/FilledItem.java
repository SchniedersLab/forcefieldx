//******************************************************************************
//
// File:    FilledItem.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.FilledItem
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

package edu.rit.draw.item;

import edu.rit.draw.Drawing;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class FilledItem is the abstract base class for a {@linkplain DrawingItem}
 * that has an outline and is filled with a paint. Specify <TT>Fill.NONE</TT>
 * for the fill paint to omit filling a drawing item's interior.
 * <P>
 * The static <TT>defaultFill()</TT> method is provided to set the default fill
 * paint. If the fill paint is not specified, the current default fill paint is
 * used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public abstract class FilledItem
	extends OutlinedItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal fill paint: White.
	 */
	public static final Fill NORMAL_FILL = ColorFill.WHITE;

// Hidden data members.

	private static final long serialVersionUID = -5025247969110396992L;

	// Fill paint, or null if the item's interior is not filled.
	Fill myFill = theDefaultFill;

	// The default fill paint.
	static Fill theDefaultFill = NORMAL_FILL;

// Exported constructors.

	/**
	 * Construct a new filled item with the default fill paint.
	 */
	public FilledItem()
		{
		super();
		}

	/**
	 * Construct a new filled item with the same outline and fill paint as the
	 * given filled item.
	 *
	 * @param  theItem  Filled item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public FilledItem
		(FilledItem theItem)
		{
		super (theItem);
		myFill = theItem.myFill;
		}

// Exported operations.

	/**
	 * Returns the default fill paint.
	 *
	 * @return  Default fill paint, or <TT>Fill.NONE</TT>.
	 */
	public static Fill defaultFill()
		{
		return theDefaultFill;
		}

	/**
	 * Set the default fill paint. Before calling this method the first time,
	 * the default fill paint is white.
	 *
	 * @param  theFill  Default fill paint, or <TT>Fill.NONE</TT>.
	 */
	public static void defaultFill
		(Fill theFill)
		{
		theDefaultFill = theFill;
		}

	/**
	 * Set this filled item's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This filled item.
	 */
	public FilledItem outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Returns this filled item's fill paint.
	 *
	 * @return  Fill paint, or <TT>Fill.NONE</TT>.
	 */
	public Fill fill()
		{
		return myFill;
		}

	/**
	 * Set this filled item's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This filled item.
	 */
	public FilledItem fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Add this filled item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This filled item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public FilledItem add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this filled item to the end of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This filled item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public FilledItem add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this filled item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This filled item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public FilledItem addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this filled item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This filled item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public FilledItem addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this filled item to the given object output stream.
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
		super.writeExternal (out);
		out.writeObject (myFill);
		}

	/**
	 * Read this filled item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this filled item cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myFill = (Fill) in.readObject();
		}

// Hidden operations.

	/**
	 * Set this filled item's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 */
	void doFill
		(Fill theFill)
		{
		myFill = theFill;
		}

	}
