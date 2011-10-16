//******************************************************************************
//
// File:    OutlinedItem.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.OutlinedItem
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
 * Class OutlinedItem is the abstract base class for a {@linkplain DrawingItem}
 * that has an outline. Specify <TT>Outline.NONE</TT> for the outline to omit
 * drawing a drawing item's outline.
 * <P>
 * The static <TT>defaultOutline()</TT> method is provided to set the default
 * outline. If the outline is not specified, the current default outline is
 * used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public abstract class OutlinedItem
	extends DrawingItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal outline: Solid, square corners, 1 point wide, black.
	 */
	public static final Outline NORMAL_OUTLINE = SolidOutline.NORMAL_OUTLINE;

// Hidden data members.

	private static final long serialVersionUID = -2644075901753850353L;

	// Outline, or null if the item's outline is not drawn.
	Outline myOutline = theDefaultOutline;

	// The default outline.
	static Outline theDefaultOutline = NORMAL_OUTLINE;

// Exported constructors.

	/**
	 * Construct a new outlined item with the default outline.
	 */
	public OutlinedItem()
		{
		super();
		}

	/**
	 * Construct a new outlined item with the same outline as the given outlined
	 * item.
	 *
	 * @param  theItem  Outlined item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public OutlinedItem
		(OutlinedItem theItem)
		{
		super (theItem);
		myOutline = theItem.myOutline;
		}

// Exported operations.

	/**
	 * Returns the default outline.
	 *
	 * @return  Default outline, or <TT>Outline.NONE</TT>.
	 */
	public static Outline defaultOutline()
		{
		return theDefaultOutline;
		}

	/**
	 * Set the default outline. Before calling this method the first time,
	 * the default outline is solid, square corners, 1 point wide, black.
	 *
	 * @param  theOutline  Default outline, or <TT>Outline.NONE</TT>.
	 */
	public static void defaultOutline
		(Outline theOutline)
		{
		theDefaultOutline = theOutline;
		}

	/**
	 * Returns this outlined item's outline.
	 *
	 * @return  Outline, or <TT>Outline.NONE</TT>.
	 */
	public Outline outline()
		{
		return myOutline;
		}

	/**
	 * Set this outlined item's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This outlined item.
	 */
	public OutlinedItem outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Add this outlined item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This outlined item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public OutlinedItem add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this outlined item to the end of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This outlined item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public OutlinedItem add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this outlined item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This outlined item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public OutlinedItem addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this outlined item to the beginning of the given drawing's sequence
	 * of drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This outlined item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public OutlinedItem addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this outlined item to the given object output stream.
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
		out.writeObject (myOutline);
		}

	/**
	 * Read this outlined item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this outlined item cannot
	 *     be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myOutline = (Outline) in.readObject();
		}

// Hidden operations.

	/**
	 * Set this outlined item's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 */
	void doOutline
		(Outline theOutline)
		{
		myOutline = theOutline;
		}

	}
