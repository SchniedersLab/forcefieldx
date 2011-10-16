//******************************************************************************
//
// File:    DisplayableList.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DisplayableList
//
// This Java source file is copyright (C) 2002-2005 by Alan Kaminsky. All rights
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

package edu.rit.swing;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;

import java.awt.geom.Rectangle2D;

import java.util.Vector;

/**
 * Class DisplayableList provides an object that combines a list of separate
 * {@link Displayable </CODE>Displayable<CODE>} objects into a single
 * displayable object.
 * <P>
 * <I>Note:</I> Class DisplayableList is multiple thread safe.
 *
 * @author  Alan Kaminsky
 * @version 20-Jul-2006
 */
public class DisplayableList
	implements Displayable
	{

// Hidden data members.

	private Vector<Displayable> myDisplayList = new Vector<Displayable>();

	private static final Rectangle2D EMPTY = new Rectangle2D.Double();
	private Rectangle2D myBoundingBox = EMPTY;

	private Paint myBackgroundPaint = Color.white;

// Exported constructors.

	/**
	 * Construct a new, empty displayable list.
	 */
	public DisplayableList()
		{
		}

// Exported operations.

	/**
	 * Clear this displayable list.
	 */
	public synchronized void clear()
		{
		myDisplayList.clear();
		myBoundingBox = EMPTY;
		myBackgroundPaint = Color.white;
		}

	/**
	 * Add the given displayable object to the end of this displayable list.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public synchronized void add
		(Displayable theDisplayable)
		{
		if (theDisplayable == null)
			{
			throw new NullPointerException();
			}
		if (myDisplayList.isEmpty())
			{
			myBackgroundPaint = theDisplayable.getBackgroundPaint();
			}
		myBoundingBox =
			myBoundingBox.createUnion (theDisplayable.getBoundingBox());
		myDisplayList.add (theDisplayable);
		}

// Exported operations implemented from interface Drawable.

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 * <P>
	 * The displayable objects in this displayable list are drawn in the order
	 * they were added.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public synchronized void draw
		(Graphics2D g2d)
		{
		for (Displayable item : myDisplayList)
			{
			item.draw (g2d);
			}
		}

// Exported operations implemented from interface Displayable.

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 * <P>
	 * This method returns the union of the bounding boxes of the displayable
	 * objects in this displayable list. If this displayable list is empty, an
	 * empty rectangle is returned.
	 */
	public synchronized Rectangle2D getBoundingBox()
		{
		return myBoundingBox;
		}

	/**
	 * Returns this displayable object's background paint.
	 * <P>
	 * This method returns the background paint of the first displayable object
	 * in this displayable list. If this displayable list is empty, the color
	 * white is returned.
	 */
	public synchronized Paint getBackgroundPaint()
		{
		return myBackgroundPaint;
		}

	}
