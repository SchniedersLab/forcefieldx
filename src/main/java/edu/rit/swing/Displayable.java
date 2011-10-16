//******************************************************************************
//
// File:    Displayable.java
// Package: edu.rit.swing
// Unit:    Interface edu.rit.swing.Displayable
//
// This Java source file is copyright (C) 2002-2004 by Alan Kaminsky. All rights
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

import java.awt.Paint;

import java.awt.geom.Rectangle2D;

/**
 * Interface Displayable specifies the interface for an object that can be
 * displayed in a graphics context. Besides the ability to draw itself
 * (inherited from interface {@link Drawable </CODE>Drawable<CODE>}), a
 * displayable object can tell the range of display coordinates it needs to
 * display itself and the background paint for the display.
 *
 * @author  Alan Kaminsky
 * @version 11-Mar-2003
 */
public interface Displayable
	extends Drawable
	{

// Exported operations.

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 */
	public Rectangle2D getBoundingBox();

	/**
	 * Returns this displayable object's background paint.
	 */
	public Paint getBackgroundPaint();

	}
