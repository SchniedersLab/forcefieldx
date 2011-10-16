//******************************************************************************
//
// File:    Drawable.java
// Package: edu.rit.swing
// Unit:    Interface edu.rit.swing.Drawable
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

import java.awt.Graphics2D;

/**
 * Interface Drawable specifies the interface for an object that can be drawn.
 * The drawable object generally draws itself at the origin of the graphics
 * context, coordinates (0, 0). It is up to the caller to set the graphics
 * context's transformation to place the drawable object at the desired
 * location.
 *
 * @author  Alan Kaminsky
 * @version 11-Mar-2003
 */
public interface Drawable
	{

// Exported operations.

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d);

	}
