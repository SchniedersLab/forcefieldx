//******************************************************************************
//
// File:    Outline.java
// Package: edu.rit.draw.item
// Unit:    Interface edu.rit.draw.item.Outline
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

import java.awt.Graphics2D;

import java.io.Externalizable;

/**
 * Interface Outline specifies the interface for an object that gives the stroke
 * and paint with which to outline an area in a {@linkplain DrawingItem}.
 * <P>
 * All outline objects must be externalizable so they can be stored in a drawing
 * file.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public interface Outline
	extends Externalizable
	{

// Exported constants.

	/**
	 * Specify <TT>Outline.NONE</TT> to omit drawing a drawing item's outline.
	 */
	public static final Outline NONE = null;

// Exported operations.

	/**
	 * Returns the stroke width of this outline.
	 *
	 * @return  Stroke width.
	 */
	public float getStrokeWidth();

	/**
	 * Set the given graphics context's stroke and paint attributes as specified
	 * by this outline object.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void setGraphicsContext
		(Graphics2D g2d);

	}
