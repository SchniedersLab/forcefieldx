//******************************************************************************
//
// File:    Fill.java
// Package: edu.rit.draw.item
// Unit:    Interface edu.rit.draw.item.Fill
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
 * Interface Fill specifies the interface for an object that gives the paint
 * with which to fill an area in a {@linkplain DrawingItem}.
 * <P>
 * All fill objects must be externalizable so they can be stored in a drawing
 * file.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public interface Fill
	extends Externalizable
	{

// Exported constants.

	/**
	 * Specify <TT>Fill.NONE</TT> to omit filling a drawing item's interior.
	 */
	public static final Fill NONE = null;

// Exported operations.

	/**
	 * Set the given graphics context's paint attribute as specified by this
	 * fill object.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void setGraphicsContext
		(Graphics2D g2d);

	}
