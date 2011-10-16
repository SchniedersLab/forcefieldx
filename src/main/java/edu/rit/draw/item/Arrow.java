//******************************************************************************
//
// File:    Arrow.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Arrow
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
import java.awt.Shape;

import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;

import java.io.Externalizable;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Arrow provides various kinds of arrowheads that can be added to the
 * ends of a {@linkplain Line}.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public class Arrow
	implements Externalizable
	{

// Exported constants.

	/**
	 * A nonexistent arrow.
	 */
	public static final Arrow NONE = new Arrow (0);

	/**
	 * A solid arrow in the shape of a narrow triangle.
	 */
	public static final Arrow SOLID = new Arrow (1);

	/**
	 * An open arrow in the shape of a narrow triangle.
	 */
	public static final Arrow OPEN = new Arrow (2);

// Hidden helper classes.

	/**
	 * Class Arrow.Info is a record containing information needed to draw an
	 * arrow. <TT>thePath</TT> draws the arrow in "standard position" -- with
	 * its tip at (0,0), pointing along the positive X axis, for a stroke width
	 * of 1. If <TT>thePath</TT> is null, the arrow is not drawn.
	 * <TT>isFilled</TT> tells whether the interior of the arrow is filled.
	 *
	 * @author  Alan Kaminsky
	 * @version 10-Jul-2006
	 */
	private static class Info
		{
		public GeneralPath thePath;
		public boolean isFilled;

		public Info
			(GeneralPath thePath,
			 boolean isFilled)
			{
			this.thePath = thePath;
			this.isFilled = isFilled;
			}
		}

// Hidden data members.

	private static final long serialVersionUID = 7782719680651946887L;

	// Table of arrow information records, indexed by arrow kind.
	private static final Info[] theArrowInfo;
	static
		{
		GeneralPath path;
		float goldenRatio = (float) ((1.0 + Math.sqrt (5.0)) / 2.0);

		theArrowInfo = new Info [3];

		theArrowInfo[0] = new Info (null, false);

		path = new GeneralPath();
		//path.moveTo (-4*goldenRatio, -2);
		path.moveTo (-8*goldenRatio, -4);
		path.lineTo (0, 0);
		//path.lineTo (-4*goldenRatio, 2);
		path.lineTo (-8*goldenRatio, 4);
		path.closePath();
		theArrowInfo[1] = new Info (path, true);

		path = new GeneralPath();
		//path.moveTo (-4*goldenRatio, -2);
		path.moveTo (-8*goldenRatio, -4);
		path.lineTo (0, 0);
		//path.lineTo (-4*goldenRatio, 2);
		path.lineTo (-8*goldenRatio, 4);
		theArrowInfo[2] = new Info (path, false);
		}

	// Kind of arrow.
	private int myKind;

// Exported constructors.

	/**
	 * Construct a new nonexistent arrow.
	 */
	public Arrow()
		{
		this (0);
		}

	/**
	 * Construct a new arrow of the given kind.
	 *
	 * @param  theKind  Kind of arrow.
	 */
	private Arrow
		(int theKind)
		{
		myKind = theKind;
		}

// Exported operations.

	/**
	 * Write this arrow to the given object output stream.
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
		out.writeInt (myKind);
		}

	/**
	 * Read this arrow from the given object input stream.
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
		myKind = in.readInt();
		if (0 > myKind || myKind >= theArrowInfo.length)
			{
			throw new InvalidObjectException
				("Invalid arrow (kind = " + myKind + ")");
			}
		}

	/**
	 * Draw this arrow in the given graphics context. It assumes the graphics
	 * context's stroke and paint are already set to the correct values. The
	 * arrow is scaled to match the <TT>width</TT>. The arrow's tip is placed at
	 * the coordinates <TT>(x,y)</TT>. The arrow is rotated so it points in the
	 * direction given by <TT>phi</TT>.
	 *
	 * @param  g2d    2-D graphics context.
	 * @param  width  Stroke width.
	 * @param  x      X coordinate of the arrow's tip.
	 * @param  y      Y coordinate of the arrow's tip.
	 * @param  phi    Angle in which the arrow points (radians).
	 */
	public void draw
		(Graphics2D g2d,
		 float width,
		 double x,
		 double y,
		 double phi)
		{
		Info info =
			0 <= myKind && myKind < theArrowInfo.length ?
				theArrowInfo[myKind] :
				null;
		if (info != null && info.thePath != null)
			{
			AffineTransform pathTransform = new AffineTransform();
			pathTransform.translate (x, y);
			pathTransform.rotate (phi);
			pathTransform.scale (width, width);
			Shape transformedPath =
				info.thePath.createTransformedShape (pathTransform);
			if (info.isFilled)
				{
				g2d.fill (transformedPath);
				}
			g2d.draw (transformedPath);
			}
		}

	}
