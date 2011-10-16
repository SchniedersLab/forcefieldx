//******************************************************************************
//
// File:    Bullet.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Bullet
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

import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Stroke;

import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Bullet provides various kinds of bullets that can be added to a
 * {@linkplain Text} item. The bullet is always added to the left of the text
 * item, with a specified distance from the left side of the bullet to the left
 * side of the text item.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public class Bullet
	implements Externalizable
	{

// Hidden helper classes.

	/**
	 * Class Bullet.Info is a record containing information needed to draw a
	 * bullet. <TT>theArea</TT> draws the bullet in "standard position" -- with
	 * its lower left corner at (0,0), for a text font size of 1. If
	 * <TT>theArea</TT> is null, the bullet is not drawn. <TT>isFilled</TT>
	 * tells whether the interior of the bullet is filled.
	 * <TT>theOutlineWidth</TT> gives the width of the outline stroke for the
	 * bullet; if 0, the outline is not drawn.
	 *
	 * @author  Alan Kaminsky
	 * @version 10-Jul-2006
	 */
	private static class Info
		{
		public Area theArea;
		public boolean isFilled;
		public float theOutlineWidth;

		public Info
			(Area theArea,
			 boolean isFilled,
			 float theOutlineWidth)
			{
			this.theArea = theArea;
			this.isFilled = isFilled;
			this.theOutlineWidth = theOutlineWidth;
			}
		}

// Hidden constants.

	/**
	 * A kind of bullet with no shape.
	 */
	private static final int BULLET_KIND_NONE = 0;

	/**
	 * A kind of bullet in the shape of a medium sized dot.
	 */
	private static final int BULLET_KIND_DOT = 1;

	/**
	 * A kind of bullet in the shape of a medium sized circle.
	 */
	private static final int BULLET_KIND_CIRCLE = 2;

// Exported constants.

	/**
	 * The normal bullet offset, 36.0 (1/2 inch).
	 */
	public static final double NORMAL_OFFSET = 36.0;

	/**
	 * A nonexistent bullet.
	 */
	public static final Bullet NONE = null;

	/**
	 * A bullet in the shape of a medium sized dot, offset the normal offset
	 * (36.0, or 1/2 inch) from the text item.
	 */
	public static final Bullet DOT =
		new Bullet (BULLET_KIND_DOT, NORMAL_OFFSET);

	/**
	 * A bullet in the shape of a medium sized circle, offset the normal offset
	 * (36.0, or 1/2 inch) from the text item.
	 */
	public static final Bullet CIRCLE =
		new Bullet (BULLET_KIND_CIRCLE, NORMAL_OFFSET);

// Hidden data members.

	private static final long serialVersionUID = 1799559688305790072L;

	// Table of bullet information records, indexed by bullet kind.
	private static final Info[] theBulletInfo;
	static
		{
		theBulletInfo = new Info [3];

		theBulletInfo[0] = new Info (null, false, 0.0f);

		theBulletInfo[1] =
			new Info
				(new Area (new Ellipse2D.Double (0.0, -0.625, 0.5, 0.5)),
				 true,
				 0.0f);

		theBulletInfo[2] =
			new Info
				(new Area (new Ellipse2D.Double (0.0, -0.625, 0.5, 0.5)),
				 false,
				 1.0f);
		}

	// Kind of bullet.
	private int myKind;

	// Offset from left end of bullet to left end of text box.
	private double myOffset = theDefaultOffset;

	// The default offset.
	private static double theDefaultOffset = NORMAL_OFFSET;

// Exported constructors.

	/**
	 * Construct a new nonexistent bullet.
	 */
	public Bullet()
		{
		}

	/**
	 * Construct a new bullet of the same kind as the given bullet, but with the
	 * default offset.
	 *
	 * @param  theBullet  Bullet.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theBullet</TT> is null.
	 */
	public Bullet
		(Bullet theBullet)
		{
		myKind = theBullet.myKind;
		myOffset = theDefaultOffset;
		}

// Hidden constructors.

	/**
	 * Construct a new bullet of the given kind with the given offset.
	 *
	 * @param  theKind    Kind of bullet (one of the <TT>BULLET_KIND_XXX</TT>
	 *                    values).
	 * @param  theOffset  Offset from left side of bullet to left side of text.
	 */
	private Bullet
		(int theKind,
		 double theOffset)
		{
		myKind = theKind;
		myOffset = theOffset;
		}

// Exported operations.

	/**
	 * Returns the default offset for bullets. The offset is the distance
	 * between the left side of the bullet and the left side of the text.
	 *
	 * @return  Default offset.
	 */
	public static double defaultOffset()
		{
		return theDefaultOffset;
		}

	/**
	 * Set the default offset for bullets. The offset is the distance between
	 * the left side of the bullet and the left side of the text.
	 *
	 * @param  theOffset  Default offset.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theOffset</TT> is less than 0.
	 */
	public static void defaultOffset
		(double theOffset)
		{
		if (theOffset < 0.0) throw new IllegalArgumentException();
		theDefaultOffset = theOffset;
		DOT.myOffset = theOffset;
		CIRCLE.myOffset = theOffset;
		}

	/**
	 * Returns this bullet's offset. The offset is the distance between the left
	 * side of the bullet and the left side of the text.
	 *
	 * @return  Offset.
	 */
	public double offset()
		{
		return myOffset;
		}

	/**
	 * Set this bullet's offset. The offset is the distance between the left
	 * side of the bullet and the left side of the text.
	 *
	 * @param  theOffset  Offset.
	 *
	 * @return  This bullet.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theOffset</TT> is less than 0.
	 */
	public Bullet offset
		(double theOffset)
		{
		if (theOffset < 0.0) throw new IllegalArgumentException();
		myOffset = theOffset;
		return this;
		}

	/**
	 * Write this bullet to the given object output stream.
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
		out.writeDouble (myOffset);
		}

	/**
	 * Read this bullet from the given object input stream.
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
		myOffset = in.readDouble();
		if (0 > myKind || myKind >= theBulletInfo.length)
			{
			throw new InvalidObjectException
				("Invalid bullet (kind = " + myKind + ")");
			}
		}

	/**
	 * Draw this bullet in the given graphics context. It assumes the graphics
	 * context's paint is already set to the correct value. The bullet is scaled
	 * to match the <TT>ascent</TT>. The bullet's left side is placed relative
	 * to the coordinates of the left end of the text baseline <TT>(x,y)</TT>.
	 *
	 * @param  g2d     2-D graphics context.
	 * @param  ascent  Text font ascent.
	 * @param  x       X coordinate of the left end of the text baseline.
	 * @param  y       Y coordinate of the left end of the text baseline.
	 */
	public void draw
		(Graphics2D g2d,
		 double ascent,
		 double x,
		 double y)
		{
		Info info =
			0 <= myKind && myKind < theBulletInfo.length ?
				theBulletInfo[myKind] :
				null;
		if (info != null && info.theArea != null)
			{
			AffineTransform areaTransform = new AffineTransform();
			areaTransform.translate (x-myOffset, y);
			areaTransform.scale (ascent, ascent);
			Area transformedArea =
				info.theArea.createTransformedArea (areaTransform);
			if (info.isFilled)
				{
				g2d.fill (transformedArea);
				}
			if (info.theOutlineWidth > 0.0f)
				{
				Stroke oldStroke = g2d.getStroke();
				g2d.setStroke (new BasicStroke (info.theOutlineWidth));
				g2d.draw (transformedArea);
				g2d.setStroke (oldStroke);
				}
			}
		}

	}
