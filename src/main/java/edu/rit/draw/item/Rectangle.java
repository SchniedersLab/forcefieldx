//******************************************************************************
//
// File:    Rectangle.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Rectangle
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

import java.awt.Shape;

import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Rectangle provides a rectangle {@linkplain DrawingItem}. The rectangle
 * may have sharp corners or round corners.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0 WIDTH=100%>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * If the rectangle has round corners, the round corner distance gives the
 * distance from where the (sharp) corner would normally be to the starting and
 * ending points of the round corner arc, as shown in the adjoining diagram.
 * (The round corner distance is <I>d</I>.) A round corner distance of 0 results
 * in a sharp corner. In a particular rectangle, the same round corner distance
 * is used for every corner; that is, either all the corners are sharp or all
 * the corners are rounded the same amount. Different rectangles can have
 * different round corner distances.
 * </TD>
 * <TD WIDTH=10> </TD>
 * <TD ALIGN="left" VALIGN="top" WIDTH=135><IMG SRC="doc-files/RectangleFig01.png"></TD>
 * </TR>
 * </TABLE>
 * <P>
 * The static <TT>defaultSize()</TT>, <TT>defaultWidth()</TT>,
 * <TT>defaultHeight()</TT>, and <TT>defaultRound()</TT> methods are provided to
 * set the default size, width, height, and round corner distance. If a
 * rectangle's size, width, height, or round corner distance is not specified,
 * the current default size, width, height, or round corner distance is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Rectangle
	extends RectangularItem
	{

// Exported constants.

	/**
	 * The normal width for rectangles (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_WIDTH = 72.0;

	/**
	 * The normal height for rectangles (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_HEIGHT = 72.0;

	/**
	 * The normal round corner distance (0), signifying sharp corners.
	 */
	public static final double NORMAL_ROUND = 0.0;

// Hidden data members.

	private static final long serialVersionUID = -1131516846619887544L;

	// Attributes.
	private double myRound = theDefaultRound;

	// The default attributes.
	private static double theDefaultWidth = NORMAL_WIDTH;
	private static double theDefaultHeight = NORMAL_HEIGHT;
	private static double theDefaultRound = NORMAL_ROUND;

// Exported constructors.

	/**
	 * Construct a new rectangle. The rectangle's northwest corner is located at
	 * (0,0). The rectangle's size is the default size. The rectangle has sharp
	 * corners.
	 */
	public Rectangle()
		{
		super();
		this.width = theDefaultWidth;
		this.height = theDefaultHeight;
		this.myRound = theDefaultRound;
		}

	/**
	 * Construct a new rectangle with the same outline, fill paint, location,
	 * and size as the given rectangle.
	 *
	 * @param  theItem  Rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Rectangle
		(Rectangle theItem)
		{
		super (theItem);
		this.myRound = theItem.myRound;
		}

// Exported operations.

	/**
	 * Returns the default size for rectangles.
	 *
	 * @return  Default size.
	 */
	public static Size defaultSize()
		{
		return new Size (theDefaultWidth, theDefaultHeight);
		}

	/**
	 * Set the default size for rectangles.
	 *
	 * @param  theSize  Default size.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public static void defaultSize
		(Size theSize)
		{
		if (theSize.width < 0.0 || theSize.height < 0.0)
			{
			throw new IllegalArgumentException();
			}
		theDefaultWidth = theSize.width;
		theDefaultHeight = theSize.height;
		}

	/**
	 * Returns the default width for rectangles.
	 *
	 * @return  Default width.
	 */
	public static double defaultWidth()
		{
		return theDefaultWidth;
		}

	/**
	 * Set the default width for rectangles.
	 *
	 * @param  theWidth  Default width.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public static void defaultWidth
		(double theWidth)
		{
		if (theWidth < 0.0) throw new IllegalArgumentException();
		theDefaultWidth = theWidth;
		}

	/**
	 * Returns the default height for rectangles.
	 *
	 * @return  Default height.
	 */
	public static double defaultHeight()
		{
		return theDefaultHeight;
		}

	/**
	 * Set the default height for rectangles.
	 *
	 * @param  theHeight  Default height.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public static void defaultHeight
		(double theHeight)
		{
		if (theHeight < 0.0) throw new IllegalArgumentException();
		theDefaultHeight = theHeight;
		}

	/**
	 * Returns the default round corner distance. A value of 0 signifies sharp
	 * corners.
	 *
	 * @return  Default round corner distance.
	 */
	public static double defaultRound()
		{
		return theDefaultRound;
		}

	/**
	 * Set the default round corner distance. A value of 0 signifies sharp
	 * corners. Before calling this method the first time, the default round
	 * corner distance is 0.
	 *
	 * @param  theRound  Default round corner distance.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRound</TT> is less than 0.
	 */
	public static void defaultRound
		(double theRound)
		{
		if (theRound < 0.0) throw new IllegalArgumentException();
		theDefaultRound = theRound;
		}

	/**
	 * Set this rectangle's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this rectangle's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Set the size of this rectangle's bounding box.
	 *
	 * @param  theSize  Size.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public Rectangle size
		(Size theSize)
		{
		doSize (theSize);
		return this;
		}

	/**
	 * Set the width of this rectangle's bounding box.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public Rectangle width
		(double theWidth)
		{
		doWidth (theWidth);
		return this;
		}

	/**
	 * Set the height of this rectangle's bounding box.
	 *
	 * @param  theHeight  Height.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public Rectangle height
		(double theHeight)
		{
		doHeight (theHeight);
		return this;
		}

	/**
	 * Set the northwest corner point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this rectangle's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this rectangle's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this rectangle's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this rectangle's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this rectangle's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this rectangle's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this rectangle's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this rectangle's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this rectangle's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This rectangle.
	 */
	public Rectangle se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this rectangle's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Rectangle se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Returns this rectangle's round corner distance. A value of 0 signifies
	 * sharp corners.
	 *
	 * @return  Round corner distance.
	 */
	public double round()
		{
		return myRound;
		}

	/**
	 * Set this rectangle's round corner distance. A value of 0 signifies sharp
	 * corners.
	 *
	 * @param  theRound  Round corner distance.
	 *
	 * @return  This line.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRound</TT> is less than 0.
	 */
	public Rectangle round
		(double theRound)
		{
		if (theRound < 0.0) throw new IllegalArgumentException();
		myRound = theRound;
		return this;
		}

	/**
	 * Add this rectangle to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Rectangle add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this rectangle to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Rectangle add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this rectangle to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Rectangle addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this rectangle to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This rectangle.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Rectangle addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this rectangle to the given object output stream.
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
		out.writeDouble (myRound);
		}

	/**
	 * Read this rectangle from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this rectangle cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myRound = in.readDouble();
		}

// Hidden operations.

	/**
	 * Determine the 2-D graphics shape that this shape object will draw.
	 *
	 * @return  Shape.
	 */
	Shape getShape()
		{
		Point nw = nw();
		return
			myRound == 0.0 ?
				new Rectangle2D.Double
					(nw.x, nw.y, width, height) :
				new RoundRectangle2D.Double
					(nw.x, nw.y, width, height, 2*myRound, 2*myRound);
		}

	}
