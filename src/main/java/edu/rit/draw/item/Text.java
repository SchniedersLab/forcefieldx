//******************************************************************************
//
// File:    Text.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Text
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

import java.awt.Font;
import java.awt.Graphics2D;

import java.awt.font.FontRenderContext;
import java.awt.font.LineMetrics;
import java.awt.font.TextAttribute;
import java.awt.font.TextLayout;

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.text.CharacterIterator;
import java.text.AttributedCharacterIterator;
import java.text.AttributedString;

import java.util.HashMap;
import java.util.Map;

/**
 * Class Text provides a textual {@linkplain DrawingItem}. A text item has these
 * attributes:
 * <UL>
 * <LI>
 * <B>Text</B> -- The text to be displayed. The text consists of one or more
 * lines. A newline character (<TT>'\n'</TT>) separates each line from the next.
 * <LI>
 * <B>Font</B> -- The font in which to display the text. All the text is
 * displayed in the same font.
 * <LI>
 * <B>Fill</B> -- The fill paint with which to display the text.
 * <LI>
 * <B>Spacing</B> -- The vertical distance between consecutive lines. This is
 * specified as a factor which is multiplied by the font size to get the actual
 * inter-line distance.
 * <LI>
 * <B>Alignment</B> -- Whether multiple lines are aligned left, center, or
 * right.
 * <LI>
 * <B>Bullet</B> -- The bullet, if any, attached to the text. The bullet is
 * always displayed to the left of the first line of the text.
 * <LI>
 * <B>Position</B> -- The text occupies a rectangular area. The text's position
 * is specified by giving the (<I>x,y</I>) coordinates of the text area's
 * northwest, north, northeast, west, center, east, southwest, south, or
 * southeast point.
 * </UL>
 * <P>
 * Class Text supports a subset of HTML for specifying attributes of characters.
 * The following HTML tags may be embedded in the text:
 * <UL>
 * <LI><TT>&lt;B&gt;Bold&lt;/B&gt;</TT>
 * <LI><TT>&lt;I&gt;Italic&lt;/I&gt;</TT>
 * <LI><TT>&lt;SUP&gt;Superscript&lt;/SUP&gt;</TT>
 * <LI><TT>&lt;SUB&gt;Subscript&lt;/SUB&gt;</TT>
 * </UL>
 * <P>
 * To get literal less-than, greater-than, and ampersand characters in the text,
 * use the following HTML character entities:
 * <UL>
 * <LI><TT>&amp;lt;</TT> -- Less-than (&lt;)
 * <LI><TT>&amp;gt;</TT> -- Greater-than (&gt;)
 * <LI><TT>&amp;amp;</TT> -- Ampersand (&amp;)
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Text
	extends DrawingItem
	implements Externalizable
	{

// Exported enumerations.

	/**
	 * Enumeration Text.Alignment specifies the alignment of multiple lines in a
	 * {@linkplain Text} item.
	 *
	 * @author  Alan Kaminsky
	 * @version 11-Jul-2006
	 */
	public static enum Alignment
		{
		/**
		 * Align the left sides of multiple text lines.
		 */
		LEFT,
		/**
		 * Align the center points of multiple text lines.
		 */
		CENTER,
		/**
		 * Align the right sides of multiple text lines.
		 */
		RIGHT,
		}

// Exported constants.

	/**
	 * The normal font (sanserif, plain, 12 point).
	 */
	public static final Font NORMAL_FONT =
		new Font ("sanserif", Font.PLAIN, 12);

	/**
	 * The normal fill paint (black).
	 */
	public static final Fill NORMAL_FILL = ColorFill.BLACK;

	/**
	 * The normal line spacing factor (7/6).
	 */
	public static final double NORMAL_SPACING = 7.0 / 6.0;

	/**
	 * The normal alignment (left).
	 */
	public static final Alignment NORMAL_ALIGNMENT = Alignment.LEFT;

	/**
	 * The normal bullet (none).
	 */
	public static final Bullet NORMAL_BULLET = Bullet.NONE;

// Hidden data members.

	private static final long serialVersionUID = -3141850073230054560L;

	// Default attributes.
	static Font theDefaultFont = NORMAL_FONT;
	static Fill theDefaultFill = NORMAL_FILL;
	static double theDefaultSpacing = NORMAL_SPACING;
	static Alignment theDefaultAlignment = NORMAL_ALIGNMENT;
	static Bullet theDefaultBullet = NORMAL_BULLET;

	// Text. If null, the text is empty.
	String myText;

	// Attributes.
	Font myFont = theDefaultFont;
	Fill myFill = theDefaultFill;
	double mySpacing = theDefaultSpacing;
	Alignment myAlignment = theDefaultAlignment;
	Bullet myBullet = theDefaultBullet;

	// Coordinates of most recently specified corner.
	double x;
	double y;

	// Factors for going from specified corner to northwest corner.
	double xFactor;
	double yFactor;

	// Size. If null, the size must be recomputed.
	Size mySize;

	// Array of attributed strings, one per line, resulting from rendering the
	// HTML text. If null, the attributed strings must be recomputed.
	AttributedString[] myLines;

// Exported constructors.

	/**
	 * Construct a new text item. The text item's text is empty. The text item's
	 * northwest corner is located at (0,0). The text item's attributes have the
	 * default values.
	 */
	public Text()
		{
		super();
		}

	/**
	 * Construct a new text item with the same text, location, and attributes
	 * as the given text item.
	 *
	 * @param  theItem  Text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Text
		(Text theItem)
		{
		super (theItem);
		this.myText = theItem.myText;
		this.myFont = theItem.myFont;
		this.myFill = theItem.myFill;
		this.mySpacing = theItem.mySpacing;
		this.myAlignment = theItem.myAlignment;
		this.myBullet = theItem.myBullet;
		this.x = theItem.x;
		this.y = theItem.y;
		this.xFactor = theItem.xFactor;
		this.yFactor = theItem.yFactor;
		this.mySize = theItem.mySize;
		if (theItem.myLines == null)
			{
			this.myLines = null;
			}
		else
			{
			int n = theItem.myLines.length;
			this.myLines = new AttributedString [n];
			for (int i = 0; i < n; ++ i)
				{
				this.myLines[n] =
					new AttributedString (theItem.myLines[n].getIterator());
				}
			}
		}

// Exported operations.

	/**
	 * Returns the default font.
	 *
	 * @return  Default font.
	 */
	public static Font defaultFont()
		{
		return theDefaultFont;
		}

	/**
	 * Set the default font. Before calling this method the first time, the
	 * default font is sanserif, plain, 12 point.
	 *
	 * @param  theFont  Default font.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 */
	public static void defaultFont
		(Font theFont)
		{
		if (theFont == null) throw new NullPointerException();
		theDefaultFont = theFont;
		}

	/**
	 * Returns the default fill paint.
	 *
	 * @return  Default fill paint.
	 */
	public static Fill defaultFill()
		{
		return theDefaultFill;
		}

	/**
	 * Set the default fill paint. Before calling this method the first time,
	 * the default fill paint is black.
	 *
	 * @param  theFill  Default fill paint.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFill</TT> is null.
	 */
	public static void defaultFill
		(Fill theFill)
		{
		if (theFill == null) throw new NullPointerException();
		theDefaultFill = theFill;
		}

	/**
	 * Returns the default line spacing factor. The line spacing factor is
	 * multiplied by the font size to get the actual inter-line distance.
	 *
	 * @return  Default line spacing factor.
	 */
	public static double defaultSpacing()
		{
		return theDefaultSpacing;
		}

	/**
	 * Set the default line spacing factor. The line spacing factor is
	 * multiplied by the font size to get the actual inter-line distance. Before
	 * calling this method the first time, the default line spacing factor is
	 * 7/6.
	 *
	 * @param  theSpacing  Default line spacing factor.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theSpacing</TT> is less than 0.
	 */
	public static void defaultSpacing
		(double theSpacing)
		{
		if (theSpacing < 0.0) throw new IllegalArgumentException();
		theDefaultSpacing = theSpacing;
		}

	/**
	 * Returns the default text alignment.
	 *
	 * @return  Default text alignment.
	 */
	public static Alignment defaultAlign()
		{
		return theDefaultAlignment;
		}

	/**
	 * Set the default text alignment. Before calling this method the first
	 * time, the default text alignment is left alignment.
	 *
	 * @param  theAlignment  Default text alignment.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theAlignment</TT> is null.
	 */
	public static void defaultAlign
		(Alignment theAlignment)
		{
		if (theAlignment == null) throw new NullPointerException();
		theDefaultAlignment = theAlignment;
		}

	/**
	 * Returns the default bullet.
	 *
	 * @return  Default bullet.
	 */
	public static Bullet defaultBullet()
		{
		return theDefaultBullet;
		}

	/**
	 * Set the default bullet. Before calling this method the first time, the
	 * default bullet is none.
	 *
	 * @param  theBullet  Default bullet.
	 */
	public static void defaultBullet
		(Bullet theBullet)
		{
		theDefaultBullet = theBullet;
		}

	/**
	 * Returns this text item's text. Newline characters (<TT>'\n'</TT>)
	 * separate multiple lines within the text.
	 *
	 * @return  Text (may be null).
	 */
	public String text()
		{
		return myText;
		}

	/**
	 * Set this text item's text. Newline characters (<TT>'\n'</TT>) separate
	 * multiple lines within the text.
	 *
	 * @param  theText  Text (may be null).
	 *
	 * @return  This text item.
	 */
	public Text text
		(String theText)
		{
		myText = theText;
		mySize = null;
		myLines = null;
		return this;
		}

	/**
	 * Returns this text item's font.
	 *
	 * @return  Font.
	 */
	public Font font()
		{
		return myFont;
		}

	/**
	 * Set this text item's font.
	 *
	 * @param  theFont  Font.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 */
	public Text font
		(Font theFont)
		{
		if (theFont == null) throw new NullPointerException();
		myFont = theFont;
		mySize = null;
		myLines = null;
		return this;
		}

	/**
	 * Returns this text item's fill paint.
	 *
	 * @return  Fill paint.
	 */
	public Fill fill()
		{
		return myFill;
		}

	/**
	 * Set this text item's fill paint.
	 *
	 * @param  theFill  Fill paint.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFill</TT> is null.
	 */
	public Text fill
		(Fill theFill)
		{
		if (theFill == null) throw new NullPointerException();
		myFill = theFill;
		return this;
		}

	/**
	 * Returns this text item's line spacing factor. The line spacing factor is
	 * multiplied by the font size to get the actual inter-line distance.
	 *
	 * @return  Line spacing factor.
	 */
	public double spacing()
		{
		return mySpacing;
		}

	/**
	 * Set this text item's line spacing factor. The line spacing factor is
	 * multiplied by the font size to get the actual inter-line distance.
	 *
	 * @param  theSpacing  Line spacing factor.
	 *
	 * @return  This text item.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theSpacing</TT> is less than 0.
	 */
	public Text spacing
		(double theSpacing)
		{
		if (theSpacing < 0.0) throw new IllegalArgumentException();
		mySpacing = theSpacing;
		mySize = null;
		return this;
		}

	/**
	 * Returns this text item's text alignment.
	 *
	 * @return  Text alignment.
	 */
	public Alignment align()
		{
		return myAlignment;
		}

	/**
	 * Set this text item's text alignment.
	 *
	 * @param  theAlignment  Text alignment.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theAlignment</TT> is null.
	 */
	public Text align
		(Alignment theAlignment)
		{
		if (theAlignment == null) throw new NullPointerException();
		myAlignment = theAlignment;
		return this;
		}

	/**
	 * Returns this text item's bullet.
	 *
	 * @return  Bullet.
	 */
	public Bullet bullet()
		{
		return myBullet;
		}

	/**
	 * Set this text item's bullet.
	 *
	 * @param  theBullet  Default bullet.
	 *
	 * @return  This text item.
	 */
	public Text bullet
		(Bullet theBullet)
		{
		myBullet = theBullet;
		return this;
		}

	/**
	 * Returns the width of this text item's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		computeSize();
		return mySize.width;
		}

	/**
	 * Returns the height of this text item's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		computeSize();
		return mySize.height;
		}

	/**
	 * Returns the northwest corner point of this text item's bounding box.
	 *
	 * @return  Northwest corner point.
	 */
	public Point nw()
		{
		computeSize();
		return new Point
			(x + xFactor * mySize.width, y + yFactor * mySize.height);
		}

	/**
	 * Set the northwest corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This text item.
	 */
	public Text nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this text item's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This text item.
	 */
	public Text n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this text item's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This text item.
	 */
	public Text ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this text item's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This text item.
	 */
	public Text w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this text item's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This text item.
	 */
	public Text c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this text item's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This text item.
	 */
	public Text e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this text item's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This text item.
	 */
	public Text sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this text item's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This text item.
	 */
	public Text s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this text item's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This text item.
	 */
	public Text se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this text item's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Text se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Add this text item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Text add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this text item to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Text add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this text item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Text addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this text item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This text item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Text addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this text item to the given object output stream.
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
		out.writeObject (myText);
		out.writeObject (myFont);
		out.writeObject (myFill);
		out.writeDouble (mySpacing);
		out.writeObject (myAlignment);
		out.writeObject (myBullet);
		out.writeDouble (x);
		out.writeDouble (y);
		out.writeDouble (xFactor);
		out.writeDouble (yFactor);
		out.writeObject (mySize);
		}

	/**
	 * Read this text item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this text item cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myText = (String) in.readObject();
		myFont = (Font) in.readObject();
		myFill = (Fill) in.readObject();
		mySpacing = in.readDouble();
		myAlignment = (Alignment) in.readObject();
		myBullet = (Bullet) in.readObject();
		x = in.readDouble();
		y = in.readDouble();
		xFactor = in.readDouble();
		yFactor = in.readDouble();
		mySize = (Size) in.readObject();
		myLines = null;
		}

	/**
	 * Draw this text item in the given graphics context. This method is
	 * allowed to change the graphics context's paint, stroke, and transform,
	 * and it doesn't have to change them back.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		super.draw (g2d);

		// Set up font and fill.
		Font oldFont = g2d.getFont();
		g2d.setFont (myFont);
		myFill.setGraphicsContext (g2d);

		// Render HTML text.
		renderLines();
		int nlines = myLines.length;

		// Convert each line to a TextLayout object.
		FontRenderContext frc = g2d.getFontRenderContext();
		TextLayout[] layoutLines = new TextLayout [nlines];
		for (int i = 0; i < nlines; ++ i)
			{
			layoutLines[i] = new TextLayout (myLines[i].getIterator(), frc);
			}

		// Determine maximum width of any line.
		double w = 0.0;
		for (TextLayout line : layoutLines)
			{
			w = Math.max (w, line.getBounds().getWidth());
			}

		// Determine height of any line.
		LineMetrics metrics = myFont.getLineMetrics ("Eg", frc);
		double a = metrics.getAscent();
		double d = metrics.getDescent();

		// Determine line spacing.
		double s = mySpacing * myFont.getSize2D();

		// Determine height of all lines.
		double h = nlines == 0 ? 0.0 : (nlines-1)*s + a + d;

		// Determine northwest corner point.
		double nw_x = x + xFactor * w;
		double nw_y = y + yFactor * h;

		// Determine X alignment factor.
		double xAlign = 0.0;
		switch (myAlignment)
			{
			case LEFT:
				xAlign = 0.0;
				break;
			case CENTER:
				xAlign = 0.5;
				break;
			case RIGHT:
				xAlign = 1.0;
				break;
			}

		// Draw all lines.
		for (int i = 0; i < nlines; ++ i)
			{
			// Draw line.
			TextLayout line = layoutLines[i];
			Rectangle2D rect = line.getBounds();
			double line_w = rect.getWidth();
			double line_x = nw_x + xAlign * (w - line_w);
			double line_y = nw_y + i*s + a;
			line.draw (g2d, (float) line_x, (float) line_y);

			// Draw bullet if any.
			if (i == 0 && myBullet != null)
				{
				myBullet.draw (g2d, a, line_x, line_y);
				}
			}

		// Restore font.
		g2d.setFont (oldFont);
		}

// Hidden operations.

	/**
	 * Compute the size of this text item's bounding box. The result is stored
	 * in mySize.
	 */
	void computeSize()
		{
		if (mySize == null)
			{
			// Render HTML text.
			renderLines();
			int nlines = myLines.length;

			// Convert each line to a TextLayout object.
			FontRenderContext frc = new FontRenderContext (null, true, true);
			TextLayout[] layoutLines = new TextLayout [nlines];
			for (int i = 0; i < nlines; ++ i)
				{
				layoutLines[i] = new TextLayout (myLines[i].getIterator(), frc);
				}

			// Determine maximum width of any line.
			double w = 0.0;
			for (TextLayout line : layoutLines)
				{
				w = Math.max (w, line.getBounds().getWidth());
				}

			// Determine height of any line.
			LineMetrics metrics = myFont.getLineMetrics ("Eg", frc);
			double a = metrics.getAscent();
			double d = metrics.getDescent();

			// Determine line spacing.
			double s = mySpacing * myFont.getSize2D();

			// Determine height of all lines.
			double h = nlines == 0 ?  0.0 : (nlines-1)*s + a + d;

			// Record bounding box size.
			mySize = new Size (w, h);
			}
		}

	/**
	 * Render this text item's lines. The result is stored in myLines.
	 */
	void renderLines()
		{
		// Early return if text has already been rendered.
		if (myLines != null) return;

		// Split HTML into lines based on newline characters.
		String[] htmlLines;
		if (myText == null)
			{
			htmlLines = new String[] {""};
			}
		else
			{
			htmlLines = myText.split ("\n");
			}
		int nlines = htmlLines.length;

		// Scan each HTML line. Eliminate HTML tags. Convert HTML character
		// entities to the corresponding characters. Store the results in
		// myLines.
		myLines = new AttributedString [nlines];
		for (int i = 0; i < nlines; ++ i)
			{
			String htmlLine = htmlLines[i];
			int nchars = htmlLine.length();
			StringBuilder buf = new StringBuilder();
			int j = 0;
			while (j < nchars)
				{
				char c = htmlLine.charAt (j);
				int k;
				if (c == '<')
					{
					k = htmlLine.indexOf ('>', j);
					if (k == -1) k = nchars - 1;
					j = k;
					}
				else if (c == '&')
					{
					k = htmlLine.indexOf (';', j);
					if (k == -1) k = nchars - 1;
					String entity = htmlLine.substring (j, k + 1);
					if (entity.equals ("&lt;"))
						{
						buf.append ('<');
						}
					else if (entity.equals ("&gt;"))
						{
						buf.append ('>');
						}
					else if (entity.equals ("&amp;"))
						{
						buf.append ('&');
						}
					j = k;
					}
				else
					{
					buf.append (c);
					}
				++ j;
				}
			myLines[i] = new AttributedString (buf.toString());
			}

		// Scan each HTML line. Render HTML tags.
		HashMap<AttributedCharacterIterator.Attribute,Object> attrMap =
			new HashMap<AttributedCharacterIterator.Attribute,Object>();
		Font currentFont = myFont.deriveFont (attrMap);
		for (int i = 0; i < nlines; ++ i)
			{
			String htmlLine = htmlLines[i];
			int nchars = htmlLine.length();
			int lineIndex = 0;
			int j = 0;
			while (j < nchars)
				{
				char c = htmlLine.charAt (j);
				int k;
				if (c == '<')
					{
					k = htmlLine.indexOf ('>', j);
					if (k == -1) k = nchars - 1;
					String tag = htmlLine.substring (j, k + 1);
					if (tag.equalsIgnoreCase ("<B>"))
						{
						attrMap.put
							(TextAttribute.WEIGHT,
							 TextAttribute.WEIGHT_BOLD);
						}
					else if (tag.equalsIgnoreCase ("</B>"))
						{
						attrMap.remove (TextAttribute.WEIGHT);
						}
					else if (tag.equalsIgnoreCase ("<I>"))
						{
						attrMap.put
							(TextAttribute.POSTURE,
							 TextAttribute.POSTURE_OBLIQUE);
						}
					else if (tag.equalsIgnoreCase ("</I>"))
						{
						attrMap.remove (TextAttribute.POSTURE);
						currentFont = myFont.deriveFont (attrMap);
						}
					else if (tag.equalsIgnoreCase ("<SUP>"))
						{
						//attrMap.put
						//	(TextAttribute.SUPERSCRIPT,
						//	 TextAttribute.SUPERSCRIPT_SUPER);
						AffineTransform t = new AffineTransform();
						t.scale (2.0/3.0, 2.0/3.0);
						t.translate (0.0, -myFont.getSize2D()*2.0/3.0);
						attrMap.put (TextAttribute.TRANSFORM, t);
						}
					else if (tag.equalsIgnoreCase ("</SUP>"))
						{
						//attrMap.remove (TextAttribute.SUPERSCRIPT);
						attrMap.remove (TextAttribute.TRANSFORM);
						}
					else if (tag.equalsIgnoreCase ("<SUB>"))
						{
						//attrMap.put
						//	(TextAttribute.SUPERSCRIPT,
						//	 TextAttribute.SUPERSCRIPT_SUB);
						AffineTransform t = new AffineTransform();
						t.scale (2.0/3.0, 2.0/3.0);
						t.translate (0.0, myFont.getSize2D()/3.0);
						attrMap.put (TextAttribute.TRANSFORM, t);
						}
					else if (tag.equalsIgnoreCase ("</SUB>"))
						{
						//attrMap.remove (TextAttribute.SUPERSCRIPT);
						attrMap.remove (TextAttribute.TRANSFORM);
						}
					currentFont = myFont.deriveFont (attrMap);
					j = k;
					}
				else if (c == '&')
					{
					k = htmlLine.indexOf (';', j);
					if (k == -1) k = nchars - 1;
					String entity = htmlLine.substring (j, k + 1);
					if (entity.equals ("&lt;") ||
							entity.equals ("&gt;") ||
							entity.equals ("&amp;"))
						{
						myLines[i].addAttribute
							(TextAttribute.FONT, currentFont,
							 lineIndex, lineIndex + 1);
						++ lineIndex;
						}
					j = k;
					}
				else
					{
					myLines[i].addAttribute
						(TextAttribute.FONT, currentFont,
						 lineIndex, lineIndex + 1);
					++ lineIndex;
					}
				++ j;
				}
			}
		}

	/**
	 * Set the northwest corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 */
	void doNw
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = 0.0;
		}

	/**
	 * Set the north middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 */
	void doN
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = 0.0;
		}

	/**
	 * Set the northeast corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 */
	void doNe
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = 0.0;
		}

	/**
	 * Set the west middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 */
	void doW
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -0.5;
		}

	/**
	 * Set the center point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 */
	void doC
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -0.5;
		}

	/**
	 * Set the east middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 */
	void doE
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -0.5;
		}

	/**
	 * Set the southwest corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 */
	void doSw
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -1.0;
		}

	/**
	 * Set the south middle point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 */
	void doS
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -1.0;
		}

	/**
	 * Set the southeast corner point of this text item's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 */
	void doSe
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -1.0;
		}

	}
