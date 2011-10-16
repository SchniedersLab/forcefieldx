//******************************************************************************
//
// File:    Label.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.Label
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.numeric.plot.impl;

import edu.rit.swing.Drawable;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;

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

import java.text.AttributedCharacterIterator;
import java.text.AttributedString;

import java.util.HashMap;

/**
 * Class Label provides an object for drawing arbitrary text on a plot. A label
 * has a <I>location,</I> a certain pair of coordinates (<I>x,y</I>) on the
 * plot. A label can be drawn at any one of nine positions relative to its
 * location: above left, above, above right, left, center, right, below left,
 * below, or below right. A label can be offset a certain distance from its
 * location. A label can be not rotated, rotated 90 degrees left, rotated 90
 * degrees right, or rotated 180 degrees. A label's location can be specified
 * either in plot coordinates or in pixel coordinates.
 * <P>
 * Class Label supports a subset of HTML for specifying attributes of
 * characters. The following HTML tags may be embedded in a label's text:
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
 * @version 11-Dec-2007
 */
public class Label
	implements Drawable, Externalizable
	{

// Exported constants.

	/**
	 * The label is positioned above and to the left of its location.
	 */
	public static final int ABOVE_LEFT = 0;

	/**
	 * The label is positioned above its location.
	 */
	public static final int ABOVE = 1;

	/**
	 * The label is positioned above and to the right of its location.
	 */
	public static final int ABOVE_RIGHT = 2;

	/**
	 * The label is positioned to the left of its location.
	 */
	public static final int LEFT = 3;

	/**
	 * The label is centered on its location.
	 */
	public static final int CENTER = 4;

	/**
	 * The label is positioned to the right of its location.
	 */
	public static final int RIGHT = 5;

	/**
	 * The label is positioned below and to the left of its location.
	 */
	public static final int BELOW_LEFT = 6;

	/**
	 * The label is positioned below its location.
	 */
	public static final int BELOW = 7;

	/**
	 * The label is positioned below and to the right of its location.
	 */
	public static final int BELOW_RIGHT = 8;

	/**
	 * The label is rotated 90 degrees left (counterclockwise).
	 */
	public static final int ROTATE_LEFT = 9;

	/**
	 * The label is rotated 90 degrees right (clockwise).
	 */
	public static final int ROTATE_RIGHT = 18;

	/**
	 * The label is rotated 180 degrees around.
	 */
	public static final int ROTATE_AROUND = 27;

	/**
	 * The label's location is specified in pixel coordinates relative to the
	 * lower left corner of the plot area, rather than in plot coordinates.
	 */
	public static final int PIXEL_COORDINATES = 36;

	/**
	 * The default label font (sanserif, plain, 12 point).
	 */
	public static final Font DEFAULT_FONT =
		new Font ("sanserif", Font.PLAIN, 12);

	/**
	 * The default label color (black).
	 */
	public static final Color DEFAULT_PAINT = Color.black;

// Hidden data members.

	private static final long serialVersionUID = -560274063744339038L;

	/**
	 * This label's text.
	 */
	protected String myText;

	/**
	 * X coordinate of this label's location.
	 */
	protected double myLocationX;

	/**
	 * Y coordinate of this label's location.
	 */
	protected double myLocationY;

	/**
	 * This label's positioning, ABOVE_LEFT through BELOW_RIGHT.
	 */
	protected int myPositioning;

	/**
	 * This label's offset from its location on the display.
	 */
	protected double myOffset;

	/**
	 * This label's font.
	 */
	protected Font myFont;

	/**
	 * The color for drawing this label on the display.
	 */
	protected Color myColor;

	/**
	 * The color for drawing this label's background on the display.
	 */
	protected Color myBackground;

// Exported constructors.

	/**
	 * Construct a new label. This constructor is intended for use only by
	 * object deserialization.
	 */
	public Label()
		{
		}

	/**
	 * Construct a new label with the given text and location. The label is
	 * centered on its location and is not rotated. The default font and color
	 * are used. The label has no background.
	 *
	 * @param  text  Label text.
	 * @param  x     X coordinate of this label's location.
	 * @param  y     Y coordinate of this label's location.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null.
	 */
	public Label
		(String text,
		 double x,
		 double y)
		{
		this (text, x, y, CENTER, 0, DEFAULT_FONT, DEFAULT_PAINT, null);
		}

	/**
	 * Construct a new label with the given text, location, and positioning. The
	 * default font and color are used. The label has no background.
	 *
	 * @param  text      Label text.
	 * @param  x         X coordinate of this label's location.
	 * @param  y         Y coordinate of this label's location.
	 * @param  position  This label's positioning, <TT>ABOVE_LEFT</TT> through
	 *                   <TT>BELOW_RIGHT</TT>. Optionally, <TT>ROTATE_LEFT</TT>
	 *                   through <TT>ROTATE_AROUND</TT> may be added; e.g.,
	 *                   <TT>LEFT+ROTATE_LEFT</TT>. Optionally,
	 *                   <TT>PIXEL_COORDINATES</TT> may be added.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>position</TT> is invalid.
	 */
	public Label
		(String text,
		 double x,
		 double y,
		 int position)
		{
		this (text, x, y, position, 0, DEFAULT_FONT, DEFAULT_PAINT, null);
		}

	/**
	 * Construct a new label with the given text, location, positioning, and
	 * offset from its location. The default font and color are used. The label
	 * has no background. If the label's positioning is <TT>CENTER</TT>, the
	 * offset is not used.
	 *
	 * @param  text      Label text.
	 * @param  x         X coordinate of this label's location.
	 * @param  y         Y coordinate of this label's location.
	 * @param  position  This label's positioning, <TT>ABOVE_LEFT</TT> through
	 *                   <TT>BELOW_RIGHT</TT>. Optionally, <TT>ROTATE_LEFT</TT>
	 *                   through <TT>ROTATE_AROUND</TT> may be added; e.g.,
	 *                   <TT>LEFT+ROTATE_LEFT</TT>. Optionally,
	 *                   <TT>PIXEL_COORDINATES</TT> may be added.
	 * @param  offset    Offset from this label's location on the display (in
	 *                   display units).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>position</TT> is invalid. Thrown
	 *     if <TT>offset</TT> is less than 0.
	 */
	public Label
		(String text,
		 double x,
		 double y,
		 int position,
		 double offset)
		{
		this (text, x, y, position, offset, DEFAULT_FONT, DEFAULT_PAINT, null);
		}

	/**
	 * Construct a new label with the given text, location, positioning, offset
	 * from its location, and font. The default color is used. The label has no
	 * background. If the label's positioning is <TT>CENTER</TT>, the offset is
	 * not used.
	 *
	 * @param  text      Label text.
	 * @param  x         X coordinate of this label's location.
	 * @param  y         Y coordinate of this label's location.
	 * @param  position  This label's positioning, <TT>ABOVE_LEFT</TT> through
	 *                   <TT>BELOW_RIGHT</TT>. Optionally, <TT>ROTATE_LEFT</TT>
	 *                   through <TT>ROTATE_AROUND</TT> may be added; e.g.,
	 *                   <TT>LEFT+ROTATE_LEFT</TT>. Optionally,
	 *                   <TT>PIXEL_COORDINATES</TT> may be added.
	 * @param  offset    Offset from this label's location on the display (in
	 *                   display units).
	 * @param  font      Font for drawing this label.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null or
	 *     <TT>font</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>position</TT> is invalid. Thrown
	 *     if <TT>offset</TT> is less than 0.
	 */
	public Label
		(String text,
		 double x,
		 double y,
		 int position,
		 double offset,
		 Font font)
		{
		this (text, x, y, position, offset, font, DEFAULT_PAINT, null);
		}

	/**
	 * Construct a new label with the given text, location, positioning, offset
	 * from its location, font, and color. The label has no background. If the
	 * label's positioning is <TT>CENTER</TT>, the offset is not used.
	 *
	 * @param  text      Label text.
	 * @param  x         X coordinate of this label's location.
	 * @param  y         Y coordinate of this label's location.
	 * @param  position  This label's positioning, <TT>ABOVE_LEFT</TT> through
	 *                   <TT>BELOW_RIGHT</TT>. Optionally, <TT>ROTATE_LEFT</TT>
	 *                   through <TT>ROTATE_AROUND</TT> may be added; e.g.,
	 *                   <TT>LEFT+ROTATE_LEFT</TT>. Optionally,
	 *                   <TT>PIXEL_COORDINATES</TT> may be added.
	 * @param  offset    Offset from this label's location on the display (in
	 *                   display units).
	 * @param  font      Font for drawing this label.
	 * @param  color     Color for drawing this label.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null, <TT>font</TT>
	 *     is null, or <TT>color</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>position</TT> is invalid. Thrown
	 *     if <TT>offset</TT> is less than 0.
	 */
	public Label
		(String text,
		 double x,
		 double y,
		 int position,
		 double offset,
		 Font font,
		 Color color)
		{
		this (text, x, y, position, offset, font, color, null);
		}

	/**
	 * Construct a new label with the given text, location, positioning, offset
	 * from its location, font, color, and background. If <TT>bg</TT> is null,
	 * the label has no background. If the label's positioning is
	 * <TT>CENTER</TT>, the offset is not used.
	 *
	 * @param  text      Label text.
	 * @param  x         X coordinate of this label's location.
	 * @param  y         Y coordinate of this label's location.
	 * @param  position  This label's positioning, <TT>ABOVE_LEFT</TT> through
	 *                   <TT>BELOW_RIGHT</TT>. Optionally, <TT>ROTATE_LEFT</TT>
	 *                   through <TT>ROTATE_AROUND</TT> may be added; e.g.,
	 *                   <TT>LEFT+ROTATE_LEFT</TT>. Optionally,
	 *                   <TT>PIXEL_COORDINATES</TT> may be added.
	 * @param  offset    Offset from this label's location on the display (in
	 *                   display units).
	 * @param  font      Font for drawing this label.
	 * @param  color     Color for drawing this label.
	 * @param  bg        Color for drawing this label's background.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>text</TT> is null, <TT>font</TT>
	 *     is null, or <TT>color</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>position</TT> is invalid. Thrown
	 *     if <TT>offset</TT> is less than 0.
	 */
	public Label
		(String text,
		 double x,
		 double y,
		 int position,
		 double offset,
		 Font font,
		 Color color,
		 Color bg)
		{
		if (text == null || font == null || color == null)
			{
			throw new NullPointerException();
			}
		if
			(ABOVE_LEFT > position ||
			 position > BELOW_RIGHT+ROTATE_AROUND+PIXEL_COORDINATES ||
			 offset < 0.0)
			{
			throw new IllegalArgumentException();
			}
		myText = text;
		myLocationX = x;
		myLocationY = y;
		myPositioning = position;
		myOffset = offset;
		myFont = font;
		myColor = color;
		myBackground = bg;
		}

// Exported operations.

	/**
	 * Returns this label's text.
	 */
	public String getText()
		{
		return myText;
		}

	/**
	 * Returns the X coordinate of this label's location.
	 */
	public double getLocationX()
		{
		return myLocationX;
		}

	/**
	 * Returns the Y coordinate of this label's location.
	 */
	public double getLocationY()
		{
		return myLocationY;
		}

	/**
	 * Returns this label's positioning. The return value is <TT>ABOVE_LEFT</TT>
	 * through <TT>BELOW_RIGHT</TT>. If the label is rotated, one of the values
	 * <TT>ROTATE_LEFT</TT> through <TT>ROTATE_AROUND</TT> is added. If the
	 * pixel's location is specified in pixel coordinates rather than plot
	 * coordinates, <TT>PIXEL_COORDINATES</TT> is added.
	 */
	public int getPositioning()
		{
		return myPositioning;
		}

	/**
	 * Returns this label's offset from its location on the display (in display
	 * units).
	 */
	public double getOffset()
		{
		return myOffset;
		}

	/**
	 * Returns this label's font.
	 */
	public Font getFont()
		{
		return myFont;
		}

	/**
	 * Returns the color for drawing this label on the display.
	 */
	public Color getColor()
		{
		return myColor;
		}

	/**
	 * Returns the color for drawing this label's background on the display. If
	 * this label has no background, null is returned.
	 */
	public Color getBackground()
		{
		return myBackground;
		}

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 * <P>
	 * The label is drawn at position (0,0). It is the caller's responsibility
	 * to transform the graphics context so the label ends up at the correct
	 * (<I>x,y</I>) coordinates.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		int j;

		// Save graphics context's state.
		Paint oldPaint = g2d.getPaint();
		Font oldFont = g2d.getFont();
		AffineTransform oldTransform = g2d.getTransform();

		// Scan the HTML text. Eliminate HTML tags. Convert HTML character
		// entities to the corresponding characters. Store the results in
		// attrString.
		int nchars = myText.length();
		StringBuilder buf = new StringBuilder();
		j = 0;
		while (j < nchars)
			{
			char c = myText.charAt (j);
			int k;
			if (c == '<')
				{
				k = myText.indexOf ('>', j);
				if (k == -1) k = nchars - 1;
				j = k;
				}
			else if (c == '&')
				{
				k = myText.indexOf (';', j);
				if (k == -1) k = nchars - 1;
				String entity = myText.substring (j, k + 1);
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
		AttributedString attrString = new AttributedString (buf.toString());

		// Scan the HTML text. Render HTML tags.
		HashMap<AttributedCharacterIterator.Attribute,Object> attrMap =
			new HashMap<AttributedCharacterIterator.Attribute,Object>();
		Font currentFont = myFont.deriveFont (attrMap);
		int lineIndex = 0;
		j = 0;
		while (j < nchars)
			{
			char c = myText.charAt (j);
			int k;
			if (c == '<')
				{
				k = myText.indexOf ('>', j);
				if (k == -1) k = nchars - 1;
				String tag = myText.substring (j, k + 1);
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
				k = myText.indexOf (';', j);
				if (k == -1) k = nchars - 1;
				String entity = myText.substring (j, k + 1);
				if (entity.equals ("&lt;") ||
						entity.equals ("&gt;") ||
						entity.equals ("&amp;"))
					{
					attrString.addAttribute
						(TextAttribute.FONT, currentFont,
						 lineIndex, lineIndex + 1);
					++ lineIndex;
					}
				j = k;
				}
			else
				{
				attrString.addAttribute
					(TextAttribute.FONT, currentFont,
					 lineIndex, lineIndex + 1);
				++ lineIndex;
				}
			++ j;
			}

		// Set up font and get metrics.
		g2d.setFont (myFont);
		FontRenderContext frc = g2d.getFontRenderContext();
		TextLayout layoutLine = new TextLayout (attrString.getIterator(), frc);
		Rectangle2D bounds = layoutLine.getBounds();
		double w = bounds.getWidth();
		LineMetrics metrics = myFont.getLineMetrics ("Eg", frc);
		double a = metrics.getAscent();
		double d = metrics.getDescent();
		double h = a + d;

		// Get positioning and rotation.
		int positioning = myPositioning % 9;
		int rotation = ((myPositioning / 9) % 4) * ROTATE_LEFT;

		// Get width and height of label after rotation.
		double rotw;
		double roth;
		switch (rotation)
			{
			case ROTATE_LEFT:
			case ROTATE_RIGHT:
				rotw = h;
				roth = w;
				break;
			default:
				rotw = w;
				roth = h;
				break;
			}

		// Translate rotated label's center to the proper spot.
		switch (positioning)
			{
			case ABOVE_LEFT:
				g2d.translate (-rotw/2-myOffset, -roth/2-myOffset);
				break;
			case ABOVE:
				g2d.translate (0, -roth/2-myOffset);
				break;
			case ABOVE_RIGHT:
				g2d.translate (+rotw/2+myOffset, -roth/2-myOffset);
				break;
			case LEFT:
				g2d.translate (-rotw/2-myOffset, 0);
				break;
			case RIGHT:
				g2d.translate (+rotw/2+myOffset, 0);
				break;
			case BELOW_LEFT:
				g2d.translate (-rotw/2-myOffset, +roth/2+myOffset);
				break;
			case BELOW:
				g2d.translate (0, +roth/2+myOffset);
				break;
			case BELOW_RIGHT:
				g2d.translate (+rotw/2+myOffset, +roth/2+myOffset);
				break;
			}

		// Rotate label around its center.
		switch (rotation)
			{
			case ROTATE_LEFT:
				g2d.rotate (-Math.PI/2);
				break;
			case ROTATE_RIGHT:
				g2d.rotate (+Math.PI/2);
				break;
			case ROTATE_AROUND:
				g2d.rotate (Math.PI);
				break;
			}

		// Paint background rectangle if necessary.
		if (myBackground != null)
			{
			g2d.setPaint (myBackground);
			g2d.fill (new Rectangle2D.Double (-w/2, -h/2, w, h));
			}

		// Draw label text.
		g2d.setPaint (myColor);
		layoutLine.draw (g2d, (float)(-w/2), (float)(-h/2+a));

		// Restore graphics context's state.
		g2d.setPaint (oldPaint);
		g2d.setFont (oldFont);
		g2d.setTransform (oldTransform);
		}

	/**
	 * Write this label to the given object output stream.
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
		out.writeUTF (myText);
		out.writeDouble (myLocationX);
		out.writeDouble (myLocationY);
		out.writeInt (myPositioning);
		out.writeDouble (myOffset);
		out.writeObject (myFont);
		out.writeObject (myColor);
		out.writeObject (myBackground);
		}

	/**
	 * Read this label from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if a class needed to deserialize this label could not be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		myText = in.readUTF();
		myLocationX = in.readDouble();
		myLocationY = in.readDouble();
		myPositioning = in.readInt();
		myOffset = in.readDouble();
		myFont = (Font) in.readObject();
		myColor = (Color) in.readObject();
		myBackground = (Color) in.readObject();
		}

	}
