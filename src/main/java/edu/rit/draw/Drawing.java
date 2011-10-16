//******************************************************************************
//
// File:    Drawing.java
// Package: edu.rit.draw
// Unit:    Class edu.rit.draw.Drawing
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

package edu.rit.draw;

import edu.rit.draw.item.ColorFill;
import edu.rit.draw.item.DrawingItem;
import edu.rit.draw.item.Group;
import edu.rit.draw.item.Point;
import edu.rit.draw.item.Size;

import edu.rit.swing.Displayable;
import edu.rit.swing.DisplayableFrame;
import edu.rit.swing.Viewable;

import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.Stroke;

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import java.io.Externalizable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;

import java.util.LinkedList;

/**
 * Class Drawing provides a drawing. A drawing consists of a sequence of
 * {@linkplain DrawingItem}s, which are displayed in order. Thus, a drawing item
 * later in the sequence will obscure a drawing item earlier in the sequence if
 * the items overlap.
 * <P>
 * A drawing contains a {@linkplain ColorFill} object that tells how to color
 * the drawing's background.
 * <P>
 * A drawing contains the size of the rectangular region within which the
 * drawing items are to be displayed. By default, the display region is just
 * large enough to contain all the drawing items, plus a border. Methods are
 * provided to specify a different display region if desired.
 * <P>
 * An instance of class Drawing may be serialized. Class Drawing provides static
 * convenience methods to write a drawing to a file and read a drawing from a
 * file. A drawing file is simply a binary file that contains a serialized
 * drawing object. The {@linkplain View} program can be used to display the
 * drawing file.
 * <P>
 * There is a default drawing. Initially, the default drawing is a new drawing
 * object. The default drawing may be changed, or set to null to signify that
 * there is no default drawing. A {@linkplain DrawingItem}'s <TT>add()</TT> and
 * <TT>addFirst()</TT> methods add the drawing item to the default drawing (if
 * there is one).
 *
 * @author  Alan Kaminsky
 * @version 18-Jul-2008
 */
public class Drawing
	implements Externalizable, Viewable
	{

// Exported enumerations.

	/**
	 * Enumeration Drawing.Alignment specifies how a {@linkplain Drawing}'s
	 * items are to be aligned within the drawing's display region.
	 *
	 * @author  Alan Kaminsky
	 * @version 10-Jul-2006
	 */
	public enum Alignment
		{
		/**
		 * No alignment. Items are drawn at their own locations and are not
		 * aligned to the display region.
		 */
		NONE,
		/**
		 * Align the northwest corner of the drawing items to the northwest
		 * corner point of the display region, inset by the border.
		 */
		NORTHWEST,
		/**
		 * Align the north middle of the drawing items to the north middle point
		 * of the display region, inset by the border.
		 */
		NORTH,
		/**
		 * Align the northeast corner of the drawing items to the northeast
		 * corner point of the display region, inset by the border.
		 */
		NORTHEAST,
		/**
		 * Align the west middle of the drawing items to the west middle point
		 * of the display region, inset by the border.
		 */
		WEST,
		/**
		 * Align the center of the drawing items to the center point of the
		 * display region, inset by the border.
		 */
		CENTER,
		/**
		 * Align the east middle of the drawing items to the east middle point
		 * of the display region, inset by the border.
		 */
		EAST,
		/**
		 * Align the southwest corner of the drawing items to the southwest
		 * corner point of the display region, inset by the border.
		 */
		SOUTHWEST,
		/**
		 * Align the south middle of the drawing items to the south middle point
		 * of the display region, inset by the border.
		 */
		SOUTH,
		/**
		 * Align the southeast corner of the drawing items to the southeast
		 * corner point of the display region, inset by the border.
		 */
		SOUTHEAST,
		}

// Exported constants.

	/**
	 * Signifies that the drawing's display region's size should be determined
	 * automatically based on the drawing items in the drawing. A {@linkplain
	 * Size} of (0,0) is used for this purpose.
	 */
	public static final Size AUTOMATIC_SIZE = new Size (0.0, 0.0);

	/**
	 * The normal display region size (automatic).
	 */
	public static final Size NORMAL_SIZE = AUTOMATIC_SIZE;

	/**
	 * The normal display region border width (2).
	 */
	public static final double NORMAL_BORDER = 2.0;

	/**
	 * The normal display region/drawing item alignment (northwest).
	 */
	public static final Alignment NORMAL_ALIGNMENT = Alignment.NORTHWEST;

	/**
	 * The normal background color (white).
	 */
	public static final ColorFill NORMAL_BACKGROUND = ColorFill.WHITE;

// Hidden data members.

	private static final long serialVersionUID = -3240330399842161140L;

	private LinkedList<DrawingItem> myItems =
		new LinkedList<DrawingItem>();

	private Size mySize = NORMAL_SIZE;
	private double myBorder = NORMAL_BORDER;
	private Alignment myAlignment = NORMAL_ALIGNMENT;
	private ColorFill myBackground = NORMAL_BACKGROUND;

	private Point myItemsNW;
	private Point myItemsSE;

	private String myTitle;

	private static Drawing theDefaultDrawing = new Drawing();

// Exported constructors.

	/**
	 * Construct a new drawing. The drawing contains no drawing items. The
	 * normal size (automatic), normal border (2), normal alignment (northwest),
	 * and normal background color (white) are used.
	 */
	public Drawing()
		{
		}

// Exported operations.

	/**
	 * Returns the default drawing.
	 *
	 * @return  Default drawing, or null if there is no default drawing.
	 */
	public static Drawing defaultDrawing()
		{
		return theDefaultDrawing;
		}

	/**
	 * Set the default drawing.
	 *
	 * @param  theDrawing  Default drawing, or null if there is no default
	 *                     drawing.
	 */
	public static void defaultDrawing
		(Drawing theDrawing)
		{
		theDefaultDrawing = theDrawing;
		}

	/**
	 * Clear this drawing's sequence of drawing items.
	 */
	public void clear()
		{
		myItems.clear();
		}

	/**
	 * Add the given drawing item to the end of this drawing's sequence of
	 * drawing items.
	 *
	 * @param  theItem  Drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public void add
		(DrawingItem theItem)
		{
		if (theItem == null) throw new NullPointerException();
		myItems.add (theItem);
		}

	/**
	 * Add the given drawing item to the beginning of this drawing's sequence of
	 * drawing items.
	 *
	 * @param  theItem  Drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public void addFirst
		(DrawingItem theItem)
		{
		if (theItem == null) throw new NullPointerException();
		myItems.addFirst (theItem);
		}

	/**
	 * Returns this drawing's display region's size. This includes the border if
	 * any. If the return value is equal to <TT>Drawing.AUTOMATIC_SIZE</TT>
	 * (0,0), it signifies that the drawing's display region's size should be
	 * determined automatically based on the drawing items in the drawing.
	 *
	 * @return  Display region size.
	 */
	public Size size()
		{
		return mySize;
		}

	/**
	 * Set this drawing's display region's size. This includes the border if
	 * any. If <TT>theSize</TT> is equal to <TT>Drawing.AUTOMATIC_SIZE</TT>
	 * (0,0), or if either the width or the height of <TT>theSize</TT> is
	 * negative, it signifies that the drawing's display region's size should be
	 * determined automatically based on the drawing items in the drawing.
	 *
	 * @param  theSize  Display region size.
	 *
	 * @return  This drawing object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null.
	 */
	public Drawing size
		(Size theSize)
		{
		if (theSize.width() <= 0.0 || theSize.height() <= 0.0)
			{
			mySize = AUTOMATIC_SIZE;
			}
		else
			{
			mySize = theSize;
			}
		return this;
		}

	/**
	 * Returns this drawing's display region's border.
	 *
	 * @return  Display region border width.
	 */
	public double border()
		{
		return myBorder;
		}

	/**
	 * Set this drawing's display region's border.
	 *
	 * @param  theBorder  Display region border width.
	 *
	 * @return  This drawing object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theBorder</TT> is less than 0.
	 */
	public Drawing border
		(double theBorder)
		{
		if (theBorder < 0.0) throw new IllegalArgumentException();
		myBorder = theBorder;
		return this;
		}

	/**
	 * Returns this drawing's alignment. The alignment determines where this
	 * drawing's drawing items are displayed relative to this drawing's display
	 * region.
	 *
	 * @return  Alignment.
	 */
	public Alignment alignment()
		{
		return myAlignment;
		}

	/**
	 * Set this drawing's alignment. The alignment determines where this
	 * drawing's drawing items are displayed relative to this drawing's display
	 * region.
	 *
	 * @param  theAlignment  Alignment.
	 *
	 * @return  This drawing object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theAlignment</TT> is null.
	 */
	public Drawing alignment
		(Alignment theAlignment)
		{
		if (theAlignment == null) throw new NullPointerException();
		myAlignment = theAlignment;
		return this;
		}

	/**
	 * Returns this drawing's background color fill.
	 *
	 * @return  Background color fill.
	 */
	public ColorFill background()
		{
		return myBackground;
		}

	/**
	 * Set this drawing's background color fill.
	 *
	 * @param  theBackground  Background color fill.
	 *
	 * @return  This drawing object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theBackground</TT> is null.
	 */
	public Drawing background
		(ColorFill theBackground)
		{
		if (theBackground == null) throw new NullPointerException();
		myBackground = theBackground;
		return this;
		}

	/**
	 * Returns a {@linkplain edu.rit.draw.item.Group Group} consisting of the
	 * drawing items in this drawing. In this way, one drawing can be embedded
	 * (as a Grouip) inside another drawing.
	 *
	 * @return  Group.
	 */
	public Group asGroup()
		{
		Group group = new Group();
		for (DrawingItem item : myItems)
			{
			group.append (item);
			}
		return group;
		}

	/**
	 * Write this drawing to the given object output stream.
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
		out.writeInt (myItems.size());
		for (DrawingItem item : myItems)
			{
			out.writeObject (item);
			}
		out.writeObject (mySize);
		out.writeDouble (myBorder);
		out.writeObject (myAlignment);
		out.writeObject (myBackground);
		out.writeObject (myTitle);
		}

	/**
	 * Read this drawing from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if a class needed to deserialize this drawing cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		int n = in.readInt();
		myItems.clear();
		for (int i = 0; i < n; ++ i)
			{
			myItems.add ((DrawingItem) in.readObject());
			}
		mySize = (Size) in.readObject();
		myBorder = in.readDouble();
		myAlignment = (Alignment) in.readObject();
		myBackground = (ColorFill) in.readObject();
		myTitle = (String) in.readObject();
		}

	/**
	 * Write the default drawing to the file with the given name. The default
	 * drawing is written in serialized form to the file. If there is no default
	 * drawing, the <TT>write()</TT> method does nothing.
	 *
	 * @param  theFileName  File name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(String theFileName)
		throws IOException
		{
		if (theDefaultDrawing != null)
			{
			write (theDefaultDrawing, new File (theFileName));
			}
		}

	/**
	 * Write the default drawing to the given file. The default drawing is
	 * written in serialized form to <TT>theFile</TT>. If there is no default
	 * drawing, the <TT>write()</TT> method does nothing.
	 *
	 * @param  theFile  File.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(File theFile)
		throws IOException
		{
		if (theDefaultDrawing != null)
			{
			write (theDefaultDrawing, theFile);
			}
		}

	/**
	 * Write the given drawing to the file with the given name.
	 * <TT>theDrawing</TT> is written in serialized form to the file.
	 *
	 * @param  theDrawing   Drawing.
	 * @param  theFileName  File name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(Drawing theDrawing,
		 String theFileName)
		throws IOException
		{
		write (theDrawing, new File (theFileName));
		}

	/**
	 * Write the given drawing to the given file. <TT>theDrawing</TT> is written
	 * in serialized form to <TT>theFile</TT>.
	 *
	 * @param  theDrawing  Drawing.
	 * @param  theFile     File.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(Drawing theDrawing,
		 File theFile)
		throws IOException
		{
		FileOutputStream fos = null;
		ObjectOutputStream oos = null;

		try
			{
			fos = new FileOutputStream (theFile);
			oos = new ObjectOutputStream (fos);
			oos.writeObject (theDrawing);
			oos.close();
			}

		catch (IOException exc)
			{
			if (fos != null)
				{
				try { fos.close(); } catch (IOException exc2) {}
				}
			throw exc;
			}
		}

	/**
	 * Read a drawing from the file with the given name. The file must contain
	 * one instance of class Drawing in serialized form; for example, as written
	 * by the static <TT>Drawing.write()</TT> method.
	 *
	 * @param  theFileName  File name.
	 *
	 * @return  Drawing.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize the drawing cannot be
	 *     found.
	 */
	public static Drawing read
		(String theFileName)
		throws IOException, ClassNotFoundException
		{
		return read (new File (theFileName));
		}

	/**
	 * Read a drawing from the given file. The file must contain one instance of
	 * class Drawing in serialized form; for example, as written by the static
	 * <TT>Drawing.write()</TT> method.
	 *
	 * @param  theFile  File.
	 *
	 * @return  Drawing.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize the drawing cannot be
	 *     found.
	 */
	public static Drawing read
		(File theFile)
		throws IOException, ClassNotFoundException
		{
		FileInputStream fis = null;
		ObjectInputStream ois = null;
		Drawing result = null;

		try
			{
			fis = new FileInputStream (theFile);
			ois = new ObjectInputStream (fis);
			result = (Drawing) ois.readObject();
			ois.close();
			return result;
			}

		catch (IOException exc)
			{
			if (fis != null)
				{
				try { fis.close(); } catch (IOException exc2) {}
				}
			throw exc;
			}
		}

// Exported operations implemented from interface Drawable.

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		// Compute translation needed to put the drawing items at the right
		// spot.
		computeItemsNWSE();
		Size itemSize = myItemsSE.difference (myItemsNW);
		Size displaySize =
			mySize.equals (AUTOMATIC_SIZE) ?
				itemSize.add (2.0 * myBorder) :
				mySize;
		double dwidth = displaySize.width() - itemSize.width();
		double dheight = displaySize.height() - itemSize.height();
		double dx = 0.0;
		double dy = 0.0;
		switch (myAlignment)
			{
			case NONE:
				break;
			case NORTHWEST:
				dx = myBorder - myItemsNW.x();
				dy = myBorder - myItemsNW.y();
				break;
			case NORTH:
				dx = 0.5 * dwidth - myItemsNW.x();
				dy = myBorder - myItemsNW.y();
				break;
			case NORTHEAST:
				dx = dwidth - myBorder - myItemsNW.x();
				dy = myBorder - myItemsNW.y();
				break;
			case WEST:
				dx = myBorder - myItemsNW.x();
				dy = 0.5 * dheight - myItemsNW.y();
				break;
			case CENTER:
				dx = 0.5 * dwidth - myItemsNW.x();
				dy = 0.5 * dheight - myItemsNW.y();
				break;
			case EAST:
				dx = dwidth - myBorder - myItemsNW.x();
				dy = 0.5 * dheight - myItemsNW.y();
				break;
			case SOUTHWEST:
				dx = myBorder - myItemsNW.x();
				dy = dheight - myBorder - myItemsNW.y();
				break;
			case SOUTH:
				dx = 0.5 * dwidth - myItemsNW.x();
				dy = dheight - myBorder - myItemsNW.y();
				break;
			case SOUTHEAST:
				dx = dwidth - myBorder - myItemsNW.x();
				dy = dheight - myBorder - myItemsNW.y();
				break;
			}

		// Save graphics context's state.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		AffineTransform oldTransform = g2d.getTransform();
		Object oldAntialiasing =
			g2d.getRenderingHint
				(RenderingHints.KEY_ANTIALIASING);
		Object oldTextAntialiasing =
			g2d.getRenderingHint
				(RenderingHints.KEY_TEXT_ANTIALIASING);
		Object oldFractionalMetrics =
			g2d.getRenderingHint
				(RenderingHints.KEY_FRACTIONALMETRICS);

		// Turn on antialiasing.
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setRenderingHint
			(RenderingHints.KEY_TEXT_ANTIALIASING,
			 RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		g2d.setRenderingHint
			(RenderingHints.KEY_FRACTIONALMETRICS,
			 RenderingHints.VALUE_FRACTIONALMETRICS_ON);

		// Iterate over all drawing items.
		for (DrawingItem item : myItems)
			{
			// Restore graphics context.
			g2d.setStroke (oldStroke);
			g2d.setPaint (oldPaint);
			g2d.setTransform (oldTransform);
			g2d.translate (dx, dy);

			// Draw item.
			item.draw (g2d);
			}

		// Restore graphics context's state.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		g2d.setTransform (oldTransform);
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 oldAntialiasing);
		g2d.setRenderingHint
			(RenderingHints.KEY_TEXT_ANTIALIASING,
			 oldTextAntialiasing);
		g2d.setRenderingHint
			(RenderingHints.KEY_FRACTIONALMETRICS,
			 oldFractionalMetrics);
		}

// Exported operations implemented from interface Displayable.

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 */
	public Rectangle2D getBoundingBox()
		{
		computeItemsNWSE();
		Size displaySize =
			mySize.equals (AUTOMATIC_SIZE) ?
				myItemsSE.difference (myItemsNW) .add (2.0 * myBorder) :
				mySize;
		return new Rectangle2D.Double
			(0.0, 0.0, displaySize.width(), displaySize.height());
		}

	/**
	 * Returns this displayable object's background paint.
	 */
	public Paint getBackgroundPaint()
		{
		return myBackground.color();
		}

// Exported operations implemented from interface Viewable.

	/**
	 * Get a displayable frame in which to view this viewable object. Initially,
	 * the returned frame is displaying this viewable object.
	 *
	 * @return  Displayable frame.
	 */
	public DisplayableFrame getFrame()
		{
		return new DrawingFrame (getTitle(), this);
		}

	/**
	 * Get the title for the frame used to view this viewable object. If the
	 * title is null, a default title is used.
	 *
	 * @return  Title.
	 */
	public String getTitle()
		{
		return myTitle;
		}

	/**
	 * Set the title for the frame used to view this viewable object. If the
	 * title is null, a default title is used.
	 *
	 * @param  theTitle  Title.
	 */
	public void setTitle
		(String theTitle)
		{
		myTitle = theTitle;
		}

// Hidden operations.

	/**
	 * Compute the furthest northwest and furthest southeast points of all this
	 * drawing's drawing items. Store the results in myItemsNW and myItemsSE.
	 */
	private void computeItemsNWSE()
		{
		Rectangle2D bbox = new Rectangle2D.Double();
		for (DrawingItem item : myItems)
			{
			bbox = bbox.createUnion (item.boundingBox());
			}
		myItemsNW = new Point (bbox.getX(), bbox.getY());
		myItemsSE = new Point
			(bbox.getX()+bbox.getWidth(), bbox.getY()+bbox.getHeight());
		}

	}
