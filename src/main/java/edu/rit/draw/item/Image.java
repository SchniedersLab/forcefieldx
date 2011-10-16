//******************************************************************************
//
// File:    Image.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Image
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

import java.awt.Graphics2D;

import java.awt.geom.AffineTransform;

import java.awt.image.BufferedImage;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.Externalizable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import javax.imageio.ImageIO;

/**
 * Class Image provides a {@linkplain DrawingItem} containing an image. The
 * image is read from a file or input stream. Image formats supported by the
 * javax.imageio package are supported; this typically includes GIF, JPEG, and
 * PNG images.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Image
	extends DrawingItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal scale factor (1).
	 */
	public static final double NORMAL_SCALE = 1.0;

// Hidden data members.

	private static final long serialVersionUID = 6921389113129571217L;

	// Default attributes.
	static double theDefaultScale = NORMAL_SCALE;

	// The contents of the image file stored in a byte array.
	byte[] myImageContents;

	// Buffered image read from the contents of the image file.
	BufferedImage myBufferedImage;

	// Attributes.
	double myScale = theDefaultScale;

	// Coordinates of most recently specified corner.
	double x;
	double y;

	// Factors for going from specified corner to northwest corner.
	double xFactor;
	double yFactor;

	// Size. If null, the size must be recomputed.
	Size mySize;

// Exported constructors.

	/**
	 * Construct a new image item. The image item's image is empty. The image
	 * item's northwest corner is located at (0,0). The image item's attributes
	 * have the default values.
	 */
	public Image()
		{
		super();
		}

	/**
	 * Construct a new image item with the same image and attributes as the
	 * given image item.
	 *
	 * @param  theItem  Image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Image
		(Image theItem)
		{
		super (theItem);
		this.myImageContents = theItem.myImageContents;
		this.myScale = theItem.myScale;
		this.x = theItem.x;
		this.y = theItem.y;
		this.xFactor = theItem.xFactor;
		this.yFactor = theItem.yFactor;
		this.mySize = theItem.mySize;
		}

// Exported operations.

	/**
	 * Returns the default scale factor. Images are displayed at a scaled size,
	 * one image pixel = (scale factor) display units.
	 *
	 * @return  Default scale factor.
	 */
	public static double defaultScale()
		{
		return theDefaultScale;
		}

	/**
	 * Set the default scale factor. Images are displayed at a scaled size, one
	 * image pixel = (scale factor) display units.
	 *
	 * @param  theScale  Default scale factor.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theScale</TT> is less than 0.
	 */
	public static void defaultScale
		(double theScale)
		{
		if (theScale < 0.0) throw new IllegalArgumentException();
		theDefaultScale = theScale;
		}

	/**
	 * Returns this image item's image.
	 *
	 * @return  Image.
	 */
	public BufferedImage image()
		{
		return myBufferedImage;
		}

	/**
	 * Read this image item's image from the file with the given name.
	 *
	 * @param  theFileName  File name.
	 *
	 * @return  This image item.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Image image
		(String theFileName)
		throws IOException
		{
		readImageContents (new FileInputStream (theFileName));
		getBufferedImage();
		return this;
		}

	/**
	 * Read this image item's image from the given file.
	 *
	 * @param  theFile  File.
	 *
	 * @return  This image item.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Image image
		(File theFile)
		throws IOException
		{
		readImageContents (new FileInputStream (theFile));
		getBufferedImage();
		return this;
		}

	/**
	 * Read this image item's image from the given input stream. The entire
	 * input stream is read until EOF is encountered, then the input stream is
	 * closed.
	 *
	 * @param  in  Input stream.
	 *
	 * @return  This image item.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Image image
		(InputStream in)
		throws IOException
		{
		readImageContents (in);
		getBufferedImage();
		return this;
		}

	/**
	 * Returns this image item's scale factor. This image is displayed at a
	 * scaled size, one image pixel = (scale factor) display units.
	 *
	 * @return  Scale factor.
	 */
	public double scale()
		{
		return myScale;
		}

	/**
	 * Set this image item's scale factor. This image is displayed at a scaled
	 * size, one image pixel = (scale factor) display units.
	 *
	 * @param  theScale  Scale factor.
	 *
	 * @return  This image.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theScale</TT> is less than 0.
	 */
	public Image scale
		(double theScale)
		{
		if (theScale < 0.0) throw new IllegalArgumentException();
		myScale = theScale;
		return this;
		}

	/**
	 * Returns the width of this image item's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		computeSize();
		return mySize.width;
		}

	/**
	 * Returns the height of this image item's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		computeSize();
		return mySize.height;
		}

	/**
	 * Returns the northwest corner point of this image item's bounding box.
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
	 * Set the northwest corner point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This image item.
	 */
	public Image nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this image item's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This image item.
	 */
	public Image n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this image item's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This image item.
	 */
	public Image ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this image item's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This image item.
	 */
	public Image w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this image item's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This image item.
	 */
	public Image c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this image item's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This image item.
	 */
	public Image e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this image item's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This image item.
	 */
	public Image sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this image item's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This image item.
	 */
	public Image s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this image item's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this image item's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This image item.
	 */
	public Image se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this image item's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Image se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Add this image item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Image add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this image item to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Image add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this image item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Image addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this image item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This image item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Image addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this image item to the given object output stream.
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
		out.writeInt (myImageContents.length);
		out.write (myImageContents);
		out.writeDouble (myScale);
		out.writeDouble (x);
		out.writeDouble (y);
		out.writeDouble (xFactor);
		out.writeDouble (yFactor);
		out.writeObject (mySize);
		}

	/**
	 * Read this image item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this image item cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		int n = in.readInt();
		myImageContents = new byte [n];
		in.readFully (myImageContents);
		getBufferedImage();
		myScale = in.readDouble();
		x = in.readDouble();
		y = in.readDouble();
		xFactor = in.readDouble();
		yFactor = in.readDouble();
		mySize = (Size) in.readObject();
		}

	/**
	 * Draw this image item in the given graphics context. This method is
	 * allowed to change the graphics context's paint, stroke, and transform,
	 * and it doesn't have to change them back.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		super.draw (g2d);

		// Determine northwest corner point.
		computeSize();
		double nw_x = x + xFactor * mySize.width();
		double nw_y = y + yFactor * mySize.height();

		// Set up drawing transform.
		AffineTransform transform = new AffineTransform();
		transform.translate (nw_x, nw_y);
		transform.scale (myScale, myScale);
		transform.translate 
			(- myBufferedImage.getMinX(),
			 - myBufferedImage.getMinY());

		// Draw image.
		g2d.drawRenderedImage (myBufferedImage, transform);
		}

// Hidden operations.

	/**
	 * Read the given input stream and store its contents in the field
	 * <TT>myImageContents</TT>. The entire input stream is read until EOF is
	 * encountered, then the input stream is closed.
	 *
	 * @param  in  Input stream from which to read the image.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void readImageContents
		(InputStream in)
		throws IOException
		{
		try
			{
			int b;
			BufferedInputStream bis = new BufferedInputStream (in);
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			while ((b = bis.read()) != -1)
				{
				baos.write (b);
				}
			myImageContents = baos.toByteArray();
			myBufferedImage = null;
			mySize = null;
			}
		finally
			{
			try { in.close(); } catch (IOException exc2) {}
			}
		}

	/**
	 * Get a buffered image from the bytes stored in the field
	 * <TT>myImageContents</TT>, and store the buffered image in the field
	 * <TT>myBufferedImage</TT>.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void getBufferedImage()
		throws IOException
		{
		try
			{
			myBufferedImage =
				ImageIO.read (new ByteArrayInputStream (myImageContents));
			mySize = null;
			}
		catch (IOException exc)
			{
			myImageContents = null;
			myBufferedImage = null;
			mySize = null;
			throw exc;
			}
		}

	/**
	 * Compute the size of this image item's bounding box. The result is stored
	 * in mySize.
	 */
	void computeSize()
		{
		if (mySize != null)
			{
			}
		else if (myBufferedImage != null)
			{
			mySize =
				new Size
					(myScale * myBufferedImage.getWidth(),
					 myScale * myBufferedImage.getHeight());
			}
		else
			{
			mySize = new Size (0.0, 0.0);
			}
		}

	/**
	 * Set the northwest corner point of this image item's bounding box.
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
	 * Set the north middle point of this image item's bounding box.
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
	 * Set the northeast corner point of this image item's bounding box.
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
	 * Set the west middle point of this image item's bounding box.
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
	 * Set the center point of this image item's bounding box.
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
	 * Set the east middle point of this image item's bounding box.
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
	 * Set the southwest corner point of this image item's bounding box.
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
	 * Set the south middle point of this image item's bounding box.
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
	 * Set the southeast corner point of this image item's bounding box.
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
