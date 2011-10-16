//******************************************************************************
//
// File:    DisplayableIO.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DisplayableIO
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

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import java.awt.geom.Rectangle2D;

import java.awt.image.BufferedImage;

import java.awt.print.PageFormat;
import java.awt.print.Printable;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.imageio.ImageIO;

import javax.print.Doc;
import javax.print.DocFlavor;
import javax.print.PrintException;
import javax.print.SimpleDoc;
import javax.print.StreamPrintService;
import javax.print.StreamPrintServiceFactory;

/**
 * Class DisplayableIO provides static methods for writing displayable objects
 * as image files. You can write PNG files and PostScript files.
 *
 * @author  Alan Kaminsky
 * @version 12-Mar-2003
 */
public class DisplayableIO
	{

// Prevent construction.

	private DisplayableIO()
		{
		}

// Exported operations.

	/**
	 * Write the given displayable object to a PNG file in the form of a bilevel
	 * (black and white) image. The default scale factor (one display unit = one
	 * pixel) is used.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeBilevelPNGFile
		(Displayable theDisplayable,
		 File theFile)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 1.0,
			 BufferedImage.TYPE_BYTE_BINARY,
			 false);
		}

	/**
	 * Write the given displayable object to a PNG file in the form of a bilevel
	 * (black and white) image with the given scale factor. One display unit =
	 * <TT>theScale</TT> pixels.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 * @param  theScale        Scale factor for converting display units to
	 *                         pixels.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeBilevelPNGFile
		(Displayable theDisplayable,
		 File theFile,
		 double theScale)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 theScale,
			 BufferedImage.TYPE_BYTE_BINARY,
			 false);
		}

	/**
	 * Write the given displayable object to a PNG file in the form of a
	 * grayscale image. The default scale factor (one display unit = one pixel)
	 * is used.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeGrayscalePNGFile
		(Displayable theDisplayable,
		 File theFile)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 1.0,
			 BufferedImage.TYPE_BYTE_GRAY,
			 true);
		}

	/**
	 * Write the given displayable object to a PNG file in the form of a
	 * grayscale image with the given scale factor. One display unit =
	 * <TT>theScale</TT> pixels.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 * @param  theScale        Scale factor for converting display units to
	 *                         pixels.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeGrayscalePNGFile
		(Displayable theDisplayable,
		 File theFile,
		 double theScale)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 theScale,
			 BufferedImage.TYPE_BYTE_GRAY,
			 true);
		}

	/**
	 * Write the given displayable object to a PNG file in the form of a 24-bit
	 * color image. The default scale factor (one display unit = one pixel) is
	 * used.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeColorPNGFile
		(Displayable theDisplayable,
		 File theFile)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 1.0,
			 BufferedImage.TYPE_INT_RGB,
			 true);
		}

	/**
	 * Write the given displayable object to a PNG file in the form of a 24-bit
	 * color image with the given scale factor. One display unit =
	 * <TT>theScale</TT> pixels.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 * @param  theScale        Scale factor for converting display units to
	 *                         pixels.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeColorPNGFile
		(Displayable theDisplayable,
		 File theFile,
		 double theScale)
		throws IOException
		{
		writePNGFile
			(theDisplayable,
			 theFile,
			 theScale,
			 BufferedImage.TYPE_INT_RGB,
			 true);
		}

	/**
	 * Write the given displayable object to a PNG file.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         Image file.
	 * @param  theScale        Scale factor for converting display units to
	 *                         pixels.
	 * @param  theImageType    Type of image.
	 * @param  antialias       True to turn antialiasing on, false to turn
	 *                         antialiasing off.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static void writePNGFile
		(Displayable theDisplayable,
		 File theFile,
		 double theScale,
		 int theImageType,
		 boolean antialias)
		throws IOException
		{
		// Get displayable object's bounds.
		Rectangle2D bounds = theDisplayable.getBoundingBox();
		int w = (int) (bounds.getWidth() * theScale + 0.5);
		int h = (int) (bounds.getHeight() * theScale + 0.5);

		// Create a buffered image of the right size and type.
		BufferedImage bi = new BufferedImage (w, h, theImageType);

		// Get graphics context to write into the buffered image.
		Graphics2D g2d = bi.createGraphics();
		if (antialias)
			{
			g2d.setRenderingHint
				(RenderingHints.KEY_ANTIALIASING,
				 RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setRenderingHint
				(RenderingHints.KEY_TEXT_ANTIALIASING,
				 RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			}

		// Fill with background paint.
		g2d.setPaint (theDisplayable.getBackgroundPaint());
		g2d.fill (new Rectangle2D.Double (0, 0, w, h));

		// Scale graphics context by scale factor.
		g2d.scale (theScale, theScale);

		// Shift graphics context's origin to bounding box's origin.
		g2d.translate
			(- bounds.getX() * theScale,
			 - bounds.getY() * theScale);

		// Draw the displayable object.
		theDisplayable.draw (g2d);

		// Write the PNG file.
		ImageIO.write (bi, "png", theFile);
		}

	/**
	 * Write the given displayable object to a PostScript file.
	 *
	 * @param  theDisplayable  Displayable object to write.
	 * @param  theFile         PostScript file.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writePostScriptFile
		(Displayable theDisplayable,
		 File theFile)
		throws IOException
		{
		// Find a suitable stream print service factory.
		DocFlavor theFlavor = DocFlavor.SERVICE_FORMATTED.PRINTABLE;
		StreamPrintServiceFactory[] factories =
			StreamPrintServiceFactory.lookupStreamPrintServiceFactories
				(theFlavor,
				 "application/postscript");
		if (factories.length == 0)
			{
			throw new IOException ("Cannot find PostScript generator");
			}

		// Get displayable object's bounds.
		final Displayable disp = theDisplayable;
		final Rectangle2D bounds = theDisplayable.getBoundingBox();
		final double w = bounds.getWidth();
		final double h = bounds.getHeight();

		// Create a printable object to print the displayable object.
		Printable thePrintData = new Printable()
			{
			public int print
				(Graphics g,
				 PageFormat pageFormat,
				 int pageIndex)
				{
				// Only print one page.
				if (pageIndex != 0)
					{
					return Printable.NO_SUCH_PAGE;
					}

				// Get graphics context to write into the page.
				Graphics2D g2d = (Graphics2D) g;

				// Shift graphics context's origin to top left of imageable
				// area.
				g2d.translate
					(pageFormat.getImageableX(),
					 pageFormat.getImageableY());

				// Scale so drawable object fits within imageable area.
				double xScale = pageFormat.getImageableWidth() / w;
				double yScale = pageFormat.getImageableHeight() / h;
				double scale = Math.min (1.0, Math.min (xScale, yScale));
				g2d.scale (scale, scale);

				// Shift graphics context's origin to bounding box's origin.
				g2d.translate
					(- bounds.getX() * scale,
					 - bounds.getY() * scale);

				// Draw the displayable object.
				disp.draw (g2d);

				return Printable.PAGE_EXISTS;
				}
			};

		// Print PostScript into file.
		FileOutputStream fos = new FileOutputStream (theFile);
		StreamPrintService sps = factories[0].getPrintService (fos);
		Doc doc = new SimpleDoc (thePrintData, theFlavor, null);
		try
			{
			sps.createPrintJob().print (doc, null);
			}
		catch (PrintException exc)
			{
			IOException exc2 = new IOException ("Cannot generate PostScript");
			exc2.initCause (exc);
			throw exc2;
			}
		}

	}
