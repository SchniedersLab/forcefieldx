//******************************************************************************
//
// File:    ImagePanel.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.ImagePanel
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

package edu.rit.swing;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.RenderingHints;

import java.awt.geom.AffineTransform;

import java.awt.image.RenderedImage;

import javax.swing.JPanel;

/**
 * Class ImagePanel provides a GUI component that displays an image. If the
 * image is larger than the panel, the image is scaled to fit completely inside
 * the panel. The image is displayed centered in the panel.
 *
 * @author  Alan Kaminsky
 * @version 30-Jun-2006
 */
public class ImagePanel
	extends JPanel
	{

// Hidden data members.

	private RenderedImage myImage;

// Exported constructors.

	/**
	 * Construct a new image panel. The image panel initially displays no image.
	 */
	public ImagePanel()
		{
		setOpaque (true);
		}

// Exported operations.

	/**
	 * Specify the image to display in this image panel.
	 *
	 * @param  theImage  Image, or null to display no image.
	 */
	public void setImage
		(RenderedImage theImage)
		{
		myImage = theImage;
		repaint();
		}

	/**
	 * Get the image being displayed in this image panel.
	 *
	 * @return  Image, or null if no image is being displayed.
	 */
	public RenderedImage getImage()
		{
		return myImage;
		}

// Hidden operations.

	/**
	 * Paint this image panel in the given graphics context.
	 *
	 * @param  g  Graphics context.
	 */
	protected void paintComponent
		(Graphics g)
		{
		super.paintComponent (g);
		if (myImage == null) return;

		// Get a new 2-D graphics context.
		Graphics2D g2d = (Graphics2D) g.create();

		// Turn on antialiasing.
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 RenderingHints.VALUE_ANTIALIAS_ON);

		// Get drawing area.
		Dimension panelSize = getSize();
		Insets panelInsets = getInsets();
		double panelWidth =
			panelSize.getWidth() - panelInsets.left - panelInsets.right;
		double panelHeight =
			panelSize.getHeight() - panelInsets.top - panelInsets.bottom;

		// Get image area.
		double imageWidth = myImage.getWidth();
		double imageHeight = myImage.getHeight();

		// Compute scale factor along each dimension.
		double xscale = Math.min (panelWidth / imageWidth, 1.0);
		double yscale = Math.min (panelHeight / imageHeight, 1.0);

		// Set up transformation so image area is centered within drawing area.
		AffineTransform transform = new AffineTransform();
		if (xscale < yscale)
			{
			transform.translate
				(panelInsets.left,
				 panelInsets.top + (panelHeight - xscale * imageHeight) / 2.0);
			transform.scale (xscale, xscale);
			}
		else if (xscale > yscale)
			{
			transform.translate
				(panelInsets.left + (panelWidth - yscale * imageWidth) / 2.0,
				 panelInsets.top);
			transform.scale (yscale, yscale);
			}
		else
			{
			transform.translate
				(panelInsets.left + (panelWidth - xscale * imageWidth) / 2.0,
				 panelInsets.top + (panelHeight - yscale * imageHeight) / 2.0);
			transform.scale (xscale, yscale);
			}

		// Draw image.
		g2d.drawRenderedImage (myImage, transform);

		// Dispose of graphics context.
		g2d.dispose();
		}

	}
