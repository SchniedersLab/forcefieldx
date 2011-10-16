//******************************************************************************
//
// File:    ColorDisplayable.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.ColorDisplayable
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

package edu.rit.image;

import edu.rit.swing.Displayable;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Rectangle;

import java.awt.geom.Rectangle2D;

/**
 * Class ColorDisplayable provides an object with which to display a color image
 * in a Swing UI.
 *
 * @author  Alan Kaminsky
 * @version 10-Nov-2008
 */
class ColorDisplayable
	implements Displayable
	{

// Hidden data members.

	private int myHeight;
	private int myWidth;
	private int[][] myMatrix;

// Exported constructors.

	/**
	 * Construct a new color displayable object.
	 *
	 * @param  theHeight  Image height.
	 * @param  theWidth   Image width.
	 * @param  theMatrix  Pixel data matrix.
	 */
	public ColorDisplayable
		(int theHeight,
		 int theWidth,
		 int[][] theMatrix)
		{
		this.myHeight = theHeight;
		this.myWidth = theWidth;
		this.myMatrix = theMatrix;
		}

// Exported operations.

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 */
	public Rectangle2D getBoundingBox()
		{
		return new Rectangle2D.Double (0, 0, myWidth, myHeight);
		}

	/**
	 * Returns this displayable object's background paint.
	 */
	public Paint getBackgroundPaint()
		{
		return Color.LIGHT_GRAY;
		}

	/**
	 * Draw this drawable object in the given graphics context. Upon return
	 * from this method, the given graphics context's state (color, font,
	 * transform, clip, and so on) is the same as it was upon entry to this
	 * method.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		// Early return if pixel matrix is not allocated.
		if (myMatrix == null) return;

		// Save graphics context state.
		Paint oldPaint = g2d.getPaint();

		// Get bounds of clipping region. We only have to draw pixels within
		// that region.
		int rowlb;
		int rowub;
		int collb;
		int colub;
		Rectangle clipBounds = g2d.getClipBounds();
		if (clipBounds == null)
			{
			rowlb = 0;
			rowub = myHeight;
			collb = 0;
			colub = myWidth;
			}
		else
			{
			rowlb = Math.max (clipBounds.y, 0);
			rowub = Math.min (clipBounds.y + clipBounds.height, myHeight);
			collb = Math.max (clipBounds.x, 0);
			colub = Math.min (clipBounds.x + clipBounds.width, myWidth);
			}

		// Set up rectangle for drawing a run of pixels.
		Rectangle rect = new Rectangle();

		// Draw all rows.
		for (int r = rowlb; r < rowub; ++ r)
			{
			// Skip row if it is not allocated.
			int[] matrix_r = myMatrix[r];
			if (matrix_r == null) continue;

			// Determine which columns to scan. If none, skip row.
			int allocatedlb = Math.min (matrix_r.length, collb);
			int allocatedub = Math.min (matrix_r.length, colub);
			if (allocatedlb >= allocatedub) continue;

			// Scan all columns looking for runs of identical pixel values.
			int runColumn = allocatedlb;
			int runLength = 1;
			int runPixel = matrix_r[runColumn];
			for (int c = allocatedlb+1; c < allocatedub; ++ c)
				{
				int pixel = matrix_r[c];
				if (runPixel != pixel)
					{
					// Draw previous run.
					rect.setBounds (runColumn, r, runLength, 1);
					g2d.setColor (new Color (runPixel));
					g2d.fill (rect);

					// Begin next run.
					runColumn = c;
					runLength = 1;
					runPixel = pixel;
					}
				else
					{
					// Continue current run.
					++ runLength;
					}
				}

			// Draw final run on row.
			rect.setBounds (runColumn, r, runLength, 1);
			g2d.setColor (new Color (runPixel));
			g2d.fill (rect);
			}

		// Restore graphics context state.
		g2d.setPaint (oldPaint);
		}

	}
