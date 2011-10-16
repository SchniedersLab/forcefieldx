//******************************************************************************
//
// File:    DisplayablePanel.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DisplayablePanel
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
import java.awt.Paint;
import java.awt.Rectangle;
import java.awt.RenderingHints;

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import javax.swing.JPanel;

/**
 * Class DisplayablePanel is a Swing {@link javax.swing.JPanel
 * </CODE>JPanel<CODE>} that displays a {@link Displayable
 * </CODE>Displayable<CODE>} object. The displayable panel can then be made part
 * of a Swing GUI.
 * <P>
 * The displayable object is specified when the displayable panel is
 * constructed. A different displayable object can be specified later using the
 * <TT>display()</TT> method. The displayable panel's preferred size and the
 * displayable panel's bounds are set to the size of the displayable object when
 * the displayable panel is constructed and when the <TT>display()</TT> method
 * is called. If the same displayable panel will be displaying different
 * displayable objects, consider embedding the displayable panel in a
 * JScrollPane so the whole displayable object can be scrolled into view if the
 * displayable object becomes larger.
 * <P>
 * The displayable panel supports zooming (scaling) the displayable object. The
 * zoom factor is set when the displayable panel is constructed and when the
 * <TT>zoom()</TT> method is called. The displayable panel's preferred size and
 * bounds are actually those of the displayable object multiplied by the zoom
 * factor.
 * <P>
 * If the autofitting flag is true (see <TT>setAutofitting()</TT>), the
 * displayable panel will display a scaled version of the displayable object,
 * such that the displayable object is completely visible within the displayable
 * panel's current size.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public class DisplayablePanel
	extends JPanel
	{

// Hidden data members.

	private Displayable myDisplayable;
	private double myZoom = 1.0;
	private Dimension mySize;
	private boolean iamAutofitting;

// Exported constructors.

	/**
	 * Construct a new displayable panel. The displayable panel's preferred size
	 * and bounds are set to the size of <TT>theDisplayable</TT>.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public DisplayablePanel
		(Displayable theDisplayable)
		{
		super();
		display (theDisplayable);
		}

	/**
	 * Construct a new displayable panel with the given zoom factor. The
	 * displayable panel's preferred size and bounds are set to the size of
	 * <TT>theDisplayable</TT> multiplied by <TT>theZoom</TT>.
	 *
	 * @param  theDisplayable  Displayable object.
	 * @param  theZoom         Zoom factor.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theZoom</TT> &lt;= 0.
	 */
	public DisplayablePanel
		(Displayable theDisplayable,
		 double theZoom)
		{
		super();
		zoom (theZoom);
		display (theDisplayable);
		}

// Exported operations.

	/**
	 * Obtain this displayable panel's displayable object.
	 *
	 * @return  Displayable object this displayable panel is displaying.
	 */
	public Displayable display()
		{
		return myDisplayable;
		}
		
	/**
	 * Tell this displayable panel to display the given displayable object.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public void display
		(Displayable theDisplayable)
		{
		if (theDisplayable == null)
			{
			throw new NullPointerException();
			}
		myDisplayable = theDisplayable;
		if (! iamAutofitting)
			{
			Rectangle2D box = theDisplayable.getBoundingBox();
			double w = box.getWidth() * myZoom;
			double h = box.getHeight() * myZoom;
			Dimension size  = new Dimension();
			size.setSize (w, h);
			setPreferredSize (size);
			Rectangle bounds = getBounds();
			setBounds (bounds.x, bounds.y, (int) w, (int) h);
			}
		repaint();
		}

	/**
	 * Obtain this displayable panel's zoom factor.
	 *
	 * @return  Zoom factor.
	 */
	public double zoom()
		{
		return myZoom;
		}

	/**
	 * Specify this displayable panel's zoom factor.
	 *
	 * @param  theZoom  Zoom factor.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theZoom</TT> &lt;= 0.
	 */
	public void zoom
		(double theZoom)
		{
		if (theZoom <= 0.0) throw new IllegalArgumentException();
		myZoom = theZoom;
		if (myDisplayable != null) display (myDisplayable);
		}

	/**
	 * Determine whether this displayable panel is autofitting.
	 *
	 * @return  True if this displayable panel is autofitting, false otherwise.
	 */
	public boolean isAutofitting()
		{
		return iamAutofitting;
		}

	/**
	 * Set this displayable panel's autofitting flag to the given value. If the
	 * autofitting flag is true, the displayable panel will display a scaled
	 * version of the displayable object, such that the displayable object is
	 * completely visible within the displayable panel's current size.
	 *
	 * @param  autofit  True to autofit, false otherwise.
	 */
	public void setAutofitting
		(boolean autofit)
		{
		iamAutofitting = autofit;
		repaint();
		}

	/**
	 * Paint this displayable panel into the given graphics context.
	 *
	 * @param  g  Graphics context.
	 */
	protected void paintComponent
		(Graphics g)
		{
		super.paintComponent (g);

		Graphics2D g2d = (Graphics2D) g;

		// Save graphics context.
		Paint oldPaint = g2d.getPaint();
		RenderingHints oldHints = g2d.getRenderingHints();
		AffineTransform oldTransform = g2d.getTransform();

		// Turn on antialiasing.
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setRenderingHint
			(RenderingHints.KEY_TEXT_ANTIALIASING,
			 RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		// Paint one way if autofitting is turned on.
		if (iamAutofitting)
			{
			// Paint background.
			Rectangle panelsize = getBounds (null);
			g2d.setPaint (myDisplayable.getBackgroundPaint());
			g2d.fill (panelsize);

			// Determine X and Y scale factors.
			double pw = panelsize.getWidth();
			double ph = panelsize.getHeight();
			Rectangle2D objectsize = myDisplayable.getBoundingBox();
			double ow = objectsize.getWidth();
			double oh = objectsize.getHeight();
			double scalex = pw / ow;
			double scaley = ph / oh;

			// Set transform if X scale factor is smaller.
			if (scalex < scaley)
				{
				g2d.translate (0.0, (ph - oh * scalex) * 0.5);
				g2d.scale (scalex, scalex);
				}

			// Set transform if Y scale factor is smaller.
			else
				{
				g2d.translate ((pw - ow * scaley) * 0.5, 0.0);
				g2d.scale (scaley, scaley);
				}

			// Draw displayable object.
			myDisplayable.draw (g2d);
			}

		// Paint another way if autofitting is turned off.
		else
			{
			// Set zoom factor.
			g2d.scale (myZoom, myZoom);

			// Paint background.
			g2d.setPaint (myDisplayable.getBackgroundPaint());
			g2d.fill (myDisplayable.getBoundingBox());

			// Draw displayable object.
			myDisplayable.draw (g2d);
			}

		// Restore graphics context.
		g2d.setPaint (oldPaint);
		g2d.setRenderingHints (oldHints);
		g2d.setTransform (oldTransform);
		}

	}
