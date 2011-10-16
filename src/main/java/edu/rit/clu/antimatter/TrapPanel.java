//******************************************************************************
//
// File:    TrapPanel.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.TrapPanel
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

package edu.rit.clu.antimatter;

import edu.rit.vector.Vector2D;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import java.awt.geom.Ellipse2D;

import javax.swing.JPanel;

/**
 * Class TrapPanel is a UI that displays the positions of antiprotons in the
 * Antimatter Simulation. It displays the particle positions stored in an array
 * of {@linkplain edu.rit.vector.Vector2D} objects.
 *
 * @author  Alan Kaminsky
 * @version 04-Feb-2008
 */
public class TrapPanel
	extends JPanel
	{

// Hidden constants.

	private static final int W = 500;
	private static final double DIAM = 5;
	private static final double RADIUS = DIAM/2;
	private static final Color PARTICLE_COLOR = Color.red;

// Hidden data members.

	// Array containing particle positions.
	private Vector2D[] pArray;

	// Scale factor for display.
	private double scale;

	// For drawing antiprotons.
	private Ellipse2D dot = new Ellipse2D.Double();

	// State of the display.
	private boolean iamWaitingForRepaint = false;

// Exported constructors.

	/**
	 * Construct a new particle panel.
	 *
	 * @param  pArray    Array containing antiprotons' positions.
	 * @param  trapSide  Length of each side of the trap.
	 */
	public TrapPanel
		(Vector2D[] pArray,
		 double trapSide)
		{
		super();
		Dimension theSize = new Dimension (W, W);
		setMinimumSize (theSize);
		setPreferredSize (theSize);
		setMaximumSize (theSize);

		this.pArray = pArray;

		scale = W / trapSide;
		}

// Exported operations.

	/**
	 * Update the display for the next time step.
	 */
	public synchronized void step()
		{
		// Repaint the display.
		iamWaitingForRepaint = true;
		repaint (0L, 0, 0, W, W);

		// Wait until repainting is complete.
		while (iamWaitingForRepaint)
			{
			try
				{
				wait();
				}
			catch (InterruptedException exc)
				{
				}
			}
		}

// Hidden operations.

	/**
	 * Paint this trap panel's display.
	 *
	 * @param  g  Graphics context.
	 */
	protected synchronized void paintComponent
		(Graphics g)
		{
		// Perform superclass behavior.
		super.paintComponent (g);

		// Set up graphics context.
		Graphics2D g2d = (Graphics2D) g.create();
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor (PARTICLE_COLOR);

		// Draw all antiprotons.
		int n = pArray.length;
		for (int i = 0; i < n; ++ i)
			{
			Vector2D p = pArray[i];
			dot.setFrame (p.x*scale-RADIUS, p.y*scale-RADIUS, DIAM, DIAM);
			g2d.fill (dot);
			}

		// Notify that repainting is complete.
		iamWaitingForRepaint = false;
		notifyAll();
		}

	}
