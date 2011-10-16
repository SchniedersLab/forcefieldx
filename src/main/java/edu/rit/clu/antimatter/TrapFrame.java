//******************************************************************************
//
// File:    TrapFrame.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.TrapFrame
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
import java.awt.Container;

import java.text.DecimalFormat;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * Class TrapFrame is a UI window that displays the positions of antiprotons and
 * the total momentum in the Antimatter Simulation.
 *
 * @author  Alan Kaminsky
 * @version 04-Feb-2008
 */
public class TrapFrame
	extends JFrame
	{

// Hidden constants.

	private static final int GAP = 3;

	private static final DecimalFormat FMT5 = new DecimalFormat ("0.0000E0");

// Hidden data members.

	// Particle position UI.
	private TrapPanel myTrapPanel;

	// Total momentum display.
	private JLabel myMomentum;

// Exported constructors.

	/**
	 * Construct a new particle window.
	 *
	 * @param  title     Window title.
	 * @param  pArray    Array containing antiprotons' positions.
	 * @param  trapSide  Length of each side of the trap.
	 */
	public TrapFrame
		(String title,
		 Vector2D[] pArray,
		 double trapSide)
		{
		super (title);

		Container pane = getContentPane();
		pane.setLayout (new BoxLayout (pane, BoxLayout.Y_AXIS));

		myTrapPanel = new TrapPanel (pArray, trapSide);
		pane.add (myTrapPanel);

		JPanel p = new JPanel();
		pane.add (p);
		p.setLayout (new BoxLayout (p, BoxLayout.X_AXIS));
		p.setBorder (BorderFactory.createEmptyBorder (GAP, GAP, GAP, GAP));
		p.setBackground (Color.white);

		myMomentum = new JLabel();
		setTotalMomentum (0.0);
		p.add (myMomentum);

		p.add (Box.createHorizontalGlue());

		pack();
		}

// Exported operations.

	/**
	 * Get this particle window's particle UI.
	 *
	 * @return  Particle UI.
	 */
	public TrapPanel getTrapPanel()
		{
		return myTrapPanel;
		}

	/**
	 * Set this particle window's total momentum.
	 *
	 * @param  mv  Total momentum.
	 */
	public void setTotalMomentum
		(double mv)
		{
		myMomentum.setText ("Total momentum = " + FMT5.format (mv));
		}

	}
