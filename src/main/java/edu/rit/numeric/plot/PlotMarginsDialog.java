//******************************************************************************
//
// File:    PlotMarginsDialog.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.PlotMarginsDialog
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

package edu.rit.numeric.plot;

import edu.rit.swing.DoubleTextField;

import java.awt.Container;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * Class PlotMarginsDialog provides a modal dialog for specifying the plot
 * margin attributes of a {@linkplain Plot}.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
class PlotMarginsDialog
	extends JDialog
	{

// Hidden data members.

	private static final int GAP = 10;

	private DoubleTextField myTopMargin;
	private DoubleTextField myLeftMargin;
	private DoubleTextField myBottomMargin;
	private DoubleTextField myRightMargin;
	private boolean okButtonClicked;

// Exported constructors.

	/**
	 * Construct a new plot margins dialog.
	 *
	 * @param  owner  Frame in which the dialog is displayed.
	 */
	public PlotMarginsDialog
		(Frame owner)
		{
		super (owner, "Format Margins", true);

		GridBagConstraints c;
		JLabel l;

		// Set up grid bag layout manager.
		Container pane = getContentPane();
		GridBagLayout layout = new GridBagLayout();
		pane.setLayout (layout);

		// Widgets for top margin.
		l = new JLabel ("Top Margin");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets (GAP, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myTopMargin = new DoubleTextField (5);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 0;
		c.insets = new Insets (GAP, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myTopMargin, c);
		pane.add (myTopMargin);

		// Widgets for left margin.
		l = new JLabel ("Left Margin");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 1;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myLeftMargin = new DoubleTextField (5);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 1;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myLeftMargin, c);
		pane.add (myLeftMargin);

		// Widgets for bottom margin.
		l = new JLabel ("Bottom Margin");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 2;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myBottomMargin = new DoubleTextField (5);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 2;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myBottomMargin, c);
		pane.add (myBottomMargin);

		// Widgets for right margin.
		l = new JLabel ("Right Margin");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 3;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myRightMargin = new DoubleTextField (5);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 3;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myRightMargin, c);
		pane.add (myRightMargin);

		// "OK" and "Cancel" buttons.
		JPanel buttons = new JPanel();
		buttons.setLayout (new BoxLayout (buttons, BoxLayout.X_AXIS));
		JButton okButton = new JButton ("OK");
		buttons.add (okButton);
		okButton.addActionListener (new ActionListener()
			{
			public void actionPerformed (ActionEvent e)
				{
				doOkay();
				}
			});
		JButton cancelButton = new JButton ("Cancel");
		buttons.add (cancelButton);
		cancelButton.addActionListener (new ActionListener()
			{
			public void actionPerformed (ActionEvent e)
				{
				doCancel();
				}
			});
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 4;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 1.0;
		c.weighty = 1.0;
		layout.setConstraints (buttons, c);
		pane.add (buttons);
		getRootPane().setDefaultButton (okButton);

		// Set up window closing actions.
		setDefaultCloseOperation (JDialog.DO_NOTHING_ON_CLOSE);
		addWindowListener (new WindowAdapter()
			{
			public void windowActivated (WindowEvent e)
				{
				myTopMargin.setSelectionStart (Integer.MAX_VALUE);
				myTopMargin.setSelectionEnd (Integer.MAX_VALUE);
				myTopMargin.requestFocusInWindow();
				}
			public void windowClosing (WindowEvent e)
				{
				doCancel();
				}
			});

		pack();
		}

// Exported operations.

	/**
	 * Set the top margin.
	 *
	 * @param  margin  Margin.
	 */
	public void setTopMargin
		(double margin)
		{
		myTopMargin.value (margin);
		}

	/**
	 * Set the left margin.
	 *
	 * @param  margin  Margin.
	 */
	public void setLeftMargin
		(double margin)
		{
		myLeftMargin.value (margin);
		}

	/**
	 * Set the bottom margin.
	 *
	 * @param  margin  Margin.
	 */
	public void setBottomMargin
		(double margin)
		{
		myBottomMargin.value (margin);
		}

	/**
	 * Set the right margin.
	 *
	 * @param  margin  Margin.
	 */
	public void setRightMargin
		(double margin)
		{
		myRightMargin.value (margin);
		}

	/**
	 * Determine if the "OK" button was clicked.
	 *
	 * @return  True if the "OK" button was clicked, false otherwise.
	 */
	public boolean isOkay()
		{
		return okButtonClicked;
		}

	/**
	 * Get the top margin.
	 *
	 * @return  Margin.
	 */
	public double getTopMargin()
		{
		return myTopMargin.value();
		}

	/**
	 * Get the left margin.
	 *
	 * @return  Margin.
	 */
	public double getLeftMargin()
		{
		return myLeftMargin.value();
		}

	/**
	 * Get the bottom margin.
	 *
	 * @return  Margin.
	 */
	public double getBottomMargin()
		{
		return myBottomMargin.value();
		}

	/**
	 * Get the right margin.
	 *
	 * @return  Margin.
	 */
	public double getRightMargin()
		{
		return myRightMargin.value();
		}

// Hidden operations.

	/**
	 * Processing when the "Okay" button is clicked.
	 */
	private void doOkay()
		{
		if (myTopMargin.isOkay (0.0, Double.MAX_VALUE) &&
				myLeftMargin.isOkay (0.0, Double.MAX_VALUE) &&
				myBottomMargin.isOkay (0.0, Double.MAX_VALUE) &&
				myRightMargin.isOkay (0.0, Double.MAX_VALUE))
			{
			okButtonClicked = true;
			setVisible (false);
			}
		}

	/**
	 * Processing when the "Cancel" button is clicked.
	 */
	private void doCancel()
		{
		okButtonClicked = false;
		setVisible (false);
		}

	}
