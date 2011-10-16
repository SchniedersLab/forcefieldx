//******************************************************************************
//
// File:    PlotAreaDialog.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.PlotAreaDialog
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

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * Class PlotAreaDialog provides a modal dialog for specifying the plot area
 * attributes of a {@linkplain Plot}.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
class PlotAreaDialog
	extends JDialog
	{

// Hidden data members.

	private static final int GAP = 10;

	private DoubleTextField myPlotAreaWidth;
	private DoubleTextField myPlotAreaHeight;
	private JRadioButton myMajorGridLinesOn;
	private JRadioButton myMajorGridLinesOff;
	private JRadioButton myMinorGridLinesOn;
	private JRadioButton myMinorGridLinesOff;
	private boolean okButtonClicked;

// Exported constructors.

	/**
	 * Construct a new plot area dialog.
	 *
	 * @param  owner  Frame in which the dialog is displayed.
	 */
	public PlotAreaDialog
		(Frame owner)
		{
		super (owner, "Format Plot Area", true);

		GridBagConstraints c;
		JLabel l;

		// Set up grid bag layout manager.
		Container pane = getContentPane();
		GridBagLayout layout = new GridBagLayout();
		pane.setLayout (layout);

		// Widgets for plot area width.
		l = new JLabel ("Plot Area Width");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets (GAP, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myPlotAreaWidth = new DoubleTextField (10);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 0;
		c.insets = new Insets (GAP, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myPlotAreaWidth, c);
		pane.add (myPlotAreaWidth);

		// Widgets for plot area height.
		l = new JLabel ("Plot Area Height");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 1;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myPlotAreaHeight = new DoubleTextField (10);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 1;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myPlotAreaHeight, c);
		pane.add (myPlotAreaHeight);

		// Widgets for major grid lines.
		l = new JLabel ("Major Grid Lines");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 2;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel majorPanel = new JPanel();
		majorPanel.setLayout (new BoxLayout (majorPanel, BoxLayout.X_AXIS));
		ButtonGroup majorGroup = new ButtonGroup();
		myMajorGridLinesOn = new JRadioButton ("On");
		majorPanel.add (myMajorGridLinesOn);
		majorGroup.add (myMajorGridLinesOn);
		majorPanel.add (Box.createHorizontalStrut (GAP));
		myMajorGridLinesOff = new JRadioButton ("Off");
		majorPanel.add (myMajorGridLinesOff);
		majorGroup.add (myMajorGridLinesOff);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 2;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (majorPanel, c);
		pane.add (majorPanel);

		// Widgets for minor grid lines.
		l = new JLabel ("Minor Grid Lines");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 3;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel minorPanel = new JPanel();
		minorPanel.setLayout (new BoxLayout (minorPanel, BoxLayout.X_AXIS));
		ButtonGroup minorGroup = new ButtonGroup();
		myMinorGridLinesOn = new JRadioButton ("On");
		minorPanel.add (myMinorGridLinesOn);
		minorGroup.add (myMinorGridLinesOn);
		minorPanel.add (Box.createHorizontalStrut (GAP));
		myMinorGridLinesOff = new JRadioButton ("Off");
		minorPanel.add (myMinorGridLinesOff);
		minorGroup.add (myMinorGridLinesOff);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 3;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (minorPanel, c);
		pane.add (minorPanel);

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
				myPlotAreaWidth.setSelectionStart (Integer.MAX_VALUE);
				myPlotAreaWidth.setSelectionEnd (Integer.MAX_VALUE);
				myPlotAreaWidth.requestFocusInWindow();
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
	 * Set the plot area width.
	 *
	 * @param  width  Width.
	 */
	public void setPlotAreaWidth
		(double width)
		{
		myPlotAreaWidth.value (width);
		}

	/**
	 * Set the plot area height.
	 *
	 * @param  height  Height.
	 */
	public void setPlotAreaHeight
		(double height)
		{
		myPlotAreaHeight.value (height);
		}

	/**
	 * Set the major grid lines on or off.
	 *
	 * @param  isOn  True if major grid lines are on, false if they are off.
	 */
	public void setMajorGridLines
		(boolean isOn)
		{
		if (isOn)
			{
			myMajorGridLinesOn.setSelected (true);
			}
		else
			{
			myMajorGridLinesOff.setSelected (true);
			}
		}

	/**
	 * Set the minor grid lines on or off.
	 *
	 * @param  isOn  True if minor grid lines are on, false if they are off.
	 */
	public void setMinorGridLines
		(boolean isOn)
		{
		if (isOn)
			{
			myMinorGridLinesOn.setSelected (true);
			}
		else
			{
			myMinorGridLinesOff.setSelected (true);
			}
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
	 * Get the plot area width.
	 *
	 * @return  Width.
	 */
	public double getPlotAreaWidth()
		{
		return myPlotAreaWidth.value();
		}

	/**
	 * Get the plot area height.
	 *
	 * @return  Height.
	 */
	public double getPlotAreaHeight()
		{
		return myPlotAreaHeight.value();
		}

	/**
	 * Determine whether the major grid lines are on or off.
	 *
	 * @return  True if major grid lines are on, false if they are off.
	 */
	public boolean getMajorGridLines()
		{
		return myMajorGridLinesOn.isSelected();
		}

	/**
	 * Determine whether the minor grid lines are on or off.
	 *
	 * @return  True if minor grid lines are on, false if they are off.
	 */
	public boolean getMinorGridLines()
		{
		return myMinorGridLinesOn.isSelected();
		}

// Hidden operations.

	/**
	 * Processing when the "Okay" button is clicked.
	 */
	private void doOkay()
		{
		if (myPlotAreaWidth.isOkay (0.0, Double.MAX_VALUE) &&
				myPlotAreaHeight.isOkay (0.0, Double.MAX_VALUE))
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
