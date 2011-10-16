//******************************************************************************
//
// File:    AxisDialog.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.AxisDialog
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

import edu.rit.swing.DecimalFormatTextField;
import edu.rit.swing.DoubleTextField;
import edu.rit.swing.FontSelector;
import edu.rit.swing.IntegerTextField;

import java.awt.Container;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.text.DecimalFormat;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * Class AxisDialog provides a modal dialog for specifying the plot axis
 * attributes of a {@linkplain Plot}.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
class AxisDialog
	extends JDialog
	{

// Hidden data members.

	private static final int GAP = 10;

	private JRadioButton myLinearButton;
	private JRadioButton myLogarithmicButton;
	private DoubleTextField myAxisStart;
	private JCheckBox myAxisStartAuto;
	private DoubleTextField myAxisEnd;
	private JCheckBox myAxisEndAuto;
	private IntegerTextField myAxisMajorDivisions;
	private JCheckBox myAxisMajorDivisionsAuto;
	private IntegerTextField myAxisMinorDivisions;
	private JCheckBox myAxisMinorDivisionsAuto;
	private DoubleTextField myAxisCrossing;
	private JCheckBox myAxisCrossingAuto;
	private FontSelector myAxisTickFont;
	private DecimalFormatTextField myAxisTickFormat;
	private DoubleTextField myAxisTickScale;
	private boolean okButtonClicked;

// Exported constructors.

	/**
	 * Construct a new axis dialog.
	 *
	 * @param  owner  Frame in which the dialog is displayed.
	 * @param  kind   Kind of axis (<TT>"X"</TT> or <TT>"Y"</TT>).
	 */
	public AxisDialog
		(Frame owner,
		 String kind)
		{
		super (owner, "Format " + kind + " Axis", true);

		GridBagConstraints c;
		JLabel l;

		// Set up grid bag layout manager.
		Container pane = getContentPane();
		GridBagLayout layout = new GridBagLayout();
		pane.setLayout (layout);

		// Widgets for axis kind.
		l = new JLabel ("Kind");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets (GAP, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel kindPanel = new JPanel();
		kindPanel.setLayout (new BoxLayout (kindPanel, BoxLayout.X_AXIS));
		ButtonGroup kindGroup = new ButtonGroup();
		myLinearButton = new JRadioButton ("Linear");
		kindPanel.add (myLinearButton);
		kindGroup.add (myLinearButton);
		kindPanel.add (Box.createHorizontalStrut (GAP));
		myLogarithmicButton = new JRadioButton ("Logarithmic");
		kindPanel.add (myLogarithmicButton);
		kindGroup.add (myLogarithmicButton);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 0;
		c.insets = new Insets (GAP, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (kindPanel, c);
		pane.add (kindPanel);

		// Widgets for axis start.
		l = new JLabel ("Start Value");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 1;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel startPanel = new JPanel();
		startPanel.setLayout (new BoxLayout (startPanel, BoxLayout.X_AXIS));
		myAxisStart = new DoubleTextField (10);
		startPanel.add (myAxisStart);
		startPanel.add (Box.createHorizontalStrut (GAP));
		myAxisStartAuto = new JCheckBox ("Automatic");
		startPanel.add (myAxisStartAuto);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 1;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (startPanel, c);
		pane.add (startPanel);
		myAxisStartAuto.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				myAxisStart.setEnabled (! myAxisStartAuto.isSelected());
				}
			});

		// Widgets for axis end.
		l = new JLabel ("End Value");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 2;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel endPanel = new JPanel();
		endPanel.setLayout (new BoxLayout (endPanel, BoxLayout.X_AXIS));
		myAxisEnd = new DoubleTextField (10);
		endPanel.add (myAxisEnd);
		endPanel.add (Box.createHorizontalStrut (GAP));
		myAxisEndAuto = new JCheckBox ("Automatic");
		endPanel.add (myAxisEndAuto);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 2;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (endPanel, c);
		pane.add (endPanel);
		myAxisEndAuto.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				myAxisEnd.setEnabled (! myAxisEndAuto.isSelected());
				}
			});

		// Widgets for major divisions.
		l = new JLabel ("Major Divisions");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 3;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel majorPanel = new JPanel();
		majorPanel.setLayout (new BoxLayout (majorPanel, BoxLayout.X_AXIS));
		myAxisMajorDivisions = new IntegerTextField (10);
		majorPanel.add (myAxisMajorDivisions);
		majorPanel.add (Box.createHorizontalStrut (GAP));
		myAxisMajorDivisionsAuto = new JCheckBox ("Automatic");
		majorPanel.add (myAxisMajorDivisionsAuto);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 3;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (majorPanel, c);
		pane.add (majorPanel);
		myAxisMajorDivisionsAuto.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				myAxisMajorDivisions.setEnabled
					(! myAxisMajorDivisionsAuto.isSelected());
				}
			});

		// Widgets for minor divisions.
		l = new JLabel ("Minor Divisions");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 4;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel minorPanel = new JPanel();
		minorPanel.setLayout (new BoxLayout (minorPanel, BoxLayout.X_AXIS));
		myAxisMinorDivisions = new IntegerTextField (10);
		minorPanel.add (myAxisMinorDivisions);
		minorPanel.add (Box.createHorizontalStrut (GAP));
		myAxisMinorDivisionsAuto = new JCheckBox ("Automatic");
		minorPanel.add (myAxisMinorDivisionsAuto);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 4;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (minorPanel, c);
		pane.add (minorPanel);
		myAxisMinorDivisionsAuto.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				myAxisMinorDivisions.setEnabled
					(! myAxisMinorDivisionsAuto.isSelected());
				}
			});

		// Widgets for crossing.
		l = new JLabel ("Crossing Value");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 5;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		JPanel crossPanel = new JPanel();
		crossPanel.setLayout (new BoxLayout (crossPanel, BoxLayout.X_AXIS));
		myAxisCrossing = new DoubleTextField (10);
		crossPanel.add (myAxisCrossing);
		crossPanel.add (Box.createHorizontalStrut (GAP));
		myAxisCrossingAuto = new JCheckBox ("Automatic");
		crossPanel.add (myAxisCrossingAuto);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 5;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (crossPanel, c);
		pane.add (crossPanel);
		myAxisCrossingAuto.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				myAxisCrossing.setEnabled (! myAxisCrossingAuto.isSelected());
				}
			});

		// Widgets for tick font.
		l = new JLabel ("Tick Font");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 6;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myAxisTickFont = new FontSelector();
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 6;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myAxisTickFont, c);
		pane.add (myAxisTickFont);

		// Widgets for tick format.
		l = new JLabel ("Tick Format");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 7;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myAxisTickFormat = new DecimalFormatTextField (10);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 7;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myAxisTickFormat, c);
		pane.add (myAxisTickFormat);

		// Widgets for tick scale.
		l = new JLabel ("Scale Factor");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 8;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myAxisTickScale = new DoubleTextField (10);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 8;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myAxisTickScale, c);
		pane.add (myAxisTickScale);

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
		c.gridy = 9;
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
				myLinearButton.requestFocusInWindow();
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
	 * Set the axis kind.
	 *
	 * @param  kind  Plot.LINEAR or Plot.LOGARITHMIC.
	 */
	public void setAxisKind
		(Plot.AxisKind kind)
		{
		switch (kind)
			{
			case LINEAR:
				myLinearButton.setSelected (true);
				break;
			case LOGARITHMIC:
				myLogarithmicButton.setSelected (true);
				break;
			}
		}

	/**
	 * Set the axis start value.
	 *
	 * @param  value  Start value, or Double.NaN if automatic.
	 */
	public void setAxisStart
		(double value)
		{
		if (Double.isNaN (value))
			{
			myAxisStart.setText ("");
			myAxisStart.setEnabled (false);
			myAxisStartAuto.setSelected (true);
			}
		else
			{
			myAxisStart.value (value);
			myAxisStart.setEnabled (true);
			myAxisStartAuto.setSelected (false);
			}
		}

	/**
	 * Set the axis end value.
	 *
	 * @param  value  End value, or Double.NaN if automatic.
	 */
	public void setAxisEnd
		(double value)
		{
		if (Double.isNaN (value))
			{
			myAxisEnd.setText ("");
			myAxisEnd.setEnabled (false);
			myAxisEndAuto.setSelected (true);
			}
		else
			{
			myAxisEnd.value (value);
			myAxisEnd.setEnabled (true);
			myAxisEndAuto.setSelected (false);
			}
		}

	/**
	 * Set the axis major divisions.
	 *
	 * @param  value  Major divisions. 10 = automatic.
	 */
	public void setAxisMajorDivisions
		(int value)
		{
		if (value == 10)
			{
			myAxisMajorDivisions.setText ("");
			myAxisMajorDivisions.setEnabled (false);
			myAxisMajorDivisionsAuto.setSelected (true);
			}
		else
			{
			myAxisMajorDivisions.value (value);
			myAxisMajorDivisions.setEnabled (true);
			myAxisMajorDivisionsAuto.setSelected (false);
			}
		}

	/**
	 * Set the axis minor divisions.
	 *
	 * @param  value  Minor divisions. 1 = automatic.
	 */
	public void setAxisMinorDivisions
		(int value)
		{
		if (value == 1)
			{
			myAxisMinorDivisions.setText ("");
			myAxisMinorDivisions.setEnabled (false);
			myAxisMinorDivisionsAuto.setSelected (true);
			}
		else
			{
			myAxisMinorDivisions.value (value);
			myAxisMinorDivisions.setEnabled (true);
			myAxisMinorDivisionsAuto.setSelected (false);
			}
		}

	/**
	 * Set the axis crossing value.
	 *
	 * @param  value  Crossing value, or Double.NaN if automatic.
	 */
	public void setAxisCrossing
		(double value)
		{
		if (Double.isNaN (value))
			{
			myAxisCrossing.setText ("");
			myAxisCrossing.setEnabled (false);
			myAxisCrossingAuto.setSelected (true);
			}
		else
			{
			myAxisCrossing.value (value);
			myAxisCrossing.setEnabled (true);
			myAxisCrossingAuto.setSelected (false);
			}
		}

	/**
	 * Set the tick font.
	 *
	 * @param  font  Font.
	 */
	public void setTickFont
		(Font font)
		{
		myAxisTickFont.setSelectedFont (font);
		}

	/**
	 * Set the tick format.
	 *
	 * @param  pattern  Pattern.
	 */
	public void setTickFormat
		(DecimalFormat pattern)
		{
		myAxisTickFormat.pattern (pattern);
		}

	/**
	 * Set the tick scale factor.
	 *
	 * @param  scale  Scale factor.
	 */
	public void setTickScale
		(double scale)
		{
		myAxisTickScale.value (scale);
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
	 * Get the axis kind.
	 *
	 * @return  Plot.LINEAR or Plot.LOGARITHMIC.
	 */
	public Plot.AxisKind getAxisKind()
		{
		return myLinearButton.isSelected() ? Plot.LINEAR : Plot.LOGARITHMIC;
		}

	/**
	 * Get the axis start value.
	 *
	 * @return  Start value, or Double.NaN if automatic.
	 */
	public double getAxisStart()
		{
		return
			myAxisStartAuto.isSelected() ?
				Double.NaN :
				myAxisStart.value();
		}

	/**
	 * Get the axis end value.
	 *
	 * @return  End value, or Double.NaN if automatic.
	 */
	public double getAxisEnd()
		{
		return
			myAxisEndAuto.isSelected() ?
				Double.NaN :
				myAxisEnd.value();
		}

	/**
	 * Get the axis major divisions.
	 *
	 * @return  Major divisions. 10 = automatic.
	 */
	public int getAxisMajorDivisions()
		{
		return
			myAxisMajorDivisionsAuto.isSelected() ?
				10 :
				myAxisMajorDivisions.value();
		}

	/**
	 * Get the axis minor divisions.
	 *
	 * @return  Minor divisions. 1 = automatic.
	 */
	public int getAxisMinorDivisions()
		{
		return
			myAxisMinorDivisionsAuto.isSelected() ?
				1 :
				myAxisMinorDivisions.value();
		}

	/**
	 * Get the axis crossing value.
	 *
	 * @return  Crossing value, or Double.NaN if automatic.
	 */
	public double getAxisCrossing()
		{
		return
			myAxisCrossingAuto.isSelected() ?
				Double.NaN :
				myAxisCrossing.value();
		}

	/**
	 * Get the tick font.
	 *
	 * @return  Font.
	 */
	public Font getTickFont()
		{
		return myAxisTickFont.getSelectedFont();
		}

	/**
	 * Get the tick format.
	 *
	 * @return  Pattern.
	 */
	public DecimalFormat getTickFormat()
		{
		return myAxisTickFormat.pattern();
		}

	/**
	 * Get the tick scale factor.
	 *
	 * @return  Scale factor.
	 */
	public double getTickScale()
		{
		return myAxisTickScale.value();
		}

// Hidden operations.

	/**
	 * Processing when the "Okay" button is clicked.
	 */
	private void doOkay()
		{
		if
			((myAxisStartAuto.isSelected() ||
				myAxisStart.isOkay (-Double.MAX_VALUE, Double.MAX_VALUE)) &&
			 (myAxisEndAuto.isSelected() ||
				myAxisEnd.isOkay (-Double.MAX_VALUE, Double.MAX_VALUE)) &&
			 (myAxisMajorDivisionsAuto.isSelected() ||
				myAxisMajorDivisions.isOkay (1, Integer.MAX_VALUE)) &&
			 (myAxisMinorDivisionsAuto.isSelected() ||
				myAxisMinorDivisions.isOkay (1, Integer.MAX_VALUE)) &&
			 (myAxisCrossingAuto.isSelected() ||
				myAxisCrossing.isOkay (-Double.MAX_VALUE, Double.MAX_VALUE)) &&
			 myAxisTickFormat.isOkay() &&
			 myAxisTickScale.isOkay (-Double.MAX_VALUE, Double.MAX_VALUE))
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
