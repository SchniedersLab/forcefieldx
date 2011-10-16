//******************************************************************************
//
// File:    TitleDialog.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.TitleDialog
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
import edu.rit.swing.FontSelector;

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

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Class TitleDialog provides a modal dialog for specifying the title attributes
 * of a {@linkplain Plot}.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
class TitleDialog
	extends JDialog
	{

// Hidden data members.

	private static final int GAP = 10;

	private JTextField myTitle;
	private FontSelector myTitleFont;
	private DoubleTextField myTitleOffset;
	private boolean okButtonClicked;

// Exported constructors.

	/**
	 * Construct a new title dialog.
	 *
	 * @param  owner  Frame in which the dialog is displayed.
	 * @param  kind   Which kind of title is being formatted.
	 */
	public TitleDialog
		(Frame owner,
		 String kind)
		{
		super (owner, "Format " + kind + " Title", true);

		GridBagConstraints c;
		JLabel l;

		// Set up grid bag layout manager.
		Container pane = getContentPane();
		GridBagLayout layout = new GridBagLayout();
		pane.setLayout (layout);

		// Widgets for title.
		l = new JLabel ("Title");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets (GAP, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myTitle = new JTextField (20);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 0;
		c.insets = new Insets (GAP, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;
		layout.setConstraints (myTitle, c);
		pane.add (myTitle);

		// Widgets for title font.
		l = new JLabel ("Font");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 1;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myTitleFont = new FontSelector();
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 1;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.weightx = 1.0;
		layout.setConstraints (myTitleFont, c);
		pane.add (myTitleFont);

		// Widgets for title offset.
		l = new JLabel ("Offset");
		c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 2;
		c.insets = new Insets (0, GAP, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		layout.setConstraints (l, c);
		pane.add (l);
		myTitleOffset = new DoubleTextField (5);
		c = new GridBagConstraints();
		c.gridx = 1;
		c.gridy = 2;
		c.insets = new Insets (0, 0, GAP, GAP);
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		layout.setConstraints (myTitleOffset, c);
		pane.add (myTitleOffset);

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
		c.gridy = 3;
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
				myTitle.setSelectionStart (Integer.MAX_VALUE);
				myTitle.setSelectionEnd (Integer.MAX_VALUE);
				myTitle.requestFocusInWindow();
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
	 * Set the title text.
	 *
	 * @param  title  Title.
	 */
	public void setTitleText
		(String title)
		{
		myTitle.setText (title);
		}

	/**
	 * Set the title font.
	 *
	 * @param  titleFont  Title font.
	 */
	public void setTitleFont
		(Font titleFont)
		{
		myTitleFont.setSelectedFont (titleFont);
		}

	/**
	 * Set the title offset.
	 *
	 * @param  titleOffset  Title offset.
	 */
	public void setTitleOffset
		(double titleOffset)
		{
		myTitleOffset.value (titleOffset);
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
	 * Get the title text.
	 *
	 * @return  Title.
	 */
	public String getTitleText()
		{
		return myTitle.getText();
		}

	/**
	 * Get the title font.
	 *
	 * @return  Title font.
	 */
	public Font getTitleFont()
		{
		return myTitleFont.getSelectedFont();
		}

	/**
	 * Get the title offset.
	 *
	 * @return  Title offset.
	 */
	public double getTitleOffset()
		{
		return myTitleOffset.value();
		}

// Hidden operations.

	/**
	 * Processing when the "Okay" button is clicked.
	 */
	private void doOkay()
		{
		if (myTitleOffset.isOkay (0.0, Double.MAX_VALUE))
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

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		TitleDialog dialog = new TitleDialog (null);
//		dialog.setTitle ("This is a test");
//		dialog.setTitleOffset (12);
//		dialog.setVisible (true);
//		System.out.println (dialog.isOkay());
//		System.out.println (dialog.getTitle());
//		System.out.println (dialog.getTitleFont());
//		System.out.println (dialog.getTitleOffset());
//		dialog.dispose();
//		}

	}
