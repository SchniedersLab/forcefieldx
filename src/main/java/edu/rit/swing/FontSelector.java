//******************************************************************************
//
// File:    FontSelector.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.FontSelector
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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
import java.awt.Font;
import java.awt.GraphicsEnvironment;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

// For unit test main program:
//import java.awt.Container;
//import java.awt.Frame;
//import java.awt.event.ActionEvent;
//import java.awt.event.ActionListener;
//import javax.swing.JButton;
//import javax.swing.JDialog;

/**
 * Class FontSelector provides a JPanel for selecting a font. The JPanel can be
 * incorporated into a window or dialog. The JPanel includes a combo box for
 * selecting the font family name, a combo box for selecting the style
 * (normal/bold/italic/bold italic), a spinner for selecting the point size, and
 * some sample text in the selected font. Only integer point sizes from 1 to 144
 * are supported.
 *
 * @author  Alan Kaminsky
 * @version 16-Oct-2007
 */
public class FontSelector
	extends JPanel
	{

// Hidden data members.

	private static final int GAP = 5;

	private String[] myFontFamilies;
	private Font mySelectedFont;

	private JComboBox myFontFamilyComboBox;
	private JComboBox myFontStyleComboBox;
	private JSpinner myFontSizeSpinner;
	private JLabel myTextSample;

// Exported constructors.

	/**
	 * Construct a new font selector.
	 */
	public FontSelector()
		{
		Dimension d;

		// Widgets are laid out horizontally.
		setLayout (new BoxLayout (this, BoxLayout.X_AXIS));

		// Set up font family names and a 12-point normal font for each family.
		myFontFamilies =
			GraphicsEnvironment.getLocalGraphicsEnvironment()
				.getAvailableFontFamilyNames();

		// Set up combo box for font family name.
		myFontFamilyComboBox = new JComboBox (myFontFamilies);
		d = myFontFamilyComboBox.getPreferredSize();
		myFontFamilyComboBox.setMinimumSize (d);
		myFontFamilyComboBox.setMaximumSize (d);
		myFontFamilyComboBox.setPreferredSize (d);
		myFontFamilyComboBox.setEditable (false);
		myFontFamilyComboBox.addItemListener (new ItemListener()
			{
			public void itemStateChanged (ItemEvent e)
				{
				if (e.getStateChange() == ItemEvent.SELECTED)
					{
					updateSelectedFont();
					}
				}
			});
		add (myFontFamilyComboBox);
		add (Box.createHorizontalStrut (GAP));

		// Set up combo box for font style.
		myFontStyleComboBox =
			new JComboBox
				(new String[] {"Plain", "Bold", "Italic", "Bold Italic"});
		d = myFontStyleComboBox.getPreferredSize();
		myFontStyleComboBox.setMinimumSize (d);
		myFontStyleComboBox.setMaximumSize (d);
		myFontStyleComboBox.setPreferredSize (d);
		myFontStyleComboBox.setEditable (false);
		myFontStyleComboBox.addItemListener (new ItemListener()
			{
			public void itemStateChanged (ItemEvent e)
				{
				if (e.getStateChange() == ItemEvent.SELECTED)
					{
					updateSelectedFont();
					}
				}
			});
		add (myFontStyleComboBox);
		add (Box.createHorizontalStrut (GAP));

		// Set up spinner for font size.
		myFontSizeSpinner =
			new JSpinner (new SpinnerNumberModel (12, 1, 144, 1));
		d = myFontSizeSpinner.getPreferredSize();
		myFontSizeSpinner.setMinimumSize (d);
		myFontSizeSpinner.setMaximumSize (d);
		myFontSizeSpinner.setPreferredSize (d);
		myFontSizeSpinner.addChangeListener (new ChangeListener()
			{
			public void stateChanged (ChangeEvent e)
				{
				updateSelectedFont();
				}
			});
		add (myFontSizeSpinner);
		add (Box.createHorizontalStrut (GAP));

		// Set up text sample.
		myTextSample = new JLabel ("Quick Brown Fox 123");
		d = myTextSample.getPreferredSize();
		myTextSample.setMinimumSize (d);
		myTextSample.setMaximumSize (d);
		myTextSample.setPreferredSize (d);
		add (myTextSample);

		// Set default selected font.
		setSelectedFont (new Font ("SansSerif", Font.PLAIN, 12));
		}

// Exported operations.

	/**
	 * Get the selected font.
	 *
	 * @return  Font.
	 */
	public Font getSelectedFont()
		{
		return mySelectedFont;
		}

	/**
	 * Set the selected font to the given font.
	 *
	 * @param  font  Font.
	 */
	public void setSelectedFont
		(Font font)
		{
		myFontFamilyComboBox.setSelectedItem (font.getFamily());
		int style = Font.PLAIN;
		if (font.isBold()) style += Font.BOLD;
		if (font.isItalic()) style += Font.ITALIC;
		myFontStyleComboBox.setSelectedIndex (style);
		myFontSizeSpinner.setValue (new Integer (font.getSize()));
		}

// Hidden operations.

	/**
	 * Update the selected font based on UI actions.
	 */
	private void updateSelectedFont()
		{
		int family = myFontFamilyComboBox.getSelectedIndex();
		int style = myFontStyleComboBox.getSelectedIndex();
		int size = ((Integer) myFontSizeSpinner.getValue()).intValue();
		mySelectedFont = new Font (myFontFamilies[family], style, size);
		myTextSample.setFont (mySelectedFont);
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		final JDialog dialog = new JDialog ((Frame) null, "FontSelector", true);
//		Container pane = dialog.getContentPane();
//		pane.setLayout (new BoxLayout (pane, BoxLayout.Y_AXIS));
//		final FontSelector fontSelector = new FontSelector();
//		pane.add (fontSelector);
//		pane.add (Box.createVerticalStrut (GAP));
//		JButton ok = new JButton ("OK");
//		pane.add (ok);
//		ok.addActionListener (new ActionListener()
//			{
//			public void actionPerformed (ActionEvent e)
//				{
//				System.out.println (fontSelector.getSelectedFont());
//				dialog.setVisible (false);
//				}
//			});
//		dialog.pack();
//		dialog.setVisible (true);
//		System.exit (0);
//		}

	}
