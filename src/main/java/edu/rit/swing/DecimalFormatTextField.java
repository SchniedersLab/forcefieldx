//******************************************************************************
//
// File:    DecimalFormatTextField.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DecimalFormatTextField
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

package edu.rit.swing;

import java.awt.Toolkit;

import java.text.DecimalFormat;

import javax.swing.JTextField;

/**
 * Class DecimalFormatTextField provides a JTextField for entering a pattern
 * string for a DecimalFormat object. The DecimalFormatTextField can be
 * incorporated into a window or dialog.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
public class DecimalFormatTextField
	extends JTextField
	{

// Exported constructors.

	/**
	 * Construct a new decimal format text field with the given number of
	 * columns. The pattern is initially <TT>"0"</TT>.
	 *
	 * @param  columns  Number of columns.
	 */
	public DecimalFormatTextField
		(int columns)
		{
		super (columns);
		pattern (new DecimalFormat ("0"));
		}

	/**
	 * Construct a new integer text field with the given pattern and number of
	 * columns.
	 *
	 * @param  pattern  Initial pattern.
	 * @param  columns  Number of columns.
	 */
	public DecimalFormatTextField
		(DecimalFormat pattern,
		 int columns)
		{
		super (columns);
		pattern (pattern);
		}

// Exported operations.

	/**
	 * Determine if this decimal format text field contains a syntactically
	 * valid pattern in the given range. If it does, true is returned. If it
	 * does not, a beep is sounded, the input focus is set to this decimal
	 * format text field, all the text is selected, and false is returned.
	 *
	 * @return  True if this decimal format text field contains a syntactically
	 *          valid pattern, false otherwise.
	 */
	public boolean isOkay()
		{
		// Check validity.
		boolean valid = false;
		try
			{
			DecimalFormat pattern = new DecimalFormat (getText());
			valid = true;
			}
		catch (Throwable exc)
			{
			}

		// Take action if invalid.
		if (! valid)
			{
			requestFocusInWindow();
			setSelectionStart (0);
			setSelectionEnd (Integer.MAX_VALUE);
			Toolkit.getDefaultToolkit().beep();
			}

		return valid;
		}

	/**
	 * Get this decimal format text field's pattern. If this decimal format text
	 * field does not contain a syntactically valid pattern, null is returned.
	 *
	 * @return  Pattern.
	 */
	public DecimalFormat pattern()
		{
		DecimalFormat pattern = null;
		try
			{
			pattern = new DecimalFormat (getText());
			}
		catch (Throwable exc)
			{
			}
		return pattern;
		}

	/**
	 * Set this decimal format text field's pattern.
	 *
	 * @param  pattern  Pattern.
	 */
	public void pattern
		(DecimalFormat pattern)
		{
		setText (pattern.toPattern());
		}

	}
