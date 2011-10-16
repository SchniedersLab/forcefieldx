//******************************************************************************
//
// File:    LongTextField.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.LongTextField
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

import javax.swing.JTextField;

/**
 * Class LongTextField provides a JTextField for entering a long value. The
 * LongTextField can be incorporated into a window or dialog.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
public class LongTextField
	extends JTextField
	{

// Exported constructors.

	/**
	 * Construct a new long text field with the given number of columns. The
	 * value is initially 0.
	 *
	 * @param  columns  Number of columns.
	 */
	public LongTextField
		(int columns)
		{
		super (columns);
		value (0);
		}

	/**
	 * Construct a new long text field with the given value and number of
	 * columns.
	 *
	 * @param  value    Initial value.
	 * @param  columns  Number of columns.
	 */
	public LongTextField
		(long value,
		 int columns)
		{
		super (columns);
		value (value);
		}

// Exported operations.

	/**
	 * Determine if this long text field contains a syntactically valid long
	 * integer in the given range. If it does, true is returned. If it does not,
	 * a beep is sounded, the input focus is set to this long text field, all
	 * the text is selected, and false is returned.
	 *
	 * @param  lb  Lower bound.
	 * @param  ub  Upper bound.
	 *
	 * @return  True if this long text field contains a syntactically valid
	 *          long integer in the given range, false otherwise.
	 */
	public boolean isOkay
		(long lb,
		 long ub)
		{
		// Check validity.
		boolean valid = false;
		try
			{
			long value = Long.parseLong (getText());
			valid = lb <= value && value <= ub;
			}
		catch (NumberFormatException exc)
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
	 * Get this long text field's value. If this long text field does not
	 * contain a syntactically valid long integer, 0 is returned.
	 *
	 * @return  Value.
	 */
	public long value()
		{
		try
			{
			return Long.parseLong (getText());
			}
		catch (NumberFormatException exc)
			{
			return 0;
			}
		}

	/**
	 * Set this long text field's value.
	 *
	 * @param  value  Value.
	 */
	public void value
		(long value)
		{
		setText (""+value);
		}

	}
