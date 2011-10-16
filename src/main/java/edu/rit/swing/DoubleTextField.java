//******************************************************************************
//
// File:    DoubleTextField.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DoubleTextField
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
 * Class DoubleTextField provides a JTextField for entering a double precision
 * floating point value. The DoubleTextField can be incorporated into a window
 * or dialog.
 *
 * @author  Alan Kaminsky
 * @version 07-Mar-2011
 */
public class DoubleTextField
	extends JTextField
	{

// Exported constructors.

	/**
	 * Construct a new double text field with the given number of columns. The
	 * value is initially 0.
	 *
	 * @param  columns  Number of columns.
	 */
	public DoubleTextField
		(int columns)
		{
		super (columns);
		value (0.0);
		}

	/**
	 * Construct a new double text field with the given value and number of
	 * columns.
	 *
	 * @param  value    Initial value.
	 * @param  columns  Number of columns.
	 */
	public DoubleTextField
		(double value,
		 int columns)
		{
		super (columns);
		value (value);
		}

// Exported operations.

	/**
	 * Determine if this double text field contains a syntactically valid double
	 * precision floating point number in the given range. If it does, true is
	 * returned. If it does not, a beep is sounded, the input focus is set to
	 * this integer text field, all the text is selected, and false is returned.
	 *
	 * @param  lb  Lower bound.
	 * @param  ub  Upper bound.
	 *
	 * @return  True if this double text field contains a syntactically valid
	 *          double precision floating point number in the given range, false
	 *          otherwise.
	 */
	public boolean isOkay
		(double lb,
		 double ub)
		{
		// Check validity.
		boolean valid = false;
		try
			{
			double value = Double.parseDouble (getText());
			valid = (lb <= value && value <= ub);
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
	 * Get this double text field's value. If this double text field does not
	 * contain a syntactically valid double precision floating point number, 0
	 * is returned.
	 *
	 * @return  Value.
	 */
	public double value()
		{
		try
			{
			return Double.parseDouble (getText());
			}
		catch (NumberFormatException exc)
			{
			return 0.0;
			}
		}

	/**
	 * Set this double text field's value.
	 *
	 * @param  value  Value.
	 */
	public void value
		(double value)
		{
		setText (""+value);
		}

	}
