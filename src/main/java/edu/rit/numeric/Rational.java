//******************************************************************************
//
// File:    Rational.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Rational
//
// This Java source file is copyright (C) 2002-2004 by Alan Kaminsky. All rights
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

package edu.rit.numeric;

import java.io.Serializable;

/**
 * Class Rational provides a 32-bit rational number. A 32-bit rational number is
 * the ratio of two 32-bit integers (type <TT>int</TT>). Operations are provided
 * for exact arithmetic and comparison with rational numbers.
 * <P>
 * Class Rational overrides the <TT>equals()</TT> and <TT>hashCode()</TT>
 * methods, making it suitable for use as a key in a hashed data structure like
 * a {@link java.util.HashMap </CODE>HashMap<CODE>} or {@link java.util.HashSet
 * </CODE>HashSet<CODE>}. However, a Rational object is mutable. Take care not
 * to change the value of a Rational object if it is used as a key in a hashed
 * data structure.
 * <P>
 * Class Rational is not multiple thread safe.
 *
 * @author  Alan Kaminsky
 * @version 01-Dec-2004
 */
public class Rational
	implements Comparable, Serializable
	{

// Hidden data members.

	private int numer;
	private int denom;

// Exported constructors.

	/**
	 * Construct a new rational number. Its value is 0.
	 */
	public Rational()
		{
		numer = 0;
		denom = 1;
		}

	/**
	 * Construct a new rational number. Its value is <TT>value</TT>.
	 *
	 * @param  value  Value.
	 */
	public Rational
		(int value)
		{
		numer = value;
		denom = 1;
		}

	/**
	 * Construct a new rational number. Its value is <TT>numer/denom</TT>.
	 *
	 * @param  numer  Numerator.
	 * @param  denom  Denominator.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>denom</TT> is 0.
	 */
	public Rational
		(int numer,
		 int denom)
		{
		if (denom == 0)
			{
			throw new ArithmeticException ("Divide by zero");
			}
		this.numer = numer;
		this.denom = denom;
		normalize();
		}

	/**
	 * Construct a new rational number. Its value is <TT>value</TT>.
	 *
	 * @param  value  Rational number.
	 */
	public Rational
		(Rational value)
		{
		this.numer = value.numer;
		this.denom = value.denom;
		}

	/**
	 * Construct a new rational number. Its value comes from the string
	 * <TT>s</TT>. The string must obey this syntax: optional minus sign
	 * (<TT>-</TT>), numerator (decimal integer), slash (<TT>/</TT>),
	 * denominator (decimal integer). The slash-and-denominator is optional. If
	 * present, the denominator must be greater than 0. There is no whitespace
	 * within the string.
	 *
	 * @param  s  String.
	 *
	 * @exception  NumberFormatException
	 *     (unchecked exception) Thrown if <TT>s</TT> cannot be parsed into a
	 *     rational number.
	 */
	public Rational
		(String s)
		{
		assign (s);
		}

// Exported operations.

	/**
	 * Set this rational number to the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x</TT>.
	 */
	public Rational assign
		(Rational x)
		{
		this.numer = x.numer;
		this.denom = x.denom;
		return this;
		}

	/**
	 * Set this rational number to the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x</TT>.
	 */
	public Rational assign
		(int x)
		{
		this.numer = x;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the value parsed from the given string. The
	 * string must obey this syntax: optional minus sign (<TT>-</TT>), numerator
	 * (decimal integer), slash (<TT>/</TT>), denominator (decimal integer). The
	 * slash-and-denominator is optional. If present, the denominator must be
	 * greater than 0. There is no whitespace within the string.
	 *
	 * @param  s  String.
	 *
	 * @return  This rational number, with its value set to <TT>s</TT>.
	 *
	 * @exception  NumberFormatException
	 *     (unchecked exception) Thrown if <TT>s</TT> cannot be parsed into a
	 *     rational number.
	 */
	public Rational assign
		(String s)
		{
		int numer, denom;
		int iSlash = s.indexOf ('/');
		if (iSlash == -1)
			{
			// No denominator.
			numer = Integer.parseInt (s);
			denom = 1;
			}
		else if (iSlash+1 < s.length())
			{
			// Slash and denominator.
			numer = Integer.parseInt (s.substring (0, iSlash));
			denom = Integer.parseInt (s.substring (iSlash+1));
			if (denom <= 0)
				{
				throw new NumberFormatException();
				}
			}
		else
			{
			// Slash but no denominator.
			throw new NumberFormatException();
			}
		this.numer = numer;
		this.denom = denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the negative of the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>-x</TT>.
	 */
	public Rational neg
		(Rational x)
		{
		this.numer = -x.numer;
		this.denom = x.denom;
		return this;
		}

	/**
	 * Set this rational number to the negative of the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>-x</TT>.
	 */
	public Rational neg
		(int x)
		{
		this.numer = -x;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the absolute value of the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>abs(x)</TT>.
	 */
	public Rational abs
		(Rational x)
		{
		this.numer = absval (x.numer);
		this.denom = x.denom;
		return this;
		}

	/**
	 * Set this rational number to the absolute value of the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>abs(x)</TT>.
	 */
	public Rational abs
		(int x)
		{
		this.numer = absval (x);
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the integer part of the given rational
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>int(x)</TT>.
	 */
	public Rational intPart
		(Rational x)
		{
		this.numer = x.numer / x.denom;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the fractional part of the given rational
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>frac(x)</TT>.
	 */
	public Rational fracPart
		(Rational x)
		{
		this.numer = x.numer % x.denom;
		this.denom = x.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the sum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x+y</TT>.
	 */
	public Rational add
		(Rational x,
		 Rational y)
		{
		int lcm = leastCommonMultiple (x.denom, y.denom);
		this.numer = (lcm / x.denom) * x.numer + (lcm / y.denom) * y.numer;
		this.denom = lcm;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the sum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x+y</TT>.
	 */
	public Rational add
		(Rational x,
		 int y)
		{
		this.numer = x.numer + y * x.denom;
		this.denom = x.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the sum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x+y</TT>.
	 */
	public Rational add
		(int x,
		 Rational y)
		{
		this.numer = x * y.denom + y.numer;
		this.denom = y.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the sum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x+y</TT>.
	 */
	public Rational add
		(int x,
		 int y)
		{
		this.numer = x + y;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the difference of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x-y</TT>.
	 */
	public Rational sub
		(Rational x,
		 Rational y)
		{
		int lcm = leastCommonMultiple (x.denom, y.denom);
		this.numer = (lcm / x.denom) * x.numer - (lcm / y.denom) * y.numer;
		this.denom = lcm;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the difference of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x-y</TT>.
	 */
	public Rational sub
		(Rational x,
		 int y)
		{
		this.numer = x.numer - y * x.denom;
		this.denom = x.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the difference of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x-y</TT>.
	 */
	public Rational sub
		(int x,
		 Rational y)
		{
		this.numer = x * y.denom - y.numer;
		this.denom = y.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the difference of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x-y</TT>.
	 */
	public Rational sub
		(int x,
		 int y)
		{
		this.numer = x - y;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the product of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x*y</TT>.
	 */
	public Rational mul
		(Rational x,
		 Rational y)
		{
		this.numer = x.numer * y.numer;
		this.denom = x.denom * y.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the product of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x*y</TT>.
	 */
	public Rational mul
		(Rational x,
		 int y)
		{
		this.numer = x.numer * y;
		this.denom = x.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the product of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x*y</TT>.
	 */
	public Rational mul
		(int x,
		 Rational y)
		{
		this.numer = x * y.numer;
		this.denom = y.denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the product of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x*y</TT>.
	 */
	public Rational mul
		(int x,
		 int y)
		{
		this.numer = x * y;
		this.denom = 1;
		return this;
		}

	/**
	 * Set this rational number to the quotient of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x/y</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>y</TT> is 0.
	 */
	public Rational div
		(Rational x,
		 Rational y)
		{
		if (y.numer == 0)
			{
			throw new ArithmeticException ("Divide by zero");
			}
		int numer = x.numer * y.denom;
		int denom = x.denom * y.numer;
		this.numer = numer;
		this.denom = denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the quotient of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x/y</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>y</TT> is 0.
	 */
	public Rational div
		(Rational x,
		 int y)
		{
		if (y == 0)
			{
			throw new ArithmeticException ("Divide by zero");
			}
		this.numer = x.numer;
		this.denom = x.denom * y;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the quotient of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x/y</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>y</TT> is 0.
	 */
	public Rational div
		(int x,
		 Rational y)
		{
		if (y.numer == 0)
			{
			throw new ArithmeticException ("Divide by zero");
			}
		int numer = x * y.denom;
		int denom = y.numer;
		this.numer = numer;
		this.denom = denom;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the quotient of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to <TT>x/y</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>y</TT> is 0.
	 */
	public Rational div
		(int x,
		 int y)
		{
		if (y == 0)
			{
			throw new ArithmeticException ("Divide by zero");
			}
		this.numer = x;
		this.denom = y;
		normalize();
		return this;
		}

	/**
	 * Set this rational number to the minimum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to the smaller of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational min
		(Rational x,
		 Rational y)
		{
		if (x.lt (y))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the minimum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to the smaller of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational min
		(Rational x,
		 int y)
		{
		if (x.lt (y))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the minimum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to the smaller of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational min
		(int x,
		 Rational y)
		{
		if (y.gt (x))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the minimum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to the smaller of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational min
		(int x,
		 int y)
		{
		this.assign (Math.min (x, y));
		return this;
		}

	/**
	 * Set this rational number to the maximum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to the larger of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational max
		(Rational x,
		 Rational y)
		{
		if (x.gt (y))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the maximum of the given numbers.
	 *
	 * @param  x  Rational number.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to the larger of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational max
		(Rational x,
		 int y)
		{
		if (x.gt (y))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the maximum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Rational number.
	 *
	 * @return  This rational number, with its value set to the larger of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational max
		(int x,
		 Rational y)
		{
		if (y.lt (x))
			{
			this.assign (x);
			}
		else
			{
			this.assign (y);
			}
		return this;
		}

	/**
	 * Set this rational number to the maximum of the given numbers.
	 *
	 * @param  x  Integer.
	 * @param  y  Integer.
	 *
	 * @return  This rational number, with its value set to the larger of
	 *          <TT>x</TT> and <TT>y</TT>.
	 */
	public Rational max
		(int x,
		 int y)
		{
		this.assign (Math.max (x, y));
		return this;
		}

	/**
	 * Determine if this rational number is equal to zero.
	 *
	 * @return  True if this rational number is equal to 0, false otherwise.
	 */
	public boolean eqZero()
		{
		return this.numer == 0;
		}

	/**
	 * Determine if this rational number is not equal to zero.
	 *
	 * @return  True if this rational number is not equal to 0, false otherwise.
	 */
	public boolean neZero()
		{
		return this.numer != 0;
		}

	/**
	 * Determine if this rational number is less than zero.
	 *
	 * @return  True if this rational number is less than 0, false otherwise.
	 */
	public boolean ltZero()
		{
		return this.numer < 0;
		}

	/**
	 * Determine if this rational number is less than or equal to zero.
	 *
	 * @return  True if this rational number is less than or equal to 0, false
	 *          otherwise.
	 */
	public boolean leZero()
		{
		return this.numer <= 0;
		}

	/**
	 * Determine if this rational number is greater than zero.
	 *
	 * @return  True if this rational number is greater than 0, false otherwise.
	 */
	public boolean gtZero()
		{
		return this.numer > 0;
		}

	/**
	 * Determine if this rational number is greater than or equal to zero.
	 *
	 * @return  True if this rational number is greater than or equal to 0,
	 *          false otherwise.
	 */
	public boolean geZero()
		{
		return this.numer >= 0;
		}

	/**
	 * Determine if this rational number is equal to the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is equal to <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean eq
		(Rational x)
		{
		return this.numer == x.numer && this.denom == x.denom;
		}

	/**
	 * Determine if this rational number is equal to the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is equal to <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean eq
		(int x)
		{
		return this.numer == x && this.denom == 1;
		}

	/**
	 * Determine if this rational number is not equal to the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is not equal to <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean ne
		(Rational x)
		{
		return this.numer != x.numer || this.denom != x.denom;
		}

	/**
	 * Determine if this rational number is not equal to the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is not equal to <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean ne
		(int x)
		{
		return this.numer != x || this.denom != 1;
		}

	/**
	 * Determine if this rational number is less than the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is less than <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean lt
		(Rational x)
		{
		return new Rational() .sub (this, x) .ltZero();
		}

	/**
	 * Determine if this rational number is less than the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is less than <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean lt
		(int x)
		{
		return new Rational() .sub (this, x) .ltZero();
		}

	/**
	 * Determine if this rational number is less than or equal to the given
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is less than or equal to
	 *          <TT>x</TT>, false otherwise.
	 */
	public boolean le
		(Rational x)
		{
		return new Rational() .sub (this, x) .leZero();
		}

	/**
	 * Determine if this rational number is less than or equal to the given
	 * number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is less than or equal to
	 *          <TT>x</TT>, false otherwise.
	 */
	public boolean le
		(int x)
		{
		return new Rational() .sub (this, x) .leZero();
		}

	/**
	 * Determine if this rational number is greater than the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is greater than <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean gt
		(Rational x)
		{
		return new Rational() .sub (this, x) .gtZero();
		}

	/**
	 * Determine if this rational number is greater than the given number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is greater than <TT>x</TT>, false
	 *          otherwise.
	 */
	public boolean gt
		(int x)
		{
		return new Rational() .sub (this, x) .gtZero();
		}

	/**
	 * Determine if this rational number is greater than or equal to the given
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  True if this rational number is greater than or equal to
	 *          <TT>x</TT>, false otherwise.
	 */
	public boolean ge
		(Rational x)
		{
		return new Rational() .sub (this, x) .geZero();
		}

	/**
	 * Determine if this rational number is greater than or equal to the given
	 * number.
	 *
	 * @param  x  Integer.
	 *
	 * @return  True if this rational number is greater than or equal to
	 *          <TT>x</TT>, false otherwise.
	 */
	public boolean ge
		(int x)
		{
		return new Rational() .sub (this, x) .geZero();
		}

	/**
	 * Converts this rational number to an integer. If this rational number is
	 * not an exact integer, any fractional part is truncated.
	 *
	 * @return  Integer value of this rational number.
	 */
	public int intValue()
		{
		return this.numer / this.denom;
		}

	/**
	 * Converts this rational number to a long integer. If this rational number
	 * is not an exact integer, any fractional part is truncated.
	 *
	 * @return  Long integer value of this rational number.
	 */
	public long longValue()
		{
		return ((long) this.numer) / ((long) this.denom);
		}

	/**
	 * Converts this rational number to a single precision floating point
	 * number.
	 *
	 * @return  Single precision floating point value of this rational number.
	 */
	public float floatValue()
		{
		return ((float) this.numer) / ((float) this.denom);
		}

	/**
	 * Converts this rational number to a double precision floating point
	 * number.
	 *
	 * @return  Double precision floating point value of this rational number.
	 */
	public double doubleValue()
		{
		return ((double) this.numer) / ((double) this.denom);
		}

	/**
	 * Returns a string version of this rational number. If this rational number
	 * is an exact integer, the string is in the form <TT>"&lt;numer&gt;"</TT>,
	 * otherwise the string is in the form
	 * <TT>"&lt;numer&gt;/&lt;denom&gt;"</TT>. If this rational number is less
	 * than 0, the string begins with a minus sign (<TT>-</TT>).
	 */
	public String toString()
		{
		if (this.denom == 1)
			{
			return Integer.toString (this.numer);
			}
		else
			{
			return this.numer + "/" + this.denom;
			}
		}

	/**
	 * Determine if this rational number is equal to the given object. To be
	 * equal, <TT>obj</TT> must be a non-null instance of class Rational whose
	 * value is the same as this rational number's value.
	 *
	 * @param  obj  Object to test.
	 *
	 * @return  True if this rational number is equal to <TT>obj</TT>, false
	 *          otherwise.
	 */
	public boolean equals
		(Object obj)
		{
		return
			obj instanceof Rational &&
			eq ((Rational) obj);
		}

	/**
	 * Returns a hash code for this rational number.
	 */
	public int hashCode()
		{
		return this.numer + this.denom;
		}

	/**
	 * Compare this rational number to the given object.
	 *
	 * @param  obj  Object to compare to.
	 *
	 * @return  A number less than, equal to, or greater than 0 if this rational
	 *          number is less than, equal to, or greater than <TT>obj</TT>,
	 *          respectively.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>obj</TT> is null.
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if <TT>obj</TT> is not an instance of
	 *     class Rational.
	 */
	public int compareTo
		(Object obj)
		{
		return new Rational() .sub (this, (Rational) obj) .numer;
		}

// Hidden operations.

	/**
	 * Normalize this rational number. Afterwards, the denominator is greater
	 * than 0, and the denominator does not divide the numerator.
	 */
	private void normalize()
		{
		int sign = signum (numer) * signum (denom);
		numer = absval (numer);
		denom = absval (denom);
		int gcd = greatestCommonDivisor (numer, denom);
		numer = sign * numer / gcd;
		denom = denom / gcd;
		}

	/**
	 * Returns the signum of the given number.
	 */
	private static int signum
		(int x)
		{
		return x < 0 ? -1 : 1;
		}

	/**
	 * Returns the absolute value of the given number.
	 */
	private static int absval
		(int x)
		{
		return x < 0 ? -x : x;
		}

	/**
	 * Returns the greatest common divisor of the given numbers. <TT>x</TT> and
	 * <TT>y</TT> are assumed to be greater than 0.
	 */
	private static int greatestCommonDivisor
		(int x,
		 int y)
		{
		int r;
		while (y != 0)
			{
			r = x % y;
			x = y;
			y = r;
			}
		return x;
		}

	/**
	 * Returns the least common multiple of the given numbers. <TT>x</TT> and
	 * <TT>y</TT> are assumed to be greater than 0.
	 */
	private static int leastCommonMultiple
		(int x,
		 int y)
		{
		return (x / greatestCommonDivisor (x, y)) * y;
		}

	}
