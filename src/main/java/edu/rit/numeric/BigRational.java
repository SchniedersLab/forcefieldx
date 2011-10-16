//******************************************************************************
//
// File:    BigRational.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.BigRational
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * Class BigRational provides an arbitrary precision rational number. An
 * arbitrary precision rational number is the ratio of two arbitrary precision
 * integers (type java.math.BigInteger). Operations are provided for exact
 * arithmetic with rational numbers.
 * <P>
 * A rational number is said to be <B>normalized</B> if
 * GCD(numerator,denominator) = 1. The methods below do <B><I>not</I></B>
 * automatically normalize the rational number. Thus, the numerator and
 * denominator tend to get larger and larger as operations are performed on the
 * rational number. To reduce the numerator and denominator to lowest terms
 * again, call the <TT>normalize()</TT> method. It is up to you to decide
 * whether to normalize the rational number after each operation, or after a
 * series of operations.
 * <P>
 * Class BigRational provides the <TT>equals()</TT> and <TT>hashCode()</TT>
 * methods, and BigRational objects can be used as keys in hashed data
 * structures. However, BigRational objects are mutable. If a BigRational object
 * is used as a hash key, be sure not to change its value.
 * <P>
 * Class BigRational is not multiple thread safe.
 *
 * @author  Alan Kaminsky
 * @version 31-Dec-2009
 */
public class BigRational
	extends Number
	implements Comparable<BigRational>
	{

// Hidden data members.

	private BigInteger numer;
	private BigInteger denom;

// Exported constructors.

	/**
	 * Construct a new rational number. Its value is 0.
	 */
	public BigRational()
		{
		numer = BigInteger.ZERO;
		denom = BigInteger.ONE;
		}

	/**
	 * Construct a new rational number. Its value is <TT>x</TT>.
	 *
	 * @param  x  Long integer.
	 */
	public BigRational
		(long x)
		{
		assign (x);
		}

	/**
	 * Construct a new rational number. Its value is <TT>n/d</TT>.
	 *
	 * @param  n  Numerator.
	 * @param  d  Denominator.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>d</TT> is 0.
	 */
	public BigRational
		(long n,
		 long d)
		{
		assign (n, d);
		}

	/**
	 * Construct a new rational number. Its value is <TT>x</TT>.
	 *
	 * @param  x  Big integer.
	 */
	public BigRational
		(BigInteger x)
		{
		assign (x);
		}

	/**
	 * Construct a new rational number. Its value is <TT>n/d</TT>.
	 *
	 * @param  n  Numerator.
	 * @param  d  Denominator.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>d</TT> is 0.
	 */
	public BigRational
		(BigInteger n,
		 BigInteger d)
		{
		assign (n, d);
		}

	/**
	 * Construct a new rational number. Its value is <TT>x</TT>.
	 *
	 * @param  x  Rational number.
	 */
	public BigRational
		(BigRational x)
		{
		assign (x);
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
	public BigRational
		(String s)
		{
		assign (s);
		}

// Exported operations.

	/**
	 * Returns this rational number's numerator.
	 *
	 * @return  Numerator.
	 */
	public BigInteger numerator()
		{
		return numer;
		}

	/**
	 * Returns this rational number's denominator.
	 *
	 * @return  Denominator.
	 */
	public BigInteger denominator()
		{
		return denom;
		}

	/**
	 * Set this rational number to the given number.
	 *
	 * @param  x  Long integer.
	 *
	 * @return  This rational number, with its value set to <TT>x</TT>.
	 */
	public BigRational assign
		(long x)
		{
		this.numer = BigInteger.valueOf (x);
		this.denom = BigInteger.ONE;
		return this;
		}

	/**
	 * Set this rational number to the given fraction.
	 *
	 * @param  n  Numerator.
	 * @param  d  Denominator.
	 *
	 * @return  This rational number, with its value set to <TT>n/d</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>d</TT> is 0.
	 */
	public BigRational assign
		(long n,
		 long d)
		{
		if (d == 0)
			{
			throw new ArithmeticException
				("BigRational.assign(): Zero denominator");
			}
		this.numer = BigInteger.valueOf (n);
		this.denom = BigInteger.valueOf (d);
		return this;
		}

	/**
	 * Set this rational number to the given number.
	 *
	 * @param  x  Big integer.
	 *
	 * @return  This rational number, with its value set to <TT>x</TT>.
	 */
	public BigRational assign
		(BigInteger x)
		{
		this.numer = x;
		this.denom = BigInteger.ONE;
		return this;
		}

	/**
	 * Set this rational number to the given fraction.
	 *
	 * @param  n  Numerator.
	 * @param  d  Denominator.
	 *
	 * @return  This rational number, with its value set to <TT>n/d</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>d</TT> is 0.
	 */
	public BigRational assign
		(BigInteger n,
		 BigInteger d)
		{
		if (d.equals (BigInteger.ZERO))
			{
			throw new ArithmeticException
				("BigRational.assign(): Zero denominator");
			}
		this.numer = n;
		this.denom = d;
		return this;
		}

	/**
	 * Set this rational number to the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>x</TT>.
	 */
	public BigRational assign
		(BigRational x)
		{
		this.numer = x.numer;
		this.denom = x.denom;
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
	public BigRational assign
		(String s)
		{
		BigInteger numer, denom;
		int iSlash = s.indexOf ('/');
		if (iSlash == -1)
			{
			// No denominator.
			numer = new BigInteger (s);
			denom = BigInteger.ONE;
			}
		else if (iSlash+1 < s.length())
			{
			// Slash and denominator.
			numer = new BigInteger (s.substring (0, iSlash));
			denom = new BigInteger (s.substring (iSlash+1));
			if (denom.compareTo (BigInteger.ZERO) <= 0)
				{
				throw new NumberFormatException
					("BigRational.assign(): Negative denominator not allowed");
				}
			}
		else
			{
			// Slash but no denominator.
			throw new NumberFormatException
				("BigRational.assign(): Missing denominator after /");
			}
		this.numer = numer;
		this.denom = denom;
		return this;
		}

	/**
	 * Set this rational number to the integer part of itself.
	 *
	 * @return  This rational number, with its value set to <TT>int(this)</TT>.
	 */
	public BigRational intPart()
		{
		numer = numer.divide (denom);
		denom = BigInteger.ONE;
		return this;
		}

	/**
	 * Set this rational number to the fractional part of itself.
	 *
	 * @return  This rational number, with its value set to <TT>frac(this)</TT>.
	 */
	public BigRational fracPart()
		{
		numer = numer.remainder (denom);
		return this;
		}

	/**
	 * Set this rational number to the sum of itself and the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this+x</TT>.
	 */
	public BigRational add
		(BigRational x)
		{
		this.numer =
			this.numer.multiply(x.denom).add
				(x.numer.multiply (this.denom));
		this.denom = this.denom.multiply (x.denom);
		return this;
		}

	/**
	 * Set this rational number to the difference of itself and the given
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this-x</TT>.
	 */
	public BigRational sub
		(BigRational x)
		{
		this.numer =
			this.numer.multiply(x.denom).subtract
				(x.numer.multiply (this.denom));
		this.denom = this.denom.multiply (x.denom);
		return this;
		}

	/**
	 * Set this rational number to the product of itself and the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this*x</TT>.
	 */
	public BigRational mul
		(BigRational x)
		{
		this.numer = this.numer.multiply (x.numer);
		this.denom = this.denom.multiply (x.denom);
		return this;
		}

	/**
	 * Set this rational number to the quotient of itself and the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this/x</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>x</TT> is 0.
	 */
	public BigRational div
		(BigRational x)
		{
		if (x.numer.equals (BigInteger.ZERO))
			{
			throw new ArithmeticException
				("BigRational.div(): Divide by zero");
			}
		this.numer = this.numer.multiply (x.denom);
		this.denom = this.denom.multiply (x.numer);
		return this;
		}

	/**
	 * Set this rational number to the remainder when divided by the given
	 * number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to
	 *          <TT>frac(this/x)*x</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if <TT>x</TT> is 0.
	 */
	public BigRational rem
		(BigRational x)
		{
		return div (x) .fracPart() .mul (x);
		}

	/**
	 * Increment this rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this+1</TT>.
	 */
	public BigRational incr()
		{
		this.numer = this.numer.add (this.denom);
		return this;
		}

	/**
	 * Decrement this rational number.
	 *
	 * @return  This rational number, with its value set to <TT>this-1</TT>.
	 */
	public BigRational decr()
		{
		this.numer = this.numer.subtract (this.denom);
		return this;
		}

	/**
	 * Set this rational number to the absolute value of itself.
	 *
	 * @return  This rational number, with its value set to <TT>abs(this)</TT>.
	 */
	public BigRational abs()
		{
		this.numer = this.numer.abs();
		this.denom = this.denom.abs();
		return this;
		}

	/**
	 * Set this rational number to the negative of itself.
	 *
	 * @return  This rational number, with its value set to <TT>-this</TT>.
	 */
	public BigRational negate()
		{
		this.numer = this.numer.negate();
		return this;
		}

	/**
	 * Set this rational number to the reciprocal of itself.
	 *
	 * @return  This rational number, with its value set to <TT>1/this</TT>.
	 *
	 * @exception  ArithmeticException
	 *     (unchecked exception) Thrown if this rational number is 0.
	 */
	public BigRational recip()
		{
		if (this.numer.equals (BigInteger.ZERO))
			{
			throw new ArithmeticException
				("BigRational.recip(): Divide by zero");
			}
		BigInteger tmp = this.numer;
		this.numer = this.denom;
		this.denom = tmp;
		return this;
		}

	/**
	 * Set this rational number to the smaller of itself and the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to
	 *          <TT>min(this,x)</TT>.
	 */
	public BigRational min
		(BigRational x)
		{
		if (this.compareTo (x) > 0) this.assign (x);
		return this;
		}

	/**
	 * Set this rational number to the larger of itself and the given number.
	 *
	 * @param  x  Rational number.
	 *
	 * @return  This rational number, with its value set to
	 *          <TT>max(this,x)</TT>.
	 */
	public BigRational max
		(BigRational x)
		{
		if (this.compareTo (x) < 0) this.assign (x);
		return this;
		}

	/**
	 * Normalize this rational number. Afterwards, the denominator is greater
	 * than or equal to 1, and the denominator does not divide the numerator; in
	 * other words, GCD(numerator,denominator) = 1.
	 *
	 * @return  This rational number, normalized.
	 */
	public BigRational normalize()
		{
		int sign = numer.signum() * denom.signum();
		numer = numer.abs();
		denom = denom.abs();
		BigInteger gcd = numer.gcd (denom);
		numer = numer.divide (gcd);
		denom = denom.divide (gcd);
		if (sign < 0) numer = numer.negate();
		return this;
		}

	/**
	 * Converts this rational number to an integer.
	 *
	 * @return  Integer value of this rational number.
	 */
	public int intValue()
		{
		BigDecimal n = new BigDecimal (numer);
		BigDecimal d = new BigDecimal (denom);
		return n.divide (d, 1, RoundingMode.HALF_UP) .intValue();
		}

	/**
	 * Converts this rational number to a long integer.
	 *
	 * @return  Long integer value of this rational number.
	 */
	public long longValue()
		{
		BigDecimal n = new BigDecimal (numer);
		BigDecimal d = new BigDecimal (denom);
		return n.divide (d, 1, RoundingMode.HALF_UP) .longValue();
		}

	/**
	 * Converts this rational number to a single precision floating point
	 * number.
	 *
	 * @return  Single precision floating point value of this rational number.
	 */
	public float floatValue()
		{
		BigDecimal n = new BigDecimal (numer);
		BigDecimal d = new BigDecimal (denom);
		return n.divide (d, MathContext.DECIMAL32) .floatValue();
		}

	/**
	 * Converts this rational number to a double precision floating point
	 * number.
	 *
	 * @return  Double precision floating point value of this rational number.
	 */
	public double doubleValue()
		{
		BigDecimal n = new BigDecimal (numer);
		BigDecimal d = new BigDecimal (denom);
		return n.divide (d, MathContext.DECIMAL64) .doubleValue();
		}

	/**
	 * Compare this rational number to the given rational number.
	 *
	 * @param  x  Rational number to compare.
	 *
	 * @return  A number less than, equal to, or greater than 0 if this rational
	 *          number is less than, equal to, or greater than <TT>x</TT>,
	 *          respectively.
	 */
	public int compareTo
		(BigRational x)
		{
		BigInteger diff =
			this.numer.multiply(x.denom).subtract
				(x.numer.multiply (this.denom));
		return diff.compareTo (BigInteger.ZERO);
		}

	/**
	 * Determine if this rational number equals the given object.
	 * <P>
	 * <I>Note:</I> Class BigRational provides the <TT>equals()</TT> and
	 * <TT>hashCode()</TT> methods, and BigRational objects can be used as keys
	 * in hashed data structures. However, BigRational objects are mutable. If a
	 * BigRational object is used as a hash key, be sure not to change its
	 * value.
	 * <P>
	 * <I>Note:</I> Two rational numbers are equal only if their numerators and
	 * denominators are equal. Thus, it is possible for <TT>equals()</TT> to
	 * return false even if the two rational numbers have the same numerical
	 * value. To ensure that <TT>equals()</TT> will return true if the two
	 * rational numbers have the same numerical value, normalize the two
	 * rational numbers first.
	 *
	 * @param  obj  Object to compare.
	 *
	 * @return  True if this rational number equals <TT>obj</TT>, false
	 *          otherwise.
	 */
	public boolean equals
		(Object obj)
		{
		return
			(obj instanceof BigRational) &&
			(((BigRational) obj).numer.equals (this.numer)) &&
			(((BigRational) obj).denom.equals (this.denom));
		}

	/**
	 * Returns a hash code for this rational number.
	 * <P>
	 * <I>Note:</I> Class BigRational provides the <TT>equals()</TT> and
	 * <TT>hashCode()</TT> methods, and BigRational objects can be used as keys
	 * in hashed data structures. However, BigRational objects are mutable. If a
	 * BigRational object is used as a hash key, be sure not to change its
	 * value.
	 * <P>
	 * <I>Note:</I> Two rational numbers have the same hash code only if their
	 * numerators and denominators are equal. Thus, it is possible for
	 * <TT>hashCode()</TT> to return different values even if the two rational
	 * numbers have the same numerical value. To ensure that <TT>hashCode()</TT>
	 * will return the same value if the two rational numbers have the same
	 * numerical value, normalize the two rational numbers first.
	 *
	 * @return  Hash code.
	 */
	public int hashCode()
		{
		return numer.hashCode()*31 + denom.hashCode();
		}

	/**
	 * Returns a string version of this rational number. If this rational
	 * number's denominator is 1, the string is in the form
	 * <TT>"&lt;numer&gt;"</TT>, otherwise the string is in the form
	 * <TT>"&lt;numer&gt;/&lt;denom&gt;"</TT>. The negative sign, if any, is
	 * attached to the numerator.
	 * <P>
	 * <I>Note:</I> The <TT>toString()</TT> method yields the numerator and
	 * denominator as they are, without normalizing this rational number.
	 */
	public String toString()
		{
		String ns = this.numer.toString();
		String ds = this.denom.toString();
		if (ns.charAt(0) == '-' && ds.charAt(0) == '-')
			{
			ns = ns.substring (1);
			ds = ds.substring (1);
			}
		else if (ns.charAt(0) != '-' && ds.charAt(0) == '-')
			{
			ns = "-" + ns;
			ds = ds.substring (1);
			}
		if (ds.equals ("1"))
			{
			return ns;
			}
		else
			{
			return ns + "/" + ds;
			}
		}

	}
