//******************************************************************************
//
// File:    Quadratic.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Quadratic
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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

/**
 * Class Quadratic solves for the real roots of a quadratic equation with real
 * coefficients. The quadratic equation is of the form
 * <P>
 * <I>ax</I><SUP>2</SUP> + <I>bx</I> + <I>c</I> = 0
 * <P>
 * To solve a quadratic equation, construct an instance of class Quadratic; call
 * the Quadratic object's <TT>solve()</TT> method, passing in the coefficients
 * <I>a</I>, <I>b</I>, and <I>c</I>; and obtain the roots from the Quadratic
 * object's fields. The number of (real) roots, either 0 or 2, is stored in
 * field <TT>nRoots</TT>. If there are no roots, fields <TT>x1</TT> and
 * <TT>x2</TT> are set to NaN. If there are two roots, they are stored in fields
 * <TT>x1</TT> and <TT>x2</TT> in descending order.
 * <P>
 * The same Quadratic object may be used to solve several quadratic equations.
 * Each time the <TT>solve()</TT> method is called, the solution is stored in
 * the Quadratic object's fields.
 * <P>
 * The formulas for the roots of a quadratic equation come from:
 * <P>
 * E. Weisstein. "Quadratic equation." From <I>MathWorld</I>--A Wolfram Web
 * Resource.
 * <A HREF="http://mathworld.wolfram.com/QuadraticEquation.html" TARGET="_top">http://mathworld.wolfram.com/QuadraticEquation.html</A>
 *
 * @author  Alan Kaminsky
 * @version 05-Mar-2008
 */
public class Quadratic
	{

// Hidden constants.

	private static final double TWO_PI = 2.0 * Math.PI;
	private static final double FOUR_PI = 4.0 * Math.PI;

// Exported fields.

	/**
	 * The number of real roots.
	 */
	public int nRoots;

	/**
	 * The first real root.
	 */
	public double x1;

	/**
	 * The second real root.
	 */
	public double x2;

// Exported constructors.

	/**
	 * Construct a new Quadratic object.
	 */
	public Quadratic()
		{
		}

// Exported operations.

	/**
	 * Solve the quadratic equation with the given coefficients. The results are
	 * stored in this Quadratic object's fields.
	 *
	 * @param  a  Coefficient of <I>x</I><SUP>2</SUP>.
	 * @param  b  Coefficient of <I>x</I>.
	 * @param  c  Constant coefficient.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if <TT>a</TT> is 0; in other words, the
	 *     coefficients do not represent a quadratic equation.
	 */
	public void solve
		(double a,
		 double b,
		 double c)
		{
		// Verify preconditions.
		if (a == 0.0)
			{
			throw new DomainException ("Quadratic.solve(): a = 0");
			}

		// Compute discriminant.
		double d = b*b - 4.0*a*c;

		if (d >= 0.0)
			{
			// Two real roots.
			nRoots = 2;
			double q = -0.5 * (b + sgn(b)*Math.sqrt(d));
			x1 = q / a;
			x2 = c / q;
			sortRoots();
			}

		else
			{
			// No real roots.
			nRoots = 0;
			x1 = Double.NaN;
			x2 = Double.NaN;
			}
		}

// Hidden operations.

	/**
	 * Returns the signum of x.
	 */
	private static double sgn
		(double x)
		{
		return x < 0.0 ? -1.0 : 1.0;
		}

	/**
	 * Sort the roots into descending order.
	 */
	private void sortRoots()
		{
		if (x1 < x2)
			{
			double tmp = x1; x1 = x2; x2 = tmp;
			}
		}

// Unit test main program.

	/**
	 * Unit test main program.
	 * <P>
	 * Usage: java edu.rit.numeric.Quadratic <I>a</I> <I>b</I> <I>c</I>
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 3) usage();
		double a = Double.parseDouble (args[0]);
		double b = Double.parseDouble (args[1]);
		double c = Double.parseDouble (args[2]);
		Quadratic quadratic = new Quadratic();
		quadratic.solve (a, b, c);
		if (quadratic.nRoots == 2)
			{
			System.out.println ("x1 = " + quadratic.x1);
			System.out.println ("x2 = " + quadratic.x2);
			}
		else
			{
			System.out.println ("No real roots");
			}
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.numeric.Quadratic <a> <b> <c>");
		System.err.println ("Solves ax^2 + bx + c = 0");
		System.exit (1);
		}

	}
