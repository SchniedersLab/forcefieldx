//******************************************************************************
//
// File:    CurveSmoothing.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.CurveSmoothing
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

package edu.rit.numeric;

import java.util.Arrays;

/**
 * Class CurveSmoothing provides operations useful for creating smooth cubic
 * spline curves.
 * <P>
 * We are given a sequence of <I>n</I> coordinate values, <TT>u[0]</TT> ..
 * <TT>u[n-1]</TT>. These may be either the X coordinates or the Y coordinates
 * of a sequence of 2-D points. We wish to join successive pairs of points with
 * a sequence of cubic Bezier curves to create an overall cubic spline curve. We
 * need to know the X or Y coordinates of the two Bezier control points for each
 * Bezier curve. The Bezier control point coordinates are designated by
 * <TT>a[i]</TT> and <TT>c[i]</TT>, such that Bezier curve <I>i</I> is defined
 * by <TT>u[i]</TT>, <TT>a[i]</TT>, <TT>c[i]</TT>, and <TT>u[i+1]</TT> in that
 * order.
 * <P>
 * There are several cases:
 * <OL TYPE=1>
 * <LI>
 * The curve is closed (there is a Bezier curve from point <I>n</I>-1 back to
 * point 0). We need the <I>n</I> Bezier control points <TT>a[0]</TT> ..
 * <TT>a[n-1]</TT> and the <I>n</I> Bezier control points <TT>c[0]</TT> ..
 * <TT>c[n-1]</TT>. The first Bezier curve will be
 * <TT>u[0]--a[0]--c[0]--u[1]</TT>, and the last Bezier curve will be
 * <TT>u[n-1]--a[n-1]--c[n-1]--u[0]</TT>.
 * <BR>&nbsp;
 * <LI>
 * The curve is open (there is no Bezier curve from point <I>n</I>-1 back to
 * point 0). We need the <I>n</I>-1 Bezier control points <TT>a[0]</TT> ..
 * <TT>a[n-2]</TT> and the <I>n</I>-1 Bezier control points <TT>c[0]</TT> ..
 * <TT>c[n-2]</TT>. The first Bezier curve will be
 * <TT>u[0]--a[0]--c[0]--u[1]</TT>, and the last Bezier curve will be
 * <TT>u[n-2]--a[n-2]--c[n-2]--u[n-1]</TT>.
 * <P>
 * For the initial point <TT>u[0]</TT>, an additional condition must be
 * specified, either:
 * <OL TYPE=a>
 * <LI>
 * Zero curvature at the initial point; or
 * <LI>
 * Initial direction, specified as a straight line from a coordinate
 * <TT>uInitial</TT> to <TT>u[0]</TT>.
 * </OL>
 * <P>
 * For the final point <TT>u[n-1]</TT>, an additional condition must be
 * specified, either:
 * <OL TYPE=a>
 * <LI>
 * Zero curvature at the final point; or
 * <LI>
 * Final direction, specified as a straight line from <TT>u[n-1]</TT> to a
 * coordinate <TT>uFinal</TT>.
 * </OL>
 * </OL>
 *
 * @author Alan Kaminsky
 * @version 07-Jul-2007
 */
public class CurveSmoothing
	{

// Prevent construction.

	private CurveSmoothing()
		{
		}

// Exported operations.

	/**
	 * Compute the Bezier control point coordinates for a closed smooth curve.
	 *
	 * @param  u  An input array of coordinates. Elements at indexes <TT>i</TT>
	 *            .. <TT>i+n-1</TT> are used.
	 * @param  a  An output array of coordinates for the first Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-1</TT> are
	 *            used.
	 * @param  c  An output array of coordinates for the second Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-1</TT> are
	 *            used.
	 * @param  i  Index of first element used in input and output arrays.
	 * @param  n  Number of elements used in input and output arrays. Must be
	 *            greater than or equal to 3.
	 *
 	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &lt; 3.
 	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the control points cannot be
	 *     calculated.
	 */
	public static void computeBezierClosed
		(double[] u,
		 double[] a,
		 double[] c,
		 int i,
		 int n)
		{
		int j;

		// First, compute vector of second derivatives of u, upp, by solving
		// this cyclic tridiagonal linear system (illustrated for n=6):
		//     [4 1       1]   [upp0]   [6 u5 - 12 u0 + 6 u1]
		//     [1 4 1      ]   [upp1]   [6 u0 - 12 u1 + 6 u2]
		//     [  1 4 1    ] x [upp2] = [6 u1 - 12 u2 + 6 u3]
		//     [    1 4 1  ]   [upp3]   [6 u2 - 12 u3 + 6 u4]
		//     [      1 4 1]   [upp4]   [6 u3 - 12 u4 + 6 u5]
		//     [1       1 4]   [upp5]   [6 u4 - 12 u5 + 6 u0]

		// Allocate storage for tridiagonal matrix, solution vector, and right-
		// hand side vector.
		double[] d = new double [n];
		double[] e = new double [n];
		double[] upp = new double [n];
		double[] rhs = new double [n];

		// Fill in tridiagonal matrix and right-hand side vector.
		Arrays.fill (d, 4.0);
		Arrays.fill (e, 1.0);
		for (j = 0; j < n; ++ j)
			{
			rhs[j] =
				   6.0 * u[i + (j+n-1)%n]
				- 12.0 * u[i + j]
				+  6.0 * u[i + (j+1)%n];
			}

		// Solve the linear system.
		Tridiagonal.solveSymmetricCyclic (d, e, rhs, upp);

		// Lastly, compute Bezier control point vectors a and c from u and upp.
		for (j = 0; j < n; ++ j)
			{
			a[i+j] = 2.0*u[i+j]/3.0 + u[i+(j+1)%n]/3.0 - upp[j]/9.0 -
						upp[(j+1)%n]/18.0;
			c[i+j] = u[i+j]/3.0 + 2.0*u[i+(j+1)%n]/3.0 - upp[j]/18.0 -
						upp[(j+1)%n]/9.0;
			}

		}

	/**
	 * Compute the Bezier control point coordinates for an open smooth curve
	 * with zero initial curvature and zero final curvature.
	 *
	 * @param  u  An input array of coordinates. Elements at indexes <TT>i</TT>
	 *            .. <TT>i+n-1</TT> are used.
	 * @param  a  An output array of coordinates for the first Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  c  An output array of coordinates for the second Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  i  Index of first element used in input and output arrays.
	 * @param  n  Number of elements used in input array; one fewer elements
	 *            used in output arrays. Must be greater than or equal to 3.
	 *
 	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &lt; 3.
 	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the control points cannot be
	 *     calculated.
	 */
	public static void computeBezierOpen
		(double[] u,
		 double[] a,
		 double[] c,
		 int i,
		 int n)
		{
		int j;

		// First, compute vector of second derivatives of u, upp, by solving
		// this tridiagonal linear system (illustrated for n=6):
		//     [1 0        ]   [upp0]   [0                  ]
		//     [1 4 1      ]   [upp1]   [6 u0 - 12 u1 + 6 u2]
		//     [  1 4 1    ] x [upp2] = [6 u1 - 12 u2 + 6 u3]
		//     [    1 4 1  ]   [upp3]   [6 u2 - 12 u3 + 6 u4]
		//     [      1 4 1]   [upp4]   [6 u3 - 12 u4 + 6 u5]
		//     [        0 1]   [upp5]   [0                  ]

		// Allocate storage for tridiagonal matrix, solution vector, and right-
		// hand side vector.
		double[] f = new double [n-1];
		double[] d = new double [n];
		double[] e = new double [n-1];
		double[] upp = new double [n];
		double[] rhs = new double [n];

		// Fill in tridiagonal matrix and right-hand side vector.
		Arrays.fill (f, 0, n-2, 1.0);
		d[0] = 1.0; Arrays.fill (d, 1, n-1, 4.0); d[n-1] = 1.0;
		Arrays.fill (e, 1, n-1, 1.0);
		for (j = 1; j < n-1; ++ j)
			{
			rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*u[i+j+1];
			}

		// Solve the linear system.
		Tridiagonal.solve (d, e, f, rhs, upp);

		// Lastly, compute Bezier control point vectors a and c from u and upp.
		for (j = 0; j < n-1; ++ j)
			{
			a[i+j] = 2.0*u[i+j]/3.0 + u[i+j+1]/3.0 - upp[j]/9.0 - upp[j+1]/18.0;
			c[i+j] = u[i+j]/3.0 + 2.0*u[i+j+1]/3.0 - upp[j]/18.0 - upp[j+1]/9.0;
			}

		}

	/**
	 * Compute the Bezier control point coordinates for an open smooth curve
	 * with a specified initial direction and zero final curvature.
	 *
	 * @param  uInitial  Specifies initial direction as a straight line from
	 *                   <TT>uInitial</TT> to <TT>u[i]</TT>.
	 * @param  u  An input array of coordinates. Elements at indexes <TT>i</TT>
	 *            .. <TT>i+n-1</TT> are used.
	 * @param  a  An output array of coordinates for the first Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  c  An output array of coordinates for the second Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  i  Index of first element used in input and output arrays.
	 * @param  n  Number of elements used in input array; one fewer elements
	 *            used in output arrays. Must be greater than or equal to 2.
	 *
 	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &lt; 2.
 	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the control points cannot be
	 *     calculated.
	 */
	public static void computeBezierOpen
		(double uInitial,
		 double[] u,
		 double[] a,
		 double[] c,
		 int i,
		 int n)
		{
		int j;

		// First, compute vector of second derivatives of u, upp, by solving
		// this tridiagonal linear system (illustrated for n=6):
		//     [2 1        ]   [upp0]   [6 uInitial - 12 u0 + 6 u1]
		//     [1 4 1      ]   [upp1]   [6 u0 - 12 u1 + 6 u2      ]
		//     [  1 4 1    ] x [upp2] = [6 u1 - 12 u2 + 6 u3      ]
		//     [    1 4 1  ]   [upp3]   [6 u2 - 12 u3 + 6 u4      ]
		//     [      1 4 1]   [upp4]   [6 u3 - 12 u4 + 6 u5      ]
		//     [        0 1]   [upp5]   [0                        ]

		// Allocate storage for tridiagonal matrix, solution vector, and right-
		// hand side vector.
		double[] f = new double [n-1];
		double[] d = new double [n];
		double[] e = new double [n-1];
		double[] upp = new double [n];
		double[] rhs = new double [n];

		// Fill in tridiagonal matrix and right-hand side vector.
		Arrays.fill (f, 0, n-2, 1.0);
		d[0] = 2.0; Arrays.fill (d, 1, n-1, 4.0); d[n-1] = 1.0;
		Arrays.fill (e, 0, n-1, 1.0);
		rhs[0] = 6.0*uInitial - 12.0*u[i] + 6.0*u[i+1];
		for (j = 1; j < n-1; ++ j)
			{
			rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*u[i+j+1];
			}

		// Solve the linear system.
		Tridiagonal.solve (d, e, f, rhs, upp);

		// Lastly, compute Bezier control point vectors a and c from u and upp.
		for (j = 0; j < n-1; ++ j)
			{
			a[i+j] = 2.0*u[i+j]/3.0 + u[i+j+1]/3.0 - upp[j]/9.0 - upp[j+1]/18.0;
			c[i+j] = u[i+j]/3.0 + 2.0*u[i+j+1]/3.0 - upp[j]/18.0 - upp[j+1]/9.0;
			}

		}

	/**
	 * Compute the Bezier control point coordinates for an open smooth curve
	 * with zero initial curvature and a specified final direction.
	 *
	 * @param  u  An input array of coordinates. Elements at indexes <TT>i</TT>
	 *            .. <TT>i+n-1</TT> are used.
	 * @param  uFinal    Specifies final direction as a straight line from
	 *                   <TT>u[i+n-1]</TT> to <TT>uFinal</TT>.
	 * @param  a  An output array of coordinates for the first Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  c  An output array of coordinates for the second Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  i  Index of first element used in input and output arrays.
	 * @param  n  Number of elements used in input array; one fewer elements
	 *            used in output arrays. Must be greater than or equal to 2.
	 *
 	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &lt; 2.
 	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the control points cannot be
	 *     calculated.
	 */
	public static void computeBezierOpen
		(double[] u,
		 double uFinal,
		 double[] a,
		 double[] c,
		 int i,
		 int n)
		{
		int j;

		// First, compute vector of second derivatives of u, upp, by solving
		// this tridiagonal linear system (illustrated for n=6):
		//     [1 0        ]   [upp0]   [0                      ]
		//     [1 4 1      ]   [upp1]   [6 u0 - 12 u1 + 6 u2    ]
		//     [  1 4 1    ] x [upp2] = [6 u1 - 12 u2 + 6 u3    ]
		//     [    1 4 1  ]   [upp3]   [6 u2 - 12 u3 + 6 u4    ]
		//     [      1 4 1]   [upp4]   [6 u3 - 12 u4 + 6 u5    ]
		//     [        1 2]   [upp5]   [6 u4 - 12 u5 + 6 uFinal]

		// Allocate storage for tridiagonal matrix, solution vector, and right-
		// hand side vector.
		double[] f = new double [n-1];
		double[] d = new double [n];
		double[] e = new double [n-1];
		double[] upp = new double [n];
		double[] rhs = new double [n];

		// Fill in tridiagonal matrix and right-hand side vector.
		Arrays.fill (f, 0, n-1, 1.0);
		d[0] = 1.0; Arrays.fill (d, 1, n-1, 4.0); d[n-1] = 2.0;
		Arrays.fill (e, 1, n-1, 1.0);
		for (j = 1; j < n-1; ++ j)
			{
			rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*u[i+j+1];
			}
		rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*uFinal;

		// Solve the linear system.
		Tridiagonal.solve (d, e, f, rhs, upp);

		// Lastly, compute Bezier control point vectors a and c from u and upp.
		for (j = 0; j < n-1; ++ j)
			{
			a[i+j] = 2.0*u[i+j]/3.0 + u[i+j+1]/3.0 - upp[j]/9.0 - upp[j+1]/18.0;
			c[i+j] = u[i+j]/3.0 + 2.0*u[i+j+1]/3.0 - upp[j]/18.0 - upp[j+1]/9.0;
			}

		}

	/**
	 * Compute the Bezier control point coordinates for an open smooth curve
	 * with a specified initial direction and a specified final direction.
	 *
	 * @param  uInitial  Specifies initial direction as a straight line from
	 *                   <TT>uInitial</TT> to <TT>u[i]</TT>.
	 * @param  u  An input array of coordinates. Elements at indexes <TT>i</TT>
	 *            .. <TT>i+n-1</TT> are used.
	 * @param  uFinal    Specifies final direction as a straight line from
	 *                   <TT>u[i+n-1]</TT> to <TT>uFinal</TT>.
	 * @param  a  An output array of coordinates for the first Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  c  An output array of coordinates for the second Bezier control
	 *            points. Elements at indexes <TT>i</TT> .. <TT>i+n-2</TT> are
	 *            used.
	 * @param  i  Index of first element used in input and output arrays.
	 * @param  n  Number of elements used in input array; one fewer elements
	 *            used in output arrays. Must be greater than or equal to 2.
	 *
 	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &lt; 2.
 	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the control points cannot be
	 *     calculated.
	 */
	public static void computeBezierOpen
		(double uInitial,
		 double[] u,
		 double uFinal,
		 double[] a,
		 double[] c,
		 int i,
		 int n)
		throws DomainException
		{
		int j;

		// First, compute vector of second derivatives of u, upp, by solving
		// this tridiagonal linear system (illustrated for n=6):
		//     [2 1        ]   [upp0]   [6 uInitial - 12 u0 + 6 u1]
		//     [1 4 1      ]   [upp1]   [6 u0 - 12 u1 + 6 u2      ]
		//     [  1 4 1    ] x [upp2] = [6 u1 - 12 u2 + 6 u3      ]
		//     [    1 4 1  ]   [upp3]   [6 u2 - 12 u3 + 6 u4      ]
		//     [      1 4 1]   [upp4]   [6 u3 - 12 u4 + 6 u5      ]
		//     [        1 2]   [upp5]   [6 u4 - 12 u5 + 6 uFinal  ]

		// Allocate storage for tridiagonal matrix, solution vector, and right-
		// hand side vector.
		double[] f = new double [n-1];
		double[] d = new double [n];
		double[] e = new double [n-1];
		double[] upp = new double [n];
		double[] rhs = new double [n];

		// Fill in tridiagonal matrix and right-hand side vector.
		Arrays.fill (f, 1.0);
		d[0] = 2.0; Arrays.fill (d, 1, n-1, 4.0); d[n-1] = 2.0;
		Arrays.fill (e, 1.0);
		rhs[0] = 6.0*uInitial - 12.0*u[i] + 6.0*u[i+1];
		for (j = 1; j < n-1; ++ j)
			{
			rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*u[i+j+1];
			}
		rhs[j] = 6.0*u[i+j-1] - 12.0*u[i+j] + 6.0*uFinal;

		// Solve the linear system.
		Tridiagonal.solve (d, e, f, rhs, upp);

		// Lastly, compute Bezier control point vectors a and c from u and upp.
		for (j = 0; j < n-1; ++ j)
			{
			a[i+j] = 2.0*u[i+j]/3.0 + u[i+j+1]/3.0 - upp[j]/9.0 - upp[j+1]/18.0;
			c[i+j] = u[i+j]/3.0 + 2.0*u[i+j+1]/3.0 - upp[j]/18.0 - upp[j+1]/9.0;
			}
		}

	}
