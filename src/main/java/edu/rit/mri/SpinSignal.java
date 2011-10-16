//******************************************************************************
//
// File:    SpinSignal.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SpinSignal
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

package edu.rit.mri;

/**
 * Class SpinSignal provides static methods for computing the spin signal model
 * function.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2008
 */
public class SpinSignal
	{

// Prevent construction.

	private SpinSignal()
		{
		}

// Exported operations.

	/**
	 * Compute the spin signal for the given spin relaxation rate and time. A
	 * spin density of <I>&rho;</I> = 1 is used. The formula is
	 * <CENTER>
	 * <I>S</I>(<I>t</I>)&emsp;=&emsp;<I>&rho;</I> [1 &minus; 2 exp(&minus;<I>R</I> <I>t</I>)]
	 * </CENTER>
	 *
	 * @param  R  Spin relaxation rate, <I>R</I>.
	 * @param  t  Time, <I>t</I>.
	 *
	 * @return  Spin signal, <I>S</I>(<I>t</I>).
	 */
	public static double S
		(double R,
		 double t)
		{
		return 1.0 - 2.0*Math.exp(-R*t);
		}

	/**
	 * Compute the spin signal for the given spin density, spin relaxation rate,
	 * and time. The formula is
	 * <CENTER>
	 * <I>S</I>(<I>t</I>)&emsp;=&emsp;<I>&rho;</I> [1 &minus; 2 exp(&minus;<I>R</I> <I>t</I>)]
	 * </CENTER>
	 *
	 * @param  rho  Spin density, <I>&rho;</I>.
	 * @param  R    Spin relaxation rate, <I>R</I>.
	 * @param  t    Time, <I>t</I>.
	 *
	 * @return  Spin signal, <I>S</I>(<I>t</I>).
	 */
	public static double S
		(double rho,
		 double R,
		 double t)
		{
		return rho*(1.0 - 2.0*Math.exp(-R*t));
		}

	/**
	 * Compute the partial derivative of the spin signal with respect to spin
	 * density for the given spin relaxation rate and time. The formula is
	 * <CENTER>
	 * &part;<I>S</I>(<I>t</I>)&nbsp;&frasl;&nbsp;&part;<I>&rho;</I>&emsp;=&emsp;1 &minus; 2 exp(&minus;<I>R</I> <I>t</I>)
	 * </CENTER>
	 *
	 * @param  R    Spin relaxation rate, <I>R</I>.
	 * @param  t    Time, <I>t</I>.
	 *
	 * @return  Partial derivative, &part;<I>S</I>(<I>t</I>)&nbsp;&frasl;&nbsp;&part;<I>&rho;</I>.
	 */
	public static double dSdrho
		(double R,
		 double t)
		{
		return 1.0 - 2.0*Math.exp(-R*t);
		}

	/**
	 * Compute the partial derivative of the spin signal with respect to spin
	 * relaxation rate for the given spin density, spin relaxation rate, and
	 * time. The formula is
	 * <CENTER>
	 * &part;<I>S</I>(<I>t</I>)&nbsp;&frasl;&nbsp;&part;<I>R</I>&emsp;=&emsp;2 <I>&rho;</I> <I>t</I> exp(&minus;<I>R</I> <I>t</I>)
	 * </CENTER>
	 *
	 * @param  rho  Spin density, <I>&rho;</I>.
	 * @param  R    Spin relaxation rate, <I>R</I>.
	 * @param  t    Time, <I>t</I>.
	 *
	 * @return  Partial derivative, &part;<I>S</I>(<I>t</I>)&nbsp;&frasl;&nbsp;&part;<I>R</I>.
	 */
	public static double dSdR
		(double rho,
		 double R,
		 double t)
		{
		return 2.0*rho*t*Math.exp(-R*t);
		}

	}
