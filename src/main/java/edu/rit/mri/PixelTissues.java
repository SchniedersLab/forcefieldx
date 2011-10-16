//******************************************************************************
//
// File:    PixelTissues.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PixelTissues
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

import edu.rit.numeric.ArraySeries;
import edu.rit.numeric.Series;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.util.List;

/**
 * Class PixelTissues encapsulates the results of a spin relaxometry analysis
 * on one pixel of a magnetic resonance image. This includes:
 * <UL>
 * <LI>
 * File index.
 * <LI>
 * Pixel index.
 * <LI>
 * Pixel's computed spin density <I>&rho;</I> for each tissue.
 * <LI>
 * Pixel's computed spin-lattice relaxation rate <I>R</I>1 for each tissue.
 * </UL>
 * <P>
 * Each spin density and spin-lattice relaxation rate is stored as a value of
 * type <TT>double</TT>.
 * <P>
 * Operations are provided to read a pixel tissues object from a DataInputStream
 * and write a pixel tissues object to a DataOutputStream. The format is:
 * <UL>
 * <LI>
 * Number of tissues, <I>L</I> (2-byte short).
 * <LI>
 * <I>L</I> spin density values, each an 8-byte double.
 * <LI>
 * <I>L</I> spin-lattice relaxation rate values, each an 8-byte double.
 * </UL>
 * <P>
 * Operations are provided to read a pixel tissues object from an
 * ObjectInputStream and write a pixel tissues object to an ObjectOutputStream.
 * The format is:
 * <UL>
 * <LI>
 * File index (4-byte int).
 * <LI>
 * Pixel index (4-byte int).
 * <LI>
 * Number of tissues, <I>L</I> (2-byte short).
 * <LI>
 * <I>L</I> spin density values, each an 8-byte double.
 * <LI>
 * <I>L</I> spin-lattice relaxation rate values, each an 8-byte double.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class PixelTissues
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 5354464511124925380L;

	// File index.
	private int fileIndex;

	// Pixel index.
	private int pixelIndex;

	// Array of computed spin densities.
	private double[] rho;

	// Array of computed spin-lattice relaxation rates.
	private double[] R1;

// Exported constructors.

	/**
	 * Construct a new, uninitialized pixel tissues object. This constructor is
	 * for use only by object deserialization.
	 */
	public PixelTissues()
		{
		}

	/**
	 * Construct a new pixel tissues object. If a solution could not be computed
	 * for this pixel, <TT>rho</TT> and <TT>R1</TT> are null.
	 *
	 * @param  fileIndex   File index.
	 * @param  pixelIndex  Pixel index.
	 * @param  rho_list    List containing computed spin densities, or null.
	 * @param  R1_list     List containing computed spin-lattice relaxation
	 *                     rates, or null.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>index</TT> &lt; 0. Thrown if
	 *     <TT>rho_list</TT> and <TT>R1_list</TT> are not the same length.
	 */
	public PixelTissues
		(int fileIndex,
		 int pixelIndex,
		 List<Double> rho_list,
		 List<Double> R1_list)
		{
		fileIndex (fileIndex);
		pixelIndex (pixelIndex);
		tissues (rho_list, R1_list);
		}

// Exported operations.

	/**
	 * Get the file index.
	 *
	 * @return  File index.
	 */
	public int fileIndex()
		{
		return fileIndex;
		}

	/**
	 * Specify the file index.
	 *
	 * @param  index  File index.
	 */
	public void fileIndex
		(int index)
		{
		fileIndex = index;
		}

	/**
	 * Get the pixel index.
	 *
	 * @return  Pixel index.
	 */
	public int pixelIndex()
		{
		return pixelIndex;
		}

	/**
	 * Specify the pixel index.
	 *
	 * @param  index  Pixel index.
	 */
	public void pixelIndex
		(int index)
		{
		pixelIndex = index;
		}

	/**
	 * Specify the tissues for this pixel. If a solution could not be computed
	 * for this pixel, <TT>rho</TT> and <TT>R1</TT> are null.
	 *
	 * @param  rho_list  List containing computed spin densities, or null.
	 * @param  R1_list   List containing computed spin-lattice relaxation rates,
	 *                   or null.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>rho_list</TT> and
	 *     <TT>R1_list</TT> are not the same length.
	 */
	public void tissues
		(List<Double> rho_list,
		 List<Double> R1_list)
		{
		// Verify preconditions.
		int rhosize = rho_list == null ? 0 : rho_list.size();
		int R1size = R1_list == null ? 0 : R1_list.size();
		if (rhosize != R1size)
			{
			throw new IllegalArgumentException
				("PixelTissues.tissues(): rho_list length (= "+rhosize+
				 ") != R1_list length (= "+R1size+"), illegal");
			}

		if (rho_list != null && R1_list != null)
			{
			int L = rho_list.size();
			this.rho = new double [L];
			this.R1 = new double [L];
			for (int i = 0; i < L; ++ i)
				{
				this.rho[i] = rho_list.get(i);
				this.R1[i] = R1_list.get(i);
				}
			}
		else
			{
			this.rho = null;
			this.R1 = null;
			}
		}

	/**
	 * Get the number of tissues computed for this pixel. If the analysis could
	 * not find a solution for this pixel, 0 is returned.
	 *
	 * @return  Number of tissues, <I>L</I>.
	 */
	public int numTissues()
		{
		return rho == null ? 0 : rho.length;
		}

	/**
	 * Get the given tissue's spin density.
	 *
	 * @param  i  Tissue index, 0 &le; <I>i</I> &le; <I>L</I>&minus;1.
	 *
	 * @return  Spin density, <I>&rho;</I><SUB><I>i</I></SUB>.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <I>i</I> is out of bounds.
	 */
	public double rho
		(int i)
		{
		return rho[i];
		}

	/**
	 * Get the given tissue's spin-lattice relaxation rate.
	 *
	 * @param  i  Tissue index, 0 &le; <I>i</I> &le; <I>L</I>&minus;1.
	 *
	 * @return  Spin-lattice relaxation rate, <I>R</I>1<SUB><I>i</I></SUB>.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <I>i</I> is out of bounds.
	 */
	public double R1
		(int i)
		{
		return R1[i];
		}

	/**
	 * Get this pixel's computed spin signal at the given time. The formula is:
	 * <CENTER>
	 * <I>S</I>(<I>t</I>)&emsp;=&emsp;<B>&Sigma;</B><SUB><I>i</I></SUB>&emsp;<I>&rho;</I><SUB><I>i</I></SUB>&nbsp;[1&nbsp;&minus;&nbsp;2&nbsp;exp(&minus;<I>R</I>1<SUB><I>i</I></SUB>&nbsp;<I>t</I>)]
	 * </CENTER>
	 * If the analysis could not find a solution for this pixel, 0 is returned.
	 *
	 * @param  t  Time.
	 *
	 * @return  Spin signal, <I>S</I>(<I>t</I>).
	 */
	public double S
		(double t)
		{
		if (rho == null) return 0.0;
		double sum = 0.0;
		for (int i = 0; i < rho.length; ++ i)
			{
			sum += SpinSignal.S (rho[i], R1[i], t);
			}
		return sum;
		}

	/**
	 * Get a series containing this pixel's computed spin signal for the given
	 * time series. If the analysis could not find a solution for this pixel,
	 * null is returned.
	 *
	 * @param  t_series  Series of times, <I>t</I><SUB><I>j</I></SUB>.
	 *
	 * @return  Series of spin signals, <I>S</I>(<I>t</I><SUB><I>j</I></SUB>),
	 *          or null.
	 */
	public Series S_series
		(Series t_series)
		{
		if (rho == null) return null;
		final int M = t_series.length();
		final double[] S_computed = new double [M];
		for (int j = 0; j < M; ++ j)
			{
			double t_j = t_series.x(j);
			double sum = 0.0;
			for (int i = 0; i < rho.length; ++ i)
				{
				sum += SpinSignal.S (rho[i], R1[i], t_j);
				}
			S_computed[j] = sum;
			}
		return new ArraySeries (S_computed);
		}

	/**
	 * Write this pixel tissues object to the given data output stream.
	 * <P>
	 * <I>Note:</I> The file index and pixel index are not written.
	 *
	 * @param  out  Data output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void write
		(DataOutput out)
		throws IOException
		{
		int L = rho == null ? 0 : rho.length;
		out.writeShort ((short) L);
		for (int i = 0; i < L; ++ i)
			{
			out.writeDouble (rho[i]);
			}
		for (int i = 0; i < L; ++ i)
			{
			out.writeDouble (R1[i]);
			}
		}

	/**
	 * Read this pixel tissues object from the given data input stream.
	 * <P>
	 * <I>Note:</I> The file index and pixel index are not read.
	 *
	 * @param  in  Data input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void read
		(DataInput in)
		throws IOException
		{
		int L = in.readShort();
		if (L == 0)
			{
			rho = null;
			R1 = null;
			}
		else
			{
			rho = new double [L];
			R1 = new double [L];
			for (int i = 0; i < L; ++ i)
				{
				rho[i] = in.readDouble();
				}
			for (int i = 0; i < L; ++ i)
				{
				R1[i] = in.readDouble();
				}
			}
		}

	/**
	 * Write this pixel tissues object to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeInt (fileIndex);
		out.writeInt (pixelIndex);
		write (out);
		}

	/**
	 * Read this pixel tissues object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		fileIndex = in.readInt();
		pixelIndex = in.readInt();
		read (in);
		}

	}
