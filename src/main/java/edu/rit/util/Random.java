//******************************************************************************
//
// File:    Random.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Random
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

package edu.rit.util;

import java.io.Serializable;
import java.lang.reflect.Constructor;

/**
 * Class Random is the abstract base class for a pseudorandom number generator
 * (PRNG) designed for use in parallel scientific programming. It differs from
 * class java.util.Random in the following ways:
 * <UL>
 * <LI>
 * Instances can be created by calling a static factory method. The factory
 * method can take the name of a subclass as an argument and return an instance
 * of that subclass of class Random. This makes it easier to write a program
 * that can substitute a different PRNG algorithm at run time. Class {@linkplain
 * DefaultRandom} provides a default PRNG algorithm.
 * <BR>&nbsp;
 * <LI>
 * Whereas class java.util.Random generates 48-bit numbers under the hood, class
 * Random generates 64-bit numbers. This makes class Random faster than class
 * java.util.Random for some operations, notably the <TT>nextDouble()</TT>
 * method.
 * <BR>&nbsp;
 * <LI>
 * Whereas the PRNG algorithm in class java.util.Random has a period of about
 * 2<SUP>48</SUP>, the default PRNG algorithm in class {@linkplain
 * DefaultRandom} has a period of about 2<SUP>64</SUP>. This lets a parallel
 * program scale up to larger problem sizes without exhausting the PRNG's
 * period.
 * <BR>&nbsp;
 * <LI>
 * To support the <I>leapfrogging</I> and <I>sequence splitting</I> techniques
 * for generating pseudorandom numbers in multiple parallel threads or
 * processes, class Random includes methods for skipping ahead in the sequence
 * of generated numbers, without having to generate all the intermediate
 * numbers.
 * <BR>&nbsp;
 * <LI>
 * To avoid unnecessary thread synchronization overhead, class Random is not
 * multiple thread safe. It assumes that the surrounding program ensures that
 * only one thread at a time calls methods on an instance.
 * </UL>
 * <P>
 * Each method for generating a number comes in two versions; for example:
 * <PRE>
 *     public double nextDouble();
 *     public double nextDouble (long skip);
 * </PRE>
 * Calling the second version with an argument <TT>skip</TT> is equivalent to
 * calling the first version <TT>skip</TT> times:
 * <PRE>
 *     Random prng1 = Random.getInstance (1234L);
 *     double x = prng1.nextDouble (1000);
 *     Random prng2 = Random.getInstance (1234L);
 *     double y;
 *     for (int i = 0; i &lt; 1000; ++ i)
 *         y = prng2.nextDouble();
 * </PRE>
 * At the end of the above code fragment, <TT>x</TT> and <TT>y</TT> will have
 * the same value. However, calling <TT>nextDouble(1000)</TT> will typically be
 * much faster than calling <TT>nextDouble()</TT> 1000 times. Conversely,
 * <TT>nextDouble(1)</TT> is equivalent to <TT>nextDouble()</TT>, but the former
 * will typically be slower than the latter.
 * <P>
 * An instance of class Random can be serialized; for example, to checkpoint the
 * state of the PRNG into a file and restore its state later.
 * <P>
 * The design of class Random is inspired in part by Coddington's JAPARA
 * library. For further information, see P. Coddington and A. Newell, JAPARA --
 * A Java random number generator library for high-performance computing, in
 * <I>Proceedings of the 18th IEEE International Parallel and Distributed
 * Processing Symposium (IPDPS'04),</I> April 26-30, 2004, page 156.
 *
 * @author  Alan Kaminsky
 * @version 17-Aug-2009
 */
public abstract class Random
	implements Serializable
	{

// Hidden constants.

	// 1/2^64, as a float and as a double.
	private static double D_2_POW_NEG_64;
	private static float  F_2_POW_NEG_64;
	static
		{
		double x = 2.0;
		x *= x; // 2^2
		x *= x; // 2^4
		x *= x; // 2^8
		x *= x; // 2^16
		x *= x; // 2^32
		x *= x; // 2^64
		D_2_POW_NEG_64 = 1.0 / x;
		F_2_POW_NEG_64 = (float) D_2_POW_NEG_64;
		}

// Hidden constructors.

	/**
	 * Construct a new PRNG.
	 */
	protected Random()
		{
		}

// Exported operations.

	/**
	 * Construct a new PRNG with the given seed using the default algorithm.
	 * <P>
	 * If the <TT>"pj.prng"</TT> Java system property is specified, it gives the
	 * fully-qualified class name of the default PRNG class that the
	 * <TT>getInstance()</TT> method will construct. Specifying the
	 * <TT>"pj.prng"</TT> property will substitute a different PRNG algorithm
	 * into a program without needing to recompile.
	 * <P>
	 * If the <TT>"pj.prng"</TT> Java system property is not specified, the
	 * <TT>getInstance()</TT> method will return an instance of class
	 * {@linkplain DefaultRandom}.
	 * <P>
	 * You can specify the <TT>"pj.prng"</TT> property on the Java command line
	 * like this:
	 * <P>
	 * <TT>&nbsp;&nbsp;&nbsp;&nbsp;java -Dpj.prng=MyOwnPrngClass . . .</TT>
	 * <P>
	 * <I>Note:</I> Depending on the PRNG algorithm, certain seed values may not
	 * be allowed. See the PRNG algorithm subclass for further information.
	 *
	 * @param  seed  Seed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the PRNG algorithm does not allow the
	 *     given seed value.
	 * @exception  TypeNotPresentException
	 *     (unchecked exception) Thrown if a PRNG instance could not be
	 *     constructed. The chained exception gives further information about
	 *     the problem.
	 */
	public static Random getInstance
		(long seed)
		{
		String algorithm = System.getProperty ("pj.prng");
		return
			algorithm == null ?
				new DefaultRandom (seed) :
				getInstance (seed, algorithm);
		}

	/**
	 * Construct a new PRNG with the given seed using the given algorithm.
	 * <P>
	 * <I>Note:</I> Depending on the PRNG algorithm, certain seed values may not
	 * be allowed. See the PRNG algorithm subclass for further information.
	 *
	 * @param  seed       Seed.
	 * @param  algorithm  Fully-qualified name of the class to construct. This
	 *                    must be a subclass of class Random.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the PRNG algorithm does not allow the
	 *     given seed value.
	 * @exception  TypeNotPresentException
	 *     (unchecked exception) Thrown if a PRNG instance could not be
	 *     constructed. The chained exception gives further information about
	 *     the problem.
	 */
	public static Random getInstance
		(long seed,
		 String algorithm)
		{
		try
			{
			Class<?> cl = Class.forName (algorithm);
			Constructor<?> ctor = cl.getConstructor (Long.TYPE);
			return (Random) ctor.newInstance (seed);
			}
		catch (Exception exc)
			{
			throw new TypeNotPresentException (algorithm, exc);
			}
		}

	/**
	 * Set this PRNG's seed.
	 * <P>
	 * <I>Note:</I> Depending on the PRNG algorithm, certain seed values may not
	 * be allowed. See the PRNG algorithm subclass for further information.
	 *
	 * @param  seed  Seed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the PRNG algorithm does not allow the
	 *     given seed value.
	 */
	public abstract void setSeed
		(long seed);

	/**
	 * Skip one position ahead in this PRNG's sequence.
	 */
	public void skip()
		{
		next();
		}

	/**
	 * Skip the given number of positions ahead in this PRNG's sequence. If
	 * <TT>skip</TT> &lt;= 0, the <TT>skip()</TT> method does nothing.
	 *
	 * @param  skip  Number of positions to skip.
	 */
	public void skip
		(long skip)
		{
		if (skip > 0L) next (skip);
		}

	/**
	 * Return the Boolean value from the next pseudorandom value in this PRNG's
	 * sequence. With a probability of 0.5 <TT>true</TT> is returned, with a
	 * probability of 0.5 <TT>false</TT> is returned.
	 *
	 * @return  Boolean value.
	 */
	public boolean nextBoolean()
		{
		// Use the high-order (sign) bit of the 64-bit random value.
		return next() >= 0L;
		}

	/**
	 * Return the Boolean value from the pseudorandom value the given number of
	 * positions ahead in this PRNG's sequence. With a probability of 0.5
	 * <TT>true</TT> is returned, with a probability of 0.5 <TT>false</TT> is
	 * returned.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Boolean value.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public boolean nextBoolean
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextBoolean(): skip = " + skip + " illegal");
			}
			
		// Use the high-order (sign) bit of the 64-bit random value.
		return next (skip) >= 0L;
		}

	/**
	 * Return the byte value from the next pseudorandom value in this PRNG's
	 * sequence. Each value in the range &minus;128 through 127 is returned with
	 * a probability of 1/2<SUP>8</SUP>.
	 *
	 * @return  Byte value in the range &minus;128 through 127 inclusive.
	 */
	public byte nextByte()
		{
		return (byte) next();
		}

	/**
	 * Return the byte value from the next pseudorandom value the given number
	 * of positions ahead in this PRNG's sequence. Each value in the range
	 * &minus;128 through 127 is returned with a probability of 1/2<SUP>8</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Byte value in the range &minus;128 through 127 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public byte nextByte
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextByte(): skip = " + skip + " illegal");
			}
		return (byte) next (skip);
		}

	/**
	 * Return the unsigned byte value from the next pseudorandom value in this
	 * PRNG's sequence. Each value in the range 0 through 255 is returned with a
	 * probability of 1/2<SUP>8</SUP>.
	 *
	 * @return  Unsigned byte value (as an <TT>int</TT>) in the range 0 through
	 *          255 inclusive.
	 */
	public int nextUnsignedByte()
		{
		return (int) (next() & 0xFFL);
		}

	/**
	 * Return the unsigned byte value from the next pseudorandom value the given
	 * number of positions ahead in this PRNG's sequence. Each value in the
	 * range 0 through 255 is returned with a probability of 1/2<SUP>8</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Unsigned byte value (as an <TT>int</TT>) in the range 0 through
	 *          255 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public int nextUnsignedByte
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextUnsignedByte(): skip = " + skip + " illegal");
			}
		return (int) (next (skip) & 0xFFL);
		}

	/**
	 * Return the character value from the next pseudorandom value in this
	 * PRNG's sequence. Each value in the range <TT>'&#92;u0000'</TT> through
	 * <TT>'&#92;uFFFF'</TT> is returned with a probability of 1/2<SUP>16</SUP>.
	 *
	 * @return  Character value in the range <TT>'&#92;u0000'</TT> through
	 *          <TT>'&#92;uFFFF'</TT> inclusive.
	 */
	public char nextCharacter()
		{
		return (char) next();
		}

	/**
	 * Return the character value from the next pseudorandom value the given
	 * number of positions ahead in this PRNG's sequence. Each value in the
	 * range <TT>'&#92;u0000'</TT> through <TT>'&#92;uFFFF'</TT> is returned
	 * with a probability of 1/2<SUP>16</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Character value in the range <TT>'&#92;u0000'</TT> through
	 *          <TT>'&#92;uFFFF'</TT> inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public int nextCharacter
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextCharacter(): skip = " + skip + " illegal");
			}
		return (char) next (skip);
		}

	/**
	 * Return the short value from the next pseudorandom value in this PRNG's
	 * sequence. Each value in the range &minus;32768 through 32767 is returned
	 * with a probability of 1/2<SUP>16</SUP>.
	 *
	 * @return  Short value in the range &minus;32768 through 32767 inclusive.
	 */
	public short nextShort()
		{
		return (short) next();
		}

	/**
	 * Return the short value from the next pseudorandom value the given number
	 * of positions ahead in this PRNG's sequence. Each value in the range
	 * &minus;32768 through 32767 is returned with a probability of
	 * 1/2<SUP>16</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Short value in the range &minus;32768 through 32767 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public short nextShort
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextShort(): skip = " + skip + " illegal");
			}
		return (short) next (skip);
		}

	/**
	 * Return the unsigned short value from the next pseudorandom value in this
	 * PRNG's sequence. Each value in the range 0 through 65535 is returned with
	 * a probability of 1/2<SUP>16</SUP>.
	 *
	 * @return  Unsigned short value (as an <TT>int</TT>) in the range 0 through
	 *          65535 inclusive.
	 */
	public int nextUnsignedShort()
		{
		return (int) (next() & 0xFFFFL);
		}

	/**
	 * Return the unsigned short value from the next pseudorandom value the
	 * given number of positions ahead in this PRNG's sequence. Each value in
	 * the range 0 through 65535 is returned with a probability of
	 * 1/2<SUP>16</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Unsigned short value (as an <TT>int</TT>) in the range 0 through
	 *          65535 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public int nextUnsignedShort
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextUnsignedShort(): skip = " + skip + " illegal");
			}
		return (int) (next (skip) & 0xFFFFL);
		}

	/**
	 * Return the integer value from the next pseudorandom value in this PRNG's
	 * sequence. Each value in the range &minus;2147483648 through 2147483647 is
	 * returned with a probability of 1/2<SUP>32</SUP>.
	 *
	 * @return  Integer value in the range &minus;2147483648 through 2147483647
	 *          inclusive.
	 */
	public int nextInteger()
		{
		return (int) next();
		}

	/**
	 * Return the integer value from the next pseudorandom value the given
	 * number of positions ahead in this PRNG's sequence. Each value in the
	 * range &minus;2147483648 through 2147483647 is returned with a probability
	 * of 1/2<SUP>32</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Integer value in the range &minus;2147483648 through 2147483647
	 *          inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public int nextInteger
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextInteger(): skip = " + skip + " illegal");
			}
		return (int) next (skip);
		}

	/**
	 * Return the integer value in the given range from the next pseudorandom
	 * value in this PRNG's sequence. Each value in the range 0 through
	 * <I>N</I>&minus;1 is returned with a probability of 1/<I>N</I>.
	 *
	 * @param  n  Range of values to return.
	 *
	 * @return  Integer value in the range 0 through <I>N</I>&minus;1 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>N</I> &lt;= 0.
	 */
	public int nextInt
		(int n)
		{
		if (n <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextInt(): n = " + n + " illegal");
			}
		return (int) Math.floor (nextDouble() * n);
		}

	/**
	 * Return the integer value in the given range from the pseudorandom value
	 * the given number of positions ahead in this PRNG's sequence. Each value
	 * in the range 0 through <I>N</I>&minus;1 is returned with a probability of
	 * 1/<I>N</I>.
	 *
	 * @param  n     Range of values to return.
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Integer value in the range 0 through <I>N</I>&minus;1 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>N</I> &lt;= 0. Thrown if
	 *     <TT>skip</TT> &lt;= 0.
	 */
	public int nextInt
		(int n,
		 long skip)
		{
		if (n <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextInt(): n = " + n + " illegal");
			}
		return (int) Math.floor (nextDouble (skip) * n);
		}

	/**
	 * Return the long value from the next pseudorandom value in this PRNG's
	 * sequence. Each value in the range &minus;9223372036854775808 through
	 * 9223372036854775807 is returned with a probability of 1/2<SUP>64</SUP>.
	 *
	 * @return  Long value in the range &minus;9223372036854775808 through 
	 *          9223372036854775807 inclusive.
	 */
	public long nextLong()
		{
		return next();
		}

	/**
	 * Return the long value from the next pseudorandom value the given number
	 * of positions ahead in this PRNG's sequence. Each value in the range
	 * &minus;9223372036854775808 through 9223372036854775807 is returned with a
	 * probability of 1/2<SUP>64</SUP>.
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Long value in the range &minus;9223372036854775808 through 
	 *          9223372036854775807 inclusive.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public long nextLong
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextLong(): skip = " + skip + " illegal");
			}
		return next (skip);
		}

	/**
	 * Return the single precision floating point value from the next
	 * pseudorandom value in this PRNG's sequence. The returned numbers have a
	 * uniform distribution in the range 0.0 (inclusive) to 1.0 (exclusive).
	 *
	 * @return  Float value.
	 */
	public float nextFloat()
		{
		// Next random number is in the range -2^63 .. +2^63 - 1.
		// Divide by 2^64 yielding a number in the range -0.5 .. +0.5.
		// Add 0.5 yielding a number in the range 0.0 .. 1.0.
		return (float) (next()) * F_2_POW_NEG_64 + 0.5f;
		}

	/**
	 * Return the single precision floating point value from the pseudorandom
	 * value the given number of positions ahead in this PRNG's sequence. The
	 * returned numbers have a uniform distribution in the range 0.0 (inclusive)
	 * to 1.0 (exclusive).
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Float value.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public float nextFloat
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextFloat(): skip = " + skip + " illegal");
			}

		// Next random number is in the range -2^63 .. +2^63 - 1.
		// Divide by 2^64 yielding a number in the range -0.5 .. +0.5.
		// Add 0.5 yielding a number in the range 0.0 .. 1.0.
		return (float) (next (skip)) * F_2_POW_NEG_64 + 0.5f;
		}

	/**
	 * Return the double precision floating point value from the next
	 * pseudorandom value in this PRNG's sequence. The returned numbers have a
	 * uniform distribution in the range 0.0 (inclusive) to 1.0 (exclusive).
	 *
	 * @return  Double value.
	 */
	public double nextDouble()
		{
		// Next random number is in the range -2^63 .. +2^63 - 1.
		// Divide by 2^64 yielding a number in the range -0.5 .. +0.5.
		// Add 0.5 yielding a number in the range 0.0 .. 1.0.
		return (double) (next()) * D_2_POW_NEG_64 + 0.5;
		}

	/**
	 * Return the double precision floating point value from the pseudorandom
	 * value the given number of positions ahead in this PRNG's sequence. The
	 * returned numbers have a uniform distribution in the range 0.0 (inclusive)
	 * to 1.0 (exclusive).
	 *
	 * @param  skip  Number of positions to skip.
	 *
	 * @return  Double value.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>skip</TT> &lt;= 0.
	 */
	public double nextDouble
		(long skip)
		{
		if (skip <= 0)
			{
			throw new IllegalArgumentException
				("Random.nextDouble(): skip = " + skip + " illegal");
			}

		// Next random number is in the range -2^63 .. +2^63 - 1.
		// Divide by 2^64 yielding a number in the range -0.5 .. +0.5.
		// Add 0.5 yielding a number in the range 0.0 .. 1.0.
		return (double) (next (skip)) * D_2_POW_NEG_64 + 0.5;
		}

// Hidden operations.

	/**
	 * Return the next 64-bit pseudorandom value in this PRNG's sequence.
	 *
	 * @return  Pseudorandom value.
	 */
	protected abstract long next();

	/**
	 * Return the 64-bit pseudorandom value the given number of positions ahead
	 * in this PRNG's sequence.
	 *
	 * @param  skip  Number of positions to skip, assumed to be &gt; 0.
	 *
	 * @return  Pseudorandom value.
	 */
	protected abstract long next
		(long skip);

	}
