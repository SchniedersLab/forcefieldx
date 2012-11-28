//******************************************************************************
//
// File:    Mcg1Random.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Mcg1Random
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

/**
 * Class Mcg1Random provides a default pseudorandom number generator (PRNG)
 * designed for use in parallel scientific programming. To create an instance of
 * class Mcg1Random, either use the <TT>Mcg1Random()</TT> constructor, or use
 * the static <TT>getInstance(long,String)</TT> method in class {@linkplain
 * Random}.
 * <P>
 * Class Mcg1Random uses L'Ecuyer's 63-bit multiplicative congruential
 * generator:
 * <PRE>
 *     seed := seed * A (mod M);
 * </PRE>
 * with <I>A</I> = 2307085864 and <I>M</I> = 2<SUP>63</SUP>-25. For further
 * information, see P. L'Ecuyer, F. Blouin, and R. Couture, A search for good
 * multiple recursive random number generators, <I>ACM Transactions on Modeling
 * and Computer Simulation,</I> 3(2):87-98, April 1993.
 *
 * @author  Alan Kaminsky
 * @version 01-Mar-2008
 */
public class Mcg1Random
	extends Random
	{

// Hidden data members.

	// Multiplicative congruential generator parameters.
	private static final long A = 2307085864L;
	private static final long M = 9223372036854775783L;

	// Table of powers of A (mod M).
	// powtable[i] = A^(2^i) (mod M), i = 0, 1, 2, ..., 62.
	private static final long[] powtable = new long[]
		{
		/* 0*/          2307085864L,
		/* 1*/ 5322645183868626496L,
		/* 2*/  983401115462215297L,
		/* 3*/ 3556108090190705823L,
		/* 4*/ 7990665143195102590L,
		/* 5*/ 2110036525984475599L,
		/* 6*/ 7043012601020815633L,
		/* 7*/ 8705155707092105232L,
		/* 8*/ 3648485552813098205L,
		/* 9*/ 3168429798853819517L,
		/*10*/ 7370936612916750461L,
		/*11*/ 7860663018156131952L,
		/*12*/ 3001105880121306407L,
		/*13*/ 2701734581708584636L,
		/*14*/   44173215984149523L,
		/*15*/ 4386281867185367357L,
		/*16*/ 6179163218358095360L,
		/*17*/ 7483044026478026567L,
		/*18*/ 3475714592143337300L,
		/*19*/ 1764426730688581302L,
		/*20*/ 3750657437672096664L,
		/*21*/  622726075290379426L,
		/*22*/ 5708473958970181660L,
		/*23*/ 4021546582722653103L,
		/*24*/ 2336213934427760687L,
		/*25*/ 1250271094601288883L,
		/*26*/ 3574383011208782094L,
		/*27*/ 8396902035548884488L,
		/*28*/ 8461483610275050157L,
		/*29*/ 4570169555765982077L,
		/*30*/ 8905831846701231221L,
		/*31*/ 8735916407118983196L,
		/*32*/ 2440495732904503112L,
		/*33*/ 1885457269016286005L,
		/*34*/ 4972446378304258072L,
		/*35*/ 5086882142287647560L,
		/*36*/ 7606891628733932672L,
		/*37*/ 1492990033908793408L,
		/*38*/ 9099993837175275499L,
		/*39*/  164616137930049276L,
		/*40*/ 5117944347055477320L,
		/*41*/ 3732738446422589684L,
		/*42*/  577797231373159603L,
		/*43*/ 2884327325873197522L,
		/*44*/ 4833803989390835826L,
		/*45*/ 7647846260763424785L,
		/*46*/ 4871120313232679781L,
		/*47*/ 2522743552130321382L,
		/*48*/ 2285147082121189109L,
		/*49*/ 3702619298913044713L,
		/*50*/ 7517285182136659617L,
		/*51*/ 1501022168611987834L,
		/*52*/ 4083684657803873370L,
		/*53*/ 1174110446001111617L,
		/*54*/   82581059520186299L,
		/*55*/ 1334190853588951475L,
		/*56*/ 3130709730706025384L,
		/*57*/ 8886205968707213290L,
		/*58*/  993283284549990895L,
		/*59*/ 3258516944203296282L,
		/*60*/ 4273233140749644635L,
		/*61*/ 7682756089153477585L,
		/*62*/ 8243539608199123644L,
		};

	// Seed for this PRNG.
	private long seed;

	// 128 bytes of extra padding to avert cache interference.
	private transient long p0, p1, p2, p3, p4, p5, p6, p7;
	private transient long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new PRNG with the given seed. The seed must not be 0.
	 *
	 * @param  seed  Seed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>seed</TT> = 0.
	 */
	public Mcg1Random
		(long seed)
		{
		setSeed (seed);
		}

// Exported operations.

	/**
	 * Set this PRNG's seed. The seed must not be 0.
	 *
	 * @param  seed  Seed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>seed</TT> = 0.
	 */
	public void setSeed
		(long seed)
		{
		if (seed == 0L)
			{
			throw new IllegalArgumentException
				("Mcg1Random.setSeed(): seed = 0 illegal");
			}

		// Make sure seed is nonnegative.
		this.seed = seed & 0x7FFFFFFFFFFFFFFFL;
		}

// Hidden operations.

	/**
	 * Return the next 64-bit pseudorandom value in this PRNG's sequence.
	 *
	 * @return  Pseudorandom value.
	 */
	protected long next()
		{
		// Multiply seed (a 64-bit number) by A (a 32-bit number) yielding x (a
		// 96-bit number). Bits 63-0 of x are stored in x_63_0. Bits 95-32 of x
		// are stored in x_95_32. Note that these overlap.
		long tmp_63_0 = (seed & 0x00000000FFFFFFFFL) * A;
		long tmp_95_32 = (seed >>> 32) * A;
		long x_63_0 = (tmp_95_32 << 32) + tmp_63_0;
		long x_95_32 = tmp_95_32 + (tmp_63_0 >>> 32);

		// Compute x mod M, where M = 2^63-25. For the algorithm, see the
		// Handbook of Applied Cryptography, Section 14.3.4.

		// q = int (x / 2^63)
		long q = x_95_32 >>> 31;

		// r = x mod 2^63
		long r = x_63_0 & 0x7FFFFFFFFFFFFFFFL;

		// r = r + (q * 25) mod 2^63
		r += q * 25L;

		// If there was a carry into the high-order bit of r, or if r >= M,
		// subtract M.
		if (r < 0L || r >= M) r -= M;

		// r = x mod M becomes the new seed.
		seed = r;

		// Since seed is only in the range 0 .. 2^63 - 1, left-shift yielding a
		// number in the range -2^63 .. 2^63 - 1.
		return seed << 1;
		}

	/**
	 * Return the 64-bit pseudorandom value the given number of positions ahead
	 * in this PRNG's sequence.
	 *
	 * @param  skip  Number of positions to skip, assumed to be &gt; 0.
	 *
	 * @return  Pseudorandom value.
	 */
	protected long next
		(long skip)
		{
		// Compute seed * A^skip (mod M).
		int i = 0;
		while (skip != 0L)
			{
			if ((skip & 1L) != 0L) seed = modMultiply (powtable[i], seed);
			skip >>>= 1;
			++ i;
			}

		// Since seed is only in the range 0 .. 2^63 - 1, left-shift yielding a
		// number in the range -2^63 .. 2^63 - 1.
		return seed << 1;
		}

	/**
	 * Returns a * b (mod M). a and b are assumed to be in the range 0 ..
	 * 2^63-1.
	 *
	 * @param  a  First number to multiply.
	 * @param  b  Second number to multiply.
	 *
	 * @return  a * b (mod M).
	 */
	private static long modMultiply
		(long a,
		 long b)
		{
		// Let a = s*2^32 + t, b = u*2^32 + v, where s, t, u, and v are 32-bit
		// numbers. Form the four 64-bit products tv, sv, tu, and su. Add these
		// up in the appropriate combinations to get x = a * b, where x =
		// x_127_64*2^64 + x_63_0:
		//                     +--------+--------+
		//                     |    s   |    t   | = a
		//                     +--------+--------+
		//                     +--------+--------+
		//                   * |    u   |    v   | = b
		//                     +--------+--------+
		// ----------------------------------------
		//                     +-----------------+
		//                     |        tv       |
		//                     +-----------------+
		//            +-----------------+
		//            |        sv       |
		//            +-----------------+
		//            +-----------------+
		//            |        tu       |
		//            +-----------------+
		//   +-----------------+
		// + |        su       |
		//   +-----------------+
		// ----------------------------------------
		//   +-----------------+-----------------+
		//   |    x_127_64     |     x_63_0      | = x = a * b
		//   +-----------------+-----------------+
		long s = a >>> 32;
		long t = a & 0xFFFFFFFFL;
		long u = b >>> 32;
		long v = b & 0xFFFFFFFFL;
		long tv = t * v;
		long sv = s * v;
		long tu = t * u;
		long su = s * u;
		long tmp = (tv >>> 32) + (sv & 0xFFFFFFFFL) + (tu & 0xFFFFFFFFL);
		long x_63_0 = (tv & 0xFFFFFFFFL) + (tmp << 32);
		long x_127_64 = (tmp >>> 32) + (sv >>> 32) + (tu >>> 32) + su;

		// Compute x mod M, where M = 2^63-25. For the algorithm, see the
		// Handbook of Applied Cryptography, Section 14.3.4.

		// q = int (x / 2^63)
		long q = (x_127_64 << 1) | (x_63_0 >>> 63);

		// r = x mod 2^63
		long r = x_63_0 & 0x7FFFFFFFFFFFFFFFL;

		while (q > 0L)
			{
			// qc = q * 25
			// Multiply q (a 64-bit number) by c (a 32-bit number) yielding qc
			// (a 96-bit number). Bits 63-0 of qc are stored in qc_63_0. Bits
			// 95-32 of qc are stored in qc_95_32. Note that these overlap.
			long tmp_63_0 = (q & 0xFFFFFFFFL) * 25L;
			long tmp_95_32 = (q >>> 32) * 25L;
			long qc_63_0 = (tmp_95_32 << 32) + tmp_63_0;
			long qc_95_32 = tmp_95_32 + (tmp_63_0 >>> 32);

			// q = int (qc / 2^63)
			q = qc_95_32 >>> 31;

			// r = r + (qc mod 2^63)
			r += qc_63_0 & 0x7FFFFFFFFFFFFFFFL;

			// If there was a carry into the high-order bit of r, or if r >= M,
			// subtract M.
			if (r < 0L || r >= M) r -= M;
			}

		// Return r = x mod M.
		return r;
		}

	}
