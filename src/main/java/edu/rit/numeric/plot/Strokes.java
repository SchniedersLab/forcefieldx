//******************************************************************************
//
// File:    Strokes.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.Strokes
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

package edu.rit.numeric.plot;

import java.awt.BasicStroke;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Strokes provides operations for creating several common patterns of
 * drawing strokes.
 *
 * @author  Alan Kaminsky
 * @version 11-Dec-2007
 */
public class Strokes
	{

// Prevent construction.

	private Strokes()
		{
		}

// Exported operations.

	/**
	 * Returns a solid stroke of the given width with square corners and ends.
	 *
	 * @param  theWidth  Stroke width.
	 */
	public static BasicStroke solid
		(double theWidth)
		{
		return new BasicStroke
			((float) theWidth,
			 BasicStroke.CAP_SQUARE,
			 BasicStroke.JOIN_MITER,
			 10.0f);
		}

	/**
	 * Returns a dotted stroke of the given width with square dots.
	 *
	 * @param  theWidth  Stroke width.
	 */
	public static BasicStroke dotted
		(double theWidth)
		{
		return dashed (theWidth, 1, 1);
		}

	/**
	 * Returns a dashed stroke of the given width with square dashes. Each
	 * dash's length is three times the stroke width. The length of the gap
	 * between dashes equals the stroke width.
	 *
	 * @param  theWidth  Stroke width.
	 */
	public static BasicStroke dashed
		(double theWidth)
		{
		return dashed (theWidth, 3, 1);
		}

	/**
	 * Returns a dashed stroke of the given width with square dashes. Each
	 * dash's length is d times the stroke width. The length of the gap between
	 * dashes equals g times the stroke width.
	 *
	 * @param  theWidth  Stroke width.
	 * @param  d         Dash length in terms of stroke width.
	 * @param  g         Gap length in terms of stroke width.
	 */
	public static BasicStroke dashed
		(double theWidth,
		 double d,
		 double g)
		{
		return new BasicStroke
			((float) theWidth,
			 BasicStroke.CAP_SQUARE,
			 BasicStroke.JOIN_MITER,
			 10.0f,
			 new float[] {(float)((d-1)*theWidth), (float)((g+1)*theWidth)},
			 0);
		}

	/**
	 * Returns a solid stroke of the given width with round corners and ends.
	 *
	 * @param  theWidth  Stroke width.
	 */
	public static BasicStroke roundSolid
		(double theWidth)
		{
		return new BasicStroke
			((float) theWidth,
			 BasicStroke.CAP_ROUND,
			 BasicStroke.JOIN_ROUND,
			 10.0f);
		}

	/**
	 * Returns a dotted stroke of the given width with round dots.
	 *
	 * @param  theWidth  Stroke width.
	 */
	public static BasicStroke roundDotted
		(double theWidth)
		{
		return roundDashed (theWidth, 1, 1);
		}

	/**
	 * Returns a dashed stroke of the given width with rounded dashes. Each
	 * dash's length is three times the stroke width. The length of the gap
	 * between dashes equals the stroke width.
	 *
	 * @param  theWidth  Stroke width in points.
	 */
	public static BasicStroke roundDashed
		(double theWidth)
		{
		return roundDashed (theWidth, 3, 1);
		}

	/**
	 * Returns a dashed stroke of the given width with rounded dashes. Each
	 * dash's length is d times the stroke width. The length of the gap between
	 * dashes equals g times the stroke width.
	 *
	 * @param  theWidth  Stroke width in points.
	 * @param  d         Dash length in terms of stroke width.
	 * @param  g         Gap length in terms of stroke width.
	 */
	public static BasicStroke roundDashed
		(double theWidth,
		 double d,
		 double g)
		{
		return new BasicStroke
			((float) theWidth,
			 BasicStroke.CAP_ROUND,
			 BasicStroke.JOIN_ROUND,
			 10.0f,
			 new float[] {(float)((d-1)*theWidth), (float)((g+1)*theWidth)},
			 0);
		}

	/**
	 * Write the given BasicStroke object to the given object output stream.
	 *
	 * @param  theStroke  Stroke.
	 * @param  out        Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void writeExternal
		(BasicStroke theStroke,
		 ObjectOutput out)
		throws IOException
		{
		out.writeObject
			(theStroke == null ? null : new StrokeWrapper (theStroke));
		}

	/**
	 * Read a BasicStroke object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @return  Stroke.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if a class needed to deserialize the stroke object could not
	 *     be found.
	 */
	public static BasicStroke readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		StrokeWrapper wrapper = (StrokeWrapper) in.readObject();
		return wrapper == null ? null : wrapper.getStroke();
		}

	/**
	 * Class StrokeWrapper is used to serialize a BasicStroke.
	 *
	 * @author  Alan Kaminsky
	 * @version 11-Dec-2007
	 */
	private static class StrokeWrapper
		implements Externalizable
		{
		private static final long serialVersionUID = 6899737592162456726L;

		private float width;
		private int cap;
		private int join;
		private float miterlimit;
		private float[] dash;
		private float dash_phase;

		public StrokeWrapper()
			{
			}

		public StrokeWrapper
			(BasicStroke theStroke)
			{
			this.width = theStroke.getLineWidth();
			this.cap = theStroke.getEndCap();
			this.join = theStroke.getLineJoin();
			this.miterlimit = theStroke.getMiterLimit();
			this.dash = theStroke.getDashArray();
			this.dash_phase = theStroke.getDashPhase();
			}

		public BasicStroke getStroke()
			{
			return new BasicStroke
				(width, cap, join, miterlimit, dash, dash_phase);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeFloat (width);
			out.writeInt (cap);
			out.writeInt (join);
			out.writeFloat (miterlimit);
			int n = dash == null ? 0 : dash.length;
			out.writeInt (n);
			for (int i = 0; i < n; ++ i)
				{
				out.writeFloat (dash[i]);
				}
			out.writeFloat (dash_phase);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException
			{
			width = in.readFloat();
			cap = in.readInt();
			join = in.readInt();
			miterlimit = in.readFloat();
			int n = in.readInt();
			dash = n == 0 ? null : new float [n];
			for (int i = 0; i < n; ++ i)
				{
				dash[i] = in.readFloat();
				}
			dash_phase = in.readFloat();
			}
		}

	}
