//******************************************************************************
//
// File:    Test01.java
// Package: edu.rit.draw.item.test
// Unit:    Class edu.rit.draw.item.test.Test01
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

package edu.rit.draw.item.test;

import edu.rit.draw.Drawing;

import edu.rit.draw.item.*;

/**
 * Class Test01 is a unit test program for class {@linkplain Drawing}. It
 * creates a drawing and stores it in the file <TT>"Test01.drw"</TT>. To see the
 * drawing, type these commands:
 * <P>
 * <TT>$ java edu.rit.draw.item.test.Test01</TT>
 * <BR><TT>$ java edu.rit.draw.View Test01.drw</TT>
 *
 * @author  Alan Kaminsky
 * @version 12-Jul-2006
 */
public class Test01
	{
	/**
	 * Prevent construction.
	 */
	private Test01()
		{
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		double W = 72.0;
		double H = 36.0;
		Size S = new Size (W, H);

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Rectangle.defaultSize (S);
		Rectangle r1 = new Rectangle().add();
		Rectangle r2 = new Rectangle().w(r1.e().e(W)).add();

		OutlinedItem.defaultOutline (DashedOutline.NORMAL_OUTLINE);
		Line.defaultEndArrow (Arrow.SOLID);
		new Line().to(r1.e()).to(r2.w()).addFirst();

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Rectangle.defaultSize (S);
		Rectangle r3 = new Rectangle().n(r1.s()).add();
		Rectangle r4 = new Rectangle().w(r3.e().e(W)).add();

		OutlinedItem.defaultOutline (SquareDottedOutline.NORMAL_OUTLINE);
		Line.defaultEndArrow (Arrow.SOLID);
		new Line().to(r3.e()).to(r4.w()).addFirst();

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Rectangle.defaultSize (S);
		Rectangle r5 = new Rectangle().n(r3.s()).add();
		Rectangle r6 = new Rectangle().w(r5.e().e(W)).add();

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Line.defaultEndArrow (Arrow.SOLID);
		new Line().to(r5.e()).to(r6.w()).addFirst();

		OutlinedItem.defaultOutline (new DottedOutline().width(3));
		Line.defaultEndArrow (Arrow.NONE);
		new Line().to(r5.c()).to(r1.c()).to(r6.c()).to(r2.c()).add();

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Line.defaultStartArrow (Arrow.SOLID);
		Line.defaultRound (18);
		new Line().to(0,144).hby(36).by(-36,36).hby(36).add();
		new Line().to(72,144).hby(36).by(0,36).hby(36).add();
		new Line().to(144,144).hby(36).by(36,36).hby(36).add();

		OutlinedItem.defaultOutline (SolidOutline.NORMAL_OUTLINE);
		Line.defaultStartArrow (Arrow.NONE);
		Line.defaultRound (4);
		new Line().to(r5.sw().s(4)).vby(4).hto(r5.s()).vby(4).add();
		new Line().to(r5.se().s(4)).vby(4).hto(r5.s()).vby(4).add();

		Drawing.write ("Test01.drw");
		}
	}
