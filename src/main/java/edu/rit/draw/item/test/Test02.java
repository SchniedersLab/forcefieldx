//******************************************************************************
//
// File:    Test02.java
// Package: edu.rit.draw.item.test
// Unit:    Class edu.rit.draw.item.test.Test02
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

import java.awt.Font;

/**
 * Class Test02 is a unit test program for class {@linkplain Drawing}. It
 * creates a drawing and stores it in the file <TT>"Test02.drw"</TT>. To see the
 * drawing, type these commands:
 * <P>
 * <TT>$ java edu.rit.draw.item.test.Test02</TT>
 * <BR><TT>$ java edu.rit.draw.View Test02.drw</TT>
 *
 * @author  Alan Kaminsky
 * @version 12-Jul-2006
 */
public class Test02
	{
	/**
	 * Prevent construction.
	 */
	private Test02()
		{
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		Outline THICK = new SolidOutline().width(3);
		Outline THIN = new SolidOutline().width(0.5f);
		Text.defaultFont (new Font ("Times New Roman", Font.ITALIC, 12));

		new Line().to(0,0).hby(72).vby(72).round(36).outline(THICK).add();

		new Line().to(36,-22).vby(18).outline(THIN).add();
		new Line().to(72,-22).vby(42).outline(THIN).add();
		new Line().to(36,-13).hby(36).outline(THIN).startArrow(Arrow.SOLID).endArrow(Arrow.SOLID).add();
		new Text().text("d").s(54,-13).add();

		new Line().to(94,36).hby(-18).outline(THIN).add();
		new Line().to(94, 0).hby(-42).outline(THIN).add();
		new Line().to(85, 0).vby(36).outline(THIN).startArrow(Arrow.SOLID).endArrow(Arrow.SOLID).add();
		new Text().text("d").w(88,18).add();

		Drawing.write ("Test02.drw");
		}
	}
