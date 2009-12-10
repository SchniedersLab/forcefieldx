/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.parameters;

import static java.lang.Math.PI;
import static java.lang.Math.pow;

/**
 * The AngleType class defines one harmonic angle bend energy term.
 */
public final class AngleType extends BaseType {
	/**
	 * Atom classes that for this Angle type.
	 */
	public final int atomClasses[];
	/**
	 * Force constant (Kcal/mole/radian^2).
	 */
	public final double forceConstant;
	/**
	 * Equilibrium angle (degrees). There can be up to three equilibrium angles,
	 * depending on the number of attached hydrogens (0, 1, or 2).
	 */
	public final double angle[];

	/**
	 * @param atomClasses
	 * @param forceConstant
	 * @param angle
	 */
	public AngleType(int atomClasses[], double forceConstant, double angle[]) {
		super(ForceField.ForceFieldType.ANGLE, sortKey(atomClasses));
		this.atomClasses = atomClasses;
		this.forceConstant = forceConstant;
		this.angle = angle;
	}

	/**
	 * This method sorts the atom classes as: min, c[1], max
	 * 
	 * @param c
	 *            atomClasses
	 * @return lookup key
	 */
	public static String sortKey(int c[]) {
		if (c == null || c.length != 3) {
			return null;
		}
		if (c[0] > c[2]) {
			int temp = c[0];
			c[0] = c[2];
			c[2] = temp;
		}
		String key = c[0] + " " + c[1] + " " + c[2];
		return key;
	}

	/**
	 * Nicely formatted Angle bending string.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		StringBuffer angleString = new StringBuffer(String.format(
				"angle  %5d  %5d  %5d  %6.2f", atomClasses[0], atomClasses[1],
				atomClasses[2], forceConstant));
		for (double eq : angle) {
			angleString.append(String.format("  %6.2f", eq));
		}
		return angleString.toString();
	}

	/**
	 * Cubic coefficient in angle bending potential.
	 */
	public static double cubic = -0.014;
	/**
	 * Quartic coefficient in angle bending potential.
	 */
	public static double quartic = 0.000056;
	/**
	 * Quintic coefficient in angle bending potential.
	 */
	public static double quintic = -0.0000007;
	/**
	 * Sextic coefficient in angle bending potential.
	 */
	public static double sextic = 0.000000022;
	/**
	 * Convert angle bending energy to kcal/mole.
	 */
	public static double units = 1.0 / pow(180.0 / PI, 2);
}
