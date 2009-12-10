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

/**
 * The BondType class defines one harmonic bond stretch energy term.
 */
public final class BondType extends BaseType {
	/**
	 * Atom classes that form this bond stretch.
	 */
	public final int atomClasses[];
	/**
	 * Force constant (Kcal/mol).
	 */
	public final double forceConstant;
	/**
	 * Equilibrium separation (Angstroms).
	 */
	public final double distance;

	/**
	 * BondType constructor.
	 * 
	 * @param atomClasses
	 *            int[]
	 * @param forceConstant
	 *            double
	 * @param distance
	 *            double
	 */
	public BondType(int atomClasses[], double forceConstant, double distance) {
		super(ForceField.ForceFieldType.BOND, sortKey(atomClasses));
		this.atomClasses = atomClasses;
		this.forceConstant = forceConstant;
		this.distance = distance;
	}

	/**
	 * Nicely formatted bond stretch string.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		return String.format("bond  %5d  %5d  %6.1f  %7.4f", atomClasses[0],
				atomClasses[1], forceConstant, distance);
	}

	/**
	 * This method sorts the atom classes as: min, max
	 * 
	 * @param c
	 *            atomClasses
	 * @return lookup key
	 */
	public static String sortKey(int c[]) {
		if (c == null || c.length != 2) {
			return null;
		}
		String key = null;
		int temp;
		if (c[1] <= c[0]) {
			temp = c[1];
			c[1] = c[0];
			c[0] = temp;
		}
		key = c[0] + " " + c[1];
		return key;
	}

	/**
	 * Convert bond stretch energy to kcal/mole.
	 */
	public static double units = 1.0;
	/**
	 * Cubic coefficient in bond stretch potential.
	 */
	public static double cubic = -2.55;
	/**
	 * Quartic coefficient in bond stretch potential.
	 */
	public static double quartic = 3.793125;
}
