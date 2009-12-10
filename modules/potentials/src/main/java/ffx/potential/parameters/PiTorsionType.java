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
 * The PiTorsionType class defines a Pi-Torsion energy term.
 */
public final class PiTorsionType extends BaseType {
	/**
	 * Atom classes that form this Pi-Torsion.
	 */
	public final int atomClasses[];
	/**
	 * Force constant.
	 */
	public final double forceConstant;

	/**
	 * PiTorsionType Constructor.
	 * 
	 * @param atomClasses
	 *            int[]
	 * @param forceConstant
	 *            double
	 */
	public PiTorsionType(int atomClasses[], double forceConstant) {
		super(ForceField.ForceFieldType.PITORS, sortKey(atomClasses));
		this.atomClasses = atomClasses;
		this.forceConstant = forceConstant;
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
	 * Nicely formatted Pi-Torsion type.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		return String.format("pitors  %5d  %5d  %4.2f", atomClasses[0],
				atomClasses[1], forceConstant);
	}

	/**
	 * Convert Pi-Torsion energy to kcal/mole.
	 */
	public static double units = 1.0;
}
