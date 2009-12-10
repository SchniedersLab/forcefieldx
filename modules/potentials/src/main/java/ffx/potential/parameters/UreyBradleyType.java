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
 * The UreyBradleyType class defines one harmonic UreyBradley cross term.
 */
public final class UreyBradleyType extends BaseType {
	/**
	 * Atom classes that form this Urey-Bradley cross term.
	 */
	public final int atomClasses[];
	/**
	 * Force constant (Kcal/mole/angstroms^2).
	 */
	public final double forceConstant;
	/**
	 * Equilibrium 1-3 separation (Angstroms).
	 */
	public final double distance;

	/**
	 * UreyBradleyType constructor.
	 * 
	 * @param atomClasses
	 *            int[]
	 * @param forceConstant
	 *            double
	 * @param distance
	 *            double
	 */
	public UreyBradleyType(int atomClasses[], double forceConstant,
			double distance) {
		super(ForceField.ForceFieldType.UREYBRAD, sortKey(atomClasses));
		this.atomClasses = atomClasses;
		this.forceConstant = forceConstant;
		this.distance = distance;
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
	 * Nicely formatted Urey-Bradley string.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		return String.format("ureybrad  %5d  %5d  %5d  %6.2f  %7.4f",
				atomClasses[0], atomClasses[1], atomClasses[2], forceConstant,
				distance);
	}

	/**
	 * Convert bond stretch energy to kcal/mole.
	 */
	public static final double units = 1.0;
	/**
	 * Cubic coefficient in bond stretch potential.
	 */
	public static final double cubic = 0.0;
	/**
	 * Quartic coefficient in bond stretch potential.
	 */
	public static final double quartic = 0.0;
}
