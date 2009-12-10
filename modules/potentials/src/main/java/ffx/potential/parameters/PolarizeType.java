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

import static java.lang.Math.pow;

/**
 * The PolarizeType class defines an isotropic atomic polarizability.
 */
public final class PolarizeType extends BaseType {
	private static final double sixth = 1.0 / 6.0;
	/**
	 * Atom type number.
	 */
	public final int atomType;
	/**
	 * Thole damping factor.
	 */
	public final double thole;
	/**
	 * Value of polarizability scale factor.
	 */
	public final double pdamp;
	/**
	 * Isotropic polarizability in units of Angstroms^3.
	 */
	public final double polarizability;
	/**
	 * Connected types in the polarization group of each atom. (may be null)
	 */
	public int[] polarizationGroup;

	/**
	 * PolarizeType Constructor.
	 * 
	 * @param atomType
	 *            int
	 * @param polarizability
	 *            double
	 * @param polarizationGroup
	 *            int[]
	 */
	public PolarizeType(int atomType, double polarizability, double thole,
			int polarizationGroup[]) {
		super(ForceField.ForceFieldType.POLARIZE, new String("" + atomType));
		this.atomType = atomType;
		this.thole = thole;
		this.polarizability = polarizability;
		this.polarizationGroup = polarizationGroup;
		if (thole == 0.0) {
			pdamp = 0.0;
		} else {
			pdamp = pow(polarizability, sixth);
		}
	}

        public void add(int key) {
            for (int i : polarizationGroup) {
                if (key == i) {
                    return;
                }
            }
            int len = polarizationGroup.length;
            int newGroup[] = new int[len+1];
            for (int i=0; i<len; i++) {
                newGroup[i] = polarizationGroup[i];
            }
            newGroup[len] = key;
            polarizationGroup = newGroup;
        }

	/**
	 * Nicely formatted polarization type.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		StringBuffer polarizeString = new StringBuffer(String.format(
				"polarize  %5d  %5.3f %5.3f", atomType, polarizability, thole));
		if (polarizationGroup != null) {
			for (int a : polarizationGroup) {
				polarizeString.append(String.format("  %5d", a));
			}
		}
		return polarizeString.toString();
	}
}
