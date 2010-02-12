/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
 * The ChargeType class defines a partial atomic charge type.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class ChargeType extends BaseType {
	/**
	 * The atom type that uses this charge parameter.
	 */
	public final int atomType;
	/**
	 * Partial atomic charge in units of electrons.
	 */
	public final double charge;

	/**
	 * ChargeType constructor.
	 * 
	 * @param atomType
	 *            int
	 * @param charge
	 *            double
	 */
	public ChargeType(int atomType, double charge) {
		super(ForceField.ForceFieldType.CHARGE, new String("" + atomType));
		this.atomType = atomType;
		this.charge = charge;
	}

	/**
	 * Nicely formatted Charge type.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		return String.format("charge  %5d  % 7.5f", atomType, charge);
	}
}
