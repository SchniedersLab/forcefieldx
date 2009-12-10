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

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.toRadians;

/**
 * The TorsionType class defines a torsional angle.
 */
public final class TorsionType extends BaseType {
	/**
	 * Atom classes that for this Torsion angle.
	 */
	public final int atomClasses[];
	/**
	 * Number of terms in the Fourier series.
	 */
	public final int terms;
	/**
	 * Amplitudes of the Fourier series.
	 */
	public final double amplitude[];
	/**
	 * Phases of the Fourier series in degrees.
	 */
	public final double phase[];
	/**
	 * Cosine of the phase angle.
	 */
	public final double cosine[];
	/**
	 * Sine of the phase angle.
	 */
	public final double sine[];
	/**
	 * Periodicity of the Fourier series.
	 */
	public final int periodicity[];

	/**
	 * TorsionType Constructor.
	 * 
	 * @param atomClasses
	 *            int[]
	 * @param amplitude
	 *            double[]
	 * @param phase
	 *            double[]
	 * @param periodicity
	 *            double[]
	 */
	public TorsionType(int atomClasses[], double amplitude[], double phase[],
			int periodicity[]) {
		super(ForceField.ForceFieldType.TORSION, sortKey(atomClasses));
		this.atomClasses = atomClasses;
		int max = 1;
		for (int i = 0; i < periodicity.length; i++) {
			if (periodicity[i] > max) {
				max = periodicity[i];
			}
		}
		terms = max;
		if (periodicity.length != max) {
			this.amplitude = new double[max];
			this.phase = new double[max];
			this.periodicity = new int[max];
			for (int i = 0; i < periodicity.length; i++) {
				this.amplitude[periodicity[i] - 1] = amplitude[i];
				this.phase[periodicity[i] - 1] = phase[i];
				this.periodicity[periodicity[i] - 1] = periodicity[i];
			}
		} else {
			this.amplitude = amplitude;
			this.phase = phase;
			this.periodicity = periodicity;
		}
		cosine = new double[terms];
		sine = new double[terms];
		for (int i = 0; i < terms; i++) {
			double angle = toRadians(this.phase[i]);
			cosine[i] = cos(angle);
			sine[i] = sin(angle);
		}
	}

	/**
	 * Nicely formatted Torsion angle.
	 * 
	 * @return String
	 */
	@Override
	public String toString() {
		StringBuffer torsionBuffer = new StringBuffer("torsion");
		for (int i : atomClasses) {
			torsionBuffer.append(String.format("  %5d", i));
		}
		for (int i = 0; i < amplitude.length; i++) {
			torsionBuffer.append(String.format("  % 5.3f  %5.3f  %1d",
					amplitude[i], phase[i], periodicity[i]));
		}
		return torsionBuffer.toString();
	}

	/**
	 * This method sorts the atom classes for the torsion.
	 * 
	 * @param c
	 *            atomClasses
	 * @return lookup key
	 */
	public static String sortKey(int c[]) {
		if (c == null || c.length != 4) {
			return null;
		}
		if (c[1] < c[2]) {
			// Do nothing.
		} else if (c[2] < c[1]) {
			// Reverse the order.
			int temp = c[0];
			c[0] = c[3];
			c[3] = temp;
			temp = c[1];
			c[1] = c[2];
			c[2] = temp;
		} else if (c[0] <= c[3]) {
			// Do nothing.
		} else {
			// Reverse the order.
			int temp = c[0];
			c[0] = c[3];
			c[3] = temp;
			temp = c[1];
			c[1] = c[2];
			c[2] = temp;
		}
		String key = c[0] + " " + c[1] + " " + c[2] + " " + c[3];
		return key;
	}

	/**
	 * Convert Torsional Angle energy to kcal/mole.
	 */
	public static final double units = 0.5;
}
