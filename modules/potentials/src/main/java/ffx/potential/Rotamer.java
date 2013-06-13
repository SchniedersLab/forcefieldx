/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential;

import ffx.potential.ResidueEnumerations.AminoAcid3;

import static java.lang.Math.max;

/**
 * The Rotamer Class represents one immutable amino acid Rotamer.
 *
 * @author Ava M. Lynn
 */
public class Rotamer {

    public final double chi1;
    public final double chi2;
    public final double chi3;
    public final double chi4;
    public final double angles[];
    public final double sigmas[];
    public final AminoAcid3 name;
    public final int length;

    public Rotamer(AminoAcid3 name, double... values) {
        length = values.length/2;
        angles = new double[max(length, 4)];
        sigmas = new double[max(length, 4)];
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        this.name = name;
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(name.toString());
        int length = angles.length;
        for (int i = 0; i < length; i++) {
            sb.append(String.format(" %6.1f %4.1f", angles[i], sigmas[i]));
        }
        return sb.toString();
    }
}
