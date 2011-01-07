/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

package ffx.xray;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 *
 * @author fenn
 */
public final class SolventBinaryFormFactor implements FormFactor {
    private final Atom atom;
    private double xyz[] = new double[3];
    private double dxyz[] = new double[3];
    private double proberad;

    public SolventBinaryFormFactor(Atom atom, double proberad) {
        this(atom, proberad, atom.getXYZ());
    }

    public SolventBinaryFormFactor(Atom atom, double proberad, double xyz[]) {
        this.atom = atom;
        this.proberad = proberad;

        update(xyz);
    }

    @Override
    public double rho(double f, double[] xyz) {
        VectorMath.diff(this.xyz, xyz, dxyz);
        return rho(f, VectorMath.r(dxyz));
    }

    public double rho(double f, double ri) {
        if (ri <= proberad) {
            return 0.0;
        } else {
            return f * 1.0;
        }
    }

    // no derivative for the binary model!!!
    @Override
    public void rho_grad(double[] xyz, double dfc, RefinementMode refinementmode) {
    }
    
    @Override
    public void update(double xyz[]) {
        update(xyz, 0.0);
    }

    @Override
    public void update(double xyz[], double badd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
    }
}
