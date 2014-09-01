/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.xray;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>
 * SolventBinaryFormFactor class.</p>
 *
 * @author fenn
 *
 */
public final class SolventBinaryFormFactor implements FormFactor {

    private final Atom atom;
    private double xyz[] = new double[3];
    private double dxyz[] = new double[3];
    private double proberad;

    /**
     * <p>
     * Constructor for SolventBinaryFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param proberad a double.
     */
    public SolventBinaryFormFactor(Atom atom, double proberad) {
        this(atom, proberad, atom.getXYZ());
    }

    /**
     * <p>
     * Constructor for SolventBinaryFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param proberad a double.
     * @param xyz an array of double.
     */
    public SolventBinaryFormFactor(Atom atom, double proberad, double xyz[]) {
        this.atom = atom;
        this.proberad = proberad;

        update(xyz);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double[] xyz) {
        VectorMath.diff(this.xyz, xyz, dxyz);
        return rho(f, lambda, VectorMath.r(dxyz));
    }

    /**
     * <p>
     * rho</p>
     *
     * @param f a double.
     * @param lambda a double.
     * @param ri a double.
     * @return a double.
     */
    public double rho(double f, double lambda, double ri) {
        if (ri <= proberad) {
            return 0.0;
        } else {
            return f * 1.0;
        }
    }

    // no derivative for the binary model!!!
    /**
     * {@inheritDoc}
     */
    @Override
    public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementmode) {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double xyz[]) {
        update(xyz, 0.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update(double xyz[], double badd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
    }
}
