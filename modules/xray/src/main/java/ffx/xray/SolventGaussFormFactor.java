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

import static org.apache.commons.math3.util.FastMath.exp;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>
 * SolventGaussFormFactor class.</p>
 *
 * @author fenn
 *
 */
public final class SolventGaussFormFactor implements FormFactor {

    private final Atom atom;
    private final double xyz[] = new double[3];
    private final double dxyz[] = new double[3];
    private final double g[] = new double[3];
    private final double isd2;

    /**
     * <p>
     * Constructor for SolventGaussFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param sd a double.
     */
    public SolventGaussFormFactor(Atom atom, double sd) {
        this(atom, sd, atom.getXYZ());
    }

    /**
     * <p>
     * Constructor for SolventGaussFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param sd a double.
     * @param xyz an array of double.
     */
    public SolventGaussFormFactor(Atom atom, double sd, double xyz[]) {
        this.atom = atom;
        isd2 = 1.0 / (sd * sd);
        update(xyz);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double xyz[]) {
        VectorMath.diff(this.xyz, xyz, dxyz);
        return rho(f, lambda, VectorMath.rsq(dxyz));
    }

    /**
     * <p>
     * rho</p>
     *
     * @param f a double.
     * @param lambda a double.
     * @param rsq a double.
     * @return a double.
     */
    public double rho(double f, double lambda, double rsq) {
        return f + exp(-rsq * isd2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementmode) {
        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES) {
            return;
        }
        VectorMath.diff(this.xyz, xyz, dxyz);
        double r2 = VectorMath.rsq(dxyz);
        double rho = exp(-r2 * isd2);
        double prefactor = -dfc * 2.0 * rho * isd2;
        g[0] = prefactor * dxyz[0];
        g[1] = prefactor * dxyz[1];
        g[2] = prefactor * dxyz[2];
        atom.addToXYZGradient(g[0], g[1], g[2]);
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
