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
 * <p>SolventGaussFormFactor class.</p>
 *
 * @author fenn
 *
 */
public final class SolventGaussFormFactor implements FormFactor {

    private final Atom atom;
    private double xyz[] = new double[3];
    private double dxyz[] = new double[3];
    private double g[] = new double[3];
    private double sd;

    /**
     * <p>Constructor for SolventGaussFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param sd a double.
     */
    public SolventGaussFormFactor(Atom atom, double sd) {
        this(atom, sd, atom.getXYZ());
    }

    /**
     * <p>Constructor for SolventGaussFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param sd a double.
     * @param xyz an array of double.
     */
    public SolventGaussFormFactor(Atom atom, double sd, double xyz[]) {
        this.atom = atom;
        this.sd = sd;

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
     * <p>rho</p>
     *
     * @param f a double.
     * @param lambda a double.
     * @param rsq a double.
     * @return a double.
     */
    public double rho(double f, double lambda, double rsq) {
        double sd2 = sd * sd;
        return f + Math.exp(-rsq / sd2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void rho_grad(double[] xyz, double dfc, RefinementMode refinementmode) {
        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES) {
            return;
        }
        VectorMath.diff(this.xyz, xyz, dxyz);
        double r2 = VectorMath.rsq(dxyz);
        double sd2 = sd * sd;

        double rho = Math.exp(-r2 / sd2);

        g[0] = dfc * (2.0 * rho * -dxyz[0] / sd2);
        g[1] = dfc * (2.0 * rho * -dxyz[1] / sd2);
        g[2] = dfc * (2.0 * rho * -dxyz[2] / sd2);

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
