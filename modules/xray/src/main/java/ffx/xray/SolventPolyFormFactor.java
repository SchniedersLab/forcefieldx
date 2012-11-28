/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 * <p>SolventPolyFormFactor class.</p>
 *
 * @author fenn
 *
 */
public final class SolventPolyFormFactor implements FormFactor {

    private final Atom atom;
    private double xyz[] = new double[3];
    private double dxyz[] = new double[3];
    private double g[] = new double[3];
    private double arad, w;

    /**
     * <p>Constructor for SolventPolyFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param arad a double.
     * @param w a double.
     */
    public SolventPolyFormFactor(Atom atom, double arad, double w) {
        this(atom, arad, w, atom.getXYZ());
    }

    /**
     * <p>Constructor for SolventPolyFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param arad a double.
     * @param w a double.
     * @param xyz an array of double.
     */
    public SolventPolyFormFactor(Atom atom, double arad, double w, double xyz[]) {
        this.atom = atom;
        this.arad = arad;
        this.w = w;

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
     * <p>rho</p>
     *
     * @param f a double.
     * @param lambda a double.
     * @param ri a double.
     * @return a double.
     */
    public double rho(double f, double lambda, double ri) {
        double bi = arad - w;
        double ei = arad + w;
        if (ri <= bi) {
            return 0.0;
        }
        if (ri >= ei) {
            return f * 1.0;
        }

        double d = ri - arad + w;
        double d2 = d * d;
        double d3 = d2 * d;
        double w2 = w * w;
        double w3 = w2 * w;

        return f * (0.75 * d2 / w2 - 0.25 * d3 / w3);
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
        double ri = VectorMath.r(dxyz);
        double bi = arad - w;
        double ei = arad + w;
        if (ri <= bi) {
            return;
        }
        if (ri >= ei) {
            return;
        }

        double d = ri - arad + w;
        double d2 = d * d;
        double d3 = d2 * d;
        double w2 = w * w;
        double w3 = w2 * w;

        double rho = 0.75 * d2 / w2 - 0.25 * d3 / w3;

        double dp = 1.5 * d / (w2 * ri) - 0.75 * d2 / (w3 * ri);

        g[0] = (dfc / rho) * (-dp * dxyz[0]);
        g[1] = (dfc / rho) * (-dp * dxyz[1]);
        g[2] = (dfc / rho) * (-dp * dxyz[2]);

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
