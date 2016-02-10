/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray;

import static java.lang.Math.sqrt;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.numerics.VectorMath.diff;

/**
 * <p>
 * SolventPolyFormFactor class.</p>
 *
 * @author fenn
 *
 */
public final class SolventPolyFormFactor implements FormFactor {

    private final Atom atom;
    private final double xyz[] = new double[3];
    private final double dxyz[] = new double[3];
    private final double g[] = new double[3];
    private final double iw;
    private final double aradMinusW, aradPlusW;
    private final double aradMinusW2, aradPlusW2;
    private final double wMinusArad;

    /**
     * <p>
     * Constructor for SolventPolyFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param arad a double.
     * @param w a double.
     */
    public SolventPolyFormFactor(Atom atom, double arad, double w) {
        this(atom, arad, w, atom.getXYZ(null));
    }

    /**
     * <p>
     * Constructor for SolventPolyFormFactor.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @param arad a double.
     * @param w a double.
     * @param xyz an array of double.
     */
    public SolventPolyFormFactor(Atom atom, double arad, double w, double xyz[]) {
        this.atom = atom;
        this.iw = 1.0 / w;
        aradMinusW = arad - w;
        aradPlusW = arad + w;
        wMinusArad = w - arad;
        aradMinusW2 = aradMinusW * aradMinusW;
        aradPlusW2 = aradPlusW * aradPlusW;
        update(xyz);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double rho(double f, double lambda, double[] xyz) {
        VectorMath.diff(this.xyz, xyz, dxyz);
        double ri2 = VectorMath.rsq(dxyz);
        if (ri2 <= aradMinusW2) {
            return 0.0;
        }
        if (ri2 >= aradPlusW2) {
            return f;
        }
        return rho(f, lambda, sqrt(ri2));
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
        if (ri <= aradMinusW) {
            return 0.0;
        }
        if (ri >= aradPlusW) {
            return f;
        }
        double d = ri + wMinusArad;
        double dw = d * iw;
        double dw2 = dw * dw;
        return f * (0.75 - 0.25 * dw) * dw2;
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
        diff(this.xyz, xyz, dxyz);
        double ri2 = VectorMath.rsq(dxyz);
        if (ri2 <= aradMinusW2 || ri2 >= aradPlusW2) {
            return;
        }
        double ri = sqrt(ri2);
        double d = ri + wMinusArad;
        double dw = d * iw;
        double dw2 = dw * dw;
        double rho = (0.75 - 0.25 * dw) * dw2;
        double iri = 1.0 / ri;
        double dp = (1.5 * dw - 0.75 * dw2) * iri * iw;
        double prefactor = -dp * (dfc / rho);
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
