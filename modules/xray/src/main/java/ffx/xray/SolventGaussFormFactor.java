/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import static org.apache.commons.math3.util.FastMath.exp;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>
 * SolventGaussFormFactor class.</p>
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
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
        this(atom, sd, atom.getXYZ(null));
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

    /** {@inheritDoc} */
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

    /** {@inheritDoc} */
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

    /** {@inheritDoc} */
    @Override
    public void update(double xyz[]) {
        update(xyz, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public void update(double xyz[], double badd) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
    }
}
