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
package ffx.crystal;

import static java.lang.Math.floor;
import static java.lang.Math.min;
import static java.lang.Math.max;

/**
 *
 * @author fennt
 */
public class ReflectionSpline {

    private final ReflectionList reflectionlist;
    private final Crystal crystal;
    private final int nparams;
    private double f;
    private int i0, i1, i2;
    private double dfi0, dfi1, dfi2;

    public ReflectionSpline(ReflectionList reflectionlist, int nparams) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.nparams = nparams;
    }

    public double f() {
        return f;
    }

    public int i0() {
        return i0;
    }

    public int i1() {
        return i1;
    }

    public int i2() {
        return i2;
    }

    public double dfi0() {
        return dfi0;
    }

    public double dfi1() {
        return dfi1;
    }

    public double dfi2() {
        return dfi2;
    }

    public double f(double invressq, double params[]) {
        double s = nparams * reflectionlist.ordinal(invressq);
        int i = (int) floor(s);
        double ds = s - i - 0.5;
        i0 = min(max(0, i - 1), nparams - 1);
        i1 = min(max(0, i), nparams - 1);
        i2 = min(max(0, i + 1), nparams - 1);

        f = params[i0] * 0.5 * (ds - 0.5) * (ds - 0.5)
                + params[i1] * (0.75 - ds * ds)
                + params[i2] * 0.5 * (ds + 0.5) * (ds + 0.5);

        dfi0 = 0.5 * (ds - 0.5) * (ds - 0.5);
        dfi1 = 0.75 - ds * ds;
        dfi2 = 0.5 * (ds + 0.5) * (ds + 0.5);

        return f;
    }
}
