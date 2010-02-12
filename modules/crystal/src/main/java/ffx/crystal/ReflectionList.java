/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.rint;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Uses some methods from:
 * Cowtan, K. 2002. Generic representation and evaluation of properties as a
 * function of position in reciprocal space. J. Appl. Cryst. 35:655-663.
 *
 * @author fennt
 */
public class ReflectionList {

    public final HashMap<String, HKL> hklmap = new HashMap<String, HKL>();
    public final ArrayList<HKL> hkllist = new ArrayList<HKL>();
    public final Crystal crystal;
    public final SpaceGroup spaceGroup;
    private final SpaceGroup.CrystalSystem crystalSystem;
    private final SpaceGroup.LaueSystem laueSystem;
    public final Resolution resolution;
    // for binning reflections based on resolution
    public double hist[] = new double[1001];
    public double minres, maxres;

    public ReflectionList(Crystal crystal, Resolution resolution) {
        this.crystal = crystal;
        spaceGroup = crystal.spaceGroup;
        crystalSystem = spaceGroup.crystalSystem;
        laueSystem = spaceGroup.laueSystem;
        this.resolution = resolution;

        int hmax = (int) (this.crystal.a / this.resolution.res_limit());
        int kmax = (int) (this.crystal.b / this.resolution.res_limit());
        int lmax = (int) (this.crystal.c / this.resolution.res_limit());

        minres = Double.POSITIVE_INFINITY;
        maxres = Double.NEGATIVE_INFINITY;
        int n = 0;

        HKL hkl = new HKL();
        for (int h = -hmax; h <= hmax; h++) {
            hkl.h(h);
            for (int k = -kmax; k <= kmax; k++) {
                hkl.k(k);
                for (int l = -lmax; l <= lmax; l++) {
                    hkl.l(l);

                    double res = Crystal.invressq(this.crystal, hkl);
                    getepsilon(hkl);
                    if (SpaceGroup.checkLaueRestrictions(laueSystem, h, k, l)
                            && res < resolution.invressq_limit()
                            && !HKL.sys_abs(hkl)) {
                        minres = min(res, minres);
                        maxres = max(res, maxres);
                        String s = ("" + h + "_" + k + "_" + l).intern();
                        hklmap.put(s, new HKL(hkl.h(), hkl.k(), hkl.l(), hkl.epsilon(), hkl.allowed));
                        n++;
                    }
                }
            }
        }

        n=0;
        for (Iterator i = hklmap.entrySet().iterator(); i.hasNext();) {
            Map.Entry ei = (Map.Entry) i.next();
            Object key = ei.getKey();
            HKL ih = (HKL) ei.getValue();

            ih.index(n);
            hkllist.add(ih);
            n++;
        }

        /*
         * set up the resolution bins
         * first build a histogram
         */
        for (HKL ih : hkllist) {
            double r = (Crystal.invressq(this.crystal, ih) - minres) / (maxres - minres);
            int i = (int) (min(r, 0.999) * 1000.0);
            hist[i + 1] += 1.0;
        }

        // convert to cumulative histogram
        for (int i = 1; i < hist.length; i++) {
            hist[i] += hist[i - 1];
        }
        for (int i = 0; i < hist.length; i++) {
            hist[i] /= hist[hist.length - 1];
        }

        // assign each reflection to a bin in the range (0-9)
        for (HKL ih : hkllist) {
            int bin = (int) (10.0 * ordinal(Crystal.invressq(this.crystal, ih)));
            ih.bin = min(bin, 9);
        }
    }

    public ReflectionList(double a, double b, double c,
            double alpha, double beta, double gamma, String sg,
            double resolution) {
        this(new Crystal(a, b, c, alpha, beta, gamma, sg),
                new Resolution(resolution));
    }

    public HKL getHKL(int h, int k, int l) {
        String s = ("" + h + "_" + k + "_" + l);
        return hklmap.get(s);
    }

    public HKL getHKL(HKL hkl) {
        return getHKL(hkl.h(), hkl.k(), hkl.l());
    }

    public boolean hasHKL(int h, int k, int l) {
        String s = ("" + h + "_" + k + "_" + l);
        if (hklmap.containsKey(s)) {
            return true;
        } else {
            return false;
        }
    }

    public boolean hasHKL(HKL hkl) {
        return hasHKL(hkl.h(), hkl.k(), hkl.l());
    }

    public void getepsilon(HKL hkl) {
        int epsilon = 1;
        int allowed = 255;

        int nsym = spaceGroup.symOps.size();
        for (int i = 1; i < nsym; i++) {
            HKL mate = new HKL();
            crystal.applySymRot(hkl, mate, spaceGroup.symOps.get(i));
            double shift = Crystal.sym_phase_shift(hkl, spaceGroup.symOps.get(i));

            if (mate.equals(hkl)) {
                if (cos(shift) > 0.999) {
                    epsilon++;
                } else {
                    allowed = 0;
                    epsilon = 0;
                }
            } else if (mate.equals(HKL.neg(hkl))) {
                // centric reflection
                allowed = (int) rint(Crystal.mod(-0.5 * shift, PI) / (PI / 12.0));
            }
        }
        if (hkl.h() == 0 && hkl.k() == 0 && hkl.l() == 0) {
            allowed = 0;
        }

        hkl.epsilon(epsilon);
        hkl.allowed(allowed);
    }

    public double ordinal(double s) {
        double r = (s - minres) / (maxres - minres);
        r = min(r, 0.999) * 1000.0;
        int i = (int) r;
        r -= floor(r);
        return ((1.0 - r) * hist[i] + r * hist[i + 1]);
    }
}
