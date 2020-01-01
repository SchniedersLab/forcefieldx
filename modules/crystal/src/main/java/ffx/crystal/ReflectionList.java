//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.crystal;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.rint;

/**
 * Class to represent a reflection list.
 *
 * @author Timothy D. Fenn
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank">
 * Cowtan, K. 2002. Generic representation and evaluation of properties as a
 * function of position in reciprocal space. J. Appl. Cryst. 35:655-663.
 * </a>
 * @since 1.0
 */
public class ReflectionList {

    /**
     * The HKL list.
     */
    public final ArrayList<HKL> hkllist = new ArrayList<>();
    /**
     * The Crystal instance.
     */
    public final Crystal crystal;
    /**
     * The space group.
     */
    public final SpaceGroup spaceGroup;
    /**
     * Resolution instance.
     */
    public final Resolution resolution;
    /**
     * For binning reflections based on resolution
     */
    public int nbins = 10;
    /**
     * String to HKL look-up.
     */
    final HashMap<String, HKL> hklmap = new HashMap<>();
    /**
     * The Laue System.
     */
    private final SpaceGroup.LaueSystem laueSystem;
    /**
     * Histogram.
     */
    private double[] hist = new double[1001];
    /**
     * Minimum resolution.
     */
    private double minResolution;
    /**
     * Maximum resolution.
     */
    private double maxResolution;

    /**
     * <p>
     * Constructor for ReflectionList.</p>
     *
     * @param crystal    a {@link ffx.crystal.Crystal} object.
     * @param resolution a {@link ffx.crystal.Resolution} object.
     */
    public ReflectionList(Crystal crystal, Resolution resolution) {
        this(crystal, resolution, null);
    }

    /**
     * <p>
     * Constructor for ReflectionList.</p>
     *
     * @param crystal    a {@link ffx.crystal.Crystal} object.
     * @param resolution a {@link ffx.crystal.Resolution} object.
     * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public ReflectionList(Crystal crystal, Resolution resolution,
                          CompositeConfiguration properties) {
        this.crystal = crystal;
        spaceGroup = crystal.spaceGroup;
        laueSystem = spaceGroup.laueSystem;
        this.resolution = resolution;

        int hmax = (int) (this.crystal.a / this.resolution.resolutionLimit());
        int kmax = (int) (this.crystal.b / this.resolution.resolutionLimit());
        int lmax = (int) (this.crystal.c / this.resolution.resolutionLimit());

        minResolution = Double.POSITIVE_INFINITY;
        maxResolution = Double.NEGATIVE_INFINITY;
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
                            && resolution.inInverseResSqRange(res)
                            && !HKL.sys_abs(hkl)) {
                        minResolution = min(res, minResolution);
                        maxResolution = max(res, maxResolution);
                        String s = ("" + h + "_" + k + "_" + l).intern();
                        hklmap.put(s, new HKL(hkl.h(), hkl.k(), hkl.l(), hkl.epsilon(), hkl.allowed));
                        n++;
                    }
                }
            }
        }

        n = 0;
        for (Map.Entry ei : hklmap.entrySet()) {
            HKL ih = (HKL) ei.getValue();
            ih.index(n);
            hkllist.add(ih);
            n++;
        }

        // Set up the resolution bins first build a histogram.
        for (HKL ih : hkllist) {
            double r = (Crystal.invressq(this.crystal, ih) - minResolution) / (maxResolution - minResolution);
            int i = (int) (min(r, 0.999) * 1000.0);
            hist[i + 1] += 1.0;
        }

        // Convert to cumulative histogram
        for (int i = 1; i < hist.length; i++) {
            hist[i] += hist[i - 1];
        }
        for (int i = 0; i < hist.length; i++) {
            hist[i] /= hist[hist.length - 1];
        }

        // Assign each reflection to a bin in the range (0-nbins)
        setResolutionBins(properties);
    }

    /**
     * <p>
     * Constructor for ReflectionList.</p>
     *
     * @param a          a double.
     * @param b          a double.
     * @param c          a double.
     * @param alpha      a double.
     * @param beta       a double.
     * @param gamma      a double.
     * @param sg         a {@link java.lang.String} object.
     * @param resolution a double.
     */
    public ReflectionList(double a, double b, double c,
                          double alpha, double beta, double gamma, String sg,
                          double resolution) {
        this(new Crystal(a, b, c, alpha, beta, gamma, sg),
                new Resolution(resolution));
    }

    /**
     * Get the maximum resolution.
     * @return Maximum resolution
     */
    public double getMaxResolution() {
        return maxResolution;
    }

    /**
     * Get the minimum resolution.
     * @return Minimum resolution
     */
    public double getMinResolution() {
        return minResolution;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return " Reflection list with " + this.hkllist.size()
                + " reflections, spacegroup " + this.spaceGroup.shortName
                + " resolution limit: " + resolution.resolutionLimit();
    }

    /**
     * <p>
     * findSymHKL</p>
     *
     * @param hkl  a {@link ffx.crystal.HKL} object.
     * @param mate a {@link ffx.crystal.HKL} object.
     * @return a boolean.
     */
    public boolean findSymHKL(HKL hkl, HKL mate) {
        return findSymHKL(hkl, mate, false);
    }

    /**
     * <p>
     * findSymHKL</p>
     *
     * @param h    a int.
     * @param k    a int.
     * @param l    a int.
     * @param mate a {@link ffx.crystal.HKL} object.
     * @return a boolean.
     */
    public boolean findSymHKL(int h, int k, int l, HKL mate) {
        return findSymHKL(new HKL(h, k, l), mate, false);
    }

    /**
     * <p>
     * findSymHKL</p>
     *
     * @param h         a int.
     * @param k         a int.
     * @param l         a int.
     * @param mate      a {@link ffx.crystal.HKL} object.
     * @param transpose a boolean.
     * @return a boolean.
     */
    public boolean findSymHKL(int h, int k, int l, HKL mate, boolean transpose) {
        return findSymHKL(new HKL(h, k, l), mate, transpose);
    }

    /**
     * <p>
     * findSymHKL</p>
     *
     * @param hkl       a {@link ffx.crystal.HKL} object.
     * @param mate      a {@link ffx.crystal.HKL} object.
     * @param transpose a boolean.
     * @return a boolean.
     */
    private boolean findSymHKL(HKL hkl, HKL mate, boolean transpose) {
        int nsym = spaceGroup.numPrimitiveSymEquiv;

        for (int i = 0; i < nsym; i++) {
            if (transpose) {
                crystal.applyTransSymRot(hkl, mate, spaceGroup.symOps.get(i));
            } else {
                crystal.applySymRot(hkl, mate, spaceGroup.symOps.get(i));
            }
            if (SpaceGroup.checkLaueRestrictions(laueSystem,
                    mate.h(), mate.k(), mate.l())) {
                return false;
            }
            if (SpaceGroup.checkLaueRestrictions(laueSystem,
                    -mate.h(), -mate.k(), -mate.l())) {
                mate.h(-mate.h());
                mate.k(-mate.k());
                mate.l(-mate.l());
                return true;
            }
        }

        mate.h(hkl.h());
        mate.k(hkl.k());
        mate.l(hkl.l());
        return false;
    }

    /**
     * <p>
     * getHKL</p>
     *
     * @param h a int.
     * @param k a int.
     * @param l a int.
     * @return a {@link ffx.crystal.HKL} object.
     */
    public HKL getHKL(int h, int k, int l) {
        String s = ("" + h + "_" + k + "_" + l);
        return hklmap.get(s);
    }

    /**
     * <p>
     * getHKL</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a {@link ffx.crystal.HKL} object.
     */
    public HKL getHKL(HKL hkl) {
        return getHKL(hkl.h(), hkl.k(), hkl.l());
    }

    /**
     * <p>
     * hasHKL</p>
     *
     * @param h a int.
     * @param k a int.
     * @param l a int.
     * @return a boolean.
     */
    private boolean hasHKL(int h, int k, int l) {
        String s = ("" + h + "_" + k + "_" + l);
        return hklmap.containsKey(s);
    }

    /**
     * <p>
     * hasHKL</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a boolean.
     */
    public boolean hasHKL(HKL hkl) {
        return hasHKL(hkl.h(), hkl.k(), hkl.l());
    }

    private void getepsilon(HKL hkl) {
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
                    break;
                }
            } else if (mate.equals(HKL.neg(hkl))) {
                // centric reflection
                allowed = (int) rint(Crystal.mod(-0.5 * shift, PI) / (PI / HKL.ndiv));
            }
        }
        if (hkl.h() == 0 && hkl.k() == 0 && hkl.l() == 0) {
            allowed = 0;
        }

        hkl.epsilon(epsilon);
        hkl.allowed(allowed);
    }

    /**
     * <p>
     * ordinal</p>
     *
     * @param s a double.
     * @return a double.
     */
    public double ordinal(double s) {
        double r = (s - minResolution) / (maxResolution - minResolution);
        r = min(r, 0.999) * 1000.0;
        int i = (int) r;
        r -= floor(r);
        return ((1.0 - r) * hist[i] + r * hist[i + 1]);
    }

    private void setResolutionBins() {
        setResolutionBins(null);
    }

    private void setResolutionBins(CompositeConfiguration properties) {
        if (properties != null) {
            nbins = properties.getInt("nbins", 10);
        }
        double nbinsd = (double) nbins;
        for (HKL ih : hkllist) {
            int bin = (int) (nbinsd * ordinal(Crystal.invressq(this.crystal, ih)));
            ih.bin = min(bin, nbins - 1);
        }
    }
}
