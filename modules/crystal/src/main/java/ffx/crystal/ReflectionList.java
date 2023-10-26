// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.crystal;

import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import static ffx.crystal.SymOp.applySymRot;
import static ffx.crystal.SymOp.applyTransSymRot;
import static ffx.numerics.math.ScalarMath.mod;
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
 * @see <a href="http://dx.doi.org/10.1107/S0021889802013420" target="_blank"> Cowtan, K. 2002.
 * Generic representation and evaluation of properties as a function of position in reciprocal
 * space. J. Appl. Cryst. 35:655-663. </a>
 * @since 1.0
 */
public class ReflectionList {

  /**
   * The HKL list.
   */
  public final ArrayList<HKL> hklList = new ArrayList<>();
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
   * String to HKL look-up.
   */
  final HashMap<String, HKL> hklMap = new HashMap<>();
  /**
   * For binning reflections based on resolution
   */
  public int nBins = 10;
  /**
   * Histogram.
   */
  private final double[] histogram = new double[1001];
  /**
   * Minimum resolution.
   */
  private double minResolution;
  /**
   * Maximum resolution.
   */
  private double maxResolution;

  /**
   * Constructor for ReflectionList.
   *
   * @param crystal    a {@link ffx.crystal.Crystal} object.
   * @param resolution a {@link ffx.crystal.Resolution} object.
   */
  public ReflectionList(Crystal crystal, Resolution resolution) {
    this(crystal, resolution, null);
  }

  /**
   * Constructor for ReflectionList.
   *
   * @param crystal    a {@link ffx.crystal.Crystal} object.
   * @param resolution a {@link ffx.crystal.Resolution} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public ReflectionList(Crystal crystal, Resolution resolution, CompositeConfiguration properties) {
    this.crystal = crystal;
    this.spaceGroup = crystal.spaceGroup;
    this.resolution = resolution;

    int hMax = (int) (this.crystal.a / this.resolution.resolutionLimit());
    int kMax = (int) (this.crystal.b / this.resolution.resolutionLimit());
    int lMax = (int) (this.crystal.c / this.resolution.resolutionLimit());

    minResolution = Double.POSITIVE_INFINITY;
    maxResolution = Double.NEGATIVE_INFINITY;
    int n = 0;

    HKL hkl = new HKL();
    for (int h = -hMax; h <= hMax; h++) {
      hkl.setH(h);
      for (int k = -kMax; k <= kMax; k++) {
        hkl.setK(k);
        for (int l = -lMax; l <= lMax; l++) {
          hkl.setL(l);

          double res = this.crystal.invressq(hkl);
          getEpsilon(hkl);
          LaueSystem laueSystem = spaceGroup.laueSystem;
          if (laueSystem.checkRestrictions(h, k, l)
              && resolution.inInverseResSqRange(res)
              && !hkl.sysAbs()) {
            minResolution = min(res, minResolution);
            maxResolution = max(res, maxResolution);
            String s = (h + "_" + k + "_" + l).intern();
            hklMap.put(s, new HKL(hkl.getH(), hkl.getK(), hkl.getL(), hkl.getEpsilon(), hkl.allowed));
            n++;
          }
        }
      }
    }

    n = 0;
    for (Entry<String, HKL> entry : hklMap.entrySet()) {
      HKL ih = entry.getValue();
      ih.setIndex(n);
      hklList.add(ih);
      n++;
    }

    // Set up the resolution bins first build a histogram.
    for (HKL ih : hklList) {
      double r = (this.crystal.invressq(ih) - minResolution) / (maxResolution - minResolution);
      int i = (int) (min(r, 0.999) * 1000.0);
      histogram[i + 1] += 1.0;
    }

    // Convert to cumulative histogram
    for (int i = 1; i < histogram.length; i++) {
      histogram[i] += histogram[i - 1];
    }
    for (int i = 0; i < histogram.length; i++) {
      histogram[i] /= histogram[histogram.length - 1];
    }

    // Assign each reflection to a bin in the range (0-nbins)
    setResolutionBins(properties);
  }

  /**
   * Constructor for ReflectionList.
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
  public ReflectionList(double a, double b, double c, double alpha, double beta, double gamma,
                        String sg, double resolution) {
    this(new Crystal(a, b, c, alpha, beta, gamma, sg), new Resolution(resolution));
  }

  /**
   * findSymHKL
   *
   * @param hkl  a {@link ffx.crystal.HKL} object.
   * @param mate a {@link ffx.crystal.HKL} object.
   * @return a boolean.
   */
  public boolean findSymHKL(HKL hkl, HKL mate) {
    return findSymHKL(hkl, mate, false);
  }

  /**
   * findSymHKL
   *
   * @param h    an int.
   * @param k    an int.
   * @param l    an int.
   * @param mate a {@link ffx.crystal.HKL} object.
   * @return a boolean.
   */
  public boolean findSymHKL(int h, int k, int l, HKL mate) {
    return findSymHKL(new HKL(h, k, l), mate, false);
  }

  /**
   * findSymHKL
   *
   * @param h         an int.
   * @param k         an int.
   * @param l         an int.
   * @param mate      a {@link ffx.crystal.HKL} object.
   * @param transpose a boolean.
   * @return a boolean.
   */
  public boolean findSymHKL(int h, int k, int l, HKL mate, boolean transpose) {
    return findSymHKL(new HKL(h, k, l), mate, transpose);
  }

  /**
   * getHKL
   *
   * @param h an int.
   * @param k an int.
   * @param l an int.
   * @return a {@link ffx.crystal.HKL} object.
   */
  public HKL getHKL(int h, int k, int l) {
    String s = (h + "_" + k + "_" + l);
    return hklMap.get(s);
  }

  /**
   * getHKL
   *
   * @param hkl a {@link ffx.crystal.HKL} object.
   * @return a {@link ffx.crystal.HKL} object.
   */
  public HKL getHKL(HKL hkl) {
    return getHKL(hkl.getH(), hkl.getK(), hkl.getL());
  }

  /**
   * Get the maximum resolution.
   *
   * @return Maximum resolution
   */
  public double getMaxResolution() {
    return maxResolution;
  }

  /**
   * Get the minimum resolution.
   *
   * @return Minimum resolution
   */
  public double getMinResolution() {
    return minResolution;
  }

  /**
   * hasHKL
   *
   * @param hkl a {@link ffx.crystal.HKL} object.
   * @return a boolean.
   */
  public boolean hasHKL(HKL hkl) {
    return hasHKL(hkl.getH(), hkl.getK(), hkl.getL());
  }

  /**
   * ordinal
   *
   * @param s a double.
   * @return a double.
   */
  public final double ordinal(double s) {
    double r = (s - minResolution) / (maxResolution - minResolution);
    r = min(r, 0.999) * 1000.0;
    int i = (int) r;
    r -= floor(r);
    return ((1.0 - r) * histogram[i] + r * histogram[i + 1]);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return " Reflection list with "
        + this.hklList.size()
        + " reflections, spacegroup "
        + this.spaceGroup.shortName
        + " resolution limit: "
        + resolution.resolutionLimit();
  }

  /**
   * findSymHKL
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
        applyTransSymRot(hkl, mate, spaceGroup.symOps.get(i));
      } else {
        applySymRot(hkl, mate, spaceGroup.symOps.get(i));
      }

      LaueSystem laueSystem = spaceGroup.laueSystem;
      if (laueSystem.checkRestrictions(mate.getH(), mate.getK(), mate.getL())) {
        return false;
      }
      if (laueSystem.checkRestrictions(-mate.getH(), -mate.getK(), -mate.getL())) {
        mate.setH(-mate.getH());
        mate.setK(-mate.getK());
        mate.setL(-mate.getL());
        return true;
      }
    }

    mate.setH(hkl.getH());
    mate.setK(hkl.getK());
    mate.setL(hkl.getL());
    return false;
  }

  /**
   * hasHKL
   *
   * @param h an int.
   * @param k an int.
   * @param l an int.
   * @return a boolean.
   */
  private boolean hasHKL(int h, int k, int l) {
    String s = (h + "_" + k + "_" + l);
    return hklMap.containsKey(s);
  }

  private void getEpsilon(HKL hkl) {
    int epsilon = 1;
    int allowed = 255;

    int nSym = spaceGroup.symOps.size();
    for (int i = 1; i < nSym; i++) {
      HKL mate = new HKL();
      SymOp symOp = spaceGroup.symOps.get(i);
      applySymRot(hkl, mate, symOp);
      double shift = symOp.symPhaseShift(hkl);
      if (mate.equals(hkl)) {
        if (cos(shift) > 0.999) {
          epsilon++;
        } else {
          allowed = 0;
          epsilon = 0;
          break;
        }
      } else if (mate.equals(hkl.neg())) {
        // centric reflection
        allowed = (int) rint(mod(-0.5 * shift, PI) / (PI / HKL.ndiv));
      }
    }
    if (hkl.getH() == 0 && hkl.getK() == 0 && hkl.getL() == 0) {
      allowed = 0;
    }

    hkl.setEpsilon(epsilon);
    hkl.setAllowed(allowed);
  }

  private void setResolutionBins(@Nullable CompositeConfiguration properties) {
    if (properties != null) {
      nBins = properties.getInt("nbins", 10);
    }
    double nBinsDouble = nBins;
    for (HKL ih : hklList) {
      int bin = (int) (nBinsDouble * ordinal(crystal.invressq(ih)));
      ih.bin = min(bin, nBins - 1);
    }
  }
}
