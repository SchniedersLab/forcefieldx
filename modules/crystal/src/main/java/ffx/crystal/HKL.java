//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import static org.apache.commons.math3.util.FastMath.PI;

import ffx.utilities.HashCodeUtil;

/**
 * <p>
 * The HKL class represents a single reflection.</p>
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 *
 * @see ReflectionList
 */
public class HKL {

    /**
     * Constant <code>ndiv=12.0</code>
     */
    public static final double ndiv = 12.0;
    protected int h;
    protected int k;
    protected int l;
    protected int epsilon;
    protected int allowed;
    protected int bin;
    protected int index;
    private int hashCode;

    /**
     * <p>
     * Constructor for HKL.</p>
     */
    public HKL() {
    }

    /**
     * <p>
     * Constructor for HKL.</p>
     *
     * @param h a int.
     * @param k a int.
     * @param l a int.
     */
    public HKL(int h, int k, int l) {
        this.h = h;
        this.k = k;
        this.l = l;
        allowed = 255;
    }

    /**
     * <p>
     * Constructor for HKL.</p>
     *
     * @param h a int.
     * @param k a int.
     * @param l a int.
     * @param eps a int.
     * @param allowed a int.
     */
    public HKL(int h, int k, int l, int eps, int allowed) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.epsilon = eps;
        this.allowed = allowed;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof HKL)) {
            return false;
        }
        if (this == obj) {
            return true;
        }
        return (this.h() == ((HKL) obj).h()
                && this.k() == ((HKL) obj).k()
                && this.l() == ((HKL) obj).l());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        if (hashCode == 0) {
            int result = HashCodeUtil.SEED;
            result = HashCodeUtil.hash(result, h);
            result = HashCodeUtil.hash(result, k);
            result = HashCodeUtil.hash(result, l);
            hashCode = result;
        }
        return hashCode;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return this.h() + " " + this.k() + " " + this.l()
                + "(allowed: " + this.allowed + " eps: " + this.epsilon + ") ";
    }

    /**
     * <p>
     * h</p>
     *
     * @return a int.
     */
    public int h() {
        return this.h;
    }

    /**
     * <p>
     * h</p>
     *
     * @param h a int.
     */
    public void h(int h) {
        this.h = h;
    }

    /**
     * <p>
     * k</p>
     *
     * @return a int.
     */
    public int k() {
        return this.k;
    }

    /**
     * <p>
     * k</p>
     *
     * @param k a int.
     */
    public void k(int k) {
        this.k = k;
    }

    /**
     * <p>
     * l</p>
     *
     * @return a int.
     */
    public int l() {
        return this.l;
    }

    /**
     * <p>
     * l</p>
     *
     * @param l a int.
     */
    public void l(int l) {
        this.l = l;
    }

    /**
     * <p>
     * bin</p>
     *
     * @return a int.
     */
    public int bin() {
        return this.bin;
    }

    /**
     * <p>
     * bin</p>
     *
     * @param bin a int.
     */
    public void bin(int bin) {
        this.bin = bin;
    }

    /**
     * <p>
     * index</p>
     *
     * @return a int.
     */
    public int index() {
        return this.index;
    }

    /**
     * <p>
     * index</p>
     *
     * @param index a int.
     */
    public void index(int index) {
        this.index = index;
    }

    /**
     * <p>
     * neg</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a {@link ffx.crystal.HKL} object.
     */
    public static HKL neg(HKL hkl) {
        return new HKL(-hkl.h(), -hkl.k(), -hkl.l());
    }

    /**
     * <p>
     * epsilon</p>
     *
     * @return a int.
     */
    public int epsilon() {
        return this.epsilon;
    }

    /**
     * <p>
     * epsilon</p>
     *
     * @param eps a int.
     */
    public void epsilon(int eps) {
        this.epsilon = eps;
    }

    /**
     * <p>
     * epsilonc</p>
     *
     * @return a int.
     */
    public int epsilonc() {
        if (this.centric()) {
            return 2 * epsilon;
        } else {
            return epsilon;
        }
    }

    /**
     * <p>
     * allowed</p>
     *
     * @return a double.
     */
    public double allowed() {
        return ((double) this.allowed) * (PI / ndiv);
    }

    /**
     * <p>
     * allowed</p>
     *
     * @param allowed a int.
     */
    public void allowed(int allowed) {
        this.allowed = allowed;
    }

    /**
     * <p>
     * centric</p>
     *
     * @return a boolean.
     */
    public boolean centric() {
        return (allowed != 255);
    }

    /**
     * <p>
     * sys_abs</p>
     *
     * @param hkl a {@link ffx.crystal.HKL} object.
     * @return a boolean.
     */
    public static boolean sys_abs(HKL hkl) {
        return (hkl.epsilon == 0);
    }
}
