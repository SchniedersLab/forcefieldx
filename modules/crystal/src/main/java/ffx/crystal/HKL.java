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

import static java.lang.Math.PI;

import ffx.utilities.HashCodeUtil;

/**
 *
 * @author fennt
 */
public class HKL {

    public int h;
    public int k;
    public int l;
    public int epsilon;
    public int allowed;
    private int hashCode;

    // null constructor
    public HKL() {
    }

    public HKL(int h, int k, int l) {
        this.h = h;
        this.k = k;
        this.l = l;
        allowed = 255;
    }

    public HKL(int h, int k, int l, int eps, int allowed) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.epsilon = eps;
        this.allowed = allowed;
    }

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

    public int h() {
        return this.h;
    }

    public void h(int h) {
        this.h = h;
    }

    public int k() {
        return this.k;
    }

    public void k(int k) {
        this.k = k;
    }

    public int l() {
        return this.l;
    }

    public void l(int l) {
        this.l = l;
    }

    public static HKL neg(HKL hkl) {
        return new HKL(-hkl.h(), -hkl.k(), -hkl.l());
    }

    public int epsilon() {
        return this.epsilon;
    }

    public void epsilon(int eps) {
        this.epsilon = eps;
    }

    public int epsilonc() {
        if (this.centric()) {
            return 2 * epsilon;
        } else {
            return epsilon;
        }
    }

    public double allowed() {
        return ((double) this.allowed) * (PI / 12.0);
    }

    public void allowed(int allowed) {
        this.allowed = allowed;
    }

    public boolean centric() {
        return (allowed != 255);
    }

    public static boolean sys_abs(HKL hkl) {
        return (hkl.epsilon == 0);
    }
}
