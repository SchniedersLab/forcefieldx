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
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The 6 coefficients of the multiplicative polynomial switch are unique given
 * the distances "off" and "cut". They are found by solving a system of 6
 * equations, which define the boundary conditions of the switch.
 * <br>
 * f(cut) = 1
 * <br>
 * f'(cut) = f"(cut) = 0
 * <br>
 * f(off) = f'(off) = f"(off) = 0
 *
 * @author Michael J. Schnieders
 */
public class MultiplicativeSwitch {

    private final double c0;
    private final double c1;
    private final double c2;
    private final double c3;
    private final double c4;
    private final double c5;
    private final double twoC2;
    private final double threeC3;
    private final double fourC4;
    private final double fiveC5;

    private final static Logger logger = Logger.getLogger(MultiplicativeSwitch.class.getName());
    
    public MultiplicativeSwitch(double off, double cut) {

        double off2 = off * off;
        double cut2 = cut * cut;

        double denom = pow(off - cut, 5.0);
        c0 = off * off2 * (off2 - 5.0 * off * cut + 10.0 * cut2) / denom;
        c1 = -30.0 * off2 * cut2 / denom;
        c2 = 30.0 * (off2 * cut + off * cut2) / denom;
        c3 = -10.0 * (off2 + 4.0 * off * cut + cut2) / denom;
        c4 = 15.0 * (off + cut) / denom;
        c5 = -6.0 / denom;
        twoC2 = 2.0 * c2;
        threeC3 = 3.0 * c3;
        fourC4 = 4.0 * c4;
        fiveC5 = 5.0 * c5;
    }
    
    public double taper(double r) {
        return taper(r, r*r, r*r*r, r*r*r*r, r*r*r*r*r);
    }
    
    public double dtaper(double r) {
        return dtaper(r, r*r, r*r*r, r*r*r*r);
    }

    public double taper(double r, double r2, double r3, double r4, double r5) {
        return c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
    }

    public double dtaper(double r, double r2, double r3, double r4) {
        return fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
    }

}
