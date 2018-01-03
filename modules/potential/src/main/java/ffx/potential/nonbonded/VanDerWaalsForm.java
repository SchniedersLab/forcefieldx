/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;
import static ffx.potential.parameters.ForceField.ForceFieldString.EPSILONRULE;
import static ffx.potential.parameters.ForceField.ForceFieldString.RADIUSRULE;
import static ffx.potential.parameters.ForceField.ForceFieldString.RADIUSSIZE;
import static ffx.potential.parameters.ForceField.ForceFieldString.RADIUSTYPE;
import static ffx.potential.parameters.ForceField.ForceFieldString.VDWTYPE;
import static ffx.potential.parameters.ForceField.toEnumForm;

/**
 * This class contains fields and methods for maintaining details of the van der
 * Waals functional form.
 *
 * @author Michael J. Schnieders
 */
public class VanDerWaalsForm {

    /**
     * The logger.
     */
    private static final Logger logger = Logger.getLogger(VanDerWaalsForm.class.getName());

    public enum VDW_TYPE {

        BUFFERED_14_7, LENNARD_JONES
    }

    public enum RADIUS_RULE {
        ARITHMETIC, CUBIC_MEAN, GEOMETRIC
    }

    public enum RADIUS_SIZE {
        DIAMETER, RADIUS
    }

    public enum RADIUS_TYPE {
        R_MIN, SIGMA
    }

    public enum EPSILON_RULE {
        GEOMETRIC, HHG
    }

    public VanDerWaalsForm(ForceField forceField) {

        /**
         * Set-up default rules.
         */
        vdwType = VDW_TYPE.BUFFERED_14_7;
        epsilonRule = EPSILON_RULE.HHG;
        radiusRule = RADIUS_RULE.CUBIC_MEAN;
        radiusSize = RADIUS_SIZE.DIAMETER;
        radiusType = RADIUS_TYPE.R_MIN;

        /**
         * Define functional form.
         */
        String value = forceField.getString(VDWTYPE, vdwType.toString());
        try {
            vdwType = VDW_TYPE.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized VDWTYPE %s; defaulting to %s.", value, vdwType));
        }

        switch (vdwType) {
            case BUFFERED_14_7:
                vdwPowers = new Buffered_14_7();
                break;
            case LENNARD_JONES:
                vdwPowers = new LJ_6_12();
                break;
            default:
                vdwPowers = new VDWPowers();
                break;
        }

        /**
         * Define epsilon combining rule.
         */
        value = forceField.getString(EPSILONRULE, epsilonRule.toString());
        try {
            epsilonRule = EPSILON_RULE.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized EPSILONRULE %s; defaulting to %s.", value, epsilonRule));
        }

        /**
         * Define radius combining rule.
         */
        value = forceField.getString(RADIUSRULE, radiusRule.toString());
        try {
            radiusRule = RADIUS_RULE.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized RADIUSRULE %s; defaulting to %s.", value, radiusRule));
        }

        /**
         * Define radius size.
         */
        value = forceField.getString(RADIUSSIZE, radiusSize.toString());
        try {
            radiusSize = RADIUS_SIZE.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized RADIUSSIZE %s; defaulting to %s.", value, radiusSize));
        }

        /**
         * Define radius type.
         */
        value = forceField.getString(RADIUSTYPE, radiusType.toString());
        try {
            radiusType = RADIUS_TYPE.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized RADIUSTYPE %s; defaulting to %s", value, radiusType));
        }

        /**
         * Configure van der Waals well shape parameters.
         */
        switch (vdwType) {
            case LENNARD_JONES:
                repulsivePower = 12;
                dispersivePower = 6;
                delta = 0.0;
                gamma = 0.0;
                break;
            case BUFFERED_14_7:
            default:
                repulsivePower = 14;
                dispersivePower = 7;
                delta = 0.07;
                gamma = 0.12;
                break;
        }

        repDispPower = repulsivePower - dispersivePower;
        dispersivePower1 = dispersivePower - 1;
        repDispPower1 = repDispPower - 1;
        delta1 = 1.0 + delta;
        t1n = pow(delta1, dispersivePower);
        gamma1 = 1.0 + gamma;

        scale12 = forceField.getDouble(ForceField.ForceFieldDouble.VDW_12_SCALE, 0.0);
        scale13 = forceField.getDouble(ForceField.ForceFieldDouble.VDW_13_SCALE, 0.0);
        scale14 = forceField.getDouble(ForceField.ForceFieldDouble.VDW_14_SCALE, 1.0);
        scale15 = forceField.getDouble(ForceField.ForceFieldDouble.VDW_15_SCALE, 1.0);

        /**
         * The convention in TINKER is a vdw-14-scale factor of 2.0 means to
         * scale by 0.5.
         */
        if (scale12 > 1.0) {
            scale12 = 1.0 / scale12;
        }
        if (scale13 > 1.0) {
            scale13 = 1.0 / scale13;
        }
        if (scale14 > 1.0) {
            scale14 = 1.0 / scale14;
        }
        if (scale15 != 1.0) {
            logger.severe(" Van Der Waals 1-5 masking rules are not supported.");
        }

        Map<String, VDWType> map = forceField.getVDWTypes();
        TreeMap<String, VDWType> vdwTypes = new TreeMap<>(map);
        maxClass = 0;
        for (VDWType currentType : vdwTypes.values()) {
            if (currentType.atomClass > maxClass) {
                maxClass = currentType.atomClass;
            }
        }
        radEps = new double[maxClass + 1][2 * (maxClass + 1)];

        /**
         * Scale factor to convert to vdW size to Rmin.
         */
        double radScale;
        switch (radiusSize) {
            case DIAMETER:
                radScale = 0.5;
                break;
            case RADIUS:
            default:
                radScale = 1.0;
                break;
        }
        switch (radiusType) {
            case SIGMA:
                radScale *= 1.122462048309372981;
                break;
            case R_MIN:
            default:
                break;

        }

        /**
         * Atom Class numbering starts at 1.
         */
        for (VDWType vdwi : vdwTypes.values()) {
            int i = vdwi.atomClass;
            double ri = radScale * vdwi.radius;
            double ri2 = ri * ri;
            double ri3 = ri * ri2;
            double e1 = vdwi.wellDepth;
            double se1 = sqrt(e1);
            for (VDWType vdwj : vdwTypes.tailMap(vdwi.getKey()).values()) {
                int j = vdwj.atomClass;
                double rj = radScale * vdwj.radius;
                double rj2 = rj * rj;
                double rj3 = rj * rj2;
                double e2 = vdwj.wellDepth;
                double se2 = sqrt(e2);
                double radmin;
                double eps;
                switch (radiusRule) {
                    case ARITHMETIC:
                        radmin = ri + rj;
                        break;
                    case GEOMETRIC:
                        radmin = 2.0 * sqrt(ri) * sqrt(rj);
                        break;
                    default:
                    case CUBIC_MEAN:
                        radmin = 2.0 * (ri3 + rj3) / (ri2 + rj2);
                }
                switch (epsilonRule) {
                    case GEOMETRIC:
                        eps = se1 * se2;
                        break;
                    default:
                    case HHG:
                        eps = 4.0 * (e1 * e2) / ((se1 + se2) * (se1 + se2));
                        break;
                }
                if (radmin > 0) {
                    radEps[i][j * 2 + RADMIN] = 1.0 / radmin;
                } else {
                    radEps[i][j * 2 + RADMIN] = 0.0;
                }

                radEps[i][j * 2 + EPS] = eps;
                if (radmin > 0) {
                    radEps[j][i * 2 + RADMIN] = 1.0 / radmin;
                } else {
                    radEps[j][i * 2 + RADMIN] = 0.0;
                }
                radEps[j][i * 2 + EPS] = eps;
            }
        }

    }

    public double rhoDisp1(double rho) {
        return vdwPowers.rhoDisp1(rho);
    }

    public double rhoDelta1(double rhoDelta) {
        return vdwPowers.rhoDisp1(rhoDelta);
    }

    /**
     * van der Waals functional form.
     */
    public VDW_TYPE vdwType;
    public EPSILON_RULE epsilonRule = EPSILON_RULE.HHG;
    public RADIUS_RULE radiusRule = RADIUS_RULE.CUBIC_MEAN;
    public RADIUS_SIZE radiusSize = RADIUS_SIZE.DIAMETER;
    public RADIUS_TYPE radiusType = RADIUS_TYPE.R_MIN;

    private final VDWPowers vdwPowers;

    /**
     * vdW Repulsive Power (e.g. 12).
     */
    public final int repulsivePower;
    /**
     * vdW Dispersive Power (e.g. 6).
     */
    public final int dispersivePower;

    /**
     * First constant suggested by Halgren for the Buffered-14-7 potential.
     */
    public final double gamma;
    public final double gamma1;

    /**
     * Second constant suggested by Halgren for the Buffered-14-7 potential.
     */
    public final double delta;
    public final double delta1;

    /**
     * Define some handy constants.
     */
    protected final double t1n;
    protected final int dispersivePower1;
    protected final int repDispPower;
    protected final int repDispPower1;

    /**
     * Define scale factors between 1-2, 1-3, etc. atoms.
     */
    protected double scale12 = 0.0;
    protected double scale13 = 0.0;
    protected double scale14 = 1.0;
    protected double scale15 = 1.0;

    /**
     * Maximum number of classes in the force field.
     */
    protected int maxClass;
    /**
     * Store combined radius and epsilon values.
     */
    protected double radEps[][];
    protected static final byte RADMIN = 0;
    protected static final byte EPS = 1;

    public double getScale14() {
        return scale14;
    }

    private class VDWPowers {

        public double rhoDisp1(double rho) {
            return pow(rho, dispersivePower1);
        }

        public double rhoDelta1(double rhoDelta) {
            return pow(rhoDelta, repDispPower1);
        }
    }

    private class LJ_6_12 extends VDWPowers {

        @Override
        public double rhoDisp1(double rho) {
            double rho2 = rho * rho;
            return rho2 * rho2 * rho;
        }

        @Override
        public double rhoDelta1(double rhoDelta) {
            double rhoDelta2 = rhoDelta * rhoDelta;
            return rhoDelta2 * rhoDelta2 * rhoDelta;
        }
    }

    private class Buffered_14_7 extends VDWPowers {

        @Override
        public double rhoDisp1(double rho) {
            double rho2 = rho * rho;
            return rho2 * rho2 * rho2;
        }

        @Override
        public double rhoDelta1(double rhoDelta) {
            double rhoDelta2 = rhoDelta * rhoDelta;
            return rhoDelta2 * rhoDelta2 * rhoDelta2;
        }
    }
}
