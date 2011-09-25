/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.math.util.MathUtils;

/**
 * The Resolution class encapsulates the sampling limits and resolution limits
 * for a given crystal and/or data set.
 *
 * @author fennt
 * @since 1.0
 * @version $Id: $
 */
public class Resolution {

    public final double sampling;
    public final double resolution;
    public final double invres;

    /**
     * <p>Constructor for Resolution.</p>
     *
     * @param resolution a double.
     * @param sampling a double.
     */
    public Resolution(double resolution, double sampling) {
        this.resolution = resolution;
        this.invres = 1.0 / (resolution * resolution);
        this.sampling = sampling;
    }

    /**
     * <p>Constructor for Resolution.</p>
     *
     * @param resolution a double.
     */
    public Resolution(double resolution) {
        this(resolution, 1.0 / 1.5);
    }

    /**
     * <p>checkProperties</p>
     *
     * @param properties a {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @return a {@link ffx.crystal.Resolution} object.
     */
    public static Resolution checkProperties(CompositeConfiguration properties) {
        double res = properties.getDouble("resolution", -1.0);
        double sampling = properties.getDouble("sampling", 1.0 / 1.5);

        if (res < 0.0) {
            return null;
        }

        return new Resolution(res, sampling);
    }

    /**
     * <p>res_limit</p>
     *
     * @return a double.
     */
    public double res_limit() {
        return resolution;
    }

    /**
     * <p>invressq_limit</p>
     *
     * @return a double.
     */
    public double invressq_limit() {
        return invres;
    }

    /**
     * <p>sampling_limit</p>
     *
     * @return a double.
     */
    public double sampling_limit() {
        return sampling;
    }

    /**
     * <p>inResolutionRange</p>
     *
     * @param res a double.
     * @return a boolean.
     */
    public boolean inResolutionRange(double res) {
        if (MathUtils.equals(res, this.resolution, 1e-8)) {
            return true;
        } else if (res > this.resolution) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * <p>inInvresolutionRange</p>
     *
     * @param res a double.
     * @return a boolean.
     */
    public boolean inInvresolutionRange(double res) {
        if (MathUtils.equals(res, this.invres, 1e-8)) {
            return true;
        } else if (res < this.invres) {
            return true;
        } else {
            return false;
        }
    }
}
