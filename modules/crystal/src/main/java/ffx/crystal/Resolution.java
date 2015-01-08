/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.crystal;

import org.apache.commons.configuration.CompositeConfiguration;

import static org.apache.commons.math3.util.FastMath.abs;

/**
 * The Resolution class encapsulates the sampling limits and resolution limits
 * for a given crystal and/or data set.
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 */
public class Resolution {

    public final double sampling;
    public final double resolution;
    public final double inverseResSq;

    /**
     * <p>
     * Constructor for Resolution.</p>
     *
     * @param resolution a double.
     * @param sampling a double.
     */
    public Resolution(double resolution, double sampling) {
        this.resolution = resolution;
        this.inverseResSq = 1.0 / (resolution * resolution);
        this.sampling = sampling;
    }

    /**
     * <p>
     * Constructor for Resolution.</p>
     *
     * @param resolution a double.
     */
    public Resolution(double resolution) {
        this(resolution, 1.0 / 1.5);
    }

    /**
     * <p>
     * checkProperties</p>
     *
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @return a {@link ffx.crystal.Resolution} object.
     */
    public static Resolution checkProperties(CompositeConfiguration properties) {
        double resolution = properties.getDouble("resolution", -1.0);
        double sampling = properties.getDouble("sampling", 0.6);

        if (resolution < 0.0) {
            return null;
        }

        return new Resolution(resolution, sampling);
    }

    /**
     * <p>
     * resolutionLimit</p>
     *
     * @return a double.
     */
    public double resolutionLimit() {
        return resolution;
    }

    /**
     * <p>
     * inverseResSqLimit</p>
     *
     * @return a double.
     */
    public double inverseResSqLimit() {
        return inverseResSq;
    }

    /**
     * <p>
     * samplingLimit</p>
     *
     * @return a double.
     */
    public double samplingLimit() {
        return sampling;
    }

    /**
     * <p>
     * inResolutionRange</p>
     *
     * @param res a double.
     * @return a boolean.
     */
    public boolean inResolutionRange(double res) {
        if (abs(res - this.resolution) < 1e-8) {
            return true;
        } else {
            return res > this.resolution;
        }
    }

    /**
     * <p>
     * inInverseResSqRange</p>
     *
     * @param res a double.
     * @return a boolean.
     */
    public boolean inInverseResSqRange(double res) {
        if (abs(res - this.inverseResSq) < 1e-8) {
            return true;
        } else {
            return res < this.inverseResSq;
        }
    }
}
