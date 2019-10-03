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
package ffx.potential.extended;

import ffx.utilities.Constants;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * <p>ExtConstants class.</p>
 *
 * @author Stephen LuCore
 *
 * @since 1.0
 */
public class ExtConstants {

    /**
     * Constant <code>titratableHydrogenNames</code>
     */
    public static final List<String> titratableHydrogenNames
            = Arrays.asList("HH", "HG", "HE2", "HD1", "HE2", "HD2", "HZ3");
    /**
     * Constant <code>backboneNames</code>
     */
    public static final List<String> backboneNames = Arrays.asList("N", "CA", "C", "O", "HA", "H");
    /**
     * Constant <code>beta=1 / Boltzmann</code>
     */
    public static final double beta = 1 / Constants.R;
    /**
     * Random force conversion to kcal/mol/A; formerly randomForce.
     */
    public static final double forceToKcal = sqrt(4.184) / 10e9;
    /**
     * Random force conversion to (kcal/mol/A)^2; formerly randomForce2.
     */
    public static final double forceToKcalSquared = forceToKcal * forceToKcal;
    /**
     * Conversion from natural to base ten.
     */
    public static final double log10 = Math.log(10);
    /**
     * Propagation occurs on master thread only; otherwise use (multiple,
     * unshared) ThreadLocalRandoms.
     */
    public static final Random RNG = new Random();

}
