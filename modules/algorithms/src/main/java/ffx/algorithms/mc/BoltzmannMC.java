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
package ffx.algorithms.mc;

import ffx.potential.AssemblyState;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import org.apache.commons.math3.util.FastMath;

/**
 * The BoltzmannMC abstract class is a skeleton for Boltzmann-weighted Metropolis
 * Monte Carlo simulations.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public abstract class BoltzmannMC  implements MetropolisMC {
    public static final double BOLTZMANN = 0.0019872041; // In kcal/(mol*K)
    protected double temperature = 298.15; // Room temperature (also SATP).
    protected double kbTinv = -1.0 / (BOLTZMANN * temperature); // Constant factor for Monte Carlo moves (-1/kbT)
    protected boolean print = true;
    
    protected double e1 = 0.0;
    protected double e2 = 0.0;
    protected double eAdjust = 0.0;
    
    /**
     * Criterion for accept/reject a move; intended to be used mostly internally.
     * @param e1 Initial energy
     * @param e2 Final energy
     * @return If move accepted
     */
    @Override
    public boolean evaluateMove(double e1, double e2) {
        if (e2 <= e1) {
            return true;
        } else {
            // p(X) = exp(-U(X)/kb*T)
            double prob = FastMath.exp(kbTinv * (e2 - e1));
            
            assert (prob >= 0.0 && prob <= 1.0);
            /* Section used for testing purposes only.
            if (prob < 0.0 || prob > 1.0) {
                throw new ArithmeticException(String.format(" Invalid Monte Carlo "
                        + "result: probability %10.6f of accepting move from "
                        + "energy %10.6f to energy %10.6f", prob, e1, e2));
            }*/
            
            double trial = ThreadLocalRandom.current().nextDouble();
            return (trial <= prob);
        }
    }

    @Override
    public void setTemperature(double temp) {
        temperature = temp;
        kbTinv = -1.0 / (BOLTZMANN * temperature);
    }
    
    @Override
    public void setPrint(boolean print) {
        this.print = print;
    }
    
    @Override
    public double getE1() {
        return e1;
    }
    
    @Override
    public double getE2() {
        return e2;
    }
    
    @Override
    public double getEAdjust() {
        return eAdjust;
    }
}
