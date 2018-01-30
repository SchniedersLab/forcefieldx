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
package ffx.algorithms;

import java.io.File;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.Keyword;

/**
 * The AlgorithmUtils class provides a local implementation, independent of the
 * User Interfaces module, of AlgorithmFunctions methods (MD, minimization).
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class AlgorithmUtils extends PotentialsUtils implements AlgorithmFunctions {

    private static final Logger logger = Logger.getLogger(AlgorithmUtils.class.getName());
    private final long initTime;
    private long interTime;

    public AlgorithmUtils() {
        initTime = System.nanoTime();
        interTime = initTime;
    }

    /**
     * Logs time since this interface was created and the last time this method
     * was called. May be more elegant to replace this by using protected
     * variables and simply inheriting the time() function.
     *
     * @return Time since last call (double).
     */
    @Override
    public double time() {
        long currTime = System.nanoTime();
        logger.info(String.format(" Time since interface established: %f", (currTime - initTime) * 1.0E-9));
        double elapsed = (currTime - interTime) * 1.0E-9;
        interTime = currTime;
        logger.info(String.format(" Time since last timer call: %f", elapsed));
        return elapsed;
    }

    /**
     * Performs molecular dynamics on a MolecularAssembly.
     *
     * @param assembly the MolecularAssembly to do MD on.
     * @param nStep the number of MD steps.
     * @param timeStep the MD time step.
     * @param printInterval the time between logging statements.
     * @param saveInterval the time between saving restart files.
     * @param temperature the target temperature.
     * @param initVelocities if true, reset velocities from a Maxwell
     * distribution.
     * @param dyn a dynamics restart.
     */
    @Override
    public void md(MolecularAssembly assembly, int nStep, double timeStep,
            double printInterval, double saveInterval, double temperature,
            boolean initVelocities, File dyn) {
        if (assembly == null) {
            logger.info(" No active system to minimize.");
        } else {
            CompositeConfiguration properties = Keyword.loadProperties(assembly.getFile());
            MolecularDynamics molecularDynamics = new MolecularDynamics(assembly,
                    assembly.getPotentialEnergy(), properties, null,
                    ThermostatEnum.BUSSI, IntegratorEnum.BEEMAN);
            molecularDynamics.dynamic(nStep, timeStep, printInterval, saveInterval,
                    temperature, initVelocities, dyn);
        }
    }

    /**
     * Minimizes a MolecularAssembly using AMOEBA potential energy.
     *
     * @param assembly Assembly to minimize
     * @param eps Convergence criterion
     * @return A Potential.
     */
    @Override
    public Potential minimize(MolecularAssembly assembly, double eps) {
        if (assembly == null) {
            logger.info(" No active system to minimize.");
            return null;
        } else {
            Minimize minimize = new Minimize(assembly, null);
            Potential potential = minimize.minimize(eps);
            return potential;
        }
    }
}
