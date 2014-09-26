/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.algorithms;

import java.io.File;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.PotentialsUtils;
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
     * was called. May be more elegant to replace this by using protected variables
     * and simply inheriting the time() function.
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
     * @param assembly
     * @param nStep
     * @param timeStep
     * @param printInterval
     * @param saveInterval
     * @param temperature
     * @param initVelocities
     * @param dyn
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
                    Thermostat.Thermostats.BUSSI, Integrator.Integrators.BEEMAN);
            molecularDynamics.dynamic(nStep, timeStep, printInterval, saveInterval,
                    temperature, initVelocities, dyn);
        }
    }

    /**
     * Minimizes a MolecularAssembly using AMOEBA potential energy.
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
