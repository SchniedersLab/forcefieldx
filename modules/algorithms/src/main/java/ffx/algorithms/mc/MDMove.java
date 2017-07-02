/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;

/**
 * Use MD as a coordinate based MC move.
 *
 * @author Mallory R. Tollefson
 */
public class MDMove implements MCMove {

    private static final Logger logger = Logger.getLogger(MDMove.class.getName());

    private int mdSteps = 50;
    private double timeStep = 1.0;
    private double printInterval = 0.05;
    private double temperature = 298.15;
    private boolean initVelocities = true;

    private final double saveInterval = 10000.0;
    private final MolecularDynamics molecularDynamics;

    public MDMove(MolecularAssembly assembly, Potential potentialEnergy,
            CompositeConfiguration properties, AlgorithmListener listener,
            Thermostats requestedThermostat, Integrators requestedIntegrator) {

        molecularDynamics = MolecularDynamics.dynamicsFactory(assembly,
                potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);

        /**
         * Ensure at least one interval is printed.
         */
        if (printInterval < mdSteps * timeStep) {
            printInterval = mdSteps * timeStep;
        }

        molecularDynamics.init(mdSteps, timeStep, printInterval, saveInterval, temperature, true, null);
        molecularDynamics.setQuiet(true);
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public void initVelocities(boolean initVelocities) {
        this.initVelocities = initVelocities;
    }

    public void setMDParameters(int mdSteps, double timeStep) {
        this.mdSteps = mdSteps;
        this.timeStep = timeStep;
        printInterval = mdSteps * timeStep;
    }

    @Override
    public void move() {
        molecularDynamics.dynamic(mdSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, null);
    }

    public double getKineticEnergy() {
        return molecularDynamics.getKineticEnergy();
    }

    @Override
    public void revertMove() {
        try {
            molecularDynamics.revertState();
        } catch (Exception ex) {
            logger.severe(" The MD state could not be reverted.");
        }
    }

}
