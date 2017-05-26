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

    private int mdSteps;
    private double timeStep;
    private double printInterval;
    private double temperature;
    private final MolecularDynamics molecularDynamics;

    public MDMove(MolecularAssembly assembly, Potential potentialEnergy,
            CompositeConfiguration properties, AlgorithmListener listener,
            Thermostats requestedThermostat, Integrators requestedIntegrator) {

        molecularDynamics = new MolecularDynamics(assembly,
                potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);

        molecularDynamics.init(mdSteps, timeStep, printInterval, printInterval, temperature, true, null);
        molecularDynamics.setQuiet(true);
    }

    @Override
    public void move() {
        mdSteps = 50;
        timeStep = 0.5;
        printInterval = 0.025;
        temperature = 298.15;
        boolean initVelocities = true;
        molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 0.025, temperature, initVelocities, null);
    }

    public double getKineticEnergy() {
        return molecularDynamics.getKineticEnergy();
    }

    @Override
    public void revertMove() {
        molecularDynamics.revertState();
    }

}
