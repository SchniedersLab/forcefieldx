/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.algorithms;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.bonded.MolecularAssembly;

/**
 * The ReplicaExchange implements various replica exchange methods including
 * temperature, Hamiltonian and lambda.
 *
 * @author Sergio Urahata
 *
 * @since 1.0
 */
public class ReplicaExchange implements Terminatable {

    private boolean done = false;
    private boolean terminate = false;
    private boolean sampling = false;
    private static final Logger logger = Logger.getLogger(ReplicaExchange.class.getName());
    private int nReplicas;
    private MolecularAssembly[] molecularAssembly;
    private MolecularDynamics[] replicas;
    private double temperatures[];
    private Thread threads[];
    private final CompositeConfiguration properties;
    private AlgorithmListener algorithmListener;

    public ReplicaExchange(MolecularAssembly assembly[],
            CompositeConfiguration properties, AlgorithmListener listener) {
        if (assembly == null) {
            nReplicas = 0;
            molecularAssembly = null;
        } else {
            nReplicas = assembly.length;
            molecularAssembly = new MolecularAssembly[nReplicas];
            replicas = new MolecularDynamics[nReplicas];
            temperatures = new double[nReplicas];
            threads = new Thread[nReplicas];
            for (int i = 0; i < nReplicas; i++) {
                molecularAssembly[i] = assembly[i];
                replicas[i] = new MolecularDynamics(molecularAssembly[i], molecularAssembly[i].getPotentialEnergy(),
                        properties, listener, Thermostats.BUSSI);
            }
        }
        this.algorithmListener = listener;
        this.properties = properties;
    }

    public void sample(double temperatures[], int cycles) {

        done = false;
        terminate = false;

        this.temperatures = temperatures;

        for (int i = 0; i < cycles; i++) {
            /**
             * Check for termination request.
             */
            if (terminate) {
                done = true;
                break;
            }

            dynamic(1000, 1.0, 10000, 1000, true);
            logger.info(String.format("Applying exchange condition for cycle %d.", i));
            exchange();

        }

    }

    private void exchange() {
        for (int i = 0; i < nReplicas - 1; i++) {

            MolecularAssembly assemblyA = molecularAssembly[i];
            MolecularAssembly assemblyB = molecularAssembly[i + 1];

            double energyA = assemblyA.getPotentialEnergy().getTotal();
            double energyB = assemblyB.getPotentialEnergy().getTotal();

            /**
             * TODO: Apply Boltzmann criteria.
             */
            boolean boltzmann = false;

            /**
             * If the Boltzmann criteria is satisfied, do the switch.
             */
            if (boltzmann) {
                MolecularAssembly systemA = molecularAssembly[i];
                MolecularDynamics replicaA = replicas[i];

                molecularAssembly[i] = molecularAssembly[i + 1];
                replicas[i] = replicas[i + 1];

                molecularAssembly[i + 1] = systemA;
                replicas[i + 1] = replicaA;

                logger.info(String.format("Exchanged systems %d and %d.", i, i + 1));
            }
        }
    }

    /**
     * Blocking molecualr dynamics; when this returns each replica has completed
     * the requested number of steps.
     *
     * @param nSteps
     * @param timeStep
     * @param printInterval
     * @param saveInterval
     * @param temperature
     * @param initVelocities
     */
    private void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final boolean initVelocities) {

        if (sampling) {
            logger.warning("Programming error - the Replicas are already being sampled.");
            return;
        }

        sampling = true;

        /**
         * Initialize each MolecularDynamics instance and start them sampling.
         */
        for (int i = 0; i < nReplicas; i++) {
            replicas[i].init(nSteps, timeStep, printInterval, saveInterval, 
                    temperatures[i], initVelocities, null);
            threads[i] = new Thread(replicas[i]);
            threads[i].start();
        }

        /**
         * Block until all Replicas have finished sampling.
         */
        synchronized (this) {
            try {
                while (sampling) {
                    sampling = false;
                    for (int i = 0; i < nReplicas; i++) {
                        if (threads[i].isAlive()) {
                            sampling = true;
                        }
                    }
                    /**
                     * Wait if at least one of the Replicas is still running.
                     */
                    if (sampling) {
                        wait(100);
                    }
                }

            } catch (Exception e) {
                String message = "Exception during ReplicaExchange sampling";
                logger.log(Level.WARNING, message, e);
                sampling = false;
            }
        }
    }

    /**
     * This should be implemented as a blocking interrupt; when the method
     * returns the <code>Terminatable</code> algorithm has reached a clean
     * termination point. For example, between minimize or molecular dynamics
     * steps.
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
                }
            }
        }
    }
}
