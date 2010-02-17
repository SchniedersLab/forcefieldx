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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.PotentialEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.XYZFilter;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class MolecularDynamics implements Terminatable, Runnable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    private final int n;
    private final int dof;
    private final double[] x;
    private final double[] xPrevious;
    private final double[] v;
    private final double[] a;
    private final double[] aPrevious;
    private final double[] grad;
    private final double[] mass;
    private double currentTemp;
    private double kinetic;
    private double potential;
    private double total;
    private double dt;
    private final MolecularAssembly molecularAssembly;
    private final PotentialEnergy potentialEnergy;
    private final CompositeConfiguration properties;
    private final Atom atoms[];
    private AlgorithmListener algorithmListener;
    private Thermostat thermostat;
    private File archive = null;
    private XYZFilter xyzFilter = null;
    private boolean done;
    private boolean terminate;
    private int nSteps = 1000;
    private int printFrequency = 100;
    private int saveFrequency = 1000;
    private double temperature = 300.0;
    private boolean initVelocities = true;

    public MolecularDynamics(MolecularAssembly assembly,
            CompositeConfiguration properties,
            AlgorithmListener listener, Thermostats requestedThermostat) {
        this.molecularAssembly = assembly;
        this.algorithmListener = listener;
        if (molecularAssembly.getPotentialEnergy() == null) {
            potentialEnergy = new PotentialEnergy(molecularAssembly);
            molecularAssembly.setPotential(potentialEnergy);
        } else {
            potentialEnergy = molecularAssembly.getPotentialEnergy();
        }
        this.properties = properties;
        ArrayList<Atom> atomList = molecularAssembly.getAtomList();
        n = atomList.size();
        atoms = atomList.toArray(new Atom[n]);
        Arrays.sort(atoms);

        mass = new double[n];
        dof = n * 3;
        x = new double[dof];
        xPrevious = new double[dof];
        v = new double[dof];
        a = new double[dof];
        aPrevious = new double[dof];
        grad = new double[dof];

        if (requestedThermostat != null) {
            switch (requestedThermostat) {
                case ISOTHERMAL:
                    thermostat = null;
                    break;
                case BERENDSEN:
                    double tau = properties.getDouble("tau-temperature", 0.2);
                    thermostat = new Berendsen(n, x, v, mass, 300.0, tau);
                    break;
                case BUSSI:
                default:
                    tau = properties.getDouble("tau-temperature", 0.2);
                    thermostat = new Bussi(n, x, v, mass, 300.0, tau);
            }
        } else {
            thermostat = null;
        }
    }

    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    public void setArchiveFile(File archive) {
        this.archive = archive;
    }

    public File getArchiveFile() {
        return archive;
    }

    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities) {

        /**
         * Return if already running.
         */
        if (!done) {
            logger.warning("Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
            return;
        }

        this.nSteps = nSteps;
        /**
         * Convert the time step from femtoseconds to picoseconds.
         */
        dt = timeStep * 1.0e-3;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 1;
        if (printInterval > dt) {
            printFrequency = (int) (printInterval / dt);
        }

        /**
         * Convert the save interval to a save frequency.
         */
        saveFrequency = -1;
        if (saveInterval > dt) {
            saveFrequency = (int) (saveInterval / dt);
            if (archive == null) {
                File file = molecularAssembly.getFile();
                String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                archive = new File(filename + ".arc");
                archive = XYZFilter.version(archive);
            }
            logger.info(" Snap shots will be written to " + archive.getAbsolutePath());
            if (xyzFilter == null) {
                xyzFilter = new XYZFilter(molecularAssembly);
            }
        }

        this.temperature = temperature;
        this.initVelocities = initVelocities;
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is done.
     *
     * @param nSteps
     * @param timeStep
     * @param printInterval
     * @param saveInterval
     * @param temperature
     * @param initVelocities
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities) {

        /**
         * Return if already running; Could happen if two threads call
         * dynamic on the same MolecularDynamics instance.
         */
        if (!done) {
            logger.warning("Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        logger.info(" Molecular dynamics starting up");
        if (thermostat != null) {
            logger.info(String.format(" Sampling the NVT ensemble via a %s thermostat", thermostat.name));
        } else {
            logger.info(String.format(" Sampling the NVE ensemble"));
        }

        init(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities);

        // Now we're committed!
        done = false;

        logger.info(String.format(" Number of steps:     %8d", nSteps));
        logger.info(String.format(" Time step:           %8.3f (fsec)", timeStep));
        logger.info(String.format(" Print interval:      %8.3f (psec)", printInterval));
        logger.info(String.format(" Save interval:       %8.3f (psec)", saveInterval));
        logger.info(String.format(" Target temperature:  %8.3f Kelvin", temperature));

        Thread dynamicThread = new Thread(this);
        dynamicThread.start();
        synchronized (this) {
            try {
                while (dynamicThread.isAlive()) {
                    wait(100);
                }
            } catch (Exception e) {
                String message = "Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
    }

    @Override
    public void run() {
        done = false;
        terminate = false;

        /**
         * Set the target temperature.
         */
        if (thermostat != null) {
            thermostat.setTargetTemperature(temperature);
        }

        /**
         * Initialize atomic positions and masses.
         */
        int j = 0;
        for (int i = 0; i < n; i++) {
            Atom atom = atoms[i];
            mass[i] = atom.getMass();
            double xyz[] = atom.getXYZ();
            x[j++] = xyz[0];
            x[j++] = xyz[1];
            x[j++] = xyz[2];
        }

        /**
         * Initialize velocities and compute the kinetic energy.
         */
        if (thermostat != null && initVelocities) {
            thermostat.maxwell();
            kinetic = thermostat.getKineticEnergy();
            currentTemp = thermostat.getCurrentTemperture();
        } else {
            for (int i = 0; i < dof; i++) {
                v[i] = 0.0;
            }
            kinetic = 0.0;
        }

        potentialEnergy.setOptimizationScaling(null);
        potential = potentialEnergy.energyAndGradient(x, grad);
        total = kinetic + potential;

        /**
         * Initialize current and previous accelerations.
         */
        int index = 0;
        for (int i = 0; i < n; i++) {
            double m = mass[i];
            for (j = 0; j < 3; j++, index++) {
                a[index] = -Thermostat.convert * grad[index] / m;
                aPrevious[index] = a[index];
            }
        }

        logger.info(String.format("\n   Step      Kinetic    Potential        Total     Temp     Time"));
        logger.info(String.format("       %13.4f%13.4f%13.4f %8.2f ", kinetic, potential, total, currentTemp));

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            beeman(dt);

            if (thermostat != null) {
                kinetic = thermostat.getKineticEnergy();
                currentTemp = thermostat.getCurrentTemperture();
            } else {
                kinetic = Double.MIN_VALUE;
                currentTemp = Double.MIN_VALUE;
            }

            total = kinetic + potential;
            if (step % printFrequency == 0) {
                time = System.nanoTime() - time;
                logger.info(String.format(" %6d%13.4f%13.4f%13.4f%9.2f%9.3f", step, kinetic, potential, total, currentTemp, time * 1.0e-9));
                time = System.nanoTime();
            }

            if (saveFrequency > 0 && step % saveFrequency == 0 && archive != null) {
                if (xyzFilter.writeFile(archive, true)) {
                    logger.info(String.format(" Appended snap shot to " + archive.getName()));
                } else {
                    logger.severe(String.format(" Appending snap shot to " + archive.getName() + " failed"));
                }
            }

            if (algorithmListener != null && step % printFrequency == 0) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }

            if (terminate) {
                logger.info(String.format("\n Terminating after %8d time steps\n", step));
                done = true;
                break;
            }
        }

        if (!terminate) {
            logger.info(String.format(" Completed %8d time steps\n", nSteps));
        }
    }

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

    /**
     * Integrate Newton's equations of motion using a Beeman multistep recursion
     * formula; the actual coefficients are Brooks' "Better Beeman" values.
     *
     * @since 1.0
     */
    private void beeman(final double dt) {
        final double dt_8 = 0.125 * dt;
        final double dt2_8 = dt * dt_8;

        if (thermostat != null) {
            thermostat.halfStep(dt);
        }

        /**
         * Store the current atom positions, then find new atom positions
         * and half-step velocities via Beeman recusion.
         */
        for (int i = 0; i < dof; i++) {
            xPrevious[i] = x[i];
            double temp = 5.0 * a[i] - aPrevious[i];
            x[i] += v[i] * dt + temp * dt2_8;
            v[i] += temp * dt_8;
        }

        /**
         * Compute the potential energy and gradients.
         */
        potential = potentialEnergy.energyAndGradient(x, grad);

        /**
         * Use Newton's second law to get the next acceleration and find
         * the full-step velocities using the Beeman recusion.
         */
        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            for (int j = 0; j < 3; j++, index++) {
                aPrevious[index] = a[index];
                a[index] = -Thermostat.convert * grad[index] / m;
                v[index] += (3.0 * a[index] + aPrevious[index]) * dt_8;
            }
        }

        if (thermostat != null) {
            thermostat.fullStep(dt);
        }
    }
}
