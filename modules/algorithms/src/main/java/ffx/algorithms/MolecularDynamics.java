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
package ffx.algorithms;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import java.util.Random;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
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
    private final Potential potentialEnergy;
    private final CompositeConfiguration properties;
    private AlgorithmListener algorithmListener;
    private Thermostat thermostat;
    private File archiveFile = null;
    private File dynFile = null;
    private File pdbFile = null;
    private XYZFilter xyzFilter = null;
    private DYNFilter dynFilter = null;
    private PDBFilter pdbFilter = null;
    private boolean done;
    private boolean terminate;
    private int nSteps = 1000;
    private int printFrequency = 100;
    private int saveFrequency = 1000;
    private int removeCOMMotionFrequency = 100;
    private double temperature = 300.0;
    private boolean initVelocities = true;
    private boolean loadRestart = false;
    private boolean initialized = false;

    /**
     * <p>Constructor for MolecularDynamics.</p>
     *
     * @param assembly a {@link ffx.potential.bonded.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a {@link ffx.algorithms.Thermostat.Thermostats} object.
     */
    public MolecularDynamics(MolecularAssembly assembly,
                             Potential potentialEnergy,
                             CompositeConfiguration properties,
                             AlgorithmListener listener,
                             Thermostats requestedThermostat) {
        this.molecularAssembly = assembly;
        this.algorithmListener = listener;
        this.potentialEnergy = potentialEnergy;

        this.properties = properties;
        mass = potentialEnergy.getMass();
        dof = potentialEnergy.getNumberOfVariables();
        x = new double[dof];
        xPrevious = new double[dof];
        v = new double[dof];
        a = new double[dof];
        aPrevious = new double[dof];
        grad = new double[dof];

        if (requestedThermostat == null) {
            requestedThermostat = Thermostats.ADIABATIC;
        }

        switch (requestedThermostat) {
            case ADIABATIC:
            default:
                thermostat = new Adiabatic(dof, x, v, mass, potentialEnergy.getVariableTypes());
                break;
            case BERENDSEN:
                double tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Berendsen(dof, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
                break;
            case BUSSI:
                tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Bussi(dof, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
        }

        if (properties.containsKey("randomseed")) {
            int randomSeed = properties.getInt("randomseed");
            thermostat.setRandomSeed(randomSeed);
        }

        done = true;
    }

    /**
     * <p>Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * <p>Setter for the field <code>x</code>.</p>
     *
     * @param x a double array to set the current parameters to.
     */
    public void setParameters(double x[]) {
        System.arraycopy(x, 0, this.x, 0, dof);
    }

    /**
     * <p>Getter for the field <code>x</code>.</p>
     *
     * @return a double array with the current parameters
     */
    public double[] getParameters() {
        return x;
    }

    /**
     * <p>Setter for the field <code>archiveFile</code>.</p>
     *
     * @param archive a {@link java.io.File} object.
     */
    public void setArchiveFile(File archive) {
        this.archiveFile = archive;
    }

    /**
     * <p>Getter for the field <code>archiveFile</code>.</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getArchiveFile() {
        return archiveFile;
    }

    /**
     * <p>init</p>
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
                     final double saveInterval, final double temperature, final boolean initVelocities,
                     final File dyn) {

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
            File file = molecularAssembly.getFile();
            String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
            if (archiveFile == null) {
                archiveFile = new File(filename + ".arc");
                archiveFile = XYZFilter.version(archiveFile);
            }
            logger.info(" Snap shots will be written to " + archiveFile.getName());

            if (dyn == null) {
                this.dynFile = new File(filename + ".dyn");
                loadRestart = false;
            } else {
                this.dynFile = dyn;
                loadRestart = true;
            }

            logger.info(" Restart file will be written to " + this.dynFile.getName());

            if (xyzFilter == null) {
                xyzFilter = new XYZFilter(file, molecularAssembly,
                                          molecularAssembly.getForceField(), properties);
            }

            if (dynFilter == null) {
                dynFilter = new DYNFilter(molecularAssembly);
            }

            if (pdbFilter == null) {
                pdbFile = new File(filename + "_dyn.pdb");
                pdbFilter = new PDBFilter(new File(filename + "_dyn.pdb"), molecularAssembly, null, null);
            }
        }

        this.temperature = temperature;
        this.initVelocities = initVelocities;


    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is done.
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
                        final double saveInterval, final double temperature, final boolean initVelocities,
                        final File dyn) {

        init(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);

        /**
         * Return if already running; Could happen if two threads call
         * dynamic on the same MolecularDynamics instance.
         */
        if (!done) {
            logger.warning("Programming error - a thread invoked dynamic when it was already running.");
            return;
        }
        done = false;

        logger.info("\n Molecular dynamics starting up");
        if (!(thermostat instanceof Adiabatic)) {
            logger.info(format(" Sampling the NVT ensemble using a %s", thermostat.toString()));
        } else {
            logger.info(format(" Sampling the NVE ensemble"));
        }

        if (dyn != null) {
            logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
        }

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
                String message = " Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
    }

    public void setRemoveCOMMotionFrequency(int removeCOMMotionFrequency) {
        this.removeCOMMotionFrequency = removeCOMMotionFrequency;
    }

    /** {@inheritDoc} */
    @Override
    public void run() {
        done = false;
        terminate = false;

        /**
         * Set the target temperature.
         */
        if (!(thermostat instanceof Adiabatic)) {
            thermostat.setTargetTemperature(temperature);
        }

        if (initialized) {
            /**
             * if we've already been here,
             * don't update coordinates or velocities or try
             * to read restart file.
             */
            if (initVelocities) {
                thermostat.maxwell();
            }
        } else {
            if (!loadRestart) {
                /**
                 * Initialize atomic coordinates.
                 */
                potentialEnergy.getCoordinates(x);
                /**
                 * Initialize atomic velocities.
                 */
                if (initVelocities) {
                    thermostat.maxwell();
                } else {
                    for (int i = 0; i < dof; i++) {
                        v[i] = 0.0;
                    }
                }
            } else {
                if (!dynFilter.readDYN(dynFile, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    return;
                }
            }
        }

        /**
         * Compute the current potential energy.
         */
        potentialEnergy.setScaling(null);
        potential = potentialEnergy.energyAndGradient(x, grad);

        /**
         * Compute the current kinetic energy.
         */
        thermostat.kineticEnergy();
        kinetic = thermostat.getKineticEnergy();
        currentTemp = thermostat.getCurrentTemperture();
        total = kinetic + potential;

        /**
         * Initialize current and previous accelerations.
         */
        if (!initialized) {
            if (!loadRestart) {
                for (int i = 0; i < dof; i++) {
                    a[i] = -Thermostat.convert * grad[i] / mass[i];
                    aPrevious[i] = a[i];
                }
            }
            initialized = true;
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

            if (step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            kinetic = thermostat.getKineticEnergy();
            currentTemp = thermostat.getCurrentTemperture();
            total = kinetic + potential;
            if (step % printFrequency == 0) {
                time = System.nanoTime() - time;
                logger.info(String.format(" %6d%13.4f%13.4f%13.4f%9.2f%9.3f", step, kinetic, potential, total, currentTemp, time * 1.0e-9));
                time = System.nanoTime();
            }

            if (saveFrequency > 0 && step % saveFrequency == 0 && archiveFile != null) {
                if (xyzFilter.writeFile(archiveFile, true)) {
                    logger.info(String.format(" Appended snap shot to " + archiveFile.getName()));
                } else {
                    logger.warning(String.format(" Appending snap shot to " + archiveFile.getName() + " failed"));
                }
                if (dynFilter.writeDYN(dynFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + dynFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + dynFile.getName() + " failed"));
                }
                if (pdbFilter.writeFile(pdbFile, false)) {
                    logger.info(String.format(" Wrote PDB file to " + pdbFile.getName()));
                }
            }

            if (algorithmListener != null && step % printFrequency == 0) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }

            if (terminate) {
                logger.info(String.format("\n Terminating after %8d time steps\n", step));
                break;
            }
        }

        if (!terminate) {
            logger.info(String.format(" Completed %8d time steps\n", nSteps));
        }

        done = true;
        terminate = false;
    }

    /** {@inheritDoc} */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating dynamics.\n", e);
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
        /**
         * Do half-step thermostat operation.
         */
        thermostat.halfStep(dt);
        /**
         * Store the current atom positions, then find new atom positions
         * and half-step velocities via Beeman recursion.
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
        for (int i = 0; i < dof; i++) {
            aPrevious[i] = a[i];
            a[i] = -Thermostat.convert * grad[i] / mass[i];
            v[i] += (3.0 * a[i] + aPrevious[i]) * dt_8;
        }

        /**
         * Compute the full-step kinetic energy.
         */
        thermostat.kineticEnergy();
        /**
         * Do full-step thermostat operation.
         */
        thermostat.fullStep(dt);
    }
    private double friction = 1.0;
    private Random random = null;
    private double pfric[] = null;
    private double vfric[] = null;
    private double afric[] = null;
    private double prand[] = null;
    private double vrand[] = null;

    private void stochastic(final double dt) {
        /**
         * Set the frictional and random coefficients.
         */
        for (int i = 0; i < dof; i++) {
            double m = mass[i];
            double fdt = friction * dt;
            /**
             * In the limit of no friction, SD recovers normal
             * molecular dynamics.
             */
            if (fdt <= 0.0) {
                pfric[i] = 1.0;
                vfric[i] = dt;
                afric[i] = 0.5 * dt * dt;
                prand[i] = 0.0;
                vrand[i] = 0.0;
            } else {
                double pterm = 0.0;
                double vterm = 0.0;
                double rho = 0.0;
                if (fdt > 0.05) {
                    /**
                     * Analytical expressions when friction coefficient is large
                     */
                    double efdt = Math.exp(-fdt);
                    pfric[i] = efdt;
                    vfric[i] = (1.0 - efdt) / friction;
                    afric[i] = (dt - vfric[i]) / friction;
                    pterm = 2.0 * fdt - 3.0 + (4.0 - efdt) * efdt;
                    vterm = 1.0 - efdt * efdt;
                    rho = (1.0 - efdt) * (1.0 - efdt) / Math.sqrt(pterm * vterm);
                } else {
                    /**
                     * Use a series expansions when friction coefficient is small.
                     */
                    double fdt2 = fdt * fdt;
                    double fdt3 = fdt * fdt2;
                    double fdt4 = fdt2 * fdt2;
                    double fdt5 = fdt2 * fdt3;
                    double fdt6 = fdt3 * fdt3;
                    double fdt7 = fdt3 * fdt4;
                    double fdt8 = fdt4 * fdt4;
                    double fdt9 = fdt4 * fdt5;
                    afric[i] = (fdt2 / 2.0 - fdt3 / 6.0 + fdt4 / 24.0
                                - fdt5 / 120.0 + fdt6 / 720.0
                                - fdt7 / 5040.0 + fdt8 / 40320.0
                                - fdt9 / 362880.0) / (friction * friction);
                    vfric[i] = dt - friction * afric[i];
                    pfric[i] = 1.0 - friction * vfric[i];
                    pterm = 2.0 * fdt3 / 3.0 - fdt4 / 2.0
                            + 7.0 * fdt5 / 30.0 - fdt6 / 12.0
                            + 31.0 * fdt7 / 1260.0 - fdt8 / 160.0
                            + 127.0 * fdt9 / 90720.0;
                    vterm = 2.0 * fdt - 2.0 * fdt2 + 4.0 * fdt3 / 3.0
                            - 2.0 * fdt4 / 3.0 + 4.0 * fdt5 / 15.0
                            - 4.0 * fdt6 / 45.0 + 8.0 * fdt7 / 315.0
                            - 2.0 * fdt8 / 315.0 + 4.0 * fdt9 / 2835.0;
                    rho = Math.sqrt(3.0) * (0.5 - fdt / 16.0
                                            - 17.0 * fdt2 / 1280.0
                                            + 17.0 * fdt3 / 6144.0
                                            + 40967.0 * fdt4 / 34406400.0
                                            - 57203.0 * fdt5 / 275251200.0
                                            - 1429487.0 * fdt6 / 13212057600.0
                                            + 1877509.0 * fdt7 / 105696460800.0);
                }
                /**
                 * Compute random terms to thermostat the nonzero friction case.
                 */
                double kB = 0;
                double kelvin = 0;
                double ktm = kB * kelvin / mass[i];
                double psig = Math.sqrt(ktm * pterm) / friction;
                double vsig = Math.sqrt(ktm * vterm);
                double rhoc = Math.sqrt(1.0 - rho * rho);
                double pnorm = random.nextGaussian();
                double vnorm = random.nextGaussian();
                prand[i] = psig * pnorm;
                vrand[i] = vsig * (rho * pnorm + rhoc * vnorm);
            }
        }
    }
}
