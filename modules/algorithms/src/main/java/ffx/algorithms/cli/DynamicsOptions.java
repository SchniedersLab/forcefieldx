//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.cli;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.Integrator;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;

import ffx.potential.cli.WriteoutOptions;
import picocli.CommandLine;
import picocli.CommandLine.Option;

import java.util.logging.Logger;

/**
 * Represents command line options for scripts that calculate thermodynamics.
 *
 * @author Michael J. Schnieders
 * @author Hernan V. Bernabe
 * @since 1.0
 */
public class DynamicsOptions {

    /**
     * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of
     * 2.0 is possible for the RESPA integrator.
     */
    @Option(names = {"-d", "--dt"}, paramLabel = "1.0",
            description = "Time discretization step in femtoseconds.")
    double dt = 1.0;

    /**
     * -b or --thermostat sets the desired thermostat: current choices are
     * Adiabatic, Berendsen, or Bussi.
     */
    @Option(names = {"-b", "--thermostat"}, paramLabel = "Bussi",
            description = "Thermostat: [Adiabatic / Berendsen / Bussi].")
    String thermostatString = "BUSSI";

    /**
     * -i or --integrator sets the desired integrator: current choices are
     * Beeman, RESPA, Stochastic (i.e. Langevin dynamics) or Verlet.
     */
    @Option(names = {"-i", "--integrator"}, paramLabel = "Verlet",
            description = "Integrator: [Beeman / Respa / Stochastic / Verlet].")
    String integratorString = "Verlet";

    /**
     * -r or --report sets the thermodynamics reporting frequency in picoseconds
     * (0.1 psec default).
     */
    @Option(names = {"-r", "--report"}, paramLabel = "0.25",
            description = "Interval to report thermodynamics (psec).")
    double report = 0.25;

    /**
     * -w or --write sets snapshot save frequency in picoseconds (1.0 psec
     * default).
     */
    @Option(names = {"-w", "--write"}, paramLabel = "10.0",
            description = "Interval to write out coordinates (psec).")
    double write = 10.0;

    /**
     * -t or --temperature sets the simulation temperature (Kelvin).
     */
    @Option(names = {"-t", "--temperature"}, paramLabel = "298.15",
            description = "Temperature (Kelvin).")
    double temp = 298.15;

    /**
     * -n or --steps sets the number of molecular dynamics steps (default is 1
     * nsec).
     */
    @Option(names = {"-n", "--numberOfSteps"}, paramLabel = "1000000",
            description = "Number of molecular dynamics steps.")
    long steps = 1000000;

    /**
     * -o or --optimize saves low-energy snapshots discovered (only for single topology simulations).
     */
    @CommandLine.Option(names = {"-o", "--optimize"},
            description = "Optimize and save low-energy snapshots.")
    private boolean optimize = false;

    /**
     * -k or --checkpoint sets the restart save frequency in picoseconds (1.0 psec default).
     */
    @CommandLine.Option(names = {"-k", "--checkpoint"}, paramLabel = "1.0",
            description = "Interval to write out restart files (.dyn, .his, etc).")
    private double checkpoint = 1.0;

    /**
     * --mdE or --molecularDynamicsEngine over-rides the default engine choice for
     * integrating the equations of motion
     */
    @CommandLine.Option(names = {"--mdE", "--molecularDynamicsEngine"}, paramLabel = "FFX",
            description = "Use FFX or OpenMM to integrate dynamics.")
    private String engineString = null;

    /**
     * -z or --trajSteps Number of steps for each OpenMM MD cycle.
     */
    @CommandLine.Option(names = {"-z", "--trajSteps"}, paramLabel = "100",
            description = "Number of steps per MD cycle (--mdE = OpenMM only).")
    int trajSteps = 100;

    private static final Logger logger = Logger.getLogger(DynamicsOptions.class.getName());

    /**
     * Thermostat.
     */
    public ThermostatEnum thermostat;

    /**
     * Integrator.
     */
    public IntegratorEnum integrator;

    private MolecularDynamics.DynamicsEngine engine = null;

    /**
     * Parse the thermostat and integrator.
     */
    public void init() {
        thermostat = Thermostat.parseThermostat(thermostatString);
        integrator = Integrator.parseIntegrator(integratorString);
        if (engineString != null) {
            try {
                engine = MolecularDynamics.DynamicsEngine.valueOf(engineString.toUpperCase());
            } catch (Exception ex) {
                logger.warning(String.format(" Could not parse %s as a valid dynamics engine! Defaulting to the Platform-recommended engine.", engineString));
                engine = null;
            }
        }
    }

    /**
     * Initialize a MolecularDynamics from the parsed options.
     *
     * @param potential      a {@link ffx.numerics.Potential} object.
     * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param sh             a {@link ffx.algorithms.AlgorithmListener} object.
     * @param writeout       a {@link WriteoutOptions} object.
     * @return a {@link MolecularDynamics} object.
     */
    public MolecularDynamics getDynamics(WriteoutOptions writeout,
                                         Potential potential,
                                         MolecularAssembly activeAssembly,
                                         AlgorithmListener sh) {
        MolecularDynamics molDyn;
        if (engine == null) {
            molDyn = MolecularDynamics.dynamicsFactory(activeAssembly, potential, activeAssembly.getProperties(), sh, thermostat, integrator);
        } else {
            molDyn = MolecularDynamics.dynamicsFactory(activeAssembly, potential, activeAssembly.getProperties(), sh, thermostat, integrator, engine);
        }
        molDyn.setFileType(writeout.getFileType());
        molDyn.setRestartFrequency(write);
        molDyn.setIntervalSteps(trajSteps);

        return molDyn;
    }

    /**
     * <p>Getter for the field <code>temp</code>.</p>
     *
     * @return a double.
     */
    public double getTemp() {
        return temp;
    }

    /**
     * <p>Getter for the field <code>dt</code>.</p>
     *
     * @return a double.
     */
    public double getDt() {
        return dt;
    }

    /**
     * <p>Getter for the field <code>report</code>.</p>
     *
     * @return a double.
     */
    public double getReport() {
        return report;
    }

    /**
     * <p>Getter for the field <code>checkpoint</code>.</p>
     *
     * @return a double.
     */
    public double getCheckpoint() {
        return checkpoint;
    }

    /**
     * <p>Getter for the field <code>optimize</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getOptimize() {
        return optimize;
    }
}
