/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.cli;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.integrators.Integrator;
import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.thermostats.Thermostat;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import picocli.CommandLine.Option;

/**
 * Dynamics options shared by Dynamics scripts that use the Pico CLI.
 */
public class DynamicsOptions {

    /**
     * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of
     * 2.0 is possible for the RESPA integrator.
     */
    @Option(names = {"-d", "--dt"}, paramLabel = "1.0", description = "Time discretization step in femtoseconds.")
    double dt = 1.0;

    /**
     * -b or --thermostat sets the desired thermostat: current choices are
     * Adiabatic, Berendsen, or Bussi.
     */
    @Option(names = {"-b", "--thermostat"},
            paramLabel = "Berendsen", description = "Thermostat: [Adiabatic / Berendsen / Bussi].")
    String thermostatString = "BERENDSEN";

    /**
     * -i or --integrator sets the desired integrator: current choices are
     * Beeman, RESPA, Velocity Verlet, or Stochastic (AKA Langevin dynamics).
     */
    @Option(names = {"-i", "--integrator"},
            paramLabel = "Beeman", description = "Integrator: [Beeman / Respa / Stochastic / VelocityVerlet ]")
    String integratorString = "VELOCITYVERLET";

    /**
     * -r or --report sets the thermodynamics reporting frequency in picoseconds
     * (0.1 psec default).
     */
    @Option(names = {"-r", "--report"}, paramLabel = "0.25", description = "Interval to report thermodynamics (psec).")
    double report = 0.25;

    /**
     * -w or --write sets snapshot save frequency in picoseconds (1.0 psec
     * default).
     */
    @Option(names = {"-w", "--write"}, paramLabel = "10.0", description = "Interval to write out coordinates (psec).")
    double write = 10.0;

    /**
     * -t or --temperature sets the simulation temperature (Kelvin).
     */
    @Option(names = {"-t", "--temperature"}, paramLabel = "298.15", description = "Temperature (Kelvin)")
    double temp = 298.15;

    /**
     * -n or --steps sets the number of molecular dynamics steps (default is 1
     * nsec).
     */
    @Option(names = {"-n"}, paramLabel = "1000000", description = "Number of molecular dynamics steps")
    int steps = 1000000;

    /**
     * -p or --npt Specify use of a MC Barostat at the given pressure (default
     * 1.0 atm).
     */
    @Option(names = {"-p", "--npt"}, paramLabel = "0",
            description = "Specify use of a MC Barostat at the given pressure (default of 0 = disabled)")
    double pressure = 0;

    /**
     * -f or --file Choose the file type to write [PDB/XYZ].
     */
    @Option(names = {"-f", "--file"}, paramLabel = "XYZ", description = "Choose file type to write [PDB/XYZ].")
    String fileType = "XYZ";

    /**
     * Thermostat.
     */
    public ThermostatEnum thermostat;
    
    /**
     * Integrator.
     */
    public IntegratorEnum integrator; 
    
    public void init() {
        thermostat = Thermostat.parseThermostat(thermostatString);
        integrator = Integrator.parseIntegrator(integratorString);        
    }
    
    public MolecularDynamics getDynamics(Potential potential,
            MolecularAssembly activeAssembly,
            AlgorithmListener sh) {
                
        MolecularDynamics molDyn = new MolecularDynamics(activeAssembly, potential,
                activeAssembly.getProperties(), sh, thermostat, integrator);
        molDyn.setFileType(fileType);
        molDyn.setRestartFrequency(write);

        return molDyn;
    }

}
