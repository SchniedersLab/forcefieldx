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

import picocli.CommandLine.Option;
import picocli.CommandLine.ITypeConverter;

import ffx.algorithms.thermostats.Thermostat;
import ffx.algorithms.thermostats.ThermostatEnum;

/**
 * Dynamics options shared by Dynamics scripts that use the Pico CLI.
 */
public class DynamicsOptions {

    /**
     * -d or --dt sets the timestep in femtoseconds (default of 1.0). A value of 2.0 is possible for the RESPA integrator.
     */
    @Option(names = {"-d", "--dt"}, paramLabel = "1.0", description = "Time discretization step in femtoseconds.")
    double dt = 1.0;

    /**
     * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
     */
    @Option(names = {"-b", "--thermostat"}, converter = ThermoConverter.class,
            paramLabel = "Berendsen", description = "Thermostat: [Adiabatic / Berendsen / Bussi].")
    ThermostatEnum tstat = ThermostatEnum.BERENDSEN;

    class ThermoConverter implements ITypeConverter<ThermostatEnum> {
        public ThermostatEnum convert(String value) {
            return Thermostat.parseThermostat(value);
        }
    }

    
}
