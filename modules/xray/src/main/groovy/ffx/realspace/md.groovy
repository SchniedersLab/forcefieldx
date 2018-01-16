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

// REALSPACE MOLECULAR DYNAMICS

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.MolecularDynamics
import ffx.algorithms.thermostats.Thermostat
import ffx.algorithms.thermostats.ThermostatEnum
import ffx.realspace.RealSpaceData
import ffx.realspace.RealSpaceFile
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 100.0;

// ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
ThermostatEnum thermostat = ThermostatEnum.BERENDSEN;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc realspace.md [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:2, valueSeparator:',', argName:'data.map,1.0', 'specify input data filename (or simply provide the datafilename argument after the PDB file) and weight to apply to the data (wA)');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.f(longOpt:'dt', args:1, argName:'1.0', 'Time step in femtoseconds.');
cli.i(longOpt:'print', args:1, argName:'0.01', 'Interval to print out thermodyanamics in picoseconds.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates in picoseconds.');
cli.t(longOpt:'temperature', args:1, argName:'100.0', 'Temperature in degrees Kelvin.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic/Berendsen/Bussi]')
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Name of the file (PDB or XYZ).
String modelfilename = arguments.get(0);

// set up real space map data (can be multiple files)
List mapfiles = new ArrayList();
if (arguments.size() > 1) {
    RealSpaceFile realspacefile = new RealSpaceFile(arguments.get(1), 1.0);
    mapfiles.add(realspacefile);
}
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=2) {
	double wA = Double.parseDouble(options.ds[i+1]);
	RealSpaceFile realspacefile = new RealSpaceFile(options.ds[i], wA);
	mapfiles.add(realspacefile);
    }
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load the time steps in femtoseconds.
if (options.f) {
    timeStep = Double.parseDouble(options.f);
}

// Print interval in picoseconds.
if (options.i) {
    printInterval = Double.parseDouble(options.i);
}

// Write interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

if (options.b) {
    thermostat = Thermostat.parseThermostat(options.b)
}


logger.info("\n Running molecular dynmaics on " + modelfilename);
open(modelfilename);

if (mapfiles.size() == 0) {
    RealSpaceFile realspacefile = new RealSpaceFile(active, 1.0);
    mapfiles.add(realspacefile);
}

RealSpaceData realspacedata = new RealSpaceData(active,
    active.getProperties(), active.getParallelTeam(),
    mapfiles.toArray(new RealSpaceFile[mapfiles.size()]));

energy();

RefinementEnergy refinementEnergy = new RefinementEnergy(realspacedata, RefinementMode.COORDINATES);

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}
MolecularDynamics molDyn = new MolecularDynamics(active, refinementEnergy, active.getProperties(), refinementEnergy, thermostat, null);
refinementEnergy.setThermostat(molDyn.getThermostat());

molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
