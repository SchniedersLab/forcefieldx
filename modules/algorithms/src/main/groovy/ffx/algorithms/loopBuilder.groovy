
/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

// LOOP BUILDER

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.algorithms.Minimize;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.SimulatedAnnealing;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;
import ffx.algorithms.RotamerOptimization
import ffx.algorithms.RotamerOptimization.Direction;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3;
import ffx.potential.bonded.Residue.ResidueType;

import edu.rit.pj.Comm
import java.util.Scanner;

// Default convergence criteria.
double eps = 1.0;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Time step in femtoseconds.
double timeStep = 2.5;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100; 

// Frequency to write out restart information in picoseconds.
double restartInterval = 1.0;  

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 10000;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;

// Integrators [ BEEMAN, RESPA, STOCHASTIC ]
Integrators integrator = Integrators.RESPA;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// File type of coordinate snapshots to write out.
String fileType = "PDB";

// Value of Lambda.
double lambda = 0.0;

// Rotamer Optimization
boolean runRotamer = false;

// Simulated Annealing
boolean runSimulatedAnnealing = false;

// Molecular Dynamics to Generate Set
boolean runMD = false;

// OSRW
boolean runOSRW = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc loopBuilder [options] <filename1>');
cli.h(longOpt:'help', 'Print this help message.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.n(longOpt:'steps', args:1, argName:'10000', 'Number of molecular dynamics steps.');
cli.d(longOpt:'dt', args:1, argName:'2.5', 'Time discretization step (fsec).');
cli.r(longOpt:'report', args:1, argName:'0.01', 'Interval to report thermodyanamics (psec).');
cli.w(longOpt:'write', args:1, argName:'100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.osrw(longOpt:'OSRW', 'Run OSRW.');
cli.sa(longOpt:'simulated annealing', 'Run simulated annealing.');
cli.md(longOpt:'molecular dynamics', args:1, argName:'10','MD generateration of sets (size).');
cli.rot(longOpt:'rotamer', 'Run rotamer optimization.');
cli.a(longOpt:'all', 'Run optimal pipeline of algorithms.');


def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

// Report interval in picoseconds.
if (options.r) {
    printInterval = Double.parseDouble(options.r);
}

// Write interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);
}

// Run OSRW
if (options.osrw){
    runOSRW = true;
}

// Run Simulated Annealing
if (options.sa){
    runSimulatedAnnealing = true;
}

// Run Rotamer Optimization
if (options.rot){
    runRotamer = true;
}

// Run MD to obtain set of possible loops
if (options.md){
    runMD = true;
}

// Default
if (!(options.osrw && options.sa)){
    runOSRW = true;
}

// Robust Default
if (options.a){
    runOSRW = true;
    runRotamer = true;
}
System.setProperty("buildLoops", "true");
System.setProperty("vdwterm", "false");

List<String> arguments = options.arguments();
String filename = null;
MolecularAssembly[] systems = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    filename = arguments.get(0);
    systems = open(filename);
} else {
    return cli.usage();
}

// Get a reference to the first system's ForceFieldEnergy.
ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
// Set built atoms active/use flags to true (false for other atoms).
Atom[] atoms = active.getAtomArray();
for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setActive(true);
        ai.setUse(true);
    } else {
        ai.setActive(false);
        ai.setUse(true);
    }
}

logger.info("\n Running minimize on built atoms of " + active.getName());
logger.info(" RMS gradient convergence criteria: " + eps);

// Minimization without vdW.
e = minimize(eps);

energy();

if(runOSRW){
    // Run OSRW with a vdW potential.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "true");
    System.setProperty("polarization", "none");
    System.setProperty("intramolecular-softcore", "true");
    System.setProperty("intermolecular-softcore", "true");
    System.setProperty("lambdaterm", "true");
    System.setProperty("ligand-vapor-elec","false");
    System.setProperty("ewald-alpha","0.0");
    System.setProperty("lambda-bias-cutoff", "3");
    System.setProperty("bias-gaussian-mag", "0.01");
    System.setProperty("lambda-bin-width", "0.01");

    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        if (ai.getBuilt()) {
            ai.setApplyLambda(true);
        } else {
            ai.setApplyLambda(false);
        }
    }

    forceFieldEnergy= new ForceFieldEnergy(active);
    forceFieldEnergy.setLambda(lambda);

    energy();

    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    forceFieldEnergy.getCrystal().setSpecialPositionCutoff(0.0);
    // OSRW will be configured for a single topology.
    File structureFile = new File(FilenameUtils.normalize(filename));
    structureFile = new File(structureFile.getAbsolutePath());
    String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
    File histogramRestart = new File(baseFilename + ".his");
    File lambdaRestart = new File(baseFilename + ".lam");

    boolean asynchronous = false;
    boolean wellTempered = false;
    OSRW osrw =  new OSRW(forceFieldEnergy, forceFieldEnergy, lambdaRestart, histogramRestart, active.getProperties(),
        temperature, timeStep, printInterval, saveInterval, asynchronous, sh, wellTempered);
    osrw.setLambda(lambda);
    osrw.setLoopBuilding(true);
    // Create the MolecularDynamics instance.
    MolecularDynamics molDyn = new MolecularDynamics(active, osrw, active.getProperties(),
        null, thermostat, integrator);
    File dyn = null;
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
        fileType, restartInterval, dyn);

    logger.info("\n Obtaining low energy coordinates");
    double[] lowEnergyCoordinates = osrw.getLowEnergyLoop();
    System.out.println(lowEnergyCoordinates);
    logger.info(" Placing low energy coordinates");
    forceFieldEnergy.setCoordinates(lowEnergyCoordinates);
    energy();
}

if (runSimulatedAnnealing) {   
    // Minimize with vdW.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "false");
    energy = new ForceFieldEnergy(active);
    e = minimize(eps);

    // SA with vdW.
    logger.info("\n Running simulated annealing on " + active.getName());
    double[] heatUpTemperatures = [150,250,400,700,1000];
    // Number of molecular dynamics steps at each temperature.
    int steps = 267; //267 at 3
    // Time step in femtoseconds.
    timeStep = 3.0;
    // Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
    thermostat = Thermostats.BERENDSEN;
    // Integrators [ BEEMAN, RESPA, STOCHASTIC]
    integrator = Integrators.RESPA;

    SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, energy, active.getProperties(), null, thermostat, integrator);
    simulatedAnnealing.annealToTargetValues(heatUpTemperatures, steps, timeStep);

    double[] annealingTargetTemperatures = [1000, 800, 600, 500, 400, 300];
    steps = 800; //800 at 3
    simulatedAnnealing.annealToTargetValues(annealingTargetTemperatures,steps,timeStep);
}

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    ai.setUse(true);
    ai.setApplyLambda(false);
}

// Optimize with the full AMOEBA potential energy.
System.setProperty("vdwterm", "true");
System.setProperty("mpoleterm", "true");
System.setProperty("polarization", "direct");
System.setProperty("intramolecularSoftcore", "false");
System.setProperty("intermolecularSoftcore", "false");
System.setProperty("lambdaterm", "false");

forceFieldEnergy = new ForceFieldEnergy(active);
e = minimize(eps);

if (runMD){
    // Number of molecular dynamics steps
    nSteps = 500000; // nSteps * timeStep = time in femtoseconds

    // Temperature in degrees Kelvin.
    temperature = 300;
    
    // Time step in femtoseconds.
    timeStep = 1;
    
    // Reset velocities (ignored if a restart file is given)
    initVelocities = true;

    // Write interval in picoseconds.
    saveInterval = 1;   //set size = (nSteps * timeStep) / (1000 * saveInterval)
    
    // Interval to write out restart file (psec)
    double restartFrequency = 1000;

    MolecularDynamics molDyn = new MolecularDynamics(active, energy, active.getProperties(), null, thermostat, integrator);
    molDyn.setFileType(fileType);
    molDyn.setRestartFrequency(restartFrequency);
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
}

if (runRotamer){
    
    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        ai.setActive(true);
        ai.setUse(true);
    }

    energy = new ForceFieldEnergy(active);
    boolean threeBodyTerm = false;
    RotamerOptimization rotamerOptimization;  


    logger.info(String.format(" Rotomer Optimization"));
    rotamerOptimization = new RotamerOptimization(active, energy, null);

    rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
    RotamerLibrary.setUseOrigCoordsRotamer(true);

    //RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);

    Polymer[] polymers = systems[0].getPolymers();
    ArrayList<Residue> fullResidueList = polymers[0].getResidues();
    ArrayList<Residue> residuesToRO = new ArrayList<>();

    //Rotomer Optimization inclusion list building (grab built residues)
    for (int i = 0; i < fullResidueList.size(); i++) {
        Residue r = fullResidueList[i];
        if (r.getBackboneAtoms().get(0).getBuilt()) {
            residuesToRO.add(fullResidueList[i]);
        } 
    }
    //Rotomer Optimization inclusion list building (grab residues within 7A of the built loop)
    boolean expandList = true
    double expansionDistance = 7.0;

    if (expandList) {
        // Do a sliding-window rotamer optimization on loop window with a radius-inclusion criterion.
        RotamerLibrary.setUseOrigCoordsRotamer(true);

       // rotamerOptimization.setForcedResidues(resID, resID);
        rotamerOptimization.setWindowSize(1);
        rotamerOptimization.setDistanceCutoff(expansionDistance);

        startResID = residuesToRO.get(0).getResidueNumber();
        finalResID = residuesToRO.get(residuesToRO.size() - 1).getResidueNumber();

        rotamerOptimization.setForcedResidues(startResID,finalResID);
        rotamerOptimization.setResidues(startResID, finalResID);
    }

    rotamerOptimization.setResiduesIgnoreNull(residuesToRO);

    residuesToRO = rotamerOptimization.getResidues();
    RotamerLibrary.measureRotamers(residuesToRO, false);
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.SLIDING_WINDOW);
}

String ext = FilenameUtils.getExtension(filename);
filename = FilenameUtils.removeExtension(filename);
saveAsPDB(systems, new File(filename + ".pdb"));


