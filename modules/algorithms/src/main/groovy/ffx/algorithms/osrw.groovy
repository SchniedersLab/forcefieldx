/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

// ORTHOGONAL SPACE RANDOM WALK

// Apache Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;
import groovy.transform.Field;

// Paralle Java Imports
import edu.rit.pj.Comm;

// Java Imports
import java.util.regex.Pattern;

// Force Field X Imports
//import ffx.algorithms.Barostat;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// Asychronous communication between walkers.
boolean asynchronous = false;

/* 
 * The @Field annotations are ncessary for variables to be used by the readFile
 * method. Without @Field, variables with explicit declaration are local to the
 * implied main method of the script, whereas @Field transforms them to a private
 * script variable... though Groovy kind of ignores the "private" access modifier.
 */

// List of topologies (MolecularAssemblies).
@Field def topologies = [];
// List of ForceFieldEnergies.
@Field def energies = [];

// First atom of the ligand.
@Field int ligandStart = 1;

// Last atom of the ligand.
@Field int ligandStop = -1;

// First atom of the ligand for the 2nd topology.
@Field int ligandStart2 = 1;

// Last atom of the ligand for the 2nd topology.
@Field int ligandStop2 = -1;

// Additional ligand atoms; intended for creating disjoint sets.
@Field List<Atom> ligandAtoms1 = new ArrayList<>();
@Field List<Atom> ligandAtoms2 = new ArrayList<>();
@Field def ranges1 = []; // Groovy mechanism for creating an untyped ArrayList.
@Field def ranges2 = [];
@Field def rangesA = [];
@Field def rangesB = [];

// First atom for no electrostatics.
@Field int noElecStart = 1;

// Last atom for no electrostatics.
@Field int noElecStop = -1;

// First atom of the 2nd topology for no electrostatics.
@Field int noElecStart2 = 1;

// Last atom of the 2nd topology for no electrostatics.
@Field int noElecStop2 = -1;

// Initial lambda value (0 is ligand in vacuum; 1 is ligand in PBC).
double lambda = 0.0;

// Friction coefficient on the Lambda particle.
double friction = 1.0e-18;

// Mass of the Lambda particle.
double mass = 1.0e-18;

// Interval between OSRW recursion kernel counts.
int count = 10;

// Magnitude of each Gaussian bias (kcal/mole).
double biasMag = 0.002;

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 100000000;

// Number of equilibration steps: default is 1 picosecond.
int eSteps = 1000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 1.0;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100.0;

// File type of coordinate snapshots to write out.
String fileType = "XYZ";

// Frequency to write out restart information in picoseconds.
double restartInterval = 0.1;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC ]
Integrators integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Reset OSRW statistics
boolean resetStatistics = false;

// Reset number of steps taken, ignoring any record in .lam file
boolean resetNumSteps = true;

// Constant pressure
boolean NPT = false;
int meanInterval = 10;
double minDensity = 0.5;
double maxDensity = 1.5;
double maxSideMove = 0.25;
double maxAngleMove = 1.0;

// Write traversal snapshots
boolean writeTraversals = false;

int numParallel = 1;
@Field int threadsAvail = edu.rit.pj.ParallelTeam.getDefaultThreadCount();
@Field int threadsPer = threadsAvail;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc osrw [options] <filename> [filename] [filename] [filename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'async', args:0, 'Walker communication is asynchronous.');
cli.n(longOpt:'steps', args:1, argName:'10000000', 'Number of molecular dynamics steps.');
cli.q(longOpt:'equilibrate', args:1, argName:'1000', 'Equilibration steps prior to OSRW counts.');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization step (fsec).');
cli.r(longOpt:'report', args:1, argName:'1.0', 'Interval to report thermodyanamics (psec).');
cli.w(longOpt:'write', args:1, argName:'100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.i(longOpt:'integrator', args:1, argName:'Beeman', 'Integrator: [Beeman / Respa / Stochastic]');
cli.s(longOpt:'start', args:1, argName:'1', 'Starting ligand atom.');
cli.s2(longOpt:'start2', args:1, argName:'1', 'Starting ligand atom for the 2nd topology.');
cli.f(longOpt:'final', args:1, argName:'-1', 'Final ligand atom.');
cli.f2(longOpt:'final2', args:1, argName:'-1', 'Final ligand atom for the 2nd topology.');
cli.la1(longOpt:'ligAtoms1', args:1, argName:'None', 'Period-separated ranges of 1st topology ligand atoms (e.g. 40-50.72-83)');
cli.la2(longOpt:'ligAtoms2', args:1, argName:'None', 'Period-separated ranges of 2nd topology ligand atoms (e.g. 40-50.72-83)');
cli.es1(longOpt:'noElecStart1', args:1, argName:'1', 'No Electrostatics Starting Atom.');
cli.es2(longOpt:'noElecStart2', args:1, argName:'1', 'No Electrostatics Starting Atom for the 2nd Topology.');
cli.ef1(longOpt:'noElecFinal1', args:1, argName:'-1', 'No Electrostatics Final Atom.');
cli.ef2(longOpt:'noElecfinal2', args:1, argName:'-1', 'No Electrostatics Final Atom for the 2nd topology.');
cli.l(longOpt:'lambda', args:1, argName:'0.0', 'Initial lambda value (> 1.0 distributes lambda across walkers)');
cli.c(longOpt:'count', args:1, argName:'10', 'Time steps between OSRW counts.');
cli.g(longOpt:'bias', args:1, argName:'0.002', 'Gaussian bias magnitude (kcal/mol).');
cli.m(longOpt:'mass', args:1, argName:'1e-18', 'Lambda particle mass.');
//cli.p(longOpt:'npt', args:0, 'Constant pressure MD (1 atm).');
cli.x(longOpt:'friction', args:1, argName:'1e-18', 'Lambda particle friction.');
cli.W(longOpt:'traversals', args:0, 'Write out lambda-traversal snapshots.');
//cli.ld(longOpt:'minDensity', args:1, argName:'0.5', 'Minimum density allowed by the barostat.');
//cli.hd(longOpt:'maxDensity', args:1, argName:'1.5', 'Maximum density allowed by the barostat.');
//cli.sm(longOpt:'maxSideMove', args:1, argName:'0.25', 'Maximum side move allowed by the barostat.');
//cli.am(longOpt:'maxAngleMove', args:1, argName:'1.0', 'Maximum angle move allowed by the barostat.');
//cli.mi(longOpt:'meanInterval', args:1, argName:'10', 'Mean number of MD steps between applications of the barostat.');
cli.rt(longOpt:'reset', args:0, 'Reset OSRW histogram once, when lambda reaches 0.99.');
cli.rn(longOpt:'resetNumSteps', args:1, argName:'true', 'Ignore prior steps logged in .lam files');
cli.uaA(longOpt:'unsharedAtomsA', args:1, argName:'None', 'Quad-Topology: Period-separated ranges of A dual-topology atoms not shared by B');
cli.uaB(longOpt:'unsharedAtomsB', args:1, argName:'None', 'Quad-Topology: Period-separated ranges of B dual-topology atoms not shared by A');
cli.np(longOpt:'numParallel', args:1, argName:'1', 'Number of topology energies to calculate in parallel');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Read in command line file.
String filename = arguments.get(0);

// Asynchronous?
if (options.a) {
    asynchronous = true;
}

// Constant Pressue?
if (options.p) {
    NPT = true;
}

// Reset OSRW statistics when L=0.99 is reached the first time.
if (options.rt) {
    resetStatistics = true;
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load the number of molecular dynamics steps.
if (options.q) {
    eSteps = Integer.parseInt(options.q);
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

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

// Starting ligand atom.
if (options.s) {
    ligandStart = Integer.parseInt(options.s);
}

// Final ligand atom.
if (options.f) {
    ligandStop = Integer.parseInt(options.f);
}

// First ligand atom from Topology 1 with no electrostatics.
if (options.es1) {
    noElecStart = Integer.parseInt(options.es1);
}

// Final ligand atom from Topology 1 with no electrostatics.
if (options.ef1) {
    noElecStop = Integer.parseInt(options.ef1);
}

// Starting ligand atom for Topology 2.
if (options.s2) {
    ligandStart2 = Integer.parseInt(options.s2);
}

// Final ligand atom for Topology 2.
if (options.f2) {
    ligandStop2 = Integer.parseInt(options.f2);
}

// First ligand atom from Topology 2 with no electrostatics.
if (options.es2) {
    noElecStart2 = Integer.parseInt(options.es2);
}

// Final ligand atom from Topology 2 with no electrostatics.
if (options.ef2) {
    noElecStop2 = Integer.parseInt(options.ef2);
}

// Starting lambda value.
if (options.l) {
    lambda = Double.parseDouble(options.l);
}

// Number of time steps between counts.
if (options.c) {
    count = Integer.parseInt(options.c);
}

// Gaussian bias magnitude (kcal/mol).
if (options.g) {
    biasMag = Double.parseDouble(options.g);
}

// Lambda particle mass.
if (options.m) {
    mass = Double.parseDouble(options.m);
}

// Lambda particle friction.
if (options.x) {
    friction = Double.parseDouble(options.x);
}

// Minimum density
if (options.ld) {
    minDensity = Double.parseDouble(options.ld);
}
// Maximum density
if (options.hd) {
    maxDensity = Double.parseDouble(options.hd);
}
// Max side move
if (options.sm) {
    maxSideMove = Double.parseDouble(options.sm);
}
// Max angle move
if (options.am) {
    maxAngleMove = Double.parseDouble(options.am);
}
// Max angle move
if (options.mi) {
    meanInterval = Integer.parseInt(options.mi);
}
// Traversal snapshots
if (options.W) {
    writeTraversals = true;
}

if (options.np) {
    numParallel = Integer.parseInt(options.np);
    if (threadsAvail % numParallel != 0) {
        logger.warning(String.format(" Number of threads available %d not evenly divisible by np %d; reverting to sequential", threadsAvail, numParallel));
        numParallel = 1;
    } else if (arguments.size() % numParallel != 0) {
        logger.warning(String.format(" Number of topologies %d not evenly divisible by np %d; reverting to sequential", arguments.size(), numParallel));
        numParallel = 1;
    } else {
        threadsPer = threadsAvail / numParallel;
    }
}

if (options.rn) {
    if (eSteps > 0) {
        println("");
        logger.warning(" Ignoring resetNumSteps input due to equilibration");
    } else if (options.rn.equalsIgnoreCase("false")) {
        resetNumSteps = false;
    }
}

if (options.la1) {
    ranges1 = options.la1.tokenize(".");
}
if (options.la2) {
    ranges2 = options.la2.tokenize(".");
}

println("\n Running Orthogonal Space Random Walk on " + filename);

File structureFile = new File(FilenameUtils.normalize(filename));
structureFile = new File(structureFile.getAbsolutePath());
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File histogramRestart = new File(baseFilename + ".his");
File lambdaOneFile = null;
File lambdaZeroFile = null;
File lambdaRestart = null;
File dyn = null;

@Field Comm world = Comm.world();
@Field int size = world.size();
int rank = 0;

// For a multi-process job, try to get the restart files from rank sub-directories.
if (size > 1) {
    rank = world.rank();
    File rankDirectory = new File(structureFile.getParent() + File.separator
        + Integer.toString(rank));
    if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
    }
    lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
    dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
    lambdaOneFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam1");
    lambdaZeroFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam0");
    structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
} else {
    // For a single process job, try to get the restart files from the current directory.
    lambdaRestart = new File(baseFilename + ".lam");
    dyn = new File(baseFilename + ".dyn");
    lambdaOneFile = new File(baseFilename + ".lam1");
    lambdaZeroFile = new File(baseFilename + ".lam0");
}

if (!dyn.exists()) {
    dyn = null;
}

// Turn on computation of lambda derivatives.
System.setProperty("lambdaterm", "true");

// Relative free energies via the DualTopologyEnergy class require different
// default OSRW parameters than absolute free energies.
if (arguments.size() >= 2) {
    // Condensed phase polarization is evaluated over the entire range.
    System.setProperty("polarization-lambda-start","0.0");
    // Polarization energy is not scaled individually by lambda, but
    // along with the overall potential energy of a topology.
    System.setProperty("polarization-lambda-exponent","0.0");
    // Ligand vapor electrostatics are not calculated. This cancels when the
    // difference between protein and water environments is considered.
    System.setProperty("ligand-vapor-elec","false");
    // Condensed phase polarization, without the ligand present, is unecessary.
    System.setProperty("no-ligand-condensed-scf","false");
}

@Field Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");

/**
 * Handles opening a file (filename), with 0-indexed number topNum.
 */
private void openFile(String toOpen, File structFile, int topNum) {
    open(toOpen, threadsPer);
    if (size > 1) {
        active.setFile(structFile);
    }
    ForceFieldEnergy energy = active.getPotentialEnergy();
    Atom[] atoms = active.getAtomArray();
    int remainder = (topNum % 2) + 1;
    switch(remainder) {
        case 1:
            for (int i = ligandStart; i <= ligandStop; i++) {
                Atom ai = atoms[i-1];
                ai.setApplyLambda(true);
                ai.print();
            }
            if (ranges1) {
                for (range in ranges1) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        // Don't need to worry about negative numbers; rangeregex just won't match.
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            Atom ai = atoms[i-1];
                            ai.setApplyLambda(true);
                            ai.print();
                        }
                    } else {
                        logger.warning(" Could not recognize ${range} as a valid range; skipping");
                    }
                }
            }

            // Apply the no electrostatics atom selection
            if (noElecStart < 1) {
                noElecStart = 1;
            }
            if (noElecStop > atoms.length) {
                noElecStop = atoms.length;
            }
            for (int i = noElecStart; i <= noElecStop; i++) {
                Atom ai = atoms[i - 1];
                ai.setElectrostatics(false);
                ai.print();
            }
            break;
        case 2:
            for (int i = ligandStart2; i <= ligandStop2; i++) {
                Atom ai = atoms[i-1];
                ai.setApplyLambda(true);
                ai.print();
            }
            if (ranges2) {
                for (range in ranges2) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        // Don't need to worry about negative numbers; rangeregex just won't match.
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            Atom ai = atoms[i-1];
                            ai.setApplyLambda(true);
                            ai.print();
                        }
                    } else {
                        logger.warning(" Could not recognize ${range} as a valid range; skipping");
                    }
                }
            }

            // Apply the no electrostatics atom selection
            if (noElecStart2 < 1) {
                noElecStart2 = 1;
            }
            if (noElecStop2 > atoms.length) {
                noElecStop2 = atoms.length;
            }
            for (int i = noElecStart2; i <= noElecStop2; i++) {
                Atom ai = atoms[i - 1];
                ai.setElectrostatics(false);
                ai.print();
            }
            break;
    }
    
    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff(0.0);
    // Save a reference to the topology.
    topologies[topNum] = active;
    energies[topNum] = energy;
}

openFile(filename, structureFile, 0);

/*
// Open the first system
open(filename);
// If this is a multi-process job, set the structure file to come from the subdirectory.
if (size > 1) {
    active.setFile(structureFile);
}

// Get a reference to the first system's ForceFieldEnergy and atom array.
ForceFieldEnergy energy = active.getPotentialEnergy();
Atom[] atoms = active.getAtomArray();
// Apply the ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
//    ai.setElectrostatics(!noElec1);
    ai.print();
}

if (ranges1) {
    for (range in ranges1) {
        def m = rangeregex.matcher(range);
        if (m.find()) {
            if (m.groupCount() > 1) {
                int rangeStart = Integer.parseInt(m.group(1));
                int rangeEnd = Integer.parseInt(m.group(2));
                if (rangeStart > rangeEnd) {
                    logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                }
                // Don't need to worry about negative numbers; rangeregex just won't match.
                for (int i = rangeStart; i <= rangeEnd; i++) {
                    Atom ai = atoms[i-1];
                    ai.setApplyLambda(true);
                    ai.print();
                }
            } else {
                int i = Integer.parseInt(m.group(1));
                Atom ai = atoms[i-1];
                ai.setApplyLambda(true);
                ai.print();
            }
        } else {
            logger.warning(" Could not recognize ${range} as a valid range; skipping");
        }
    }
}

// Apply the no electrostatics atom selection
if (noElecStart < 1) {
    noElecStart = 1;
}
if (noElecStop > atoms.length) {
    noElecStop = atoms.length;
}
for (int i = noElecStart; i <= noElecStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setElectrostatics(false);
    ai.print();
}

// Turn off checks for overlapping atoms, which is expected for lambda=0.
energy.getCrystal().setSpecialPositionCutoff(0.0);
*/

// OSRW will be configured for the appropriate number of topologies.
OSRW osrw = null;
// Save a reference to the first topology.
//topologies[0] = active;

if (arguments.size() == 1) {
    // Check for constant pressure
    if (NPT) {
        //        // Create a barostat.
        //        Barostat barostat = new Barostat(active);
        //        barostat.setMaxdUdL(1000.0);
        //        barostat.setMaxDensity(maxDensity);
        //        barostat.setMinDensity(minDensity);
        //        barostat.setMaxSideMove(maxSideMove);
        //        barostat.setMaxAngleMove(maxAngleMove);
        //        barostat.setMeanBarostatInterval(meanInterval);
        //
        //        // Create the OSRW instance.
        //        osrw = new OSRW(energy, barostat, lambdaRestart, histogramRestart, active.getProperties(),
        //            temperature, timeStep, printInterval, saveInterval, asynchronous, sh);
        //        osrw.setResetStatistics(resetStatistics);
        if (writeTraversals) {
            osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[0]);
        }
    } else {
        // Wrap the single topology ForceFieldEnergy inside an OSRW instance.
        osrw = new OSRW(energies[0], energies[0], lambdaRestart, histogramRestart, active.getProperties(),
            temperature, timeStep, printInterval, saveInterval, asynchronous, resetNumSteps, sh);
        osrw.setResetStatistics(resetStatistics);
        if (writeTraversals) {
            osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[0]);
        }
    }
} else if (arguments.size() == 2) {
    // Open the 2nd topology.
    filename = arguments.get(1);
    openFile(filename, structureFile, 1);
    
    /*
    // If this is a multi-process job, set the structure file to come from the subdirectory.
    if (size > 1) {
        active.setFile(structureFile);
    }
    energy = active.getPotentialEnergy();
    atoms = active.getAtomArray();
    // Apply the ligand atom selection for the 2nd topology.
    for (int i = ligandStart2; i <= ligandStop2; i++) {
        Atom ai = atoms[i - 1];
        ai.setApplyLambda(true);
        ai.print();
    }
    
    if (ranges2) {
        for (range in ranges2) {
            def m = rangeregex.matcher(range);
            if (m.find()) {
                if (m.groupCount() > 1) {
                    int rangeStart = Integer.parseInt(m.group(1));
                    int rangeEnd = Integer.parseInt(m.group(2));
                    if (rangeStart > rangeEnd) {
                        logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                    }
                    // Don't need to worry about negative numbers; rangeregex just won't match.
                    for (int i = rangeStart; i <= rangeEnd; i++) {
                        Atom ai = atoms[i-1];
                        ai.setApplyLambda(true);
                        ai.print();
                    }
                } else {
                    int i = Integer.parseInt(m.group(1));
                    Atom ai = atoms[i-1];
                    ai.setApplyLambda(true);
                    ai.print();
                }
            } else {
                logger.warning(" Could not recognize ${range} as a valid range; skipping");
            }
        }
    }

    // Apply the no electrostatics atom selection
    if (noElecStart2 < 1) {
        noElecStart2 = 1;
    }
    if (noElecStop2 > atoms.length) {
        noElecStop2 = atoms.length;
    }
    for (int i = noElecStart2; i <= noElecStop2; i++) {
        Atom ai = atoms[i - 1];
        ai.setElectrostatics(false);
        ai.print();
    }

    // Save a reference to the second topology.
    topologies[1] = active;
    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff(0.0);
    */
   
    // Create the DualTopology potential energy.
    DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topologies[0], topologies[1]);
    if (numParallel == 2) {
        dualTopologyEnergy.setParallel(true);
    }
    // Wrap the DualTopology potential energy inside an OSRW instance.
    osrw = new OSRW(dualTopologyEnergy, dualTopologyEnergy, lambdaRestart,
        histogramRestart, active.getProperties(), temperature, timeStep, printInterval,
        saveInterval, asynchronous, resetNumSteps, sh);
    osrw.setResetStatistics(resetStatistics);
    if (writeTraversals) {
        osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[1]);
    }
} else if (arguments.size() == 4) {
    logger.info(" For quad-topology calculations, systems should be set up as follows:");
    logger.info(" The first two define the first dual topology, the second two should define the second dual topology in the same order.");
    logger.info(" For example, for a dual force field correction on decharging sodium, topologies should be in this order:");
    logger.info(" Sodium-AMOEBA, sodium-AMBER, decharged Na-AMOEBA, decharged Na-AMBER");
    
    openFile(arguments.get(1), structureFile, 1);
    openFile(arguments.get(2), structureFile, 2);
    openFile(arguments.get(3), structureFile, 3);
    DualTopologyEnergy dtA = new DualTopologyEnergy(topologies[0], topologies[1]);
    // Intentionally reversed order.
    DualTopologyEnergy dtB = new DualTopologyEnergy(topologies[3], topologies[2]);
    List<Integer> uniqueA = new ArrayList<>();
    List<Integer> uniqueB = new ArrayList<>();
    
    if (options.uaA) {
        rangesA = options.uaA.tokenize(".");
        def ra = [] as Set;
        for (range in rangesA) {
            def m = rangeregex.matcher(range);
            if (m.find()) {
                int rangeStart = Integer.parseInt(m.group(1));
                int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                if (rangeStart > rangeEnd) {
                    logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
}
                for (int i = rangeStart; i <= rangeEnd; i++) {
                    ra.add(i-1);
                }
            }
        }
        Atom[] atA1 = topologies[0].getAtomArray();
        int counter = 0;
        def raAdj = [] as Set; // Indexed by common variables in dtA.
        for (int i = 0; i < atA1.length; i++) {
            Atom ai = atA1[i];
            if (i in ra) {
                if (ai.applyLambda()) {
                    logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                } else {
                    logger.info(String.format(" Unshared A: %d variables %d-%d", i, counter, counter+2));
                    for (int j = 0; j < 3; j++) {
                        raAdj.add(new Integer(counter + j));
                    }
                }
            }
            if (! ai.applyLambda()) {
                counter += 3;
            }
        }
        if (raAdj) {
            uniqueA.addAll(raAdj);
        }
    }
    if (options.uaB) {
        rangesB = options.uaB.tokenize(".");
        def rb = [] as Set;
        for (range in rangesB) {
            def m = rangeregex.matcher(range);
            if (m.find()) {
                int rangeStart = Integer.parseInt(m.group(1));
                int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                if (rangeStart > rangeEnd) {
                    logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                }
                for (int i = rangeStart; i <= rangeEnd; i++) {
                    rb.add(i-1);
                }
            }
        }
        Atom[] atB1 = topologies[2].getAtomArray();
        int counter = 0;
        def rbAdj = [] as Set; // Indexed by common variables in dtA.
        for (int i = 0; i < atB1.length; i++) {
            Atom bi = atB1[i];
            if (i in rb) {
                if (bi.applyLambda()) {
                    logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                } else {
                    logger.info(String.format(" Unshared B: %d variables %d-%d", i, counter, counter+2));
                    for (int j = 0; j < 3; j++) {
                        rbAdj.add(counter + j);
                    }
                }
            }
            if (! bi.applyLambda()) {
                counter += 3;
            }
        }
        if (rbAdj) {
            uniqueB.addAll(rbAdj);
        }
    }

    QuadTopologyEnergy qte = new QuadTopologyEnergy(dtA, dtB, uniqueA, uniqueB);
    if (numParallel >= 2) {
        qte.setParallel(true);
        if (numParallel == 2) {
            dtA.setParallel(true);
            dtB.setParallel(true);
        }
    }
    osrw = new OSRW(qte, qte, lambdaRestart, histogramRestart, 
        active.getProperties(), temperature, timeStep, printInterval, 
        saveInterval, asynchronous, resetNumSteps, sh);
} else {
    logger.severe(" Must have 1, 2, or 4 topologies to test.");
}

// Apply the command line lambda value if a lambda restart file does not exist.
if (!lambdaRestart.exists()) {
    if (lambda >= 0.0 && lambda <= 1.0) {
        osrw.setLambda(lambda);
        logger.info(String.format(" Setting lambda to %5.3f.", lambda));
    } else {
        if (size > 1) {
            dL = 1.0 / (size - 1.0);
            lambda = rank * dL;
            if (lambda > 1.0) {
                lambda = 1.0;
            }
            if (lambda < 0.0) {
                lambda = 0.0;
            }
            logger.info(String.format(" Setting lambda to %5.3f.", lambda));
            osrw.setLambda(lambda);
        } else {
            lambda = 0.5;
            logger.info(String.format(" Setting lambda to %5.3f", lambda));
            osrw.setLambda(lambda);
        }
    }
}

// Apply the command line OSRW values if a histogram restart file does not exist.
if (!histogramRestart.exists()) {
    osrw.setThetaFrication(friction);
    osrw.setThetaMass(mass);
    osrw.setCountInterval(count);
    osrw.setBiasMagnitude(biasMag);
}

// Create the MolecularDynamics instance.
MolecularDynamics molDyn = new MolecularDynamics(topologies[0], osrw, topologies[0].getProperties(), null, thermostat, integrator);

// Start sampling.
if (eSteps > 0) {
    logger.info(" Beginning equilibration");
    osrw.setPropagateLambda(false);
    molDyn.dynamic(eSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
    logger.info(" Beginning OSRW sampling");
    osrw.setPropagateLambda(true);
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, false,
        fileType, restartInterval, dyn);
} else {
    logger.info(" Beginning OSRW sampling without equilibration");
    if (!resetNumSteps) {
        int nEnergyCount = osrw.getEnergyCount();
        if (nEnergyCount > 0) {
            nSteps -= nEnergyCount;
            logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
        }
    }
    if (nSteps > 0) {
    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
        fileType, restartInterval, dyn);
    } else {
        logger.info(" No steps remaining for this process!");
    }
}

