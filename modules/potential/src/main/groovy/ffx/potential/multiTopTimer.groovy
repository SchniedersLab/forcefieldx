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

// TIMER

// Groovy Imports
import groovy.util.CliBuilder;
import groovy.transform.Field;

// Java imports
import java.util.regex.Pattern;

// FFX Imports
import ffx.potential.ForceFieldEnergy;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.Atom;
import ffx.numerics.Potential;

// The number of iterations.
int nEvals = 5;

// Compute the atomic coordinate gradient.
boolean gradient = true;

// Print the energy for each iteraction.
boolean print = true;

@Field def topologies = [];
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
double lambda = 0.5;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc timer [options] <filename> [filename] [filename] [filename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.n(longOpt:'iterations', args:1, argName:'5', 'Number of iterations');
cli.c(longOpt:'threads', args:1, argName:'all', 'Number of SMP threads (ie. default uses all CPU cores)');
cli.g(longOpt:'gradient', args:1, argName:'true', 'Compute the atomic coordinates gradient');
cli.v(longOpt:'verbose', args:1, argName:'true', 'Print out the energy for each step');
cli.l(longOpt:'lambda', args:1, argName:'0.5', 'Lambda value to test.');
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
cli.uaA(longOpt:'unsharedAtomsA', args:1, argName:'None', 'Quad-Topology: Period-separated ranges of A dual-topology atoms not shared by B');
cli.uaB(longOpt:'unsharedAtomsB', args:1, argName:'None', 'Quad-Topology: Period-separated ranges of B dual-topology atoms not shared by A');
cli.np(longOpt:'numParallel', args:1, argName:'1', 'Number of topology energies to calculate in parallel');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Load the number iterations.
if (options.n) {
    nEvals = Integer.parseInt(options.n);
}

@Field int threadsAvail = edu.rit.pj.ParallelTeam.getDefaultThreadCount();
// Load the number of threads.
if (options.c) {
    System.setProperty("pj.nt", options.c);
    threadsAvail = Integer.parseInt(options.c);
}

int numParallel = 1;
@Field int threadsPer = threadsAvail;

// Compute the gradient for each step.
if (options.g) {
    gradient = Boolean.parseBoolean(options.g);
}

// Print the energy for each step.
if (options.v) {
    print = Boolean.parseBoolean(options.v);
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
    //logger.severe("" + threadsPer);
}

if (options.la1) {
    ranges1 = options.la1.tokenize(".");
}
if (options.la2) {
    ranges2 = options.la2.tokenize(".");
}

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
 * Handles opening a file (filenmae), with 0-indexed number topNum.
 */
private void openFile(String toOpen, int topNum) {
    open(toOpen, threadsPer);
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



String filename = arguments.get(0);
StringBuilder sb = new StringBuilder("\n Timing energy and gradient for " + filename);
for (int i = 1; i < arguments.size(); i++) {
    sb.append(",").append(arguments.get(i));
}
logger.info(sb.toString());

for (int i = 0; i < arguments.size(); i++) {
    openFile(arguments.get(i), i);
}
/*for (topFile in arguments) {
    open(topFile);
    topologies.add(active);
    energies.add(active.getPotentialEnergy());
}*/
//ForceFieldEnergy energy = active.getPotentialEnergy();
Potential energy;
switch (arguments.size()) {
    case 1:
        energy = energies[0];
        break;
    case 2:
        energy = new DualTopologyEnergy(topologies[0], topologies[1]);
        energy.setLamda(lambda);
        if (numParallel == 2) {
            energy.setParallel(true);
        }
        break;
    case 4:
        DualTopologyEnergy dta = new DualTopologyEnergy(topologies[0], topologies[1]);
        DualTopologyEnergy dtb = new DualTopologyEnergy(topologies[3], topologies[2]);
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
        energy = new QuadTopologyEnergy(dta, dtb, uniqueA, uniqueB);
        if (numParallel >= 2) {
            energy.setParallel(true);
            if (numParallel == 4) {
                dta.setParallel(true);
                dtb.setParallel(true);
            }
        }
        energy.setLambda(lambda);
        break;
    default:
        logger.severe(String.format(" Must have 1, 2, or 4 file arguments; found %d", arguments.size()));
        break;
}

long minTime = Long.MAX_VALUE;
double sumTime2 = 0.0;
int halfnEvals = (nEvals % 2 == 1) ? (nEvals/2) : (nEvals/2) - 1; // Halfway point
int nVars = energy.getNumberOfVariables();
double[] x = new double[nVars];
energy.getCoordinates(x);
double[] g = gradient ? new double[nVars] : null;
def eCall = gradient ? { energy.energyAndGradient(x, g, print) } : { energy.energy(x, print) };

for (int i=0; i<nEvals; i++) {
    long time = -System.nanoTime();
    //energy.energy(gradient, print);
    eCall();
    time += System.nanoTime();
    minTime = time < minTime ? time : minTime;
    if (i >= (int) (nEvals/2)) {
        double time2 = time * 1.0E-9;
        sumTime2 += (time2*time2);
    }
}
++halfnEvals;
double rmsTime = Math.sqrt(sumTime2/halfnEvals);
logger.info(String.format(" Minimum time: %14.5f (sec)", minTime * 1.0E-9));
logger.info(String.format(" RMS time (latter half): %14.5f (sec)", rmsTime));
for (int i = 0; i < energies.size(); i++) {
    int numt = ((ForceFieldEnergy) energies[i]).parallelTeam.getThreadCount();
    logger.info(String.format(" Number of threads for topology %d: %d", i, numt));
}
