/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

// MINIMIZE

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports

import ffx.algorithms.Minimize;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;

// Default convergence criteria.
double eps = 1.0;

// First softcore atom for topology 1.
double s = -1;

// First softcore atom for topology 2.
double s2 = -1;

// Last softcore atom for topology 1.
double f = -2;

// Last softcore atom for topology 2.
double f2 = -2;

// First atom for no electrostatics.
int noElecStart = 1;

// Last atom for no electrostatics.
int noElecStop = -2;

// First atom of the 2nd topology for no electrostatics.
int noElecStart2 = 1;

// Last atom of the 2nd topology for no electrostatics.
int noElecStop2 = -2;

// Fixed lambda value.
double lambda = -1;

// Things below this line normally do not need to be changed.
// ===============================================================================================
for (String arg : args) {
    logger.info(arg);
}

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc minimize [options] <filename1> [filename2]');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName:'0', 'Starting ligand atom.');
cli.s2(longOpt:'start2', args:1, argName:'0', 'Starting ligand atom for the 2nd topology.');
cli.f(longOpt:'final', args:1, argName:'0', 'Final ligand atom.');
cli.f2(longOpt:'final2', args:1, argName:'0', 'Final ligand atom for the 2nd topology.');
cli.es1(longOpt:'noElecStart1', args:1, argName:'1', 'No Electrostatics Starting Atom.');
cli.es2(longOpt:'noElecStart2', args:1, argName:'1', 'No Electrostatics Starting Atom for the 2nd Topology.');
cli.ef(longOpt:'noElecFinal', args:1, argName:'-1', 'No Electrostatics Final Atom.');
cli.ef2(longOpt:'noElecfinal2', args:1, argName:'-1', 'No Electrostatics Final Atom for the 2nd topology.');
cli.l(longOpt:'lambda', args:1, argName:'0.0', 'Initial lambda value.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.p(longOpt:'polarization', args:1, 'polarization model: [none / direct / mutual]');
def options = cli.parse(args); 
    
if (options.h) {
    return cli.usage();
}

// Starting ligand atom.
if (options.s) {
    s = Integer.parseInt(options.s);
}

// Final ligand atom.
if (options.f) {
    f = Integer.parseInt(options.f);
}

// Starting ligand atom for the 2nd topology.
if (options.s2) {
    s2 = Integer.parseInt(options.s2);
}

// Final ligand atom for the 2nd topology.
if (options.f2) {
    f2 = Integer.parseInt(options.f2);
}

// No electrostatics on Topology 1.
if (options.es1) {
    noElecStart = Integer.parseInt(options.es1);
}

// First atom from Topology 1 with no electrostatics.
if (options.ef) {
    noElecStop = Integer.parseInt(options.ef);
}

// No electrostatics on Topology 2.
if (options.es2) {
    noElecStart2 = Integer.parseInt(options.es2);
}

// First atom from Topology 2 with no electrostatics.
if (options.ef2) {
    noElecStop2 = Integer.parseInt(options.ef2);
}

// Starting lambda value.
if (options.l) {
    lambda = Double.parseDouble(options.l);
}

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);    
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

List<String> arguments = options.arguments();
String filename = null;
String filename2 = null;
MolecularAssembly[] systems = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    filename = arguments.get(0);
    systems = open(filename);
    if (arguments.size() > 1) {
        filename2 = arguments.get(1);
    }
} else if (active == null) {
    return cli.usage();
} else {
    filename = active.getFile();
}

boolean lambdaTerm = false;
if ((lambda >= 0.0 && lambda <= 1.0) || filename2 != null) {
    lambdaTerm = true;
    System.setProperty("lambdaterm", 'true');
    // Check that lambda is within the limit 0..1
    if (lambda < 0.0 || lambda > 1.0) {
        lambda = 0.0;
    }
}

if (lambdaTerm) {
    // Get a reference to the first system's ForceFieldEnergy.
    ForceFieldEnergy energy = active.getPotentialEnergy();
    // Set the lambda value.
    energy.setLambda(lambda);
    // Apply the ligand atom selection
    Atom[] atoms = active.getAtomArray();
    for (int i = s; i <= f; i++) {
        Atom ai = atoms[i - 1];
        ai.setApplyLambda(true);
        ai.print();
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
}

if (filename2 == null) {
    logger.info("\n Running minimize on " + active.getName());
    logger.info(" RMS gradient convergence criteria: " + eps);

    // Do the minimization
    energy();
    e = minimize(eps);
    energy();
    
    String ext = FilenameUtils.getExtension(filename);
    filename = FilenameUtils.removeExtension(filename);
    if (ext.toUpperCase().contains("XYZ")) {
        saveAsXYZ(new File(filename + ".xyz"));
    } else {
        saveAsPDB(systems, new File(filename + ".pdb"));
    }
} else {
    logger.info("\n Running minimize on a dual topology from " + filename + " and " + filename2);
    logger.info(" RMS gradient convergence criteria: " + eps);

    // Save a reference to the first topology.
    MolecularAssembly topology1 = active;

    systems2 = open(filename2);

    ForceFieldEnergy energy = active.getPotentialEnergy();
    atoms = active.getAtomArray();
    // Apply the ligand atom selection for the 2nd topology.
    for (int i = s2; i <= f2; i++) {
        Atom ai = atoms[i - 1];
        ai.setApplyLambda(true);
        ai.print();
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

    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff(0.0);

    // Create the DualTopology potential energy.
    DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topology1, active);
    dualTopologyEnergy.setLambda(lambda);

    Minimize minimize = new Minimize(topology1, dualTopologyEnergy, sh);

    minimize.minimize(eps);

    // Save topology 1
    String ext = FilenameUtils.getExtension(filename);
    filename = FilenameUtils.removeExtension(filename);

    if (ext.toUpperCase().contains("XYZ")) {
        dat.setActive(0);
        saveAsXYZ(new File(filename + ".xyz"));
    } else {
        saveAsPDB(systems, new File(filename + ".pdb"));
    }

    // Save topology 2
    ext = FilenameUtils.getExtension(filename2);
    filename2 = FilenameUtils.removeExtension(filename2);

    if (ext.toUpperCase().contains("XYZ")) {
        dat.setActive(1);
        saveAsXYZ(new File(filename2 + ".xyz"));
    } else {
        saveAsPDB(systems2, new File(filename2 + ".pdb"));
    }
}

