
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

// MOLECULAR & STOCHASTIC DYNAMICS

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// PJ Imports
import edu.rit.pj.Comm;

// FFX Imports
import ffx.algorithms.mc.MolecularMC;
import ffx.algorithms.mc.CoordShakeMove;
import ffx.algorithms.mc.MCMove;
import ffx.potential.AssemblyState;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Residue;

// Number of Monte Carlo steps
int nSteps = 1000000;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Shake coordinates every n steps.
int xShakeFreq = 1;

// Normal distribution width for coordinate shakes.
double xShakeWidth = 0.1;

// File type of snapshots.
String fileType = "PDB";

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc md [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of Monte Carlo steps.');
// cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
cli.xf(longOpt:'xShakeFreq', args:1, argName:'1', 'Shake coordinates every n steps.');
cli.xw(longOpt:'xShakeWidth', args:1, argName:'0.1', 'Normal distribution width for coordinate shakes.');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Load the number of Monte Carlo steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Set file type to save to.
if (options.f) {
    fileType = options.f.toUpperCase();
}

// Write snapshot interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

if (options.xw) {
    xShakeWidth = Double.parseDouble(options.xw);
}

if (options.xf) {
    xShakeFreq = Integer.parseInt(options.xf);
}

List<String> arguments = options.arguments();
String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    modelfilename = arguments.get(0);
    open(modelfilename);
} else if (active == null) {
    return cli.usage();
} else {
    modelfilename = active.getFile();
}

boolean master = true;
if (Comm.world().size() > 1) {
    int rank = Comm.world().rank();
    if (rank != 0) {
        master = false;
    }
}

logger.info("\n Running Metropolis Monte Carlo on " + modelfilename);

MolecularMC mc = new MolecularMC(active);
mc.setTemperature(temperature);

// Preliminary testing implementation.
List<MCMove> moveList = new ArrayList<>();
CoordShakeMove move = new CoordShakeMove(active);
move.setSigma(xShakeWidth);
moveList.add(move);

// Run the first step to get it running.
double eOpt = mc.tryMove(moveList) ? mc.getE2() : mc.getE1();
AssemblyState optimum = new AssemblyState(active);

for (int i = 1; i < nSteps; i++) {
    if (mc.tryMove(moveList) && mc.getE2() < eOpt) {
        eOpt = mc.getE2();
        optimum = new AssemblyState(active);
    }
}

String ext = FilenameUtils.getExtension(modelfilename);
String basename = FilenameUtils.removeExtension(modelfilename);
if (master) {
    String finalFileName = basename + "_fin.pdb";
    
    if (fileType.equalsIgnoreCase("PDB")) {
        logger.info(String.format(" Printing out final structure to file %s", finalFileName));
        saveAsPDB(new File(finalFileName));
    } else if (fileType.equalsIgnoreCase("XYZ")) {
        finalFileName = basename + "_fin.xyz";
        logger.info(String.format(" Printing out final structure to file %s", finalFileName));
        saveAsXYZ(new File(finalFileName));
    } else {
        logger.info(String.format(" File extension type %s not recognized, saving as PDB", fileType));
        logger.info(String.format(" Printing out final structure to file %s", finalFileName));
        saveAsPDB(new File(finalFileName));
    }
}

logger.info(String.format(" Optimum energy %10.6f", eOpt));
optimum.revertState();

if (master) {
    energy();
    String finalFileName = basename + "_opt.pdb";
    
    if (fileType.equalsIgnoreCase("PDB")) {
        logger.info(String.format(" Printing out optimum structure found to file %s", finalFileName));
        saveAsPDB(new File(finalFileName));
    } else if (fileType.equalsIgnoreCase("XYZ")) {
        finalFileName = basename + "_opt.xyz";
        logger.info(String.format(" Printing out optimum structure found to file %s", finalFileName));
        saveAsXYZ(new File(finalFileName));
    } else {
        logger.info(String.format(" File extension type %s not recognized, saving as PDB", fileType));
        logger.info(String.format(" Printing out optimum structure found to file %s", finalFileName));
        saveAsPDB(new File(finalFileName));
    }
}