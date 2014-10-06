/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */

package ffx.xray

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Java Imports
import java.nio.file.Paths;
import java.nio.file.Path;

// Force Field X X-ray Imports
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

// Force Field X Minimize imports
import ffx.algorithms.Minimize;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.CoordinateFileFilter;

// Force Field X real-space Imports
import ffx.xray.RealSpaceData;
import ffx.xray.RealSpaceFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.RescoreAndCluster

// Other Force Field X imports.
import ffx.utilities.DoubleIndexPair;

double boltzRestraint = 0.0;
String targetFunct = "xray_energy";
boolean useFree = false; // I very strongly suggest this remains false.
String suffix = "_rw";
double eps = 0.1;
int maxiter = 1000;
List<DiffractionFile> diffractionFiles = new ArrayList<>();
List<File> modelFiles = new ArrayList<>();
double[] weights;
int numModels = 0;


// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc rescore [options] <ensemblefilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.br(longOpt:'boltzmannRestraint', args:1, argName:'0.0', 'Coefficient of harmonic restraint to initial Boltzmann weights for reweighting.');
cli.t(longOpt:'targetFunction', args:1, argName:'xray_energy', 'Sets target function for initial Boltzmann weight generation (none, energy, xray_energy, or xray). If no target function, all start with equal weights.');
cli.f(longOpt:'freeData', args:1, argName:'false', 'Include free reflections while reweighting.');
cli.s(longOpt:'suffix', args:1, argName:'_rw', 'Suffix for reweighted .mtz and weights text file');
cli.e(longOpt:'eps', args:1, argName:'0.1', 'RMS gradient convergence criteria for reweighting');
cli.xd(longOpt:'xrayData', args:3, valueSeparator:',', argName:'data,1.0,false', 'specify input xray data filename or filenames (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment.');
cli.i(longOpt:'maxiter', args:1, argName:'1000', 'maximum number of allowed reweighting iterations');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual]');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1 || arguments.size > 2) {
    logger.info(" Must have an ensemble file name (or directory) and a data file name (or data specified by xd flag).");
    return cli.usage();
}

if (options.br) {
    boltzRestraint = Double.parseDouble(options.br);
}

if (options.f) {
    useFree = Boolean.parseBoolean(options.f);
}

if (options.s) {
    suffix = options.s;
}

if (options.e) {
    eps = Double.parseDouble(options.e);
}

if (options.i) {
    maxiter = Integer.parseInt(options.i);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

if (options.t) {
    String targF = options.t.toLowerCase();
    switch (targF) {
    case "xray_energy":
        targetFunct = "xray_energy";
        break;
    case "energy":
        targetFunct = "energy";
        break;
    case "xray":
        targetFunct = "xray";
        break;
    default: // Bit of a weird positioning: default falls through into the "none" case.
        logger.warning(" No valid target function entered: must be none, energy, xray_energy, or xray. Setting to \"none\"");
    case "none":
        targetFunct = "none";
        break;
    }
}

