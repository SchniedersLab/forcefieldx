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

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.ForceFieldFilter
import ffx.potential.parsers.PDBFilter
import org.apache.commons.configuration.CompositeConfiguration

boolean writemaps = false;

boolean writemtz = false;

boolean timings = false;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.scaleBulk [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual ]');
cli.m(longOpt:'maps', 'set to output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps');
cli.t(longOpt:'timings', 'set to perform FFT test timings');
cli.w(longOpt:'mtz', 'write out MTZ containing structure factor coefficients');
//cli.o(longOpt:'omit-mode', 'Treat all rotamer-optimizable residues as Alanine during map generation.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h) {
    return cli.usage();
}

String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    modelfilename = arguments.get(0);
} else {
    return cli.usage();
}


// set up diffraction data (can be multiple files)
List diffractionfiles = new ArrayList();
if (arguments.size() > 1) {
    DiffractionFile diffractionfile = new DiffractionFile(arguments.get(1), 1.0, false);
    diffractionfiles.add(diffractionfile);
}
if (options.d) {
    for (int i=0; i<options.ds.size(); i+=3) {
	double wA = Double.parseDouble(options.ds[i+1]);
	boolean neutron = Boolean.parseBoolean(options.ds[i+2]);
	DiffractionFile diffractionfile = new DiffractionFile(options.ds[i], wA, neutron);
	diffractionfiles.add(diffractionfile);
    }
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

if (options.m) {
    writemaps = true;
}

if (options.t) {
    timings = true;
}

if (options.w) {
    writemtz = true;
}

systems = open(modelfilename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

if (writemtz) {
    diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_ffx.mtz");
}

if (writemaps) {
//    if (options.o) {
//        def polymers = active.getChains();
//        int nPolymers = polymers.length;
//        for (int p=0; p < nPolymers; p++) {
//            Polymer polymer = polymers[p];
//            ArrayList<Residue> residues = polymer.getResidues();
//            for (int i=0; i < residues.size(); i++) {
//                Residue residue = residues.get(i);
//                Rotamer[] rotamers = RotamerLibrary.getRotamers(residue);
//                if (rotamers != null) {
//                    logger.info("Turning off atoms of " + residue.toString());
//                    RotamerOptimization.turnOffAtoms(residue);
//                    logger.info("Turning on CBeta of " + residue.toString());
//                    RotamerOptimization.turnOnCBeta(residue);
//                }
//            }
//        }
//    }
    diffractiondata.writeMaps(FilenameUtils.removeExtension(modelfilename) + "_ffx");
}

if (timings) {
    diffractiondata.timings();
}
