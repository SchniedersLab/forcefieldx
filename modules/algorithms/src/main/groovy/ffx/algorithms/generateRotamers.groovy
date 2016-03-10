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

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// ENERGY
import ffx.potential.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;

import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.algorithms.GenerateRotamers;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports

// Things below this line normally do not need to be changed.
// ===============================================================================================
Residue residue = null;
int nChi;
File outFile;
AminoAcid3 baseResidue = AminoAcid3.UNK;
boolean useBase = false;
int initDepth = 0;
int finalDepth = -1;
double incr = 10.0;
double width = 180.0; // Internally half of overall search width (so -180 and +180 to get 360 degrees search).
int library = 2;
boolean verbose = false;
String videoFile = null;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc generateRotamers [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.r(longOpt:'residue', args:1, argName:'resid', 'Required: Specify residue as ChainNumber (e.g. B71).');
cli.n(longOpt:'nChi', args:1, argName:'1-7', 'Required: Specify number of torsions per rotamer.');
cli.b(longOpt:'baseResidue', args:1, argName:'aa3', '3-letter code for standard amino acid to start from (or UNK for patch-defined).');
cli.i(longOpt:'initialDepth', args:1, argName:'1', 'First chi to operate on.');
cli.f(longOpt:'finalDepth', args:1, argName:'-1', 'Last chi to operate on (or -1 to optimize to nChi).');
cli.a(longOpt:'angleIncrement', args:1, argName:'10', 'Angle in degrees to spin torsions by.');
cli.w(longOpt:'width', args:1, argName:'360.0', 'Total width to scan in degrees.');
cli.o(longOpt:'outFile', args:1, argName:'file.tor.csv', 'File to print output to.');
cli.l(longOpt:'library', args:1, argName:'2', 'Available rotamer libraries are Ponder and Richards (1) or Richardson (2).');
cli.v(longOpt:'verbose', args:1, argName:'false', 'Log rotamer energies to console.');
cli.vw(longOpt:'videoWriter', args: 1, argName:'file', 'Writes video to a file.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);
open(filename);

nChi = Integer.parseInt(options.n);
String resn = options.r;
char chain = resn.charAt(0);
int resIndex = Integer.parseInt(resn.substring(1));
for (Polymer polymer : active.getChains()) {
    if (polymer.getChainID() == chain) {
        for (Residue pres : polymer.getResidues()) {
            if (pres.getResidueNumber() == resIndex) {
                residue = pres;
            }
        }
    }
}

if (residue == null) {
    logger.severe(String.format(" Could not find residue %s", resn));
}

if (options.b) {
    try {
        baseResidue = AminoAcid3.valueOf(options.b.toUpperCase());
        useBase = true;
    } catch (Exception ex) {
        logger.warning(String.format(" Could not find residue matching %s", options.b));
    }
}

if (options.i) {
    initDepth = Integer.parseInt(options.i) - 1;
}

if (options.f) {
    finalDepth = Integer.parseInt(options.f) - 1;
}

if (options.w) {
    width = Double.parseDouble(options.w) / 2.0;
}

if (options.a) {
    incr = Double.parseDouble(options.a);
}

if (options.o) {
    outFile = new File(options.o);
} else {
    String newName = FilenameUtils.getBaseName(filename);
    outFile = new File(newName + ".rotout");
}

if (options.l) {
    library = Integer.parseInt(options.l);
}

if (options.v) {
    verbose = Boolean.parseBoolean(options.v);
}

if (options.vw) {
    videoFile = options.vw;
}

/**
 * This needs to come before setting the baseline residue.
 */
if (library == 1) {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
} else {
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
}

GenerateRotamers genr = new GenerateRotamers(active, active.getPotentialEnergy(), residue, outFile, nChi, sh);
genr.setDepth(initDepth, finalDepth);
if (useBase) {
    genr.setBaselineAARes(baseResidue);
}
genr.setIncrement(incr);
genr.setSearchWidth(width);
genr.setPrint(verbose);
genr.setVideo(videoFile);
energy();
genr.tryRotamers();