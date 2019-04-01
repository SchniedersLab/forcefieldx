//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.groovy.test

import org.apache.commons.io.FilenameUtils

import groovy.cli.picocli.CliBuilder

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.AminoAcidUtils
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.parsers.PDBFilter

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc exchangeCoordinates -i [start] -f [end] <toFilename> <fromFilename>');
cli.h(longOpt: 'help', 'Print this help message.');
cli.i(longOpt: 'start', args: 1, argName: '-1', 'Start residue.');
cli.f(longOpt: 'finish', args: 1, argName: '-1', 'End residue.');
cli.x(longOpt: 'start-2', args: 1, argName: '-1', 'Alternative start residue.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 2) {
    return cli.usage();
}

// Read in command line.
String toFile = arguments.get(0);
String fromFile = arguments.get(1);


int startPosition = 0;
int endPosition = 0;
// Sets where we will start and stop
if (options.i) {
    startPosition = Integer.parseInt(options.i);
    endPosition = Integer.parseInt(options.f);
} else {
    logger.info(String.format("Error reading start and end positions"));
}

int startPosition2 = startPosition;
int endPosition2 = endPosition;
// Sets where second start and stop occur for pdb files with different residue numbers
if (options.x) {
    startPosition2 = Integer.parseInt(options.x);
}

open(toFile);
MolecularAssembly original = active;
Polymer[] toPolymers = original.getChains();

// Begin Loop Here...
// For 
open(fromFile);
Polymer[] fromPolymers = active.getChains();

for (int i = startPosition; i <= endPosition; i++) {
    Residue fromResidue = fromPolymers[0].getResidue(i - startPosition + startPosition2);
    Residue toResidue = toPolymers[0].getResidue(i);

    AminoAcidUtils.copyResidue(fromResidue, toResidue);
}


File file = original.getFile();
String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
modifiedFile = new File(filename + "fullLoop.pdb");

modFilter = new PDBFilter(new File(modifiedFile.getName()), original, null, null);
modFilter.writeFile(modifiedFile, false);
