
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

// LOOPCLOSURE

// Groovy Imports
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Loop
import ffx.potential.parsers.PDBFilter
import groovy.util.CliBuilder;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc loopClosure [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName: '-1', 'First residue to draw. Indexed from residue 0.');
cli.f(longOpt:'finish', args:1, argName: '-1', 'Last residue to draw. Indexed from residue 0.');

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

int stt_res = -1;
int end_res = -1;

 //Start point
if (options.s) {
    stt_res = Integer.parseInt(options.s);
}
else
{
    stt_res = 1;
}

// Finish point
if (options.f) {
    end_res = Integer.parseInt(options.f);
} 
else
{
    end_res = 5;
}


if (end_res - stt_res != 4)
{
    logger.info("\n\nInvalid residue range. Residue range must consist of only 5 residues.\n\n");
}

List<String> arguments = options.arguments();
String modelFilename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    modelFilename = arguments.get(0);
    open(modelFilename);
} else if (active == null) {
    return cli.usage();
} else {
    modelFilename = active.getFile();
}

logger.info("\n Running loopClosure on " + modelFilename);
boolean writeFile = true;
ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
//Loop loop = new Loop(active, stt_res, end_res, writeFile);

Loop loop = new Loop(active);

List<double[]> loopSolutions = loop.generateLoops(stt_res, end_res);
for(int i = 0; i < loopSolutions.size(); i++){
    forceFieldEnergy.setCoordinates(loopSolutions.get(i));
    File modifiedFile = new File("loop_"+modelFilename + "_"+i);
    PDBFilter modFilter = new PDBFilter(modifiedFile, active, null, null);
    if (writeFile) {
        modFilter.writeFile(modifiedFile, true);
    }        
}

/*
 * Example of using the new PotentialsFunctions interface instead of Groovy method
 * closures:

logger.info("\n Running energy on " + modelfilename);
PotentialsFunctions functions; // This is an interface specifying the closure-like methods.
try {
    // Use a method closure to try to get an instance of UIUtils (the User Interfaces
    // implementation, which interfaces with the GUI, etc.).
    functions = getPotentialsFunctions();
} catch (MissingMethodException ex) {
    // If Groovy can't find the appropriate closure, catch the exception and build
    // an instance of the local implementation.
    functions = new PotentialsUtils();
}

// Use PotentialsFunctions methods instead of Groovy method closures to do work.
MolecularAssembly[] assemblies = functions.open(modelfilename);
MolecularAssembly activeAssembly = assemblies[0];
functions.energy(activeAssembly);

*/