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

// ENERGY

// Groovy Imports
import groovy.util.CliBuilder;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc energy [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
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

logger.info("\n Running energy on " + modelfilename);

energy();

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