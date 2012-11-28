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

//POTENTIAL

import groovy.util.CliBuilder

//Supply a choice (1 - 4)
//1. Get QM Potential from a Gaussian CUBE File
//2. Calculate the Model Potential for a System
//3. Compare a Model Potential to a Target Grid
//4. Fit Electrostatic Parameters to Target Grid
//if choice == 1, then supply the cube filename
//if choice == 2 || choice == 3, then supply the xyz filename
//if choice == 4, then supply the xyz filename and the eps

def cli = new CliBuilder(usage:' ffxc potential <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}

//Get XYZ File
Double eps = null;
Integer choice = Integer.parseInt(arguments.get(0));
String filename = arguments.get(1);
int out_type = 5;
if(choice == 4){
    eps = Double.parseDouble(arguments.get(2));
}
potential(choice, filename, eps);
