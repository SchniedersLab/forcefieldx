/**
 * Title: Force Field X 
 * Description: Force Field X - Software for Molecular Biophysics 
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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

// MINIMIZE 2

import groovy.util.CliBuilder

def cli = new CliBuilder(usage:' ffxc minimize_2 [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}

String xyzfilename = arguments.get(0);
double eps = Double.parseDouble(arguments.get(1));
if(arguments.size() > 3 && arguments.get(2).equals("-k")){
    String keyfname = arguments.get(3);
    minimize_2(xyzfilename,eps,keyfname);
}
else{
    minimize_2(xyzfilename,eps,null);
}
