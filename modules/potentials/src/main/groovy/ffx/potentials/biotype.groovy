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

// BIOTYPE

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc biotype [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String xyzname = arguments.get(0);

systems = open(xyzname);
energy();

String mol = FilenameUtils.getBaseName(xyzname);
List atoms = active.getAtomList();
int index = 1;
for (Atom atom : atoms) {
    StringBuilder sb = new StringBuilder();
    sb.append(String.format(" biotype %3d %4s \"%s\" %3d", index++, atom.getName(), mol, atom.getAtomType().type));
    List bonds = atom.getBonds();
    if (bonds != null) {
        for (Bond bond : bonds) {
            sb.append(String.format(" %4s", bond.get1_2(atom).getName()));
        }
    }
    logger.info(sb.toString());
}

return;

