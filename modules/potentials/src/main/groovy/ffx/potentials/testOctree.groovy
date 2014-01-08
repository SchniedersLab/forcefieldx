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
import ffx.potential.nonbonded.Octree;
import javax.vecmath.Point3d

int maxTreeDepth = 10;
int maxAtomsPerVolume = 20;
boolean writePartitionFile = false;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc biotype [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.d(longOpt:'maxDepth', args:1, argName:'10', 'Maximum octree depth.');
cli.m(longOpt:'maxOccupancy', args:1, argName:'20', 'Maximum number of atoms per sub-volume.');
cli.w(longOpt:'partFile', args:0, 'Write atom partitioning file.');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

if (options.d) {
    maxTreeDepth = Integer.parseInt(options.d);
}
if (options.m) {
    maxAtomsPerVolume = Integer.parseInt(options.m);
}
if (options.w) {
    writePartitionFile = true;
}

// Read in command line.
String pdbName = arguments.get(0);
systems = open(pdbName);

File structureFile = new File(FilenameUtils.normalize(pdbName));
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File partFile = new File(baseFilename + ".prt");

double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE, minZ = Double.MAX_VALUE;
double maxX = Double.MIN_VALUE, maxY = Double.MIN_VALUE, maxZ = Double.MIN_VALUE;
List<Atom> atoms = active.getAtomList();
for (Atom atom : atoms) {
    double atomX = atom.getX();
    double atomY = atom.getY();
    double atomZ = atom.getZ();
    double vdwr = atom.getVDWR();
    minX = (atomX - vdwr < minX) ? atomX - vdwr : minX;
    minY = (atomY - vdwr < minY) ? atomY - vdwr : minY;
    minZ = (atomZ - vdwr < minZ) ? atomZ - vdwr : minZ;
    maxX = (atomX + vdwr > maxX) ? atomX + vdwr : maxX;
    maxY = (atomY + vdwr > maxY) ? atomY + vdwr : maxY;
    maxZ = (atomZ + vdwr > maxZ) ? atomZ + vdwr : maxZ;
}
logger.info(String.format("\n\n System corners:"));
logger.info(String.format("    (%4.2f, %4.2f, %4.2f)", minX, minY, minZ));
logger.info(String.format("    (%4.2f, %4.2f, %4.2f)\n", maxX, maxY, maxZ));

Point3d corner = new Point3d(minX, minY, minZ);
double edgeLength = Math.max(Math.max(maxX - minX, maxY - minY), maxZ - minZ);
logger.info(String.format(" Total cube volume: %10.2g\n", Math.pow(edgeLength,3)));
    
Octree octree = new Octree(0, corner, edgeLength);
octree.setMaxTreeDepth(maxTreeDepth);
octree.setMaxAtomsPerVolume(maxAtomsPerVolume);

octree.addAtoms(atoms);
octree.debugPrintStats(writePartitionFile, partFile);

return;

//String mol = FilenameUtils.getBaseName(xyzname);
//List atoms = active.getAtomList();
//int index = 1;
//for (Atom atom : atoms) {
//    StringBuilder sb = new StringBuilder();
//    sb.append(String.format(" biotype %3d %4s \"%s\" %3d", index++, atom.getName(), mol, atom.getAtomType().type));
//    List bonds = atom.getBonds();
//    if (bonds != null) {
//        for (Bond bond : bonds) {
//            sb.append(String.format(" %4s", bond.get1_2(atom).getName()));
//        }
//    }
//    logger.info(sb.toString());
//}
