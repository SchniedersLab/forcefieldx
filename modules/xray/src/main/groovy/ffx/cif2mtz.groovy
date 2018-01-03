/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import org.apache.commons.io.FilenameUtils

import ffx.crystal.Crystal
import ffx.crystal.ReflectionList
import ffx.crystal.Resolution
import ffx.xray.DiffractionRefinementData
import ffx.xray.parsers.CIFFilter
import ffx.xray.parsers.MTZWriter
import ffx.xray.parsers.MTZWriter.MTZType

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc cif2mtz [options] <pdbfilename> <ciffilename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 2) {
    return cli.usage();
}

// Name of the PDB with crystal header information
String modelfilename = arguments.get(0);

// input CIF file
String datafilename = arguments.get(1);

systems = open(modelfilename);

CIFFilter ciffilter = new CIFFilter();
ReflectionList reflectionlist = ciffilter.getReflectionList(new File(datafilename), systems[0].getProperties());

if (reflectionlist == null) {
    println("Using crystal information from PDB to generate MTZ file");

    Crystal crystal = systems[0].getCrystal().getUnitCell();
    double res = ciffilter.getResolution(new File(datafilename), crystal);
    if (res < 0.0) {
        println("resolution could not be determined from PDB and CIF file");
        return;
    }

    Resolution resolution = new Resolution(res);
    reflectionlist = new ReflectionList(crystal, resolution, systems[0].getProperties());
}

DiffractionRefinementData refinementdata = new DiffractionRefinementData(systems[0].getProperties(), reflectionlist);
ciffilter.readFile(new File(datafilename), reflectionlist, refinementdata, systems[0].getProperties());

MTZWriter mtzwriter = new MTZWriter(reflectionlist, refinementdata, FilenameUtils.removeExtension(datafilename) + "_cif.mtz", MTZType.DATAONLY);
mtzwriter.write();
