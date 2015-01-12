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

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;

// First atom to test.
int atomID = 0;

// atomic number to test
int atomicNum = -1;

// Finite-difference step size in Angstroms.
double step = 1.0e-4;

// Print out the energy for each step.
boolean print = false;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc xray.testGradient [options] <pdbfilename> [datafilename]');
cli.h(longOpt:'help', 'Print this help message.');
cli.a(longOpt:'atomID', args:1, argName:'1', 'Number of the first atom to test');
cli.n(longOpt:'num', args:1, argName:'-1', 'only test gradient on a given atomic number');
cli.s(longOpt:'dx', args:1, argName:'1.0e-4', 'Finite-difference step size (Angstroms)');
cli.d(longOpt:'data', args:3, valueSeparator:',', argName:'data.mtz,1.0,false', 'specify input data filename (or simply provide the datafilename argument after the PDB file), weight applied to the data (wA) and if the data is from a neutron experiment');
cli.p(longOpt:'polarization', args:1, argName:'mutual', 'polarization model: [none / direct / mutual ]');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Name of the file (PDB or XYZ).
String modelfilename = arguments.get(0);

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

// First atom to test. Subtract 1 for Java array indexing.
if (options.a) {
    atomID = Integer.parseInt(options.a) - 1;
}

if (options.n) {
    atomicNum = Integer.parseInt(options.n);
}

// Load the finite-difference step size in Angstroms.
if (options.s) {
    step = Double.parseDouble(options.s);
}

systems = open(modelfilename);

if (diffractionfiles.size() == 0) {
    DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
    diffractionfiles.add(diffractionfile);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]));

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
// energy();


diffractiondata.computeAtomicGradients(RefinementMode.COORDINATES_AND_BFACTORS);

double llk0 = diffractiondata.computeLikelihood();
double llk1, llk2, fdb;
double[] gxyz = new double[3];
double[] fd = new double[3];
double bg;
double[] ganisou = null;
double[] fdanisou = new double[6];
double[] danisou = new double[6];
double gradientTolerance = 1.0e-2;

Atom[] atoms = diffractiondata.getAtomArray();
int natoms = atoms.length;

if (atomID >= natoms) {
    atomID = 0;
}

logger.info("atom list:");
for (int i = atomID; i < natoms; i++) {
    Atom atom = atoms[i];

    logger.info("atom " + i + ": (" + atom.getXYZIndex() + ") " + atom.toString());
}

for (int i = atomID; i < natoms; i++) {
    Atom atom = atoms[i];

    if (atomicNum > 0 && atom.getAtomType().atomicNumber != atomicNum) {
        continue;
    }

    atom.getXYZGradient(gxyz);
    bg = atom.getTempFactorGradient();
    ganisou = null;
    if (atom.getAnisou() != null) {
	ganisou = atom.getAnisouGradient();
    }

    // X
    diffractiondata.crs_fc[0].deltaX(i, step);
    diffractiondata.crs_fs[0].deltaX(i, step);
    diffractiondata.computeAtomicDensity();
    llk1 = diffractiondata.computeLikelihood();
    diffractiondata.crs_fc[0].deltaX(i, -step);
    diffractiondata.crs_fs[0].deltaX(i, -step);
    diffractiondata.computeAtomicDensity();
    llk2 = diffractiondata.computeLikelihood();
    fd[0] = (llk1 - llk2) / (2.0 * step);

    double dx = gxyz[0] - fd[0];
    diffractiondata.crs_fc[0].deltaX(i, 0.0);
    diffractiondata.crs_fs[0].deltaX(i, 0.0);

    // Y
    diffractiondata.crs_fc[0].deltaY(i, step);
    diffractiondata.crs_fs[0].deltaY(i, step);
    diffractiondata.computeAtomicDensity();
    llk1 = diffractiondata.computeLikelihood();
    diffractiondata.crs_fc[0].deltaY(i, -step);
    diffractiondata.crs_fs[0].deltaY(i, -step);
    diffractiondata.computeAtomicDensity();
    llk2 = diffractiondata.computeLikelihood();
    fd[1] = (llk1 - llk2) / (2.0 * step);

    double dy = gxyz[1] - fd[1];
    diffractiondata.crs_fc[0].deltaY(i, 0.0);
    diffractiondata.crs_fs[0].deltaY(i, 0.0);

    // Z
    diffractiondata.crs_fc[0].deltaZ(i, step);
    diffractiondata.crs_fs[0].deltaZ(i, step);
    diffractiondata.computeAtomicDensity();
    llk1 = diffractiondata.computeLikelihood();
    diffractiondata.crs_fc[0].deltaZ(i, -step);
    diffractiondata.crs_fs[0].deltaZ(i, -step);
    diffractiondata.computeAtomicDensity();
    llk2 = diffractiondata.computeLikelihood();
    fd[2] = (llk1 - llk2) / (2.0 * step);

    double dz = gxyz[2] - fd[2];
    diffractiondata.crs_fc[0].deltaZ(i, 0.0);
    diffractiondata.crs_fs[0].deltaZ(i, 0.0);

    double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
    logger.info("atom " + i + ": (" + atom.getXYZIndex() + ") " + atom.toString()
        + String.format(" %s: %10.6f.", len > gradientTolerance ? "FAILED" : "passed", len)
        + String.format("\n Analytic XYZ: (%12.4f, %12.4f, %12.4f)", gxyz[0], gxyz[1], gxyz[2])
        + String.format("\n  Numeric XYZ: (%12.4f, %12.4f, %12.4f)", fd[0], fd[1], fd[2]));

    // B
    if (atom.getAnisou() == null) {
	double b = atom.getTempFactor();
	atom.setTempFactor(b + step);
	diffractiondata.computeAtomicDensity();
	llk1 = diffractiondata.computeLikelihood();
	atom.setTempFactor(b - step);
	diffractiondata.computeAtomicDensity();
	llk2 = diffractiondata.computeLikelihood();
	fdb = (llk1 - llk2) / (2.0 * step);

	double db = bg - fdb;
	atom.setTempFactor(b);

	double lenb = Math.sqrt(db * db);
	logger.info("atom " + i + ": (" + atom.getXYZIndex() + ") " + atom.toString()
            + String.format(" %s: %10.6f.", lenb > gradientTolerance ? "FAILED" : "passed", lenb)
            + String.format("\n Analytic B: (%12.4f)", bg)
            + String.format("\n  Numeric B: (%12.4f)\n", fdb));
    } else {
	double[] anisou = atom.getAnisou();
	for (int j = 0; j < 6; j++) {
	    double tmpu = anisou[j];
	    anisou[j] = tmpu + VectorMath.b2u(step);
	    diffractiondata.computeAtomicDensity();
	    llk1 = diffractiondata.computeLikelihood();
	    anisou[j] = tmpu - VectorMath.b2u(step);
	    diffractiondata.computeAtomicDensity();
	    llk2 = diffractiondata.computeLikelihood();
	    fdanisou[j] = (llk1 - llk2) / (2.0 * step);

	    danisou[j] = ganisou[j] - fdanisou[j];
	    anisou[j] = tmpu;
	}

	double lenb = Math.sqrt(danisou[0] * danisou[0] + danisou[1] * danisou[1]
            + danisou[2] * danisou[2] + danisou[3] * danisou[3]
            + danisou[4] * danisou[4] + danisou[5] * danisou[5]);
	logger.info("atom " + i + ": (" + atom.getXYZIndex() + ") " + atom.toString()
            + String.format(" %s: %10.6f.", lenb > gradientTolerance ? "FAILED" : "passed", lenb)
            + String.format("\n Analytic B(anis): (%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f)",
                ganisou[0], ganisou[1], ganisou[2], ganisou[3], ganisou[4], ganisou[5])
            + String.format("\n  Numeric B(anis): (%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f)\n",
                fdanisou[0], fdanisou[1], fdanisou[2], fdanisou[3], fdanisou[4], fdanisou[5]));
    }
}
