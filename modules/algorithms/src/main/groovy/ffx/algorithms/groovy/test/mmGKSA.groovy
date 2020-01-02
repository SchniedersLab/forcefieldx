//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.groovy.test;

import groovy.cli.picocli.CliBuilder

import ffx.algorithms.misc.MMgksa
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

// Default weight parameters
double elecWt = 1.0;
double solvWt = 1.0;
double vdwWt = 1.0;

int maxFrames = -1;
int freq = 1;
boolean decompose = false;

// Things below this line normally do not need to be changed.
// ===============================================================================================
for (String arg : args) {
    logger.info(arg);
}

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc mmGKSA [options] <filename>');
cli.h(longOpt: 'help', 'Print this help message.');
cli.l(longOpt: 'ligand', args: 1, argName: '0', 'Ligand atoms');
cli.i(longOpt: 'ignore', args: 1, argName: '0', 'Atoms to ignore (not in ligand or protein)');
cli.f(longOpt: 'frequency', args: 1, argName: '1', 'Evaluate every nth frame');
cli.m(longOpt: 'maxFrames', args: 1, argName: '-1', 'Evaluate at most this many frames (-1 to evaluate to end of file)');
cli.e(longOpt: 'electrostaticWeight', args: 1, argName: '1.0', 'Weight to electrostatic interactions.');
cli.s(longOpt: 'solvationWeight', args: 1, argName: '1.0', 'Weight to solvation interactions.');
cli.v(longOpt: 'vanDerWaalsWeight', args: 1, argName: '1.0', 'Weight to van der Waals interactions');
cli.d(longOpt: 'decompose', args: 1, argName: 'false', 'Decompose electrostatics and solvation to components.');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

String filename = null;

filename = arguments.get(0);
//MolecularAssembly[] systems = open(filename);
PotentialsFunctions functs;
try {
    functs = getPotentialsUtils();
} catch (MissingMethodException ex) {
    functs = new PotentialsUtils();
}

MolecularAssembly mola = functs.open(filename)[0];
SystemFilter filter = functs.getFilter();
Atom[] atoms = mola.getAtomArray();
int nAtoms = atoms.length;

Set<Atom> iAtoms = new HashSet<>();
if (options.i) {
    String[] optiToks = options.i.split(",");
    for (String itok : optiToks) {
        try {
            int[] inactAtoms = SystemFilter.parseAtNumArg("", itok, nAtoms);
            for (int i = inactAtoms[0]; i <= inactAtoms[1]; i++) {
                iAtoms.add(atoms[i]);
            }
        } catch (IllegalArgumentException ex) {
            logger.warning(ex.toString());
        }
    }
}

if (options.e) {
    elecWt = Double.parseDouble(options.e);
}
if (options.s) {
    solvWt = Double.parseDouble(options.s);
}
if (options.v) {
    vdwWt = Double.parseDouble(options.v);
}
if (options.f) {
    freq = Integer.parseInt(options.f);
}
if (options.m) {
    maxFrames = Integer.parseInt(options.m);
}
if (options.d) {
    decompose = Boolean.parseBoolean(options.d);
}

boolean useLigand = false;
Atom[] ligAtoms;
Atom[] protAtoms = atoms; // And yes, I know this is an array reference.
if (options.l) {
    Set ligAtomSet = new HashSet<>();
    String[] optltoks = options.l.split(",");
    for (String ltok : optltoks) {
        try {
            int[] lAtoms = SystemFilter.parseAtNumArg("", ltok, nAtoms);
            for (int i = lAtoms[0]; i <= lAtoms[1]; i++) {
                ligAtomSet.add(atoms[i]);
            }
        } catch (IllegalArgumentException ex) {
            logger.warning(ex.toString());
        }
    }
    if (ligAtomSet.isEmpty()) {
        logger.warning(" No valid input for ligand atoms! Defaulting to scoring trajectory without binding");
    } else {
        useLigand = true;
        int nLig = ligAtomSet.size();
        ligAtoms = new Atom[nLig];
        ligAtomSet.toArray(ligAtoms);
        protAtoms = new Atom[nAtoms - nLig];
        List<Atom> patomList = new ArrayList<>(nAtoms - nLig);
        for (Atom atom : atoms) {
            if (!ligAtoms.contains(atom)) {
                patomList.add(atom);
            }
        }
        patomList.toArray(protAtoms);
    }
} else {
    logger.info(" Scoring trajectory without binding calculations");
}

MMgksa mmGKSA;
if (useLigand) {
    mmGKSA = new MMgksa(mola, functs, filter, protAtoms, ligAtoms);
} else {
    mmGKSA = new MMgksa(mola, functs, filter);
}

if (!iAtoms.isEmpty()) {
    Atom[] iAtArray = new Atom[iAtoms.size()];
    iAtoms.toArray(iAtArray);
    mmGKSA.setIgnoredAtoms(iAtArray);
}

mmGKSA.setElectrostaticsWeight(elecWt);
mmGKSA.setSolvationWeight(solvWt);
mmGKSA.setVdwWeight(vdwWt);
mmGKSA.setDecompose(decompose);
mmGKSA.runMMgksa(freq, maxFrames);

