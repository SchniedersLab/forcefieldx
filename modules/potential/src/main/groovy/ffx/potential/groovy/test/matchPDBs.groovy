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
package ffx.potential.groovy.test

import groovy.cli.picocli.CliBuilder

import ffx.potential.parsers.PDBFileFilter
import ffx.potential.parsers.PDBFileMatcher
import ffx.potential.parsers.SimplePDBMatcher

File sourceFileSource; // The "good" structures: for our purposes, FFX structures.
File matchFileSource; // The structures to back-correlate and fix: MSMBuilder structures.
List<File> matchStructures;
List<File> sourceStructures;
String suffix = "_match";
boolean parseDeep = true;
boolean verbose = false;
boolean headerLink = true;
boolean parallel = true;
boolean bFactors = false;
boolean ssBonds = false;
boolean superpose = true;
boolean robustMatch = false;
boolean fixCryst = false;
boolean simple = false;
int atomsUsed = 1;
PDBFileMatcher fileMatcher;

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc test.matchPDBs [options] <sourcefiles> <filestomatch>');
cli.h(longOpt: 'help', 'Print this message.');
cli.d(longOpt: 'parseDeep', args: 1, argName: 'true', 'Accept files with valid coordinate file format (true) or look only at extension (false).');
cli.v(longOpt: 'verbose', 'Flag to print progress messages.');
cli.hL(longOpt: 'headerLink', args: 1, argName: 'true', 'Will print a HEADER line in matched PDBs with filename of source file.');
cli.b(longOpt: 'bFactors', args: 1, argName: 'false', 'Will replace B-factors in matched PDBs with B-factors from source.');
cli.sB(longOpt: 'ssBondRecords', args: 1, argName: 'false', 'Will re-introduce SSBOND records from source to matched PDBs.');
cli.s(longOpt: 'superpose', args: 1, argName: 'true', 'Will use BioJava alignment tool with default settings to calculate RMSDs.');
cli.a(longOpt: 'atomsUsed', args: 1, argName: '1', 'Over-ridden by superpose. Will use only CA (protein) and N1/9 (nucleic acid) for RMSD (1), all protein/nucleic acid atoms (2), all non-water atoms (3), or all atoms including water (4).');
cli.p(longOpt: 'parallel', args: 1, argName: 'true', 'Will parallelize comparisons');
cli.x(longOpt: 'suffix', args: 1, argName: '_match', 'Suffix for modified copies of matched files.');
cli.r(longOpt: 'robustAtomMatch', args: 1, argName: 'false', 'If true, if an atom\'s counterpart from another file cannot be easily found, will scan over atoms to find the counterpart.');
cli.c(longOpt: 'crystRecord', args: 1, argName: 'false', 'If true, uses source CRYST1 record if available to replace or add to matched file.');
cli.sm(longOpt: 'simple', 'Set to use a simple match (ignores all flags except parallel, simply prints output of RMSD matching)');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 2) {
    return cli.usage();
}

if (options.d) {
    parseDeep = Boolean.parseBoolean(options.d);
}
if (options.v) {
    verbose = true;
}
if (options.hL) {
    headerLink = Boolean.parseBoolean(options.hL);
}
if (options.b) {
    bFactors = Boolean.parseBoolean(options.b);
}
if (options.sB) {
    ssBonds = Boolean.parseBoolean(options.sB);
}
if (options.s) {
    superpose = Boolean.parseBoolean(options.s);
}
if (options.a) {
    atomsUsed = Integer.parseInt(options.a);
    if (atomsUsed < 1 || atomsUsed > 4) {
        logger.warning(String.format(" Incorrect usage of atomsUsed %d: must be in range 1-4. Using default of 1.", atomsUsed));
        atomsUsed = 1;
    }
}
if (options.p) {
    parallel = Boolean.parseBoolean(options.p);
}
if (options.x) {
    suffix = options.x;
}
if (options.r) {
    robustMatch = Boolean.parseBoolean(options.r);
}

if (options.c) {
    fixCryst = Boolean.parseBoolean(options.c);
}

sourceFileSource = new File(arguments.get(0));
if (!sourceFileSource.exists()) {
    logger.info(" Source file list or directory does not exist.");
    return cli.usage();
}
matchFileSource = new File(arguments.get(1));
if (!matchFileSource.exists()) {
    logger.info(" Match file list or directory does not exist.");
    return cli.usage();
}

sourceStructures = new ArrayList<>();
if (sourceFileSource.isDirectory()) {
    if (parseDeep) {
        PDBFileFilter filter = new PDBFileFilter();
        File[] fileArray = sourceFileSource.listFiles();
        for (File file : fileArray) {
            if (filter.acceptDeep(file)) {
                sourceStructures.add(file);
            }
        }
    } else {
        sourceStructures.addAll(Arrays.asList(sourceFileSource.listFiles(new PDBFileFilter())));
    }
} else {
    BufferedReader bw = new BufferedReader(new FileReader(sourceFileSource));
    try {
        PDBFileFilter filter = new PDBFileFilter();
        String line = bw.readLine();
        while (line != null) {
            line = line.trim();
            File file = new File(line);
            if (file.exists() && file.isFile()) {
                if (parseDeep) {
                    if (filter.acceptDeep(file)) {
                        sourceStructures.add(file);
                    }
                } else {
                    if (filter.accept(file)) {
                        sourceStructures.add(file);
                    }
                }
            }
            line = bw.readLine();
        }
    } catch (IOException ex) {
        logger.severe(String.format(" Exception in reading source files %s", ex.toString()));
    }
    bw.close();
}
File[] sourceFileArray = sourceStructures.toArray(new File[sourceStructures.size()]);

matchStructures = new ArrayList<>();
if (matchFileSource.isDirectory()) {
    if (parseDeep) {
        PDBFileFilter filter = new PDBFileFilter();
        File[] fileArray = matchFileSource.listFiles();
        for (File file : fileArray) {
            if (filter.acceptDeep(file)) {
                matchStructures.add(file);
            }
        }
    } else {
        matchStructures.addAll(Arrays.asList(matchFileSource.listFiles(new PDBFileFilter())));
    }
} else {
    BufferedReader br = new BufferedReader(new FileReader(matchFileSource));
    try {
        PDBFileFilter filter = new PDBFileFilter();
        String line = br.readLine();
        while (line != null) {
            line = line.trim();
            File file = new File(line);
            if (file.exists() && file.isFile()) {
                if (parseDeep) {
                    if (filter.acceptDeep(file)) {
                        matchStructures.add(file);
                    }
                } else {
                    if (filter.accept(file)) {
                        matchStructures.add(file);
                    }
                }
            }
            line = br.readLine();
        }
    } catch (IOException ex) {
        logger.severe(String.format(" Exception in reading files to align %s", ex.toString()));
    }
    br.close();
}
File[] matchFileArray = matchStructures.toArray(new File[matchStructures.size()]);

if (options.sm) {
    SimplePDBMatcher simpleMatcher = new SimplePDBMatcher(matchFileArray, sourceFileArray);
    if (parallel) {
        simpleMatcher.matchParallel();
    } else {
        simpleMatcher.match();
    }
} else {
    fileMatcher = new PDBFileMatcher(sourceFileArray, matchFileArray);
    fileMatcher.setVerbose(verbose);
    fileMatcher.setParallel(parallel);
    fileMatcher.setHeaderLink(headerLink);
    fileMatcher.setFixBFactors(bFactors);
    fileMatcher.setFixSSBonds(ssBonds);
    fileMatcher.setSuperpose(superpose);
    fileMatcher.setAtomsUsed(atomsUsed);
    fileMatcher.setSuffix(suffix);
    fileMatcher.setRobustMatch(robustMatch);
    fileMatcher.setFixCryst(fixCryst);
    fileMatcher.match();
}
