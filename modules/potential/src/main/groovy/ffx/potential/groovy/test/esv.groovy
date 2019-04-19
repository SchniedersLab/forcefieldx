//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import java.util.logging.Logger
import static java.lang.String.format

import org.apache.commons.lang3.StringUtils

import groovy.cli.picocli.CliBuilder

import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.PotentialComponent
import ffx.potential.extended.ExtendedSystem
import ffx.potential.utils.PotentialsUtils
import static ffx.potential.PotentialComponent.Bias
import static ffx.potential.PotentialComponent.Bonded
import static ffx.potential.PotentialComponent.InducedRealSpace
import static ffx.potential.PotentialComponent.InducedReciprocal
import static ffx.potential.PotentialComponent.InducedSelf
import static ffx.potential.PotentialComponent.PermanentRealSpace
import static ffx.potential.PotentialComponent.PermanentReciprocal
import static ffx.potential.PotentialComponent.PermanentSelf
import static ffx.potential.PotentialComponent.Topology
import static ffx.potential.PotentialComponent.VanDerWaals

// finite-difference parameters
double lambda = 0.5;
double step = 0.0001;

// ESV discretization bias height
double biasMag = 1.0;
double pH = 7.4;
boolean verbose = false;
int testOneIterations = 1;
boolean usePersistentIndex = true;
PotentialsUtils utils = new PotentialsUtils();
double cutoff;
boolean tablesToFile = false;
int tableDimensions = 10;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage: ' ffxc test.esv [options]');
cli.h(longOpt: 'help', 'Print this help message.');
cli.l(longOpt: 'lambda', args: 1, argName: 'random', 'Initial lambda value (for all ESVs).');
cli.t1(longOpt: 'test1', 'Test 1: Lambda derivatives by finite difference.');
cli.t2(longOpt: 'test2', 'Test 2: End state energies verification.');
cli.t3(longOpt: 'test3', 'Test 3: Switching function and path smoothness.');
cli.i(longOpt: 'iterations', args: 1, argName: '1', 'Repeat FDs to verify threaded replicability.');
cli.v(longOpt: 'verbose', 'Print out all the ForceFieldEnergy decompositions.');
cli.cut(longOpt: 'cutoff', args: 1, argName: '1000', 'Value of vdw-cutoff and pme-cutoff.');
cli.nocryst('Test a vacuum system.');
cli.pH(args: 1, argName: '7.4', 'Constant system pH.');
cli.ttf(longOpt: 'tables-to-file', args: 1, argName: '100', 'Export Utot and dUdEsv to files with this granularity (third arg to MATLAB\'s linspace).');

def options = cli.parse(args);
if (options.h) {
    return cli.usage();
}

String filename = args;

boolean useCrystal = (Boolean) !(options.nocryst);
boolean testVdw = true;
boolean testPerm = true;
boolean testInduced = true;
boolean testPolarization = true;

if (options.ttf) {
    tablesToFile = true;
    tableDimensions = Integer.parseInt(options.ttf);
}

String esvFilename;         // 1x ESV (fully protonated)
String[] endStates;         // 4x Mutated Endpoints

if (useCrystal) {
    esvFilename = "lys-lys-cryst.pdb";
    endStates = ["lyd-lyd-cryst", "lyd-lys-cryst", "lys-lyd-cryst", "lys-lys-cryst"];
    cutoff = 5.0;
} else {
    esvFilename = "lys-lys.pdb";
    endStates = ["lyd-lyd", "lyd-lys", "lys-lyd", "lys-lys"];
    cutoff = 1000.0;
}

String[] rlTokens = new String[2];
rlTokens[0] = "A2";
rlTokens[1] = "A4";

/* Auto-find the examples directory */
String filePrefix = "";
if (!(new File("lys-lys.pdb")).exists()) {
    logger.info("Couldn't find lys-lys.pdb, trying examples/lys-lys.pdb...");
    String fn = "examples/lys-lys.pdb";
    if (!(new File(fn)).exists()) {
        logger.severe("Couldn't find input files. Try examples directory.");
    } else {
        filePrefix = "examples/";
        esvFilename = filePrefix + esvFilename;
        for (int i = 0; i < endStates.length; i++) {
            endStates[i] = filePrefix + endStates[i];
        }
    }
}
logger.info(format(" Loading Files: %s, %s", esvFilename, Arrays.toString(endStates)));

/* Default Tests */
boolean test1 = true;
boolean test2 = false;
boolean test3 = false;

if (options.t1 || options.t2 || options.t3) {
    test1 = false; test2 = false; test3 = false;
    if (options.t1) {
        test1 = true;
    }
    if (options.t2) {
        test2 = true;
    }
    if (options.t3) {
        test3 = true;
    }
}

if (options.l) {
    lambda = Double.parseDouble(options.l);
}

if (options.i) {
    testOneIterations = Integer.parseInt(options.i);
}

if (options.cut) {
    cutoff = Double.parseDouble(options.cut);
}

if (options.v) {
    verbose = true;
}

if (options.pH) {
    pH = Double.parseDouble(options.pH);
} else {
    pH = 7.4;
}

// Active Potential
System.setProperty("forcefield", "AMOEBA_PROTEIN_2013");
System.setProperty("esvterm", "true");
System.setProperty("lambdaterm", "false");
System.setProperty("pme-qi", "true");

// Optional Potential
System.setProperty("vdwterm", String.valueOf(testVdw));                 // van der Waals
System.setProperty("mpoleterm", String.valueOf(testPerm));          // permanent real space
System.setProperty("esv-vdw", String.valueOf(testVdw));
System.setProperty("esv-pme", String.valueOf(testPerm));
System.setProperty("esv-propagation", "false");         // don't allow ESV particle to undergo dynamics
System.setProperty("esv-biasTerm", "true");             // include discretization and pH biases
System.setProperty("esv-scaleBonded", "true");          // include effects on bonded terms
System.setProperty("esv-backgroundBonded", "true");     // hook up BG bonded terms to FG node
System.setProperty("esv-scaleUnshared", "true");        // use multipole scaling in all cases (eliminates softcoring)

// Inactive Potential
System.setProperty("gkterm", "false");
System.setProperty("restrainterm", "false");
System.setProperty("comrestrainterm", "false");
System.setProperty("torsion-lambdaterm", "false");

// Potential Settings
System.setProperty("permanent-lambda-alpha", "2.0");
System.setProperty("permanent-lambda-exponent", "3.0");
System.setProperty("polarization-lambda-start", "0.0");      // polarize on the whole range [0,1]
System.setProperty("polarization-lambda-exponent", "0.0");   // polarization not softcored, only prefactored
System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization
System.setProperty("intramolecular-softcore", "true");
System.setProperty("intermolecular-softcore", "true");
System.setProperty("vdw-cutoff", String.valueOf(cutoff));
System.setProperty("ewald-cutoff", String.valueOf(cutoff));

// Logging Settings
System.setProperty("ffe-combineBonded", "true");        // fold all bonded contributions into a single line
System.setProperty("ffe-decomposePme", "true");         // print the six components of PME: {real,recip}*{self,perm,ind}

// Open fully protonated and create TitrationESV objects.
MolecularAssembly mola = utils.open(esvFilename);
ForceFieldEnergy ffe = mola.getPotentialEnergy();
ExtendedSystem esvSystem = new ExtendedSystem(mola);
esvSystem.setConstantPh(pH);
esvSystem.populate(rlTokens);
int numESVs = esvSystem.size();
// Attach populated esvSystem to the potential.
ffe.attachExtendedSystem(esvSystem);
// Turn off checks for overlapping atoms, which is expected for lambda=0.
ffe.getCrystal().setSpecialPositionCutoff(0.0);

/* Try block to request crash-dump on error. */
try {
    /*******************************************************************************
     * Analytic Derivative vs Finite Difference: VdW,PermReal,PermRecip,Total
     *******************************************************************************/
    if (test1) {
        for (int iter = 0; iter < testOneIterations; iter++) {
            for (int i = 0; i < numESVs; i++) {
                for (int k = 0; k < numESVs; k++) {
                    esvSystem.setLambda(k, lambda);
                }

                esvSystem.setLambda(i, lambda - step);
                ffe.energy(true, false);
                final double vdwLow = ffe.getEnergyComponent(VanDerWaals);
                final double bondedLow = ffe.getEnergyComponent(Bonded);
                final double biasLow = ffe.getEnergyComponent(Bias);
                final double permRealLow = ffe.getEnergyComponent(PermanentRealSpace);
                final double permSelfLow = ffe.getEnergyComponent(PermanentSelf);
                final double permRecipLow = ffe.getEnergyComponent(PermanentReciprocal);
                final double indRealLow = ffe.getEnergyComponent(InducedRealSpace);
                final double indSelfLow = ffe.getEnergyComponent(InducedSelf);
                final double indRecipLow = ffe.getEnergyComponent(InducedReciprocal);
                final double totalLow = ffe.getEnergyComponent(Topology);

                esvSystem.setLambda(i, lambda + step);
                ffe.energy(true, false);
                final double vdwHigh = ffe.getEnergyComponent(VanDerWaals);
                final double bondedHigh = ffe.getEnergyComponent(Bonded);
                final double biasHigh = ffe.getEnergyComponent(Bias);
                final double permRealHigh = ffe.getEnergyComponent(PermanentRealSpace);
                final double permSelfHigh = ffe.getEnergyComponent(PermanentSelf);
                final double permRecipHigh = ffe.getEnergyComponent(PermanentReciprocal);
                final double indRealHigh = ffe.getEnergyComponent(InducedRealSpace);
                final double indSelfHigh = ffe.getEnergyComponent(InducedSelf);
                final double indRecipHigh = ffe.getEnergyComponent(InducedReciprocal);
                final double totalHigh = ffe.getEnergyComponent(Topology);

                // Get analytic derivatives from the center.
                esvSystem.setLambda(i, lambda);
                ffe.energy(true, false);
                final double vdwAna = esvSystem.getDerivativeComponent(VanDerWaals, i);
                final double bondedAna = esvSystem.getDerivativeComponent(Bonded, i);
                final double biasAna = esvSystem.getDerivativeComponent(Bias, i);
                final double permRealAna = esvSystem.getDerivativeComponent(PermanentRealSpace, i);
                final double permSelfAna = esvSystem.getDerivativeComponent(PermanentSelf, i);
                final double permRecipAna = esvSystem.getDerivativeComponent(PermanentReciprocal, i);
                final double indRealAna = esvSystem.getDerivativeComponent(InducedRealSpace, i);
                final double indSelfAna = esvSystem.getDerivativeComponent(InducedSelf, i);
                final double indRecipAna = esvSystem.getDerivativeComponent(InducedReciprocal, i);
                final double totalAna = esvSystem.getDerivativeComponent(Topology, i);

                // Calculate numeric derivatives and error.
                final double vdwNum = (vdwHigh - vdwLow) / (2 * step);
                final double vdwErr = Math.abs(vdwNum - vdwAna);
                final double bondedNum = (bondedHigh - bondedLow) / (2 * step);
                final double bondedErr = Math.abs(bondedNum - bondedAna);
                final double biasNum = (biasHigh - biasLow) / (2 * step);
                final double biasErr = Math.abs(biasNum - biasAna);
                final double permRealNum = (permRealHigh - permRealLow) / (2 * step);
                final double permRealErr = Math.abs(permRealNum - permRealAna);
                final double permSelfNum = (permSelfHigh - permSelfLow) / (2 * step);
                final double permSelfErr = Math.abs(permSelfNum - permSelfAna);
                final double permRecipNum = (permRecipHigh - permRecipLow) / (2 * step);
                final double permRecipErr = Math.abs(permRecipNum - permRecipAna);
                final double indRealNum = (indRealHigh - indRealLow) / (2 * step);
                final double indRealErr = Math.abs(indRealNum - indRealAna);
                final double indSelfNum = (indSelfHigh - indSelfLow) / (2 * step);
                final double indSelfErr = Math.abs(indSelfNum - indSelfAna);
                final double indRecipNum = (indRecipHigh - indRecipLow) / (2 * step);
                final double indRecipErr = Math.abs(indRecipNum - indRecipAna);
                final double totalNum = (totalHigh - totalLow) / (2 * step);
                final double totalErr = Math.abs(totalNum - totalAna);

                String esvName = esvSystem.getEsv(i).getName();
                final int columns = 10;
                final String headFormat = " %9s";
                final String dataFormat = " %+9.4f";
                logger.info(format(" %-31s", format(" Finite Difference Test: %s", esvName)));
                logger.info(format(" %-31s %s",
                        StringUtils.repeat("*", esvName.length() + 26),
                        format(new String(new char[columns]).replace('\0', headFormat),
                                "vdw ", "bonded", "bias", "permReal", "permSelf", "permRecip", "indReal", "indSelf", "indRecip", "total")));
                logger.info(format(" %-31s %s",
                        format("Numeric  Derivatives (@L %4.2f):", lambda),
                        format(new String(new char[columns]).replace('\0', dataFormat),
                                vdwNum, bondedNum, biasNum, permRealNum, permSelfNum, permRecipNum, indRealNum, indSelfNum, indRecipNum, totalNum)));
                logger.info(format(" %-31s %s",
                        format("Analytic Derivatives (@L %4.2f):", lambda),
                        format(new String(new char[columns]).replace('\0', dataFormat),
                                vdwAna, bondedAna, biasAna, permRealAna, permSelfAna, permRecipAna, indRealAna, indSelfAna, indRecipAna, totalAna)));
                logger.info(format(" %-31s %s",
                        "Error:",
                        format(new String(new char[columns]).replace('\0', dataFormat),
                                vdwErr, bondedErr, biasErr, permRealErr, permSelfErr, permRecipErr, indRealErr, indSelfErr, indRecipErr, totalErr)));
            }
        }
    }

    /*******************************************************************************
     * End-State Verification
     Verify that a lys-lys system with two ESVs can exactly reproduce the VdW
     energy yielded by vanilla energy() calls on mutated PDB files.
     *******************************************************************************/
    if (test2) {
        esvSystem.setLambda((char) 'A', 0.0);
        esvSystem.setLambda((char) 'B', 0.0);
        final double esvTot00 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvVdw00 = ffe.getVanDerWaalsEnergy();
        final double esvPme00 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp00 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda((char) 'A', 0.0);
        esvSystem.setLambda((char) 'B', 1.0);
        final double esvTot01 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvVdw01 = ffe.getVanDerWaalsEnergy();
        final double esvPme01 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp01 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda((char) 'A', 1.0);
        esvSystem.setLambda((char) 'B', 0.0);
        final double esvTot10 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvVdw10 = ffe.getVanDerWaalsEnergy();
        final double esvPme10 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp10 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda((char) 'A', 1.0);
        esvSystem.setLambda((char) 'B', 1.0);
        final double esvTot11 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvVdw11 = ffe.getVanDerWaalsEnergy();
        final double esvPme11 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp11 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        StringBuilder ext = new StringBuilder();
        ext.append(format("%-9s %-7s %10.6f\n", "Total", "0-0", esvTot00));
        ext.append(format("%-9s %-7s %10.6f\n", "Total", "0-1", esvTot01));
        ext.append(format("%-9s %-7s %10.6f\n", "Total", "1-0", esvTot10));
        ext.append(format("%-9s %-7s %10.6f\n", "Total", "1-1", esvTot11));
        if (testVdw) {
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "0-0", esvVdw00));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "0-1", esvVdw01));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "1-0", esvVdw10));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "1-1", esvVdw11));
        }
        if (testPerm) {
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "0-0", esvPme00));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "0-1", esvPme01));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "1-0", esvPme10));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "1-1", esvPme11));
        }
        if (testInduced) {
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "0-0", esvRcp00));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "0-1", esvRcp01));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "1-0", esvRcp10));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "1-1", esvRcp11));
        }

        StringBuilder tot = new StringBuilder();
        StringBuilder vdw = new StringBuilder();
        StringBuilder pme = new StringBuilder();
        StringBuilder rcp = new StringBuilder();
        for (String state : endStates) {
            MolecularAssembly vanilla = utils.openQuietly(state + ".pdb");
            ForceFieldEnergy vanillaFFE = vanilla.getPotentialEnergy();
            final double Utot = vanillaFFE.energy(true, false);
            tot.append(format("%-9s %-7s %10.6f\n", "Total", state, Utot));
            if (testVdw) {
                final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
                vdw.append(format("%-9s %-7s %10.6f\n", "VanWaals", state, vanWaals));
            }
            if (testPerm) {
                final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
                pme.append(format("%-9s %-7s %10.6f\n", "PermReal", state, permReal));
            }
            if (testInduced) {
                final double reciprocal = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
                rcp.append(format("%-9s %-7s %10.6f\n", "PermRecip", state, reciprocal));
            }
            utils.close(vanilla);
        }
        StringBuilder vanilla = new StringBuilder();
        vanilla.append(tot.toString());
        if (testVdw) vanilla.append(vdw.toString());
        if (testPerm) vanilla.append(pme.toString());
        if (testInduced) vanilla.append(rcp.toString());
        String[] vanLines = vanilla.toString().split("\\n");
        String[] extLines = ext.toString().split("\\n");
        StringBuilder sb = new StringBuilder();
        sb.append(format("  Two-site End State Analysis: \n"));
        sb.append(format(" ****************************** \n"));
        sb.append(format(" %-30s     %-30s\n", "Extended System (1 Assembly)", "Mutated Files (4 Assemblies)"));
        for (int i = 0; i < extLines.length; i++) {
            sb.append(format(" %-30s     %-30s\n", extLines[i], vanLines[i]));
        }
        logger.info(sb.toString());
    }

    /*******************************************************************************
     * Switching and Smoothness
     Numerically ensure that the VdW energy and lambda derivatives are smooth all
     along both ESV coordinates in the dilysine system.
     *******************************************************************************/
    // TODO Make this loop more granular for surface plot.
    if (test3) {
        final int dimPlusOne = tableDimensions + 1;
        double[][] totalEnergies, totalDerivsA, totalDerivsB;
        double[][] vdwEnergies, vdwDerivsA, vdwDerivsB;
        double[][] permEnergies, permDerivsA, permDerivsB;
        double[][] polEnergies, polDerivsA, polDerivsB;
        totalEnergies = new double[dimPlusOne][dimPlusOne];
        totalDerivsA = new double[dimPlusOne][dimPlusOne];
        totalDerivsB = new double[dimPlusOne][dimPlusOne];
        if (testVdw) {
            vdwEnergies = new double[dimPlusOne][dimPlusOne];
            vdwDerivsA = new double[dimPlusOne][dimPlusOne];
            vdwDerivsB = new double[dimPlusOne][dimPlusOne];
        }
        if (testPerm) {
            permEnergies = new double[dimPlusOne][dimPlusOne];
            permDerivsA = new double[dimPlusOne][dimPlusOne];
            permDerivsB = new double[dimPlusOne][dimPlusOne];
        }
        if (testInduced) {
            polEnergies = new double[dimPlusOne][dimPlusOne];
            polDerivsA = new double[dimPlusOne][dimPlusOne];
            polDerivsB = new double[dimPlusOne][dimPlusOne];
        }
        logger.info(format(" %30s %4d/%d", "Progress, outer loop:", 0, tableDimensions));
        for (int idxA = 0; idxA <= tableDimensions; idxA++) {
            double evA = idxA.div(tableDimensions);
            logger.info(format(" %30s %4d/%d", "", idxA, tableDimensions));
            for (int idxB = 0; idxB <= tableDimensions; idxB++) {
                final double evB = idxB.div(tableDimensions);
                esvSystem.setLambda((char) 'A', evA);
                esvSystem.setLambda((char) 'B', evB);
                final double totalEnergy = ffe.energy(true, false);
                totalEnergies[idxA][idxB] = totalEnergy;
                totalDerivsA[idxA][idxB] = esvSystem.getDerivative(0);
                totalDerivsB[idxA][idxB] = esvSystem.getDerivative(1);
                if (testVdw) {
                    vdwEnergies[idxA][idxB] = ffe.getVanDerWaalsEnergy();
                    vdwDerivsA[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.VanDerWaals, 0);
                    vdwDerivsB[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.VanDerWaals, 1);
                }
                if (testPerm) {
                    permEnergies[idxA][idxB] = ffe.getPermanentMultipoleEnergy();
                    permDerivsA[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.Permanent, 0);
                    permDerivsB[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.Permanent, 1);
                }
                if (testInduced) {
                    polEnergies[idxA][idxB] = ffe.getPolarizationEnergy();
                    polDerivsA[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.Induced, 0);
                    polDerivsB[idxA][idxB] = esvSystem.getDerivativeComponent(PotentialComponent.Induced, 1);
                }
            }
        }

        /* Printing loop */
        logger.info(format("  Smoothness Verification: Total "));
        logger.info(format(" ******************************** "));
        if (tablesToFile) {
            tableToFile(totalEnergies, "U_Total", new File("esv.Utot"));
            tableToFile(totalDerivsA, "dU_dEsvA", new File("esv.dUdEvA"));
            tableToFile(totalDerivsB, "dU_dEsvB", new File("esv.dUdEvB"));
        } else {
            printAsTable(totalEnergies, "U_Total");
            printAsTable(totalDerivsA, "dU_dEsvA");
            printAsTable(totalDerivsB, "dU_dEsvB");
        }
        if (testVdw && verbose) {
            logger.info(format("  Smoothness Verification: van der Waals "));
            logger.info(format(" ****************************** "));
            printAsTable(vdwEnergies, "van der Waals");
            printAsTable(vdwDerivsA, "dVdw_dA");
            printAsTable(vdwDerivsB, "dVdw_dB");
        }
        if (testPerm && verbose) {
            logger.info(format("  Smoothness Verification: Permanent "));
            logger.info(format(" *********************************** "));
            printAsTable(permEnergies, "Permanent");
            printAsTable(permDerivsA, "dPerm_dA");
            printAsTable(permDerivsB, "dPerm_dB");
        }
        if (testInduced && verbose) {
            logger.info(format("  Smoothness Verification: Polarization "));
            logger.info(format(" ************************************ "));
            printAsTable(polEnergies, "Polarization");
            printAsTable(polDerivsA, "dPol_dA");
            printAsTable(polDerivsB, "dPol_dB");
        }
    }   // Test3
} catch (RuntimeException ex) {
    /* Dump system info on any unhandled exception. */
    esvSystem.crashDump();
    throw ex;
}

def void tableToFile(double[][] values, String title, File out) {
    printAsTable(values, title, out);
}

def void printAsTable(double[][] values, String title) {
    printAsTable(values, title, null);
}

def void printAsTable(double[][] values, String title, File out) {
    boolean tablesToFile = (out != null);
    int tableDimensions = values.length - 1;
    if (title == null) title = "Title?";
    StringBuilder sb = new StringBuilder();
    if (!tablesToFile) {
        sb.append(format(" %8s %8.1f %8.1f %8.1f %8.1f %8s %8s %8s %8.1f %8.1f %8.1f %8.1f",
                title, 0.0, 0.1, 0.2, 0.3, "<   ", " EvA  ", ">   ", 0.7, 0.8, 0.9, 1.0));
    }
    for (int idxA = 0; idxA <= tableDimensions; idxA++) {
        double evA = idxA.div(tableDimensions);
        if (tablesToFile) {
            sb.append(format("\n"));
        } else {
            switch (evA) {
                case 0.4:
                    sb.append(format("\n %8s", " ^ "));
                    break;
                case 0.5:
                    sb.append(format("\n %8s", "EvB"));
                    break;
                case 0.6:
                    sb.append(format("\n %8s", " v "));
                    break;
                default:
                    sb.append(format("\n %8.1f", evA));
            }
        }
        for (int idxB = 0; idxB <= tableDimensions; idxB++) {
            double evB = idxB.div(tableDimensions);
            if (tablesToFile) {
                sb.append(format(",%.4f", values[idxA][idxB]));
            } else {
                String value = format("%8.4f", values[idxA][idxB]);
                int max = (value.length() < 8) ? value.length() : 8;
                sb.append(format(" %8s", value.substring(0, max)));
            }
        }
    }
    sb.append(format("\n"));
    if (tablesToFile) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
            bw.write(sb.toString());
            bw.close();
        } catch (IOException ex) {
        }
    } else {
        Logger.getLogger("ffx").info(sb.toString());
    }
}

