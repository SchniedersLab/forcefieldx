
/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.MolecularAssembly;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.extended.TitrationESV;
import ffx.potential.extended.TitrationUtils;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.extended.ExtUtils.SB_Instance
import ffx.potential.utils.PotentialsUtils;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.logging.Level.FINE;
import static java.lang.String.format;
import org.apache.commons.lang.StringUtils;

Logger ffxlog = Logger.getLogger("ffx");

// finite-difference parameters
double lambda = 0.28;
double step = 0.0001;

// ESV discretization bias height
double biasMag = 1.0;
double pH;
boolean verbose = false;
int testOneIterations = 1;
boolean usePersistentIndex = true;
PotentialsUtils utils = new PotentialsUtils();
double cutoff;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testESVs [options]');
cli.h(longOpt:'help', 'Print this help message.');
cli.l(longOpt:'lambda', args:1, argName:'rng', 'Initial lambda value (for all ESVs).');
cli.t1(longOpt:'test1', 'Test 1: Lambda derivatives by finite difference.');
cli.t2(longOpt:'test2', 'Test 2: End state energies verification.');
cli.t3(longOpt:'test3', 'Test 3: Switching function and path smoothness.');
cli.i(longOpt:'iterations', args:1, argName:'1', 'Repeat FDs to verify threaded replicability.');
cli.v(longOpt:'verbose', 'Print out all the ForceFieldEnergy decompositions.');
cli.xyz(longOpt:'xyzIndex', 'Forego use of atom-indexing=PERSIST.');
cli.cut(longOpt:'cutoff', args:1, argName:'1000', 'Value of vdw-cutoff and pme-cutoff.');
cli.c(longOpt:'cryst', 'Use a faked crystal lattice (necessary for reciprocal space).');
cli.norecip('Use non-crystal inputs and do not include reciprocal space contributions.');

def options = cli.parse(args);
if (options.h) {
    return cli.usage();
}

boolean testVdw = true;
boolean testPermReal = true;
boolean testReciprocal = (Boolean) !(options.norecip);
boolean testPolarization = false;
boolean useCrystal = (Boolean) options.c || testReciprocal;

String esvFilename;         // 1x ESV (fully protonated)
String[] endStates;         // 4x Mutated Endpoints
String[] rlTokens = new String[2];
if (useCrystal) {
    esvFilename = "lys-lys-cryst.pdb";
    endStates = ["lyd-lyd-cryst", "lyd-lys-cryst", "lys-lyd-cryst", "lys-lys-cryst"];
    cutoff = 5.0;
} else {
    esvFilename = "lys-lys.pdb";
    endStates = ["lyd-lyd", "lyd-lys", "lys-lyd", "lys-lys"];
    cutoff = 1000.0;
}
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
boolean test2 = true;
boolean test3 = true;

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

if (!options.xyz) {
    usePersistentIndex = true;
    System.setProperty("atom-indexing", "PERSIST");
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

if (testReciprocal && cutoff > 20.0) {
    ffxlog.warning("Large cutoff for periodic system! Resetting to 10.0 Ang.");
    cutoff = 10.0;
}

// Active Potential
System.setProperty("forcefield", "AMOEBA_PROTEIN_2013");
System.setProperty("esvterm", "true");
System.setProperty("lambdaterm", "true");
System.setProperty("bondterm", "true");
System.setProperty("angleterm", "true");
System.setProperty("strbndterm", "true");
System.setProperty("ureyterm", "true");
System.setProperty("opbendterm", "true");
System.setProperty("torsionterm", "true");
System.setProperty("pitorsterm", "true");
System.setProperty("tortorterm", "true");
System.setProperty("improperterm", "true");

// Optional Potential
System.setProperty("vdwterm", String.valueOf(testVdw));                 // van der Waals
System.setProperty("esv-useVdw", String.valueOf(testVdw));
System.setProperty("mpoleterm", String.valueOf(testPermReal));          // permanent real space
System.setProperty("pme-qi", String.valueOf(testPermReal));
System.setProperty("esv-usePme", String.valueOf(testPermReal));
System.setProperty("recipterm", String.valueOf(testReciprocal));         // permanent reciprocal space
System.setProperty("polarizeterm", String.valueOf(testPolarization));   // polarization
if (testPolarization) {
    System.setProperty("polarization", "DIRECT");
} else {
    System.setProperty("polarization", "NONE");
}

// Inactive Potential
System.setProperty("gkterm", "false");
System.setProperty("restrainterm", "false");
System.setProperty("comrestrainterm", "false");
System.setProperty("lambda_torsions", "false");

// Potential Settings
System.setProperty("permanent-lambda-alpha","2.0");
System.setProperty("permanent-lambda-exponent","3.0");
System.setProperty("polarization-lambda-start","0.0");      // polarize on the whole range [0,1]
System.setProperty("polarization-lambda-exponent","0.0");   // polarization not softcored, only prefactored
System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization
System.setProperty("intramolecular-softcore", "true");
System.setProperty("intermolecular-softcore", "true");
System.setProperty("vdw-cutoff", String.valueOf(cutoff));
System.setProperty("ewald-cutoff", String.valueOf(cutoff));

// ESV Settings
System.setProperty("esv-propagation", "false");         // don't allow ESV particle to undergo dynamics
System.setProperty("esv-biasTerm", "true");             // include discretization and pH biases
System.setProperty("esv-scaleBonded", "true");          // include effects on bonded terms
System.setProperty("esv-backgroundBonded", "true");     // hook up BG bonded terms to FG node
System.setProperty("esv-scaleUnshared", "true");        // use multipole scaling in all cases (eliminates softcoring)

// Logging Settings
System.setProperty("ffe-combineBonded", "true");        // fold all bonded contributions into a single line
System.setProperty("ffe-decomposePme", "true");         // print the six components of PME: {real,recip}*{self,perm,ind}

// Open fully protonated and create TitrationESV objects.
MolecularAssembly mola = utils.open(esvFilename);
ForceFieldEnergy ffe = mola.getPotentialEnergy();
ExtendedSystem esvSystem = new ExtendedSystem(mola);

Polymer[] polymers = mola.getChains();
double temperature = 298.15;
int numESVs = rlTokens.length;
for (int i = 0; i < numESVs; i++) {
    Character chainID = rlTokens[i].charAt(0);
    int resNum = Integer.parseInt(rlTokens[i].substring(1));
    Optional<Residue> target = new Optional<>();
    for (Polymer p : polymers) {
        if (p.getChainID().equals(chainID)) {
            target = p.getResidues().stream()
                .filter {res -> res.getResidueNumber() == resNum}
                .findFirst();
            break;
        }
    }
    if (!target.isPresent()) {
        ffxlog.severe("Couldn't find target residue " + rlTokens[i]);
    }

    MultiResidue titrating = TitrationUtils.titrationFactory(mola, target.get());
    TitrationESV esv = new TitrationESV(esvSystem.getConfig(), titrating, pH, biasMag);
    esvSystem.addVariable(esv);
}

// Attach populated esvSystem to the potential.
ffe.attachExtendedSystem(esvSystem);

// Turn off checks for overlapping atoms, which is expected for lambda=0.
ffe.getCrystal().setSpecialPositionCutoff(0.0);

// Declare arrays and get coordinates.
n = ffe.getNumberOfVariables();
assert(n % 3 == 0);
int nAtoms = n.intdiv(3);
double[] xyz = new double[n];
double[] gradient = new double[n];
double[][] esvLambdaGrad = new double[numESVs][n];
double[][][] esvLambdaGradFD = new double[2][numESVs][n];
ffe.getCoordinates(xyz);

double width = 2.0 * step;      // default step = 0.0001
double tolerance = 1.0e-3;      // pass-fail margin

/* Try block to request crash-dump on error. */
try {
    /*******************************************************************************
     * Analytic Derivative vs Finite Difference: VdW,PermReal,PermRecip,Total
    *******************************************************************************/
    if (test1) {
        for (int iter = 0; iter < testOneIterations; iter++) {
            for (int i = 0; i < numESVs; i++) {
                esvSystem.setLambda(i, lambda);
            }
            for (int i = 0; i < numESVs; i++) {
                esvSystem.setLambda(i, 0.0);
                final double e0 = ffe.energyAndGradient(xyz,gradient);
                esvSystem.setLambda(i, 1.0);
                final double e1 = ffe.energyAndGradient(xyz,gradient);
                StringBuilder sb = new StringBuilder();
                sb.append(format(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f\n", i, e1, e0, e1-e0));

                esvSystem.setLambda(i, lambda - step);
                final double totalLow = ffe.energy(true, false);
                final double vdwLow = ffe.getVanDerWaalsEnergy();
                final double realLow = ffe.getPermanentRealSpaceEnergy();
                final double selfLow = ffe.getPermanentReciprocalSelfEnergy();
                final double recipLow = ffe.getPermanentReciprocalMpoleEnergy();

                esvSystem.setLambda(i, lambda + step);
                final double totalHigh = ffe.energy(true, false);
                final double vdwHigh = ffe.getVanDerWaalsEnergy();
                final double realHigh = ffe.getPermanentRealSpaceEnergy();
                final double selfHigh = ffe.getPermanentReciprocalSelfEnergy();
                final double recipHigh = ffe.getPermanentReciprocalMpoleEnergy();

                esvSystem.setLambda(i, lambda);
                final double totalCenter = ffe.energy(true, true);
                final double vdwCenter = ffe.getVanDerWaalsEnergy();
                final double realCenter = ffe.getPermanentRealSpaceEnergy();
                final double selfCenter = ffe.getPermanentReciprocalSelfEnergy();
                final double recipCenter = ffe.getPermanentReciprocalMpoleEnergy();

                final double vdwAnalytic = esvSystem.getdVdwdL(i);
                final double vdwNumeric = (vdwHigh - vdwLow).div(2 * step);
                final double vdwError = Math.abs(vdwNumeric - vdwAnalytic);
                final double realAnalytic = esvSystem.getdPermRealdL(i);
                final double realNumeric = (realHigh - realLow).div(2 * step);
                final double realError = Math.abs(realNumeric - realAnalytic);
                final double selfAnalytic = esvSystem.getdPermRecipSelf_dL(i);
                final double selfNumeric = (selfHigh - selfLow).div(2 * step);
                final double selfError = Math.abs(selfNumeric - selfAnalytic);
                final double recipAnalytic = esvSystem.getdPermRecipMpole_dL(i);
                final double recipNumeric = (recipHigh - recipLow).div(2 * step);
                final double recipError = Math.abs(recipNumeric - recipAnalytic);
                final double totalAnalytic = esvSystem.getdEdL(i, true);
                final double totalNumeric = (totalHigh - totalLow).div(2 * step);
                final double totalError = Math.abs(totalNumeric - totalAnalytic);

                String esvName = esvSystem.getEsv(i).getName();
                sb.append(format("\n  Finite Difference Test: %s", esvName));
                sb.append(format("\n *************************"));
                sb.append(format(StringUtils.repeat("*", esvName.length()))).append("*");
                sb.append(format("   %10s   %10s    %10s    %10s    %10s\n", "vdw", "permReal", "permSelf", "permRecip", "total"));
                sb.append(format(" Numeric  Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f\n",
                        lambda, vdwNumeric, realNumeric, selfNumeric, recipNumeric, totalNumeric));
                sb.append(format(" Analytic Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f\n",
                        lambda, vdwAnalytic, realAnalytic, selfAnalytic, recipAnalytic, totalAnalytic));
                sb.append(format(" Error:                               %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f\n",
                        vdwError, realError, selfError, recipError, totalError));
                ffxlog.info(sb.toString());
            }
        }
    }

    /*******************************************************************************
     * End-State Verification
    Verify that a lys-lys system with two ESVs can exactly reproduce the VdW
    energy yielded by vanilla energy() calls on mutated PDB files.
    *******************************************************************************/
    if (test2) {
        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 0.0);
        ffe.energy(true, false);
        final double esvVdw00 = ffe.getVanDerWaalsEnergy();
        final double esvPme00 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp00 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 1.0);
        ffe.energy(true, false);
        final double esvVdw01 = ffe.getVanDerWaalsEnergy();
        final double esvPme01 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp01 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 0.0);
        ffe.energy(true, false);
        final double esvVdw10 = ffe.getVanDerWaalsEnergy();
        final double esvPme10 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp10 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 1.0);
        ffe.energy(true, false);
        final double esvVdw11 = ffe.getVanDerWaalsEnergy();
        final double esvPme11 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp11 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        StringBuilder ext = new StringBuilder();
        if (testVdw) {
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "0-0", esvVdw00));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "0-1", esvVdw01));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "1-0", esvVdw10));
            ext.append(format("%-9s %-7s %10.6f\n", "VanWaals", "1-1", esvVdw11));
        }
        if (testPermReal) {
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "0-0", esvPme00));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "0-1", esvPme01));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "1-0", esvPme10));
            ext.append(format("%-9s %-7s %10.6f\n", "PermReal", "1-1", esvPme11));
        }
        if (testReciprocal) {
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "0-0", esvRcp00));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "0-1", esvRcp01));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "1-0", esvRcp10));
            ext.append(format("%-9s %-7s %10.6f\n", "PermRecip", "1-1", esvRcp11));
        }

        StringBuilder vdw = new StringBuilder();
        StringBuilder pme = new StringBuilder();
        StringBuilder rcp = new StringBuilder();
        for (String state : endStates) {
            MolecularAssembly vanilla = utils.openQuietly(state + ".pdb");
            ForceFieldEnergy vanillaFFE = vanilla.getPotentialEnergy();
            vanillaFFE.energy(true, false);
            if (testVdw) {
                final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
                vdw.append(format("%-9s %-7s %10.6f\n", "VanWaals", state, vanWaals));
            }
            if (testPermReal) {
                final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
                pme.append(format("%-9s %-7s %10.6f\n", "PermReal", state, permReal));
            }
            if (testReciprocal) {
                final double reciprocal = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
                rcp.append(format("%-9s %-7s %10.6f\n", "PermRecip", state, reciprocal));
            }
            utils.close(vanilla);
        }
        StringBuilder vanilla = new StringBuilder();
        if (testVdw) vanilla.append(vdw.toString());
        if (testPermReal) vanilla.append(pme.toString());
        if (testReciprocal) vanilla.append(rcp.toString());
        String[] vanLines = vanilla.toString().split("\\n");
        String[] extLines = ext.toString().split("\\n");
        StringBuilder sb = new StringBuilder();
        sb.append(format("  Two-site End State Analysis: \n"));
        sb.append(format(" ****************************** \n"));
        sb.append(format(" %-30s     %-30s\n","Extended System (1 Assembly)","Mutated Files (4 Assemblies)"));
        for (int i = 0; i < extLines.length; i++) {
            sb.append(format(" %-30s     %-30s\n", extLines[i], vanLines[i]));
        }
        ffxlog.info(sb.toString());
    }

    /*******************************************************************************
     * Switching and Smoothness
    Numerically ensure that the VdW energy and lambda derivatives are smooth all 
    along both ESV coordinates in the dilysine system.
    *******************************************************************************/
    if (test3) {
        double[][] totalEnergies, totalDerivsA, totalDerivsB;
        double[][] vdwEnergies, vdwDerivsA, vdwDerivsB;
        double[][] permRealEnergies, permRealDerivsA, permRealDerivsB;
        double[][] permRecipEnergies, permRecipDerivsA, permRecipDerivsB;
        totalEnergies = new double[11][11];
        totalDerivsA = new double[11][11];
        totalDerivsB = new double[11][11];
        if (testVdw) {
            vdwEnergies = new double[11][11];
            vdwDerivsA = new double[11][11];
            vdwDerivsB = new double[11][11];
        }
        if (testPermReal) {
            permRealEnergies = new double[11][11];
            permRealDerivsA = new double[11][11];
            permRealDerivsB = new double[11][11];
        }
        if (testReciprocal) {
            permRecipEnergies = new double[11][11];
            permRecipDerivsA = new double[11][11];
            permRecipDerivsB = new double[11][11];
        }
        for (int idxA = 0; idxA <= 10; idxA++) {
            double evA = idxA.div(10.0);
            for (int idxB = 0; idxB <= 10; idxB++) {
                final double evB = idxB.div(10.0);
                esvSystem.setLambda('A', evA);
                esvSystem.setLambda('B', evB);
                final double totalEnergy = ffe.energy(true, false);
                totalEnergies[idxA][idxB] = totalEnergy;
                totalDerivsA[idxA][idxB] = esvSystem.getdEdL(0);
                totalDerivsB[idxA][idxB] = esvSystem.getdEdL(1);
                if (testVdw) {
                    vdwEnergies[idxA][idxB] = ffe.getVanDerWaalsEnergy();
                    vdwDerivsA[idxA][idxB] = esvSystem.getdVdwdL(0);
                    vdwDerivsB[idxA][idxB] = esvSystem.getdVdwdL(1);
                }
                if (testPermReal) {
                    permRealEnergies[idxA][idxB] = ffe.getPermanentRealSpaceEnergy();
                    permRealDerivsA[idxA][idxB] = esvSystem.getdPermRealdL(0);
                    permRealDerivsB[idxA][idxB] = esvSystem.getdPermRealdL(1);
                }
                if (testReciprocal) {
                    permRecipEnergies[idxA][idxB] = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();
                    permRecipDerivsA[idxA][idxB] = esvSystem.getdPermRecipSelf_dL(0) + esvSystem.getdPermRecipMpole_dL(0);
                    permRecipDerivsB[idxA][idxB] = esvSystem.getdPermRecipSelf_dL(1) + esvSystem.getdPermRecipMpole_dL(1);
                }
            }
        }
        
        /* Printing loop */
        ffxlog.info(format("  Smoothness Verification: Total "));
        ffxlog.info(format(" ******************************** "));
        printAsTable(totalEnergies, "U_Total");
        printAsTable(totalDerivsA, "dU_dEsvA");
        printAsTable(totalDerivsB, "dU_dEsvB");
        if (testVdw && verbose) {
            ffxlog.info(format("  Smoothness Verification: VdW "));
            ffxlog.info(format(" ****************************** "));
            printAsTable(vdwEnergies, "vanWaals");
            printAsTable(vdwDerivsA, "dVdw_dA");
            printAsTable(vdwDerivsB, "dVdw_dB");
        }
        if (testPermReal && verbose) {
            ffxlog.info(format("  Smoothness Verification: PermReal "));
            ffxlog.info(format(" *********************************** "));
            printAsTable(permRealEnergies, "permReal");
            printAsTable(permRealDerivsA, "dPRealdA");
            printAsTable(permRealDerivsB, "dPRealdB");
        }
        if (testReciprocal && verbose) {
            ffxlog.info(format("  Smoothness Verification: PermRecip "));
            ffxlog.info(format(" ************************************ "));
            printAsTable(permRecipEnergies, "permRcp");
            printAsTable(permRecipDerivsA, "dPRcp_dA");
            printAsTable(permRecipDerivsB, "dPRcp_dB");
        }
    }   // Test3
} catch (RuntimeException ex) {
    /* Dump system info on any unhandled exception. */
    esvSystem.crashDump();
    throw ex;
}

def void printAsTable(double[][] values, String title) {
    if (title == null) title = "Title?";
    StringBuilder sb = new StringBuilder();
    sb.append(format(" %8s %8.1f %8.1f %8.1f %8.1f %8s %8s %8s %8.1f %8.1f %8.1f %8.1f", 
            title, 0.0, 0.1, 0.2, 0.3, "<   ", " EvA  ", ">   ", 0.7, 0.8, 0.9, 1.0));
    for (int idxA = 0; idxA <= 10; idxA++) {
        double evA = idxA.div(10.0);
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
        for (int idxB = 0; idxB <= 10; idxB++) {
            double evB = idxB.div(10.0);
            String value = format("%8.4f", values[idxA][idxB]);
            int max = (value.length() < 8) ? value.length() : 8;
            sb.append(format(" %8s", value.substring(0,max)));
        }
    }
    sb.append(format("\n"));
    Logger.getLogger("ffx").info(sb.toString());
}

