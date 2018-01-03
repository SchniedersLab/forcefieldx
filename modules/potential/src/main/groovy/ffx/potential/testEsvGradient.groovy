
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

// TEST LAMBDA GRADIENT

import java.util.logging.Level
import java.util.logging.Logger
import static java.lang.String.format

import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.extended.ExtendedSystem
import ffx.potential.extended.TitrationESV
import ffx.potential.extended.TitrationUtils

// finite-difference parameters
double lambda = 0.5;
double step = 0.0001;

// ESV discretization bias height
double biasMag = 1.0;
double pH;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testEsvGradient [options] <XYZ|PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.l(longOpt:'lambda', args:1, argName:'1.0', 'Initial lambda value (for all ESVs).');
cli.dx(longOpt:'stepSize', args:1, argName:'0.0001', 'Finite difference step size.');
cli.rl(longOpt:'resList', required:true, args:1, 'List of residues to titrate (e.g. A4.A8.B2.B34)');
cli.bm(longOpt:'biasMag', args:1, argName:'1.0', 'ESV discretization bias height.');
cli.st(longOpt:'singleThread', args:0, 'Limit PJ concurrency.');
cli.fv(longOpt:'ffe-verbose', args:0, 'Always print FFE decomposition.');
cli.v(longOpt:'verbose', 'Print out every atomic interaction.');
cli.pH(longOpt:'pH', args:1, argName:'7.4', 'Constant simulation pH.');
cli.t1(longOpt:'test1', 'Test 1: Lambda derivatives by finite difference.');
cli.t2(longOpt:'test2', 'Test 2: End state energies verification.');
cli.t3(longOpt:'test3', 'Test 3: Switching function and path smoothness.');
cli.i(longOpt:'iterations', args:1, argName:'n', 'Repeat Test1 several times to verify threaded replicability.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h) {
    return cli.usage();
} else if (arguments == null || arguments.size() != 1) {
    cli.usage();
    logger.info("Requires exactly one filename argument.");
    return;
} else if (!options.rl) {
    cli.usage();
    logger.info(" Specify titratable residues with --resList\n"
              + "     e.g. to titrate chainA res4 and chainB res6: -rl A4.B6");
    return;
}

// Read in command line file.
String filename = arguments.get(0);

boolean test1 = true, test2 = true, test3 = true;
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

if (options.pH) {
    pH = Double.parseDouble(options.pH);
} else {
    pH = 7.4;
}

if (options.l) {
    lambda = Double.parseDouble(options.l);
}

if (options.dx) {
    step = Double.parseDouble(options.dx);
}

if (options.bm) {
    biasMag = Double.parseDouble(options.bm);
}

if (options.st) {
    logger.info(" Launching in single-threaded mode.");
    System.setProperty("pj.nt", "1");
}

if (options.v) {
    System.setProperty("esv-verbose", "true");
}

if (options.fv) {
    System.setProperty("ffe-printOverride", "true");
}

if (arguments.size() == 0) {
    filename = "examples/lys-lys.pdb";
} else if (arguments.size() == 1) {
    logger.info("\n Testing lambda derivatives for " + filename);
} else {
    // No dual-topology mode for ESV derivatives.
    return cli.usage();
}

// ForceField
System.setProperty("forcefield", "AMOEBA_PROTEIN_2013");
// ForceField: Active
System.setProperty("esvterm", "true");
System.setProperty("lambdaterm", "true");
System.setProperty("vdwterm", "true");
System.setProperty("bondterm", "true");
System.setProperty("angleterm", "true");
System.setProperty("strbndterm", "true");
System.setProperty("ureyterm", "true");
System.setProperty("opbendterm", "true");
System.setProperty("torsionterm", "true");
System.setProperty("pitorsterm", "true");
System.setProperty("tortorterm", "true");
System.setProperty("improperterm", "true");
// ForceField: Inactive
System.setProperty("mpoleterm", "false");
System.setProperty("polarizeterm", "false");
System.setProperty("gkterm", "false");
System.setProperty("restrainterm", "false");
System.setProperty("comrestrainterm", "false");
System.setProperty("lambda_torsions", "false");

// Potential Settings
System.setProperty("polarization", "NONE");
System.setProperty("polarization-lambda-start","0.0");      // polarize on the whole range [0,1]
System.setProperty("polarization-lambda-exponent","0.0");   // polarization not softcored, only prefactored
System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization
System.setProperty("vdw-cutoff", "1000");

// ESV Settings
//System.setProperty("ffe-combineBonded", "true");
System.setProperty("esv-propagation", "false");             // don't allow ESV particle to undergo dynamics
//System.setProperty("esv-biasTerm", "true");                 // include discretization and pH biases
//System.setProperty("esv-scaleBonded", "true");             // include effects on bonded terms

Logger ffxlog = Logger.getLogger("ffx");

// Open the first topology.
open(filename);

// Parse the required lambda arguments and create TitrationESV objects.
String[] rlTokens = (options.rl).tokenize(',');

MolecularAssembly mola = (MolecularAssembly) active;
ForceFieldEnergy ffe = mola.getPotentialEnergy();
ExtendedSystem esvSystem = new ExtendedSystem(mola);

Polymer[] polymers = active.getChains();
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
        logger.severe("Couldn't find target residue " + rlTokens[i]);
    }

    MultiResidue titrating = TitrationUtils.titrationFactory(mola, target.get());
    TitrationESV esv = new TitrationESV(titrating, pH, biasMag);
    esvSystem.addVariable(esv);
}

// Attach populated esvSystem to the potential.
ffe.attachExtendedSystem(esvSystem);

// Turn off checks for overlapping atoms, which is expected for lambda=0.
ffe.getCrystal().setSpecialPositionCutoff(0.0);

// Declare arrays and get coordinates.
n = ffe.getNumberOfVariables();
assert(n % 3 == 0);
int nAtoms = n / 3;
double[] xyz = new double[n];
double[] gradient = new double[n];
double[][] esvLambdaGrad = new double[numESVs][n];
double[][][] esvLambdaGradFD = new double[2][numESVs][n];
ffe.getCoordinates(xyz);

double width = 2.0 * step;      // default step = 0.0001
double tolerance = 1.0e-3;      // pass-fail margin

/*******************************************************************************
 * Finite Difference Tests of the VdW Lambda Derivative
Check that the analytic lambda derivatives reported by van der Waals agree with
the central finite difference.
*******************************************************************************/
if (test1) {
    for (int iter = 0; iter < 10; iter++) {
        for (int i = 0; i < numESVs; i++) {
            ffe.setPrintOverride(false);
            esvSystem.setLambda(i, 0.0);
            double e0 = ffe.energyAndGradient(xyz,gradient);
            esvSystem.setLambda(i, 1.0);
            double e1 = ffe.energyAndGradient(xyz,gradient);
            logger.info(String.format(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f\n", i, e1, e0, e1-e0));
            if (iter == 0) {
                ffe.setPrintOverride(true);
            }
            
            esvSystem.setLambda(i, lambda - step);
            double eLower = ffe.energyAndGradient(xyz,esvLambdaGradFD[1][i]);

            esvSystem.setLambda(i, lambda + step);
            double eUpper = ffe.energyAndGradient(xyz,esvLambdaGradFD[0][i]);

            esvSystem.setLambda(i, lambda);
            double center = ffe.energyAndGradient(xyz,esvLambdaGrad[i]);

            double analytic = esvSystem.getdEdL(i, null, false, true);
            double numeric = (eUpper - eLower) / (2 * step);
            double error = Math.abs(numeric - analytic);

            StringBuilder sb = new StringBuilder();
            sb.append(format("\n  ESV_%d Derivative Test \n", i));
            sb.append(format(  " *********************** \n"));
            sb.append(String.format(" Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", eUpper, eLower, width, numeric));
            sb.append(String.format(" Analytic Derivative @ lambda %4.2f  >  dEdL: %+9.6g\n", lambda, analytic));
            String passFail = (error < tolerance) ? "passed" : "failed";
            sb.append(String.format(" dE/dL %6s:   %+10.6f\n", passFail, error));
            sb.append(String.format(" Numeric:        %+10.6f\n", numeric));
            sb.append(String.format(" Analytic:       %+10.6f\n", analytic));
            logger.info(sb.toString());
        }
    }
}

/*******************************************************************************
 * Dilysine End-State Verification
Verify that a lys-lys system with two ESVs can exactly reproduce the VdW
energy yielded by vanilla energy() calls on mutated PDB files.
*******************************************************************************/
if (test2) {
    esvSystem.setEsvBiasTerm(false);
    ffe.setPrintOverride(true);
    StringBuilder sb = new StringBuilder();
    sb.append(format("\n  Two-site ESV Analysis: \n"));
    sb.append(format(" ************************ \n"));

    esvSystem.setLambda(0, 0.0);
    esvSystem.setLambda(1, 0.0);
    //    System.setProperty("vdw-printInteractions", "inter-ddesv");
    ffe.energy(false, false);
    double esvVdw00 = ffe.getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "0-0", esvVdw00));
    esvSystem.setLambda(0, 1.0);
    esvSystem.setLambda(1, 0.0);
    //    System.setProperty("vdw-printInteractions", "inter-sdesv");
    ffe.energy(false, false);
    double esvVdw10 = ffe.getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "1-0", esvVdw10));
    esvSystem.setLambda(0, 0.0);
    esvSystem.setLambda(1, 1.0);
    //    System.setProperty("vdw-printInteractions", "inter-dsesv");
    ffe.energy(false, false);
    double esvVdw01 = ffe.getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "0-1", esvVdw01));
    esvSystem.setLambda(0, 1.0);
    esvSystem.setLambda(1, 1.0);
    //    System.setProperty("vdw-printInteractions", "inter-ssesv");
    ffe.energy(false, false);
    double esvVdw11 = ffe.getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "1-1", esvVdw11));
    logger.info(sb.toString());

    sb = new StringBuilder();
    sb.append(format(" Vanilla End-States: \n"));
    ffxlog.setLevel(Level.OFF);
    //    System.setProperty("vdw-printInteractions", "inter-dd");
    open("lyd-lyd.pdb");
    double lydlyd = energy().getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "lyd-lyd", lydlyd));
    //    System.setProperty("vdw-printInteractions", "inter-sd");
    open("lys-lyd.pdb");
    double lyslyd = energy().getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "lys-lyd", lyslyd));
    //    System.setProperty("vdw-printInteractions", "inter-ds");
    open("lyd-lys.pdb");
    double lydlys = energy().getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "lyd-lys", lydlys));
    //    System.setProperty("vdw-printInteractions", "inter-ss");
    open("lys-lys.pdb");
    double lyslys = energy().getVanDerWaalsEnergy();
    sb.append(format("   vdw %-7s %10.6f\n", "lys-lys", lyslys));
    ffxlog.setLevel(Level.INFO);
    logger.info(sb.toString());
}

/*******************************************************************************
 * Switching and Smoothness
Numerically ensure that the VdW energy and lambda derivatives are smooth all 
along both ESV coordinates in the dilysine system.
*******************************************************************************/
if (test3) {
    logger.info(format("  Smoothness Verification: "));
    logger.info(format(" ************************** "));
    ffxlog.setLevel(Level.OFF);
    StringBuilder energyTable = new StringBuilder(format(" %6s", "Energy"));
    StringBuilder dedlaTable = new StringBuilder(format("\n %6s", "dEdEv1"));
    StringBuilder dedlbTable = new StringBuilder(format("\n %6s", "dEdEv2"));
    energyTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %6s    %6s    %6s    %8.1f  %8.1f  %8.1f  %8.1f", 
            0.0, 0.1, 0.2, 0.3, "<", " Ev1", ">", 0.7, 0.8, 0.9, 1.0));
    dedlaTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
            0.0, 0.1, 0.2, 0.3, "<", " Ev1", ">", 0.7, 0.8, 0.9, 1.0));
    dedlbTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
            0.0, 0.1, 0.2, 0.3, "<", " Ev1", ">", 0.7, 0.8, 0.9, 1.0));
    for (double la = 0.0; la <= 1.0; la += 0.1) {
        switch (la) {
            case 0.4:
                energyTable.append(format("\n    ^ "));
                dedlaTable.append(format("\n    ^ "));
                dedlbTable.append(format("\n    ^ "));
                break;
            case 0.5:
                energyTable.append(format("\n   Ev2"));
                dedlaTable.append(format("\n   Ev2"));
                dedlbTable.append(format("\n   Ev2"));
                break;
            case 0.6:
                energyTable.append(format("\n    v "));
                dedlaTable.append(format("\n    v "));
                dedlbTable.append(format("\n    v "));
                break;
            default:
                energyTable.append(format("\n  %4.1f", la));
                dedlaTable.append(format("\n  %4.1f", la));
                dedlbTable.append(format("\n  %4.1f", la));
        }
        for (double lb = 0.0; lb <= 1.0; lb += 0.1) {
            esvSystem.setLambda(0, la);
            esvSystem.setLambda(1, lb);
            ffe.energy(true, false);
            double evdw = ffe.getVanDerWaalsEnergy();
            double dvdwdla = esvSystem.getdVdwdL(0);
            double dvdwdlb = esvSystem.getdVdwdL(1);
            energyTable.append(format("  %8.5f", evdw));
            dedlaTable.append(format("  %8.6f", dvdwdla));
            dedlbTable.append(format("  %8.6f", dvdwdlb));
        }
    }
    ffxlog.setLevel(Level.INFO);
    logger.info(energyTable.toString());
    logger.info(dedlaTable.toString());
    logger.info(dedlbTable.toString());
}
