    
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

// TEST LAMBDA GRADIENT

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
import ffx.potential.extended.TitrationESV.TitrationUtils;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
//import edu.rit.pj.ParallelTeam;

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

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.s) {
    return cli.usage();
} else if (arguments == null || arguments.size() != 1) {
    logger.info("Requires exactly one filename argument.");
    return cli.usage();
} else if (!options.rl) {
    logger.info(" Specify titratable residues with --resList\n"
              + "     e.g. to titrate chainA res4 and chainB res6: -rl A4.B6");
    return cli.usage();
} else if (options.lbd && options.fl) {
    logger.info("Contradictory: --lambda and --fixLambda.");
    return cli.usage();
}

// Read in command line file.
String filename = arguments.get(0);

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

if (arguments.size() == 1) {
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
System.setProperty("ffe-combineBonded", "true");
System.setProperty("esv-propagation", "false");             // don't allow ESV particle to undergo dynamics
//System.setProperty("esv-biasTerm", "true");                 // include discretization and pH biases
//System.setProperty("esv-scaleBonded", "false");             // include effects on bonded terms

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
Finite Difference tests of each ESV.
First derivatives only; dEdLdX and dEdL2 not yet implemented in PME.
*******************************************************************************/
for (int i = 0; i < numESVs; i++) {
    esvSystem.setLambda(i, 0.0);
    double e0 = ffe.energyAndGradient(xyz,gradient);
    esvSystem.setLambda(i, 1.0);
    double e1 = ffe.energyAndGradient(xyz,gradient);
    logger.info(String.format(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f\n", i, e1, e0, e1-e0));
    
    esvSystem.setLambda(i, lambda - step);
    double eLower = ffe.energyAndGradient(xyz,esvLambdaGradFD[1][i]);

    esvSystem.setLambda(i, lambda + step);
    double eUpper = ffe.energyAndGradient(xyz,esvLambdaGradFD[0][i]);

    esvSystem.setLambda(i, lambda);
    double center = ffe.energyAndGradient(xyz,esvLambdaGrad[i]);
    
    double analytic = esvSystem.getdEdL(i, null, false, false);
    double numeric = (eUpper - eLower) / (2 * step);
    double error = Math.abs(numeric - analytic);

    StringBuilder sb = new StringBuilder();
    sb.append(String.format(" ESV%d> Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", i, eUpper, eLower, width, numeric));
    sb.append(String.format(" ESV%d> Analytic Derivative @ lambda %4.2f  >  dEdL: %+9.6g\n", i, lambda, analytic));
    String passFail = (error < tolerance) ? "passed" : "failed";
    sb.append(String.format(" ESV%d> dE/dL %6s:   %+10.6f\n", i, passFail, error));
    sb.append(String.format(" ESV%d> Numeric:        %+10.6f\n", i, numeric));
    sb.append(String.format(" ESV%d> Analytic:       %+10.6f\n", i, analytic));
    logger.info(sb.toString());
}

/*******************************************************************************
Dilysine end-state analysis.
*******************************************************************************/
if (esvSystem.count() != 2
        || !files.contains("lys-lys.pdb") || !files.contains("lyd-lyd.pdb")
        || !files.contains("lys-lyd.pdb") || !files.contains("lyd-lys.pdb")) {
    logger.info("Run me from examples directory.");
} else {
//    String fn = "ffx/potential/structures/dilysine-lydlyd.pdb";
//    List<String> files = Arrays.asList(File.list());

    esvSystem.setEsvBiasTerm(false);
    ffe.setPrintOverride(true);
    StringBuilder sb = new StringBuilder();
    sb.append(format(" Two-site ESV Analysis: \n"));

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
    Logger.getLogger("ffx").setLevel(Level.OFF);
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
    Logger.getLogger("ffx").setLevel(Level.INFO);
    logger.info(sb.toString());
}
