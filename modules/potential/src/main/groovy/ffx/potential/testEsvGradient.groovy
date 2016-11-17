    
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
//import edu.rit.pj.ParallelTeam;

// First ligand atom.
int ligandStart = 1;

// Last ligand atom.
int ligandStop = -1;

// First active atom.
int activeStart = 1;

// Last active atom.
int activeStop = -1;

// First ligand atom of the 2nd topology.
int ligandStart2 = 1;

// Last ligand atom of the 2nd topology.
int ligandStop2 = -1;

// First atom for no electrostatics.
int noElecStart = 1;

// Last atom for no electrostatics.
int noElecStop = -1;

// First atom of the 2nd topology for no electrostatics.
int noElecStart2 = 1;

// Last atom of the 2nd topology for no electrostatics.
int noElecStop2 = -1;

// Lambda value.
double lambda = 0.5;

// Print out the energy for each step.
boolean print = false;
boolean singleThread = true;

// FD step size.
double step = 0.0001;

// ESV discretization bias height
double biasMag = 1.0;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testEsvGradient [options] <XYZ|PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName:'1', '(Lbd) Starting ligand atom.');
cli.f(longOpt:'final', args:1, argName:'n', '(Lbd) Final ligand atom.');
//cli.es(longOpt:'noElecStart', args:1, argName:'1', 'No Electrostatics Starting Atom.');
//cli.ef(longOpt:'noElecFinal', args:1, argName:'-1', 'No Electrostatics Final Atom.');
cli.l(longOpt:'lambda', args:1, argName:'1.0', 'Initial lambda value (for all ESVs).');
cli.rl(longOpt:'resList', required:true, args:1, 'List of residues to titrate; eg A4.A8.B2.B34');
cli.dx(longOpt:'stepSize', args:1, argName:'0.0001', 'Finite difference step size.');
cli.v(longOpt:'verbose', 'Print out the energy for each step.');
cli.bm(longOpt:'biasMag', args:1, argName:'1.0', 'ESV discretization bias height.');
cli.st(longOpt:'singleThread', args:0, 'Limit PJ concurrency.')

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

// Ligand definition.
if (options.s) {
    ligandStart = Integer.parseInt(options.s);
}
if (options.f) {
    ligandStop = Integer.parseInt(options.f);
}

// No-ES definition.
if (options.es) {
    noElecStart = Integer.parseInt(options.es);
}
if (options.ef) {
    noElecStop = Integer.parseInt(options.ef);
}

// Starting lambda value.
if (options.l) {
    lambda = Double.parseDouble(options.l);
}

// Print the energy for each step.
if (options.v) {
    print = true;
}

if (options.dx) {
    step = Double.parseDouble(options.dx);
}

if (options.bm) {
    biasMag = Double.parseDouble(options.bm);
}

if (options.st || singleThread) {
    System.setProperty("pj.nt", "1");
}

if (arguments.size() == 1) {
    logger.info("\n Testing lambda derivatives for " + filename);
} else {
    return cli.usage();
}

/*
        bondTerm = forceField.getBoolean(               ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(              ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(        ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(        ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(     ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(            ForceFieldBoolean.TORSIONTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(   ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(     ForceFieldBoolean.TORTORTERM, true);
        improperTorsionTerm = forceField.getBoolean(    ForceFieldBoolean.IMPROPERTERM, true);
        vanderWaalsTerm = forceField.getBoolean(        ForceFieldBoolean.VDWTERM, true);
        multipoleTerm = forceField.getBoolean(          ForceFieldBoolean.MPOLETERM, true);
        polarizationTerm = forceField.getBoolean(       ForceFieldBoolean.POLARIZETERM, true);
        generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
        lambdaTerm = forceField.getBoolean(             ForceField.ForceFieldBoolean.LAMBDATERM, false);
        restrainTerm = forceField.getBoolean(           ForceFieldBoolean.RESTRAINTERM, false);
        comTerm = forceField.getBoolean(                ForceFieldBoolean.COMRESTRAINTERM, false);
        lambdaTorsions = forceField.getBoolean(         ForceFieldBoolean.LAMBDA_TORSIONS, false);
        esvTerm = forceField.getBoolean(                ForceFieldBoolean.ESVTERM, false);
        printOnFailure = forceField.getBoolean(         ForceFieldBoolean.PRINT_ON_FAILURE, true);
*/

// Active terms
System.setProperty("ESVTERM", "true");
System.setProperty("LAMBDATERM", "true");
System.setProperty("VDWTERM", "true");

// Inactive terms
System.setProperty("MPOLETERM", "false");
System.setProperty("BONDTERM", "false");
System.setProperty("ANGLETERM", "false");
System.setProperty("STRBNDTERM", "false");
System.setProperty("UREYTERM", "false");
System.setProperty("OPBENDTERM", "false");
System.setProperty("TORSIONTERM", "false");
System.setProperty("PITORSTERM", "false");
System.setProperty("TORTORTERM", "false");
System.setProperty("IMPROPERTERM", "false");
System.setProperty("POLARIZETERM", "false");
System.setProperty("GKTERM", "false");
System.setProperty("RESTRAINTERM", "false");
System.setProperty("COMRESTRAINTERM", "false");
System.setProperty("LAMBDA_TORSIONS", "false");

// Potential Settings
System.setProperty("polarization", "NONE");             // !! TODO
System.setProperty("polarization-lambda-start","0.0");      // polarize on the whole range [0,1]
System.setProperty("polarization-lambda-exponent","0.0");   // polarization not softcored, only prefactored
System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization
System.setProperty("vdw-cutoff", "1000");

// Debug options
//System.setProperty("debug", "true");
System.setProperty("verbose", "true");
System.setProperty("esv-propagation", "false");
System.setProperty("pj.nt", "1");

// Open the first topology.
open(filename);

// Parse the required lambda arguments and create TitrationESV objects.
String[] rlTokens = (options.rl).tokenize(',');

MolecularAssembly mola = (MolecularAssembly) active;
ForceFieldEnergy ffe = mola.getPotentialEnergy();
ExtendedSystem esvSystem = new ExtendedSystem(mola);

List<ExtendedVariable> esvList = new ArrayList<>();
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
//    titrating.finalize(false, mola.getForceField());
    TitrationESV esv = new TitrationESV(7.4, titrating, biasMag);
    esvSystem.addVariable(esv);
}

ffe.attachExtendedSystem(esvSystem);

// Select ligand atoms
Atom[] atoms = active.getAtomArray();
int n = atoms.length;

// Turn off checks for overlapping atoms, which is expected for lambda=0.
ffe.getCrystal().setSpecialPositionCutoff(0.0);

// Reset the number of variables for the case of dual topology.
n = ffe.getNumberOfVariables();
double[] xyz = new double[n];
double[] gradient = new double[n];
double[][] esvLambdaGrad = new double[numESVs][n];
double[][][] esvLambdaGradFD = new double[2][numESVs][n];

// Number of independent atoms.
assert(n % 3 == 0);
int nAtoms = n / 3;

ffe.getCoordinates(xyz);
ffe.setPrintOverride(true);

// Compute the L = 0.0 energy.
esvSystem.setAllLambdas(0.0);
//ffe.reInit();
double e0 = ffe.energyAndGradient(xyz,gradient);

// Compute the L = 1.0 energy.
esvSystem.setAllLambdas(1.0);
//ffe.reInit();
double e1 = ffe.energyAndGradient(xyz,gradient);

logger.info(String.format(" E(0):      %20.8f.", e0));
logger.info(String.format(" E(1):      %20.8f.", e1));
logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1-e0));

// Finite-difference step size.
double width = 2.0 * step;    // default step = 0.0

// Error tolerence
double tolerance = 1.0e-3;
// Upper bound for typical gradient sizes (expected gradient)
double expGrad = 1000.0;

/******************************************************************************/
// Finite Difference tests of each ESV.
/******************************************************************************/
for (int i = 0; i < numESVs; i++) {
    double lambdaLower = lambda - step;
    double lambdaUpper = lambda + step;
    
    esvSystem.setAllLambdas(lambdaLower);
//    ffe.reInit();
    double energyLower = ffe.energyAndGradient(xyz,esvLambdaGradFD[1][i]);

    esvSystem.setAllLambdas(lambdaUpper);
//    ffe.reInit();
    double energyUpper = ffe.energyAndGradient(xyz,esvLambdaGradFD[0][i]);

    esvSystem.setAllLambdas(lambda);
//    ffe.reInit();
    double center = ffe.energyAndGradient(xyz,esvLambdaGrad[i]);
    double[] esvAnalytics = esvSystem.getdEdL(false, 0.0);
    double analytic = esvAnalytics[i];

    double numeric = (energyUpper - energyLower) / (lambdaUpper - lambdaLower);
    double error = Math.abs(numeric - analytic);

    StringBuilder sb = new StringBuilder();
    sb.append(String.format(" (ESV%d) Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", i, energyUpper, energyLower, width, numeric));
    sb.append(String.format(" (ESV%d) Analytic Derivative @ lambda %4.2f  >  dEdL: %+9.6g\n\n", i, lambda, analytic));
    String passFail = (error < tolerance) ? "passed" : "failed";
    sb.append(String.format(" (ESV%d) dE/dL %6s:   %+10.6f\n", i, passFail, error));
    sb.append(String.format(" (ESV%d) Numeric:        %+10.6f\n", i, numeric));
    sb.append(String.format(" (ESV%d) Analytic:       %+10.6f\n", i, analytic));
    logger.info(sb.toString());
}

/******************************************************************************/
// Skip second derivatives w'r't esvLambda.
return;
/******************************************************************************/

// Test Lambda gradient in the neighborhood of the lambda variable.
for (int j=0; j<3; j++) {
    lambda = initialLambda - 0.01 + 0.01 * j;
    if (lambda - step < 0.0) {
        continue;
    }
    if (lambda + step > 1.0) {
        continue;
    }

    logger.info(String.format(" Current lambda value:  %6.4f", lambda));

    // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
    ffe.setAllLambdas(lambda);
    double e = potential.energyAndGradient(xyz,gradient);

    // Analytic dEdL, d2E/dL2 and dE/dL/dX
    double dEdL = esvSystem.getdEdL();
    for (int i = 0; i < n; i++) {
        lambdaGrad[i] = 0.0;
    }
    lambdaGrad = esvSystem.getdEdXdL();

    // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
    // Plus step
    ffe.setAllLambdas(lambda + step);
    double lp = potential.energyAndGradient(xyz,lambdaGradFD[0]);
    
    // Minus step
    ffe.setAllLambdas(lambda - step);
    double lm = potential.energyAndGradient(xyz,lambdaGradFD[1]);
    
    double dEdLFD = (lp - lm) / width;
//    logger.info(String.format("ldhPlus, ldhMinus, width: %.6g, %.6g, %.6g", ldhPlus, ldhMinus, width));
//    logger.info(String.format("dEdLdhPlusSum, dEdLdhMinusSum: %.6g, %.6g", dEdLdhPlusSum, dEdLdhMinusSum));
//    logger.info(String.format("numLdh, dEdLdh[0], dEdLdhSum: %d, %.8f, %.8f", numLdh, dEdLdh[0], dEdLdhSum));
    
    double err = Math.abs(dEdLFD - dEdL);
    if (err < tolerance) {
        logger.info(String.format(" dE/dL passed:   %10.6f", err));
    } else {
        logger.info(String.format(" dE/dL failed: %10.6f", err));
    }
    logger.info(String.format(" Numeric:   %15.8f", dEdLFD));
    logger.info(String.format(" Analytic:  %15.8f", dEdL));
    
    boolean passed = true;
    for (int i = 0; i < nAtoms; i++) {
        int ii = i * 3;
        double dX = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
        double dXa = lambdaGrad[ii];
        double eX = dX - dXa;
            
        double dXdLdhFD = (esvLambdaGradFD[0][ii] - esvLambdaGradFD[1][ii]) / width;
        double dXdLdhAna = dEdXdL[ii];
        double eXLdh = dXdLdhFD - dXdLdhAna;
        
        ii++;
        double dYdLdhFD = (esvLambdaGradFD[0][ii] - esvLambdaGradFD[1][ii]) / width;
        double dYdLdhAna = dEdXdL[ii];
        double eYLdh = dYdLdhFD - dYdLdhAna;
        
        ii++;
        double dZdLdhFD = (esvLambdaGradFD[0][ii] - esvLambdaGradFD[1][ii]) / width;
        double dZdLdhAna = dEdXdL[ii];
        double eZLdh = dZdLdhFD - dZdLdhAna;
        
        double errorLdh = Math.sqrt(eXLdh*eXLdh + eYLdh*eYLdh + eZLdh*eZLdh);
        if (errorLdh < tolerance) {
            logger.fine(String.format(" dE/dX/dLdh for Atom %d passed: %10.6f", i + 1, errorLdh));
        } else {
            logger.info(String.format(" dE/dX/dLdh for Atom %d failed: %10.6f", i + 1, errorLdh));
            logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXdLdhAna,dYdLdhAna,dZdLdhAna));
            logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dXdLdhFD, dYdLdhFD, dZdLdhFD));
            passed = false;
        }
    }
    if (passed) {
        logger.info(String.format(" dE/dX/dLdh passed for all atoms"));
    }
    
    logger.info("");
}

/*
if (!fixedLambda) {
    lambdaInterface.setLambda(initialLambda);
}
for (ExtendedVariable esv : esvList) {
    esv.setLamedh(initialLamedh[esv.index]);
}
potential.getCoordinates(xyz);
potential.energyAndGradient(xyz,gradient);

logger.info(String.format(" Checking Cartesian coordinate gradient"));

double[] numeric = new double[3];
double avLen = 0.0;
int nFailures = 0;
double avGrad = 0.0;
for (int i = 0; i < nAtoms; i++) {
    int i3 = i*3;
    int i0 = i3 + 0;
    int i1 = i3 + 1;
    int i2 = i3 + 2;

    // Find numeric dX
    double orig = xyz[i0];
    xyz[i0] = xyz[i0] + step;
    double e = potential.energyAndGradient(xyz,esvLambdaGradFD[0]);
    xyz[i0] = orig - step;
    e -= potential.energyAndGradient(xyz,esvLambdaGradFD[1]);
    xyz[i0] = orig;
    numeric[0] = e / width;

    // Find numeric dY
    orig = xyz[i1];
    xyz[i1] = xyz[i1] + step;
    e = potential.energyAndGradient(xyz,esvLambdaGradFD[0]);
    xyz[i1] = orig - step;
    e -= potential.energyAndGradient(xyz,esvLambdaGradFD[1]);
    xyz[i1] = orig;
    numeric[1] = e / width;

    // Find numeric dZ
    orig = xyz[i2];
    xyz[i2] = xyz[i2] + step;
    e = potential.energyAndGradient(xyz,esvLambdaGradFD[0]);
    xyz[i2] = orig - step;
    e -= potential.energyAndGradient(xyz,esvLambdaGradFD[1]);
    xyz[i2] = orig;
    numeric[2] = e / width;

    double dx = gradient[i0] - numeric[0];
    double dy = gradient[i1] - numeric[1];
    double dz = gradient[i2] - numeric[2];
    double len = dx * dx + dy * dy + dz * dz;
    avLen += len;
    len = Math.sqrt(len);

    double grad2 = gradient[i0] * gradient[i0] + gradient[i1] * gradient[i1] + gradient[i2] * gradient[i2];
    avGrad += grad2;

    if (len > errTol) {
        logger.info(String.format(" Atom %d failed: %10.6f.",i+1,len)
            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
        ++nFailures;
        //return;
    } else {
        logger.info(String.format(" Atom %d passed: %10.6f.",i+1,len)
            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]));
    }

    if (grad2 > expGrad) {
        logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i+1, grad2));
    }
    logger.info("\n");
}

avLen = avLen / nAtoms;
avLen = Math.sqrt(avLen);
if (avLen > errTol) {
    logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol));
} else {
    logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol));
}
logger.info(String.format(" Number of atoms failing gradient test: %d", nFailures));

avGrad = avGrad / nAtoms;
avGrad = Math.sqrt(avGrad);
if (avGrad > expGrad) {
    logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad));
} else {
    logger.info(String.format(" RMS gradient: %10.6f", avGrad));
}

*/