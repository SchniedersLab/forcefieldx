    
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
double lambda = 1.0;

// Print out the energy for each step.
boolean print = false;

// FD step size.
double lbdStep = 0.0;
double ldhStep = 0.0001;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc testEsvGradient [options] <XYZ|PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.s(longOpt:'start', args:1, argName:'1', '(Lbd) Starting ligand atom.');
cli.f(longOpt:'final', args:1, argName:'n', '(Lbd) Final ligand atom.');
cli.es(longOpt:'noElecStart', args:1, argName:'1', '(Lbd) No Electrostatics Starting Atom.');
cli.ef(longOpt:'noElecFinal', args:1, argName:'-1', '(Lbd) No Electrostatics Final Atom.');
cli.lbd(longOpt:'lambda', args:1, argName:'1.0', 'Lambda (Lbd) value; multiplicative* with ldh.');
cli.ldh(longOpt:'lamedh', args:1, 'Array of lamedh (Ldh) values to test.');
cli.rl(longOpt:'resList', required:true, args:1, '(Ldh) Titrate a list of residues (eg A4.A8.B2.B34)');
cli.dx(longOpt:'stepSize', args:1, argName:'0.0001', 'Finite difference step size.');
cli.v(longOpt:'verbose', 'Print out the energy for each step.');

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

// Fixed lambda value.
if (options.lbd) {
    lambda = Double.parseDouble(options.lbd);
}

// Print the energy for each step.
if (options.v) {
    print = true;
}

if (options.dx) {
    ldhStep = Double.parseDouble(options.dx);
}

if (arguments.size() == 1) {
    logger.info("\n Testing lamedh derivatives for " + filename);
} else {
    return cli.usage();
}

// Stuff that's OFF
System.setProperty("strbndterm", "false");
System.setProperty("opbendterm", "false");
System.setProperty("torsionterm", "false");
System.setProperty("tortorterm", "false");
System.setProperty("pitorsterm", "false");
System.setProperty("mpoleterm", "false");               // !! TODO

// Polarization keys
System.setProperty("polarization", "NONE");             // !! TODO
System.setProperty("polarization-lambda-start","0.0");      // polarize on the whole range [0,1]
System.setProperty("polarization-lambda-exponent","0.0");   // polarization not softcored, only prefactored
System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization

// Stuff that's ON
System.setProperty("esvterm", "true");
System.setProperty("lambdaterm", "true");
System.setProperty("bondterm", "true");
System.setProperty("angleterm", "true");
System.setProperty("vdwterm", "true");

// Test parameters
System.setProperty("vdw-cutoff", "1000");

// Open the first topology.
open(filename);
MolecularAssembly mola = active;

// Parse the required lamedh arguments and create TitrationESV objects.
String[] rlTokens = (options.rl).tokenize(',');
final int numESVs = rlTokens.length;
String[] ldhTokens;
if (options.ldh) {
    ldhTokens = (options.ldh).tokenize(',');
    if (ldhTokens.length != numLdh) {
        logger.warning("Number of --lamedh inputs must match --resList.");
    }
} else {
    ldhTokens = new String[numESVs];
    for (int i = 0; i < numESVs; i++) {
        ldhTokens[i] = 0.5;
    }
}

for (int i = 0; i < numESVs; i++) {
    logger.info(" (Groovy) Ldh: " + rlTokens[i] + ", " + ldhTokens[i]);
}

List<ExtendedVariable> esvList = new ArrayList<>();
Polymer[] polymers = active.getChains();
double[] lamedh = new double[numESVs];
double temperature = 298.15;
double dt = 1.0;
for (int i = 0; i < numESVs; i++) {
    if (ldhTokens != null) {
        lamedh[i] = Double.parseDouble(ldhTokens[i]);
    } else {
        lamedh[i] = 0.5;
    }
    
    Character chainID = rlTokens[i].charAt(0);
    int resNum = Integer.parseInt(rlTokens[i].substring(1));
    Optional<Residue> target = new Optional<>();
    for (Polymer p : polymers) {
        if (p.getChainID().equals(chainID)) {
            target = p.getResidues().parallelStream()
                .filter {res -> res.getResidueNumber() == resNum}
                .findFirst();
            break;
        }
    }
    if (!target.isPresent()) {
        logger.severe("Couldn't find target residue " + rlTokens[i]);
    }
    
    
    
//        AminoAcid3 aa = res.getAminoAcid3();
//        if (aa == AminoAcid3.HIE) {
//            throw new UnsupportedOperationException("3-way HIS Ldh under construction.");
//        }
//        MultiResidue titrating = TitrationUtils.titrationFactory(mola, mola.getForceField(), ffe, res);
    MultiResidue titrating = TitrationUtils.titrationFactory(mola, target.get());
    TitrationESV esv = new TitrationESV(titrating, temperature, dt);
    active.getPotentialEnergy().addExtendedVariable(esv);
    esvList.add(esv);
}

// Select ligand atoms
Atom[] atoms = active.getAtomArray();
int n = atoms.length;

// Apply ligand atom selection
for (int i = ligandStart; i <= ligandStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setApplyLambda(true);
    ai.print();
}

// Apply the no electrostatics atom selection
if (noElecStart < 1) {
    noElecStart = 1;
}
if (noElecStop > atoms.length) {
    noElecStop = atoms.length;
}
for (int i = noElecStart; i <= noElecStop; i++) {
    Atom ai = atoms[i - 1];
    ai.setElectrostatics(false);
    ai.print();
}

ForceFieldEnergy ffe = active.getPotentialEnergy();
// Turn off checks for overlapping atoms, which is expected for lambda=0.
ffe.getCrystal().setSpecialPositionCutoff(0.0);

// Reset the number of variables for the case of dual topology.
n = ffe.getNumberOfVariables();
double[] x = new double[n];
double[] gradient = new double[n];
double[] lambdaGrad = new double[n];
double[][] lamedhGrad = new double[numESVs][n];
double[][] lambdaGradFD = new double[2][n];
double[][][] lamedhGradFD = new double[2][numESVs][n];

// Number of independent atoms.
assert(n % 3 == 0);
int nAtoms = n / 3;

ffe.setLambda(lambda);

// Compute the L = 0.0 energy.
for (ExtendedVariable esv : esvList) {
    esv.setLamedh(0.0);
}
ffe.getCoordinates(x);
double e0 = ffe.energyAndGradient(x,gradient);

// Compute the L = 1.0 energy.
for (ExtendedVariable esv : esvList) {
    esv.setLamedh(1.0);
}
double e1 = ffe.energyAndGradient(x,gradient);

logger.info(String.format(" E(0):      %20.8f.", e0));
logger.info(String.format(" E(1):      %20.8f.", e1));
logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1-e0));

// Finite-difference step size.
double lbdWidth = 2.0 * lbdStep;    // default lbdStep = 0.0
double ldhWidth = 2.0 * ldhStep;

// Error tolerence
double errTol = 1.0e-3;
// Upper bound for typical gradient sizes (expected gradient)
double expGrad = 1000.0;

/******************************************************************************/
// Finite Difference tests of each ESV individually.
/******************************************************************************/
for (int iESV = 0; iESV < numESVs; iESV++) {
    esvList.get(iESV).setLamedh(lamedh[iESV] - ldhStep);
    double lower = ffe.energyAndGradient(x,lamedhGradFD[1][iESV]);

    esvList.get(iESV).setLamedh(lamedh[iESV] + ldhStep);
    double upper = ffe.energyAndGradient(x,lamedhGradFD[0][iESV]);

    esvList.get(iESV).setLamedh(lamedh[iESV]);
    double center = ffe.energyAndGradient(x,lamedhGrad[iESV]);
    double[] esvAnalytics = ffe.getdEdLdh();
    double alytic1 = esvAnalytics[iESV];

    double fd1 = (upper - lower) / ldhWidth;
    double err1 = Math.abs(fd1 - alytic1);

    StringBuilder sb = new StringBuilder();
    sb.append(String.format(" (ESV%d) Numeric FD @ LdhPlus,LdhMinus,width,dEdLdh: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", iESV, upper, lower, ldhWidth, fd1));
    sb.append(String.format(" (ESV%d) Analytic Derivative @ Ldh %4.2f  >  dEdLdh: %+9.6g\n\n", iESV, lamedh[iESV], alytic1));
    String passFail1 = (err1 < errTol) ? "passed" : "failed";
    sb.append(String.format(" (ESV%d) dE/dL %6s:   %+10.6f\n", iESV, passFail1, err1));
    sb.append(String.format(" (ESV%d) Numeric:        %+10.6f\n", iESV, fd1));
    sb.append(String.format(" (ESV%d) Analytic:       %+10.6f\n", iESV, alytic1));
    logger.info(sb.toString());
}

/******************************************************************************/
// Skip the loop?
return;
/******************************************************************************/

/*

// Test Lambda gradient in the neighborhood of the lambda variable.
for (int j=0; j<3; j++) {
    if (!fixedLambda) {
        lambda = initialLambda - 0.01 + 0.01 * j;
        if (lambda - lbdStep < 0.0) {
            continue;
        }
        if (lambda + lbdStep > 1.0) {
            continue;
        }
    } else {
        lambda = 1.0;
        lbdStep = 0.0;
        lbdWidth = 0.0;
    }
    
    // Initialize lamedh arrays
    double[] lamedh = new double[numLdh];
    double[] lamedhPlus = new double[numLdh];
    double[] lamedhMinus = new double[numLdh];
    boolean oob = false;
    for (int i = 0; i < numLdh; i++) {
        lamedh[i] = initialLamedh[i] - 0.01 + 0.01 * j;
        if (lamedh[i] - lbdStep < 0.0 || lamedh[i] + lbdStep > 1.0) {
            oob = true;
            break;
        }
    }
    if (oob) {
        continue;
    }
    for (int i = 0; i < numLdh; i++) {
        lamedhPlus[i] = lamedh[i] + ldhStep;
        lamedhMinus[i] = lamedh[i] - ldhStep;
    }

    logger.info(String.format(" Current lambda value:  %6.4f", lambda));
    logger.info(String.format(" Current lamedh values: %s", Arrays.toString(lamedh)));
    if (!fixedLambda) {
        lambdaInterface.setLambda(lambda);
    }
    for (ExtendedVariable esv : esvList) {
        esv.setLamedh(lamedh[esv.index]);
    }

    // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
    double e = potential.energyAndGradient(x,gradient);

    // Analytic dEdL, d2E/dL2 and dE/dL/dX
    if (!fixedLambda) {
        double dEdL = lambdaInterface.getdEdL();
        double d2EdL2 = lambdaInterface.getd2EdL2();
        for (int i = 0; i < n; i++) {
            lambdaGrad[i] = 0.0;
        }
        potential.getdEdXdL(lambdaGrad);
    }

    // Analytic dEdLdh, d2E/dLdh2 and dE/dX/dLdh
    double[] dEdLdh = ffe.getdEdLdh();
    double[] d2EdLdh2 = ffe.getd2EdLdh2();
    double[][] dEdXdLdh = ffe.getdEdXdLdh();

    // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
    // Plus step
    if (!fixedLambda) {
        lambdaInterface.setLambda(lambda + lbdStep);
    }
    for (ExtendedVariable esv : esvList) {
        esv.setLamedh(lamedhPlus[esv.index]);
    }
    if (!fixedLambda) {
        double lp = potential.energyAndGradient(x,lambdaGradFD[0]);
        double dedlp = lambdaInterface.getdEdL();
    }
    double ldhPlus = potential.energyAndGradient(x,lamedhGradFD[0]);
    double[] dEdLdhPlus = forceFieldEnergy.getdEdLdh();
    double dEdLdhPlusSum = 0.0;
    for (int i = 0; i < numLdh; i++) {
        dEdLdhPlusSum += dEdLdhPlus[i];
    }
    
    // Minus step
    if (!fixedLambda) {
        lambdaInterface.setLambda(lambda - lbdStep);
    }
    for (ExtendedVariable esv : esvList) {
        esv.setLamedh(lamedhMinus[esv.index]);
    }
    if (!fixedLambda) {
        double lm = potential.energyAndGradient(x,lambdaGradFD[1]);
        double dedlm = lambdaInterface.getdEdL();
    }
    double ldhMinus = potential.energyAndGradient(x,lamedhGradFD[1]);
    double[] dEdLdhMinus = forceFieldEnergy.getdEdLdh();
    double dEdLdhMinusSum = 0.0;
    for (int i = 0; i < numLdh; i++) {
        dEdLdhMinusSum += dEdLdhMinus[i];
    }

    double dEdLdhSum = 0.0;
    double d2EdLdh2Sum = 0.0;
    double[] dEdXdLdhSum = new double[n];
    for (int iESV = 0; iESV < numLdh; iESV++) {
        dEdLdhSum += dEdLdh[iESV];
        d2EdLdh2Sum += d2EdLdh2[iESV];
        for (int iAtom = 0; iAtom < n; iAtom++) {
            dEdXdLdhSum[iAtom] += dEdXdLdh[iESV][iAtom];
        }
    }
    
    if (!fixedLambda) {
        double dEdLFD = (lp - lm) / lbdWidth;
        double d2EdL2FD = (dedlp - dedlm) / lbdWidth;
    }
    double dEdLdhFD = (ldhPlus - ldhMinus) / ldhWidth;
    double d2EdLdh2FD = (dEdLdhPlusSum - dEdLdhMinusSum) / ldhWidth;
//    logger.info(String.format("ldhPlus, ldhMinus, ldhWidth: %.6g, %.6g, %.6g", ldhPlus, ldhMinus, ldhWidth));
//    logger.info(String.format("dEdLdhPlusSum, dEdLdhMinusSum: %.6g, %.6g", dEdLdhPlusSum, dEdLdhMinusSum));
//    logger.info(String.format("numLdh, dEdLdh[0], dEdLdhSum: %d, %.8f, %.8f", numLdh, dEdLdh[0], dEdLdhSum));
    
    if (!fixedLambda) {
        double err = Math.abs(dEdLFD - dEdL);
        if (err < errTol) {
            logger.info(String.format(" dE/dL passed:   %10.6f", err));
        } else {
            logger.info(String.format(" dE/dL failed: %10.6f", err));
        }
        logger.info(String.format(" Numeric:   %15.8f", dEdLFD));
        logger.info(String.format(" Analytic:  %15.8f", dEdL));
    }
    double errLdh = Math.abs(dEdLdhFD - dEdLdhSum);
    if (errLdh < errTol) {
        logger.info(String.format(" dE/dLdh passed:   %10.6f", errLdh));
    } else {
        logger.info(String.format(" dE/dLdh failed: %10.6f", errLdh));
    }
    logger.info(String.format(" Numeric:   %15.8f", dEdLdhFD));
    logger.info(String.format(" Analytic:  %15.8f", dEdLdhSum));

    if (!fixedLambda) {
        err = Math.abs(d2EdL2FD - d2EdL2);
        if (err < errTol) {
            logger.info(String.format(" d2E/dL2 passed: %10.6f", err));
        } else {
            logger.info(String.format(" d2E/dL2 failed: %10.6f", err));
        }
        logger.info(String.format(" Numeric:   %15.8f", d2EdL2FD));
        logger.info(String.format(" Analytic:  %15.8f", d2EdL2));
    }
    errLdh = Math.abs(d2EdLdh2FD - d2EdLdh2Sum);
    if (errLdh < errTol) {
        logger.info(String.format(" d2E/dLdh2 passed: %10.6f", errLdh));
    } else {
        logger.info(String.format(" d2E/dLdh2 failed: %10.6f", errLdh));
    }
    logger.info(String.format(" Numeric:   %15.8f", d2EdLdh2FD));
    logger.info(String.format(" Analytic:  %15.8f", d2EdLdh2Sum));
    
    boolean passed = true;

    for (int i = 0; i < nAtoms; i++) {   // TODO should be nAtoms
        int ii = i * 3;
        if (!fixedLambda) {
            double dX = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / lbdWidth;
            double dXa = lambdaGrad[ii];
            double eX = dX - dXa;
        }
        double dXdLdhFD = (lamedhGradFD[0][ii] - lamedhGradFD[1][ii]) / ldhWidth;
        double dXdLdhAna = dEdXdLdhSum[ii];
        double eXLdh = dXdLdhFD - dXdLdhAna;
        ii++;
        if (!fixedLambda) {
            double dY = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / lbdWidth;
            double dYa = lambdaGrad[ii];
            double eY = dY - dYa;
        }
        double dYdLdhFD = (lamedhGradFD[0][ii] - lamedhGradFD[1][ii]) / ldhWidth;
        double dYdLdhAna = dEdXdLdhSum[ii];
        double eYLdh = dYdLdhFD - dYdLdhAna;
        ii++;
        if (!fixedLambda) {
            double dZ = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / lbdWidth;
            double dZa = lambdaGrad[ii];
            double eZ = dZ - dZa;
        }
        double dZdLdhFD = (lamedhGradFD[0][ii] - lamedhGradFD[1][ii]) / ldhWidth;
        double dZdLdhAna = dEdXdLdhSum[ii];
        double eZLdh = dZdLdhFD - dZdLdhAna;
//        logger.info(String.format("i,dXYZ,gradFD[x,y,z][0,1]: %d, (%.6g %.6g %.6g), %.6g %.6g %.6g %.6g %.6g %.6g", 
//                i, dXdLdhFD, dYdLdhFD, dZdLdhFD, 
//                lamedhGradFD[0][ii-2], lamedhGradFD[1][ii-2], 
//                lamedhGradFD[0][ii-1], lamedhGradFD[1][ii-1], 
//                lamedhGradFD[0][ii],   lamedhGradFD[1][ii]  ));
//        logger.info(String.format("i,dXYZ,gradFD[x,y,z][0,1]: %d, (%10.8f %10.8f %10.8f), (%10.8f %10.8f %10.8f)", 
//                i, dXdLdhFD, dYdLdhFD, dZdLdhFD, dXdLdhAna, dYdLdhAna, dZdLdhAna));

        if (!fixedLambda) {
            double error = Math.sqrt(eX * eX + eY * eY + eZ * eZ);
            if (error < errTol) {
                logger.fine(String.format(" dE/dX/dL for Atom %d passed: %10.6f", i + 1, error));
            } else {
                logger.info(String.format(" dE/dX/dL for Atom %d failed: %10.6f", i + 1, error));
                logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa,dYa,dZa));
                logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ));
                passed = false;
            }
        }
        
        double errorLdh = Math.sqrt(eXLdh*eXLdh + eYLdh*eYLdh + eZLdh*eZLdh);
        if (errorLdh < errTol) {
            logger.fine(String.format(" dE/dX/dLdh for Atom %d passed: %10.6f", i + 1, errorLdh));
        } else {
            logger.info(String.format(" dE/dX/dLdh for Atom %d failed: %10.6f", i + 1, errorLdh));
            logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXdLdhAna,dYdLdhAna,dZdLdhAna));
            logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dXdLdhFD, dYdLdhFD, dZdLdhFD));
            passed = false;
        }
    }
    if (passed) {
        if (!fixedLambda) {
            logger.info(String.format(" dE/dX/dL passed for all atoms"));
        }
        logger.info(String.format(" dE/dX/dLdh passed for all atoms"));
    }

    logger.info("");
}

if (!fixedLambda) {
    lambdaInterface.setLambda(initialLambda);
}
for (ExtendedVariable esv : esvList) {
    esv.setLamedh(initialLamedh[esv.index]);
}
potential.getCoordinates(x);
potential.energyAndGradient(x,gradient);

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
    double orig = x[i0];
    x[i0] = x[i0] + ldhStep;
    double e = potential.energyAndGradient(x,lamedhGradFD[0]);
    x[i0] = orig - ldhStep;
    e -= potential.energyAndGradient(x,lamedhGradFD[1]);
    x[i0] = orig;
    numeric[0] = e / ldhWidth;

    // Find numeric dY
    orig = x[i1];
    x[i1] = x[i1] + ldhStep;
    e = potential.energyAndGradient(x,lamedhGradFD[0]);
    x[i1] = orig - ldhStep;
    e -= potential.energyAndGradient(x,lamedhGradFD[1]);
    x[i1] = orig;
    numeric[1] = e / ldhWidth;

    // Find numeric dZ
    orig = x[i2];
    x[i2] = x[i2] + ldhStep;
    e = potential.energyAndGradient(x,lamedhGradFD[0]);
    x[i2] = orig - ldhStep;
    e -= potential.energyAndGradient(x,lamedhGradFD[1]);
    x[i2] = orig;
    numeric[2] = e / ldhWidth;

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