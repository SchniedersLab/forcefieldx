
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
double cutoff = 1000.0;

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

def options = cli.parse(args);
if (options.h) {
    return cli.usage();
}

String filename = "lys-lys.pdb";
String[] rlTokens = new String[2];
rlTokens[0] = "A2";
rlTokens[1] = "A4";

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

boolean testVdw = true;
boolean testPermReal = true;
boolean testReciprocal = true;
boolean testPolarization = false;

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
Logger ffxlog = Logger.getLogger("ffx");
System.setProperty("ffe-combineBonded", "true");        // fold all bonded contributions into a single line
System.setProperty("ffe-decomposePme", "true");         // print the six components of PME: {real,recip}*{self,perm,ind}

// Open fully protonated and create TitrationESV objects.
MolecularAssembly mola = utils.open(filename);
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
int nAtoms = n / 3;
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
     * Analytic Derivative vs Finite Difference: VdW
    *******************************************************************************/
   /*
    if (test1 && testVdw) {
        for (int iter = 0; iter < testOneIterations; iter++) {
            for (int i = 0; i < numESVs; i++) {
                esvSystem.setLambda(i, 0.0);
                double e0 = ffe.energyAndGradient(xyz,gradient);
                esvSystem.setLambda(i, 1.0);
                double e1 = ffe.energyAndGradient(xyz,gradient);
                StringBuilder sb = new StringBuilder();
                sb.append(format(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f\n", i, e1, e0, e1-e0));

                esvSystem.setLambda(i, lambda - step);
                ffe.energy(true, true);
                double vdwLow = ffe.getVanDerWaalsEnergy();

                esvSystem.setLambda(i, lambda + step);
                ffe.energy(true, true);
                double vdwHigh = ffe.getVanDerWaalsEnergy();

                esvSystem.setLambda(i, lambda);
                ffe.energy(true, true);
                double vdwCenter = ffe.getVanDerWaalsEnergy();
                
                double vdwAnalytic = esvSystem.getdVdwdL(i);
                double vdwNumeric = (vdwHigh - vdwLow) / (2 * step);
                double vdwError = Math.abs(vdwNumeric - vdwAnalytic);

                sb.append(format("\n  ESV_%d Derivative: VdW \n", i));
                sb.append(format(" *********************** \n"));
                sb.append(format(" Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", vdwHigh, vdwLow, width, vdwNumeric));
                sb.append(format(" Analytic Derivative @ lambda %4.2f  >  dEdL: %+9.6g %+9.6g %+9.6g %+9.6g\n",
                        lambda, vdwAnalytic, realAnalytic, recipAnalytic, esvSystem.getdEdL(i)));
                String vdwPassFail = (error < tolerance) ? "passed" : "failed";
                sb.append(format(" dE/dL %6s:   %+10.6f\n", vdwPassFail, vdwError));
                sb.append(format(" Numeric:        %+10.6f\n", vdwNumeric));
                sb.append(format(" Analytic:       %+10.6f\n", vdwAnalytic));
                ffxlog.info(sb.toString());
            }
        }
    }   */

    /*******************************************************************************
     * Analytic Derivative vs Finite Difference: PME
    *******************************************************************************/
    if (test1) {
        for (int iter = 0; iter < testOneIterations; iter++) {
            for (int i = 0; i < numESVs; i++) {
                esvSystem.setLambda(i, 0.0);
                double e0 = ffe.energyAndGradient(xyz,gradient);
                esvSystem.setLambda(i, 1.0);
                double e1 = ffe.energyAndGradient(xyz,gradient);
                StringBuilder sb = new StringBuilder();
                sb.append(format(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f\n", i, e1, e0, e1-e0));

                esvSystem.setLambda(i, lambda - step);
                ffe.energy(true, true);
                double vdwLow = ffe.getVanDerWaalsEnergy();
                double realLow = ffe.getPermanentRealSpaceEnergy();
                double recipLow = ffe.getPermanentReciprocalEnergy();

                esvSystem.setLambda(i, lambda + step);
                ffe.energy(true, true);
                double vdwHigh = ffe.getVanDerWaalsEnergy();
                double realHigh = ffe.getPermanentRealSpaceEnergy();
                double recipHigh = ffe.getPermanentReciprocalEnergy();

                esvSystem.setLambda(i, lambda);
                ffe.energy(true, true);
                double vdwCenter = ffe.getVanDerWaalsEnergy();
                double realCenter = ffe.getPermanentRealSpaceEnergy();
                double recipCenter = ffe.getPermanentReciprocalEnergy();

                double vdwAnalytic = esvSystem.getdVdwdL(i);
                double vdwNumeric = (vdwHigh - vdwLow) / (2 * step);
                double vdwError = Math.abs(vdwNumeric - vdwAnalytic);
                double realAnalytic = esvSystem.getdPermRealdL(i);
                double realNumeric = (realHigh - realLow) / (2 * step);
                double realError = Math.abs(realNumeric - realAnalytic);
                double recipAnalytic = esvSystem.getdReciprocaldL(i);
                double recipNumeric = (recipHigh - recipLow) / (2 * step);
                double recipError = Math.abs(recipNumeric - recipAnalytic);

                sb.append(format("\n Finite Difference Tests: %s", esvSystem.getEsv(i).toString()));
                sb.append(format(" > Van der Waals \n"));
                sb.append(format(" Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", vdwHigh, vdwLow, width, vdwNumeric));
                sb.append(format(" Analytic Derivatives @ lambda %4.2f  >  dEdL: %+9.6g %+9.6g %+9.6g %+9.6g\n",
                        lambda, vdwAnalytic, realAnalytic, recipAnalytic, esvSystem.getdEdL(i)));
                String vdwPassFail = (vdwError < tolerance) ? "passed" : "failed";
                sb.append(format(" dE/dL %6s:   %+10.6f\n", vdwPassFail, vdwError));
                sb.append(format(" Numeric:        %+10.6f\n", vdwNumeric));
                sb.append(format(" Analytic:       %+10.6f\n", vdwAnalytic));
                
                sb.append(format(" > Permanent Real Space \n"));
                sb.append(format(" Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", realHigh, realLow, width, realNumeric));
                sb.append(format(" Analytic Derivatives @ lambda %4.2f  >  dEdL: %+9.6g %+9.6g %+9.6g %+9.6g\n",
                        lambda, vdwAnalytic, realAnalytic, recipAnalytic, esvSystem.getdEdL(i)));
                sb.append(format(" Analytic Vdw,Perm,Total:  %g %g %g %g\n",
                        esvSystem.getdVdwdL(i), esvSystem.getdPermRealdL(i), esvSystem.getdReciprocaldL(i), esvSystem.getdEdL(i)));
                String realPassFail = (realError < tolerance) ? "passed" : "failed";
                sb.append(format(" dE/dL %6s:   %+10.6f\n", realPassFail, realError));
                sb.append(format(" Numeric:        %+10.6f\n", realNumeric));
                sb.append(format(" Analytic:       %+10.6f\n", realAnalytic));
                
                sb.append(format(" > Permanent Reciprocal \n"));
                sb.append(format(" Numeric FD @ upper,lower,width,dEdL: %+9.6g %+9.6g %4.2g  >  %+9.6g\n", realHigh, realLow, width, recipNumeric));
                sb.append(format(" Analytic Derivatives @ lambda %4.2f  >  dEdL: %+9.6g %+9.6g %+9.6g %+9.6g\n",
                        lambda, vdwAnalytic, realAnalytic, recipAnalytic, esvSystem.getdEdL(i)));
                sb.append(format(" Analytic Vdw,Perm,Total:  %g %g %g %g\n",
                        esvSystem.getdVdwdL(i), esvSystem.getdPermRealdL(i), esvSystem.getdReciprocaldL(i), esvSystem.getdEdL(i)));
                String recipPassFail = (recipError < tolerance) ? "passed" : "failed";
                sb.append(format(" dE/dL %6s:   %+10.6f\n", recipPassFail, recipError));
                sb.append(format(" Numeric:        %+10.6f\n", recipNumeric));
                sb.append(format(" Analytic:       %+10.6f\n", recipAnalytic));
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
        esvSystem.setEsvBiasTerm(false);

        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 0.0);
        ffe.energy(true, true);
        double esvVdw00 = ffe.getVanDerWaalsEnergy();
        double esvPme00 = ffe.getPermanentRealSpaceEnergy();

        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 1.0);
        ffe.energy(true, true);
        double esvVdw01 = ffe.getVanDerWaalsEnergy();
        double esvPme01 = ffe.getPermanentRealSpaceEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 0.0);
        ffe.energy(true, true);
        double esvVdw10 = ffe.getVanDerWaalsEnergy();
        double esvPme10 = ffe.getPermanentRealSpaceEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 1.0);
        ffe.energy(true, true);
        double esvVdw11 = ffe.getVanDerWaalsEnergy();
        double esvPme11 = ffe.getPermanentRealSpaceEnergy();

        StringBuilder sb = new StringBuilder();
        sb.append(format("\n  Two-site ESV Analysis: \n"));
        sb.append(format(" ************************ \n"));
        sb.append(format(" Extended System (1 Assembly) \n"));
        if (testVdw) {
            sb.append(format("   vanWaals %-7s %10.6f\n", "0-0", esvVdw00));
            sb.append(format("   vanWaals %-7s %10.6f\n", "0-1", esvVdw01));
            sb.append(format("   vanWaals %-7s %10.6f\n", "1-0", esvVdw10));
            sb.append(format("   vanWaals %-7s %10.6f\n", "1-1", esvVdw11));
        }
        if (testPermReal) {
            sb.append(format("   permReal %-7s %10.6f\n", "0-0", esvPme00));
            sb.append(format("   permReal %-7s %10.6f\n", "0-1", esvPme01));
            sb.append(format("   permReal %-7s %10.6f\n", "1-0", esvPme10));
            sb.append(format("   permReal %-7s %10.6f\n", "1-1", esvPme11));
        }

        ForceFieldEnergy vanillaFFE;
        sb.append(format(" Mutated Files (4 Assemblies) \n"));
        StringBuilder vdw = new StringBuilder();
        StringBuilder pme = new StringBuilder();
        List<String> endStates = Arrays.asList("lyd-lyd", "lyd-lys", "lys-lyd", "lys-lys");
        for (String state : endStates) {
            MolecularAssembly vanilla = utils.openQuietly(state + ".pdb");
            vanillaFFE = vanilla.getPotentialEnergy();
            vanillaFFE.energy(true, true);
            if (testVdw) {
                double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
                vdw.append(format("   vanWaals %-7s %10.6f\n", state, vanWaals));
            }
            if (testPermReal) {
                double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
                pme.append(format("   permReal %-7s %10.6f\n", state, permReal));
            }
            utils.close(vanilla);
        }
        if (testVdw) sb.append(vdw.toString());
        if (testPermReal) sb.append(pme.toString());
        ffxlog.info(sb.toString());
        
        esvSystem.setEsvBiasTerm(true);
    }

    /*******************************************************************************
     * Switching and Smoothness
    Numerically ensure that the VdW energy and lambda derivatives are smooth all 
    along both ESV coordinates in the dilysine system.
    *******************************************************************************/
    if (test3) {
        StringBuilder vdwTable = new StringBuilder(format(" %6s", "Energy"));
        StringBuilder dvdwdaTable = new StringBuilder(format("\n %6s", "dEdEvA"));
        StringBuilder dvdwdbTable = new StringBuilder(format("\n %6s", "dEdEvB"));
        StringBuilder pmeTable = new StringBuilder(format(" %6s", "Energy"));
        StringBuilder dpmedaTable = new StringBuilder(format("\n %6s", "dEdEvA"));
        StringBuilder dpmedbTable = new StringBuilder(format("\n %6s", "dEdEvB"));
        if (testVdw) {
            vdwTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %6s    %6s    %6s    %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
            dvdwdaTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
            dvdwdbTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
        }
        if (testPermReal) {
            pmeTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %6s    %6s    %6s    %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
            dpmedaTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
            dpmedbTable.append(format(" %8.1f  %8.1f  %8.1f  %8.1f  %8s  %8s  %8s  %8.1f  %8.1f  %8.1f  %8.1f", 
                    0.0, 0.1, 0.2, 0.3, "<", " EvA", ">", 0.7, 0.8, 0.9, 1.0));
        }
        for (double evA = 0.0; evA <= 1.0; evA += 0.1) {
            switch (evA) {
                case 0.4:
                    vdwTable.append(format("\n    ^ "));
                    dvdwdaTable.append(format("\n    ^ "));
                    dvdwdbTable.append(format("\n    ^ "));
                    pmeTable.append(format("\n    ^ "));
                    dpmedaTable.append(format("\n    ^ "));
                    dpmedbTable.append(format("\n    ^ "));
                    break;
                case 0.5:
                    vdwTable.append(format("\n   EvB"));
                    dvdwdaTable.append(format("\n   EvB"));
                    dvdwdbTable.append(format("\n   EvB"));
                    pmeTable.append(format("\n   EvB"));
                    dpmedaTable.append(format("\n   EvB"));
                    dpmedbTable.append(format("\n   EvB"));
                    break;
                case 0.6:
                    vdwTable.append(format("\n    v "));
                    dvdwdaTable.append(format("\n    v "));
                    dvdwdbTable.append(format("\n    v "));
                    pmeTable.append(format("\n    v "));
                    dpmedaTable.append(format("\n    v "));
                    dpmedbTable.append(format("\n    v "));
                    break;
                default:
                    vdwTable.append(format("\n  %4.1f", evA));
                    dvdwdaTable.append(format("\n  %4.1f", evA));
                    dvdwdbTable.append(format("\n  %4.1f", evA));
                    pmeTable.append(format("\n  %4.1f", evA));
                    dpmedaTable.append(format("\n  %4.1f", evA));
                    dpmedbTable.append(format("\n  %4.1f", evA));
            }
            for (double evB = 0.0; evB <= 1.0; evB += 0.1) {
                esvSystem.setLambda('A', evA);
                esvSystem.setLambda('B', evB);
                ffe.energy(true, false);
                if (testVdw) {
                    double evdw = ffe.getVanDerWaalsEnergy();
                    double dvdwdevA = esvSystem.getdVdwdL(0);
                    double dvdwdevB = esvSystem.getdVdwdL(1);
                    vdwTable.append(format("  %8.5f", evdw));
                    dvdwdaTable.append(format("  %8.6f", dvdwdevA));
                    dvdwdbTable.append(format("  %8.6f", dvdwdevB));
                }
                if (testPermReal) {
                    double ePermReal = ffe.getPermanentRealSpaceEnergy();
                    double dPermRealdEvA = esvSystem.getdPermRealdL(0);
                    double dPermRealdEvB = esvSystem.getdPermRealdL(1);
                    pmeTable.append(format(" %9.4f", ePermReal));
                    dpmedaTable.append(format(" %9.4f", dPermRealdEvA));
                    dpmedbTable.append(format(" %9.4f", dPermRealdEvB));
                }
            }
        }
        if (testVdw) {
            ffxlog.info(format("  Smoothness Verification: VDW "));
            ffxlog.info(format(" ****************************** "));
            ffxlog.info(vdwTable.toString());
            ffxlog.info(dvdwdaTable.toString());
            ffxlog.info(dvdwdbTable.toString());
            ffxlog.info("");
        }
        if (testPermReal) {
            ffxlog.info(format("\n  Smoothness Verification: PermReal "));
            ffxlog.info(format(" *********************************** "));
            ffxlog.info(pmeTable.toString());
            ffxlog.info(dpmedaTable.toString());
            ffxlog.info(dpmedbTable.toString());
            ffxlog.info("");
        }    
    }
} catch (RuntimeException ex) {
    /* Dump system info on any unhandled exception. */
    esvSystem.crashDump();
    throw ex;
}
