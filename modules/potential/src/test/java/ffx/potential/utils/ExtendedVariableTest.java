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
package ffx.potential.utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.lang.StringUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.SBLogger;

import static ffx.potential.extended.SBLogger.SB;

/**
 *
 * @author slucore
 */
@RunWith(Parameterized.class)
public class ExtendedVariableTest {
    /**
     * Handle on the root log allows easy global silencing (like during multi-file setup).
     */
    private static final Logger logMaster = Logger.getLogger("ffx");
    private static final PotentialsUtils utils = new PotentialsUtils();
    
    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            /* Arguments to the class constructor. */
            {"A2", "A4", 1.0e-4, 
             "lyd-lyd-cryst.pdb", 65.339655,
             "lyd-lys-cryst.pdb", 129.231033,
             "lys-lyd-cryst.pdb", 127.384741,
             "lys-lys-cryst.pdb", 197.920812}
        });
    }
    
    private static final boolean testVdw = true;
    private static final boolean testPermReal = true;
    private static final boolean testReciprocal = true;
    private static final boolean testPolarization = false;
    private static final double cutoff = 5.0;
    private static final boolean verbose = false;
    private static final String resourcePrefix = "ffx/potential/structures/";
    
    private final String residA, residB;
    private final String filename00;
    private final double energy00;
    private final String filename01;
    private final double energy01;
    private final String filename10;
    private final double energy10;
    private final String filename11;
    private final double energy11;
    private final double tolerance;
    
    private static final double[] initialLambda = new double[]{0.28, 0.55};
    private static final double pH = 7.4;
    private static final int numESVs = 2;
    
    private MolecularAssembly mola;
    private ExtendedSystem esvSystem;
    private ForceFieldEnergy ffe;
    
    public ExtendedVariableTest(
            String residA, String residB, double tolerance,
            String filename1, double energy1,
            String filename2, double energy2,
            String filename3, double energy3,
            String filename4, double energy4) {
        this.tolerance = tolerance;
        this.residA = residA;
        this.residB = residB;
        this.filename00 = filename1;
        this.filename01 = filename2;
        this.filename10 = filename3;
        this.filename11 = filename4;
        this.energy00 = energy1;
        this.energy01 = energy2;
        this.energy10 = energy3;
        this.energy11 = energy4;
    }
    
    /**
     * Set system properties prior to class construction.
     */
    @BeforeClass
    public static void setUpClass() {
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
        System.setProperty("esv-vdw", String.valueOf(testVdw));
        System.setProperty("mpoleterm", String.valueOf(testPermReal));          // permanent real space
        System.setProperty("pme-qi", String.valueOf(testPermReal));
        System.setProperty("esv-pme", String.valueOf(testPermReal));
        System.setProperty("recipterm", String.valueOf(testReciprocal));        // permanent reciprocal space
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
    }
    
    /**
     * Unset all modified system properties.
     */
    @AfterClass
    public static void tearDownClass() {
        // Active Potential
        System.clearProperty("forcefield");
        System.clearProperty("esvterm");
        System.clearProperty("lambdaterm");
        System.clearProperty("bondterm");
        System.clearProperty("angleterm");
        System.clearProperty("strbndterm");
        System.clearProperty("ureyterm");
        System.clearProperty("opbendterm");
        System.clearProperty("torsionterm");
        System.clearProperty("pitorsterm");
        System.clearProperty("tortorterm");
        System.clearProperty("improperterm");

        // Optional Potential
        System.clearProperty("vdwterm");
        System.clearProperty("esv-useVdw");
        System.clearProperty("mpoleterm");
        System.clearProperty("pme-qi");
        System.clearProperty("esv-usePme");
        System.clearProperty("recipterm");
        System.clearProperty("polarizeterm");
        System.clearProperty("polarization");

        // Inactive Potential
        System.clearProperty("gkterm");
        System.clearProperty("restrainterm");
        System.clearProperty("comrestrainterm");
        System.clearProperty("lambda_torsions");

        // Potential Settings
        System.clearProperty("permanent-lambda-alpha");
        System.clearProperty("permanent-lambda-exponent");
        System.clearProperty("polarization-lambda-start");
        System.clearProperty("polarization-lambda-exponent");
        System.clearProperty("ligand-vapor-elec");
        System.clearProperty("no-ligand-condensed-scf");
        System.clearProperty("intramolecular-softcore");
        System.clearProperty("intermolecular-softcore");
        System.clearProperty("vdw-cutoff");
        System.clearProperty("ewald-cutoff");

        // ESV Settings
        System.clearProperty("esv-propagation");
        System.clearProperty("esv-biasTerm");
        System.clearProperty("esv-scaleBonded");
        System.clearProperty("esv-backgroundBonded");
        System.clearProperty("esv-scaleUnshared");

        // Logging Settings
        System.clearProperty("ffe-combineBonded");
        System.clearProperty("ffe-decomposePme");
    }
    
    @After
    public void tearDown() {
        mola = null;
        ffe = null;
        esvSystem = null;
        logMaster.setLevel(Level.INFO);
    }
    
    public void setup(boolean quietly) {
        mola = openResource(filename11, quietly);
        ffe = mola.getPotentialEnergy();
        esvSystem = new ExtendedSystem(mola, pH);
        // Create extended variables and add them to the system.
        String[] resIDs = new String[]{residA, residB};
        esvSystem.populate(resIDs);
        ffe.attachExtendedSystem(esvSystem);
        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        ffe.getCrystal().setSpecialPositionCutoff(0.0);
        for (int i = 0; i < numESVs; i++) {
            esvSystem.setLambda(i, initialLambda[i]);
        }
    }
    
    /**
     * Locates files packed into the uberjar by Maven.
     */
    private MolecularAssembly openResource(String filename, boolean quietly) {
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(resourcePrefix + filename).getPath());
        return (quietly) ? utils.openQuietly(structure) : utils.open(structure);
    }
    private MolecularAssembly openResource(String filename) {
        return openResource(filename, false);
    }
    private MolecularAssembly openResourceQuietly(String filename) {
        return openResource(filename, true);
    }
        
    /**
     * Analytic Derivative vs Finite Difference: VdW,PermReal,PermRecip,Total
     */
    @Test
    public void testDerivatives() {
        logMaster.setLevel(Level.WARNING);
        setup(false);
        final double step = 0.000001;
        final double width = 2.0 * step;
        for (int i = 0; i < numESVs; i++) {
            esvSystem.setLambda(i, initialLambda[i]);
        }
        for (int i = 0; i < numESVs; i++) {
            esvSystem.setLambda(i, 0.0);
            final double e0 = ffe.energy(true, false);
            esvSystem.setLambda(i, 1.0);
            final double e1 = ffe.energy(true, false);
            SB.logfn(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f", i, e1, e0, e1-e0);
            
            esvSystem.setLambda(i, initialLambda[i] - step);
            final double totalLow = ffe.energy(true, false);
            final double vdwLow = ffe.getVanDerWaalsEnergy();
            final double realLow = ffe.getPermanentRealSpaceEnergy();
            final double selfLow = ffe.getPermanentReciprocalSelfEnergy();
            final double recipLow = ffe.getPermanentReciprocalMpoleEnergy();
            
            esvSystem.setLambda(i, initialLambda[i] + step);
            final double totalHigh = ffe.energy(true, false);
            final double vdwHigh = ffe.getVanDerWaalsEnergy();
            final double realHigh = ffe.getPermanentRealSpaceEnergy();
            final double selfHigh = ffe.getPermanentReciprocalSelfEnergy();
            final double recipHigh = ffe.getPermanentReciprocalMpoleEnergy();
            
            esvSystem.setLambda(i, initialLambda[i]);
            final double vdwAnalytic = esvSystem.getdVdwdL(i);
            final double realAnalytic = esvSystem.getdPermRealdL(i);
            final double selfAnalytic = esvSystem.getdPermRecipSelf_dL(i);
            final double recipAnalytic = esvSystem.getdPermRecipMpole_dL(i);
            final double totalAnalytic = esvSystem.getdEdL(i, true);
            
            final double vdwNumeric = (vdwHigh - vdwLow) / (2 * step);
            final double vdwError = Math.abs(vdwNumeric - vdwAnalytic);            
            final double realNumeric = (realHigh - realLow) / width;
            final double realError = Math.abs(realNumeric - realAnalytic);            
            final double selfNumeric = (selfHigh - selfLow) / width;
            final double selfError = Math.abs(selfNumeric - selfAnalytic);            
            final double recipNumeric = (recipHigh - recipLow) / width;
            final double recipError = Math.abs(recipNumeric - recipAnalytic);            
            final double totalNumeric = (totalHigh - totalLow) / width;
            final double totalError = Math.abs(totalNumeric - totalAnalytic);
            
            String esvName = esvSystem.getEsv(i).getName();
            SB.nlogf("  Finite Difference Test: %s", esvName);
            SB.nlogf(" *************************");
            SB.logf(StringUtils.repeat("*", esvName.length() + 1));
            SB.logfn("   %10s   %10s    %10s    %10s    %10s", "vdw", "permReal", "permSelf", "permRecip", "total");
            SB.logfn(" Numeric  Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    initialLambda[i], vdwNumeric, realNumeric, selfNumeric, recipNumeric, totalNumeric);
            SB.logfn(" Analytic Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    initialLambda[i], vdwAnalytic, realAnalytic, selfAnalytic, recipAnalytic, totalAnalytic);
            SB.logfn(" Error:                               %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    vdwError, realError, selfError, recipError, totalError);
            SB.force();
            
            assertEquals("VdW Deriv", 0.0, vdwError, tolerance);
            assertEquals("PermReal Deriv", 0.0, realError, tolerance);
            assertEquals("PermSelf Deriv", 0.0, selfError, tolerance);
            assertEquals("PermMpole Deriv", 0.0, recipError, tolerance);
            assertEquals("Total Deriv", 0.0, totalError, tolerance);
        }
    }
    
    /**
     * Verify that a lys-lys system with two ESVs can exactly reproduce the energy
     * yielded by vanilla energy() calls on mutated PDB files.
     */
    @Test
    public void testEndStates() {
        logMaster.setLevel(Level.WARNING);
        setup(true);        
        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 0.0);
        final double esvTot00 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvBias00 = esvSystem.getBiasEnergy();
        SB.force("Tot,Bias: %g %g", esvTot00, esvBias00);
        final double esvVdw00 = ffe.getVanDerWaalsEnergy();
        final double esvPme00 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp00 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();
        
        esvSystem.setLambda("A", 0.0);
        esvSystem.setLambda("B", 1.0);
        final double esvTot01 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvBias01 = esvSystem.getBiasEnergy();
        SB.force("Tot,Bias: %g %g", esvTot01, esvBias01);
        final double esvVdw01 = ffe.getVanDerWaalsEnergy();
        final double esvPme01 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp01 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 0.0);
        final double esvTot10 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvBias10 = esvSystem.getBiasEnergy();
        SB.force("Tot,Bias: %g %g", esvTot10, esvBias10);
        final double esvVdw10 = ffe.getVanDerWaalsEnergy();
        final double esvPme10 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp10 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        esvSystem.setLambda("A", 1.0);
        esvSystem.setLambda("B", 1.0);
        final double esvTot11 = ffe.energy(true, false) - esvSystem.getBiasEnergy();
        final double esvBias11 = esvSystem.getBiasEnergy();
        SB.force("Tot,Bias: %g %g", esvTot11, esvBias11);
        final double esvVdw11 = ffe.getVanDerWaalsEnergy();
        final double esvPme11 = ffe.getPermanentRealSpaceEnergy();
        final double esvRcp11 = ffe.getPermanentReciprocalSelfEnergy() + ffe.getPermanentReciprocalMpoleEnergy();

        StringBuilder ext = new StringBuilder();
        ext.append(format("%-9s %-7s %10.5f\n", "Total", "0-0", esvTot00));
        ext.append(format("%-9s %-7s %10.5f\n", "Total", "0-1", esvTot01));
        ext.append(format("%-9s %-7s %10.5f\n", "Total", "1-0", esvTot10));
        ext.append(format("%-9s %-7s %10.5f\n", "Total", "1-1", esvTot11));
        if (testVdw) {
            ext.append(format("%-9s %-7s %10.5f\n", "VanWaals", "0-0", esvVdw00));
            ext.append(format("%-9s %-7s %10.5f\n", "VanWaals", "0-1", esvVdw01));
            ext.append(format("%-9s %-7s %10.5f\n", "VanWaals", "1-0", esvVdw10));
            ext.append(format("%-9s %-7s %10.5f\n", "VanWaals", "1-1", esvVdw11));
        }
        if (testPermReal) {
            ext.append(format("%-9s %-7s %10.5f\n", "PermReal", "0-0", esvPme00));
            ext.append(format("%-9s %-7s %10.5f\n", "PermReal", "0-1", esvPme01));
            ext.append(format("%-9s %-7s %10.5f\n", "PermReal", "1-0", esvPme10));
            ext.append(format("%-9s %-7s %10.5f\n", "PermReal", "1-1", esvPme11));
        }
        if (testReciprocal) {
            ext.append(format("%-9s %-7s %10.5f\n", "PermRecip", "0-0", esvRcp00));
            ext.append(format("%-9s %-7s %10.5f\n", "PermRecip", "0-1", esvRcp01));
            ext.append(format("%-9s %-7s %10.5f\n", "PermRecip", "1-0", esvRcp10));
            ext.append(format("%-9s %-7s %10.5f\n", "PermRecip", "1-1", esvRcp11));
        }

        StringBuilder tot = new StringBuilder();
        StringBuilder vdw = new StringBuilder();
        StringBuilder pme = new StringBuilder();
        StringBuilder rcp = new StringBuilder();

        /* Open vanilla end states. */
        MolecularAssembly vanilla;
        ForceFieldEnergy vanillaFFE;
        String state;
        
        /* lyd-lyd */
        vanilla = openResource(filename00, true);
        vanillaFFE = vanilla.getPotentialEnergy();
        state = filename00;
        {
            final double Utot = vanillaFFE.energy(true, false);
            assertEquals("total00", energy00, Utot, tolerance);
            assertEquals("total00", Utot, esvTot00, tolerance);
            tot.append(format("%-9s %-7.7s %10.5f\n", "Total", state, Utot));
        }
        if (testVdw) {
            final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
            assertEquals("vdw00", vanWaals, esvVdw00, tolerance);
            vdw.append(format("%-9s %-7.7s %10.5f\n", "VanWaals", state, vanWaals));
        }
        if (testPermReal) {
            final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
            assertEquals("permReal00", permReal, esvPme00, tolerance);
            pme.append(format("%-9s %-7s %10.5f\n", "PermReal", state, permReal));
        }
        if (testReciprocal) {
            final double permRecip = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
            assertEquals("permRecip00", permRecip, esvRcp00, tolerance);
            rcp.append(format("%-9s %-7s %10.5f\n", "PermRecip", state, permRecip));
        }
        utils.close(vanilla);
        
        /* lyd-lys */
        vanilla = openResource(filename01, true);
        vanillaFFE = vanilla.getPotentialEnergy();
        state = filename01;
        {
            final double Utot = vanillaFFE.energy(true, false);
            assertEquals("total01", energy01, Utot, tolerance);
            assertEquals("total01", Utot, esvTot01, tolerance);
            tot.append(format("%-9s %-7.7s %10.5f\n", "Total", state, Utot));
        }
        if (testVdw) {
            final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
            assertEquals("vdw01", vanWaals, esvVdw01, tolerance);
            vdw.append(format("%-9s %-7.7s %10.5f\n", "VanWaals", state, vanWaals));
        }
        if (testPermReal) {
            final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
            assertEquals("permReal01", permReal, esvPme01, tolerance);
            pme.append(format("%-9s %-7s %10.5f\n", "PermReal", state, permReal));
        }
        if (testReciprocal) {
            final double permRecip = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
            assertEquals("permRecip01", permRecip, esvRcp01, tolerance);
            rcp.append(format("%-9s %-7s %10.5f\n", "PermRecip", state, permRecip));
        }
        utils.close(vanilla);
        
        /* lys-lyd */
        vanilla = openResource(filename10, true);
        vanillaFFE = vanilla.getPotentialEnergy();
        state = filename10;
        {
            final double Utot = vanillaFFE.energy(true, false);
            assertEquals("total10", energy10, Utot, tolerance);
            assertEquals("total10", Utot, esvTot10, tolerance);
            tot.append(format("%-9s %-7.7s %10.5f\n", "Total", state, Utot));
        }
        if (testVdw) {
            final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
            assertEquals("vdw10", vanWaals, esvVdw10, tolerance);
            vdw.append(format("%-9s %-7.7s %10.5f\n", "VanWaals", state, vanWaals));
        }
        if (testPermReal) {
            final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
            assertEquals("permReal10", permReal, esvPme10, tolerance);
            pme.append(format("%-9s %-7s %10.5f\n", "PermReal", state, permReal));
        }
        if (testReciprocal) {
            final double permRecip = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
            assertEquals("permRecip10", permRecip, esvRcp10, tolerance);
            rcp.append(format("%-9s %-7s %10.5f\n", "PermRecip", state, permRecip));
        }
        utils.close(vanilla);
        
        /* lys-lys */
        vanilla = openResource(filename11, true);
        vanillaFFE = vanilla.getPotentialEnergy();
        state = filename11;
        {
            final double Utot = vanillaFFE.energy(true, false);
            assertEquals("total11", energy11, Utot, tolerance);
            assertEquals("total11", Utot, esvTot11, tolerance);
            tot.append(format("%-9s %-7.7s %10.5f\n", "Total", state, Utot));
        }
        if (testVdw) {
            final double vanWaals = vanillaFFE.getVanDerWaalsEnergy();
            assertEquals("vdw11", vanWaals, esvVdw11, tolerance);
            vdw.append(format("%-9s %-7.7s %10.5f\n", "VanWaals", state, vanWaals));
        }
        if (testPermReal) {
            final double permReal = vanillaFFE.getPermanentRealSpaceEnergy();
            assertEquals("permReal11", permReal, esvPme11, tolerance);
            pme.append(format("%-9s %-7s %10.5f\n", "PermReal", state, permReal));
        }
        if (testReciprocal) {
            final double permRecip = vanillaFFE.getPermanentReciprocalSelfEnergy() + vanillaFFE.getPermanentReciprocalMpoleEnergy();
            assertEquals("permRecip11", permRecip, esvRcp11, tolerance);
            rcp.append(format("%-9s %-7s %10.5f\n", "PermRecip", state, permRecip));
        }
        utils.close(vanilla);

        StringBuilder van = new StringBuilder();
        van.append(tot.toString());
        if (testVdw) van.append(vdw.toString());
        if (testPermReal) van.append(pme.toString());
        if (testReciprocal) van.append(rcp.toString());
        String[] vanLines = van.toString().split("\\n");
        String[] extLines = ext.toString().split("\\n");
        StringBuilder sb = new StringBuilder();
        SB.nlogfn("  Two-site End State Analysis ");
        SB.logfn(" ***************************** ");
        SB.logfn(" %-30s     %-30s","Extended System (1 Assembly)","Mutated Files (4 Assemblies)");
        for (int i = 0; i < extLines.length; i++) {
            SB.logfn(" %-30s     %-30s", extLines[i], vanLines[i]);
        }
        SB.force();
    }
    
    /**
     * Numerically ensure that the energy and lambda derivatives are smooth all 
     * along both ESV coordinates in the dilysine system.
     */
    @Test
    public void testSmoothness() {
        logMaster.setLevel(Level.WARNING);
        setup(true);        
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
            double evA = idxA / 10.0;
            for (int idxB = 0; idxB <= 10; idxB++) {
                final double evB = idxB / 10.0;
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
        StringBuilder sb = new StringBuilder();
        SB.nlogfn("  Smoothness Verification: Total ");
        SB.logfn(" ******************************** ");
        printAsTable(totalEnergies, "U_Total", SB);
        printAsTable(totalDerivsA, "dU_dEsvA", SB);
        printAsTable(totalDerivsB, "dU_dEsvB", SB);
        SB.force();
        
        /* TODO: improve upon the following arbitrary definition of approximate smoothness. */
        final int testRow = 2;
        final double min = Arrays.stream(totalEnergies[testRow]).min().getAsDouble();
        final double max = Arrays.stream(totalEnergies[testRow]).max().getAsDouble();
        final double maxChange = (max - min) * 0.5;
        final int maxDirectionSignChanges = 2;
        int directionSignChanges = 0;
        for (int idxB = 0; idxB < 10; idxB++) {
            final double here = totalEnergies[testRow][idxB];
            final double next = totalEnergies[testRow][idxB+1];
            final double absChange = Math.abs(next - here);
            if (idxB > 0) {
                final double prev = totalEnergies[testRow][idxB-1];
                if (Math.signum(next - here) != Math.signum(here - prev)) {
                    directionSignChanges++;
                }
            }
            if (absChange > maxChange) {
                Assert.fail(format("Failed max_change smoothness criterion: %g > %g. (change:%g->%g, range:{%g,%g})",
                        absChange, maxChange, min, max,
                        totalEnergies[testRow][idxB], totalEnergies[testRow][idxB+1]));
            }
        }
        if (directionSignChanges > maxDirectionSignChanges) {
            Assert.fail(format("Failed direction_changes smoothness criterion: %d > %d.",
                    directionSignChanges, maxDirectionSignChanges));
        }
        
        if (testVdw && verbose) {
            sb = new StringBuilder();
            SB.nlogfn("  Smoothness Verification: VdW ");
            SB.logfn(" ****************************** ");
            printAsTable(vdwEnergies, "vanWaals", SB);
            printAsTable(vdwDerivsA, "dVdw_dA", SB);
            printAsTable(vdwDerivsB, "dVdw_dB", SB);
            SB.force();
        }
        if (testPermReal && verbose) {
            sb = new StringBuilder();
            SB.nlogfn("  Smoothness Verification: PermReal ");
            SB.logfn(" *********************************** ");
            printAsTable(permRealEnergies, "permReal", SB);
            printAsTable(permRealDerivsA, "dPRealdA", SB);
            printAsTable(permRealDerivsB, "dPRealdB", SB);
            SB.force();
        }
        if (testReciprocal && verbose) {
            sb = new StringBuilder();
            SB.nlogfn("  Smoothness Verification: PermRecip ");
            SB.logfn(" ************************************ ");
            printAsTable(permRecipEnergies, "permRcp", SB);
            printAsTable(permRecipDerivsA, "dPRcp_dA", SB);
            printAsTable(permRecipDerivsB, "dPRcp_dB", SB);
            SB.force();
        }
    }
    
    private void printAsTable(double[][] values, String title) {
        printAsTable(values, title, null);
    }
    
    private SBLogger printAsTable(double[][] values, String title, SBLogger sb) {
        if (SB == null) {
            sb = new SBLogger();
        }
        if (title == null) title = "Title?";
        sb.logf(" %8s %8.1f %8.1f %8.1f %8.1f %8s %8s %8s %8.1f %8.1f %8.1f %8.1f", 
                title, 0.0, 0.1, 0.2, 0.3, "<   ", " EvA  ", ">   ", 0.7, 0.8, 0.9, 1.0);
        for (int idxA = 0; idxA <= 10; idxA++) {
            double evA = idxA / 10.0;
            switch (idxA) {
                case 4:
                    sb.nlogf(" %8s", " ^ ");
                    break;
                case 5:
                    sb.nlogf(" %8s", "EvB");
                    break;
                case 6:
                    sb.nlogf(" %8s", " v ");
                    break;
                default:
                    sb.nlogf(" %8.1f", evA);
            }
            for (int idxB = 0; idxB <= 10; idxB++) {
                double evB = idxB / 10.0;
                String value = format("%8.4f", values[idxA][idxB]);
                int max = (value.length() < 8) ? value.length() : 8;
                sb.logf(" %8s", value.substring(0,max));
            }
        }
        sb.nl();
        return sb;
    }
}
