package ffx.potential.utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Optional;
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
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.TitrationESV;
import ffx.potential.extended.TitrationUtils;

/**
 *
 * @author slucore
 */
@RunWith(Parameterized.class)
public class ExtendedVariableTest {
    
    private static final Logger logger = Logger.getLogger(ExtendedVariableTest.class.getName());
    private final PotentialsUtils utils = new PotentialsUtils();
    
    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            /* Arguments to the class constructor. */
            {"A2", "A4", 1.0e-5, 
             "lyd-lyd-cryst.pdb", 65.3397,
             "lyd-lys-cryst.pdb", 78.7910,
             "lys-lyd-cryst.pdb", 76.9447,
             "lys-lys-cryst.pdb", 97.0408}
        });
    }
    
    private static final boolean testVdw = true;
    private static final boolean testPermReal = true;
    private static final boolean testReciprocal = true;
    private static final boolean testPolarization = false;
    private static final double cutoff = 5.0;
    private static final boolean verbose = false;
    private static String resourcePrefix = "ffx/potential/structures/";
    
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
    
    private final double lambda = 0.28;
    private final double step = 0.0001;
    private final double biasMag = 1.0;
    private final double pH = 7.4;
    private final int numESVs = 2;
    
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
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}

    private MolecularAssembly openResource(String filename) {
        return openResource(filename, false);
    }
    
    private MolecularAssembly openResource(String filename, boolean quietly) {
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(resourcePrefix + filename).getPath());
        return (quietly) ? utils.openQuietly(structure) : utils.open(structure);
    }
    
    @Test
    public void test() {
        mola = openResource(filename11, false);
        ffe = mola.getPotentialEnergy();
        esvSystem = new ExtendedSystem(mola);
        
        Polymer[] polymers = mola.getChains();
        String[] rlTokens = {residA, residB};
        double temperature = 298.15;
        for (int i = 0; i < numESVs; i++) {
            Character chainID = rlTokens[i].charAt(0);
            int resNum = Integer.parseInt(rlTokens[i].substring(1));
            Optional<Residue> target = Optional.empty();
            for (Polymer p : polymers) {
                if (p.getChainID().equals(chainID)) {
                    target = p.getResidues().stream()
                            .filter((Residue res) -> res.getResidueNumber() == resNum)
                            .findFirst();
                    break;
                }
            }
            if (!target.isPresent()) {
                logger.log(Level.SEVERE, "Couldn't find target residue {0}", rlTokens[i]);
            }

            SB.logfp("Residue: %s", target.get());
            MultiResidue titrating = TitrationUtils.titrationFactory(mola, target.get());
            
            TitrationESV esv = new TitrationESV(esvSystem.getConfig(), titrating, pH, biasMag);
            esvSystem.addVariable(esv);
        }
        
        // Attach populated esvSystem to the potential.
        ffe.attachExtendedSystem(esvSystem);
        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        ffe.getCrystal().setSpecialPositionCutoff(0.0);
        
        testDerivatives();
        testEndStates();
        testSmoothness();
    }
    
    /**
     * Analytic Derivative vs Finite Difference: VdW,PermReal,PermRecip,Total
     */
    private void testDerivatives() {
        for (int i = 0; i < numESVs; i++) {
            esvSystem.setLambda(i, lambda);
        }
        for (int i = 0; i < numESVs; i++) {
            esvSystem.setLambda(i, 0.0);
            final double e0 = ffe.energy(true, false);
            esvSystem.setLambda(i, 1.0);
            final double e1 = ffe.energy(true, false);
            SB.logfn(" ESV%d> E(1),E(0),diff:  %14.6f - %14.6f = %14.6f", i, e1, e0, e1-e0);

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
            final double vdwNumeric = (vdwHigh - vdwLow) / (2 * step);
            final double vdwError = Math.abs(vdwNumeric - vdwAnalytic);
            final double realAnalytic = esvSystem.getdPermRealdL(i);
            final double realNumeric = (realHigh - realLow) / (2 * step);
            final double realError = Math.abs(realNumeric - realAnalytic);
            final double selfAnalytic = esvSystem.getdPermRecipSelf_dL(i);
            final double selfNumeric = (selfHigh - selfLow) / (2 * step);
            final double selfError = Math.abs(selfNumeric - selfAnalytic);
            final double recipAnalytic = esvSystem.getdPermRecipMpole_dL(i);
            final double recipNumeric = (recipHigh - recipLow) / (2 * step);
            final double recipError = Math.abs(recipNumeric - recipAnalytic);
            final double totalAnalytic = esvSystem.getdEdL(i, true);
            final double totalNumeric = (totalHigh - totalLow) / (2 * step);
            final double totalError = Math.abs(totalNumeric - totalAnalytic);

            String esvName = esvSystem.getEsv(i).getName();
            SB.nlogf("  Finite Difference Test: %s", esvName);
            SB.nlogf(" *************************");
            SB.logf(StringUtils.repeat("*", esvName.length() + 1));
            SB.logfn("   %10s   %10s    %10s    %10s    %10s", "vdw", "permReal", "permSelf", "permRecip", "total");
            SB.logfn(" Numeric  Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    lambda, vdwNumeric, realNumeric, selfNumeric, recipNumeric, totalNumeric);
            SB.logfn(" Analytic Derivatives (@lambda %4.2f): %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    lambda, vdwAnalytic, realAnalytic, selfAnalytic, recipAnalytic, totalAnalytic);
            SB.logfn(" Error:                               %+12.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f",
                    vdwError, realError, selfError, recipError, totalError);
            SB.print();
            
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
    private void testEndStates() {
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
        vanillaFFE.energy(true, false);
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
        vanillaFFE.energy(true, false);
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
        vanillaFFE.energy(true, false);
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
        vanillaFFE.energy(true, false);
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
        if (testVdw) van.append(vdw.toString());
        if (testPermReal) van.append(pme.toString());
        if (testReciprocal) van.append(rcp.toString());
        String[] vanLines = van.toString().split("\\n");
        String[] extLines = ext.toString().split("\\n");
        StringBuilder sb = new StringBuilder();
        sb.append(format("\n  Two-site End State Analysis \n"));
        sb.append(format(" ***************************** \n"));
        sb.append(format(" %-30s     %-30s\n","Extended System (1 Assembly)","Mutated Files (4 Assemblies)"));
        for (int i = 0; i < extLines.length; i++) {
            sb.append(format(" %-30s     %-30s\n", extLines[i], vanLines[i]));
        }
        logger.info(sb.toString());
    }
    
    /**
     * Numerically ensure that the energy and lambda derivatives are smooth all 
     * along both ESV coordinates in the dilysine system.
     */
    private void testSmoothness() {
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
        sb.append(format("\n  Smoothness Verification: Total \n"));
        sb.append(format(" ******************************** \n"));
        printAsTable(totalEnergies, "U_Total", sb);
        printAsTable(totalDerivsA, "dU_dEsvA", sb);
        printAsTable(totalDerivsB, "dU_dEsvB", sb);
        logger.info(sb.toString());
        
        /* TODO: improve upon the following arbitrary definition of approximate smoothness. */
        final int testRow = 2;
        final double min = Arrays.stream(totalEnergies[testRow]).min().getAsDouble();
        final double max = Arrays.stream(totalEnergies[testRow]).max().getAsDouble();
        final double maxChange = (max - min) * 0.5;
        final double maxDirectionSignChanges = 2;
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
            sb.append(format("\n  Smoothness Verification: VdW "));
            sb.append(format(" ****************************** "));
            printAsTable(vdwEnergies, "vanWaals", sb);
            printAsTable(vdwDerivsA, "dVdw_dA", sb);
            printAsTable(vdwDerivsB, "dVdw_dB", sb);
            logger.info(sb.toString());
        }
        if (testPermReal && verbose) {
            sb = new StringBuilder();
            sb.append(format("\n  Smoothness Verification: PermReal "));
            sb.append(format(" *********************************** "));
            printAsTable(permRealEnergies, "permReal", sb);
            printAsTable(permRealDerivsA, "dPRealdA", sb);
            printAsTable(permRealDerivsB, "dPRealdB", sb);
            logger.info(sb.toString());
        }
        if (testReciprocal && verbose) {
            sb = new StringBuilder();
            sb.append(format("\n  Smoothness Verification: PermRecip "));
            sb.append(format(" ************************************ "));
            printAsTable(permRecipEnergies, "permRcp", sb);
            printAsTable(permRecipDerivsA, "dPRcp_dA", sb);
            printAsTable(permRecipDerivsB, "dPRcp_dB", sb);
            logger.info(sb.toString());
        }
    }
    
    private void printAsTable(double[][] values, String title) {
        printAsTable(values, title, null);
    }
    
    private StringBuilder printAsTable(double[][] values, String title, StringBuilder sb) {
        if (sb == null) {
            sb = new StringBuilder();
        }
        if (title == null) title = "Title?";
        sb.append(format(" %8s %8.1f %8.1f %8.1f %8.1f %8s %8s %8s %8.1f %8.1f %8.1f %8.1f", 
                title, 0.0, 0.1, 0.2, 0.3, "<   ", " EvA  ", ">   ", 0.7, 0.8, 0.9, 1.0));
        for (int idxA = 0; idxA <= 10; idxA++) {
            double evA = idxA / 10.0;
            switch (idxA) {
                case 4:
                    sb.append(format("\n %8s", " ^ "));
                    break;
                case 5:
                    sb.append(format("\n %8s", "EvB"));
                    break;
                case 6:
                    sb.append(format("\n %8s", " v "));
                    break;
                default:
                    sb.append(format("\n %8.1f", evA));
            }
            for (int idxB = 0; idxB <= 10; idxB++) {
                double evB = idxB / 10.0;
                String value = format("%8.4f", values[idxA][idxB]);
                int max = (value.length() < 8) ? value.length() : 8;
                sb.append(format(" %8s", value.substring(0,max)));
            }
        }
        sb.append(format("\n"));
        return sb;
    }
}
