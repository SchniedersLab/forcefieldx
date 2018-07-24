/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import static org.junit.Assert.assertEquals;

import edu.rit.pj.Comm;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.utils.PotentialsUtils;

/**
 * Test the Goldstein elimination criteria for both self and pair eliminations.
 *
 * @author Mallory R. Tollefson
 * @author Claire E. OConnell
 */
@RunWith(Parameterized.class)
public class RotamerOptimizationTest {

    private static final Logger logger = Logger.getLogger(RotamerOptimizationTest.class.getName());

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {"Chignolin Direct with Orig Rot - No Pruning (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun0.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        0, // Pruning Level.
                        true, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492687, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        1018.8163922682276, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - Singles Pruning (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        1, // Pruning Level.
                        true, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492684, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - Pairs Pruning (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun2.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        2, // Pruning Level.
                        true, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492684, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - 3-body (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.3body.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        1, // Pruning Level.
                        true, // Goldstein Elimination.
                        true, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492687, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        true, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        -197.51830687083316, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - No Pruning (DEE)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun0.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        0, // Pruning Level.
                        false, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492687, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        1018.8163922682276, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - Singles Pruning (DEE)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        1, // Pruning Level.
                        false, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492684, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - Pairs Pruning (DEE)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun2.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        2, // Pruning Level.
                        false, // Goldstein Elimination.
                        false, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492684, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        false, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        0.0, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                },
                {"Chignolin Direct with Orig Rot - 3-body (DEE)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.3body.residues1-4.restart",
                        1, // Start Residue.
                        4, // End Residue.
                        1, // Pruning Level.
                        false, // Goldstein Elimination.
                        true, // Use 3-body Energies.
                        true, // Use Original Rotamers.
                        true, // Do Overall Opt.
                        -203.10632181492687, // Expected Energy.
                        true, // Do Self-Energy Opt.
                        -203.10632181492687, // Expected Self-Energy.
                        true, // Do Pair-Energy Opt.
                        1, // Pair residue
                        -186.71062663844978, // Expected Pair-Energy.
                        true, // Do Trimer-Energy Opt.
                        1, // Trimer residue 1.
                        2, // Trimer residue 2.
                        -197.51830687083316, // Expected trimer energy.
                        1.0e-3 // Energy Tolerance.
                }

                /**
                 * {
                 * "Trpcage Direct with Orig Rot (Goldstein)",
                 * "ffx/algorithms/structures/trpcage.pdb",
                 * "ffx/algorithms/structures/trpcage.direct.orig.restart", 0, //
                 * Pruning Level. true, // Goldstein Elimination. false, // Use 3-body
                 * Energies. true, // Use Original Rotamers. false, // Do Overall Opt.
                 * -628.143879, // Expected Energy. true, // Do Self-Energy Opt.
                 * 22086.367273, // Expected Self-Energy. false, // Do Pair-Energy Opt.
                 * 5, // Pair residue 1280.865248, // Expected Pair-Energy. false, // Do
                 * Trimer-Energy Opt. 5, // Trimer residue 1. 6, // Trimer residue 2.
                 * 0.0, // Expected trimer energy. 1.0e-3 // Energy Tolerance. }, {
                 * "Trpcage Direct with Orig Rot (DEE, Prune=1)",
                 * "ffx/algorithms/structures/trpcage.pdb",
                 * "ffx/algorithms/structures/trpcage.direct.orig.restart", 1, //
                 * Pruning Level. false, // Goldstein Elimination. false, // Use 3-body
                 * Energies. true, // Use Original Rotamers. false, // Do Overall Opt.
                 * (false because DEE leaves too many permutations) -628.143879, //
                 * Expected Energy. true, // Do Self-Energy Opt. 22086.367273, //
                 * Expected Self-Energy. true, // Do Pair-Energy Opt. 5, // Pair residue
                 * 1280.865248, // Expected Pair-Energy. false, // Do Trimer-Energy Opt.
                 * 5, // Trimer residue 1. 6, // Trimer residue 2. 0.0, // Expected
                 * trimer energy. 1.0e-3 // Energy Tolerance. }, { "Trpcage Direct with
                 * Orig Rot (DEE, Prune=2)", "ffx/algorithms/structures/trpcage.pdb",
                 * "ffx/algorithms/structures/trpcage.direct.orig.restart", 2, //
                 * Pruning Level. false, // Goldstein Elimination. false, // Use 3-body
                 * Energies. true, // Use Original Rotamers. true, // Do Overall Opt.
                 * (false because DEE leaves too many permutations) -628.143879, //
                 * Expected Energy. true, // Do Self-Energy Opt. -616.237522, //
                 * Expected Self-Energy. true, // Do Pair-Energy Opt. 5, // Pair residue
                 * -563.929018, // Expected Pair-Energy. false, // Do Trimer-Energy Opt.
                 * 5, // Trimer residue 1. 6, // Trimer residue 2. 0.0, // Expected
                 * trimer energy. 1.0e-3 // Energy Tolerance. }, { "Trpcage Direct
                 * (Goldstein)", "ffx/algorithms/structures/trpcage.pdb",
                 * "ffx/algorithms/structures/trpcage.direct.3body.restart", 1, //
                 * Pruning Level. true, // Goldstein Elimination. true, // Use 3-body
                 * Energies. true, // Use Original Rotamers. true, // Do Overall Opt.
                 * -628.143879, // Expected Energy. true, // Do Self-Energy Opt.
                 * 22086.367273, // Expected Self-Energy. true, // Do Pair-Energy Opt.
                 * 5, // Pair residue. 1280.865247, // Expected Pair-Energy. true, // Do
                 * Trimer-Energy Opt. 5, // Trimer residue 1. 6, // Trimer residue 2.
                 * 1276.608834, // Expected trimer energy. 1.0e-3 // Energy Tolerance. }
                */
        });
    }

    String info;
    String filename;
    String restartName;
    int startResID;
    int endResID;
    int pruningLevel;
    boolean useGoldstein;
    boolean useThreeBody;
    boolean useOriginalRotamers;
    boolean doOverallOpt;
    double expectedEnergy;
    boolean doSelfOpt;
    double expectedSelfEnergy;
    boolean doPairOpt;
    int pairResidue;
    double expectedPairEnergy;
    boolean doTripleOpt;
    int tripleResidue1;
    int tripleResidue2;
    double expectedTripleEnergy;
    double tolerance;
    File structure;
    File restartFile;
    MolecularAssembly molecularAssembly;
    ForceFieldEnergy forceFieldEnergy;

    public RotamerOptimizationTest(String info,
                                   String filename,
                                   String restartName,
                                   int startResID,
                                   int endResID,
                                   int pruningLevel,
                                   boolean useGoldstein,
                                   boolean useThreeBody,
                                   boolean useOriginalRotamers,
                                   boolean doOverallOpt,
                                   double expectedEnergy,
                                   boolean doSelfOpt,
                                   double expectedSelfEnergy,
                                   boolean doPairOpt,
                                   int pairResidue,
                                   double expectedPairEnergy,
                                   boolean doTripleOpt,
                                   int tripleResidue1,
                                   int tripleResidue2,
                                   double expectedTripleEnergy,
                                   double tolerance) {
        this.info = info;
        this.filename = filename;
        this.restartName = restartName;
        this.startResID = startResID;
        this.endResID = endResID;
        this.pruningLevel = pruningLevel;
        this.useGoldstein = useGoldstein;
        this.useThreeBody = useThreeBody;
        this.useOriginalRotamers = useOriginalRotamers;
        this.doOverallOpt = doOverallOpt;
        this.expectedEnergy = expectedEnergy;
        this.doSelfOpt = doSelfOpt;
        this.expectedSelfEnergy = expectedSelfEnergy;
        this.doPairOpt = doPairOpt;
        this.pairResidue = pairResidue;
        this.expectedPairEnergy = expectedPairEnergy;
        this.doTripleOpt = doTripleOpt;
        this.tripleResidue1 = tripleResidue1;
        this.tripleResidue2 = tripleResidue2;
        this.expectedTripleEnergy = expectedTripleEnergy;
        this.tolerance = tolerance;
    }

    public void load() {
        /**
         * Load the test system.
         */
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());
        restartFile = new File(cl.getResource(restartName).getPath());
        PotentialsUtils potentialUtils = new PotentialsUtils();
        molecularAssembly = potentialUtils.openQuietly(structure.getAbsolutePath());
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    }

    @BeforeClass
    public static void beforeClass() {
        // Initialize Parallel Java
        try {
            Comm.world();
        } catch (IllegalStateException ise) {
            try {
                String args[] = new String[0];
                Comm.init(args);
            } catch (Exception e) {
                String message = " Exception starting up the Parallel Java communication layer.";
                logger.log(Level.WARNING, message, e.toString());
                message = " Skipping rotamer optimization test.";
                logger.log(Level.WARNING, message, e.toString());
            }
        }
    }

    @After
    public void after() {
        forceFieldEnergy.destroy();
        System.gc();
    }

    @Test
    public void testSelfEnergyElimination() {
        // Load the test system.
        load();

        // Run the optimization.
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(useOriginalRotamers);

        int counter = 1;
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = molecularAssembly.getChains();
        int nPolymers = polymers.length;
        for (int p = 0; p < nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            for (int i = 0; i < endResID; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers(rLib);
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    }
                    if (counter >= startResID) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }
        Residue[] residues = residueList.toArray(new Residue[residueList.size()]);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(useThreeBody);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setPruning(pruningLevel);
        rotamerOptimization.setEnergyRestartFile(restartFile);
        rotamerOptimization.setResidues(residueList);

        double energy;
        int nRes = residueList.size();
        if (doOverallOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestOverallOpt(true);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            //System.out.println("The expected overall energy is: " + energy);
            assertEquals(info + " Total Energy", expectedEnergy, energy, tolerance);
        }

        if (doSelfOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestSelfEnergyEliminations(true);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            //System.out.println("The expected self is: " + energy);
            assertEquals(info + " Self-Energy", expectedSelfEnergy, energy, tolerance);

            // Check that optimized rotamers are equivalent to the lowest self-energy of each residue.
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            // Loop over all residues
            for (int i = 0; i < nRes; i++) {
                Residue res = residueList.get(i);
                Rotamer[] rotI = res.getRotamers(rLib);
                int nRot = rotI.length;

                int rotCounter = 0;
                while (rotCounter < nRot && rotamerOptimization.checkPrunedSingles(i, rotCounter)) {
                    rotCounter++;
                }

                double lowEnergy = rotamerOptimization.getSelf(i, rotCounter);
                int bestRot = rotCounter;
                for (int ri = 1; ri < nRot; ri++) {
                    if (rotamerOptimization.checkPrunedSingles(i, ri)) {
                        continue;
                    } else {
                        double selfEnergy = rotamerOptimization.getSelf(i, ri);
                        if (selfEnergy < lowEnergy) {
                            lowEnergy = selfEnergy;
                            bestRot = ri;
                        }
                    }
                }
                assertEquals(String.format(" %s Self-Energy of residue %d", info, i), optimum[i], bestRot);
            }
        }

        if (doPairOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestPairEnergyEliminations(pairResidue);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            assertEquals(info + " Pair-Energy", expectedPairEnergy, energy, tolerance);

            // Check that optimized rotamers are equivalent to the lowest 2-Body energy sum for the "pairResidue".
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            Residue resI = residueList.get(pairResidue);
            Rotamer rotI[] = resI.getRotamers(rLib);
            int ni = rotI.length;

            double minEnergy = Double.POSITIVE_INFINITY;
            int bestRotI = -1;

            // Loop over the pairResidue rotamers to find its lowest energy rotamer.
            for (int ri = 0; ri < ni; ri++) {
                double energyForRi = 0.0;
                if (rotamerOptimization.checkPrunedSingles(pairResidue, ri)) {
                    continue;
                }
                // Loop over residue J
                for (int j = 0; j < nRes; j++) {
                    if (j == pairResidue) {
                        continue;
                    }
                    Residue resJ = residueList.get(j);
                    Rotamer[] rotJ = resJ.getRotamers(rLib);
                    int nRot = rotJ.length;

                    int rj = 0;
                    while (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, ri, j, rj)) {
                        if (++rj >= nRot) {
                            logger.warning("RJ is too large.");
                        }
                    }

                    double lowEnergy = rotamerOptimization.get2Body(pairResidue, ri, j, rj);

                    for (rj = 1; rj < nRot; rj++) {
                        if (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, ri, j, rj)) {
                            continue;
                        } else {
                            double pairEnergy = rotamerOptimization.get2Body(pairResidue, ri, j, rj);
                            if (pairEnergy < lowEnergy) {
                                lowEnergy = pairEnergy;
                            }
                        }
                    }
                    energyForRi += lowEnergy;
                }
                if (energyForRi < minEnergy) {
                    minEnergy = energyForRi;
                    bestRotI = ri;
                }
            }

            assertEquals(String.format(" %s Best 2-body energy sum for residue %d is with rotamer %d at %10.4f.", info, pairResidue, bestRotI, minEnergy),
                    optimum[pairResidue], bestRotI);

            // Given the minimum energy rotamer for "pairResidue" is "bestRotI", we can check selected rotamers for all other residues.
            for (int j = 0; j < nRes; j++) {
                if (j == pairResidue) {
                    continue;
                }
                Residue resJ = residueList.get(j);
                Rotamer[] rotJ = resJ.getRotamers(rLib);
                int nRotJ = rotJ.length;

                int rotCounter = 0;
                while (rotamerOptimization.checkPrunedPairs(pairResidue, bestRotI, j, rotCounter) && rotCounter < nRotJ) {
                    rotCounter++;
                }

                double lowEnergy = rotamerOptimization.get2Body(pairResidue, bestRotI, j, rotCounter);
                int bestRotJ = rotCounter;
                for (int rj = 1; rj < nRotJ; rj++) {
                    if (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, bestRotI, j, rj)) {
                        continue;
                    } else {
                        double pairEnergy = rotamerOptimization.get2Body(pairResidue, bestRotI, j, rj);
                        if (pairEnergy < lowEnergy) {
                            lowEnergy = pairEnergy;
                            bestRotJ = rj;
                        }
                    }
                }
                assertEquals(String.format(" %s Pair-Energy of residue (%d,%d) with residue %d", info, pairResidue, bestRotI, j), optimum[j], bestRotJ);
            }
        }

        //Test 3-Body Energy Eliminations.
        if (doTripleOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestTripleEnergyEliminations(tripleResidue1, tripleResidue2);
            try {
                energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
                assertEquals(info + " Triple-Energy", expectedTripleEnergy, energy, tolerance);
            } catch (Exception e) {
                e.fillInStackTrace();
                e.printStackTrace();
                logger.log(java.util.logging.Level.INFO, "Error in doTripleOpt", e);
            }

            //Check that optimized rotamers are equivalent to the lowest 3-body energy of each residue with the tripleResidue1 and 2.
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            //fix residue 1 and gets its rotamers
            Residue resI = residueList.get(tripleResidue1);
            Rotamer rotI[] = resI.getRotamers(rLib);
            int ni = rotI.length;

            //fix residue 2 and get its rotamers
            Residue resJ = residueList.get(tripleResidue2);
            Rotamer rotJ[] = resJ.getRotamers(rLib);
            int nj = rotJ.length;

            double minEnergyIJ = Double.POSITIVE_INFINITY;
            int bestRotI = -1;
            int bestRotJ = -1;

            for (int ri = 0; ri < ni; ri++) { //loop through rot I
                if (rotamerOptimization.check(tripleResidue1, ri)) {
                    continue;
                }
                for (int rj = 0; rj < nj; rj++) { //loop through rot J
                    if (rotamerOptimization.checkPrunedSingles(tripleResidue2, rj) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, tripleResidue2, rj)) {
                        continue;
                    }
                    double currentEnergy = 0.0;
                    for (int k = 0; k < nRes; k++) { //loop through all other residues
                        if (k == tripleResidue1 || k == tripleResidue2) {
                            continue;
                        }
                        Residue resK = residueList.get(k);
                        Rotamer rotK[] = resK.getRotamers(rLib);
                        int nk = rotK.length;

                        int rkStart = 0;
                        while (rotamerOptimization.checkPrunedSingles(k, rkStart) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, k, rkStart) || rotamerOptimization.checkPrunedPairs(tripleResidue2, rj, k, rkStart)) {
                            if (++rkStart >= nk) {
                                logger.warning("RJ is too large.");
                            }
                        }

                        double lowEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, ri, tripleResidue2, rj, k, rkStart);
                        for (int rk = rkStart; rk < nk; rk++) {
                            if (rotamerOptimization.checkPrunedSingles(k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue2, rj, k, rk)) {
                                continue;
                            } else {
                                double tripleEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, ri, tripleResidue2, rj, k, rk);
                                if (tripleEnergy < lowEnergy) {
                                    lowEnergy = tripleEnergy;
                                }
                            }
                        }
                        currentEnergy += lowEnergy; //adds lowest energy conformation of residue k to that of the rotamer I
                    }
                    if (currentEnergy < minEnergyIJ) {
                        minEnergyIJ = currentEnergy;
                        bestRotI = ri;
                        bestRotJ = rj;
                    }
                }
            }

            assertEquals(String.format(" %s Best three-body energy sum for residue %d is with rotamer %d at %10.4f.",
                    info, tripleResidue1, bestRotI, minEnergyIJ),
                    optimum[tripleResidue1], bestRotI);

            assertEquals(String.format(" %s Best three-body energy sum for residue %d is with rotamer %d at %10.4f.",
                    info, tripleResidue2, bestRotJ, minEnergyIJ),
                    optimum[tripleResidue2], bestRotJ);

            //loop over the residues to find the best rotamer per residue given bestRotI and bestRotJ
            for (int k = 0; k < nRes; k++) {
                if (k == tripleResidue1 || k == tripleResidue2) {
                    continue;
                }
                Residue resK = residueList.get(k);
                Rotamer rotK[] = resK.getRotamers(rLib);
                int nk = rotK.length;

                int rotCounter = 0;
                while (rotamerOptimization.checkPrunedPairs(tripleResidue1, bestRotI, k, rotCounter) && rotamerOptimization.checkPrunedPairs(tripleResidue2, bestRotJ, k, rotCounter) && rotCounter < nk) {
                    rotCounter++;
                }
                double lowEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k, rotCounter);
                int bestRotK = rotCounter;
                for (int rk = 1; rk < nk; rk++) {
                    if (rotamerOptimization.checkPrunedSingles(k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue1, bestRotI, k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue2, bestRotJ, k, rk)) {
                        continue;
                    } else {
                        double tripleEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k, rk);
                        if (tripleEnergy < lowEnergy) {
                            lowEnergy = tripleEnergy;
                            bestRotK = rk;
                        }
                    }
                }
                assertEquals(String.format(" %s Triple-Energy of residue (%d,%d) and residue (%d,%d) with residue %d",
                        info, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k), optimum[k], bestRotK);
            }
        }

    }

    @Test
    public void testPairEnergyElimination() {
        // Load the test system.
        load();

        // Run the optimization.
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(useOriginalRotamers);

        int counter = 1;
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = molecularAssembly.getChains();
        int nPolymers = polymers.length;
        for (int p = 0; p < nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            for (int i = 0; i < endResID; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers(rLib);
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    }
                    if (counter >= startResID) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }
        Residue[] residues = residueList.toArray(new Residue[residueList.size()]);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(useThreeBody);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setPruning(pruningLevel);
        rotamerOptimization.setEnergyRestartFile(restartFile);
        rotamerOptimization.setResidues(residueList);

        double energy;
        int nRes = residueList.size();
        if (doOverallOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            //System.out.println("The expected overall energy is: " + energy);
            assertEquals(info + " Total Energy", expectedEnergy, energy, tolerance);
        }

        // ToDo: Test self-energy use for rotamer 2-body eliminations.
        if (doSelfOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestSelfEnergyEliminations(true);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            //System.out.println("The expected self energy is: " + energy);
            assertEquals(info + " Self-Energy", expectedSelfEnergy, energy, tolerance);

            // Check that optimized rotamers are equivalent to the lowest self-energy of each residue.
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            // Loop over all residues
            for (int i = 0; i < nRes; i++) {
                Residue res = residueList.get(i);
                Rotamer[] rotI = res.getRotamers(rLib);
                int nRot = rotI.length;

                int rotCounter = 0;
                while (rotCounter < nRot && rotamerOptimization.checkPrunedSingles(i, rotCounter)) {
                    rotCounter++;
                }

                double lowEnergy = rotamerOptimization.getSelf(i, rotCounter);
                int bestRot = rotCounter;
                for (int ri = 1; ri < nRot; ri++) {
                    if (rotamerOptimization.checkPrunedSingles(i, ri)) {
                        continue;
                    } else {
                        double selfEnergy = rotamerOptimization.getSelf(i, ri);
                        if (selfEnergy < lowEnergy) {
                            lowEnergy = selfEnergy;
                            bestRot = ri;
                        }
                    }
                }
                assertEquals(String.format(" %s Self-Energy of residue %d", info, i), optimum[i], bestRot);
            }
        }

        // ToDo: Test 2-body energy use for rotamer pair eliminations.
        if (doPairOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestPairEnergyEliminations(pairResidue);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            //System.out.println("The expected 2-body energy is: " + energy);
            assertEquals(info + " Pair-Energy", expectedPairEnergy, energy, tolerance);

            // Check that optimized rotamers are equivalent to the lowest 2-body energy sum for the "pairResidue".
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            Residue resI = residueList.get(pairResidue);
            Rotamer rotI[] = resI.getRotamers(rLib);
            int ni = rotI.length;

            double minEnergy = Double.POSITIVE_INFINITY;
            int bestRotI = -1;

            // Loop over the pairResidue rotamers to find its lowest energy rotamer.
            for (int ri = 0; ri < ni; ri++) {
                double energyForRi = 0.0;
                if (rotamerOptimization.checkPrunedSingles(pairResidue, ri)) {
                    continue;
                }
                // Loop over residue J
                for (int j = 0; j < nRes; j++) {
                    if (j == pairResidue) {
                        continue;
                    }
                    Residue resJ = residueList.get(j);
                    Rotamer[] rotJ = resJ.getRotamers(rLib);
                    int nRot = rotJ.length;

                    int rj = 0;
                    while (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, ri, j, rj)) {
                        if (++rj >= nRot) {
                            logger.warning("RJ is too large.");
                        }
                    }

                    double lowEnergy = rotamerOptimization.get2Body(pairResidue, ri, j, rj);

                    for (rj = 1; rj < nRot; rj++) {
                        if (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, ri, j, rj)) {
                            continue;
                        } else {
                            double pairEnergy = rotamerOptimization.get2Body(pairResidue, ri, j, rj);
                            if (pairEnergy < lowEnergy) {
                                lowEnergy = pairEnergy;
                            }
                        }
                    }
                    energyForRi += lowEnergy;
                }
                if (energyForRi < minEnergy) {
                    minEnergy = energyForRi;
                    bestRotI = ri;
                }
            }

            assertEquals(String.format(" %s Best 2-body energy sum for residue %d is with rotamer %d at %10.4f.", info, pairResidue, bestRotI, minEnergy),
                    optimum[pairResidue], bestRotI);

            // Given the minimum energy rotamer for "pairResidue" is "bestRotI", we can check selected rotamers for all other residues.
            for (int j = 0; j < nRes; j++) {
                if (j == pairResidue) {
                    continue;
                }
                Residue resJ = residueList.get(j);
                Rotamer[] rotJ = resJ.getRotamers(rLib);
                int nRotJ = rotJ.length;

                int rotCounter = 0;
                while (rotamerOptimization.checkPrunedPairs(pairResidue, bestRotI, j, rotCounter) && rotCounter < nRotJ) {
                    rotCounter++;
                }

                double lowEnergy = rotamerOptimization.get2Body(pairResidue, bestRotI, j, rotCounter);
                int bestRotJ = rotCounter;
                for (int rj = 1; rj < nRotJ; rj++) {
                    if (rotamerOptimization.checkPrunedSingles(j, rj) || rotamerOptimization.checkPrunedPairs(pairResidue, bestRotI, j, rj)) {
                        continue;
                    } else {
                        double pairEnergy = rotamerOptimization.get2Body(pairResidue, bestRotI, j, rj);
                        if (pairEnergy < lowEnergy) {
                            lowEnergy = pairEnergy;
                            bestRotJ = rj;
                        }
                    }
                }
                if (bestRotJ != optimum[j]) {
                    // Check if 2-body energies are equal.
                    if (lowEnergy == rotamerOptimization.get2Body(pairResidue, bestRotI, j, optimum[j])) {
                        logger.warning(String.format(" Identical 2-body energies for %s: resi %d-%d, resj %d, best rotamer J %d, optimum J %d, 2-body energy (both) %10.6f", info, pairResidue, bestRotI, j, bestRotJ, optimum[j], lowEnergy));
                    } else {
                        assertEquals(String.format(" %s Pair-Energy of residue (%d,%d) with residue %d", info, pairResidue, bestRotI, j), optimum[j], bestRotJ);
                    }
                }
            }
        }

        // ToDo: Test 3-Body use for rotamer pair eliminations.
        if (doTripleOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestTripleEnergyEliminations(tripleResidue1, tripleResidue2);
            try {
                energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
                assertEquals(info + " Triple-Energy", expectedTripleEnergy, energy, tolerance);
            } catch (Exception e) {
                e.fillInStackTrace();
                e.printStackTrace();
                logger.log(java.util.logging.Level.INFO, "Error in doTripleOpt", e);
            }

            //Check that optimized rotamers are equivalent to the lowest 3-body energy of each residue with the tripleResidue1 and 2.
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            //fix residue 1 and gets its rotamers
            Residue resI = residueList.get(tripleResidue1);
            Rotamer rotI[] = resI.getRotamers(rLib);
            int ni = rotI.length;

            //fix residue 2 and get its rotamers
            Residue resJ = residueList.get(tripleResidue2);
            Rotamer rotJ[] = resJ.getRotamers(rLib);
            int nj = rotJ.length;

            double minEnergyIJ = Double.POSITIVE_INFINITY;
            int bestRotI = -1;
            int bestRotJ = -1;

            for (int ri = 0; ri < ni; ri++) { //loop through rot I
                if (rotamerOptimization.check(tripleResidue1, ri)) {
                    continue;
                }
                for (int rj = 0; rj < nj; rj++) { //loop through rot J
                    if (rotamerOptimization.checkPrunedSingles(tripleResidue2, rj) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, tripleResidue2, rj)) {
                        continue;
                    }
                    double currentEnergy = 0.0;
                    for (int k = 0; k < nRes; k++) { //loop through all other residues
                        if (k == tripleResidue1 || k == tripleResidue2) {
                            continue;
                        }
                        Residue resK = residueList.get(k);
                        Rotamer rotK[] = resK.getRotamers(rLib);
                        int nk = rotK.length;

                        int rkStart = 0;
                        while (rotamerOptimization.checkPrunedSingles(k, rkStart) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, k, rkStart) || rotamerOptimization.checkPrunedPairs(tripleResidue2, rj, k, rkStart)) {
                            if (++rkStart >= nk) {
                                logger.warning("RJ is too large.");
                            }
                        }

                        double lowEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, ri, tripleResidue2, rj, k, rkStart);
                        for (int rk = rkStart; rk < nk; rk++) {
                            if (rotamerOptimization.checkPrunedSingles(k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue1, ri, k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue2, rj, k, rk)) {
                                continue;
                            } else {
                                double tripleEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, ri, tripleResidue2, rj, k, rk);
                                if (tripleEnergy < lowEnergy) {
                                    lowEnergy = tripleEnergy;
                                }
                            }
                        }
                        currentEnergy += lowEnergy; //adds lowest energy conformation of residue k to that of the rotamer I
                    }
                    if (currentEnergy < minEnergyIJ) {
                        minEnergyIJ = currentEnergy;
                        bestRotI = ri;
                        bestRotJ = rj;
                    }
                }
            }

            assertEquals(String.format(" %s Best three-body energy sum for residue %d is with rotamer %d at %10.4f.",
                    info, tripleResidue1, bestRotI, minEnergyIJ),
                    optimum[tripleResidue1], bestRotI);

            assertEquals(String.format(" %s Best three-body energy sum for residue %d is with rotamer %d at %10.4f.",
                    info, tripleResidue2, bestRotJ, minEnergyIJ),
                    optimum[tripleResidue2], bestRotJ);

            //loop over the residues to find the best rotamer per residue given bestRotI and bestRotJ
            for (int k = 0; k < nRes; k++) {
                if (k == tripleResidue1 || k == tripleResidue2) {
                    continue;
                }
                Residue resK = residueList.get(k);
                Rotamer rotK[] = resK.getRotamers(rLib);
                int nk = rotK.length;

                int rotCounter = 0;
                while (rotamerOptimization.checkPrunedPairs(tripleResidue1, bestRotI, k, rotCounter) && rotamerOptimization.checkPrunedPairs(tripleResidue2, bestRotJ, k, rotCounter) && rotCounter < nk) {
                    rotCounter++;
                }
                double lowEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k, rotCounter);
                int bestRotK = rotCounter;
                for (int rk = 1; rk < nk; rk++) {
                    if (rotamerOptimization.checkPrunedSingles(k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue1, bestRotI, k, rk) || rotamerOptimization.checkPrunedPairs(tripleResidue2, bestRotJ, k, rk)) {
                        continue;
                    } else {
                        double tripleEnergy = rotamerOptimization.get3Body(residues, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k, rk);
                        if (tripleEnergy < lowEnergy) {
                            lowEnergy = tripleEnergy;
                            bestRotK = rk;
                        }
                    }
                }
                assertEquals(String.format(" %s Triple-Energy of residue (%d,%d) and residue (%d,%d) with residue %d",
                        info, tripleResidue1, bestRotI, tripleResidue2, bestRotJ, k), optimum[k], bestRotK);
            }
        }
    }

}
