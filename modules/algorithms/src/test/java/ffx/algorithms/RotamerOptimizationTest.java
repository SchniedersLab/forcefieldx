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

@RunWith(Parameterized.class)
public class RotamerOptimizationTest {

    private static final Logger logger = Logger.getLogger(RotamerOptimizationTest.class.getName());

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                    {
                        "Chignolin Direct with Orig Rot (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun0.restart",
                        0,                      // Pruning Level.
                        true,                   // Goldstein Elimination.
                        false,                  // Use 3-body Energies.
                        true,                   // Use Original Rotamers.
                        true,                   // Do Overall Opt.
                        -208.7299423233922,     // Expected Energy.
                        true,                   // Do Self-Energy Opt.
                        1.0E100,                // Expected Self-Energy.
                        false,                  // Do Pair-Energy Opt.
                        5,                      // Pair residue
                        0.0,                    // Expected Pair-Energy.
                        false,                  // Do Trimer-Energy Opt.
                        5,                      // Trimer residue 1.
                        6,                      // Trimer residue 2.
                        0.0,                    // Expected trimer energy.
                        1.0e-3                  // Energy Tolerance.
                    },
                    {
                        "Chignolin Direct with Orig Rot (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.restart",
                        1,                      // Pruning Level.
                        true,                   // Goldstein Elimination.
                        false,                  // Use 3-body Energies.
                        true,                   // Use Original Rotamers.
                        true,                   // Do Overall Opt.
                        -208.72994232339215,    // Expected Energy.
                        true,                   // Do Self-Energy Opt.
                        1.0E100,                // Expected Self-Energy.
                        false,                  // Do Pair-Energy Opt.
                        5,                      // Pair residue
                        0.0,                    // Expected Pair-Energy.
                        false,                  // Do Trimer-Energy Opt.
                        5,                      // Trimer residue 1.
                        6,                      // Trimer residue 2.
                        0.0,                    // Expected trimer energy.
                        1.0e-3                  // Energy Tolerance.
                    },
                    {
                        "Chignolin Direct with Orig Rot (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun2.restart",
                        2,                      // Pruning Level.
                        true,                   // Goldstein Elimination.
                        false,                  // Use 3-body Energies.
                        true,                   // Use Original Rotamers.
                        true,                   // Do Overall Opt.
                        -208.72994232339218,    // Expected Energy.
                        true,                   // Do Self-Energy Opt.
                        -201.19341910515593,    // Expected Self-Energy.
                        false,                   // Do Pair-Energy Opt.
                        5,                      // Pair residue
                        0.0,                    // Expected Pair-Energy.
                        false,                  // Do Trimer-Energy Opt.
                        5,                      // Trimer residue 1.
                        6,                      // Trimer residue 2.
                        0.0,                    // Expected trimer energy.
                        1.0e-3                  // Energy Tolerance.
                    },
                   /* {
                        "Chignolin Direct with Orig Rot (Goldstein)",
                        "ffx/algorithms/structures/5awl.pdb",
                        "ffx/algorithms/structures/5awl.direct.orig.prun1.3body.restart",
                        2,                      // Pruning Level.
                        true,                   // Goldstein Elimination.
                        true,                   // Use 3-body Energies.
                        true,                   // Use Original Rotamers.
                        true,                   // Do Overall Opt.
                        -208.72994232339215,    // Expected Energy.
                        true,                   // Do Self-Energy Opt.
                        1.0E100,                // Expected Self-Energy.
                        false,                  // Do Pair-Energy Opt.
                        5,                      // Pair residue
                        0.0,                    // Expected Pair-Energy.
                        false,                   // Do Trimer-Energy Opt.
                        5,                      // Trimer residue 1.
                        6,                      // Trimer residue 2.
                        0.0,                    // Expected trimer energy.
                        1.0e-3                  // Energy Tolerance.
                    }*/
               
 
                    
                /** {
                        "Trpcage Direct with Orig Rot (Goldstein)",
                        "ffx/algorithms/structures/trpcage.pdb",
                        "ffx/algorithms/structures/trpcage.direct.orig.restart",
                        0,              // Pruning Level.
                        true,           // Goldstein Elimination.
                        false,          // Use 3-body Energies.
                        true,           // Use Original Rotamers.
                        false,          // Do Overall Opt.
                        -628.143879,    // Expected Energy.
                        true,           // Do Self-Energy Opt.
                        22086.367273,   // Expected Self-Energy.
                        false,          // Do Pair-Energy Opt.
                        5,              // Pair residue
                        1280.865248,    // Expected Pair-Energy.
                        false,          // Do Trimer-Energy Opt.
                        5,              // Trimer residue 1.
                        6,              // Trimer residue 2.
                        0.0,            // Expected trimer energy.
                        1.0e-3          // Energy Tolerance.
                }, 
        {
        "Trpcage Direct with Orig Rot (DEE, Prune=1)",
        "ffx/algorithms/structures/trpcage.pdb",
        "ffx/algorithms/structures/trpcage.direct.orig.restart",
        1,              // Pruning Level.
        false,          // Goldstein Elimination.
        false,          // Use 3-body Energies.
        true,           // Use Original Rotamers.
        false,          // Do Overall Opt. (false because DEE leaves too many permutations)
        -628.143879,    // Expected Energy.
        true,           // Do Self-Energy Opt.
        22086.367273,   // Expected Self-Energy.
        true,           // Do Pair-Energy Opt.
        5,              // Pair residue
        1280.865248,    // Expected Pair-Energy.
        false,          // Do Trimer-Energy Opt.
        5,              // Trimer residue 1.
        6,              // Trimer residue 2.
        0.0,            // Expected trimer energy.
        1.0e-3          // Energy Tolerance.
        },
        {
        "Trpcage Direct with Orig Rot (DEE, Prune=2)",
        "ffx/algorithms/structures/trpcage.pdb",
        "ffx/algorithms/structures/trpcage.direct.orig.restart",
        2,              // Pruning Level.
        false,          // Goldstein Elimination.
        false,          // Use 3-body Energies.
        true,           // Use Original Rotamers.
        true,           // Do Overall Opt. (false because DEE leaves too many permutations)
        -628.143879,    // Expected Energy.
        true,           // Do Self-Energy Opt.
        -616.237522,    // Expected Self-Energy.
        true,           // Do Pair-Energy Opt.
        5,              // Pair residue
        -563.929018,    // Expected Pair-Energy.
        false,          // Do Trimer-Energy Opt.
        5,              // Trimer residue 1.
        6,              // Trimer residue 2.
        0.0,            // Expected trimer energy.
        1.0e-3          // Energy Tolerance.
        },
        {
        "Trpcage Direct (Goldstein)",
        "ffx/algorithms/structures/trpcage.pdb",
        "ffx/algorithms/structures/trpcage.direct.3body.restart",
        1,              // Pruning Level.
        true,           // Goldstein Elimination.
        true,           // Use 3-body Energies.
        true,           // Use Original Rotamers.
        true,           // Do Overall Opt.
        -628.143879,    // Expected Energy.
        true,           // Do Self-Energy Opt.
        22086.367273,   // Expected Self-Energy.
        true,           // Do Pair-Energy Opt.
        5,              // Pair residue.
        1280.865247,    // Expected Pair-Energy.
        true,           // Do Trimer-Energy Opt.
        5,              // Trimer residue 1.
        6,              // Trimer residue 2.
        1276.608834,    // Expected trimer energy.
        1.0e-3          // Energy Tolerance.
        } */
        });
    }

    String info;
    String filename;
    String restartName;
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

    @Test
    public void testSelfEnergyElimination() {
        // Load the test system.
        load();

        // Initialize Parallel Java
        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
            message = String.format(" Skipping rotamer optimization test.");
            logger.log(Level.WARNING, message, e.toString());
            return;
        }

        // Run the optimization.
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(useOriginalRotamers);

        int counter = 1;
        int allStartResID = 1;
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = molecularAssembly.getChains();
        int nPolymers = polymers.length;
        for (int p = 0; p < nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            int nResidues = residues.size();
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers(rLib);
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    }
                    if (counter >= allStartResID) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }

        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(useThreeBody);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setPruning(pruningLevel);
        rotamerOptimization.setEnergyRestartFile(restartFile);
        rotamerOptimization.setResidues(residueList);

        double energy = 0.0;
        if (doOverallOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected overall energy is: " + energy);
            assertEquals(info + " Total Energy", expectedEnergy, energy, tolerance);
        }

        int nRes = residueList.size();
        if (doSelfOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestSelfEnergyEliminations(true);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected self energy is: " + energy);
            assertEquals(info + " Self-Energy", expectedSelfEnergy, energy, tolerance);

            // Check that optimized rotamers are equivalent to the lowest self-energy of each residue.
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            // Loop over all residues
            for (int i = 0; i < nRes; i++) {
                Residue res = residueList.get(i);
                Rotamer[] rotI = res.getRotamers(rLib);
                int nRot = rotI.length;
                double lowEnergy = rotamerOptimization.self(i, 0);
                int bestRot = 0;
                for (int ri = 1; ri < nRot; ri++) {
                    double selfEnergy = rotamerOptimization.self(i, ri);
                    if (selfEnergy < lowEnergy) {
                        lowEnergy = selfEnergy;
                        bestRot = ri;
                    }
                }
                //assertEquals(String.format(" %s Self-Energy of residue %d", info, i), optimum[i], bestRot);
            }
        }

        if (doPairOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestPairEnergyEliminations(pairResidue);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected pair energy is: " + energy);

            // Check that optimized rotamers are equivalent to the lowest pair-energy sum for the "pairResidue".
            int optimum[] = rotamerOptimization.getOptimumRotamers();

            Residue resI = residueList.get(pairResidue);
            Rotamer rotI[] = resI.getRotamers(rLib);
            int ni = rotI.length;

            double minEnergy = Double.POSITIVE_INFINITY;
            int bestRotI = -1;

            // Loop over the pairResidue rotamers to find its lowest energy rotamer.
            for (int ri = 0; ri < ni; ri++) {
                double energyForRi = 0.0;
                // Loop over residue J
                for (int j = 0; j < nRes; j++) {
                    if (j == pairResidue) {
                        continue;
                    }
                    Residue resJ = residueList.get(j);
                    Rotamer[] rotJ = resJ.getRotamers(rLib);
                    int nRot = rotJ.length;
                    double lowEnergy = rotamerOptimization.pair(pairResidue, ri, j, 0);
                    for (int rj = 1; rj < nRot; rj++) {
                        double pairEnergy = rotamerOptimization.pair(pairResidue, ri, j, rj);
                        if (pairEnergy < lowEnergy) {
                            lowEnergy = pairEnergy;
                        }
                    }
                    energyForRi += lowEnergy;
                }
                if (energyForRi < minEnergy) {
                    minEnergy = energyForRi;
                    bestRotI = ri;
                }
            }

            assertEquals(String.format(" %s Best pair-energy sum for residue %d is with rotamer %d at %10.4f.", info, pairResidue, bestRotI, minEnergy),
                    optimum[pairResidue], bestRotI);

            // Given the minimum energy rotamer for "pairResidue" is "bestRotI", we can check selected rotamers for all other residues.
            for (int j = 0; j < nRes; j++) {
                if (j == pairResidue) {
                    continue;
                }
                Residue resJ = residueList.get(j);
                Rotamer[] rotJ = resJ.getRotamers(rLib);
                int nRotJ = rotJ.length;
                double lowEnergy = rotamerOptimization.pair(pairResidue, bestRotI, j, 0);
                int bestRotJ = 0;
                for (int rj = 1; rj < nRotJ; rj++) {
                    double pairEnergy = rotamerOptimization.pair(pairResidue, bestRotI, j, rj);
                    if (pairEnergy < lowEnergy) {
                        lowEnergy = pairEnergy;
                        bestRotJ = rj;
                    }
                }
                assertEquals(String.format(" %s Pair-Energy of residue (%d,%d) with residue %d", info, pairResidue, bestRotI, j), optimum[j], bestRotJ);
            }
        }

        // ToDo: Test 3-Body Energy Eliminations.
        if (doTripleOpt) {
            rotamerOptimization.turnRotamerPairEliminationOff();
            rotamerOptimization.setTestTripleEnergyEliminations(tripleResidue1, tripleResidue2);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected pair energy is: " + energy);
            assertEquals(info + " Triple-Energy", expectedTripleEnergy, energy, tolerance);

            // ToDo: Check that optimized rotamers are equivalent to the lowest 3-body energy of each residue with the tripleResidue1 and 2.
        }
    } 
    
    
    
    
    
 
    /*
    @Test
    public void testPairEnergyElimination() {
        // Load the test system.
        load();

        // Initialize Parallel Java
        try {
            String args[] = new String[0];
            Comm.init(args);
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
            message = String.format(" Skipping rotamer optimization test.");
            logger.log(Level.WARNING, message, e.toString());
            return;
        }

        // Run the optimization.
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(useOriginalRotamers);

        int counter = 1;
        int allStartResID = 1;
        ArrayList<Residue> residueList = new ArrayList<Residue>();
        Polymer[] polymers = molecularAssembly.getChains();
        int nPolymers = polymers.length;
        for (int p = 0; p < nPolymers; p++) {
            Polymer polymer = polymers[p];
            ArrayList<Residue> residues = polymer.getResidues();
            int nResidues = residues.size();
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues.get(i);
                Rotamer[] rotamers = residue.getRotamers(rLib);
                if (rotamers != null) {
                    int nrot = rotamers.length;
                    if (nrot == 1) {
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    }
                    if (counter >= allStartResID) {
                        residueList.add(residue);
                    }
                }
                counter++;
            }
        }

        RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
        rotamerOptimization.setThreeBodyEnergy(useThreeBody);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setPruning(pruningLevel);
        rotamerOptimization.setEnergyRestartFile(restartFile);
        rotamerOptimization.setResidues(residueList);

        System.out.println("BEGINNING PAIR ELIMINATION TEST");
        double energy = 0.0;
        if (doOverallOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected overall energy is: " + energy);
            //assertEquals(info + " Total Energy", expectedEnergy, energy, tolerance);
        }


        // ToDo: Test self-energy use for rotamer pair eliminations.
        if (doSelfOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestSelfEnergyEliminations(true);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected self energy is: " + energy);
            //assertEquals(info + " Self-Energy", expectedSelfEnergy, energy, tolerance);
        } 

        // ToDo: Test pair-energy use for rotamer pair eliminations.
        if (doPairOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestPairEnergyEliminations(pairResidue);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected pair energy is: " + energy);
            //assertEquals(info + " Pair-Energy", expectedPairEnergy, energy, tolerance);
        }
        
        // ToDo: Test 3-Body use for rotamer pair eliminations.
       if (doTripleOpt) {
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.setTestTripleEnergyEliminations(tripleResidue1, tripleResidue2);
            energy = rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
            System.out.println("The expected trimer energy is: " + energy);
            //assertEquals(info + " Triple-Energy", expectedTripleEnergy, energy, tolerance);
        }
    } */

}
