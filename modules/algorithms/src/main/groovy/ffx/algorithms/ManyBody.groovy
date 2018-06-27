package ffx.algorithms

import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Residue.ResidueType
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.MolecularAssembly


import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Run ManyBody algorithm on a system.", name = "ffxc ManyBody")
class ManyBody extends AlgorithmsScript {

    @Mixin
    ManyBodyOptions manyBody;

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "XYZ or PDB input files.")
    private List<String> filenames
    def run() {
        
        if (!init()) {
            return
        }
        
        String filename = null;
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            active = assemblies[0]
            filename = filenames.get(0)
        } else if (active == null) {
            logger.info(helpString())
            return
        } else {
            filename = active.getFile().getAbsolutePath();
        }
    
        /**
         * Algorithm options.
         */
        boolean useGoldstein = !manyBody.dee
        boolean revert = !manyBody.noRevert
        boolean useOrigCoordsRotamer = !manyBody.original
        boolean verboseEnergies = manyBody.verbose
        boolean decomposeOriginal = manyBody.decompose
        if (decomposeOriginal) {
            useOrigCoordsRotamer = true;
        }

        // By default, rotamer optimization should silence GK warnings, because occasionally we will have unreasonable configurations.
        if (System.getProperty("gk-suppressWarnings") == null) {
            System.setProperty("gk-suppressWarnings", "true");
        }

        /**
         * Chain, Residue and/or Box selections.
         */
        int allStartResID = manyBody.all
        // Internal machinery indexed 0 to (n-1)
        int boxStart = manyBody.start - 1
        int boxEnd = manyBody.finish - 1

        if (manyBody.algorithm != 5) {
            // Not Box optimization.
            if (allStartResID < 1 && manyBody.listResidues.equalsIgnoreCase("none")) {
                if (manyBody.finish < manyBody.start || manyBody.start < 0 || manyBody.finish < 0) {
                    logger.warning(" FFX shutting down: no residues specified for optimization.")
                    return;
                }
            }
        } else {
            // Box optimization.
            if (allStartResID > 0) {
                boxStart = allStartResID - 1 // Internal machinery indexed 0 to (n-1)
                if (boxStart < 0) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).")
                    return;
                }
            } else {
                if (boxStart < 0 || (boxEnd > -1 && boxEnd < boxStart)) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
                    return;
                }
            }
        }


       /**
         * Box optimization options.
         */
        int[] numXYZBoxes = new int[3];
        if (manyBody.algorithm == 5) {
            String input = manyBody.numBoxes
            Scanner boxNumInput = new Scanner(input);
            boxNumInput.useDelimiter(",");
            int inputLoopCounter = 0;
            numXYZBoxes[0] = 3; // Default
            while (inputLoopCounter < numXYZBoxes.length) {
                if (boxNumInput.hasNextInt()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt();
                    inputLoopCounter++;
                } else if (boxNumInput.hasNextDouble()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextDouble();
                    inputLoopCounter++;
                    logger.info("Double input to nB truncated to integer.");
                } else if (boxNumInput.hasNext()) {
                    logger.info("Non-numeric input to nB discarded");
                    boxNumInput.next();
                } else {
                    logger.info("Insufficient input to nB. Non-input values assumed either equal to X or default to 3");
                    break;
                }
            }
            boxNumInput.close();
            for (int i = inputLoopCounter; i < numXYZBoxes.length; i++) {
                numXYZBoxes[i] = numXYZBoxes[0];
            }
            for (int i = 0; i < numXYZBoxes.length; i++) {
                if (numXYZBoxes[i] == 0) {
                    numXYZBoxes[i] = 3;
                    logger.info("Input of zero to nB reset to default of three.");
                } else if (numXYZBoxes[i] < 0) {
                    numXYZBoxes[i] = -1 * numXYZBoxes[i];
                    logger.info("Input of negative number to nB reset to positive number");
                }
            }
        }

        double boxBorderSize = manyBody.boxBorderSize
        double approxBoxLength = manyBody.approxBoxLength
        if (approxBoxLength < 0) {
            logger.info(" Negative box length value changed to -1 * input.")
            approxBoxLength *= -1;
        }
        int boxInclusionCriterion = manyBody.boxInclusionCriterion

        /**
         * Minimize the final structure.
         */
        boolean min = false;
        double eps = manyBody.minimize
        if (eps > 0.0) {
            min = true
        }

        /**
         * Monte Carlo search of remaining permutations.
         */
        boolean monteCarlo = false
        int nMCSteps = manyBody.monteCarlo;
        if (nMCSteps > 1) {
            monteCarlo = true
        }

        /**
         * Nucleic acid options.
         */
        double nucleicCorrectionThreshold = manyBody.nucleicCorrectionThreshold
        int minimumNumberAcceptedNARotamers = manyBody.minimumAcceptedNARotamers

        /**
         * Energy restart.
         */
        boolean useEnergyRestart = false;
        File energyRestartFile = null;
        if (!manyBody.energyRestart.equalsIgnoreCase("none")) {
            useEnergyRestart = true;
            energyRestartFile = new File(manyBody.energyRestart)
        }

        /**
         * Force residues.
         */
        String forceResidues = manyBody.forceResidues
        int forceResiduesStart = -1
        int forceResiduesEnd = -1;

        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();

        List<String> resList = new ArrayList<>();
        if (!manyBody.listResidues.equalsIgnoreCase("none")) {
            def tok = (manyBody.listResidues).tokenize('.');
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }

        List<String> sequenceOptimizationList = new ArrayList<>();
        if (!manyBody.sequence.equalsIgnoreCase("none")) {
            def tok = manyBody.sequence.tokenize('.');
            for (String t : tok) {
                logger.info(" Sequence optimizing " + t);
                sequenceOptimizationList.add(t);
            }
            if (System.getProperty("relative-solvation") == null) {
                System.setProperty("relative-solvation", "AUTO");
            }
        }

        List<String> titrationOptimizationList = new ArrayList<>();
        if (!manyBody.titrationOptimization.equalsIgnoreCase("none")) {
            def tok = manyBody.titrationOptimization.tokenize('.');
            for (String t : tok) {
                logger.info(" Optimizing protonation states of " + t);
                titrationOptimizationList.add(t);
            }
        }

        /**
         * Evaluate forced residues for the sliding window algorithm
         */
        if (manyBody.algorithm == 4 && !forceResidues.equalsIgnoreCase('-1,-1')) {
            String input = forceResidues;
            Scanner frScan = new Scanner(input);
            frScan.useDelimiter(",");
            try {
                if (!frScan.hasNextInt()) {
                    frScan.next(); // Discards extra input to indicate a negative value of frStart.
                }
                forceResiduesStart = frScan.nextInt();
                forceResiduesEnd = frScan.nextInt();
            } catch (Exception ex) {
                logger.severe(String.format(" FFX shutting down: input to -fR could not be parsed as a pair of integers: %s", input));
            }
            if (forceResiduesStart > forceResiduesEnd) {
                logger.info(" Start of range higher than ending: start flipped with end.");
                int temp = forceResiduesStart;
                forceResiduesStart = forceResiduesEnd;
                forceResiduesEnd = temp;
            }
            if (forceResiduesEnd < 1) {
                logger.severe(String.format(" FFX shutting down: end range for -fR must be at least 1; input range %d to %d",
                        forceResiduesStart, forceResiduesEnd));
            }
        }

        if (manyBody.algorithm != 5) {
            if (!manyBody.listResidues.equalsIgnoreCase("none")) {
                String info = "\n Evaluating rotamers for residues ";
                for (String i : resList) {
                    info += String.format("%s, ", i);
                }
                logger.info(info);
            } else if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for residues " + manyBody.start + " to " + manyBody.finish);
            } else {
                logger.info("\n Evaluating rotamers for all residues beginning at " + allStartResID);
            }
        } else {
            if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1));
            } else {
                logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1));
            }
        }

        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(
                activeAssembly, activeAssembly.getPotentialEnergy(), algorithmListener);

        rotamerOptimization.setTwoBodyCutoff(manyBody.twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(manyBody.threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(manyBody.threeBody);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setRevert(revert);
        rotamerOptimization.setPruning(manyBody.prune);
        rotamerOptimization.setWindowSize(manyBody.window);
        rotamerOptimization.setIncrement(manyBody.increment);
        rotamerOptimization.setDistanceCutoff(manyBody.cutoff);
        rotamerOptimization.setNucleicCorrectionThreshold(manyBody.nucleicCorrectionThreshold);
        rotamerOptimization.setMinimumNumberAcceptedNARotamers(manyBody.minimumAcceptedNARotamers);
        rotamerOptimization.setVerboseEnergies(verboseEnergies);
        rotamerOptimization.setBoxBorderSize(manyBody.boxBorderSize);
        rotamerOptimization.setApproxBoxLength(manyBody.approxBoxLength);
        rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
        rotamerOptimization.setBoxInclusionCriterion(manyBody.boxInclusionCriterion);
        rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
        rotamerOptimization.setMonteCarlo(monteCarlo, nMCSteps);

        if (useEnergyRestart) {
            rotamerOptimization.setEnergyRestartFile(energyRestartFile);
        }

        if (manyBody.library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        }

        if (useOrigCoordsRotamer) {
            rLib.setUseOrigCoordsRotamer(true);
        }

        int counter = 1;
        if (manyBody.algorithm != 5) {
            if (allStartResID > 0) {
                ArrayList<Residue> residueList = new ArrayList<Residue>();
                Polymer[] polymers = active.getChains();
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
                        } else if (!manyBody.forceResidues.equalsIgnoreCase("none")) {
                            if (counter >= allStartResID && counter >= forceResiduesStart
                                && counter <= forceResiduesEnd) {
                                residueList.add(residue);
                            }
                        }
                        counter++;
                    }
                }
                rotamerOptimization.setResidues(residueList);
            } else if (!manyBody.listResidues.equalsIgnoreCase("none")) {
                ArrayList<Residue> residueList = new ArrayList<>();
                Polymer[] polymers = active.getChains();
                int n = 0;
                for (String s : resList) {
                    Character chainID = s.charAt(0);
                    int i = Integer.parseInt(s.substring(1));
                    for (Polymer p : polymers) {
                        if (p.getChainID() == chainID) {
                            List<Residue> rs = p.getResidues();
                            for (Residue r : rs) {
                                if (r.getResidueNumber() == i) {
                                    residueList.add(r);
                                    Rotamer[] rotamers = r.getRotamers(rLib);
                                    if (rotamers != null) {
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }
                rotamerOptimization.setResiduesIgnoreNull(residueList);
                if (n < 1) {
                    return;
                }
            } else if (!manyBody.chain.equalsIgnoreCase('-1')) {
                rotamerOptimization.setResidues(manyBody.chain, manyBody.start, manyBody.finish);
            } else {
                rotamerOptimization.setResidues(manyBody.start, manyBody.finish);
            }
        } else {
            boolean ignoreNA = false;
            String ignoreNAProp = System.getProperty("ignoreNA");
            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
                ignoreNA = true;
            }
            ArrayList<Residue> residueList = new ArrayList<Residue>();
            Polymer[] polymers = active.getChains();
            int nPolymers = polymers.length;
            for (int p = 0; p < nPolymers; p++) {
                Polymer polymer = polymers[p];
                ArrayList<Residue> residues = polymer.getResidues();
                int nResidues = residues.size();
                for (int i = 0; i < nResidues; i++) {
                    Residue residue = residues.get(i);
                    if (ignoreNA && residue.getResidueType() == ResidueType.NA) {
                        continue;
                    }
                    Rotamer[] rotamers = residue.getRotamers(rLib);
                    if (rotamers != null) {
                        int nrot = rotamers.length;
                        if (nrot == 1) {
                            RotamerLibrary.applyRotamer(residue, rotamers[0]);
                        } else if (nrot > 1) {
                            residueList.add(residue);
                        }
                    }
                    counter++;
                }
            }
            rotamerOptimization.setResidues(residueList);
            rotamerOptimization.setBoxStart(boxStart);
            if (allStartResID != -1) {
                rotamerOptimization.setBoxEnd(boxEnd);
            }
        }

        ArrayList<Residue> residueList = rotamerOptimization.getResidues();

        if (!manyBody.sequence.equalsIgnoreCase("none")) {
            for (String s : sequenceOptimizationList) {
                Character chainID = s.charAt(0);
                int num = Integer.parseInt(s.substring(1));
                for (int i = 0; i < residueList.size(); i++) {
                    Residue res = residueList.get(i);
                    if (res.getChainID() == chainID && res.getResidueNumber() == num) {
                        MultiResidue multiRes = new MultiResidue(res, active.getForceField(), active.getPotentialEnergy());
                        for (Polymer polymer : active.getChains()) {
                            if (polymer.getChainID() == chainID) {
                                logger.info(String.format(" Adding multiresidue %s to chain %c.", multiRes, chainID));
                                polymer.addMultiResidue(multiRes);
                            }
                        }
                        for (CommonAminoAcid3 aa : CommonAminoAcid3.values()) {
                            if (aa.toString().equals("PRO") || aa.toString().equals("GLY")) {
                                continue;
                            }
                            if (!aa.toString().equalsIgnoreCase(res.getName())) {
                                logger.info(String.format(" Adding %s to residue %s.", aa.toString(), multiRes.toString()));
                                multiRes.addResidue(new Residue(aa.toString(), res.getResidueNumber(), ResidueType.AA));
                            }
                        }
                        //multiRes.requestSetActiveResidue(ResidueEnumerations.AminoAcid3.valueOf(res.getName()));
                        multiRes.setActiveResidue(res);
                        active.getPotentialEnergy().reInit();
                        residueList.remove(i);
                        residueList.add(i, multiRes);
                    }
                }
            }
        }

        if (!manyBody.titrationOptimization.equalsIgnoreCase("none")) {
            ArrayList<Residue> titrating = new ArrayList<>();
            for (String s : titrationOptimizationList) {
                Character chainID = s.charAt(0);
                int num = Integer.parseInt(s.substring(1));
                for (int i = 0; i < residueList.size(); i++) {
                    Residue res = residueList.get(i);
                    if (res.getChainID() == chainID && res.getResidueNumber() == num) {
                        titrating.add(res);
                    }
                }
            }
            rotamerOptimization.titrationSetResidues(titrating);
        }

        boolean master = true;
        if (Comm.world().size() > 1) {
            int rank = Comm.world().rank();
            if (rank != 0) {
                master = false;
            }
        }

        energy();
        RotamerLibrary.measureRotamers(residueList, false);

        if (decomposeOriginal) {
            rLib.setUseOrigCoordsRotamer(true);
            rotamerOptimization.setDecomposeOriginal(true);
        }

        if (manyBody.algorithm == 1) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT);
        } else if (manyBody.algorithm == 2) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
        } else if (manyBody.algorithm == 3) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BRUTE_FORCE);
        } else if (manyBody.algorithm == 4) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW);
        } else if (manyBody.algorithm == 5) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX);
        }

        if (master) {
            if (min) {
                minimize(eps);
            }
            logger.info(" Final Minimum Energy");
            energy();
            String ext = FilenameUtils.getExtension(filename);
            filename = FilenameUtils.removeExtension(filename);
            if (ext.toUpperCase().contains("XYZ")) {
                algorithmFunctions.saveAsXYZ(active, new File(filename + ".xyz"));
            } else {
                algorithmFunctions.saveAsPDB(active, new File(filename + ".pdb"));
            }
        }
        
        if (manyBody.saveOutput) {
            rotamerOptimization.outputEliminated();
        }

    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
