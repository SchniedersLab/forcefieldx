package realspace

import ffx.algorithms.AlgorithmFunctions
import ffx.potential.MolecularAssembly
import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.algorithms.RotamerOptimization
import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Residue.ResidueType
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.realspace.RealSpaceData
import ffx.realspace.parsers.RealSpaceFile
import ffx.xray.CrystalReciprocalSpace
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize
import ffx.xray.parsers.DiffractionFile

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.realspace.cli.RealSpaceOptions

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "ffxc realspace.ManyBody")
class ManyBody extends AlgorithmsScript {

    @Mixin
    RealSpaceOptions realSpaceOptions

    @Mixin
    ManyBodyOptions manyBodyOptions

    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames

    def run() {

        if (!init()) {
            return
        }

        realSpaceOptions.init()

        String modelfilename
        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
            assemblies = { activeAssembly }
        }
        /**
         * Algorithm options.
         */
        int library = manyBodyOptions.library
        int algorithm = manyBodyOptions.algorithm
        boolean useGoldstein = !manyBodyOptions.dee
        boolean threeBodyTerm = manyBodyOptions.threeBody
        int pruningType = manyBodyOptions.prune
        boolean revert = !manyBodyOptions.noRevert
        boolean useOrigCoordsRotamer = !manyBodyOptions.original
        boolean verboseEnergies = manyBodyOptions.verbose
        boolean decomposeOriginal = manyBodyOptions.decompose
        double twoBodyCutoff = manyBodyOptions.twoBodyCutoff
        double threeBodyCutoff = manyBodyOptions.threeBodyCutoff
        if (decomposeOriginal) {
            useOrigCoordsRotamer = true
        }

        // By default, rotamer optimization should silence GK warnings, because occasionally we will have unreasonable configurations.
        if (System.getProperty("gk-suppressWarnings") == null) {
            System.setProperty("gk-suppressWarnings", "true")
        }

        /**
         * Chain, Residue and/or Box selections.
         */
        String chain = manyBodyOptions.chain
        int startResID = manyBodyOptions.start
        int allStartResID = manyBodyOptions.all
        int finalResID = manyBodyOptions.finish
        // Internal machinery indexed 0 to (n-1)
        int boxStart = manyBodyOptions.start - 1
        int boxEnd = manyBodyOptions.finish - 1

        if (algorithm != 5) {
            // Not Box optimization.
            if (allStartResID < 1 && manyBodyOptions.listResidues.equalsIgnoreCase('none')) {
                if (finalResID < startResID || startResID < 0 || finalResID < 0) {
                    logger.warning(" FFX shutting down: no residues specified for optimization.")
                    return
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
         * Sliding window options.
         */
        int windowSize = manyBodyOptions.window
        int increment = manyBodyOptions.increment
        double distance = manyBodyOptions.cutoff

        /**
         * Box optimization options.
         */
        int[] numXYZBoxes = new int[3];
        if (algorithm == 5) {
            String input = manyBodyOptions.numBoxes
            Scanner boxNumInput = new Scanner(input)
            boxNumInput.useDelimiter(",")
            int inputLoopCounter = 0
            numXYZBoxes[0] = 3 // Default
            while (inputLoopCounter < numXYZBoxes.length) {
                if (boxNumInput.hasNextInt()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt()
                    inputLoopCounter++
                } else if (boxNumInput.hasNextDouble()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextDouble()
                    inputLoopCounter++
                    logger.info("Double input to nB truncated to integer.")
                } else if (boxNumInput.hasNext()) {
                    logger.info("Non-numeric input to nB discarded")
                    boxNumInput.next()
                } else {
                    logger.info("Insufficient input to nB. Non-input values assumed either equal to X or default to 3")
                    break
                }
            }
            boxNumInput.close()
            for (int i = inputLoopCounter; i < numXYZBoxes.length; i++) {
                numXYZBoxes[i] = numXYZBoxes[0]
            }
            for (int i = 0; i < numXYZBoxes.length; i++) {
                if (numXYZBoxes[i] == 0) {
                    numXYZBoxes[i] = 3
                    logger.info("Input of zero to nB reset to default of three.")
                } else if (numXYZBoxes[i] < 0) {
                    numXYZBoxes[i] = -1 * numXYZBoxes[i]
                    logger.info("Input of negative number to nB reset to positive number")
                }
            }
        }

        double boxBorderSize = manyBodyOptions.boxBorderSize
        double approxBoxLength = manyBodyOptions.approxBoxLength
        if (approxBoxLength < 0) {
            logger.info(" Negative box length value changed to -1 * input.")
            approxBoxLength *= -1
        }
        int boxInclusionCriterion = manyBodyOptions.boxInclusionCriterion

        /**
         * Minimize the final structure.
         */
        boolean min = false
        double eps = manyBodyOptions.minimize
        if (eps > 0.0) {
            min = true
        }

        /**
         * Monte Carlo search of remaining permutations.
         */
        boolean monteCarlo = false
        int nMCSteps = monteCarlo
        if (nMCSteps > 1) {
            monteCarlo = true
        }

        /**
         * Nucleic acid options.
         */
        double nucleicCorrectionThreshold = manyBodyOptions.nucleicCorrectionThreshold
        int minimumNumberAcceptedNARotamers = manyBodyOptions.minimumAcceptedNARotamers

        /**
         * Energy restart.
         */
        boolean useEnergyRestart = false
        File energyRestartFile = null
        if (!manyBodyOptions.energyRestart.equalsIgnoreCase('none')) {
            useEnergyRestart = true
            energyRestartFile = new File(manyBodyOptions.energyRestart)
        }

        /**
         * Force residues.
         */
        String forceResidues = manyBodyOptions.forceResidues
        int forceResiduesStart = -1
        int forceResiduesEnd = -1

        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();

        List<String> resList = new ArrayList<>()
        if (!manyBodyOptions.listResidues.equalsIgnoreCase('none')) {
            def tok = (manyBodyOptions.listResidues).tokenize('.')
            for (String t : tok) {
                logger.info(" Adding " + t)
                resList.add(t)
            }
        }

        List<String> sequenceOptimizationList = new ArrayList<>()
        if (!manyBodyOptions.sequence.equalsIgnoreCase('none')) {
            def tok = manyBodyOptions.sequence.tokenize('.')
            for (String t : tok) {
                logger.info(" Sequence optimizing " + t)
                sequenceOptimizationList.add(t)
            }
            if (System.getProperty("relative-solvation") == null) {
                System.setProperty("relative-solvation", "AUTO")
            }
        }

        List<String> titrationOptimizationList = new ArrayList<>()
        if (!manyBodyOptions.titrationOptimization.equalsIgnoreCase('none')) {
            def tok = manyBodyOptions.titrationOptimization.tokenize('.')
            for (String t : tok) {
                logger.info(" Protonation state optimizing " + t)
                titrationOptimizationList.add(t)
            }
        }

        /**
         * Evaluate forced residues for the sliding window algorithm
         */
        if (algorithm == 4 && !forceResidues.equalsIgnoreCase('-1,-1')) {
            String input = forceResidues
            Scanner frScan = new Scanner(input)
            frScan.useDelimiter(",")
            try {
                if (!frScan.hasNextInt()) {
                    frScan.next() // Discards extra input to indicate a negative value of frStart.
                }
                forceResiduesStart = frScan.nextInt()
                forceResiduesEnd = frScan.nextInt()
            } catch (Exception ex) {
                logger.severe(String.format(" FFX shutting down: input to -fR could not be parsed as a pair of integers: %s", input))
            }
            if (forceResiduesStart > forceResiduesEnd) {
                logger.info(" Start of range higher than ending: start flipped with end.")
                int temp = forceResiduesStart
                forceResiduesStart = forceResiduesEnd
                forceResiduesEnd = temp
            }
            if (forceResiduesEnd < 1) {
                logger.severe(String.format(" FFX shutting down: end range for -fR must be at least 1; input range %d to %d",
                        forceResiduesStart, forceResiduesEnd))
            }
        }

        if (algorithm != 5) {
            if (!manyBodyOptions.listResidues.equalsIgnoreCase('none')) {
                String info = "\n Evaluating rotamers for residues "
                for (String i : resList) {
                    info += String.format("%s, ", i)
                }
                logger.info(info);
            } else if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for residues " + startResID + " to " + finalResID)
            } else {
                logger.info("\n Evaluating rotamers for all residues beginning at " + allStartResID)
            }
        } else {
            if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1))
            } else {
                logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1))
            }
        }

        algorithmFunctions.open(modelfilename)
        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);

        def mapFiles = []
        int nDiffractionData = 0
        if (filenames.size() > 1) {
            String dataFileName = filenames.get(1);
            if (FilenameUtils.isExtension(dataFileName, "map")) {
                RealSpaceFile realspacefile = new RealSpaceFile(dataFileName, 1.0);
                mapFiles.add(realspacefile);
            } else {
                DiffractionFile diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
                DiffractionData diffractionData = new DiffractionData(activeAssembly, activeAssembly.getProperties(), CrystalReciprocalSpace.SolventModel.POLYNOMIAL, diffractionFile);
                diffractionData.scaleBulkFit();
                diffractionData.printStats();
                String mapFileName = String.format("%s_ffx_%d", FilenameUtils.removeExtension(dataFileName), ++nDiffractionData);
                diffractionData.writeMaps(mapFileName)
                mapFiles.add(new RealSpaceFile(mapFileName + "_2fofc.map", 1.0));
            }
        } else {
            RealSpaceFile realspacefile = new RealSpaceFile(activeAssembly);
            mapFiles.add(realspacefile);
        }


        RealSpaceData realSpaceData = new RealSpaceData(activeAssembly,
                activeAssembly.getProperties(), activeAssembly.getParallelTeam(),
                mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
        RefinementEnergy refinementEnergy = new RefinementEnergy(realSpaceData, RefinementMinimize.RefinementMode.COORDINATES, null);

        double[] x = new double[refinementEnergy.getNumberOfVariables()];
        x = refinementEnergy.getCoordinates(x);
        refinementEnergy.energy(x, true);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, refinementEnergy, sh);
        rotamerOptimization.setTwoBodyCutoff(twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setRevert(revert);
        rotamerOptimization.setPruning(pruningType);
        rotamerOptimization.setWindowSize(windowSize);
        rotamerOptimization.setIncrement(increment);
        rotamerOptimization.setDistanceCutoff(distance);
        rotamerOptimization.setNucleicCorrectionThreshold(nucleicCorrectionThreshold);
        rotamerOptimization.setMinimumNumberAcceptedNARotamers(minimumNumberAcceptedNARotamers);
        rotamerOptimization.setVerboseEnergies(verboseEnergies);
        rotamerOptimization.setBoxBorderSize(boxBorderSize);
        rotamerOptimization.setApproxBoxLength(approxBoxLength);
        rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
        rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
        rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
        rotamerOptimization.setMonteCarlo(monteCarlo, nMCSteps);

        if (useEnergyRestart) {
            rotamerOptimization.setEnergyRestartFile(energyRestartFile)
        }

        if (library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards)
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson)
        }

        if (useOrigCoordsRotamer) {
            rLib.setUseOrigCoordsRotamer(true);
        }

        int counter = 1;
        if (algorithm != 5) {
            if (allStartResID > 0) {
                ArrayList<Residue> residueList = new ArrayList<Residue>();
                Polymer[] polymers = activeAssembly.getChains();
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
                        } else if (!manyBodyOptions.forceResidues.equalsIgnoreCase('none')) {
                            if (counter >= allStartResID && counter >= forceResiduesStart
                                    && counter <= forceResiduesEnd) {
                                residueList.add(residue);
                            }
                        }
                        counter++
                    }
                }
                rotamerOptimization.setResidues(residueList);
            } else if (!manyBodyOptions.listResidues.equalsIgnoreCase('none')) {
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
            } else if (!chain.equalsIgnoreCase('-1')) {
                rotamerOptimization.setResidues(chain, startResID, finalResID);
            } else {
                rotamerOptimization.setResidues(startResID, finalResID);
            }
        } else {
            boolean ignoreNA = false;
            String ignoreNAProp = System.getProperty("ignoreNA");
            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
                ignoreNA = true;
            }
            ArrayList<Residue> residueList = new ArrayList<Residue>();
            Polymer[] polymers = activeAssembly.getChains();
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

        if (!manyBodyOptions.sequence.equalsIgnoreCase('none')) {
            for (String s : sequenceOptimizationList) {
                Character chainID = s.charAt(0);
                int num = Integer.parseInt(s.substring(1));
                for (int i = 0; i < residueList.size(); i++) {
                    Residue res = residueList.get(i);
                    if (res.getChainID() == chainID && res.getResidueNumber() == num) {
                        MultiResidue multiRes = new MultiResidue(res, activeAssembly.getForceField(), activeAssembly.getPotentialEnergy());
                        for (Polymer polymer : activeAssembly.getChains()) {
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
                        activeAssembly.getPotentialEnergy().reInit();
                        residueList.remove(i);
                        residueList.add(i, multiRes);
                    }
                }
            }
        }

        if (!manyBodyOptions.titrationOptimization.equalsIgnoreCase('none')) {
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

        x = refinementEnergy.getCoordinates(x);
        refinementEnergy.energy(x, true);

        RotamerLibrary.measureRotamers(residueList, false);

        if (decomposeOriginal) {
            rLib.setUseOrigCoordsRotamer(true);
            rotamerOptimization.setDecomposeOriginal(true);
        }

        if (algorithm == 1) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT);
        } else if (algorithm == 2) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
        } else if (algorithm == 3) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BRUTE_FORCE);
        } else if (algorithm == 4) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW);
        } else if (algorithm == 5) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX);
        }

        if (master) {
            if (min) {
                manyBodyOptions.minimize(eps);
            }
            logger.info(" Final Minimum Energy");
            algorithmFunctions.energy(activeAssembly);
            String ext = FilenameUtils.getExtension(modelfilename);
            modelfilename = FilenameUtils.removeExtension(modelfilename);
            if (ext.toUpperCase().contains("XYZ")) {
                algorithmFunctions.saveAsXYZ(activeAssembly, new File(modelfilename + ".xyz"));
            } else {
                algorithmFunctions.saveAsPDB(activeAssembly, new File(modelfilename + ".pdb"));
            }
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

