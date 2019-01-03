/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.cli;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.RotamerOptimization;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;

import picocli.CommandLine;

/**
 * Represents command line options for scripts that can create multiple walkers,
 * such as multi-walker OSRW. Should be kept agnostic to whether it is an MD-based
 * algorithm, or some other flavor of Monte Carlo.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MultiDynamicsOptions {

    private static final Logger logger = Logger.getLogger(MultiDynamicsOptions.class.getName());

    /**
     * -y or --synchronous sets synchronous walker communication (not recommended)
     */
    @CommandLine.Option(names = {"-y", "--synchronous"},
            description = "Walker communication is synchronous")
    private boolean synchronous = false;

    /**
     * -dw or --distributeWalkers allows walkers to start from multiple
     * conformations; AUTO picks up per-walker conformations as
     * filename.pdb_(walker number), and specifying a residue starts a
     * rotamer optimization to generate side-chain configurations to start
     * from.
     */
    @CommandLine.Option(names = {"--dw", "--distributeWalkers"}, paramLabel = "OFF", description = "AUTO: Pick up per-walker configurations as [filename.pdb]_[num], or specify a residue to distribute on.")
    private String distributeWalkersString = "OFF";

    /**
     * <p>isSynchronous.</p>
     *
     * @return a boolean.
     */
    public boolean isSynchronous() {
        return synchronous;
    }

    /**
     * Opens a file and processes it. Extends the behavior of AlchemicalOptions.openFile
     * by permitting use of a rank-dependent File.
     *
     * @param afuncts       AlgorithmFunctions object.
     * @param topOptions    Optional Topology Options.
     * @param threadsPer    Threads to use per system.
     * @param toOpen        Filename to open.
     * @param topNum        Number of the topology to open.
     * @param alchemy       Alchemical Options.
     * @param rank          Rank in the world communicator.
     * @param structureFile a {@link java.io.File} object.
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    public MolecularAssembly openFile(AlgorithmFunctions afuncts, Optional<TopologyOptions> topOptions, int threadsPer, String toOpen, int topNum, AlchemicalOptions alchemy, File structureFile, int rank) {
        boolean autoDist = distributeWalkersString.equalsIgnoreCase("AUTO");

        if (autoDist) {
            String openName = String.format("%s_%d", toOpen, rank + 1);
            File testFile = new File(openName);
            if (testFile.exists()) {
                toOpen = openName;
            } else {
                logger.warning(String.format(" File %s does not exist; using default %s", openName, toOpen));
            }
        }
        MolecularAssembly assembly = alchemy.openFile(afuncts, topOptions, threadsPer, toOpen, topNum);
        assembly.setFile(structureFile);
        return assembly;
    }

    /**
     * Parses --dw into optimization tokens if it's not "OFF", "AUTO", or null.
     *
     * @return
     */
    private String[] parseDistributed() {
        if (distributeWalkersString.equalsIgnoreCase("OFF") ||
                distributeWalkersString.equalsIgnoreCase("AUTO") ||
                distributeWalkersString.isEmpty()) {
            return null;
        }
        return distributeWalkersString.split("\\.");
    }

    /**
     * If residues selected for distributing initial configurations, performs many-body optimization for this distribution.
     *
     * @param topologies an array of {@link ffx.potential.MolecularAssembly} objects.
     * @param cpot       Overall CrystalPotential in use.
     * @param afuncts    a {@link ffx.algorithms.AlgorithmFunctions} object.
     * @param rank       a int.
     * @param worldSize  a int.
     */
    public void distribute(MolecularAssembly[] topologies, CrystalPotential cpot, AlgorithmFunctions afuncts, int rank, int worldSize) {
        int ntops = topologies.length;
        Potential[] energies = new Potential[ntops];
        for (int i = 0; i < ntops; i++) {
            energies[i] = topologies[i].getPotentialEnergy();
        }
        distribute(topologies, energies, cpot, afuncts, rank, worldSize);
    }

    /**
     * If residues selected for distributing initial configurations, performs many-body optimization for this distribution.
     *
     * @param topologies an array of {@link ffx.potential.MolecularAssembly} objects.
     * @param energies   ForceFieldEnergy for each topology.
     * @param cpot       Overall CrystalPotential in use.
     * @param afuncts    a {@link ffx.algorithms.AlgorithmFunctions} object.
     * @param rank       a int.
     * @param worldSize  a int.
     */
    public void distribute(MolecularAssembly[] topologies, Potential[] energies, CrystalPotential cpot, AlgorithmFunctions afuncts, int rank, int worldSize) {
        if (!distributeWalkersString.equalsIgnoreCase("AUTO") && !distributeWalkersString.equalsIgnoreCase("OFF")) {
            logger.info(" Distributing walker conformations.");
            int nSys = topologies.length;
            assert nSys == energies.length;
            switch (nSys) {
                case 1: {
                    optStructure(topologies[0], cpot, afuncts, rank, worldSize);
                }
                break;

                case 2: {
                    DualTopologyEnergy dte = (DualTopologyEnergy) cpot;
                    if (dte.getNumSharedVariables() == dte.getNumberOfVariables()) {
                        logger.info(" Generating starting structures based on dual-topology:");
                        optStructure(topologies[0], dte, afuncts, rank, worldSize);
                    } else {
                        logger.info(" Generating separate starting structures for each topology of the dual toplogy:");
                        optStructure(topologies[0], energies[0], afuncts, rank, worldSize);
                        optStructure(topologies[1], energies[1], afuncts, rank, worldSize);
                    }
                }
                break;

                case 4: {
                    QuadTopologyEnergy qte = (QuadTopologyEnergy) cpot;
                    optStructure(topologies[0], qte.getDualTopA(), afuncts, rank, worldSize);
                    optStructure(topologies[3], qte.getDualTopB(), afuncts, rank, worldSize);
                }
                break;

                // Oct-topology is deprecated on account of not working as intended.
                default: {
                    logger.severe(" First: must have 1, 2, or 4 topologies.");
                }
                break;
            }
        } else {
            logger.finer(" Skipping RO-based distribution of initial configurations.");
        }
    }

    /**
     * Distribute side-chain conformations of mola.
     *
     * @param mola To distribute
     * @param pot  Potential to use
     */
    private void optStructure(MolecularAssembly mola, Potential pot, AlgorithmFunctions aFuncts, int rank, int worldSize) {
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
        String[] distribRes = parseDistributed();

        if (distribRes == null || distribRes.length == 0) {
            throw new IllegalArgumentException(" Programming error: Must have list of residues to split on!");
        }

        LambdaInterface linter = (pot instanceof LambdaInterface) ? (LambdaInterface) pot : null;
        double initLam = -1.0;
        if (linter != null) {
            initLam = linter.getLambda();
            linter.setLambda(0.5);
        }

        Pattern chainMatcher = Pattern.compile("^([a-zA-Z])?([0-9]+)$");

        List<Residue> residueList = new ArrayList<>(distribRes.length);

        for (String ts : distribRes) {
            Matcher m = chainMatcher.matcher(ts);
            Character chainID;
            int resNum;
            if (m.find()) {
                if (m.groupCount() == 2) {
                    chainID = m.group(1).charAt(0);
                    resNum = Integer.parseInt(m.group(2));
                } else {
                    chainID = ' ';
                    resNum = Integer.parseInt(m.group(1));
                }
            } else {
                logger.warning(String.format(" Could not parse %s as a valid residue!", ts));
                continue;
            }
            logger.info(String.format(" Looking for chain %c residue %d", chainID, resNum));

            for (Polymer p : mola.getChains()) {
                if (p.getChainID() == chainID) {
                    for (Residue r : p.getResidues()) {
                        if (r.getResidueNumber() == resNum && r.getRotamers(rLib) != null) {
                            residueList.add(r);
                        }
                    }
                }
            }
        }

        if (residueList.isEmpty()) {
            throw new IllegalArgumentException(" No valid entries for distWalkers!");
        }

        AlgorithmListener alist = aFuncts.getDefaultListener();
        RotamerOptimization ropt = new RotamerOptimization(mola, pot, alist);

        ropt.setThreeBodyEnergy(false);
        if (System.getProperty("ro-ensembleNumber") == null && System.getProperty("ro-ensembleEnergy") == null) {
            logger.info(String.format(" Setting ensemble to default of number of walkers %d", worldSize));
            ropt.setEnsemble(worldSize);
        }
        ropt.setPrintFiles(false);
        ropt.setResiduesIgnoreNull(residueList);

        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(false);
        RotamerLibrary.measureRotamers(residueList, false);

        String oldLazyMat = System.getProperty("ro-lazyMatrix");
        System.setProperty("ro-lazyMatrix", "true");

        ropt.optimize(RotamerOptimization.Algorithm.ALL);
        ropt.setCoordinatesToEnsemble(rank);

        // One final energy call to ensure the coordinates are properly set at the
        // end of rotamer optimization.
        double[] xyz = new double[pot.getNumberOfVariables()];
        pot.getCoordinates(xyz);
        logger.info(" Final Optimized Energy:");
        pot.energy(xyz, true);

        if (linter != null) {
            linter.setLambda(initLam);
        }

        if (oldLazyMat != null) {
            System.setProperty("ro-lazyMatrix", oldLazyMat);
        } else {
            System.clearProperty("ro-lazyMatrix");
        }
    }
}
