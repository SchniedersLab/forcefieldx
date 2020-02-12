//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.groovy

import edu.rit.pj.Comm
import ffx.algorithms.cli.RepexOSTOptions
import ffx.algorithms.cli.ThermodynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.thermodynamics.LambdaReader
import ffx.algorithms.thermodynamics.MonteCarloOST
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.algorithms.thermodynamics.RepExOST
import ffx.algorithms.thermodynamics.SynchronousSend
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine
import ffx.algorithms.groovy.Thermodynamics

import java.util.stream.Collectors;

/**
 * The Thermodynamics script uses the Transition-Tempered Orthogonal Space Random Walk
 * algorithm to estimate a free energy.
 * <br>
 * Usage:
 * <br>
 * ffxc Thermodynamics [options] &lt;filename&gt [file2...];
 */
@CommandLine.Command(description = " Use Orthogonal Space Tempering with histogram replica exchange to estimate a free energy.", name = "ffxc RepexThermo")
class RepexThermo extends Thermodynamics {

    @CommandLine.Mixin
    RepexOSTOptions repex;

    private RepExOST repExOST;
    private OrthogonalSpaceTempering[] allOST;
    private MonteCarloOST[] allMC;
    private MolecularDynamics[] allMD;
    // Can be either identical to allOST, or a set of Barostats built on top of them.
    private CrystalPotential[] allPotentials;

    @Override
    RepexThermo run() {

        // Begin boilerplate "make a topology" code.
        if (!init()) {
            return
        }

        boolean fromActive;
        List<String> arguments;
        if (filenames) {
            arguments = filenames;
            fromActive = false;
        } else {
            logger.warning(" Untested: use of active assembly instead of provided filenames!");
            MolecularAssembly mola = algorithmFunctions.getActiveAssembly();
            if (mola == null) {
                return helpString()
            }
            arguments = Collections.singletonList(mola.getFile().getName())
            fromActive = true;
        }

        int nArgs = arguments.size()

        topologies = new MolecularAssembly[nArgs]

        int numParallel = topology.getNumParallel(threadsAvail, nArgs)
        threadsPer = (int) (threadsAvail / numParallel)

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
        The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nArgs == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OST parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        List<MolecularAssembly> topologyList = new ArrayList<>(nArgs)

        Comm world = Comm.world()
        int size = world.size()
        if (size < 2) {
            logger.severe(" RepexThermo requires multiple processes, found only one!");
        }
        int rank = (size > 1) ? world.rank() : 0

        double initLambda = alchemical.getInitialLambda(size, rank)

        // Segment of code for MultiDynamics and OST.
        List<File> structureFiles = arguments.stream().
                map { fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath()) }.
                collect(Collectors.toList())

        File firstStructure = structureFiles.get(0)
        String filePathNoExtension = firstStructure.getAbsolutePath().replaceFirst(~/\.[^.]+$/, "")

        // SEGMENT DIFFERS FROM THERMODYNAMICS.

        String filepath = FilenameUtils.getFullPath(filePathNoExtension);
        String fileBase = FilenameUtils.getBaseName(FilenameUtils.getName(filePathNoExtension));
        String rankDirName = String.format("%s%d", filepath, rank);
        File rankDirectory = new File(rankDirName);
        if (!rankDirectory.exists()) {
            rankDirectory.mkdir();
        }
        rankDirName = "${rankDirName}${File.separator}";
        String withRankName = "${rankDirName}${fileBase}";

        File lambdaRestart = new File("${withRankName}.lam");
        boolean lamExists = lambdaRestart.exists();

        int histoRank = rank;
        if (lamExists) {
            LambdaReader lr = new LambdaReader(new BufferedReader(new FileReader(lambdaRestart)));
            lr.readLambdaFile(false);
            histoRank = lr.getHistogramIndex();
            lr.close();
        }
        File dynFile = new File("${withRankName}.dyn");

        String histoDirName = "${filepath}${histoRank}";
        File histoDir = new File(histoDirName);
        if (!histoDir.exists()) {
            histoDir.mkdir();
        }
        histoDirName = "${histoDirName}${File.separator}";
        String withHistoName = "${histoDirName}${fileBase}";

        File histogramRestart = new File("${withHistoName}.his");
        File[] trajectoryFiles = new File[nArgs];
        for (int i = 0; i < nArgs; i++) {
            trajectoryFiles[i] = new File(String.format("%s%s.arc", histoDirName, FilenameUtils.getBaseName(args[i])));
        }

        // Read in files.
        logger.info(String.format(" Initializing %d topologies...", nArgs))
        if (fromActive) {
            topologyList.add(alchemical.processFile(topology, mola, 0))
        } else {
            for (int i = 0; i < nArgs; i++) {
                topologyList.add(multidynamics.openFile(algorithmFunctions, topology,
                        threadsPer, arguments.get(i), i, alchemical, structureFiles.get(i), rank))
            }
        }

        MolecularAssembly[] topologies = topologyList.toArray(new MolecularAssembly[topologyList.size()])

        StringBuilder sb = new StringBuilder("\n Running ");
        switch (thermodynamics.getAlgorithm()) {
        // Labeled case blocks needed because Groovy (can't tell the difference between a closure and an anonymous code block).
            case ThermodynamicsOptions.ThermodynamicsAlgorithm.OST:
                ostAlg: {
                    sb.append("Orthogonal Space Tempering");
                }
                break;
            default:
                defAlg: {
                    throw new IllegalArgumentException(" RepexThermo currently does not support fixed-lambda alchemy!");
                }
                break;
        }
        sb.append(" with histogram replica exchange for ");

        potential = (CrystalPotential) topology.assemblePotential(topologies, threadsAvail, sb)
        LambdaInterface linter = (LambdaInterface) potential
        logger.info(sb.toString())

        double[] x = new double[potential.getNumberOfVariables()]
        potential.getCoordinates(x)
        linter.setLambda(initLambda)
        potential.energy(x, true)

        if (nArgs == 1) {
            randomSymop.randomize(topologies[0], potential)
        }

        multidynamics.distribute(topologies, potential, algorithmFunctions, rank, size)

        allOST = new OrthogonalSpaceTempering[size];
        allMD = new MolecularDynamics[size];
        allPotentials = new CrystalPotential[size];
        boolean isMC = ostOptions.isMc();
        boolean twoStep = ostOptions.isTwoStep();
        allMC = isMC ? new MonteCarloOST[size] : null;

        if (thermodynamics.getAlgorithm() == ThermodynamicsOptions.ThermodynamicsAlgorithm.OST) {
            for (int i = 0; i < size; i++) {
                File rankIHisto = new File("${filepath}${i}${File.separator}${fileBase}.his");
                boolean hisExists = rankIHisto.exists();

                orthogonalSpaceTempering = ostOptions.constructOST(potential, lambdaRestart, rankIHisto, topologies[0],
                        additionalProperties, dynamics, thermodynamics, algorithmListener, false, i);
                allOST[i] = orthogonalSpaceTempering;

                if (!lamExists) {
                    orthogonalSpaceTempering.setLambda(initLambda);
                }
                CrystalPotential ostPotential = ostOptions.applyAllOSTOptions(orthogonalSpaceTempering, topologies[0],
                        dynamics, lambdaParticle, barostat, hisExists)
                allPotentials[i] = ostPotential;
                if (isMC) {
                    allMC[i] = ostOptions.setupMCOST(orthogonalSpaceTempering, topologies, dynamics, thermodynamics, verbose, algorithmListener);
                    allMD[i] = allMC[i].getMD();
                } else {
                    allMD[i] = ostOptions.assembleMolecularDynamics(topologies, ostPotential, dynamics, algorithmListener);
                }
            }

            if (isMC) {
                repExOST = RepExOST.repexMC(allOST, allMC, dynamics, writeout.getFileType(), twoStep, repex.getRepexFrequency());
            } else {
                repExOST = RepExOST.repexMD(allOST, allMD, dynamics, writeout.getFileType(), repex.getRepexFrequency());
            }

            repExOST.mainLoop(thermodynamics.getEquilSteps(), true);
            repExOST.mainLoop(dynamics.getNumSteps(), false);
        } else {
            logger.severe(" RepexThermo currently does not support fixed-lambda alchemy!")
        }

        logger.info(" Algorithm with Histogram Replica Exchange Done: " + thermodynamics.getAlgorithm());

        return this

    }


    @Override
    OrthogonalSpaceTempering getOST() {
        return repex == null ? null : repex.getCurrentOST();
    }

    @Override
    CrystalPotential getPotential() {
        return (repex == null) ? potential : repex.getCurrentOST();
    }

    @Override
    List<Potential> getPotentials() {
        if (repex == null) {
            return potential == null ? Collections.emptyList() : Collections.singletonList(potential);
        }
        return Arrays.stream(repex.getAllOST()).collect(Collectors.toList());
    }
}
