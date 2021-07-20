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
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.BarostatOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.cli.RepExOptions
import ffx.algorithms.dynamics.Barostat
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.ReplicaExchange
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.WriteoutOptions
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

import static java.lang.String.format
import static java.lang.String.format
import static java.lang.String.format

/**
 * The Dynamics script implements molecular and stochastic dynamics algorithms.
 * <br>
 * Usage:
 * <br>
 * ffxc Dynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run dynamics on a system.", name = "ffxc Dynamics")
class Dynamics extends AlgorithmsScript {

    @Mixin
    AtomSelectionOptions atomSelectionOptions

    @Mixin
    DynamicsOptions dynamicsOptions

    @Mixin
    BarostatOptions barostatOptions

    @Mixin
    WriteoutOptions writeOut

    @Mixin
    RepExOptions repEx


    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    /**
     * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
     * Originally it is declared in the run method
     */
    public Potential potential = null
    public MolecularDynamics molDyn = null

    MolecularDynamics getMolecularDynamics() {
        return molDyn
    }

    Potential getPotentialObject() {
        return potential
    }

    /**
     * Dynamics Constructor.
     */
    Dynamics() {
        this(new Binding())
    }

    /**
     * Dynamics Constructor.
     * @param binding The Groovy Binding to use.
     */
    Dynamics(Binding binding) {
        super(binding)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    Dynamics run() {

        if (!init()) {
            return this
        }

        dynamicsOptions.init()

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = [algorithmFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        atomSelectionOptions.setActiveAtoms(activeAssembly)

        /* File structureFile = new File(FilenameUtils.normalize(modelFilename))
         structureFile = new File(structureFile.getAbsolutePath())
         String baseFilename = FilenameUtils.removeExtension(structureFile.getName())*/

        potential = activeAssembly.getPotentialEnergy()
        double[] x = new double[potential.getNumberOfVariables()]
        potential.getCoordinates(x)
        potential.energy(x, true)

        if (barostatOptions.pressure > 0) {
            CrystalPotential crystalPotential = (CrystalPotential) potential
            Barostat barostat = barostatOptions.createBarostat(activeAssembly, crystalPotential)
            potential = barostat
        }

        Comm world = Comm.world()
        int size = world.size()
        List<String> files = filenames
        int nFiles = filenames.size()
        nFiles = (nFiles < 1) ? 1 : nFiles

        if (!repEx.repEx || size < 2) {
            logger.info("\n Running molecular dynamics on " + modelFilename)
            // Restart File
            File dyn = new File(FilenameUtils.removeExtension(modelFilename) + ".dyn")
            if (!dyn.exists()) {
                dyn = null
            }

            molDyn = dynamicsOptions.getDynamics(writeOut, potential, activeAssembly, algorithmListener)

            molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt,
                    dynamicsOptions.report, dynamicsOptions.write, dynamicsOptions.temperature, true, dyn)

        } else {
            logger.info("\n Running replica exchange molecular dynamics on " + modelFilename)
            int rank = (size > 1) ? world.rank() : 0
            logger.info("Rank:" + rank.toString())

            List<File> structureFiles = files.stream().
                    map { fn -> new File(new File(FilenameUtils.normalize(fn)).getAbsolutePath())
                    }.
                    collect(Collectors.toList())

            File firstStructure = structureFiles.get(0)
            String filePathNoExtension = firstStructure.getAbsolutePath().replaceFirst(~/\.[^.]+$/, "")
            File histogramRestart = new File(filePathNoExtension + ".his")

            String withRankName = filePathNoExtension
            String baseName = FilenameUtils.getBaseName(filenames.get(0))

            List<File> rankedFiles = new ArrayList<>(nFiles)
            String rankDirName = FilenameUtils.getFullPath(filePathNoExtension)
            rankDirName = format("%s%d", rankDirName, rank)
            File rankDirectory = new File(rankDirName)
            if (!rankDirectory.exists()) {
                rankDirectory.mkdir()
            }
            rankDirName = rankDirName + File.separator
            withRankName = format("%s%s", rankDirName, FilenameUtils.getName(filePathNoExtension))
            logger.info("With Rank Name:" + withRankName)

            for (File structureFile : structureFiles) {
                rankedFiles.add(new File(format("%s%s", rankDirName,
                        FilenameUtils.getName(structureFile.getName()))))
            }
            structureFiles = rankedFiles
            logger.info("ranked Files:" + rankedFiles.get(0))
            File dyn = new File(withRankName + ".dyn")

            /*File rankDirectory = new File(structureFile.getParent() + File.separator
                + Integer.toString(rank))
            if (!rankDirectory.exists()) {
              rankDirectory.mkdir()
            }
            String withRankName = rankDirectory.getPath() + File.separator + baseFilename

            if (!dyn.exists()) {
              dyn = null
            }*/
            activeAssembly.setFile(rankedFiles.get(0))
            molDyn = dynamicsOptions.getDynamics(writeOut, potential, activeAssembly, algorithmListener)
            ReplicaExchange replicaExchange = new ReplicaExchange(molDyn, algorithmListener,
                    dynamicsOptions.temperature, repEx.exponent)

            long totalSteps = dynamicsOptions.steps
            int nSteps = repEx.replicaSteps
            int cycles = (int) (totalSteps / nSteps)
            if (cycles <= 0) {
                cycles = 1
            }

            replicaExchange.
                    sample(cycles, nSteps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write)
        }

        return this
    }

    /**
     * {@inheritDoc}
     */
    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (potential == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList(potential)
        }
        return potentials
    }

}
