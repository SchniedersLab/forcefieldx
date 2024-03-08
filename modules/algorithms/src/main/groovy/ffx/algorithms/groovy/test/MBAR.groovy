//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.algorithms.groovy.test

import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.EstimateBootstrapper
import ffx.potential.parsers.MBARFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option
import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.estimator.MultistateBennettAcceptanceRatio
import ffx.numerics.estimator.MultistateBennettAcceptanceRatio.*

/**
 * Simple wrapper for the MBAR class and does not support energy evaluations, which need to be precomputed.
 * <br>
 * Usage:
 * <br>
 * ffxc test.MBAR [options] &lt;path&gt
 */
@Command(description = " Evaluates a free energy change with the Multistate Bennett Acceptance Ratio algorithm using pregenerated snapshot energy evaluations.",
        name = "test.MBAR")
class MBAR extends AlgorithmsScript {

    @Option(names = ["--bar"], paramLabel = "false",
            description = "Run BAR calculation as well.")
    boolean bar = false

    @Option(names = ["--numBootstrap", "--nb"], paramLabel = "50",
            description = "Number of bootstrap samples to use.")
    int numBootstrap = 50

    @Option(names = ["--seed"], paramLabel = "BAR",
            description = "Seed MBAR calculation with this: ZEROS, ZWANZIG, BAR.")
    String seedWith = "BAR"

    /**
     * The path to MBAR/BAR files.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'Path to MBAR/BAR files.')
    String path = null

    public MultistateBennettAcceptanceRatio mbar = null

    /**
     * MBAR Constructor.
     */
    MBAR() {
        this(new Binding())
    }

    /**
     * MBAR Constructor.
     * @param binding The Groovy Binding to use.
     */
    MBAR(Binding binding) {
        super(binding)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    MBAR run() {
        // Begin boilerplate code.
        if (!init()) {
            return this
        }
        if (path == null) {
            logger.severe("No path to MBAR/BAR files specified.")
            return this
        }
        File path = new File(this.path)
        if (!path.exists()) {
            logger.severe("Path to MBAR/BAR files does not exist: " + path)
            return this
        }
        if (!path.isDirectory()) {
            logger.severe("Path to MBAR/BAR files is not a directory: " + path)
            return this
        }

        MBARFilter filter = new MBARFilter(path);
        seedWith = seedWith.toUpperCase()
        SeedType seed = SeedType.valueOf(seedWith) as MultistateBennettAcceptanceRatio.SeedType
        if (seed == null) {
            logger.severe("Invalid seed type: " + seedWith)
            return this
        }

        MultistateBennettAcceptanceRatio mbar = filter.getMBAR(seed)
        this.mbar = mbar
        if (mbar == null) {
            logger.severe("Could not create MBAR object.")
            return this
        }
        logger.info("\n MBAR Results:")
        logger.info(String.format(" Total dG = %10.4f +/- %10.4f kcal/mol", mbar.getFreeEnergy(),
                mbar.getUncertainty()))
        double[] dGs = mbar.getBinEnergies()
        double[] uncertainties = mbar.getBinUncertainties()
        double[][] uncertaintyMatrix = mbar.getDiffMatrix()
        for (int i = 0; i < dGs.length; i++) {
            logger.info(String.format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
        }
        logger.info("\n MBAR uncertainty between all i & j: ")
        for(int i = 0; i < uncertaintyMatrix.length; i++) {
            StringBuilder sb = new StringBuilder()
            sb.append("    [")
            for(int j = 0; j < uncertaintyMatrix[i].length; j++) {
                sb.append(String.format(" %6.5f ", uncertaintyMatrix[i][j]))
            }
            sb.append("]")
            logger.info(sb.toString())
        }
        logger.info("\n")
        if(bar){
            try {
                logger.info("\n BAR Results:")
                BennettAcceptanceRatio bar = mbar.getBAR()
                logger.info(String.format(" Total dG = %10.4f +/- %10.4f kcal/mol", bar.getFreeEnergy(),
                        bar.getUncertainty()))
                dGs = bar.getBinEnergies()
                uncertainties = bar.getBinUncertainties()
                for (int i = 0; i < dGs.length; i++) {
                    logger.info(String.format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                }
            } catch (Exception e) {
                logger.warning(" BAR calculation failed to converge.")
            }
        }
        logger.info("\n")
        if(numBootstrap != 0) {
            EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar)
            bootstrapper.bootstrap(numBootstrap)
            logger.info("\n MBAR Bootstrap Results from " + numBootstrap + " Samples:")
            logger.info(String.format(" Total dG = %10.4f +/- %10.4f kcal/mol", bootstrapper.getTotalFE(),
                    bootstrapper.getTotalUncertainty()))
            dGs = bootstrapper.getFE()
            uncertainties = bootstrapper.getUncertainty()
            for (int i = 0; i < dGs.length; i++) {
                logger.info(String.format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
            }
            logger.info("\n")
            if (bar) {
                try {
                    logger.info("\n BAR Bootstrap Results:")
                    bootstrapper = new EstimateBootstrapper(mbar.getBAR())
                    bootstrapper.bootstrap(numBootstrap)
                    logger.info(String.format(" Total dG = %10.4f +/- %10.4f kcal/mol", bootstrapper.getTotalFE(),
                            bootstrapper.getTotalUncertainty()))
                    dGs = bootstrapper.getFE()
                    uncertainties = bootstrapper.getUncertainty()
                    for (int i = 0; i < dGs.length; i++) {
                        logger.info(String.format("    dG_%d = %10.4f +/- %10.4f kcal/mol", i, dGs[i], uncertainties[i]))
                    }
                } catch (Exception e) {
                    logger.warning(" BAR calculation failed to converge.")
                }
            }
        }
        return this
    }

    MultistateBennettAcceptanceRatio getMBAR() {
        return
    }
}
