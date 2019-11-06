//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.BarostatOptions
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.SystemFilter
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
import org.apache.commons.math3.util.FastMath
import picocli.CommandLine

/**
 * The BAR script find the free energy difference across a lambda window. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt &lt;structures2&gt;
 */
@CommandLine.Command(description = " Evaluates a free energy change with the Bennett Acceptance Ratio algorithm using pregenerated snapshots.", name = "ffxc BAR")
class BAR extends AlgorithmsScript {

    @CommandLine.Mixin
    private AlchemicalOptions alchemical

    @CommandLine.Mixin
    private TopologyOptions topology

    @CommandLine.Mixin
    private BarostatOptions barostat

    @CommandLine.Option(names = ["--l2", "--lambdaTwo"], paramLabel = "1.0",
            description = "Lambda value for the upper edge of the window")
    private double lambda2 = 1.0;

    @CommandLine.Option(names = ["--t1", "--temperature1"], paramLabel = "298.15",
            description = "Temperature for system 1")
    private double temp1 = 298.15;

    @CommandLine.Option(names = ["--t2", "--temperature2"], paramLabel = "298.15",
            description = "Temperature for system 2")
    private double temp2 = 298.15;

    /**
     * The final argument(s) should be filenames for lambda windows in order..
     */
    @CommandLine.Parameters(arity = "1..*", paramLabel = "files", description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
    List<String> filenames = null

    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail
    MolecularAssembly[] topologies1
    MolecularAssembly[] topologies2
    SystemFilter[] openers1;
    SystemFilter[] openers2;

    CrystalPotential potential1
    CrystalPotential potential2

    private Configuration additionalProperties1;
    private Configuration additionalProperties2;

    /**
     * Sets an optional Configuration with additional properties.
     * @param additionalProps
     */
    void setProperties(Configuration additionalProps1, Configuration additionalProps2) {
        this.additionalProperties1 = additionalProps1
        this.additionalProperties2 = additionalProps2;
    }

    @Override
    BAR run() {
        // Begin boilerplate code.
        if (!init()) {
            return
        }
        if (filenames == null || filenames.size() % 2 != 0) {
            return;
        }

        int nFiles = filenames.size();
        int nPer = nFiles / 2;

        topologies1 = new MolecularAssembly[nPer]
        topologies2 = new MolecularAssembly[nPer]
        openers1 = new SystemFilter[nPer]
        openers2 = new SystemFilter[nPer]

        int numParallel = topology.getNumParallel(threadsAvail, nPer)
        threadsPer = threadsAvail / numParallel

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
        The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nPer == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nPer >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        logger.info(String.format(" Initializing %d topologies for each end", nPer));
        for (int i = 0; i < nPer; i++) {
            MolecularAssembly ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i);
            topologies1[i] = ma;
            openers1[i] = algorithmFunctions.getFilter();
            int iSecond = i + nPer;
            ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[iSecond], i);
            topologies2[i] = ma;
            openers2[i] = algorithmFunctions.getFilter();
        }

        double lambda1 = alchemical.initialLambda;

        StringBuilder sb = new StringBuilder(String.format("\n Using BAR to analyze a free energy change between L=%.5f and L=%.5f for systems ", lambda1, lambda2));
        potential1 = (CrystalPotential) topology.assemblePotential(topologies1, threadsAvail, sb);
        potential1 = barostat.checkNPT(topologies1[0], potential1);
        sb.append(" and systems ");
        potential2 = (CrystalPotential) topology.assemblePotential(topologies2, threadsAvail, sb);
        potential2 = barostat.checkNPT(topologies2[0], potential2);

        LambdaInterface linter1 = (LambdaInterface) potential1;
        LambdaInterface linter2 = (LambdaInterface) potential2;
        logger.info(sb.toString());

        int nSnapshots1 = openers1[0].countNumModels();
        int nSnapshots2 = openers2[0].countNumModels();
        double[] e1L1 = new double[nSnapshots1];
        double[] e1L2 = new double[nSnapshots1];
        double[] e2L1 = new double[nSnapshots2];
        double[] e2L2 = new double[nSnapshots2];
        double[] eDiff1 = new double[nSnapshots1];
        double[] eDiff2 = new double[nSnapshots2];

        logger.info(" Preliminary energy evaluation for first end of the window.");
        double[] x1 = new double[potential1.getNumberOfVariables()];
        x1 = potential1.getCoordinates(x1);
        linter1.setLambda(lambda1);
        e1L1[0] = potential1.energy(x1, true);
        linter1.setLambda(lambda2);
        e1L2[0] = potential1.energy(x1, false);
        eDiff1[0] = e1L2[0] - e1L1[0];

        logger.info(" Preliminary energy evaluation for second end of the window.");
        double[] x2 = new double[potential2.getNumberOfVariables()];
        potential2.getCoordinates(x2);
        linter2.setLambda(lambda2);
        e2L2[0] = potential2.energy(x2, true);
        e2L1[0] = potential2.energy(x2, false);
        eDiff2[0] = e2L2[0] - e2L1[0];

        String lamString1 = String.format("%.3f", lambda1);
        String lamString2 = String.format("%.3f", lambda2);

        // TODO: Increment by stride instead of always 1.
        for (int i = 1; i < nSnapshots1; i++) {
            linter1.setLambda(lambda1);
            for (int j = 0; j < nPer; j++) {
                openers1[j].readNext(false, false);
            }
            // TODO: Make repeat measurements (esp. for mixed-precision OMM).
            x1 = potential1.getCoordinates(x1);
            e1L1[i] = potential1.energy(x1, false);
            linter1.setLambda(lambda2);
            e1L2[i] = potential1.energy(x1, false);

            eDiff1[i] = e1L2[i] - e1L1[i];
            logger.info(String.format(" Snapshot %d of system 1: E(L=%s) = %14.7f, E(L=%s) = %14.7f, difference = %14.7f",
                    i+1, lamString1, e1L1[i], lamString2, e1L2[i], eDiff1[i]));
        }

        for (int i = 1; i < nSnapshots2; i++) {
            linter2.setLambda(lambda1);
            for (int j = 0; j < nPer; j++) {
                openers2[j].readNext(false, false);
            }
            x2 = potential2.getCoordinates(x2);
            e2L1[i] = potential2.energy(x2, false);
            linter2.setLambda(lambda2);
            e2L2[i] = potential2.energy(x2, false);

            eDiff2[i] = e2L2[i] - e2L1[i];
            logger.info(String.format(" Snapshot %d of system 2: E(L=%s) = %14.7f, E(L=%s) = %14.7f, difference = %14.7f",
                    i+1, lamString1, e1L1[i], lamString2, e1L2[i], eDiff2[i]));
        }

        double mean1 = 0;
        double var1 = 0;
        double max1 = Integer.MIN_VALUE;
        double min1 = Integer.MAX_VALUE;
        for (int i = 1; i <= nSnapshots1; i++) {
            double priorMean = mean1;
            double val = eDiff1[i-1];
            mean1 += ((val - mean1) / i);
            var1 += ((val - priorMean)*(val - mean1));
            min1 = Math.min(min1, val);
            max1 = Math.max(max1, val);
        }
        double sd1 = FastMath.sqrt(var1 / (nSnapshots1 - 1));
        logger.info(String.format(" System 1 differences: mean %14.7f, sample standard deviation %14.7f, min %14.7f, max %14.7f over %d samples", mean1, sd1, min1, max1, nSnapshots1));

        double mean2 = 0;
        double var2 = 0;
        double max2 = Integer.MIN_VALUE;
        double min2 = Integer.MAX_VALUE;
        for (int i = 1; i <= nSnapshots2; i++) {
            double priorMean = mean2;
            double val = eDiff2[i-1];
            mean2 += ((val - mean2) / i);
            var2 += ((val - priorMean)*(val - mean2));
            min2 = Math.min(min2, val);
            max2 = Math.max(max2, val);
        }
        double sd2 = FastMath.sqrt(var2 / (nSnapshots2 - 1));
        logger.info(String.format(" System 2 differences: mean %14.7f, sample standard deviation %14.7f, min %14.7f, max %14.7f over %d samples", mean2, sd2, min2, max2, nSnapshots2));

        String barFileName = FilenameUtils.removeExtension(filenames.get(0)) + ".bar";
        logger.info(" Writing Tinker-compatible .bar file to ${barFileName}. For now: use Tinker's bar command; built-in FFX calculations not yet implemented.");
        File barFile = new File(barFileName);
        BufferedWriter bw;
        try {
            bw = new BufferedWriter(new FileWriter(barFile));

            StringBuilder fnames1 = new StringBuilder();
            for (int i = 0; i < nPer; i++) {
                fnames1.append("  ").append(filenames.get(i));
            }
            bw.write(String.format("%8d %9.3f%s\n", nSnapshots1, temp1, fnames1.toString()));
            for (int i = 0; i < nSnapshots1; i++) {
                bw.write(String.format("%8d %20.10f %20.10f\n", i+1, e1L1[i], e1L2[i]));
            }

            StringBuilder fnames2 = new StringBuilder();
            for (int i = nPer; i < nFiles; i++) {
                fnames2.append("  ").append(filenames.get(i));
            }
            bw.write(String.format("%8d %9.3f  %s\n", nSnapshots2, temp2, fnames2.toString()));
            for (int i = 0; i < nSnapshots2; i++) {
                bw.write(String.format("%8d %20.10f %20.10f\n", i+1, e2L1[i], e2L2[i]));
            }
        } finally {
            bw?.close();
        }

        return this;
    }

    @Override
    List<Potential> getPotentials() {
        ArrayList<Potential> pots = new ArrayList<>(2);
        pots.add(potential1);
        pots.add(potential2);
        return pots;
    }
}
