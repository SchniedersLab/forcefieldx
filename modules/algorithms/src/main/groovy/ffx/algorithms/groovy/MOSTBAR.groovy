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

import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.BarostatOptions
import ffx.algorithms.thermodynamics.HistogramReader
import ffx.crystal.CrystalPotential
import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.SequentialEstimator
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.SystemFilter
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine

/**
 * The MOSTBAR script uses a single set of archive file(s) from a Metropolized
 * Orthogonal Space Tempering run to evaluate free energy via the Bennett Acceptance Ratio
 * <br>
 * Usage:
 * <br>
 * ffxc MOSTBAR [options] &lt;structures1&gt
 */
@CommandLine.Command(description = " Evaluates free energy of an M-OST run using the BAR estimator.", name = "ffxc MOSTBAR")
class MOSTBAR extends AlgorithmsScript {

    @CommandLine.Mixin
    private AlchemicalOptions alchemical

    @CommandLine.Mixin
    private TopologyOptions topology

    @CommandLine.Mixin
    private BarostatOptions barostat

    @CommandLine.Option(names = ["-t", "--temperature"], paramLabel = "298.15",
            description = "Temperature in Kelvins")
    private double temp = 298.15

    @CommandLine.Option(names = ["--his", "--histogram"], paramLabel = "file.his",
            description = "Manually provided path to a histogram file (otherwise, attempts to autodetect from same directory as input files).")
    private String histogramName = "";

    @CommandLine.Option(names = ["--lb", "--lambdaBins"], paramLabel = "autodetected",
            description = "Manually specified number of lambda bins (else auto-detected from histogram")
    private int lamBins = -1;

    /**
     * The final argument(s) should be filenames for lambda windows in order..
     */
    @CommandLine.Parameters(arity = "1..*", paramLabel = "files",
            description = 'Trajectory files for the first end of the window, followed by trajectories for the other end')
    List<String> filenames = null

    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail
    private MolecularAssembly[] topologies
    private SystemFilter[] openers

    private CrystalPotential potential
    private LambdaInterface linter;

    private Configuration additionalProperties

    private List<List<Double>> energiesL;
    private List<List<Double>> energiesUp;
    private List<List<Double>> energiesDown;

    private double[] lamPoints;
    private double lamSep;
    private double halfLamSep;
    private double[] x;
    private final double[] lastEntries = new double[3];
    private static final String energyFormat = "%11.4f kcal/mol";
    private static final String nanFormat = String.format("%20s", "N/A");

    void setProperties(CompositeConfiguration addedProperties) {
        additionalProperties = addedProperties;
    }

    @Override
    MOSTBAR run() {
        // Begin boilerplate code.
        if (!init()) {
            return null
        }
        if (filenames == null || filenames.isEmpty()) {
            return null
        }

        int nFiles = filenames.size()

        topologies = new MolecularAssembly[nFiles]
        openers = new SystemFilter[nFiles]

        int numParallel = topology.getNumParallel(threadsAvail, nFiles)
        threadsPer = (int) (threadsAvail / numParallel)

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
        The Minimize script, for example, may be running on a single, unscaled physical topology. */
        boolean lambdaTerm = (nFiles == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OST parameters than absolute free energies.
        if (nFiles >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }


        logger.info(String.format(" Initializing %d topologies", nFiles))
        for (int i = 0; i < nFiles; i++) {
            MolecularAssembly ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i)
            topologies[i] = ma
            openers[i] = algorithmFunctions.getFilter()
        }

        StringBuilder sb = new StringBuilder("\n Using BAR to analyze an M-OST free energy change between for systems ")
        potential = (CrystalPotential) topology.assemblePotential(topologies, threadsAvail, sb)
        potential = barostat.checkNPT(topologies[0], potential)
        linter = (LambdaInterface) potential;
        logger.info(sb.toString())

        int nSnapshots = openers[0].countNumModels();

        if (histogramName.isEmpty()) {
            histogramName = FilenameUtils.removeExtension(filenames.get(0)) + ".his";
        }

        if (lamBins < 1) {
            File histogramFile = new File(histogramName);
            if (!histogramFile.exists() || !histogramFile.canRead()) {
                logger.severe(" Histogram file ${histogramName} does not exist or could not be read!");
            }

            HistogramReader hr = null;
            try {
                hr = new HistogramReader(new BufferedReader(new FileReader(histogramFile)));
                hr.readHistogramFile();
                lamBins = hr.getLambdaBins();
                logger.info(" Autodetected ${lamBins} from histogram file.");
            } finally {
                hr?.close();
            }
        }

        energiesL = new ArrayList<>(lamBins);
        energiesUp = new ArrayList<>(lamBins);
        energiesDown = new ArrayList<>(lamBins);
        for (int i = 0; i < lamBins; i++) {
            energiesL.add(new ArrayList<Double>());
            energiesUp.add(new ArrayList<Double>());
            energiesDown.add(new ArrayList<Double>());
        }

        lamSep = 1.0 / (lamBins - 1);
        halfLamSep = 0.5 * lamSep;
        lamPoints = new double[lamBins];
        // TODO: Remove assumption that it's using discrete lambda bins.
        for (int i = 0; i < (lamBins - 1); i++) {
            lamPoints[i] = i * lamSep;
        }
        lamPoints[lamBins - 1] = 1.0; // Eliminate machine precision error.

        OptionalDouble optLam = openers[0].getLastReadLambda();
        if (optLam.isEmpty()) {
            throw new IllegalArgumentException(" No lambda records found in the first header of archive file ${filenames[0]}");
        }

        double lambda = optLam.getAsDouble();
        int nVar = potential.getNumberOfVariables();
        x = new double[nVar];

        logger.info(String.format(" Evaluating snapshot     1 of %5d", nSnapshots));
        addEntries(lambda, 0)

        for (int i = 1; i < nSnapshots; i++) {
            logger.info(String.format(" Evaluating snapshot %5d of %5d", (i+1), nSnapshots));
            for (int j = 0; j < nFiles; j++) {
                openers[j].readNext();
            }
            lambda = openers[0].getLastReadLambda().getAsDouble();
            addEntries(lambda, i);
        }

        double[][] eLow = new double[lamBins][];
        double[][] eAt = new double[lamBins][];
        double[][] eHigh = new double[lamBins][];
        for (int i = 0; i < lamBins; i++) {
            eLow[i] = energiesDown.get(i).stream().mapToDouble(Double::doubleValue).toArray();
            eAt[i] = energiesL.get(i).stream().mapToDouble(Double::doubleValue).toArray();
            eHigh[i] = energiesUp.get(i).stream().mapToDouble(Double::doubleValue).toArray();
        }

        SequentialEstimator bar = new BennettAcceptanceRatio(lamPoints, eLow, eAt, eHigh, new double[]{temp});
        SequentialEstimator forwards = bar.getInitialForwardsGuess();
        SequentialEstimator backwards = bar.getInitialBackwardsGuess();

        logger.info(String.format(" Free energy via BAR:           %15.9f +/- %.9f kcal/mol.", bar.getFreeEnergy(), bar.getUncertainty()))
        logger.warning(" FEP uncertainties in FFX are currently underestimated and unreliable!");
        logger.info(String.format(" Free energy via forwards FEP:  %15.9f +/- %.9f kcal/mol.", forwards.getFreeEnergy(), forwards.getUncertainty()));
        logger.info(String.format(" Free energy via backwards FEP: %15.9f +/- %.9f kcal/mol.", backwards.getFreeEnergy(), backwards.getUncertainty()));

        double[] barFE = bar.getWindowEnergies();
        double[] barVar = bar.getWindowUncertainties();
        double[] forwardsFE = forwards.getWindowEnergies();
        double[] forwardsVar = forwards.getWindowUncertainties();
        double[] backwardsFE = backwards.getWindowEnergies();
        double[] backwardsVar = backwards.getWindowUncertainties();
        
        sb = new StringBuilder(" Free Energy Profile\n Min_Lambda Max_Lambda          BAR_dG      BAR_Var          FEP_dG      FEP_Var     FEP_Back_dG FEP_Back_Var\n");
        for (int i = 0; i < (lamBins - 1); i++) {
            sb.append(String.format(" %-10.8f %-10.8f %15.9f %12.9f %15.9f %12.9f %15.9f %12.9f\n",
                    lamPoints[i], lamPoints[i+1], barFE[i], barVar[i], forwardsFE[i], forwardsVar[i], backwardsFE[i], backwardsVar[i]));
        }
        logger.info(sb.toString());

        return this;
    }

    private void addEntries(double lambda, int index) {
        x = potential.getCoordinates(x);
        lastEntries[0] = addLambdaDown(lambda);
        lastEntries[1] = addAtLambda(lambda);
        lastEntries[2] = addLambdaUp(lambda);

        String low = Double.isNaN(lastEntries[0]) ? nanFormat : String.format(energyFormat, lastEntries[0]);
        String high = Double.isNaN(lastEntries[2]) ? nanFormat : String.format(energyFormat, lastEntries[2]);
        logger.info(String.format(" Energies for snapshot %5d: " +
                "%s, %s, %s", (index+1), low, String.format(energyFormat, lastEntries[1]), high));
    }

    private double addAtLambda(double lambda) {
        assert lambda >= 0.0 && lambda <= 1.0;
        linter.setLambda(lambda);
        double e = potential.energy(x, false);
        int bin = binForLambda(lambda);
        energiesL.get(bin).add(e);
        return e;
    }

    private double addLambdaUp(double lambda) {
        int bin = binForLambda(lambda);
        double modLambda = lambda + lamSep;
        // DISCRETE ONLY: assert lambda == 1.0d || modLambda < (1.0 + 1.0E-6);
        modLambda = Math.min(1.0d, modLambda);
        if (bin == (lamBins - 1)) {
            energiesUp.get(bin).add(Double.NaN);
            return Double.NaN;
        } else {
            linter.setLambda(modLambda);
            double e = potential.energy(x, false);
            energiesUp.get(bin).add(e);
            linter.setLambda(lambda);
            return e;
        }
    }

    private double addLambdaDown(double lambda) {
        int bin = binForLambda(lambda);
        double modLambda = lambda - lamSep;
        // DISCRETE ONLY: assert lambda == 0.0d || modLambda > -1.0E-6;
        modLambda = Math.max(0.0d, modLambda);

        if (bin == 0) {
            energiesDown.get(0).add(Double.NaN);
            return Double.NaN;
        } else {
            linter.setLambda(modLambda);
            double e = potential.energy(x, false);
            energiesDown.get(bin).add(e);
            linter.setLambda(lambda);
            return e;
        }
    }

    /**
     * <p>binForLambda.</p>
     *
     * @param lambda a double.
     * @return a int.
     */
    private int binForLambda(double lambda) {
        // TODO: Robust to existence of half-bins.
        return (int) Math.round(lambda / lamSep);
    }
}
