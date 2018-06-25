package ffx.potential.test

import edu.rit.pj.ParallelTeam

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.PotentialScript;
import ffx.potential.cli.TopologyOptions;

import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

/**
 * The LambdaGradient script tests numeric gradients w.r.t. lambda against
 * numeric, finite-difference gradients
 * <br>
 * Usage:
 * <br>
 * ffxc test.LambdaGradient [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Test potential energy derivatives with respect to Lambda.", name = "ffxc test.LambdaGradient")
class LambdaGradient extends PotentialScript {

    @Mixin
    private AlchemicalOptions alchemical;

    @Mixin
    private TopologyOptions topology;

    /**
     * -v or --verbose is a flag to print out energy at each step.
     */
    @Option(names = ['-v', '--verbose'], paramLabel = 'false', description = 'Print out the energy for each step.')
    boolean print = false

    /**
     * -qi or --quasi-internal sets use of the quasi-internal multipole tensor formulation; not presently recommended for production use.
     */
    @Option(names = ['--qi', '--quasi-internal'], paramLabel = 'false', description = 'Use quasi-internal multipole tensors.')
    boolean qi = false

    /**
     * -d or --dx sets the finite-difference step size (in Angstroms).
     */
    @Option(names = ['-d', '--dx'], paramLabel = '1.0e-5', description = 'Finite-difference step size (Angstroms).')
    double step = 1.0e-5

    /**
     * -sdX or --skipdX disables the very long atomic gradient tests and
     * only performs dU/dL, dU2/dL2, and dU2/dLdX tests. Useful if the
     * script is being used to generate a dU/dL profile.
     */
    @Option(names = ['--sdX', '--skipdX'], paramLabel = 'false',
            description = 'Skip calculating per-atom dUdX values and only test lambda gradients.')
    boolean skipAtomGradients = false

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private int threadsAvail = ParallelTeam.getDefaultThreadCount();
    private int threadsPer = threadsAvail;
    MolecularAssembly[] topologies;

    /**
     * Script run method.
     */
    def run() {

        if (!init()) {
            return
        }

        List<String> arguments = filenames;
        // Check nArgs should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1;
        nArgs = (nArgs < 1) ? 1 : nArgs;

        topologies = new MolecularAssembly[nArgs];

        int numParallel = topology.getNumParallel(threadsAvail, nArgs);
        threadsPer = threadsAvail / numParallel;

        if (qi) {
            System.setProperty("pme.qi", "true");
            System.setProperty("no-ligand-condensed-scf", "false");
            System.setProperty("ligand-vapor-elec", "false");
            System.setProperty("polarization", "NONE");
        }
        boolean debug = (System.getProperty("debug") != null);

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        boolean lambdaTerm = (nArgs == 1 || alchemical.hasSoftcore() || topology.hasSoftcore());

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        alchemical.getInitialLambda();

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false");
        }

        List<MolecularAssembly> topologyList = new ArrayList<>(4);

        /**
         * Read in files.
         */
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = potentialFunctions.getActiveAssembly();
            if (mola == null) {
                return helpString();
            }
            arguments = new ArrayList<>();
            arguments.add(mola.getFile().getName());
            topologyList.add(alchemical.processFile(Optional.of(topology), mola, 0));
        } else {
            logger.info(String.format(" Initializing %d topologies...", nArgs));
            for (int i = 0; i < nArgs; i++) {
                topologyList.add(alchemical.openFile(potentialFunctions, Optional.of(topology), threadsPer, arguments.get(i), i));
            }
        }

        MolecularAssembly[] topologies = topologyList.toArray(new MolecularAssembly[topologyList.size()]);

        /**
         * Configure the potential to test.
         */
        StringBuilder sb = new StringBuilder("\n Testing lambda derivatives for ");
        Potential potential = topology.assemblePotential(topologies, threadsAvail, sb);

        logger.info(sb.toString());

        LambdaInterface linter = (LambdaInterface) potential;

        // End boilerplate open-topologies code.

        // Reset the number of variables for the case of dual topology.
        int n = potential.getNumberOfVariables();
        double[] x = new double[n];
        double[] gradient = new double[n];
        double[] lambdaGrad = new double[n];
        double[][] lambdaGradFD = new double[2][n];

        // Number of independent atoms.
        assert (n % 3 == 0);
        int nAtoms = n / 3;

        // Compute the Lambda = 0.0 energy.
        double lambda = 0.0;
        linter.setLambda(lambda);
        potential.getCoordinates(x);
        double e0 = potential.energyAndGradient(x, gradient);

        // Compute the Lambda = 1.0 energy.
        lambda = 1.0;
        linter.setLambda(lambda);
        double e1 = potential.energyAndGradient(x, gradient);

        logger.info(String.format(" E(0):      %20.8f.", e0));
        logger.info(String.format(" E(1):      %20.8f.", e1));
        logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1 - e0));

        // Finite-difference step size.
        double width = 2.0 * step;

        // Error tolerence
        double errTol = 1.0e-3;
        // Upper bound for typical gradient sizes (expected gradient)
        double expGrad = 1000.0;

        // Test Lambda gradient in the neighborhood of the lambda variable.
        for (int j = 0; j < 3; j++) {
            lambda = alchemical.initialLambda - 0.01 + 0.01 * j;

            if (lambda - step < 0.0) {
                continue;
            }
            if (lambda + step > 1.0) {
                continue;
            }

            logger.info(String.format(" Current lambda value %6.4f", lambda));
            linter.setLambda(lambda);

            // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
            double e = potential.energyAndGradient(x, gradient);

            // Analytic dEdL, d2E/dL2 and dE/dL/dX
            double dEdL = linter.getdEdL();
            double d2EdL2 = linter.getd2EdL2();
            for (int i = 0; i < n; i++) {
                lambdaGrad[i] = 0.0;
            }
            linter.getdEdXdL(lambdaGrad);

            // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
            linter.setLambda(lambda + step);
            double lp = potential.energyAndGradient(x, lambdaGradFD[0]);
            double dedlp = linter.getdEdL();
            linter.setLambda(lambda - step);
            double lm = potential.energyAndGradient(x, lambdaGradFD[1]);
            double dedlm = linter.getdEdL();

            double dEdLFD = (lp - lm) / width;
            double d2EdL2FD = (dedlp - dedlm) / width;

            if (debug) {
                logger.info(String.format(" db> Lambda FD Test  lower,center,upper: %g %g %g", lambda - step, lambda, lambda + step));
                logger.info(String.format(" db> dE/dL   Numeric lp,lm,width,val: %+9.6g %+9.6g %4.2g (%+9.6g)", lp, lm, width, dEdLFD));
                logger.info(String.format(" db> d2E/dL2 Numeric lp,lm,width,val: %+9.6g %+9.6g %4.2g [%+9.6g]", dedlp, dedlm, width, d2EdL2FD));
                logger.info(String.format(" db> Analytic vers   l,dEdL,d2EdL2: %4.2f (%+9.6g) [%+9.6g]", lambda, dEdL, d2EdL2));
            }

            double err = Math.abs(dEdLFD - dEdL);
            if (err < errTol) {
                logger.info(String.format(" dE/dL passed:   %10.6f", err));
            } else {
                logger.info(String.format(" dE/dL failed: %10.6f", err));
            }
            logger.info(String.format(" Numeric:   %15.8f", dEdLFD));
            logger.info(String.format(" Analytic:  %15.8f", dEdL));

            err = Math.abs(d2EdL2FD - d2EdL2);
            if (err < errTol) {
                logger.info(String.format(" d2E/dL2 passed: %10.6f", err));
            } else {
                logger.info(String.format(" d2E/dL2 failed: %10.6f", err));
            }
            logger.info(String.format(" Numeric:   %15.8f", d2EdL2FD));
            logger.info(String.format(" Analytic:  %15.8f", d2EdL2));

            //boolean passed = true
            int ndEdXdLFails = 0;

            double rmsError = 0;
            for (int i = 0; i < nAtoms; i++) {
                int ii = i * 3;
                double dX = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dXa = lambdaGrad[ii];
                double eX = dX - dXa;
                ii++;
                double dY = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dYa = lambdaGrad[ii];
                double eY = dY - dYa;
                ii++;
                double dZ = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dZa = lambdaGrad[ii];
                double eZ = dZ - dZa;

                double error = eX * eX + eY * eY + eZ * eZ;
                rmsError += error;
                error = Math.sqrt(error);
                if (error < errTol) {
                    logger.fine(String.format(" dE/dX/dL for Atom %d passed: %10.6f", i + 1, error));
                } else {
                    logger.info(String.format(" dE/dX/dL for Atom %d failed: %10.6f", i + 1, error));
                    logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa, dYa, dZa));
                    logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ));
                    ndEdXdLFails++;
                }
            }
            rmsError = Math.sqrt(rmsError / nAtoms);
            if (ndEdXdLFails == 0) {
                logger.info(String.format(" dE/dX/dL passed for all atoms: RMS error %15.8f", rmsError));
            } else {
                logger.info(String.format(" dE/dX/dL failed for %d of %d atoms: RMS error %15.8f", ndEdXdLFails, nAtoms, rmsError));
            }

            logger.info("");
        }

        linter.setLambda(alchemical.initialLambda);
        potential.getCoordinates(x);
        potential.energyAndGradient(x, gradient, print);

        if (!skipAtomGradients) {
            logger.info(String.format(" Checking Cartesian coordinate gradient"));
            double[] numeric = new double[3];
            double avLen = 0.0;
            int nFailures = 0;
            double avGrad = 0.0;
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                int i0 = i3 + 0;
                int i1 = i3 + 1;
                int i2 = i3 + 2;

                // Find numeric dX
                double orig = x[i0];
                x[i0] = x[i0] + step;
                double e = potential.energyAndGradient(x, lambdaGradFD[0], print);
                x[i0] = orig - step;
                e -= potential.energyAndGradient(x, lambdaGradFD[1], print);
                x[i0] = orig;
                numeric[0] = e / width;

                // Find numeric dY
                orig = x[i1];
                x[i1] = x[i1] + step;
                e = potential.energyAndGradient(x, lambdaGradFD[0], print);
                x[i1] = orig - step;
                e -= potential.energyAndGradient(x, lambdaGradFD[1], print);
                x[i1] = orig;
                numeric[1] = e / width;

                // Find numeric dZ
                orig = x[i2];
                x[i2] = x[i2] + step;
                e = potential.energyAndGradient(x, lambdaGradFD[0], print);
                x[i2] = orig - step;
                e -= potential.energyAndGradient(x, lambdaGradFD[1], print);
                x[i2] = orig;
                numeric[2] = e / width;

                double dx = gradient[i0] - numeric[0];
                double dy = gradient[i1] - numeric[1];
                double dz = gradient[i2] - numeric[2];
                double len = dx * dx + dy * dy + dz * dz;
                avLen += len;
                len = Math.sqrt(len);

                double grad2 = gradient[i0] * gradient[i0] + gradient[i1] * gradient[i1] + gradient[i2] * gradient[i2];
                avGrad += grad2;

                if (len > errTol) {
                    logger.info(String.format(" Atom %d failed: %10.6f.", i + 1, len)
                            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
                            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
                    ++nFailures
                    //return
                } else {
                    logger.info(String.format(" Atom %d passed: %10.6f.", i + 1, len)
                            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
                            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
                }

                if (grad2 > expGrad) {
                    logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i + 1, grad2))
                }
                logger.info("\n")
            }

            avLen = avLen / nAtoms
            avLen = Math.sqrt(avLen)
            if (avLen > errTol) {
                logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol))
            } else {
                logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol))
            }
            logger.info(String.format(" Number of atoms failing gradient test: %d", nFailures))

            avGrad = avGrad / nAtoms
            avGrad = Math.sqrt(avGrad)
            if (avGrad > expGrad) {
                logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
            } else {
                logger.info(String.format(" RMS gradient: %10.6f", avGrad))
            }
        } else {
            logger.info(" Skipping atomic dU/dX gradients.")
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
