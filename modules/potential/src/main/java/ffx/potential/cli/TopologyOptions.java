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
package ffx.potential.cli;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.stream.Collectors;

import ffx.numerics.Potential;
import ffx.numerics.switching.PowerSwitch;
import ffx.numerics.switching.SquaredTrigSwitch;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.MultiplicativeSwitch;

import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize multiple physical topologies.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class TopologyOptions {

    /**
     * The logger for this class.
     */
    public static final Logger logger = Logger.getLogger(TopologyOptions.class.getName());

    /**
     * -s2 or --start2 defines the first softcored atom for the second topology.
     */
    @Option(names = {"--s2", "--start2"}, paramLabel = "0",
            description = "Starting ligand atom for 2nd topology")
    int s2 = 0;

    /**
     * -f2 or --final2 defines the last softcored atom for the second topology.
     */
    @Option(names = {"--f2", "--final2"}, paramLabel = "-1",
            description = "Final ligand atom for the 2nd topology")
    int f2 = -1;

    /**
     * --la2 or -ligAtoms2 allows for multiple ranges and/or singletons of ligand atoms in the second topology, separated by periods.
     */
    @Option(names = {"--la2", "--ligAtoms2"},
            description = "Period-separated ranges of 2nd toplogy ligand atoms (e.g. 40-50.72-83)")
    String ligAt2 = null;

    /**
     * -es2 or --noElecStart2 defines the first atom of the second topology to have no electrostatics.
     */
    @Option(names = {"--es2", "--noElecStart2"}, paramLabel = "1",
            description = "Starting no-electrostatics atom for 2nd topology")
    int es2 = 1;

    /**
     * -ef2 or --noElecFinal2 defines the last atom of the second topology to have no electrostatics.
     */
    @Option(names = {"--ef2", "--noElecFinal2"}, paramLabel = "-1",
            description = "Final no-electrostatics atom for 2nd topology")
    int ef2 = -1;

    /**
     * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, or 4.
     */
    @Option(names = {"--np", "--nParallel"}, paramLabel = "1",
            description = "Number of topologies to evaluate in parallel")
    int nPar = 1;

    /**
     * -uaA or -unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated ranges or singletons.
     */
    @Option(names = {"--uaA", "--unsharedA"},
            description = "Unshared atoms in the A dual topology (period-separated hyphenated ranges)")
    String unsharedA = null;

    /**
     * -uaB or -unsharedB sets atoms unique to the B dual-topology, as period-separated hyphenated ranges or singletons.
     */
    @Option(names = {"--uaB", "--unsharedB"},
            description = "Unshared atoms in the B dual topology (period-separated hyphenated ranges)")
    String unsharedB = null;

    /**
     * -sf or --switchingFunction sets the switching function to be used by
     * dual topologies; TRIG produces the function sin^2(pi/2*lambda)*E1(lambda)
     * + cos^2(pi/2*lambda)*E2(1-lambda), MULT uses a 5"th-order polynomial
     * switching function with zero first and second derivatives at the end
     * (same function as used for van der Waals switch), and a number uses
     * the original function, of l^beta*E1(lambda) + (1-lambda)^beta*E2(1-lambda).
     * <p>
     * All of these are generalizations of Udt = f(l)*E1(l) +
     * f(1-l)*E2(1-lambda), where f(l) is a continuous switching function
     * such that f(0) = 0, f(1) = 1, and 0 <= f(l) <= 1 for lambda 0-1.
     * The trigonometric switch can be restated thusly, since
     * cos^2(pi/2*lambda) is identical to sin^2(pi/2*(1-lambda)), f(1-l).
     */
    @Option(names = {"--sf", "--switchingFunction"}, paramLabel = "1.0",
            description = "Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)")
    String lambdaFunction = "1.0";

    /**
     * Return the switching function between topology energies.
     *
     * @return the switching function.
     */
    public UnivariateSwitchingFunction getSwitchingFunction() {
        UnivariateSwitchingFunction sf;
        if (!lambdaFunction.equalsIgnoreCase("1.0")) {
            String lf = lambdaFunction.toUpperCase();
            switch (lf) {
                // TODO implement the following RegEx case.
//                case ~/^-?[0-9]*\.?[0-9]+/:
//                    double exp = Double.parseDouble(lf);
//                    sf = new ffx.numerics.PowerSwitch(1.0, exp);
//                    break;
                case "TRIG":
                    sf = new SquaredTrigSwitch(false);
                    break;
                case "MULT":
                    sf = new MultiplicativeSwitch(0.0, 1.0);
                    break;
                default:
                    try {
                        double beta = Double.parseDouble(lf);
                        sf = new PowerSwitch(1.0, beta);
                    } catch (NumberFormatException ex) {
                        logger.warning(String.format("Argument to option -sf %s could not be properly parsed; using default linear switch",
                                lambdaFunction));
                        sf = new PowerSwitch(1.0, 1.0);
                    }
            }
        } else {
            sf = new PowerSwitch(1.0, 1.0);
        }
        return sf;
    }

    /**
     * The number of topologies to run in parallel.
     *
     * @param threadsAvail a int.
     * @param nArgs        a int.
     * @return a int.
     */
    public int getNumParallel(int threadsAvail, int nArgs) {
        int numParallel = nPar;
        if (threadsAvail % numParallel != 0) {
            logger.warning(String.format(" Number of threads available %d not evenly divisible by np %d reverting to sequential",
                    threadsAvail, numParallel));
            numParallel = 1;
        } else if (nArgs % numParallel != 0) {
            logger.warning(String.format(" Number of topologies %d not evenly divisible by np %d reverting to sequential",
                    nArgs, numParallel));
            numParallel = 1;
        }
        return numParallel;
    }

    /**
     * Set the alchemical atoms for this topology.
     *
     * @param topology a {@link ffx.potential.MolecularAssembly} object.
     */
    public void setSecondSystemAlchemistry(MolecularAssembly topology) {
        AlchemicalOptions.setAlchemicalAtoms(topology, s2, f2, ligAt2);
    }

    /**
     * Set uncharged atoms for this topology.
     *
     * @param topology a {@link ffx.potential.MolecularAssembly} object.
     */
    public void setSecondSystemUnchargedAtoms(MolecularAssembly topology) {
        AlchemicalOptions.setUnchargedAtoms(topology, es2, ef2);
    }

    /**
     * Collect unique atoms for the A dual-topology.
     *
     * @param topology A MolecularAssembly from dual-topology A.
     * @return A List of Integers.
     */
    public List<Integer> getUniqueAtomsA(MolecularAssembly topology) {
        return getUniqueAtoms(topology, "A", unsharedA);
    }

    /**
     * Collect unique atoms for the B dual-topology.
     *
     * @param topology A MolecularAssembly from dual-topology B.
     * @return A List of Integers.
     */
    public List<Integer> getUniqueAtomsB(MolecularAssembly topology) {
        return getUniqueAtoms(topology, "B", unsharedB);
    }

    /**
     * Collect unique atoms for a dual-topology. List MUST be sorted at the end.
     *
     * @param assembly A MolecularAssembly from the dual topology.
     * @param label    Either 'A' or 'B'.
     * @param unshared Atoms this dual topology isn't sharing.
     * @return A sorted List of Integers.
     */
    public List<Integer> getUniqueAtoms(MolecularAssembly assembly, String label, String unshared) {
        if (unshared != null && !unshared.isEmpty()) {
            logger.info(" Finding unique atoms for dual topology " + label);
            Set<Integer> indices = new HashSet<>();
            String[] toks = unshared.split("\\.");
            Atom[] atoms1 = assembly.getAtomArray();
            for (String range : toks) {
                Matcher m = AlchemicalOptions.rangeregex.matcher(range);
                if (m.find()) {
                    int rangeStart = Integer.parseInt(m.group(1));
                    int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                    if (rangeStart > rangeEnd) {
                        logger.severe(String.format(" Range %s was invalid start was greater than end", range));
                    }
                    logger.info(String.format(" Range %s for %s, start %d end %d", range, label, rangeStart, rangeEnd));
                    logger.fine(String.format(" First atom in range: %s", atoms1[rangeStart - 1]));
                    if (rangeEnd > rangeStart) {
                        logger.fine(String.format(" Last atom in range: %s", atoms1[rangeEnd - 1]));
                    }
                    for (int i = rangeStart; i <= rangeEnd; i++) {
                        indices.add(i - 1);
                    }
                }
            }
            int counter = 0;
            Set<Integer> adjustedIndices = new HashSet<>(); // Indexed by common variables in dtA.
            for (int i = 0; i < atoms1.length; i++) {
                Atom ai = atoms1[i];
                if (indices.contains(i)) {
                    if (ai.applyLambda()) {
                        logger.warning(String.format(
                                " Ranges defined in %s should not overlap with ligand atoms they are assumed to not be shared.", label));
                    } else {
                        logger.fine(String.format(" Unshared %s: %d variables %d-%d", label, i, counter, counter + 2));
                        for (int j = 0; j < 3; j++) {
                            adjustedIndices.add(counter + j);
                        }
                    }
                }
                if (!ai.applyLambda()) {
                    counter += 3;
                }
            }
            return adjustedIndices.stream().sorted().collect(Collectors.toList());
        } else {
            return Collections.emptyList();
        }
    }

    /**
     * Performs the bulk of the work of setting up a multi-topology system.
     * <p>
     * The sb StringBuilder is often something like "Timing energy and gradients for". The
     * method will append the exact type of Potential being assembled.
     *
     * @param assemblies   Opened MolecularAssembly(s).
     * @param threadsAvail Number of available threads.
     * @param sb           A StringBuilder describing what is to be done.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential assemblePotential(MolecularAssembly[] assemblies, int threadsAvail, StringBuilder sb) {
        int nargs = assemblies.length;
        int numPar = getNumParallel(threadsAvail, nargs);
        UnivariateSwitchingFunction sf = nargs > 1 ? getSwitchingFunction() : null;
        List<Integer> uniqueA;
        List<Integer> uniqueB;
        if (assemblies.length >= 4) {
            uniqueA = getUniqueAtomsA(assemblies[0]);
            uniqueB = getUniqueAtomsB(assemblies[2]);
        } else {
            uniqueA = Collections.emptyList();
            uniqueB = Collections.emptyList();
        }

        return getTopology(assemblies, sf, uniqueA, uniqueB, numPar, sb);
    }

    /**
     * Configure a Dual-, Quad- or Oct- Topology.
     *
     * @param topologies  The topologies.
     * @param sf          The Potential switching function.
     * @param uniqueA     The unique atoms of topology A.
     * @param uniqueB     The unique atoms of topology B.
     * @param numParallel The number of energies to evaluate in parallel.
     * @param sb          A StringBuilder for logging.
     * @return The Potential for the Topology.
     */
    public Potential getTopology(MolecularAssembly[] topologies,
                                 UnivariateSwitchingFunction sf,
                                 List<Integer> uniqueA, List<Integer> uniqueB,
                                 int numParallel, StringBuilder sb) {
        Potential potential = null;

        switch (topologies.length) {
            case 1:
                sb.append("single topology ");
                potential = topologies[0].getPotentialEnergy();
                break;
            case 2:
                sb.append("dual topology ");
                ffx.potential.DualTopologyEnergy dte = new ffx.potential.DualTopologyEnergy(topologies[0], topologies[1], sf);
                if (numParallel == 2) {
                    dte.setParallel(true);
                }
                potential = dte;
                break;
            case 4:
                sb.append("quad topology ");
                ffx.potential.DualTopologyEnergy dta = new ffx.potential.DualTopologyEnergy(topologies[0], topologies[1], sf);
                ffx.potential.DualTopologyEnergy dtb = new ffx.potential.DualTopologyEnergy(topologies[3], topologies[2], sf);
                ffx.potential.QuadTopologyEnergy qte = new ffx.potential.QuadTopologyEnergy(dta, dtb, uniqueA, uniqueB);
                if (numParallel >= 2) {
                    qte.setParallel(true);
                    if (numParallel == 4) {
                        dta.setParallel(true);
                        dtb.setParallel(true);
                    }
                }
                potential = qte;
                break;
            default:
                logger.severe(" Must have 2, 4, or 8 topologies!");
                break;
        }
        sb.append(Arrays.stream(topologies).
                map(MolecularAssembly::toString).
                collect(Collectors.joining(", ", " [", "] ")));
        return potential;
    }

    /**
     * If any softcore Atoms have been detected.
     *
     * @return Presence of softcore Atoms.
     */
    public boolean hasSoftcore() {
        return ((ligAt2 != null && ligAt2.length() > 0) || s2 > 0);
    }
}
