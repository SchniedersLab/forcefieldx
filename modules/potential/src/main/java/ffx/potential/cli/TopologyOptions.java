/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential.cli;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ffx.numerics.Potential;
import ffx.numerics.PowerSwitch;
import ffx.numerics.SquaredTrigSwitch;
import ffx.numerics.UnivariateSwitchingFunction;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.MultiplicativeSwitch;

import picocli.CommandLine.Option;

/**
 * Alchemical options shared by scripts that use Dual or Quad Topologies and use the Pico CLI.
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
     * @param threadsAvail
     * @param nArgs
     * @return
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
     * @param topology
     */
    public void setAlchemicalAtoms(MolecularAssembly topology) {
        Atom atoms[] = topology.getAtomArray();
        if (s2 > 0) {
            for (int i = s2; i <= f2; i++) {
                Atom ai = atoms[i - 1];
                ai.setApplyLambda(true);
                ai.print();
            }
        }

        if (ligAt2 != null) {
            String ranges[] = ligAt2.split(".");
            Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");
            for (String range : ranges) {
                Matcher m = rangeregex.matcher(range);
                if (m.find()) {
                    int rangeStart = Integer.parseInt(m.group(1));
                    int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                    if (rangeStart > rangeEnd) {
                        logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                    }
                    // Don't need to worry about negative numbers; rangeregex just won't match.
                    for (int i = rangeStart; i <= rangeEnd; i++) {
                        Atom ai = atoms[i - 1];
                        ai.setApplyLambda(true);
                        ai.print();
                    }
                } else {
                    logger.warning(" Could not recognize ${range} as a valid range; skipping");
                }
            }
        }
    }

    /**
     * Set uncharged atoms for this topology.
     *
     * @param topology
     */
    public void setUnchargedAtoms(MolecularAssembly topology) {
        Atom atoms[] = topology.getAtomArray();
        // Apply the no electrostatics atom selection
        int noElecStart = es2;
        noElecStart = (noElecStart < 1) ? 1 : noElecStart;

        int noElecStop = ef2;
        noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop;

        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1];
            ai.setElectrostatics(false);
            ai.print();
        }
    }

    /**
     * Collect unique atoms for a topology.
     *
     * @param topology The MolecularAsembly.
     * @param label    A label for logging (i.e. A or B).
     * @return A List of Integers.
     */
    public List<Integer> getUniqueAtoms(MolecularAssembly topology, String label) {
        List<Integer> unique = new ArrayList<>();
        Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");

        if (unsharedA != null) {
            Set<Integer> ra = new HashSet<>();
            String[] toksA = unsharedA.split(".");
            Atom[] atA1 = topology.getAtomArray();
            for (String range : toksA) {
                Matcher m = rangeregex.matcher(range);
                if (m.find()) {
                    int rangeStart = Integer.parseInt(m.group(1));
                    int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                    if (rangeStart > rangeEnd) {
                        logger.severe(String.format(" Range %s was invalid start was greater than end", range));
                    }
                    logger.info(String.format(" Range %s for %s, start %d end %d", range, label, rangeStart, rangeEnd));
                    logger.info(String.format(" First atom in range: %s", atA1[rangeStart - 1]));
                    if (rangeEnd > rangeStart) {
                        logger.info(String.format(" Last atom in range: %s", atA1[rangeEnd - 1]));
                    }
                    for (int i = rangeStart; i <= rangeEnd; i++) {
                        ra.add(i - 1);
                    }
                }
            }
            int counter = 0;
            Set<Integer> raAdj = new HashSet<>(); // Indexed by common variables in dtA.
            for (int i = 0; i < atA1.length; i++) {
                Atom ai = atA1[i];
                if (ra.contains(i)) {
                    if (ai.applyLambda()) {
                        logger.warning(String.format(
                                " Ranges defined in %s should not overlap with ligand atoms they are assumed to not be shared.", label));
                    } else {
                        logger.fine(String.format(" Unshared %s: %d variables %d-%d", label, i, counter, counter + 2));
                        for (int j = 0; j < 3; j++) {
                            raAdj.add(Integer.valueOf(counter + j));
                        }
                    }
                }
                if (!ai.applyLambda()) {
                    counter += 3;
                }
            }
            unique.addAll(raAdj);
        }

        return unique;
    }

    /**
     * Configure a Dual-, Quad- or Oct- Topology.
     *
     * @param num         The number of topologies.
     * @param topologies  The topologies.
     * @param sf          The Potential switching function.
     * @param uniqueA     The unique atoms of topology A.
     * @param uniqueB     The unique atoms of topology B.
     * @param numParallel The number of energies to evaluate in paralle.
     * @param sb          A StringBuilder for logging.
     * @return The Potential for the Topology.
     */
    public Potential getTopology(int num, MolecularAssembly[] topologies,
                                 UnivariateSwitchingFunction sf,
                                 List<Integer> uniqueA, List<Integer> uniqueB,
                                 int numParallel, StringBuilder sb) {
        Potential potential = null;

        switch (num) {
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
            case 8:
                sb.append("oct-topology ");
                ffx.potential.DualTopologyEnergy dtga = new ffx.potential.DualTopologyEnergy(topologies[0], topologies[1], sf);
                ffx.potential.DualTopologyEnergy dtgb = new ffx.potential.DualTopologyEnergy(topologies[3], topologies[2], sf);
                ffx.potential.QuadTopologyEnergy qtg = new ffx.potential.QuadTopologyEnergy(dtga, dtgb, uniqueA, uniqueB);
                ffx.potential.DualTopologyEnergy dtda = new ffx.potential.DualTopologyEnergy(topologies[4], topologies[5], sf);
                ffx.potential.DualTopologyEnergy dtdb = new ffx.potential.DualTopologyEnergy(topologies[7], topologies[6], sf);
                ffx.potential.QuadTopologyEnergy qtd = new ffx.potential.QuadTopologyEnergy(dtda, dtdb, uniqueA, uniqueB);
                ffx.potential.OctTopologyEnergy ote = new ffx.potential.OctTopologyEnergy(qtg, qtd, true);
                if (numParallel >= 2) {
                    ote.setParallel(true);
                    if (numParallel >= 4) {
                        qtg.setParallel(true);
                        qtd.setParallel(true);
                        if (numParallel == 8) {
                            dtga.setParallel(true);
                            dtgb.setParallel(true);
                            dtda.setParallel(true);
                            dtdb.setParallel(true);
                        }
                    }
                }
                potential = ote;
                break;
            default:
                logger.severe(" Must have 2, 4, or 8 topologies!");
                break;
        }
        return potential;
    }

}
