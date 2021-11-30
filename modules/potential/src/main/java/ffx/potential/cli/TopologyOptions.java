// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.potential.cli;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.numerics.switching.PowerSwitch;
import ffx.numerics.switching.SquaredTrigSwitch;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.Atom;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.stream.Collectors;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize multiple physical topologies.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class TopologyOptions {

  /** The logger for this class. */
  public static final Logger logger = Logger.getLogger(TopologyOptions.class.getName());

  /**
   * The ArgGroup keeps the TopologyOptions together when printing help.
   */
  @ArgGroup(heading = "%n Alchemical Options for Dual and Quad Topologies%n", validate = false)
  public TopologyOptionGroup group = new TopologyOptionGroup();

  /**
   * --ac2 or --alchemicalAtoms2 Specify alchemical atoms [ALL, NONE, Range(s): 1-3,6-N].
   *
   * @return Returns alchemical atoms for the 2nd topology.
   */
  public String getAlchemicalAtoms2() {
    return group.alchemicalAtoms2;
  }

  /**
   * --uc2 or --unchargedAtoms2 Specify atoms without electrostatics [ALL, NONE, Range(s): 1-3,6-N].
   *
   * @return Returns atoms without electrostatics for the 2nd topology.
   */
  public String getUnchargedAtoms2() {
    return group.unchargedAtoms2;
  }

  /**
   * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, or 4.
   *
   * @return Returns Number of topologies to evaluate in parallel.
   */
  public int getNPar() {
    return group.nPar;
  }

  /**
   * --uaA or --unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated
   * ranges or singletons.
   *
   * @return Return atoms unique to the A dual-topology.
   */
  public String getUnsharedA() {
    return group.unsharedA;
  }

  /**
   * --uaB or --unsharedB sets atoms unique to the A dual-topology, as period-separated hyphenated
   * ranges or singletons.
   *
   * @return Return atoms unique to the B dual-topology.
   */
  public String getUnsharedB() {
    return group.unsharedB;
  }

  /**
   * -sf or --switchingFunction
   *
   * <p>Sets the switching function to be used by dual topologies.
   *
   * <p>
   *
   * <ul>
   *   TRIG produces the function sin^2(pi/2*lambda)*E1(lambda) + cos^2(pi/2*lambda)*E2(1-lambda)
   * </ul>
   *
   * <ul>
   *   MULT uses a 5th-order polynomial switching function with zero first and second derivatives at
   *   the end (same function as used for van der Waals switch)
   * </ul>
   *
   * <ul>
   *   A number uses the original function, of l^beta*E1(lambda) + (1-lambda)^beta*E2(1-lambda).
   * </ul>
   *
   * <p>All of these are generalizations of <code>Udt = f(l)*E1(l) + f(1-l)*E2(1-lambda)</code>,
   * where f(l) is a continuous switching function such that f(0) = 0, f(1) = 1, and 0 <= f(l) <= 1
   * for lambda 0-1. The trigonometric switch can be restated thusly, since cos^2(pi/2*lambda) is
   * identical to sin^2(pi/2*(1-lambda)), f(1-l).
   *
   * @return Returns the lambda function.
   */
  public String getLambdaFunction() {
    return group.lambdaFunction;
  }

  /**
   * The number of topologies to run in parallel.
   *
   * @param threadsAvail a int.
   * @param nArgs a int.
   * @return a int.
   */
  public int getNumParallel(int threadsAvail, int nArgs) {
    int numParallel = getNPar();
    if (threadsAvail % numParallel != 0) {
      logger.warning(
          format(
              " Number of threads available %d not evenly divisible by np %d reverting to sequential",
              threadsAvail, numParallel));
      numParallel = 1;
    } else if (nArgs % numParallel != 0) {
      logger.warning(
          format(
              " Number of topologies %d not evenly divisible by np %d reverting to sequential",
              nArgs, numParallel));
      numParallel = 1;
    }
    return numParallel;
  }

  /**
   * Performs the bulk of the work of setting up a multi-topology system.
   *
   * <p>The sb StringBuilder is often something like "Timing energy and gradients for". The method
   * will append the exact type of Potential being assembled.
   *
   * @param assemblies Opened MolecularAssembly(s).
   * @param threadsAvail Number of available threads.
   * @param sb A StringBuilder describing what is to be done.
   * @return a {@link ffx.numerics.Potential} object.
   */
  public Potential assemblePotential(
      MolecularAssembly[] assemblies, int threadsAvail, StringBuilder sb) {
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
   * Return the switching function between topology energies.
   *
   * @return the switching function.
   */
  public UnivariateSwitchingFunction getSwitchingFunction() {
    UnivariateSwitchingFunction sf;
    String lambdaFunction = getLambdaFunction();
    if (!lambdaFunction.equalsIgnoreCase("1.0")) {
      String lf = lambdaFunction.toUpperCase();
      switch (lf) {
        case "TRIG":
          sf = new SquaredTrigSwitch(false);
          break;
        case "MULT":
          sf = new MultiplicativeSwitch(1.0, 0.0);
          break;
        default:
          try {
            double beta = parseDouble(lf);
            sf = new PowerSwitch(1.0, beta);
          } catch (NumberFormatException ex) {
            logger.warning(
                format(
                    "Argument to option -sf %s could not be properly parsed; using default linear switch",
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
   * Collect unique atoms for the A dual-topology.
   *
   * @param topology A MolecularAssembly from dual-topology A.
   * @return A List of Integers.
   */
  public List<Integer> getUniqueAtomsA(MolecularAssembly topology) {
    return getUniqueAtoms(topology, "A", getUnsharedA());
  }

  /**
   * Configure a Dual-, Quad- or Oct- Topology.
   *
   * @param topologies The topologies.
   * @param sf The Potential switching function.
   * @param uniqueA The unique atoms of topology A.
   * @param uniqueB The unique atoms of topology B.
   * @param numParallel The number of energies to evaluate in parallel.
   * @param sb A StringBuilder for logging.
   * @return The Potential for the Topology.
   */
  public Potential getTopology(
      MolecularAssembly[] topologies,
      UnivariateSwitchingFunction sf,
      List<Integer> uniqueA,
      List<Integer> uniqueB,
      int numParallel,
      StringBuilder sb) {
    Potential potential = null;
    switch (topologies.length) {
      case 1:
        sb.append("Single Topology ");
        potential = topologies[0].getPotentialEnergy();
        break;
      case 2:
        sb.append("Dual Topology ");
        DualTopologyEnergy dte = new DualTopologyEnergy(topologies[0], topologies[1], sf);
        if (numParallel == 2) {
          dte.setParallel(true);
        }
        potential = dte;
        break;
      case 4:
        sb.append("Quad Topology ");
        DualTopologyEnergy dta = new DualTopologyEnergy(topologies[0], topologies[1], sf);
        DualTopologyEnergy dtb = new DualTopologyEnergy(topologies[3], topologies[2], sf);
        QuadTopologyEnergy qte = new QuadTopologyEnergy(dta, dtb, uniqueA, uniqueB);
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
    sb.append(
        Arrays.stream(topologies)
            .map(MolecularAssembly::toString)
            .collect(Collectors.joining(", ", " [", "] ")));
    return potential;
  }

  /**
   * Collect unique atoms for a dual-topology. List MUST be sorted at the end.
   *
   * @param assembly A MolecularAssembly from the dual topology.
   * @param label Either 'A' or 'B'.
   * @param unshared Atoms this dual topology isn't sharing.
   * @return A sorted List of Integers.
   */
  public List<Integer> getUniqueAtoms(MolecularAssembly assembly, String label, String unshared) {
    if (!unshared.isEmpty()) {
      logger.info(" Finding unique atoms for dual topology " + label);
      Set<Integer> indices = new HashSet<>();
      String[] toks = unshared.split("\\.");
      Atom[] atoms1 = assembly.getAtomArray();
      for (String range : toks) {
        Matcher m = AlchemicalOptions.rangeRegEx.matcher(range);
        if (m.find()) {
          int rangeStart = Integer.parseInt(m.group(1));
          int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
          if (rangeStart > rangeEnd) {
            logger.severe(format(" Range %s was invalid start was greater than end", range));
          }
          logger.info(
              format(" Range %s for %s, start %d end %d", range, label, rangeStart, rangeEnd));
          logger.fine(format(" First atom in range: %s", atoms1[rangeStart - 1]));
          if (rangeEnd > rangeStart) {
            logger.fine(format(" Last atom in range: %s", atoms1[rangeEnd - 1]));
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
            logger.warning(
                format(
                    " Ranges defined in %s should not overlap with ligand atoms they are assumed to not be shared.",
                    label));
          } else {
            logger.fine(format(" Unshared %s: %d variables %d-%d", label, i, counter, counter + 2));
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
   * Collect unique atoms for the B dual-topology.
   *
   * @param topology A MolecularAssembly from dual-topology B.
   * @return A List of Integers.
   */
  public List<Integer> getUniqueAtomsB(MolecularAssembly topology) {
    return getUniqueAtoms(topology, "B", getUnsharedB());
  }

  /**
   * If any softcore Atoms have been detected.
   *
   * @return Presence of softcore Atoms.
   */
  public boolean hasSoftcore() {
    String alchemicalAtoms2 = getAlchemicalAtoms2();
    return (alchemicalAtoms2 != null
        && !alchemicalAtoms2.equalsIgnoreCase("NONE")
        && !alchemicalAtoms2.equalsIgnoreCase(""));
  }

  /**
   * Set the alchemical atoms for this topology.
   *
   * @param topology a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setSecondSystemAlchemistry(MolecularAssembly topology) {
    String alchemicalAtoms2 = getAlchemicalAtoms2();
    AlchemicalOptions.setAlchemicalAtoms(topology, alchemicalAtoms2);
  }

  /**
   * Set uncharged atoms for this topology.
   *
   * @param topology a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setSecondSystemUnchargedAtoms(MolecularAssembly topology) {
    String unchargedAtoms2 = getUnchargedAtoms2();
    AlchemicalOptions.setUnchargedAtoms(topology, unchargedAtoms2);
  }

  /**
   * Collection of Topology Options.
   */
  private static class TopologyOptionGroup {

    /** --ac2 or --alchemicalAtoms2 Specify alchemical atoms [ALL, NONE, Range(s): 1-3,6-N]. */
    @Option(
        names = {"--ac2", "--alchemicalAtoms2"},
        paramLabel = "<selection>",
        defaultValue = "",
        description = "Specify alchemical atoms for the 2nd topology [ALL, NONE, Range(s): 1-3,6-N].")
    String alchemicalAtoms2;

    /**
     * --uc2 or --unchargedAtoms2 Specify atoms without electrostatics [ALL, NONE, Range(s):
     * 1-3,6-N]."
     */
    @Option(
        names = {"--uc2", "--unchargedAtoms2"},
        paramLabel = "<selection>",
        defaultValue = "",
        description =
            "Specify atoms without electrostatics for the 2nd topology [ALL, NONE, Range(s): 1-3,6-N].")
    String unchargedAtoms2;

    /**
     * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, or
     * 4.
     */
    @Option(
        names = {"--np", "--nParallel"},
        paramLabel = "1",
        description = "Number of topologies to evaluate in parallel")
    int nPar = 1;

    /**
     * --uaA or --unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated
     * ranges or singletons.
     */
    @Option(
        names = {"--uaA", "--unsharedA"},
        paramLabel = "-1",
        description = "Unshared atoms in the A dual topology (e.g. 1-24.32-65).")
    String unsharedA = null;

    /**
     * --uaB or --unsharedB sets atoms unique to the B dual-topology, as period-separated hyphenated
     * ranges or singletons.
     */
    @Option(
        names = {"--uaB", "--unsharedB"},
        paramLabel = "-1",
        description = "Unshared atoms in the B dual topology (e.g. 1-24.32-65).")
    String unsharedB = null;

    /**
     * -sf or --switchingFunction
     *
     * <p>Sets the switching function to be used by dual topologies.
     *
     * <p>
     *
     * <ul>
     *   TRIG produces the function sin^2(pi/2*lambda)*E1(lambda) + cos^2(pi/2*lambda)*E2(1-lambda)
     * </ul>
     *
     * <ul>
     *   MULT uses a 5th-order polynomial switching function with zero first and second derivatives at
     *   the end (same function as used for van der Waals switch)
     * </ul>
     *
     * <ul>
     *   A number uses the original function, of l^beta*E1(lambda) + (1-lambda)^beta*E2(1-lambda).
     * </ul>
     *
     * <p>All of these are generalizations of <code>Udt = f(l)*E1(l) + f(1-l)*E2(1-lambda)</code>,
     * where f(l) is a continuous switching function such that f(0) = 0, f(1) = 1, and 0 <= f(l) <= 1
     * for lambda 0-1. The trigonometric switch can be restated thusly, since cos^2(pi/2*lambda) is
     * identical to sin^2(pi/2*(1-lambda)), f(1-l).
     */
    @Option(
        names = {"--sf", "--switchingFunction"},
        paramLabel = "1.0",
        description =
            "Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)")
    String lambdaFunction = "1.0";
  }
}
