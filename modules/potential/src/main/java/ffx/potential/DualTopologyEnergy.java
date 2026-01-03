// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential;

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.openmm.OpenMMDualTopologyEnergy;
import ffx.potential.openmm.OpenMMEnergy;
import ffx.potential.parameters.ForceField;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.crystal.SymOp.applyCartesianSymOp;
import static ffx.crystal.SymOp.applyCartesianSymRot;
import static ffx.crystal.SymOp.invertSymOp;
import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.utilities.StringUtils.parseAtomRanges;
import static ffx.utilities.StringUtils.writeAtomRanges;
import static java.lang.Double.parseDouble;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Compute the potential energy and derivatives for a dual-topology system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DualTopologyEnergy implements CrystalPotential, LambdaInterface {

  /**
   * Logger for the DualTopologyEnergy class.
   */
  private static final Logger logger = Logger.getLogger(DualTopologyEnergy.class.getName());
  /**
   * Shared atoms between topologies 1 and 2.
   */
  private final int nShared;
  /**
   * Total number of variables
   */
  private final int nVariables;
  /**
   * Region for computing energy/gradient in parallel.
   */
  private final EnergyRegion energyRegion;
  /**
   * Mass array for shared and softcore atoms.
   */
  private final double[] mass;
  /**
   * VARIABLE_TYPE array for shared and softcore atoms.
   */
  private final VARIABLE_TYPE[] variableTypes;
  /**
   * Coordinates for topology 1.
   */
  private final double[] x1;
  /**
   * Coordinates for topology 2.
   */
  private final double[] x2;
  /**
   * Topology 1 coordinate gradient.
   */
  private final double[] g1;
  /**
   * Topology 2 coordinate gradient.
   */
  private final double[] g2;
  /**
   * Topology 1 restraint gradient for end state bonded terms.
   */
  private final double[] rg1;
  /**
   * Topology 2 restraint gradient for end state bonded terms.
   */
  private final double[] rg2;
  /**
   * Topology 1 derivative of the coordinate gradient with respect to lambda.
   */
  private final double[] gl1;
  /**
   * Topology 2 derivative of the coordinate gradient with respect to lambda.
   */
  private final double[] gl2;
  /**
   * Topology 1 derivative of the coordinate gradient with respect to lambda for end state bonded
   * terms
   */
  private final double[] rgl1;
  /**
   * Topology 2 derivative of the coordinate gradient with respect to lambda for end state bonded
   * terms
   */
  private final double[] rgl2;
  /**
   * Topology 1 Potential.
   */
  private final CrystalPotential potential1;
  /**
   * Topology 2 Potential.
   */
  private final CrystalPotential potential2;
  /**
   * Topology 1 LambdaInterface.
   */
  private final LambdaInterface lambdaInterface1;
  /**
   * Topology 2 LambdaInterface.
   */
  private final LambdaInterface lambdaInterface2;
  /**
   * Topology 1 ForceFieldEnergy.
   */
  private final ForceFieldEnergy forceFieldEnergy1;
  /**
   * Topology 2 ForceFieldEnergy.
   */
  private final ForceFieldEnergy forceFieldEnergy2;
  /**
   * Number of active atoms in Topology 1.
   */
  private final int nActive1;
  /**
   * Number of active atoms in Topology 2.
   */
  private final int nActive2;
  /**
   * Array of active atoms in Topology 1.
   */
  private final Atom[] activeAtoms1;
  /**
   * Array of active atoms in Topology 2.
   */
  private final Atom[] activeAtoms2;
  /**
   * Flag is true if the atom is shared (not alchemical) for Topology 1.
   */
  private final boolean[] sharedAtoms1;
  /**
   * Flag is true if the atom is shared (not alchemical) for Topology 2.
   */
  private final boolean[] sharedAtoms2;
  /**
   * Index for each atom in topology 1 into the overall dual-topology atom array.
   * This includes all atoms, including (in)active, shared atoms and alchemical atoms.
   */
  private final int[] system1AtomIndex;
  /**
   * Index for each atom in topology 2 into the overall dual-topology atom array.
   * This includes all atoms, including (in)active, shared atoms and alchemical atoms.
   */
  private final int[] system2AtomIndex;
  /**
   * An array of shared and alchemical atoms for the dual topology system.
   * <p>
   * The shared atoms from Topology 1.
   */
  private final Atom[] dualTopologyAtoms;
  /**
   * An array of shared and alchemical atoms for the dual topology system.
   * <p>
   * The shared atoms from Topology 2.
   */
  private final Atom[] dualTopologyAtoms2;
  /**
   * Will default to a power-1 PowerSwitch function.
   */
  private final UnivariateSwitchingFunction switchFunction;
  /**
   * Current potential energy of topology 1 (kcal/mol).
   */
  private double energy1 = 0;
  /**
   * Current potential energy of topology 2 (kcal/mol).
   */
  private double energy2 = 0;
  /**
   * ParallelTeam to execute the EnergyRegion.
   */
  private ParallelTeam parallelTeam;
  /**
   * Include a valence restraint energy for atoms being "disappeared."
   */
  private final boolean doValenceRestraint1;
  /**
   * Include a valence restraint energy for atoms being "disappeared."
   */
  private final boolean doValenceRestraint2;
  /**
   * Use System 1 Bonded Energy.
   */
  private final boolean useFirstSystemBondedEnergy;
  /**
   * End-state restraint energy of topology 1 (kcal/mol).
   */
  private double restraintEnergy1 = 0;
  /**
   * End-state restraint energy of topology 2 (kcal/mol).
   */
  private double restraintEnergy2 = 0;
  /**
   * dEdL of topology 1 (kcal/mol).
   */
  private double dEdL_1 = 0;
  /**
   * dEdL of topology 2 (kcal/mol).
   */
  private double dEdL_2 = 0;
  /**
   * End-state restraint dEdL of topology 1 (kcal/mol).
   */
  private double restraintdEdL_1 = 0;
  /**
   * End-state restraint dEdL of topology 2 (kcal/mol).
   */
  private double restraintdEdL_2 = 0;
  /**
   * d2EdL2 of topology 1 (kcal/mol).
   */
  private double d2EdL2_1 = 0;
  /**
   * d2EdL2 of topology 2 (kcal/mol).
   */
  private double d2EdL2_2 = 0;
  /**
   * End-state restraint d2EdL2 of topology 1 (kcal/mol).
   */
  private double restraintd2EdL2_1 = 0;
  /**
   * End-state restraint d2EdL2 of topology 2 (kcal/mol).
   */
  private double restraintd2EdL2_2 = 0;
  /**
   * Total energy of the dual topology, including lambda scaling.
   */
  private double totalEnergy = 0;
  /**
   * Current lambda value.
   */
  private double lambda = 1.0;
  /**
   * Value of the switching function at lambda.
   */
  private double f1L = 1.0;
  /**
   * Value of the switching function at (1-lambda).
   */
  private double f2L = 0.0;
  /**
   * First derivative with respect to lambda of the switching function.
   */
  private double dF1dL = 0.0;
  /**
   * First derivative with respect to (1-lambda) of the switching function.
   */
  private double dF2dL = 0.0;
  /**
   * Second derivative with respect to lambda of the switching function.
   */
  private double d2F1dL2 = 0.0;
  /**
   * Second derivative with respect to (1-lambda) of the switching function.
   */
  private double d2F2dL2 = 0.0;
  /**
   * Scaling array for shared and softcore atoms.
   */
  private double[] scaling = null;
  /**
   * State for calculating energies.
   */
  private STATE state = STATE.BOTH;
  /**
   * Symmetry operator to superimpose system 2 onto system 1
   */
  private SymOp[] symOp;
  /**
   * Symmetry operator to move system 2 back to original frame
   */
  private SymOp[] inverse;
  /**
   * Atom selection for each symmetry operator.
   */
  private final ArrayList<List<Integer>> symOpAtoms = new ArrayList<>();
  /**
   * Mask to determine which atoms a sym op will be applied.
   */
  private boolean[][] mask = null;
  /**
   * Utilize a provided SymOp
   */
  private boolean useSymOp = false;

  /**
   * Constructor for DualTopologyEnergy.
   *
   * @param topology1      a {@link ffx.potential.MolecularAssembly} object.
   * @param topology2      a {@link ffx.potential.MolecularAssembly} object.
   * @param switchFunction a {@link UnivariateSwitchingFunction} object.
   */
  public DualTopologyEnergy(
      MolecularAssembly topology1,
      MolecularAssembly topology2,
      UnivariateSwitchingFunction switchFunction) {
    forceFieldEnergy1 = topology1.getPotentialEnergy();
    forceFieldEnergy2 = topology2.getPotentialEnergy();
    potential1 = forceFieldEnergy1;
    potential2 = forceFieldEnergy2;
    lambdaInterface1 = forceFieldEnergy1;
    lambdaInterface2 = forceFieldEnergy2;

    /* Atom array for topology 1. */
    Atom[] atoms1 = topology1.getAtomArray();
    /* Atom array for topology 2. */
    Atom[] atoms2 = topology2.getAtomArray();

    ForceField forceField1 = topology1.getForceField();
    doValenceRestraint1 = forceField1.getBoolean("LAMBDA_VALENCE_RESTRAINTS", true);
    ForceField forceField2 = topology2.getForceField();
    doValenceRestraint2 = forceField2.getBoolean("LAMBDA_VALENCE_RESTRAINTS", true);

    useFirstSystemBondedEnergy = forceField2.getBoolean("USE_FIRST_SYSTEM_BONDED_ENERGY", false);

    // Check that all atoms that are not undergoing alchemy are common to both topologies.
    int shared1 = 0;
    int shared2 = 0;
    int activeCount1 = 0;
    int activeCount2 = 0;
    for (Atom a1 : atoms1) {
      if (a1.isActive()) {
        activeCount1++;
        if (!a1.applyLambda()) {
          shared1++;
        }
      }
    }
    for (Atom a2 : atoms2) {
      if (a2.isActive()) {
        activeCount2++;
        if (!a2.applyLambda()) {
          shared2++;
        }
      }
    }

    if (shared1 != shared2) {
      logger.severe(" Shared active atoms are not equal between topologies:\n"
          + " Topology 1 has " + shared1 + " shared atoms.\n"
          + " Topology 2 has " + shared2 + " shared atoms.\n"
          + " Please check the input files or force field parameters.");
    }

    assert (shared1 == shared2);
    nActive1 = activeCount1;
    nActive2 = activeCount2;
    activeAtoms1 = new Atom[nActive1];
    activeAtoms2 = new Atom[nActive2];
    sharedAtoms1 = new boolean[nActive1];
    sharedAtoms2 = new boolean[nActive2];
    Arrays.fill(sharedAtoms1, true);
    Arrays.fill(sharedAtoms2, true);

    // Create a dual topology list that includes all atoms from both topologies (even inactive atoms).
    int nAlchemical1 = 0;
    int nShared1 = 0;
    int nAlchemical2 = 0;
    int nShared2 = 0;
    for (Atom a : atoms1) {
      if (a.applyLambda()) {
        nAlchemical1++;
      } else {
        nShared1++;
      }
    }
    for (Atom a : atoms2) {
      if (a.applyLambda()) {
        nAlchemical2++;
      } else {
        nShared2++;
      }
    }
    if (nShared1 != nShared2) {
      logger.severe(" Shared atoms are not equal between topologies:\n"
          + " Topology 1 has " + nShared1 + " shared atoms.\n"
          + " Topology 2 has " + nShared2 + " shared atoms.\n"
          + " Please check the input files or force field parameters.");
    }
    dualTopologyAtoms = new Atom[nShared1 + nAlchemical1 + nAlchemical2];
    dualTopologyAtoms2 = new Atom[nShared1 + nAlchemical1 + nAlchemical2];
    int indexShared = 0;
    int indexAlchemical = nShared1;
    int index1 = 0;
    system1AtomIndex = new int[atoms1.length];
    for (Atom a : atoms1) {
      a.setTopologyIndex(0);
      if (a.applyLambda()) {
        system1AtomIndex[index1] = indexAlchemical;
        a.setTopologyAtomIndex(indexAlchemical);
        dualTopologyAtoms[indexAlchemical] = a;
        dualTopologyAtoms2[indexAlchemical] = a;
        indexAlchemical++;
      } else {
        system1AtomIndex[index1] = indexShared;
        a.setTopologyAtomIndex(indexShared);
        dualTopologyAtoms[indexShared] = a;
        indexShared++;
      }
      index1++;
    }
    // Reset the shared index.
    indexShared = 0;
    int index2 = 0;
    system2AtomIndex = new int[atoms2.length];
    for (Atom a : atoms2) {
      a.setTopologyIndex(1);
      if (a.applyLambda()) {
        system2AtomIndex[index2] = indexAlchemical;
        a.setTopologyAtomIndex(indexAlchemical);
        dualTopologyAtoms[indexAlchemical] = a;
        dualTopologyAtoms2[indexAlchemical] = a;
        indexAlchemical++;
      } else {
        system2AtomIndex[index2] = indexShared;
        a.setTopologyAtomIndex(indexShared);
        dualTopologyAtoms2[indexShared] = a;
        indexShared++;
      }
      index2++;
    }

    // Fill the active atom arrays and shared atom flags.
    int index = 0;
    for (Atom a1 : atoms1) {
      if (a1.isActive()) {
        activeAtoms1[index] = a1;
        if (a1.applyLambda()) {
          // This atom is softcore with independent coordinates.
          sharedAtoms1[index] = false;
        }
        index++;
      }
    }
    index = 0;
    for (Atom a2 : atoms2) {
      if (a2.isActive()) {
        activeAtoms2[index] = a2;
        if (a2.applyLambda()) {
          // This atom is softcore with independent coordinates.
          sharedAtoms2[index] = false;
        }
        index++;
      }
    }

    nShared = shared1;
    /* Topology 1 number of softcore atoms. */
    int nSoftCore1 = nActive1 - nShared;
    /* Topology 2 number of softcore atoms. */
    int nSoftCore2 = nActive2 - nShared;
    /* Total number of softcore and shared atoms: nTotal = nShared + nSoftcore1 + nSoftcore2 */
    int nTotal = nShared + nSoftCore1 + nSoftCore2;
    nVariables = 3 * nTotal;

    // Allocate memory for coordinates and derivatives.
    x1 = new double[nActive1 * 3];
    x2 = new double[nActive2 * 3];
    g1 = new double[nActive1 * 3];
    g2 = new double[nActive2 * 3];
    rg1 = new double[nActive1 * 3];
    rg2 = new double[nActive2 * 3];
    gl1 = new double[nActive1 * 3];
    gl2 = new double[nActive2 * 3];
    rgl1 = new double[nActive1 * 3];
    rgl2 = new double[nActive2 * 3];

    // All variables are coordinates.
    index = 0;
    variableTypes = new VARIABLE_TYPE[nVariables];
    for (int i = 0; i < nTotal; i++) {
      variableTypes[index++] = VARIABLE_TYPE.X;
      variableTypes[index++] = VARIABLE_TYPE.Y;
      variableTypes[index++] = VARIABLE_TYPE.Z;
    }

    // Fill the mass array.
    int commonIndex = 0;
    int softcoreIndex = 3 * nShared;
    mass = new double[nVariables];
    for (int i = 0; i < nActive1; i++) {
      Atom a = activeAtoms1[i];
      double m = a.getMass();
      if (sharedAtoms1[i]) {
        mass[commonIndex++] = m;
        mass[commonIndex++] = m;
        mass[commonIndex++] = m;
      } else {
        mass[softcoreIndex++] = m;
        mass[softcoreIndex++] = m;
        mass[softcoreIndex++] = m;
      }
    }
    commonIndex = 0;
    for (int i = 0; i < nActive2; i++) {
      Atom a = activeAtoms2[i];
      double m = a.getMass();
      if (sharedAtoms2[i]) {
        for (int j = 0; j < 3; j++) {
          double massCommon = mass[commonIndex];
          massCommon = Math.max(massCommon, m);
          mass[commonIndex++] = massCommon;
        }
      } else {
        mass[softcoreIndex++] = m;
        mass[softcoreIndex++] = m;
        mass[softcoreIndex++] = m;
      }
    }
    energyRegion = new EnergyRegion();
    parallelTeam = new ParallelTeam(1);

    this.switchFunction = switchFunction;
    logger.info(format("\n Dual topology using switching function:\n  %s", switchFunction));
    logger.info(format(" Shared atoms: %d (1: %d of %d; 2: %d of %d)\n ", nShared, shared1, nActive1, shared2, nActive2));
    String symOpString = forceField2.getString("symop", null);
    if (symOpString != null) {
      readSymOp(symOpString);
    } else {
      // Both end-states will use the same crystal object.
      potential2.setCrystal(potential1.getCrystal());
      useSymOp = false;
      symOp = null;
      inverse = null;
    }

    // Check that all Dual-Topology atoms start with identical coordinates.
    int i1 = 0;
    int i2 = 0;
    // Not true if mapping atoms via symmetry operator.
    if (!useSymOp) {
      for (int i = 0; i < nShared; i++) {
        Atom a1 = atoms1[i1++];
        while (a1.applyLambda()) {
          a1 = atoms1[i1++];
        }
        Atom a2 = atoms2[i2++];
        while (a2.applyLambda()) {
          a2 = atoms2[i2++];
        }
        assert (a1.getX() == a2.getX());
        assert (a1.getY() == a2.getY());
        assert (a1.getZ() == a2.getZ());
        reconcileAtoms(a1, a2, Level.INFO);
      }
    }
  }

  /**
   * Static factory method to create a DualTopologyEnergy, possibly via FFX or OpenMM implementations.
   *
   * @param molecularAssembly1 MolecularAssembly for topology 1.
   * @param molecularAssembly2 MolecularAssembly for topology 2.
   * @param switchFunction     The UnivariateSwitchingFunction to use for the dual topology energy.
   * @return A DualTopologyEnergy on some Platform
   */
  @SuppressWarnings("fallthrough")
  public static DualTopologyEnergy energyFactory(MolecularAssembly molecularAssembly1,
                                                 MolecularAssembly molecularAssembly2,
                                                 UnivariateSwitchingFunction switchFunction) {
    ForceField forceField = molecularAssembly1.getForceField();
    String platformString = toEnumForm(forceField.getString("PLATFORM-DT", "FFX"));
    try {
      Platform platform = Platform.valueOf(platformString);
      switch (platform) {
        case OMM, OMM_REF, OMM_CUDA, OMM_OPENCL:
          try {
            return new OpenMMDualTopologyEnergy(molecularAssembly1, molecularAssembly2, switchFunction, platform);
          } catch (Exception ex) {
            logger.warning(format(" Exception creating OpenMMDualTopologyEnergy: %s", ex));
            return new DualTopologyEnergy(molecularAssembly1, molecularAssembly2, switchFunction);
          }
        case OMM_CPU:
          logger.warning(format(" Platform %s not supported; defaulting to FFX", platform));
        default:
          return new DualTopologyEnergy(molecularAssembly1, molecularAssembly2, switchFunction);
      }
    } catch (IllegalArgumentException | NullPointerException ex) {
      logger.warning(
          format(" String %s did not match a known energy implementation", platformString));
      return new DualTopologyEnergy(molecularAssembly1, molecularAssembly2, switchFunction);
    }
  }


  /**
   * Get the number of dual topology atoms.
   *
   * @return The number of dual topology atoms.
   */
  public int getNumberOfAtoms() {
    return dualTopologyAtoms.length;
  }

  /**
   * Get the dual topology atoms for the specified topology.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @return An array of Atoms for dual-topology with shared atoms from the specified topology.
   */
  public Atom[] getDualTopologyAtoms(int topology) {
    if (topology == 0) {
      return dualTopologyAtoms;
    } else if (topology == 1) {
      return dualTopologyAtoms2;
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Get the atom from the dual-topology atom array corresponding to the specified topology and index.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @param index    The index of the atom in the dual topology atom array.
   * @return The Atom from the dual-topology atom array corresponding to the specified topology and index.
   */
  public Atom getDualTopologyAtom(int topology, int index) {
    if (topology == 0) {
      return dualTopologyAtoms[index];
    } else if (topology == 1) {
      return dualTopologyAtoms2[index];
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Map an atomic index from a single topology to the overall dual-topology atom array.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @param index    The index of the atom in the specified topology.
   * @return The index of the atom in the overall dual-topology atom array.
   */
  public int mapToDualTopologyIndex(int topology, int index) {
    if (topology == 0) {
      return system1AtomIndex[index];
    } else if (topology == 1) {
      return system2AtomIndex[index];
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Get the ForceFieldEnergy for topology 1.
   *
   * @return The ForceFieldEnergy for topology 1.
   */
  public ForceFieldEnergy getForceFieldEnergy1() {
    return forceFieldEnergy1;
  }

  /**
   * Get the ForceFieldEnergy for topology 2.
   *
   * @return The ForceFieldEnergy for topology 2.
   */
  public ForceFieldEnergy getForceFieldEnergy2() {
    return forceFieldEnergy2;
  }

  /**
   * Get the atom index for topology 1 into the overall dual-topology atom array.
   *
   * @return An array of indices for topology 1 atoms into the dual-topology atom array.
   */
  public int[] getDualTopologyIndex1() {
    return system1AtomIndex;
  }

  /**
   * Get the atom index for topology 2 into the overall dual-topology atom array.
   *
   * @return An array of indices for topology 2 atoms into the dual-topology atom array.
   */
  public int[] getDualTopologyIndex2() {
    return system2AtomIndex;
  }

  /**
   * The scale factor for the specified topology.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @return The scale factor for the topology.
   */
  public double getTopologyScale(int topology) {
    if (topology == 0) {
      return f1L;
    } else if (topology == 1) {
      return f2L;
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Temporary method to be used as a benchmark. Move to SymOp?
   *
   * @param symOpString String containing sym op.
   */
  private void readSymOp(String symOpString) {
    int numSymOps;
    try {
      String[] tokens = symOpString.split(" +");
      int numTokens = tokens.length;
      // Only need to move one system onto the other.
      if (numTokens == 12) {
        numSymOps = 1;
        symOp = new SymOp[numSymOps];
        inverse = new SymOp[numSymOps];
        // Assume applies to all atoms...
        symOpAtoms.add(parseAtomRanges("symmetryAtoms", "1-N", nActive2));
        // Twelve values makes 1 sym op. so assign to [0][0].
        symOp[0] = new SymOp(new double[][]{
            {parseDouble(tokens[0]), parseDouble(tokens[1]), parseDouble(tokens[2])},
            {parseDouble(tokens[3]), parseDouble(tokens[4]), parseDouble(tokens[5])},
            {parseDouble(tokens[6]), parseDouble(tokens[7]), parseDouble(tokens[8])}},
            new double[]{parseDouble(tokens[9]), parseDouble(tokens[10]), parseDouble(tokens[11])});
      } else if (numTokens % 14 == 0) {
        numSymOps = numTokens / 14;
        symOp = new SymOp[numSymOps];
        inverse = new SymOp[numSymOps];
        if (logger.isLoggable(Level.FINER)) {
          logger.finer(format(" Number of Tokens: %3d Number of Sym Ops: %2d", numTokens, numSymOps));
        }
        for (int i = 0; i < numSymOps; i++) {
          int j = i * 14;
          // Should be tokens[j] or tokens[j+1]?? Switches forward vs backward...
          symOpAtoms.add(parseAtomRanges("symmetryAtoms", tokens[j], nActive2));
          symOp[i] = new SymOp(new double[][]{
              {parseDouble(tokens[j + 2]), parseDouble(tokens[j + 3]), parseDouble(tokens[j + 4])},
              {parseDouble(tokens[j + 5]), parseDouble(tokens[j + 6]), parseDouble(tokens[j + 7])},
              {parseDouble(tokens[j + 8]), parseDouble(tokens[j + 9]), parseDouble(tokens[j + 10])}},
              new double[]{parseDouble(tokens[j + 11]), parseDouble(tokens[j + 12]), parseDouble(tokens[j + 13])});
        }
      } else {
        logger.warning(format(" Provided symmetry operator is formatted incorrectly. Should have 12 or 14 entries, but found %3d.", numTokens));
        return;
      }
      useSymOp = true;
      // Define mask(s) needed to apply sym ops.
      mask = new boolean[numSymOps][nActive2];
      // Mask(s) used locally to check symmetry operator overlap.
      boolean[][] sharedMask = new boolean[numSymOps][nShared];
      for (int i = 0; i < numSymOps; i++) {
        List<Integer> symAtoms = symOpAtoms.get(i);
        for (int atom : symAtoms) {
          if (sharedAtoms2[atom]) {
            mask[i][atom] = true;
          }
        }
        // TODO: Condense with loop above by replacing nActive2 with incremental index?
        int sharedIndex = 0;
        for (int j = 0; j < nActive2; j++) {
          if (sharedAtoms2[j]) {
            if (mask[i][j]) {
              sharedMask[i][sharedIndex] = true;
            }
            sharedIndex++;
          }
        }
      }

      // Get the coordinates for active atoms.
      potential1.getCoordinates(x1);
      potential2.getCoordinates(x2);

      // Collect coordinates for shared atoms (remove alchemical atoms from the superposition).
      double[] x1Shared = new double[nShared * 3];
      double[] x2Shared = new double[nShared * 3];
      int sharedIndex = 0;
      for (int i = 0; i < nActive1; i++) {
        if (sharedAtoms1[i]) {
          int index = i * 3;
          x1Shared[sharedIndex++] = x1[index];
          x1Shared[sharedIndex++] = x1[index + 1];
          x1Shared[sharedIndex++] = x1[index + 2];
        }
      }
      sharedIndex = 0;
      for (int i = 0; i < nActive2; i++) {
        if (sharedAtoms2[i]) {
          int index = i * 3;
          x2Shared[sharedIndex++] = x2[index];
          x2Shared[sharedIndex++] = x2[index + 1];
          x2Shared[sharedIndex++] = x2[index + 2];
        }
      }

      double origRMSD = rmsd(x1Shared, x2Shared, mass);
      logger.info("\n RMSD to topology 2 coordinates from input file:");
      logger.info(format(" RMSD: %12.3f A", origRMSD));

      // Get a fresh copy of the original coordinates.
      double[] origX1 = Arrays.copyOf(x1Shared, nShared * 3);
      double[] origX2 = Arrays.copyOf(x2Shared, nShared * 3);

      // Determine if first topology or second?
      // Explicitly written (as opposed to applySymOp()) since only applied to shared atoms.
      for (int i = 0; i < numSymOps; i++) {
        // Transpose == inverse for orthogonal matrices (value needed later).
        inverse[i] = invertSymOp(symOp[i]);
        // Apply symOp to appropriate coords (validity check).
        // Applied to topology 1 to generate approximate coordinates for topology 2.
        applyCartesianSymOp(origX1, origX1, symOp[i], sharedMask[i]);
        logger.info(format("\n SymOp %3d of %3d between topologies:\n Applied to atoms: %s\n %s\n Inverse SymOp:\n %s",
            i + 1, numSymOps, writeAtomRanges(symOpAtoms.get(i).stream().mapToInt(j -> j).toArray()), symOp[i].toString(), inverse[i].toString()));
      }
      double symOpRMSD = rmsd(origX1, origX2, mass);
      logger.info("\n RMSD to topology 2 coordinates via loaded symop:");
      logger.info(format(" RMSD: %12.3f A", symOpRMSD));

      // TODO: Validate 0th order check for alchemical atoms (i.e., not bonded between sym ops).
      logger.info("\n Checking alchemical atoms for system 1:");
      validateAlchemicalAtoms(nActive1, sharedAtoms1, activeAtoms1);
      logger.info(" Checking alchemical atoms for system 2:");
      validateAlchemicalAtoms(nActive2, sharedAtoms2, activeAtoms2);

      if (logger.isLoggable(Level.FINE)) {
        for (int i = 0; i < nShared; i++) {
          int i3 = i * 3;
          double rmsd = rmsd(new double[]{origX1[i3], origX1[i3 + 1], origX1[i3 + 2]}, new double[]{origX2[i3], origX2[i3 + 1], origX2[i3 + 2]}, new double[]{mass[i]});
          if (rmsd > 0.5) {
            logger.fine(format(" Shared atom %3d has a distance of %8.3f", i + 1, rmsd));
          }
        }
      }
    } catch (Exception ex) {
      logger.severe(" Error parsing SymOp for Dual Topology:\n (" + symOpString + ")" + Utilities.stackTraceToString(ex));
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns true if both force field energies are zero at the ends, and the switching function's
   * first derivative is zero at the end bounds.
   */
  @Override
  public boolean dEdLZeroAtEnds() {
    if (!forceFieldEnergy1.dEdLZeroAtEnds() || !forceFieldEnergy2.dEdLZeroAtEnds()) {
      return false;
    }
    return (switchFunction.highestOrderZeroDerivativeAtZeroBound() > 0);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    boolean ffe1Destroy = forceFieldEnergy1.destroy();
    boolean ffe2Destroy = forceFieldEnergy2.destroy();
    try {
      if (parallelTeam != null) {
        parallelTeam.shutdown();
      }
      return ffe1Destroy && ffe2Destroy;
    } catch (Exception ex) {
      logger.warning(format(" Exception in shutting down DualTopologyEnergy: %s", ex));
      logger.info(Utilities.stackTraceToString(ex));
      return false;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x, boolean verbose) {
    try {
      energyRegion.setX(x);
      energyRegion.setVerbose(verbose);
      parallelTeam.execute(energyRegion);
    } catch (Exception ex) {
      throw new EnergyException(format(" Exception in calculating dual-topology energy: %s", ex));
    }
    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   *
   * <p>The coordinate and gradient arrays are unpacked/packed based on the dual topology.
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    return energyAndGradient(x, g, false);
  }

  /**
   * {@inheritDoc}
   *
   * <p>The coordinate and gradient arrays are unpacked/packed based on the dual topology.
   */
  @Override
  public double energyAndGradient(double[] x, double[] g, boolean verbose) {
    assert Arrays.stream(x).allMatch(Double::isFinite);
    try {
      energyRegion.setX(x);
      energyRegion.setG(g);
      energyRegion.setVerbose(verbose);
      parallelTeam.execute(energyRegion);
    } catch (Exception ex) {
      throw new EnergyException(format(" Exception in calculating dual-topology energy: %s", ex));
    }
    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    if (acceleration == null || acceleration.length < nVariables) {
      acceleration = new double[nVariables];
    }
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    double[] accel = new double[3];
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      atom.getAcceleration(accel);
      if (sharedAtoms1[i]) {
        acceleration[indexCommon++] = accel[0];
        acceleration[indexCommon++] = accel[1];
        acceleration[indexCommon++] = accel[2];
      } else {
        acceleration[indexUnique++] = accel[0];
        acceleration[indexUnique++] = accel[1];
        acceleration[indexUnique++] = accel[2];
      }
    }
    for (int i = 0; i < nActive2; i++) {
      if (!sharedAtoms2[i]) {
        Atom atom = activeAtoms2[i];
        atom.getAcceleration(accel);
        acceleration[indexUnique++] = accel[0];
        acceleration[indexUnique++] = accel[1];
        acceleration[indexUnique++] = accel[2];
      }
    }

    return acceleration;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] x) {
    if (x == null) {
      x = new double[nVariables];
    }
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    for (int i = 0; i < nActive1; i++) {
      Atom a = activeAtoms1[i];
      if (sharedAtoms1[i]) {
        x[indexCommon++] = a.getX();
        x[indexCommon++] = a.getY();
        x[indexCommon++] = a.getZ();
      } else {
        x[indexUnique++] = a.getX();
        x[indexUnique++] = a.getY();
        x[indexUnique++] = a.getZ();
      }
    }
    for (int i = 0; i < nActive2; i++) {
      if (!sharedAtoms2[i]) {
        Atom a = activeAtoms2[i];
        x[indexUnique++] = a.getX();
        x[indexUnique++] = a.getY();
        x[indexUnique++] = a.getZ();
      }
    }
    return x;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCoordinates(@Nullable double[] x) {
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    double[] xyz = new double[3];
    for (int i = 0; i < nActive1; i++) {
      if (sharedAtoms1[i]) {
        xyz[0] = x[indexCommon++];
        xyz[1] = x[indexCommon++];
        xyz[2] = x[indexCommon++];
      } else {
        xyz[0] = x[indexUnique++];
        xyz[1] = x[indexUnique++];
        xyz[2] = x[indexUnique++];
      }
      Atom a = activeAtoms1[i];
      a.setXYZ(xyz);
    }
    // Reset the common index.
    indexCommon = 0;
    for (int i = 0; i < nActive2; i++) {
      if (sharedAtoms2[i]) {
        xyz[0] = x[indexCommon++];
        xyz[1] = x[indexCommon++];
        xyz[2] = x[indexCommon++];
      } else {
        xyz[0] = x[indexUnique++];
        xyz[1] = x[indexUnique++];
        xyz[2] = x[indexUnique++];
      }
      Atom a = activeAtoms2[i];
      a.setXYZ(xyz);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Crystal getCrystal() {
    // TODO: Handle the case where each system has a unique crystal instance (e.g., different space groups).
    if (useSymOp && logger.isLoggable(Level.FINE)) {
      logger.fine(" Get method only returned first crystal.");
    }
    return potential1.getCrystal();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    // TODO: Handle the case where each system has a unique crystal instance (e.g., different space groups).
    if (useSymOp && logger.isLoggable(Level.FINE)) {
      logger.fine(" Both systems set to the same crystal.");
    }
    potential1.setCrystal(crystal);
    potential2.setCrystal(crystal);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public STATE getEnergyTermState() {
    return state;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
    potential1.setEnergyTermState(state);
    potential2.setEnergyTermState(state);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getLambda() {
    return lambda;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setLambda(double lambda) {
    if (lambda <= 1.0 && lambda >= 0.0) {
      this.lambda = lambda;
      double oneMinusLambda = 1.0 - lambda;
      lambdaInterface1.setLambda(lambda);
      lambdaInterface2.setLambda(oneMinusLambda);

      f1L = switchFunction.valueAt(lambda);
      dF1dL = switchFunction.firstDerivative(lambda);
      d2F1dL2 = switchFunction.secondDerivative(lambda);

      f2L = switchFunction.valueAt(oneMinusLambda);
      dF2dL = -1.0 * switchFunction.firstDerivative(oneMinusLambda);
      d2F2dL2 = switchFunction.secondDerivative(oneMinusLambda);
    } else {
      String message = format("Lambda value %8.3f is not in the range [0..1].", lambda);
      throw new IllegalArgumentException(message);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    return mass;
  }

  /**
   * Returns the number of shared variables (3 * number of shared atoms).
   *
   * @return An int
   */
  public int getNumSharedVariables() {
    return 3 * nShared;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return nVariables;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    if (previousAcceleration == null || previousAcceleration.length < nVariables) {
      previousAcceleration = new double[nVariables];
    }
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    double[] prev = new double[3];
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      atom.getPreviousAcceleration(prev);
      if (sharedAtoms1[i]) {
        previousAcceleration[indexCommon++] = prev[0];
        previousAcceleration[indexCommon++] = prev[1];
        previousAcceleration[indexCommon++] = prev[2];
      } else {
        previousAcceleration[indexUnique++] = prev[0];
        previousAcceleration[indexUnique++] = prev[1];
        previousAcceleration[indexUnique++] = prev[2];
      }
    }
    for (int i = 0; i < nActive2; i++) {
      Atom atom = activeAtoms2[i];
      if (!sharedAtoms2[i]) {
        atom.getPreviousAcceleration(prev);
        previousAcceleration[indexUnique++] = prev[0];
        previousAcceleration[indexUnique++] = prev[1];
        previousAcceleration[indexUnique++] = prev[2];
      }
    }

    return previousAcceleration;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return scaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(double[] scaling) {
    this.scaling = scaling;
  }

  /**
   * Returns the switching function used by this DualTopologyEnergy; presently, switching functions
   * are immutable, and cannot be changed once a DualTopologyEnergy is constructed.
   *
   * @return The switching function.
   */
  public UnivariateSwitchingFunction getSwitch() {
    return switchFunction;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  @Override
  public List<Potential> getUnderlyingPotentials() {
    List<Potential> under = new ArrayList<>(2);
    under.add(forceFieldEnergy1);
    under.add(forceFieldEnergy2);
    under.addAll(forceFieldEnergy1.getUnderlyingPotentials());
    under.addAll(forceFieldEnergy2.getUnderlyingPotentials());
    return under;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    return variableTypes;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    if (velocity == null || velocity.length < nVariables) {
      velocity = new double[nVariables];
    }
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    double[] vel = new double[3];
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      atom.getVelocity(vel);
      if (sharedAtoms1[i]) {
        velocity[indexCommon++] = vel[0];
        velocity[indexCommon++] = vel[1];
        velocity[indexCommon++] = vel[2];
      } else {
        velocity[indexUnique++] = vel[0];
        velocity[indexUnique++] = vel[1];
        velocity[indexUnique++] = vel[2];
      }
    }
    for (int i = 0; i < nActive2; i++) {
      if (!sharedAtoms2[i]) {
        Atom atom = activeAtoms2[i];
        atom.getVelocity(vel);
        velocity[indexUnique++] = vel[0];
        velocity[indexUnique++] = vel[1];
        velocity[indexUnique++] = vel[2];
      }
    }

    return velocity;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    double e1 =
        f1L * d2EdL2_1
            + 2.0 * dF1dL * dEdL_1
            + d2F1dL2 * energy1
            + f2L * restraintd2EdL2_1
            + 2.0 * dF2dL * restraintdEdL_1
            + d2F2dL2 * restraintEnergy1;
    double e2;
    if (!useFirstSystemBondedEnergy) {
      e2 = f2L * d2EdL2_2
          + 2.0 * dF2dL * dEdL_2
          + d2F2dL2 * energy2
          + f1L * restraintd2EdL2_2
          + 2.0 * dF1dL * restraintdEdL_2
          + d2F1dL2 * restraintEnergy2;
    } else {
      e2 = f2L * d2EdL2_2
          + 2.0 * dF2dL * dEdL_2
          + d2F2dL2 * energy2
          - f2L * restraintd2EdL2_2
          - 2.0 * dF2dL * restraintdEdL_2
          - d2F2dL2 * restraintEnergy2;
    }
    return e1 + e2;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    // Subtraction is implicit, as dEdL_2 and dF2dL are both multiplied by
    // negative 1 when set.
    double e1 = f1L * dEdL_1 + dF1dL * energy1 + f2L * restraintdEdL_1 + dF2dL * restraintEnergy1;
    double e2;
    if (!useFirstSystemBondedEnergy) {
      e2 = f2L * dEdL_2 + dF2dL * energy2 + f1L * restraintdEdL_2 + dF1dL * restraintEnergy2;
    } else {
      e2 = f2L * dEdL_2 + dF2dL * energy2 - f2L * restraintdEdL_2 - dF2dL * restraintEnergy2;
    }
    return e1 + e2;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] g) {
    if (g == null) {
      g = new double[nVariables];
    }

    int index = 0;
    int indexCommon = 0;
    int indexUnique = nShared * 3;
    // Coordinate Gradient from Topology 1.
    for (int i = 0; i < nActive1; i++) {
      if (sharedAtoms1[i]) {
        g[indexCommon++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
        g[indexCommon++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
        g[indexCommon++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
      } else {
        g[indexUnique++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
        g[indexUnique++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
        g[indexUnique++] =
            f1L * gl1[index] + dF1dL * g1[index] + f2L * rgl1[index] + dF2dL * rg1[index++];
      }
    }

    // Coordinate Gradient from Topology 2.
    if (!useFirstSystemBondedEnergy) {
      index = 0;
      indexCommon = 0;
      for (int i = 0; i < nActive2; i++) {
        if (sharedAtoms2[i]) {
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
        } else {
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] - f1L * rgl2[index]
              + dF1dL * rg2[index++]);
        }
      }
    } else {
      index = 0;
      indexCommon = 0;
      for (int i = 0; i < nActive2; i++) {
        if (sharedAtoms2[i]) {
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
          g[indexCommon++] += (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
        } else {
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
          g[indexUnique++] = (-f2L * gl2[index] + dF2dL * g2[index] + f2L * rgl2[index]
              - dF2dL * rg2[index++]);
        }
      }

    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    double[] accel = new double[3];
    int indexCommon = 0;
    int indexUnique = 3 * nShared;
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      if (sharedAtoms1[i]) {
        accel[0] = acceleration[indexCommon++];
        accel[1] = acceleration[indexCommon++];
        accel[2] = acceleration[indexCommon++];
      } else {
        accel[0] = acceleration[indexUnique++];
        accel[1] = acceleration[indexUnique++];
        accel[2] = acceleration[indexUnique++];
      }
      atom.setAcceleration(accel);
    }
    indexCommon = 0;
    for (int i = 0; i < nActive2; i++) {
      Atom atom = activeAtoms2[i];
      if (sharedAtoms2[i]) {
        accel[0] = acceleration[indexCommon++];
        accel[1] = acceleration[indexCommon++];
        accel[2] = acceleration[indexCommon++];
      } else {
        accel[0] = acceleration[indexUnique++];
        accel[1] = acceleration[indexUnique++];
        accel[2] = acceleration[indexUnique++];
      }
      atom.setAcceleration(accel);
    }
  }

  /**
   * setParallel.
   *
   * @param parallel a boolean.
   */
  public void setParallel(boolean parallel) {
    if (parallelTeam != null) {
      try {
        parallelTeam.shutdown();
      } catch (Exception e) {
        logger.severe(format(" Exception in shutting down old ParallelTeam for DualTopologyEnergy: %s", e));
      }
    }
    parallelTeam = parallel ? new ParallelTeam(2) : new ParallelTeam(1);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    double[] prev = new double[3];
    int indexCommon = 0;
    int indexUnique = 3 * nShared;
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      if (sharedAtoms1[i]) {
        prev[0] = previousAcceleration[indexCommon++];
        prev[1] = previousAcceleration[indexCommon++];
        prev[2] = previousAcceleration[indexCommon++];
      } else {
        prev[0] = previousAcceleration[indexUnique++];
        prev[1] = previousAcceleration[indexUnique++];
        prev[2] = previousAcceleration[indexUnique++];
      }
      atom.setPreviousAcceleration(prev);
    }
    indexCommon = 0;
    for (int i = 0; i < nActive2; i++) {
      Atom atom = activeAtoms2[i];
      if (sharedAtoms2[i]) {
        prev[0] = previousAcceleration[indexCommon++];
        prev[1] = previousAcceleration[indexCommon++];
        prev[2] = previousAcceleration[indexCommon++];
      } else {
        prev[0] = previousAcceleration[indexUnique++];
        prev[1] = previousAcceleration[indexUnique++];
        prev[2] = previousAcceleration[indexUnique++];
      }
      atom.setPreviousAcceleration(prev);
    }
  }

  /**
   * Sets the printOnFailure flag; if override is true, over-rides any existing property. Essentially
   * sets the default value of printOnFailure for an algorithm. For example, rotamer optimization
   * will generally run into force field issues in the normal course of execution as it tries
   * unphysical self and pair configurations, so the algorithm should not print out a large number of
   * error PDBs.
   *
   * @param onFail   To set
   * @param override Override properties
   */
  public void setPrintOnFailure(boolean onFail, boolean override) {
    forceFieldEnergy1.setPrintOnFailure(onFail, override);
    forceFieldEnergy2.setPrintOnFailure(onFail, override);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    double[] vel = new double[3];
    int indexCommon = 0;
    int indexUnique = 3 * nShared;
    for (int i = 0; i < nActive1; i++) {
      Atom atom = activeAtoms1[i];
      if (sharedAtoms1[i]) {
        vel[0] = velocity[indexCommon++];
        vel[1] = velocity[indexCommon++];
        vel[2] = velocity[indexCommon++];
      } else {
        vel[0] = velocity[indexUnique++];
        vel[1] = velocity[indexUnique++];
        vel[2] = velocity[indexUnique++];
      }
      atom.setVelocity(vel);
    }

    indexCommon = 0;
    for (int i = 0; i < nActive2; i++) {
      Atom atom = activeAtoms2[i];
      if (sharedAtoms2[i]) {
        vel[0] = velocity[indexCommon++];
        vel[1] = velocity[indexCommon++];
        vel[2] = velocity[indexCommon++];
      } else {
        vel[0] = velocity[indexUnique++];
        vel[1] = velocity[indexUnique++];
        vel[2] = velocity[indexUnique++];
      }
      atom.setVelocity(vel);
    }
  }

  /**
   * Moves two shared atoms together if there is a small discrepancy (such as that caused by the
   * mutator script).
   *
   * @param a1      Atom from topology 1
   * @param a2      Atom from topology 2
   * @param warnlev Logging level to use when warning about small movements
   */
  private void reconcileAtoms(Atom a1, Atom a2, Level warnlev) {
    double dist = 0;
    double[] xyz1 = a1.getXYZ(null);
    double[] xyz2 = a2.getXYZ(null);
    double[] xyzAv = new double[3];
    for (int i = 0; i < 3; i++) {
      double dx = xyz1[i] - xyz2[i];
      dist += (dx * dx);
      xyzAv[i] = xyz1[i] + (0.5 * dx);
    }

    // Square of the maximum distance permissible between two shared atoms.
    double maxDist2 = 0.09;

    /*
     Square of the minimum distance between shared atoms which will cause a
     warning. Intended to cover anything larger than a rounding error.
    */
    double minDistWarn2 = 0.00001;

    if (dist > maxDist2) {
      logger.log(
          Level.SEVERE,
          format(
              " Distance between atoms %s " + "and %s is %7.4f >> maximum allowed %7.4f",
              a1, a2, sqrt(dist), sqrt(maxDist2)));
    } else if (dist > minDistWarn2) {
      logger.log(
          warnlev,
          format(
              " Distance between atoms %s " + "and %s is %7.4f; moving atoms together.",
              a1, a2, sqrt(dist)));
      a1.setXYZ(xyzAv);
      a2.setXYZ(xyzAv);
    } else if (dist > 0) {
      // Silently move them together; probably just a rounding error.
      a1.setXYZ(xyzAv);
      a2.setXYZ(xyzAv);
    }
  }

  private void packGradient(double[] x, double[] g) {
    if (g == null) {
      g = new double[nVariables];
    }
    int indexCommon = 0;
    int indexUnique = nShared * 3;

    // Coordinate Gradient from Topology 1.
    int index = 0;
    for (int i = 0; i < nActive1; i++) {
      if (sharedAtoms1[i]) {
        g[indexCommon++] = f1L * g1[index] + f2L * rg1[index++];
        g[indexCommon++] = f1L * g1[index] + f2L * rg1[index++];
        g[indexCommon++] = f1L * g1[index] + f2L * rg1[index++];
      } else {
        g[indexUnique++] = f1L * g1[index] + f2L * rg1[index++];
        g[indexUnique++] = f1L * g1[index] + f2L * rg1[index++];
        g[indexUnique++] = f1L * g1[index] + f2L * rg1[index++];
      }
    }

    if (!useFirstSystemBondedEnergy) {
      // Coordinate Gradient from Topology 2.
      indexCommon = 0;
      index = 0;
      for (int i = 0; i < nActive2; i++) {
        if (sharedAtoms2[i]) {
          g[indexCommon++] += f2L * g2[index] + f1L * rg2[index++];
          g[indexCommon++] += f2L * g2[index] + f1L * rg2[index++];
          g[indexCommon++] += f2L * g2[index] + f1L * rg2[index++];
        } else {
          g[indexUnique++] = f2L * g2[index] + f1L * rg2[index++];
          g[indexUnique++] = f2L * g2[index] + f1L * rg2[index++];
          g[indexUnique++] = f2L * g2[index] + f1L * rg2[index++];
        }
      }
    } else {
      // Coordinate Gradient from Topology 2.
      indexCommon = 0;
      index = 0;
      for (int i = 0; i < nActive2; i++) {
        if (sharedAtoms2[i]) {
          g[indexCommon++] += f2L * g2[index] - f2L * rg2[index++];
          g[indexCommon++] += f2L * g2[index] - f2L * rg2[index++];
          g[indexCommon++] += f2L * g2[index] - f2L * rg2[index++];
        } else {
          g[indexUnique++] = f2L * g2[index] - f2L * rg2[index++];
          g[indexUnique++] = f2L * g2[index] - f2L * rg2[index++];
          g[indexUnique++] = f2L * g2[index] - f2L * rg2[index++];
        }
      }
    }

    scaleCoordinatesAndGradient(x, g);
  }

  private void unpackCoordinates(double[] x) {

    // Unscale the coordinates.
    unscaleCoordinates(x);

    int index = 0;
    int indexCommon = 0;
    int indexUnique = 3 * nShared;
    for (int i = 0; i < nActive1; i++) {
      if (sharedAtoms1[i]) {
        x1[index++] = x[indexCommon++];
        x1[index++] = x[indexCommon++];
        x1[index++] = x[indexCommon++];
      } else {
        x1[index++] = x[indexUnique++];
        x1[index++] = x[indexUnique++];
        x1[index++] = x[indexUnique++];
      }
    }

    index = 0;
    indexCommon = 0;
    for (int i = 0; i < nActive2; i++) {
      if (sharedAtoms2[i]) {
        x2[index++] = x[indexCommon++];
        x2[index++] = x[indexCommon++];
        x2[index++] = x[indexCommon++];
      } else {
        x2[index++] = x[indexUnique++];
        x2[index++] = x[indexUnique++];
        x2[index++] = x[indexUnique++];
      }
    }
  }

  /**
   * Reload the common atomic masses. Intended for quad-topology dual force field corrections, where
   * common atoms may have slightly different masses; this was found to be the case between AMOEBA
   * 2013 and AMBER99SB carbons.
   *
   * @param secondTopology Load from second topology
   */
  void reloadCommonMasses(boolean secondTopology) {
    int commonIndex = 0;
    if (secondTopology) {
      for (int i = 0; i < nActive2; i++) {
        Atom a = activeAtoms2[i];
        double m = a.getMass();
        if (sharedAtoms2[i]) {
          mass[commonIndex++] = m;
          mass[commonIndex++] = m;
          mass[commonIndex++] = m;
        }
      }
    } else {
      for (int i = 0; i < nActive1; i++) {
        Atom a = activeAtoms1[i];
        double m = a.getMass();
        if (sharedAtoms1[i]) {
          mass[commonIndex++] = m;
          mass[commonIndex++] = m;
          mass[commonIndex++] = m;
        }
      }
    }
  }

  /**
   * Bonded terms for alchemical atoms are used to prevent them from being unconstrained.
   * These restraints have the potential to interact with symmetry operators applied in dual-topology simulations.
   *
   * @param nActive     Number of active atoms (shared and alchemical) in simulation for this system.
   * @param sharedAtoms Mask to determine which atoms are shared vs alchemical.
   * @param activeAtoms Active atoms for this simulation.
   */
  private void validateAlchemicalAtoms(int nActive, boolean[] sharedAtoms, Atom[] activeAtoms) {
    int numSymOps = symOpAtoms.size();
    // Loop through atoms in simulation.
    for (int i = 0; i < nActive; i++) {
      // Determine if this atom is alchemical.
      if (!sharedAtoms[i]) {
        Atom alchAtom = activeAtoms[i];
        int symOpForAlchemicalAtom = -1;
        for (int j = 0; j < numSymOps; j++) {
          List<Integer> symOpAtomGroup = symOpAtoms.get(j);
          if (symOpAtomGroup.contains(i)) {
            symOpForAlchemicalAtom = j;
            break;
          }
        }

        StringBuilder sb = new StringBuilder();
        sb.append(format("\n Alchemical atom %s\n Symmetry group %s:\n", alchAtom, symOpForAlchemicalAtom + 1));
        List<Integer> symOpAtomGroup = symOpAtoms.get(symOpForAlchemicalAtom);
        for (Integer j : symOpAtomGroup) {
          if (j == i) {
            continue;
          }
          sb.append(format("  %d %s\n", j + 1, activeAtoms[j]));
        }
        sb.append(format("\n 1-2, 1-3 or 1-4 atoms in other symmetry groups:\n"));

        // Collect atoms that are 1-2, 1-3 or 1-4 to the alchemical atom.
        ArrayList<Integer> bondedAtoms = new ArrayList<>();
        addUniqueBondedIndices(bondedAtoms, alchAtom.get12List());
        addUniqueBondedIndices(bondedAtoms, alchAtom.get13List());
        addUniqueBondedIndices(bondedAtoms, alchAtom.get14List());
        // Loop through atoms bonded to alchemical atom to determine if they're part of
        // a different SymOp group.
        boolean conflict = false;
        for (int j = 0; j < numSymOps; j++) {
          // If the bonded atom is in the same symmetry group as the alchemical atom, then there is no conflict.
          if (j == symOpForAlchemicalAtom) {
            continue;
          }
          symOpAtomGroup = symOpAtoms.get(j);
          for (Integer integer : bondedAtoms) {
            if (symOpAtomGroup.contains(integer)) {
              conflict = true;
              sb.append(format("  %d %s uses symmetry group %2d.\n", integer + 1, activeAtoms[integer], j + 1));
            }
          }
        }
        if (conflict) {
          logger.info(sb.toString());
        }
      }
    }
  }

  /**
   * Take list of aggregated unique indices based on input list of bonded atoms..
   *
   * @param uniqueIndices List of unique indices of bonded atoms.
   * @param bondedAtoms   Atoms that are bonded to the atom of interest.
   */
  private void addUniqueBondedIndices(ArrayList<Integer> uniqueIndices, List<Atom> bondedAtoms) {
    // Loop through bonded atoms
    for (Atom a : bondedAtoms) {
      int index = a.getIndex() - 1;
      // Determine if bonded atom's index is already in unique list.
      if (!uniqueIndices.contains(index)) {
        // Add to list if it is not already included.
        uniqueIndices.add(index);
      }
    }
  }

  private class EnergyRegion extends ParallelRegion {

    private final Energy1Section e1sect;
    private final Energy2Section e2sect;
    private double[] x;
    private double[] g;
    private boolean gradient = false;
    private boolean verbose = false;

    EnergyRegion() {
      e1sect = new Energy1Section();
      e2sect = new Energy2Section();
    }

    @Override
    public void finish() {
      // Apply the dual-topology scaling for the total energy.
      if (!useFirstSystemBondedEnergy) {
        totalEnergy =
            f1L * energy1 + f2L * restraintEnergy1 + f2L * energy2 + f1L * restraintEnergy2;
      } else {
        totalEnergy =
            f1L * energy1 + f2L * restraintEnergy1 + f2L * energy2 - f2L * restraintEnergy2;
      }

      if (gradient) {
        packGradient(x, g);
      } else {
        scaleCoordinates(x);
      }

      if (verbose) {
        logger.info(format(" Total dual-topology energy: %12.4f", totalEnergy));
      }
      setVerbose(false);
      setGradient(false);
    }

    @Override
    public void run() throws Exception {
      execute(e1sect, e2sect);
    }

    public void setG(double[] g) {
      this.g = g;
      setGradient(true);
    }

    public void setGradient(boolean grad) {
      this.gradient = grad;
      e1sect.setGradient(grad);
      e2sect.setGradient(grad);
    }

    public void setVerbose(boolean verb) {
      this.verbose = verb;
      e1sect.setVerbose(verb);
      e2sect.setVerbose(verb);
    }

    public void setX(double[] x) {
      this.x = x;
    }

    @Override
    public void start() {
      unpackCoordinates(x);
    }
  }

  private class Energy1Section extends ParallelSection {

    private boolean gradient = false;
    private boolean verbose = false;

    @Override
    public void run() {
      if (gradient) {
        fill(gl1, 0.0);
        fill(rgl1, 0.0);
        energy1 = potential1.energyAndGradient(x1, g1, verbose);
        dEdL_1 = lambdaInterface1.getdEdL();
        d2EdL2_1 = lambdaInterface1.getd2EdL2();
        lambdaInterface1.getdEdXdL(gl1);
        if (doValenceRestraint1) {
          if (verbose) {
            logger.info(" Calculating lambda bonded terms for topology 1");
          }
          forceFieldEnergy1.setLambdaBondedTerms(true, useFirstSystemBondedEnergy);
          if (forceFieldEnergy1 instanceof OpenMMEnergy openMMEnergy) {
            // Use FFX to calculate the energy and gradient for the restraints.
            restraintEnergy1 = openMMEnergy.energyAndGradientFFX(x1, rg1, verbose);
          } else {
            restraintEnergy1 = forceFieldEnergy1.energyAndGradient(x1, rg1, verbose);
          }
          restraintdEdL_1 = forceFieldEnergy1.getdEdL();
          restraintd2EdL2_1 = forceFieldEnergy1.getd2EdL2();
          forceFieldEnergy1.getdEdXdL(rgl1);
          forceFieldEnergy1.setLambdaBondedTerms(false, false);
        } else {
          restraintEnergy1 = 0.0;
          restraintdEdL_1 = 0.0;
          restraintd2EdL2_1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n", f1L * energy1, f2L * restraintEnergy1));
          logger.fine(format(" Topology 1:    %15.8f * (%.2f)", energy1, f1L));
          logger.fine(format(" T1 Restraints: %15.8f * (%.2f)", restraintEnergy1, f2L));
        }
      } else {
        energy1 = potential1.energy(x1, verbose);
        if (doValenceRestraint1) {
          if (verbose) {
            logger.info(" Calculating lambda bonded terms for topology 1");
          }
          forceFieldEnergy1.setLambdaBondedTerms(true, useFirstSystemBondedEnergy);
          if (forceFieldEnergy1 instanceof OpenMMEnergy openMMEnergy) {
            // Use FFX to calculate the energy and gradient for the restraints.
            restraintEnergy1 = openMMEnergy.energyFFX(x1, verbose);
          } else {
            restraintEnergy1 = potential1.energy(x1, verbose);
          }
          forceFieldEnergy1.setLambdaBondedTerms(false, false);
        } else {
          restraintEnergy1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n", f1L * energy1, f2L * restraintEnergy1));
        }
      }
    }

    public void setGradient(boolean grad) {
      this.gradient = grad;
    }

    public void setVerbose(boolean verb) {
      this.verbose = verb;
    }
  }

  private class Energy2Section extends ParallelSection {

    private boolean gradient = false;
    private boolean verbose = false;

    @Override
    public void run() {
      // Number of symmetry operators to be applied to the system.
      int numSymOps = 0;
      if (useSymOp) {
        numSymOps = symOpAtoms.size();
        for (int i = 0; i < numSymOps; i++) {
          // Apply sym ops to system 1 coords (duplicate of x1) to get x2.
          applyCartesianSymOp(x2, x2, symOp[i], mask[i]);
        }
      }

      if (gradient) {
        fill(gl2, 0.0);
        fill(rgl2, 0.0);

        // Compute the energy and gradient of topology 2.
        energy2 = potential2.energyAndGradient(x2, g2, verbose);
        dEdL_2 = -lambdaInterface2.getdEdL();
        d2EdL2_2 = lambdaInterface2.getd2EdL2();
        lambdaInterface2.getdEdXdL(gl2);
        if (useSymOp) {
          // Rotate the gradient back for appropriate atoms.
          for (int i = 0; i < numSymOps; i++) {
            applyCartesianSymRot(g2, g2, inverse[i], mask[i]);
            applyCartesianSymRot(gl2, gl2, inverse[i], mask[i]);
          }
        }
        if (doValenceRestraint2) {
          if (verbose) {
            logger.info(" Calculating lambda bonded terms for topology 2");
          }
          forceFieldEnergy2.setLambdaBondedTerms(true, useFirstSystemBondedEnergy);
          if (forceFieldEnergy2 instanceof OpenMMEnergy openMMEnergy) {
            // Use FFX to calculate the energy and gradient for the restraints.
            restraintEnergy2 = openMMEnergy.energyAndGradientFFX(x2, rg2, verbose);
          } else {
            restraintEnergy2 = forceFieldEnergy2.energyAndGradient(x2, rg2, verbose);
          }
          restraintdEdL_2 = -forceFieldEnergy2.getdEdL();
          restraintd2EdL2_2 = forceFieldEnergy2.getd2EdL2();
          forceFieldEnergy2.getdEdXdL(rgl2);
          if (useSymOp) {
            // Rotate the gradient back for appropriate atoms.
            for (int i = 0; i < numSymOps; i++) {
              applyCartesianSymRot(rg2, rg2, inverse[i], mask[i]);
              applyCartesianSymRot(rgl2, rgl2, inverse[i], mask[i]);
            }
          }
          forceFieldEnergy2.setLambdaBondedTerms(false, false);
        } else {
          restraintEnergy2 = 0.0;
          restraintdEdL_2 = 0.0;
          restraintd2EdL2_2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Topology 2 Energy & Restraints: %15.8f %15.8f", f2L * energy2, f1L * restraintEnergy2));
          logger.fine(format(" Topology 2:    %15.8f * (%.2f)", energy2, f2L));
          logger.fine(format(" T2 Restraints: %15.8f * (%.2f)", restraintEnergy2, f1L));
        }
      } else {
        energy2 = potential2.energy(x2, verbose);
        if (doValenceRestraint2) {
          if (verbose) {
            logger.info(" Calculating lambda bonded terms for topology 2");
          }
          forceFieldEnergy2.setLambdaBondedTerms(true, useFirstSystemBondedEnergy);
          if (forceFieldEnergy2 instanceof OpenMMEnergy openMMEnergy) {
            // Use FFX to calculate the energy and gradient for the restraints.
            restraintEnergy2 = openMMEnergy.energyFFX(x2, verbose);
          } else {
            restraintEnergy2 = potential2.energy(x2, verbose);
          }
          forceFieldEnergy2.setLambdaBondedTerms(false, false);
        } else {
          restraintEnergy2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Topology 2 Energy & Restraints: %15.8f %15.8f", f2L * energy2, f1L * restraintEnergy2));
        }
      }
    }

    public void setGradient(boolean grad) {
      this.gradient = grad;
    }

    public void setVerbose(boolean verb) {
      this.verbose = verb;
    }
  }
}
