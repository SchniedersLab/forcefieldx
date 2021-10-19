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
package ffx.potential.nonbonded;

import static ffx.potential.parameters.ForceField.toEnumForm;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedInteger;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The Van der Waals class computes Van der Waals interaction in parallel using a {@link
 * ffx.potential.nonbonded.NeighborList} for any {@link ffx.crystal.Crystal}. The repulsive power
 * (e.g. 12), attractive power (e.g. 6) and buffering (e.g. for the AMOEBA buffered-14-7) can all be
 * specified such that both Lennard-Jones and AMOEBA are supported.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class VanDerWaals implements MaskingInterface, LambdaInterface {

  private static final Logger logger = Logger.getLogger(VanDerWaals.class.getName());
  private static final byte HARD = 0;
  private static final byte SOFT = 1;
  private static final byte XX = 0;
  private static final byte YY = 1;
  private static final byte ZZ = 2;
  private final boolean doLongRangeCorrection;
  // *************************************************************************
  // Parallel variables.
  private final ParallelTeam parallelTeam;
  private final int threadCount;
  private final IntegerSchedule pairwiseSchedule;
  private final SharedInteger sharedInteractions;
  private final SharedDouble sharedEnergy;
  private final SharedDouble shareddEdL;
  private final SharedDouble sharedd2EdL2;
  private final VanDerWaalsRegion vanDerWaalsRegion;
  /** Timing variables. */
  private final long[] initializationTime;

  private final long[] vdwTime;
  private final long[] reductionTime;
  private final VanDerWaalsForm vdwForm;
  private final NonbondedCutoff nonbondedCutoff;
  private final MultiplicativeSwitch multiplicativeSwitch;
  /** This field specifies resolution for multi-scale modeling. */
  private Resolution resolution = null;
  /** Boundary conditions and crystal symmetry. */
  private Crystal crystal;
  /** An array of all atoms in the system. */
  private Atom[] atoms;
  /** Previous atom array. */
  private Atom[] previousAtoms;
  /** Specification of the molecular index for each atom. */
  private int[] molecule;
  /** The Force Field that defines the Van der Waals interactions. */
  private ForceField forceField;
  /** An array of whether each atom in the system should be used in the calculations. */
  private boolean[] use = null;
  /** A local convenience variable equal to atoms.length. */
  private int nAtoms;
  /** A local convenience variable equal to the number of crystal symmetry operators. */
  private int nSymm;
  /** ********************************************************************** Lambda variables. */
  private boolean gradient;

  private boolean lambdaTerm;
  private boolean esvTerm;
  private boolean[] isSoft;
  /**
   * There are 2 softCore arrays of length nAtoms.
   *
   * <p>The first is used for atoms in the outer loop that are hard. This mask equals: false for
   * inner loop hard atoms true for inner loop soft atoms
   *
   * <p>The second is used for atoms in the outer loop that are soft. This mask equals: true for
   * inner loop hard atoms false for inner loop soft atoms
   */
  private boolean[][] softCore;
  /** Turn on inter-molecular softcore interactions using molecular index. */
  private boolean intermolecularSoftcore = false;

  // *************************************************************************
  // Coordinate arrays.
  /** Turn on intra-molecular softcore interactions using molecular index. */
  private boolean intramolecularSoftcore = false;
  /** Current value of the lambda state variable. */
  private double lambda = 1.0;
  /** Exponent on lambda (beta). */
  private double vdwLambdaExponent = 3.0;
  /** Offset in Angstroms (alpha). */
  private double vdwLambdaAlpha = 0.25;
  /** Polymorphic inner class to set sc1,sc2,dsc1,etc only when necessary. [nThreads] */
  private LambdaFactors[] lambdaFactors = null;

  private double sc1 = 0.0; // alpha * (1 - lambdaProduct)^2
  private double sc2 = 1.0; // lambdaProduct

  // *************************************************************************
  // Force field parameters and constants for the Buffered-14-7 potential.
  private double dsc1dL = 0.0;
  private double dsc2dL = 0.0;
  private double d2sc1dL2 = 0.0;
  private double d2sc2dL2 = 0.0;
  /** Generalized extended system variables. */
  private ExtendedSystem esvSystem;
  private int numESVs = 0;

  /** A local copy of atomic coordinates, including reductions on the hydrogen atoms. */
  private double[] coordinates;
  /** Reduced coordinates of size: [nSymm][nAtoms * 3] */
  private double[][] reduced;

  private double[] reducedXYZ;
  /** Neighbor lists for each atom. Size: [nSymm][nAtoms][nNeighbors] */
  private int[][][] neighborLists;
  /** A local reference to the atom class of each atom in the system. */
  private int[] atomClass;
  /**
   * Hydrogen atom vdW sites are located toward their heavy atom relative to their nucleus. This is a
   * look-up that gives the heavy atom index for each hydrogen.
   */
  private int[] reductionIndex;

  private int[][] mask12;
  private int[][] mask13;
  private int[][] mask14;
  /**
   * Each hydrogen vdW site is located a fraction of the way from the heavy atom nucleus to the
   * hydrogen nucleus (~0.9).
   */
  private double[] reductionValue;

  private boolean reducedHydrogens;
  private double longRangeCorrection;
  private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
  /** Cartesian coordinate gradient. */
  private AtomicDoubleArray3D grad;
  /** Lambda derivative of the Cartesian coordinate gradient. */
  private AtomicDoubleArray3D lambdaGrad;
  /**
   * The neighbor-list includes 1-2 and 1-3 interactions, which are masked out in the Van der Waals
   * energy code. The AMOEBA force field includes 1-4 interactions fully.
   */
  private NeighborList neighborList;

  private boolean neighborListOnly = true;
  private long initializationTotal, vdwTotal, reductionTotal;

  public VanDerWaals() {
    // Empty constructor for use with VanDerWaalsTornado
    doLongRangeCorrection = false;
    parallelTeam = null;
    threadCount = 0;
    pairwiseSchedule = null;
    sharedInteractions = null;
    sharedEnergy = null;
    shareddEdL = null;
    sharedd2EdL2 = null;
    vanDerWaalsRegion = null;
    initializationTime = null;
    vdwTime = null;
    reductionTime = null;
    vdwForm = null;
    nonbondedCutoff = null;
    multiplicativeSwitch = null;
  }

  /**
   * The VanDerWaals class constructor.
   *
   * @param atoms the Atom array to do Van Der Waals calculations on.
   * @param molecule the molecule number for each atom.
   * @param crystal The boundary conditions.
   * @param forceField the ForceField parameters to apply.
   * @param parallelTeam The parallel environment.
   * @param vdwCutoff a double.
   * @param neighborListCutoff a double.
   * @since 1.0
   */
  public VanDerWaals(
      Atom[] atoms,
      int[] molecule,
      Crystal crystal,
      ForceField forceField,
      ParallelTeam parallelTeam,
      double vdwCutoff,
      double neighborListCutoff) {
    this.atoms = atoms;
    this.molecule = molecule;
    this.crystal = crystal;
    this.parallelTeam = parallelTeam;
    this.forceField = forceField;

    nAtoms = atoms.length;
    nSymm = crystal.spaceGroup.getNumberOfSymOps();

    vdwForm = new VanDerWaalsForm(forceField);
    reducedHydrogens = forceField.getBoolean("REDUCE_HYDROGENS", true);

    // Lambda parameters.
    lambdaTerm =
        forceField.getBoolean("VDW_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false));
    if (lambdaTerm) {
      shareddEdL = new SharedDouble();
      sharedd2EdL2 = new SharedDouble();
      vdwLambdaAlpha = forceField.getDouble("VDW_LAMBDA_ALPHA", 0.25);
      vdwLambdaExponent = forceField.getDouble("VDW_LAMBDA_EXPONENT", 3.0);
      if (vdwLambdaAlpha < 0.0) {
        logger.warning(
            format(
                " Invalid value %8.3g for vdw-lambda-alpha; must be greater than or equal to 0. Resetting to 0.25.",
                vdwLambdaAlpha));
        vdwLambdaAlpha = 0.25;
      }
      if (vdwLambdaExponent < 1.0) {
        logger.warning(
            format(
                " Invalid value %8.3g for vdw-lambda-exponent; must be greater than or equal to 1. Resetting to 3.",
                vdwLambdaExponent));
        vdwLambdaExponent = 3.0;
      }
      intermolecularSoftcore = forceField.getBoolean("INTERMOLECULAR_SOFTCORE", false);
      intramolecularSoftcore = forceField.getBoolean("INTRAMOLECULAR_SOFTCORE", false);
    } else {
      shareddEdL = null;
      sharedd2EdL2 = null;
    }

    // Parallel constructs.
    threadCount = parallelTeam.getThreadCount();
    sharedInteractions = new SharedInteger();
    sharedEnergy = new SharedDouble();
    doLongRangeCorrection = forceField.getBoolean("VDW_CORRECTION", false);
    vanDerWaalsRegion = new VanDerWaalsRegion();
    initializationTime = new long[threadCount];
    vdwTime = new long[threadCount];
    reductionTime = new long[threadCount];

    // Define how force arrays will be accumulated.
    atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
    String value = forceField.getString("ARRAY_REDUCTION", "MULTI");
    try {
      atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(
          format(
              " Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value, atomicDoubleArrayImpl));
    }

    // Allocate coordinate arrays and set up reduction indices and values.
    initAtomArrays();

    /*
     Define the multiplicative switch, which sets vdW energy to zero
     at the cutoff distance using a window that begin at 90% of the
     cutoff distance.
    */
    double vdwTaper = 0.9 * vdwCutoff;
    double buff = 2.0;
    nonbondedCutoff = new NonbondedCutoff(vdwCutoff, vdwTaper, buff);
    multiplicativeSwitch = new MultiplicativeSwitch(vdwTaper, vdwCutoff);
    neighborList =
        new NeighborList(null, this.crystal, atoms, neighborListCutoff, buff, parallelTeam);
    pairwiseSchedule = neighborList.getPairwiseSchedule();
    neighborLists = new int[nSymm][][];

    // Reduce and expand the coordinates of the asymmetric unit. Then build the first neighbor-list.
    buildNeighborList(atoms);

    // Then, optionally, prevent that neighbor list from ever updating.
    neighborList.setDisableUpdates(forceField.getBoolean("DISABLE_NEIGHBOR_UPDATES", false));

    logger.info(toString());
  }

  /**
   * {@inheritDoc}
   *
   * <p>Apply masking rules for 1-2, 1-3 and 1-4 interactions.
   */
  @Override
  public void applyMask(final int i, final boolean[] is14, final double[]... masks) {
    double[] mask = masks[0];
    var m12 = mask12[i];
    for (int value : m12) {
      mask[value] = vdwForm.scale12;
    }
    var m13 = mask13[i];
    for (int value : m13) {
      mask[value] = vdwForm.scale13;
    }
    var m14 = mask14[i];
    for (int value : m14) {
      mask[value] = vdwForm.scale14;
      is14[value] = true;
    }
  }

  /**
   * attachExtendedSystem.
   *
   * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public void attachExtendedSystem(ExtendedSystem system) {
    if (system == null) {
      logger.severe("Tried to attach null extended system.");
    }
    esvTerm = true;
    esvSystem = system;
    numESVs = esvSystem.getExtendedResidueList().size();

    // Launch shared lambda/esvLambda initializers if missed (ie. !lambdaTerm) in constructor.
    vdwLambdaAlpha = forceField.getDouble("VDW_LAMBDA_ALPHA", 0.05);
    vdwLambdaExponent = forceField.getDouble("VDW_LAMBDA_EXPONENT", 1.0);
    if (vdwLambdaExponent != 1.0) {
      logger.warning(
          format(
              "ESVs are compatible only with a vdwLambdaExponent of unity!"
                  + " (found %g, resetting to 1.0)",
              vdwLambdaExponent));
      vdwLambdaExponent = 1.0;
    }
    if (vdwLambdaAlpha < 0.0) {
      vdwLambdaAlpha = 0.05;
    }

    intermolecularSoftcore = forceField.getBoolean("INTERMOLECULAR_SOFTCORE", false);
    intramolecularSoftcore = forceField.getBoolean("INTRAMOLECULAR_SOFTCORE", false);

    previousAtoms = atoms;
    Atom[] atomsExt = esvSystem.getExtendedAtoms();
    int[] moleculeExt = esvSystem.getExtendedMolecule();
    setAtoms(atomsExt, moleculeExt);
  }

  /**
   * destroy.
   *
   * @throws java.lang.Exception if any.
   */
  public void destroy() throws Exception {
    if (neighborList != null) {
      neighborList.destroy();
    }
  }

  /** detachExtendedSystem. */
  public void detachExtendedSystem() {
    setAtoms(previousAtoms, molecule);
    esvTerm = false;
    esvSystem = null;
    numESVs = 0;
    initSoftCore(); // To remove entries from isSoft[] that were due only to ESVs.
  }

  /**
   * The energy routine may be called repeatedly.
   *
   * @param gradient If true, gradients with respect to atomic coordinates are computed.
   * @param print If true, there is verbose printing.
   * @return The energy.
   * @since 1.0
   */
  public double energy(boolean gradient, boolean print) {
    this.gradient = gradient;

    try {
      parallelTeam.execute(vanDerWaalsRegion);
    } catch (Exception e) {
      String message = " Fatal exception expanding coordinates.\n";
      logger.log(Level.SEVERE, message, e);
    }
    return sharedEnergy.get();
  }

  /**
   * getAlpha.
   *
   * @return a double.
   */
  public double getAlpha() {
    return vdwLambdaAlpha;
  }

  /**
   * getBeta.
   *
   * @return a double.
   */
  public double getBeta() {
    return vdwLambdaExponent;
  }

  /**
   * Get the buffer size.
   *
   * @return The buffer.
   * @since 1.0
   */
  public double getBuffer() {
    return nonbondedCutoff.buff;
  }

  /**
   * Return use of the long-range vdW correction.
   *
   * @return True if it is on.
   */
  public boolean getDoLongRangeCorrection() {
    return doLongRangeCorrection;
  }

  /**
   * Get the total Van der Waals potential energy.
   *
   * @return The energy.
   * @since 1.0
   */
  public double getEnergy() {
    return sharedEnergy.get();
  }

  /**
   * Get the number of interacting pairs.
   *
   * @return The interaction count.
   * @since 1.0
   */
  public int getInteractions() {
    return sharedInteractions.get();
  }

  /** {@inheritDoc} */
  @Override
  public double getLambda() {
    return lambda;
  }

  /** {@inheritDoc} */
  @Override
  public void setLambda(double lambda) {
    assert (lambda >= 0.0 && lambda <= 1.0);
    if (!lambdaTerm) {
      return;
    }

    this.lambda = lambda;
    sc1 = vdwLambdaAlpha * (1.0 - lambda) * (1.0 - lambda);
    dsc1dL = -2.0 * vdwLambdaAlpha * (1.0 - lambda);
    d2sc1dL2 = 2.0 * vdwLambdaAlpha;
    sc2 = pow(lambda, vdwLambdaExponent);
    dsc2dL = vdwLambdaExponent * pow(lambda, vdwLambdaExponent - 1.0);
    if (vdwLambdaExponent >= 2.0) {
      d2sc2dL2 =
          vdwLambdaExponent * (vdwLambdaExponent - 1.0) * pow(lambda, vdwLambdaExponent - 2.0);
    } else {
      d2sc2dL2 = 0.0;
    }

    // If LambdaFactors are in OST mode, update them now.
    if (!esvTerm) {
      for (LambdaFactors lf : lambdaFactors) {
        lf.setFactors();
      }
    }

    initSoftCore();

    // Redo the long range correction.
    if (doLongRangeCorrection) {
      longRangeCorrection = getLongRangeCorrection();
      logger.info(format(" Long-range VdW correction %12.8f (kcal/mole).", longRangeCorrection));
    } else {
      longRangeCorrection = 0.0;
    }
  }

  /**
   * Getter for the field <code>bondMask</code>.
   *
   * @return an array of {@link int} objects.
   */
  public int[][] getMask12() {
    return mask12;
  }

  /**
   * Getter for the field <code>angleMask</code>.
   *
   * @return an array of {@link int} objects.
   */
  public int[][] getMask13() {
    return mask13;
  }

  /**
   * Getter for the field <code>torsionMask</code>.
   *
   * @return an array of {@link int} objects.
   */
  public int[][] getMask14() {
    return mask14;
  }

  /**
   * Allow sharing the of the VanDerWaals NeighborList with ParticleMeshEwald.
   *
   * @return The NeighborList.
   */
  public NeighborList getNeighborList() {
    return neighborList;
  }

  /**
   * Getter for the field <code>neighborLists</code>.
   *
   * @return an array of int.
   */
  public int[][][] getNeighborLists() {
    return neighborLists;
  }

  /**
   * Get details of the non-bonded cutoff.
   *
   * @return a {@link ffx.potential.nonbonded.NonbondedCutoff} object.
   */
  public NonbondedCutoff getNonbondedCutoff() {
    return nonbondedCutoff;
  }

  /**
   * Getter for the field <code>pairwiseSchedule</code>.
   *
   * @return a {@link edu.rit.pj.IntegerSchedule} object.
   */
  public IntegerSchedule getPairwiseSchedule() {
    return pairwiseSchedule;
  }

  /**
   * Get the reduction index.
   *
   * @return an array of {@link int} objects.
   */
  public int[] getReductionIndex() {
    return reductionIndex;
  }

  /**
   * getVDWForm.
   *
   * @return a {@link ffx.potential.nonbonded.VanDerWaalsForm} object.
   */
  public VanDerWaalsForm getVDWForm() {
    return vdwForm;
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    if (sharedd2EdL2 == null || !lambdaTerm) {
      return 0.0;
    }
    return sharedd2EdL2.get();
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {
    if (shareddEdL == null || !lambdaTerm) {
      return 0.0;
    }
    return shareddEdL.get();
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] lambdaGradient) {
    if (lambdaGrad == null || !lambdaTerm) {
      return;
    }
    int index = 0;
    for (int i = 0; i < nAtoms; i++) {
      if (atoms[i].isActive()) {
        lambdaGradient[index++] += lambdaGrad.getX(i);
        lambdaGradient[index++] += lambdaGrad.getY(i);
        lambdaGradient[index++] += lambdaGrad.getZ(i);
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Remove the masking rules for 1-2, 1-3 and 1-4 interactions.
   */
  @Override
  public void removeMask(final int i, final boolean[] is14, final double[]... masks) {
    double[] mask = masks[0];
    var m12 = mask12[i];
    for (int value : m12) {
      mask[value] = 1.0;
    }
    var m13 = mask13[i];
    for (int value : m13) {
      mask[value] = 1.0;
    }
    var m14 = mask14[i];
    for (int value : m14) {
      mask[value] = 1.0;
      is14[value] = false;
    }
  }

  /**
   * Setter for the field <code>atoms</code>.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param molecule an array of {@link int} objects.
   */
  public void setAtoms(Atom[] atoms, int[] molecule) {
    this.atoms = atoms;
    this.nAtoms = atoms.length;
    this.molecule = molecule;

    if (nAtoms != molecule.length) {
      logger.warning("Atom and molecule arrays are of different lengths.");
      throw new IllegalArgumentException();
    }
    initAtomArrays();
    buildNeighborList(atoms);
  }

  /**
   * If the crystal being passed in is not equal to the current crystal, then some Van der Waals data
   * structures may need to updated. If <code>nSymm</code> has changed, update arrays dimensioned by
   * nSymm. Finally, rebuild the neighbor-lists.
   *
   * @param crystal The new crystal instance defining the symmetry and boundary conditions.
   */
  public void setCrystal(Crystal crystal) {
    this.crystal = crystal;
    int newNSymm = crystal.spaceGroup.getNumberOfSymOps();
    if (nSymm != newNSymm) {
      nSymm = newNSymm;

      // Allocate memory if necessary.
      if (reduced == null || reduced.length < nSymm) {
        reduced = new double[nSymm][nAtoms * 3];
        reducedXYZ = reduced[0];
        neighborLists = new int[nSymm][][];
      }
    }
    neighborList.setCrystal(crystal);
    neighborListOnly = true;
    try {
      parallelTeam.execute(vanDerWaalsRegion);
    } catch (Exception e) {
      String message = " Fatal exception expanding coordinates.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Setter for the field <code>intermolecularSoftcore</code>.
   *
   * @param intermolecularSoftcore a boolean.
   */
  public void setIntermolecularSoftcore(boolean intermolecularSoftcore) {
    if (!(lambdaTerm || esvTerm)) {
      logger.warning("Illegal softcoring.");
      throw new IllegalArgumentException();
    }
    this.intermolecularSoftcore = intermolecularSoftcore;
  }

  /**
   * Setter for the field <code>intramolecularSoftcore</code>.
   *
   * @param intramolecularSoftcore a boolean.
   */
  public void setIntramolecularSoftcore(boolean intramolecularSoftcore) {
    if (!(lambdaTerm || esvTerm)) {
      logger.warning("Illegal softcoring.");
      throw new IllegalArgumentException();
    }
    this.intramolecularSoftcore = intramolecularSoftcore;
  }

  /**
   * Setter for the field <code>resolution</code>.
   *
   * @param resolution a {@link ffx.potential.bonded.Atom.Resolution} object.
   */
  public void setResolution(Resolution resolution) {
    this.resolution = resolution;
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer("\n  Van der Waals\n");
    sb.append(
        format(
            "   Switch Start:                         %6.3f (A)\n",
            multiplicativeSwitch.getSwitchStart()));
    sb.append(
        format(
            "   Cut-Off:                              %6.3f (A)\n",
            multiplicativeSwitch.getSwitchEnd()));
    sb.append(format("   Long-Range Correction:                %b\n", doLongRangeCorrection));
    if (!reducedHydrogens) {
      sb.append(format("   Reduce Hydrogens:                     %b\n", reducedHydrogens));
    }
    if (lambdaTerm) {
      sb.append("   Alchemical Parameters\n");
      sb.append(format("    Softcore Alpha:                       %5.3f\n", vdwLambdaAlpha));
      sb.append(format("    Lambda Exponent:                      %5.3f\n", vdwLambdaExponent));
    }
    return sb.toString();
  }

  /** Allocate coordinate arrays and set up reduction indices and values. */
  private void initAtomArrays() {
    if (esvTerm) {
      atoms = esvSystem.getExtendedAtoms();
      nAtoms = atoms.length;
    }
    if (atomClass == null || nAtoms > atomClass.length || lambdaTerm || esvTerm) {
      atomClass = new int[nAtoms];
      coordinates = new double[nAtoms * 3];
      reduced = new double[nSymm][nAtoms * 3];
      reducedXYZ = reduced[0];
      reductionIndex = new int[nAtoms];
      reductionValue = new double[nAtoms];
      mask12 = new int[nAtoms][];
      mask13 = new int[nAtoms][];
      mask14 = new int[nAtoms][];
      use = new boolean[nAtoms];
      isSoft = new boolean[nAtoms];
      softCore = new boolean[2][nAtoms];
      grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
      if (lambdaTerm) {
        lambdaGrad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
      } else {
        lambdaGrad = null;
      }
    }

    // Initialize all atoms to be used in the energy.
    fill(use, true);
    fill(isSoft, false);
    fill(softCore[HARD], false);
    fill(softCore[SOFT], false);

    lambdaFactors = new LambdaFactors[threadCount];
    for (int i = 0; i < threadCount; i++) {
     if (lambdaTerm) {
        lambdaFactors[i] = new LambdaFactorsOST();
      } else {
        lambdaFactors[i] = new LambdaFactors();
      }
    }

    String vdwIndex = forceField.getString("VDWINDEX", "Class");

    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];
      assert (i == ai.getXyzIndex() - 1);
      double[] xyz = ai.getXYZ(null);
      int i3 = i * 3;
      coordinates[i3 + XX] = xyz[XX];
      coordinates[i3 + YY] = xyz[YY];
      coordinates[i3 + ZZ] = xyz[ZZ];

      VDWType vdwType = ai.getVDWType();
      if (vdwType == null) {
        // Find the vdW parameters from the AtomType.
        AtomType atomType = ai.getAtomType();
        if (atomType == null) {
          logger.severe(ai.toString());
          return;
        }

        if (vdwIndex.equalsIgnoreCase("Type")) {
          atomClass[i] = atomType.type;
        } else {
          atomClass[i] = atomType.atomClass;
        }
        vdwType = forceField.getVDWType(Integer.toString(atomClass[i]));

        if (vdwType == null) {
          logger.info(" No VdW type for atom class " + atomClass[i]);
          logger.severe(" No VdW type for atom " + ai);
          return;
        }
        ai.setVDWType(vdwType);
      }

      atomClass[i] = vdwType.atomClass;

      List<Bond> bonds = ai.getBonds();
      int numBonds = bonds.size();
      if (reducedHydrogens && vdwType.reductionFactor > 0.0 && numBonds == 1) {
        Bond bond = bonds.get(0);
        Atom heavyAtom = bond.get1_2(ai);
        // Atom indexes start at 1
        reductionIndex[i] = heavyAtom.getIndex() - 1;
        reductionValue[i] = vdwType.reductionFactor;
      } else {
        reductionIndex[i] = i;
        reductionValue[i] = 0.0;
      }

      // Collect 1-2 interactions.
      List<Atom> n12 = ai.get12List();
      mask12[i] = new int[n12.size()];
      int j = 0;
      for (Atom a12 : n12) {
        mask12[i][j++] = a12.getIndex() - 1;
      }

      // Collect 1-3 interactions.
      List<Atom> n13 = ai.get13List();
      mask13[i] = new int[n13.size()];
      j = 0;
      for (Atom a13 : n13) {
        mask13[i][j++] = a13.getIndex() - 1;
      }

      // Collect 1-4 interactions.
      List<Atom> n14 = ai.get14List();
      mask14[i] = new int[n14.size()];
      j = 0;
      for (Atom a14 : n14) {
        mask14[i][j++] = a14.getIndex() - 1;
      }
    }
  }

  /**
   * buildNeighborList.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   */
  private void buildNeighborList(Atom[] atoms) {
    neighborList.setAtoms(atoms);
    if (esvTerm) {
      neighborList.buildList(reduced, neighborLists, null, neighborListOnly, true);
    } else {
      neighborListOnly = true;
      try {
        parallelTeam.execute(vanDerWaalsRegion);
      } catch (Exception e) {
        String message = " Fatal exception expanding coordinates.\n";
        logger.log(Level.SEVERE, message, e);
      }
      neighborListOnly = false;
    }
  }

  /**
   * Computes the long range van der Waals correction to the energy via numerical integration
   *
   * @return Long range correction (kcal/mol)
   * @see "M. P. Allen and D. J. Tildesley, "Computer Simulation of Liquids, 2nd Ed.", Oxford
   *     University Press, 2017, Section 2.8"
   */
  private double getLongRangeCorrection() {
    // Need to treat esvLambda chain terms below before you can do this.
    if (esvTerm) {
      throw new UnsupportedOperationException();
    }

    // Only applicable under periodic boundary conditions.
    if (crystal.aperiodic()) {
      return 0.0;
    }

    /*
     Integrate to maxR = 100 Angstroms or ~33 sigma.
     Integration step size of delR to be 0.01 Angstroms.
    */
    int step = 2;
    double range = 100.0;
    int nDelta = (int) ((double) step * (range - nonbondedCutoff.cut));
    double rDelta = (range - nonbondedCutoff.cut) / (double) nDelta;
    double offset = nonbondedCutoff.cut - 0.5 * rDelta;
    double oneMinLambda = 1.0 - lambda;

    // Count the number of classes and their frequencies
    int maxClass = vdwForm.maxClass;
    int[] radCount = new int[maxClass + 1];
    int[] softRadCount = new int[maxClass + 1];
    fill(radCount, 0);
    fill(softRadCount, 0);
    for (int i = 0; i < nAtoms; i++) {
      radCount[atomClass[i]]++;
      if (isSoft[i]) {
        softRadCount[atomClass[i]]++;
      }
    }

    double total = 0.0;
    // Loop over vdW classes.
    for (int i = 1; i < maxClass + 1; i++) {
      for (int j = i; j < maxClass + 1; j++) {
        if (radCount[i] == 0 || radCount[j] == 0) {
          continue;
        }
        double irv = vdwForm.getCombinedInverseRmin(i, j);
        double ev = vdwForm.getCombinedEps(i, j);
        if (isNaN(irv) || irv == 0 || isNaN(ev)) {
          continue;
        }
        double sume = 0.0;
        for (int k = 1; k <= nDelta; k++) {
          double r = offset + (double) k * rDelta;
          double r2 = r * r;
          final double rho = r * irv;
          final double rho3 = rho * rho * rho;
          final double rhod = rho + vdwForm.delta;
          final double rhod3 = rhod * rhod * rhod;
          double t1 = 0, t2 = 0;
          switch (vdwForm.vdwType) {
            case BUFFERED_14_7:
              final double rho7 = rho3 * rho3 * rho;
              final double rhod7 = rhod3 * rhod3 * rhod;
              t1 = vdwForm.t1n / rhod7;
              t2 = vdwForm.gamma1 / (rho7 + vdwForm.gamma) - 2.0;
              break;
            case LENNARD_JONES:
              final double rho6 = rho3 * rho3;
              final double rhod6 = rhod3 * rhod3;
              t1 = vdwForm.t1n / rhod6;
              t2 = vdwForm.gamma1 / (rho6 + vdwForm.gamma) - 2.0;
              break;
          }
          final double eij = ev * t1 * t2;
          /*
           Apply one minus the multiplicative switch if the
           interaction distance is less than the end of the
           switching window.
          */
          double taper = 1.0;
          if (r2 < nonbondedCutoff.off2) {
            double r3 = r * r2;
            double r4 = r2 * r2;
            double r5 = r2 * r3;
            taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
            taper = 1.0 - taper;
          }

          double jacobian = 4.0 * PI * r2;
          double e = jacobian * eij * taper;
          if (j != i) {
            sume += e;
          } else {
            sume += 0.5 * e;
          }
        }
        double trapezoid = rDelta * sume;

        // Normal correction
        total += radCount[i] * radCount[j] * trapezoid;
        // Correct for softCore vdW that are being turned off.
        if (lambda < 1.0) {
          total -=
              (softRadCount[i] * radCount[j] + (radCount[i] - softRadCount[i]) * softRadCount[j])
                  * oneMinLambda
                  * trapezoid;
        }
      }
    }

    // Divide by the volume of the asymmetric unit.
    Crystal unitCell = crystal.getUnitCell();
    total = total / (unitCell.volume / unitCell.spaceGroup.getNumberOfSymOps());
    return total;
  }

  /**
   * Log the Van der Waals interaction.
   *
   * @param i Atom i.
   * @param k Atom j.
   * @param minr The minimum vdW separation distance.
   * @param r The distance rij.
   * @param eij The interaction energy.
   * @since 1.0
   */
  private void log(int i, int k, double minr, double r, double eij) {
    logger.info(
        format(
            "VDW %6d-%s %d %6d-%s %d %10.4f  %10.4f  %10.4f",
            atoms[i].getIndex(),
            atoms[i].getAtomType().name,
            atoms[i].getVDWType().atomClass,
            atoms[k].getIndex(),
            atoms[k].getAtomType().name,
            atoms[k].getVDWType().atomClass,
            minr, r, eij));
  }

  private void initSoftCore() {

    boolean rebuild = false;

    for (int i = 0; i < nAtoms; i++) {
      boolean soft = atoms[i].applyLambda();
      if (soft != isSoft[i]) {
        isSoft[i] = soft;
        rebuild = true;
      }
    }

    if (!rebuild) {
      return;
    }

    // Initialize the softcore atom masks.
    for (int i = 0; i < nAtoms; i++) {
      if (isSoft[i]) {
        // Outer loop atom hard, inner loop atom soft.
        softCore[HARD][i] = true;
        // Both soft: full intramolecular ligand interactions.
        softCore[SOFT][i] = false;
      } else {
        // Both hard: full interaction between atoms.
        softCore[HARD][i] = false;
        // Outer loop atom soft, inner loop atom hard.
        softCore[SOFT][i] = true;
      }
    }
  }

  /** Test if both atoms match the set Resolution (or true when unset). */
  private boolean include(Atom atom1, Atom atom2) {
    return ((resolution == null)
        || (atom1.getResolution() == resolution && atom2.getResolution() == resolution));
  }

  /**
   * The trick: The setFactors(i,k) method is called every time through the inner VdW loop, avoiding
   * an "if (esv)" branch statement. A plain OST run will have an object of type LambdaFactorsOST
   * instead, which contains an empty version of setFactors(i,k). The OST version sets new factors
   * only on lambda updates, in setLambda().
   */
  public class LambdaFactors {

    double sc1 = 0.0;
    double dsc1dL = 0.0;
    double d2sc1dL2 = 0.0;
    double sc2 = 1.0;
    double dsc2dL = 0.0;
    double d2sc2dL2 = 0.0;

    /** Overriden by the OST version which updates only during setLambda(). */
    public void setFactors() {
    }
  }

  public class LambdaFactorsOST extends LambdaFactors {

    @Override
    public void setFactors() {
      sc1 = VanDerWaals.this.sc1;
      dsc1dL = VanDerWaals.this.dsc1dL;
      d2sc1dL2 = VanDerWaals.this.d2sc1dL2;
      sc2 = VanDerWaals.this.sc2;
      dsc2dL = VanDerWaals.this.dsc2dL;
      d2sc2dL2 = VanDerWaals.this.d2sc2dL2;
    }
  }

  private class VanDerWaalsRegion extends ParallelRegion {

    private final InitializationLoop[] initializationLoop;
    private final ExpandLoop[] expandLoop;
    private final VanDerWaalsLoop[] vanDerWaalsLoop;
    private final ReductionLoop[] reductionLoop;

    VanDerWaalsRegion() {
      initializationLoop = new InitializationLoop[threadCount];
      expandLoop = new ExpandLoop[threadCount];
      vanDerWaalsLoop = new VanDerWaalsLoop[threadCount];
      reductionLoop = new ReductionLoop[threadCount];
    }

    @Override
    public void finish() {
      neighborListOnly = false;
    }

    @Override
    public void run() throws Exception {
      int threadIndex = getThreadIndex();

      // Locally initialize the Loops to help with NUMA?
      if (initializationLoop[threadIndex] == null) {
        initializationLoop[threadIndex] = new InitializationLoop();
        expandLoop[threadIndex] = new ExpandLoop();
        vanDerWaalsLoop[threadIndex] = new VanDerWaalsLoop();
        reductionLoop[threadIndex] = new ReductionLoop();
      }

      // Initialize and expand coordinates.
      try {
        if (threadIndex == 0) {
          initializationTotal = -System.nanoTime();
        }
        execute(0, nAtoms - 1, initializationLoop[threadIndex]);
        execute(0, nAtoms - 1, expandLoop[threadIndex]);
        if (threadIndex == 0) {
          initializationTotal += System.nanoTime();
        }
      } catch (RuntimeException ex) {
        logger.warning("Runtime exception expanding coordinates in thread: " + threadIndex);
        throw ex;
      } catch (Exception e) {
        String message = "Fatal exception expanding coordinates in thread: " + threadIndex + "\n";
        logger.log(Level.SEVERE, message, e);
      }

      // Build the neighbor-list (if necessary) using reduced coordinates.
      if (threadIndex == 0) {
        neighborList.buildList(reduced, neighborLists, null, neighborListOnly, false);
      }
      barrier();

      if (neighborListOnly) {
        return;
      }

      // Compute Van der Waals energy and gradient.
      try {
        if (threadIndex == 0) {
          vdwTotal = -System.nanoTime();
        }
        execute(0, nAtoms - 1, vanDerWaalsLoop[threadIndex]);
        if (threadIndex == 0) {
          vdwTotal += System.nanoTime();
        }
      } catch (RuntimeException ex) {
        logger.warning(
            "Runtime exception evaluating van der Waals energy in thread: " + threadIndex);
        throw ex;
      } catch (Exception e) {
        String message =
            "Fatal exception evaluating van der Waals energy in thread: " + threadIndex + "\n";
        logger.log(Level.SEVERE, message, e);
      }

      // Reduce derivatives.
      if (gradient || lambdaTerm) {
        try {
          if (threadIndex == 0) {
            reductionTotal = -System.nanoTime();
          }
          execute(0, nAtoms - 1, reductionLoop[threadIndex]);
          if (threadIndex == 0) {
            reductionTotal += System.nanoTime();
          }
        } catch (RuntimeException ex) {
          logger.warning(
              "Runtime exception reducing van der Waals gradient in thread: " + threadIndex);
          throw ex;
        } catch (Exception e) {
          String message =
              "Fatal exception reducing van der Waals gradient in thread: " + threadIndex + "\n";
          logger.log(Level.SEVERE, message, e);
        }
      }

      // Log timings.
      if (threadIndex == 0 && logger.isLoggable(Level.FINE)) {
        double total = (initializationTotal + vdwTotal + reductionTotal) * 1e-9;
        logger.fine(format("\n Van der Waals: %7.4f (sec)", total));
        logger.fine(" Thread    Init    Energy  Reduce  Total     Counts");
        long initMax = 0;
        long vdwMax = 0;
        long reductionMax = 0;
        long initMin = Long.MAX_VALUE;
        long vdwMin = Long.MAX_VALUE;
        long reductionMin = Long.MAX_VALUE;
        int countMin = Integer.MAX_VALUE;
        int countMax = 0;
        for (int i = 0; i < threadCount; i++) {
          int count = vanDerWaalsLoop[i].getCount();
          long totalTime = initializationTime[i] + vdwTime[i] + reductionTime[i];
          logger.fine(
              format(
                  "    %3d   %7.4f %7.4f %7.4f %7.4f %10d",
                  i,
                  initializationTime[i] * 1e-9,
                  vdwTime[i] * 1e-9,
                  reductionTime[i] * 1e-9,
                  totalTime * 1e-9,
                  count));
          initMax = max(initializationTime[i], initMax);
          vdwMax = max(vdwTime[i], vdwMax);
          reductionMax = max(reductionTime[i], reductionMax);
          countMax = max(countMax, count);
          initMin = min(initializationTime[i], initMin);
          vdwMin = min(vdwTime[i], vdwMin);
          reductionMin = min(reductionTime[i], reductionMin);
          countMin = min(countMin, count);
        }
        long totalMin = initMin + vdwMin + reductionMin;
        long totalMax = initMax + vdwMax + reductionMax;
        long totalActual = initializationTotal + vdwTotal + reductionTotal;
        logger.fine(
            format(
                " Min      %7.4f %7.4f %7.4f %7.4f %10d",
                initMin * 1e-9, vdwMin * 1e-9, reductionMin * 1e-9, totalMin * 1e-9, countMin));
        logger.fine(
            format(
                " Max      %7.4f %7.4f %7.4f %7.4f %10d",
                initMax * 1e-9, vdwMax * 1e-9, reductionMax * 1e-9, totalMax * 1e-9, countMax));
        logger.fine(
            format(
                " Delta    %7.4f %7.4f %7.4f %7.4f %10d",
                (initMax - initMin) * 1e-9,
                (vdwMax - vdwMin) * 1e-9,
                (reductionMax - reductionMin) * 1e-9,
                (totalMax - totalMin) * 1e-9,
                (countMax - countMin)));
        logger.fine(
            format(
                " Actual   %7.4f %7.4f %7.4f %7.4f %10d\n",
                initializationTotal * 1e-9,
                vdwTotal * 1e-9,
                reductionTotal * 1e-9,
                totalActual * 1e-9,
                sharedInteractions.get()));
      }
    }

    /**
     * {@inheritDoc}
     *
     * <p>This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void start() {

      // Initialize the shared variables.
      if (doLongRangeCorrection) {
        longRangeCorrection = getLongRangeCorrection();
        sharedEnergy.set(longRangeCorrection);
      } else {
        sharedEnergy.set(0.0);
      }
      sharedInteractions.set(0);
      if (lambdaTerm) {
        shareddEdL.set(0.0);
        sharedd2EdL2.set(0.0);
      }
      if (esvTerm) {
        esvSystem.initEsvVdw();
        lambdaFactors = new LambdaFactors[threadCount];
        for (int i = 0; i < threadCount; i++) {
          lambdaFactors[i] = new LambdaFactors();
        }
      }

      grad.alloc(nAtoms);
      if (lambdaTerm) {
        lambdaGrad.alloc(nAtoms);
      }
    }

    /** Update the local coordinate array and initialize reduction variables. */
    private class InitializationLoop extends IntegerForLoop {

      private int threadID;

      @Override
      public void run(int lb, int ub) {
        for (int i = lb, i3 = 3 * lb; i <= ub; i++, i3 += 3) {
          Atom atom = atoms[i];
          coordinates[i3 + XX] = atom.getX();
          coordinates[i3 + YY] = atom.getY();
          coordinates[i3 + ZZ] = atom.getZ();
          use[i] = atom.getUse();

          VDWType vdwType = atom.getVDWType();

          if (vdwType == null) {
            logger.info(" No VdW type for atom " + atom);
            logger.severe(" No VdW type for atom " + atom);
            return;
          }

          atomClass[i] = vdwType.atomClass;

          // Set reduction values.
          List<Bond> bonds = atom.getBonds();
          int numBonds = bonds.size();
          if (reducedHydrogens && vdwType.reductionFactor > 0.0 && numBonds == 1) {
            Bond bond = bonds.get(0);
            Atom heavyAtom = bond.get1_2(atom);
            // Atom indexes start at 1
            reductionIndex[i] = heavyAtom.getIndex() - 1;
            reductionValue[i] = vdwType.reductionFactor;
          } else {
            reductionIndex[i] = i;
            reductionValue[i] = 0.0;
          }

          // Collect 1-2 interactions.
          List<Atom> n12 = atom.get12List();
          mask12[i] = new int[n12.size()];
          int j = 0;
          for (Atom a12 : n12) {
            mask12[i][j++] = a12.getIndex() - 1;
          }

          // Collect 1-3 interactions.
          List<Atom> n13 = atom.get13List();
          mask13[i] = new int[n13.size()];
          j = 0;
          for (Atom a13 : n13) {
            mask13[i][j++] = a13.getIndex() - 1;
          }

          // Collect 1-4 interactions.
          List<Atom> n14 = atom.get14List();
          mask14[i] = new int[n14.size()];
          j = 0;
          for (Atom a14 : n14) {
            mask14[i][j++] = a14.getIndex() - 1;
          }
        }
        if (gradient) {
          grad.reset(threadID, lb, ub);
        }
        if (lambdaTerm) {
          lambdaGrad.reset(threadID, lb, ub);
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
        initializationTime[threadID] = -System.nanoTime();
      }
    }

    private class ExpandLoop extends IntegerForLoop {

      private final double[] in = new double[3];
      private final double[] out = new double[3];
      private int threadID;
      // Extra padding to avert cache interference.
      private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
      private long pad8, pad9, pada, padb, padc, padd, pade, padf;

      @Override
      public void finish() {
        initializationTime[threadID] += System.nanoTime();
      }

      @Override
      public void run(int lb, int ub) {
        /*
         Set up the local coordinate array for the asymmetric unit,
         applying reduction factors to the hydrogen atoms.
        */
        for (int i = lb; i <= ub; i++) {
          int i3 = i * 3;
          int iX = i3 + XX;
          int iY = i3 + YY;
          int iZ = i3 + ZZ;
          double x = coordinates[iX];
          double y = coordinates[iY];
          double z = coordinates[iZ];
          int redIndex = reductionIndex[i];
          if (redIndex >= 0) {
            int r3 = redIndex * 3;
            double rx = coordinates[r3++];
            double ry = coordinates[r3++];
            double rz = coordinates[r3];
            double a = reductionValue[i];
            reducedXYZ[iX] = a * (x - rx) + rx;
            reducedXYZ[iY] = a * (y - ry) + ry;
            reducedXYZ[iZ] = a * (z - rz) + rz;
            double[] rxyz = {reducedXYZ[iX], reducedXYZ[iY], reducedXYZ[iZ]};
            atoms[i].setRedXYZ(rxyz);
          } else {
            reducedXYZ[iX] = x;
            reducedXYZ[iY] = y;
            reducedXYZ[iZ] = z;
          }
        }

        List<SymOp> symOps = crystal.spaceGroup.symOps;

        if (symOps.size() != nSymm) {
          String message =
              format(" Programming Error: nSymm %d != symOps.size %d", nSymm, symOps.size());
          logger.log(Level.WARNING, message);
          logger.log(Level.WARNING, " Replicates\n{0}", crystal.toString());
          logger.log(Level.WARNING, " Unit Cell\n{0}", crystal.getUnitCell().toString());
        }

        double sp2 = crystal.getSpecialPositionCutoff();
        sp2 *= sp2;
        for (int iSymOp = 1; iSymOp < nSymm; iSymOp++) {
          SymOp symOp = symOps.get(iSymOp);
          double[] xyz = reduced[iSymOp];
          for (int i = lb; i <= ub; i++) {
            int i3 = i * 3;
            int iX = i3 + XX;
            int iY = i3 + YY;
            int iZ = i3 + ZZ;
            in[0] = reducedXYZ[iX];
            in[1] = reducedXYZ[iY];
            in[2] = reducedXYZ[iZ];
            crystal.applySymOp(in, out, symOp);
            xyz[iX] = out[0];
            xyz[iY] = out[1];
            xyz[iZ] = out[2];

            // Check if the atom is at a special position.
            double dx = in[0] - out[0];
            double dy = in[1] - out[1];
            double dz = in[2] - out[2];
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < sp2) {
              logger.log(
                  Level.WARNING, " Atom may be at a special position: {0}", atoms[i].toString());
            }
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
      }
    }

    /**
     * The Van der Waals loop class contains methods and thread local variables used to evaluate the
     * Van der Waals energy and gradients with respect to atomic coordinates.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class VanDerWaalsLoop extends IntegerForLoop {

      private final double[] dx_local;
      private final double[][] transOp;
      private int count;
      private double energy;
      private int threadID;
      private double dEdL;
      private double d2EdL2;
      private double[] mask;
      private boolean[] vdw14;
      private LambdaFactors lambdaFactorsLocal;

      // Extra padding to avert cache interference.
      private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
      private long pad8, pad9, pada, padb, padc, padd, pade, padf;

      VanDerWaalsLoop() {
        super();
        dx_local = new double[3];
        transOp = new double[3][3];
      }

      @Override
      public void finish() {
        /*
         Reduce the energy, interaction count and gradients from this
         thread into the shared variables.
        */
        sharedEnergy.addAndGet(energy);
        sharedInteractions.addAndGet(count);
        if (lambdaTerm) {
          shareddEdL.addAndGet(dEdL);
          sharedd2EdL2.addAndGet(d2EdL2);
        }
        vdwTime[threadID] += System.nanoTime();
      }

      public int getCount() {
        return count;
      }

      @Override
      public void run(int lb, int ub) {
        double e = 0.0;
        double[] xyzS = reduced[0];
        // neighborLists array: [nSymm][nAtoms][nNeighbors]
        int[][] list = neighborLists[0];
        double[] esvVdwPrefactori = new double[3];
        double[] esvVdwPrefactork = new double[3];
        for (int i = lb; i <= ub; i++) {
          if (!use[i]) {
            continue;
          }
          Atom atomi = atoms[i];
          // Flag to indicate if atom i is effected by an extended system variable.
          final boolean esvi = esvTerm && esvSystem.isTitratingHydrogen(i);
          if(esvTerm){ esvSystem.getVdwPrefactor(i, esvVdwPrefactori); }
          int i3 = i * 3;
          // Atomic coordinates of atom i, including application of reduction factors.
          final double xi = reducedXYZ[i3++];
          final double yi = reducedXYZ[i3++];
          final double zi = reducedXYZ[i3];
          // If this hydrogen atomic position is reduced, this is index of its heavy atom.
          final int redi = reductionIndex[i];
          // If this hydrogen atomic position is reduced, this is the reduction factor.
          final double redv = reductionValue[i];
          final double rediv = 1.0 - redv;
          final int classI = atomClass[i];
          // Molecule index for atom i.
          final int moleculei = molecule[i];
          // Zero out local gradient accumulation variables.
          double gxi = 0.0;
          double gyi = 0.0;
          double gzi = 0.0;
          // Zero out local gradient accumulation variables for reduced atoms.
          double gxredi = 0.0;
          double gyredi = 0.0;
          double gzredi = 0.0;
          // Zero out local variables for accumulation of dU/dX/dL.
          double lxi = 0.0;
          double lyi = 0.0;
          double lzi = 0.0;
          // Zero out local variables for accumulation of dU/dX/dL for reduced atoms.
          double lxredi = 0.0;
          double lyredi = 0.0;
          double lzredi = 0.0;
          // Collect information about 1-4 interactions, and fill the mask array due to
          // application of 1-2, 1-3 and 1-4 vdW scale factors.
          applyMask(i, vdw14, mask);
          // Default is that the outer loop atom is hard.
          boolean[] softCorei = softCore[HARD];
          if (isSoft[i]) {
            softCorei = softCore[SOFT];
          }
          // Loop over the neighbor list.
          final int[] neighbors = list[i];
          for (final int k : neighbors) {
            Atom atomk = atoms[k];
            // Check if atom k is in use, and if we're including the pairwise interaction.
            if (!use[k] || !include(atomi, atomk)) {
              continue;
            }
            // Flag to indicate if atom k is effected by an extended system variable.
            final boolean esvk = esvTerm && esvSystem.isTitratingHydrogen(k);
            if(esvTerm){esvSystem.getVdwPrefactor(k, esvVdwPrefactork);}
            // Hide these global variable names for thread safety.
            final double sc1, dsc1dL, d2sc1dL2;
            final double sc2, dsc2dL, d2sc2dL2;
            int k3 = k * 3;
            // Atomic coordinates of atom k, including application of reduction factors.
            final double xk = xyzS[k3++];
            final double yk = xyzS[k3++];
            final double zk = xyzS[k3];
            // Compute the separation vector between atom i and k.
            dx_local[0] = xi - xk;
            dx_local[1] = yi - yk;
            dx_local[2] = zi - zk;
            // Apply the minimum image convention (if periodic).
            final double r2 = crystal.image(dx_local);
            int classK = atomClass[k];
            double irv = vdwForm.getCombinedInverseRmin(classI, classK);
            if (vdw14[k]) {
              // Some force fields such as CHARMM have special 1-4 vdW parameters.
              irv = vdwForm.getCombinedInverseRmin14(classI, classK);
            }
            if (r2 <= nonbondedCutoff.off2 && mask[k] > 0 && irv > 0) {
              // Compute the separation distance.
              final double r = sqrt(r2);
              // Check if i and k are part of the same molecule.
              boolean sameMolecule = (moleculei == molecule[k]);
              // This a soft interaction based either 1) the softCore flag or 2) the two ESV flags.
              boolean soft = softCorei[k];
              // If both atoms are softcore, respect the intermolecularSoftcore and
              // intramolecularSoftcore flags.
              if (isSoft[i] && isSoft[k]) {
                if (intermolecularSoftcore && !sameMolecule) {
                  soft = true;
                } else if (intramolecularSoftcore && sameMolecule) {
                  soft = true;
                }
              }
              if (soft) {
                /*
                The setFactors(i,k) method is empty unless ESVs are present.
                If OST lambda present, lambdaFactors will already have been updated during setLambda().
                */
                sc1 = lambdaFactorsLocal.sc1;
                dsc1dL = lambdaFactorsLocal.dsc1dL;
                d2sc1dL2 = lambdaFactorsLocal.d2sc1dL2;
                sc2 = lambdaFactorsLocal.sc2;
                dsc2dL = lambdaFactorsLocal.dsc2dL;
                d2sc2dL2 = lambdaFactorsLocal.d2sc2dL2;
              } else {
                sc1 = 0.0;
                dsc1dL = 0.0;
                d2sc1dL2 = 0.0;
                sc2 = 1.0;
                dsc2dL = 0.0;
                d2sc2dL2 = 0.0;
              }
              final double alpha = sc1;
              final double lambda5 = sc2;
              /*
               Calculate Van der Waals interaction energy.
               Notation of Schnieders et al. The structure,
               thermodynamics, and solubility of organic
               crystals from simulation with a polarizable force
               field. J. Chem. Theory Comput. 8, 1721–1736 (2012).
              */
              double ev = mask[k] * vdwForm.getCombinedEps(classI, classK);
              if (vdw14[k]) {
                ev = mask[k] * vdwForm.getCombinedEps14(classI, classK);
              }
              final double eps_lambda = ev * lambda5;
              final double rho = r * irv;
              final double rhoDisp1 = vdwForm.rhoDisp1(rho);
              final double rhoDisp = rhoDisp1 * rho;
              final double rhoDelta1 = vdwForm.rhoDelta1(rho + vdwForm.delta);
              final double rhoDelta = rhoDelta1 * (rho + vdwForm.delta);
              final double alphaRhoDelta = alpha + rhoDelta;
              final double alphaRhoDispGamma = alpha + rhoDisp + vdwForm.gamma;
              final double t1d = 1.0 / alphaRhoDelta;
              final double t2d = 1.0 / alphaRhoDispGamma;
              final double t1 = vdwForm.t1n * t1d;
              final double t2a = vdwForm.gamma1 * t2d;
              final double t2 = t2a - 2.0;
              double eik = eps_lambda * t1 * t2;
              /*
               Apply a multiplicative switch if the interaction
               distance is greater than the beginning of the taper.
              */
              double taper = 1.0;
              double dtaper = 0.0;
              if (r2 > nonbondedCutoff.cut2) {
                final double r3 = r2 * r;
                final double r4 = r2 * r2;
                final double r5 = r2 * r3;
                taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
              }
              double esvik = 1.0;
              if (esvi || esvk){
                esvik = esvVdwPrefactori[0] * esvVdwPrefactork[0];
              }
              e += eik * taper * esvik;
              count++;
              if (!gradient && !soft) {
                continue;
              }
              final int redk = reductionIndex[k];
              final double red = reductionValue[k];
              final double redkv = 1.0 - red;
              final double dt1d_dr = vdwForm.repDispPower * rhoDelta1 * irv;
              final double dt2d_dr = vdwForm.dispersivePower * rhoDisp1 * irv;
              final double dt1_dr = t1 * dt1d_dr * t1d;
              final double dt2_dr = t2a * dt2d_dr * t2d;
              double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
              final double ir = 1.0 / r;
              final double drdx = dx_local[0] * ir;
              final double drdy = dx_local[1] * ir;
              final double drdz = dx_local[2] * ir;
              if (gradient) {
                final double dswitch = esvik * (eik * dtaper + dedr * taper);
                final double dedx = dswitch * drdx;
                final double dedy = dswitch * drdy;
                final double dedz = dswitch * drdz;
                gxi += dedx * redv;
                gyi += dedy * redv;
                gzi += dedz * redv;
                gxredi += dedx * rediv;
                gyredi += dedy * rediv;
                gzredi += dedz * rediv;
                grad.sub(threadID, k, red * dedx, red * dedy, red * dedz);
                grad.sub(threadID, redk, redkv * dedx, redkv * dedy, redkv * dedz);
                // Assign this gradient to attached ESVs.
                if (esvi) {
                  esvSystem.addVdwDeriv(i, eik * taper, esvVdwPrefactori, esvVdwPrefactork[0]);
                }
                if (esvk) {
                  esvSystem.addVdwDeriv(k, eik * taper, esvVdwPrefactork, esvVdwPrefactori[0]);
                }
              }
              if (gradient && soft) {
                final double dt1 = -t1 * t1d * dsc1dL;
                final double dt2 = -t2a * t2d * dsc1dL;
                final double f1 = dsc2dL * t1 * t2;
                final double f2 = sc2 * dt1 * t2;
                final double f3 = sc2 * t1 * dt2;
                final double dedl = ev * (f1 + f2 + f3);
                dEdL += dedl * taper * esvik;
                final double t1d2 = -dsc1dL * t1d * t1d;
                final double t2d2 = -dsc1dL * t2d * t2d;
                final double d2t1 = -dt1 * t1d * dsc1dL - t1 * t1d * d2sc1dL2 - t1 * t1d2 * dsc1dL;
                final double d2t2 =
                    -dt2 * t2d * dsc1dL - t2a * t2d * d2sc1dL2 - t2a * t2d2 * dsc1dL;
                final double df1 = d2sc2dL2 * t1 * t2 + dsc2dL * dt1 * t2 + dsc2dL * t1 * dt2;
                final double df2 = dsc2dL * dt1 * t2 + sc2 * d2t1 * t2 + sc2 * dt1 * dt2;
                final double df3 = dsc2dL * t1 * dt2 + sc2 * dt1 * dt2 + sc2 * t1 * d2t2;
                final double de2dl2 = ev * (df1 + df2 + df3);
                d2EdL2 += de2dl2 * taper * esvik;
                final double t11 = -dsc2dL * t2 * dt1_dr;
                final double t21 = -dsc2dL * t1 * dt2_dr;
                final double t13 = 2.0 * sc2 * t2 * dt1_dr * dsc1dL * t1d;
                final double t23 = 2.0 * sc2 * t1 * dt2_dr * dsc1dL * t2d;
                final double t12 = -sc2 * dt2 * dt1_dr;
                final double t22 = -sc2 * dt1 * dt2_dr;
                final double dedldr = ev * (t11 + t12 + t13 + t21 + t22 + t23);
                final double dswitch = esvik * (dedl * dtaper + dedldr * taper);
                final double dedldx = dswitch * drdx;
                final double dedldy = dswitch * drdy;
                final double dedldz = dswitch * drdz;
                lxi += dedldx * redv;
                lyi += dedldy * redv;
                lzi += dedldz * redv;
                lxredi += dedldx * rediv;
                lyredi += dedldy * rediv;
                lzredi += dedldz * rediv;
                if (lambdaTerm) {
                  lambdaGrad.sub(threadID, k, red * dedldx, red * dedldy, red * dedldz);
                  lambdaGrad.sub(threadID, redk, redkv * dedldx, redkv * dedldy, redkv * dedldz);
                }

              }
            }
          }
          if (gradient) {
            grad.add(threadID, i, gxi, gyi, gzi);
            grad.add(threadID, redi, gxredi, gyredi, gzredi);
            if (lambdaTerm) {
              lambdaGrad.add(threadID, i, lxi, lyi, lzi);
              lambdaGrad.add(threadID, redi, lxredi, lyredi, lzredi);
            }
          }
          removeMask(i, vdw14, mask);
        }
        energy += e;
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        for (int iSymOp = 1; iSymOp < nSymm; iSymOp++) {
          e = 0.0;
          SymOp symOp = symOps.get(iSymOp);
          // Compute the total transformation operator: R = ToCart * Rot * ToFrac.
          crystal.getTransformationOperator(symOp, transOp);
          xyzS = reduced[iSymOp];
          list = neighborLists[iSymOp];

          for (int i = lb; i <= ub; i++) {
            int i3 = i * 3;
            if (!use[i]) {
              continue;
            }
            Atom atomi = atoms[i];
            final boolean esvi = esvTerm && esvSystem.isTitratingHydrogen(i);
            if(esvTerm){ esvSystem.getVdwPrefactor(i, esvVdwPrefactori); }
            final double xi = reducedXYZ[i3++];
            final double yi = reducedXYZ[i3++];
            final double zi = reducedXYZ[i3];
            final int redi = reductionIndex[i];
            final double redv = reductionValue[i];
            final double rediv = 1.0 - redv;
            final int classI = atomClass[i];
            double gxi = 0.0;
            double gyi = 0.0;
            double gzi = 0.0;
            double gxredi = 0.0;
            double gyredi = 0.0;
            double gzredi = 0.0;
            double lxi = 0.0;
            double lyi = 0.0;
            double lzi = 0.0;
            double lxredi = 0.0;
            double lyredi = 0.0;
            double lzredi = 0.0;
            // Default is that the outer loop atom is hard.
            boolean[] softCorei = softCore[HARD];
            if (isSoft[i]) {
              softCorei = softCore[SOFT];
            }

            // Loop over the neighbor list.
            final int[] neighbors = list[i];
            for (final int k : neighbors) {
              Atom atomk = atoms[k];
              if (!use[k] || !include(atomi, atomk)) {
                continue;
              }
              final boolean esvk = esvTerm && esvSystem.isTitratingHydrogen(k);
              if(esvTerm){ esvSystem.getVdwPrefactor(k, esvVdwPrefactork); }
              // Hide these global variable names for thread safety.
              final double sc1, dsc1dL, d2sc1dL2;
              final double sc2, dsc2dL, d2sc2dL2;
              int k3 = k * 3;
              final double xk = xyzS[k3++];
              final double yk = xyzS[k3++];
              final double zk = xyzS[k3];
              dx_local[0] = xi - xk;
              dx_local[1] = yi - yk;
              dx_local[2] = zi - zk;
              final double r2 = crystal.image(dx_local);
              int classK = atomClass[k];
              final double irv = vdwForm.getCombinedInverseRmin(classI, classK);
              if (r2 <= nonbondedCutoff.off2 && irv > 0) {
                final double selfScale = (i == k) ? 0.5 : 1.0;
                final double r = sqrt(r2);
                boolean soft = isSoft[i] || softCorei[k];
                if (soft) {
                  sc1 = lambdaFactorsLocal.sc1;
                  dsc1dL = lambdaFactorsLocal.dsc1dL;
                  d2sc1dL2 = lambdaFactorsLocal.d2sc1dL2;
                  sc2 = lambdaFactorsLocal.sc2;
                  dsc2dL = lambdaFactorsLocal.dsc2dL;
                  d2sc2dL2 = lambdaFactorsLocal.d2sc2dL2;
                } else {
                  sc1 = 0.0;
                  dsc1dL = 0.0;
                  d2sc1dL2 = 0.0;
                  sc2 = 1.0;
                  dsc2dL = 0.0;
                  d2sc2dL2 = 0.0;
                }
                final double alpha = sc1;
                final double lambdaN = sc2;
                final double ev = vdwForm.getCombinedEps(classI, classK);
                final double eps_lambda = ev * lambdaN;
                final double rho = r * irv;
                final double rhoDisp1 = vdwForm.rhoDisp1(rho);
                final double rhoDisp = rhoDisp1 * rho;
                final double rhoDelta1 = vdwForm.rhoDelta1(rho + vdwForm.delta);
                final double rhoDelta = rhoDelta1 * (rho + vdwForm.delta);
                final double alphaRhoDelta = alpha + rhoDelta;
                final double alphaRhoDispGamma = alpha + rhoDisp + vdwForm.gamma;
                final double t1d = 1.0 / alphaRhoDelta;
                final double t2d = 1.0 / alphaRhoDispGamma;
                final double t1 = vdwForm.t1n * t1d;
                final double t2a = vdwForm.gamma1 * t2d;
                final double t2 = t2a - 2.0;
                double eik = eps_lambda * t1 * t2;
                // Apply a multiplicative switch if the interaction distance is greater than the
                // beginning of the taper.
                double taper = 1.0;
                double dtaper = 0.0;

                if (r2 > nonbondedCutoff.cut2) {
                  final double r3 = r2 * r;
                  final double r4 = r2 * r2;
                  final double r5 = r2 * r3;
                  taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                  dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
                }

                double esvik = 1.0;
                if (esvi || esvk){
                  esvik = esvVdwPrefactori[0] * esvVdwPrefactork[0];
                }
                e += selfScale * eik * taper * esvik;

                count++;
                if (!gradient && !soft) {
                  continue;
                }

                final int redk = reductionIndex[k];
                final double red = reductionValue[k];
                final double redkv = 1.0 - red;
                final double dt1d_dr = vdwForm.repDispPower * rhoDelta1 * irv;
                final double dt2d_dr = vdwForm.dispersivePower * rhoDisp1 * irv;
                final double dt1_dr = t1 * dt1d_dr * t1d;
                final double dt2_dr = t2a * dt2d_dr * t2d;
                double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);

                final double ir = 1.0 / r;
                double drdx = dx_local[0] * ir;
                double drdy = dx_local[1] * ir;
                double drdz = dx_local[2] * ir;
                dedr = esvik * (eik * dtaper + dedr * taper);
                if (gradient) {
                  double dedx = selfScale * dedr * drdx;
                  double dedy = selfScale * dedr * drdy;
                  double dedz = selfScale * dedr * drdz;
                  gxi += dedx * redv;
                  gyi += dedy * redv;
                  gzi += dedz * redv;
                  gxredi += dedx * rediv;
                  gyredi += dedy * rediv;
                  gzredi += dedz * rediv;

                  // Apply the transpose of the transformation operator.
                  final double dedxk =
                      dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                  final double dedyk =
                      dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                  final double dedzk =
                      dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
                  grad.sub(threadID, k, red * dedxk, red * dedyk, red * dedzk);
                  grad.sub(threadID, redk, redkv * dedxk, redkv * dedyk, redkv * dedzk);

                  // Assign this gradient to attached ESVs.
                  if (esvi) {
                    esvSystem.addVdwDeriv(i, selfScale * eik * taper, esvVdwPrefactori, esvVdwPrefactork[0]);
                  }
                  if (esvk) {
                    esvSystem.addVdwDeriv(k, selfScale * eik * taper, esvVdwPrefactork, esvVdwPrefactori[0]);
                  }
                }
                if (gradient && lambdaTerm && soft) {
                  double dt1 = -t1 * t1d * dsc1dL;
                  double dt2 = -t2a * t2d * dsc1dL;
                  double f1 = dsc2dL * t1 * t2;
                  double f2 = sc2 * dt1 * t2;
                  double f3 = sc2 * t1 * dt2;
                  final double dedl = ev * (f1 + f2 + f3);
                  dEdL += selfScale * dedl * taper * esvik;
                  double t1d2 = -dsc1dL * t1d * t1d;
                  double t2d2 = -dsc1dL * t2d * t2d;
                  double d2t1 = -dt1 * t1d * dsc1dL - t1 * t1d * d2sc1dL2 - t1 * t1d2 * dsc1dL;
                  double d2t2 = -dt2 * t2d * dsc1dL - t2a * t2d * d2sc1dL2 - t2a * t2d2 * dsc1dL;
                  double df1 = d2sc2dL2 * t1 * t2 + dsc2dL * dt1 * t2 + dsc2dL * t1 * dt2;
                  double df2 = dsc2dL * dt1 * t2 + sc2 * d2t1 * t2 + sc2 * dt1 * dt2;
                  double df3 = dsc2dL * t1 * dt2 + sc2 * dt1 * dt2 + sc2 * t1 * d2t2;
                  double de2dl2 = ev * (df1 + df2 + df3);
                  d2EdL2 += selfScale * de2dl2 * taper * esvik;
                  double t11 = -dsc2dL * t2 * dt1_dr;
                  double t12 = -sc2 * dt2 * dt1_dr;
                  double t13 = 2.0 * sc2 * t2 * dt1_dr * dsc1dL * t1d;
                  double t21 = -dsc2dL * t1 * dt2_dr;
                  double t22 = -sc2 * dt1 * dt2_dr;
                  double t23 = 2.0 * sc2 * t1 * dt2_dr * dsc1dL * t2d;
                  double dedldr = ev * (t11 + t12 + t13 + t21 + t22 + t23);
                  dedldr = esvik *(dedl * dtaper + dedldr * taper);
                  double dedldx = selfScale * dedldr * drdx;
                  double dedldy = selfScale * dedldr * drdy;
                  double dedldz = selfScale * dedldr * drdz;
                  lxi += dedldx * redv;
                  lyi += dedldy * redv;
                  lzi += dedldz * redv;
                  lxredi += dedldx * rediv;
                  lyredi += dedldy * rediv;
                  lzredi += dedldz * rediv;

                  // Apply the transpose of the transformation operator.
                  final double dedldxk =
                      dedldx * transOp[0][0] + dedldy * transOp[1][0] + dedldz * transOp[2][0];
                  final double dedldyk =
                      dedldx * transOp[0][1] + dedldy * transOp[1][1] + dedldz * transOp[2][1];
                  final double dedldzk =
                      dedldx * transOp[0][2] + dedldy * transOp[1][2] + dedldz * transOp[2][2];
                  lambdaGrad.sub(threadID, k, red * dedldxk, red * dedldyk, red * dedldzk);
                  lambdaGrad.sub(threadID, redk, redkv * dedldxk, redkv * dedldyk, redkv * dedldzk);
                }
              }
            }
            if (gradient) {
              grad.add(threadID, i, gxi, gyi, gzi);
              grad.add(threadID, redi, gxredi, gyredi, gzredi);
              if (lambdaTerm) {
                lambdaGrad.add(threadID, i, lxi, lyi, lzi);
                lambdaGrad.add(threadID, redi, lxredi, lyredi, lzredi);
              }
            }
          }
          energy += e;
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return pairwiseSchedule;
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
        vdwTime[threadID] = -System.nanoTime();
        energy = 0.0;
        count = 0;
        if (lambdaTerm) {
          dEdL = 0.0;
          d2EdL2 = 0.0;
        }
        lambdaFactorsLocal = lambdaFactors[threadID];
        if (lambdaFactorsLocal == null) {
          System.exit(1);
        }

        if (mask == null || mask.length < nAtoms) {
          mask = new double[nAtoms];
          fill(mask, 1.0);
          vdw14 = new boolean[nAtoms];
          fill(vdw14, false);
        }
      }
    }

    /** Reduce Van der Waals gradient. */
    private class ReductionLoop extends IntegerForLoop {

      int threadID;

      @Override
      public void finish() {
        reductionTime[threadID] += System.nanoTime();
      }

      @Override
      public void run(int lb, int ub) {
        if (gradient) {
          grad.reduce(lb, ub);
          for (int i = lb; i <= ub; i++) {
            Atom ai = atoms[i];
            ai.addToXYZGradient(grad.getX(i), grad.getY(i), grad.getZ(i));
          }
        }
        if (lambdaTerm) {
          lambdaGrad.reduce(lb, ub);
        }
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
        reductionTime[threadID] = -System.nanoTime();
      }
    }
  }
}
