package ffx.potential.constraint;

import ffx.numerics.Constraint;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVectorPreservingVisitor;
import org.apache.commons.math3.util.FastMath;

public class CcmaConstraint implements Constraint {
  public static final double DEFAULT_CCMA_NONZERO_CUTOFF = 0.01;
  private static final Logger logger = Logger.getLogger(CcmaConstraint.class.getName());
  private static final int DEFAULT_MAX_ITERS = 150;
  // Might be slightly more elegant to have "component constraint" objects.
  // This is how OpenMM does it, though, and arrays are likely more performative.
  private final int[] atoms1;
  private final int[] atoms2;
  private final double[] lengths;
  private final int nConstraints;
  private final int[] uniqueIndices;
  private final int maxIters = DEFAULT_MAX_ITERS;
  private final double elementCutoff = DEFAULT_CCMA_NONZERO_CUTOFF;
  private final double[] reducedMasses;
  private final RealMatrix kInvSparse;

  /**
   * Constructs a set of bond length Constraints to be satisfied using the Constaint Constraint
   * Matrix Approximation, a parallelizable stable numeric method.
   *
   * <p>The nonzeroCutoff field specifies how large an element of K-1 needs to be to be kept;
   * smaller elements (indicating weak constraint-constraint couplings) are set to zero to make K-1
   * sparse and thus keep CCMA tractable.
   *
   * <p>Due to leaking-this issues, a thin Factory method is provided to build this.
   *
   * <p>The constrainedAngles will constrain both component bonds and the triangle-closing distance.
   * The constrainedBonds list should not include any bonds included in the constrained angles.
   *
   * @param constrainedBonds Bonds to be constrained, not included in constrainedAngles
   * @param constrainedAngles Angles to be constrained indirectly (by constructing a rigid
   *     triangle).
   * @param allAtoms All Atoms of the system, including unconstrained Atoms.
   * @param masses All masses of the system, including unconstrained atom masses.
   * @param nonzeroCutoff CCMA parameter defining how sparse/dense K-1 should be.
   */
  private CcmaConstraint(
      List<Bond> constrainedBonds,
      List<Angle> constrainedAngles,
      final Atom[] allAtoms,
      final double[] masses,
      double nonzeroCutoff) {
    long time = -System.nanoTime();
    int nBonds = constrainedBonds.size();
    int nAngles = constrainedAngles.size();
    assert constrainedAngles.stream()
        .flatMap((Angle a) -> a.getBondList().stream())
        .noneMatch(constrainedBonds::contains);
    nConstraints = nBonds + (3 * nAngles);
    atoms1 = new int[nConstraints];
    atoms2 = new int[nConstraints];
    lengths = new double[nConstraints];

    for (int i = 0; i < nBonds; i++) {
      Bond bi = constrainedBonds.get(i);
      atoms1[i] = bi.getAtom(0).getXyzIndex() - 1;
      atoms2[i] = bi.getAtom(1).getXyzIndex() - 1;
      lengths[i] = bi.bondType.distance;
    }

    /*logger.info(" Atoms 1: " + Arrays.toString(atoms1));
    logger.info(" Atoms 2: " + Arrays.toString(atoms2));
    logger.info(" Lengths: " + Arrays.toString(lengths));*/

    for (int i = 0; i < nAngles; i++) {
      int iAng = nBonds + (3 * i);
      Angle ai = constrainedAngles.get(i);
      Atom center = ai.getCentralAtom();
      Bond b1 = ai.getBond(0);
      Bond b2 = ai.getBond(1);
      Atom at0 = b1.get1_2(center);
      Atom at2 = b2.get1_2(center);

      double angVal = ai.angleType.angle[ai.nh];
      double dist1 = b1.bondType.distance;
      double dist2 = b2.bondType.distance;
      double dist3 = SettleConstraint.lawOfCosines(dist1, dist2, angVal);

      int index0 = at0.getXyzIndex() - 1;
      int index1 = center.getXyzIndex() - 1;
      int index2 = at2.getXyzIndex() - 1;

      atoms1[iAng] = index0;
      atoms2[iAng] = index1;
      lengths[iAng] = dist1;

      ++iAng;
      atoms1[iAng] = index1;
      atoms2[iAng] = index2;
      lengths[iAng] = dist2;

      ++iAng;
      atoms1[iAng] = index0;
      atoms2[iAng] = index2;
      lengths[iAng] = dist3;
    }

    uniqueIndices =
        IntStream.concat(Arrays.stream(atoms1), Arrays.stream(atoms2))
            .sorted()
            .distinct()
            .toArray();

    OpenMapRealMatrix kSparse = new OpenMapRealMatrix(nConstraints, nConstraints);

    int nAtoms = allAtoms.length;

    // Points from an Atom index to the Set of all Constraint indices it's involved in.
    List<Set<Integer>> atomsToConstraints = new ArrayList<>(nAtoms);
    for (int i = 0; i < nAtoms; i++) {
      // Can assume 4 constraints max, since CHONPS is rarely bonded to 5+ atoms.
      atomsToConstraints.add(new HashSet<>(4));
    }

    for (int i = 0; i < nConstraints; i++) {
      atomsToConstraints.get(atoms1[i]).add(i);
      atomsToConstraints.get(atoms2[i]).add(i);
    }

    logger.info(
        String.format(" Initial CCMA setup: %10.6g sec", 1.0E-9 * (time + System.nanoTime())));
    long subTime = -System.nanoTime();

    for (int i = 0; i < nConstraints; i++) {
      int atomi0 = atoms1[i];
      int atomi1 = atoms2[i];
      // DO YOU HAVE ANY IDEA HOW LONG IT TOOK FOR ME TO REALIZE MASSES WERE XYZ-INDEXED?
      double invMassI0 = 1.0 / masses[atomi0 * 3];
      double invMassI1 = 1.0 / masses[atomi1 * 3];
      double sumInv = invMassI0 + invMassI1;

      // Indices of all Constraints that involve either atomi0 or atomi1.
      Set<Integer> coupledConstraints = new HashSet<>(atomsToConstraints.get(atomi0));
      coupledConstraints.addAll(atomsToConstraints.get(atomi1));

      // Iterate over all coupled Constraints, sharing at least one common Atom with constraint i.
      for (int j : coupledConstraints) {
        // Diagonal element, coupling is obviously 1.0.
        if (i == j) {
          kSparse.setEntry(i, j, 1.0);
          continue;
        }

        double scale;
        int atomj0 = atoms1[j];
        int atomj1 = atoms2[j];
        int atoma; // Atom unique to constraint i.
        int atomb; // Atom shared between both constraints.
        int atomc; // Atom unique to constraint j.

        if (atomi0 == atomj0) {
          assert atomi1 != atomj1;
          atoma = atomi1;
          atomb = atomi0;
          atomc = atomj1;
          scale = invMassI0 / sumInv;
        } else if (atomi0 == atomj1) {
          assert atomi1 != atomj0;
          atoma = atomi1;
          atomb = atomi0;
          atomc = atomj0;
          scale = invMassI0 / sumInv;
        } else if (atomi1 == atomj0) {
          // Yes IntelliJ, it should always be true. That's why I'm asserting it.
          assert atomi0 != atomj1;
          atoma = atomi0;
          atomb = atomi1;
          atomc = atomj1;
          scale = invMassI1 / sumInv;
        } else if (atomi1 == atomj1) {
          assert atomi0 != atomj0;
          atoma = atomi1;
          atomb = atomi0;
          atomc = atomj1;
          scale = invMassI1 / sumInv;
        } else {
          throw new IllegalArgumentException(
              " Despite by necessity sharing an atom, these constraints don't share an atom.");
        }

        // We now have a pair of constraints a-b b-c. Find the a-b-c angle.

        // Search for a constraint a-c that closes the triangle.
        boolean foundAngle = false;
        for (int constraintK : atomsToConstraints.get(atoma)) {
          if (atoms1[constraintK] == atomc || atoms2[constraintK] == atomc) {
            double dab = lengths[i];
            double dbc = lengths[j];
            double dac = lengths[constraintK];

            // Apply the Law of Cosines to find the angle defined by the a-b b-c a-c lengths.
            double angle = (dab * dab) + (dbc * dbc) - (dac * dac);
            angle /= (2 * dab * dbc);
            // The angle is formally the cosine of its current value, but all we need is that
            // cosine.
            double coupling = scale * angle;
            kSparse.setEntry(i, j, coupling);
            foundAngle = true;
            break;
          }
        }

        // Otherwise, look for a force field term. In this case, the approximate matrix K
        // will not quite match the true Jacobian matrix J.
        if (!foundAngle) {
          Atom atA = allAtoms[atoma];
          Atom atB = allAtoms[atomb];
          Atom atC = allAtoms[atomc];
          Angle angleB = atA.getAngle(atB, atC);
          double angVal = angleB.angleType.angle[angleB.nh];
          double coupling = scale * FastMath.cos(FastMath.toRadians(angVal));
          kSparse.setEntry(i, j, coupling);
          foundAngle = true;
        }

        if (!foundAngle) {
          logger.severe(
              String.format(" Could not find the angle between coupled constraints %d, %d", i, j));
        }
      }
    }

    subTime += System.nanoTime();
    logger.info(
        String.format(" Time to construct K as a sparse matrix: %10.6g sec", 1.0E-9 * subTime));

    // K is constructed. Invert it.
    // TODO: Find a fast way to invert sparse matrices, since this is at least O(n^2).
    subTime = -System.nanoTime();
    QRDecomposition qrd = new QRDecomposition(kSparse);
    RealMatrix kInvDense = qrd.getSolver().getInverse();
    subTime += System.nanoTime();
    logger.info(String.format(" Time to invert K: %10.6g sec", 1.0E-9 * subTime));
    subTime = -System.nanoTime();

    kInvSparse = new OpenMapRealMatrix(nConstraints, nConstraints);
    IntStream.range(0, nConstraints)
        .
        // parallel(). // I have no idea if OpenMapRealMatrix is thread-safe.
        forEach(
            (int i) -> {
              double[] rowI = kInvDense.getRow(i);
              for (int j = 0; j < nConstraints; j++) {
                double val = rowI[j];
                if (Math.abs(val) > elementCutoff) {
                  kInvSparse.setEntry(i, j, val);
                }
              }
            });

    // TODO: Delete this block.
    double[][] dataDense = kInvDense.getData();
    double[][] dataSparse = kInvSparse.getData();
    logger.fine(" Dense array:");
    for (int i = 0; i < nConstraints; i++) {
      logger.fine(Arrays.toString(dataDense[i]));
    }
    logger.fine(" Sparse array:");
    for (int i = 0; i < nConstraints; i++) {
      logger.fine(Arrays.toString(dataSparse[i]));
    }

    // TODO: Actually do this.
    reducedMasses = new double[nConstraints];
    for (int i = 0; i < nConstraints; i++) {
      int atI = atoms1[i];
      atI *= 3; // The mass array is XYZ-indexed, not atom-indexed.
      int atJ = atoms2[i];
      atJ *= 3;
      double invMassI = 1.0 / masses[atI];
      double invMassJ = 1.0 / masses[atJ];
      reducedMasses[i] = 0.5 / (invMassI + invMassJ);
    }
  }

  /**
   * Constructs a set of bond length Constraints to be satisfied using the Constaint Constraint
   * Matrix Approximation, a parallelizable stable numeric method.
   *
   * <p>The nonzeroCutoff field specifies how large an element of K-1 needs to be to be kept;
   * smaller elements (indicating weak constraint-constraint couplings) are set to zero to make K-1
   * sparse and thus keep CCMA tractable.
   *
   * <p>The constrainedAngles will constrain both component bonds and the triangle-closing distance.
   * The constrainedBonds list should not include any bonds included in the constrained angles.
   *
   * @param constrainedBonds Bonds to be constrained, not included in constrainedAngles
   * @param constrainedAngles Angles to be constrained indirectly (by constructing a rigid
   *     triangle).
   * @param allAtoms All Atoms of the system, including unconstrained Atoms.
   * @param masses All masses of the system, including unconstrained atom masses.
   * @param nonzeroCutoff CCMA parameter defining how sparse/dense K-1 should be.
   * @return Returns a new CcmaConstraint instance.
   */
  public static CcmaConstraint ccmaFactory(
      List<Bond> constrainedBonds,
      List<Angle> constrainedAngles,
      final Atom[] allAtoms,
      final double[] masses,
      double nonzeroCutoff) {
    CcmaConstraint newC =
        new CcmaConstraint(constrainedBonds, constrainedAngles, allAtoms, masses, nonzeroCutoff);
    constrainedBonds.forEach((Bond b) -> b.setConstraint(newC));
    constrainedAngles.forEach((Angle a) -> a.setConstraint(newC));
    return newC;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void applyConstraintToStep(double[] xPrior, double[] xNew, double[] masses, double tol) {
    applyConstraints(xPrior, xNew, masses, false, tol);
  }

  /** {@inheritDoc} */
  @Override
  public void applyConstraintToVelocities(double[] x, double[] v, double[] masses, double tol) {
    applyConstraints(x, v, masses, true, tol);
  }

  /** {@inheritDoc} */
  @Override
  public int[] constrainedAtomIndices() {
    return Arrays.copyOf(uniqueIndices, uniqueIndices.length);
  }

  /** {@inheritDoc} */
  @Override
  public boolean constraintSatisfied(double[] x, double tol) {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public boolean constraintSatisfied(double[] x, double[] v, double xTol, double vTol) {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public int getNumDegreesFrozen() {
    return nConstraints;
  }

  /**
   * Much of the math for applying to coordinates/velocities is the same. As such, OpenMM just uses
   * a single driver method with a flag to indicate velocities or positions.
   *
   * @param xPrior Either pre-step coordinates (positions) or current coordinates (velocities)
   * @param output Output vector, either constrained coordinates or velocities
   * @param masses Masses of all particles
   * @param constrainV Whether to apply constraints to velocities, or to coordinates.
   * @param tol Numerical tolerance.
   */
  private void applyConstraints(
      double[] xPrior, double[] output, double[] masses, boolean constrainV, double tol) {
    if (xPrior == output) {
      throw new IllegalArgumentException(" xPrior and output must be different arrays!");
    }
    long time = -System.nanoTime();

    /*// Distance vector for all constraints. TODO: Flatten to 1D array.
    double[][] r_ij = new double[nConstraints][3];
    // Distance magnitude for all constraints.
    double[] d_ij = new double[nConstraints];
    // Squared distance magnitude for all constraints.
    double[] d_ij2 = new double[nConstraints];

    for (int i = 0; i < nConstraints; i++) {
        int atom1 = 3 * atoms1[i];
        int atom2 = 3 * atoms2[i];
        for (int j = 0; j < 3; j++) {
            r_ij[i][j] = xPrior[atom1 + j] - xPrior[atom2 + j];
        }
        d_ij2[i] = VectorMath.dot(r_ij[i], r_ij[i]);
    }

    double lowerTol = 1 - 2*tol + tol*tol;
    // TODO: Determine if adding tol+tol to upperTol is a typo in OpenMM. Also why they do it in the first place.
    double upperTol = 1 + 2*tol + tol*tol;

    int nConverged = 0;
    double[] constraintDelta = new double[nConstraints];
    double[] tempDelta = new double[nConstraints];

    // Main CCMA loop.
    for (int constraintIter = 0; constraintIter <= maxIters; constraintIter++) {
        if (constraintIter >= maxIters) {
            throw new IllegalArgumentException(String.format(" CCMA constraint failed to converge in %d iterations!", maxIters));
        }
        nConverged = 0;
        double rmsd = 0;
        for (int i = 0; i < nConstraints; i++) {
            int atom1 = atoms1[i] * 3;
            int atom2 = atoms2[i] * 3;

            // Separation vector I-J at this iteration.
            double[] rp_ij = new double[3];
            for (int j = 0; j < 3; j++) {
                rp_ij[j] = output[atom1 + j] - output[atom2 + j];
            }
            if (constrainV) {
                double rrpr = VectorMath.dot(rp_ij, r_ij[i]);
                constraintDelta[i] = -2 * reducedMasses[i] * rrpr / d_ij2[i];
                if (Math.abs(constraintDelta[i]) <= tol) {
                    ++nConverged;
                }
            } else {
                double rp2 = VectorMath.dot(rp_ij, rp_ij);
                double dist2 = lengths[i] * lengths[i];
                double diff = dist2 - rp2;
                double rrpr = VectorMath.dot(rp_ij, r_ij[i]);
                constraintDelta[i] = reducedMasses[i] * diff / rrpr;
                if (rp2 >= lowerTol * dist2 && rp2 <= upperTol * dist2) {
                    ++nConverged;
                }
            }
        }

        // Test if the last iteration satisfied all constraints.
        if (nConverged == nConstraints) {
            break;
        }

        if (kInvSparse != null) {
            for (int i = 0; i < nConstraints; i++) {
                RealVector row = kInvSparse.getRowVector(i);
                RealVectorPreservingVisitor neo = new MatrixWalker(constraintDelta);
                double sum = row.walkInOptimizedOrder(neo);
                tempDelta[i] = sum;
            }
            System.arraycopy(tempDelta, 0, constraintDelta, 0, nConstraints);
        }
        for (int i = 0; i < nConstraints; i++) {
            int atom1 = atoms1[i];
            int atom2 = atoms2[i];
            int at13 = 3 * atom1;
            int at23 = 3 * atom2;

            double[] dr = new double[3];
            VectorMath.scalar(r_ij[i], constraintDelta[i], dr);
            for (int j = 0; j < 3; j++) {
                output[at13 + j] = dr[j] / masses[at13];
                output[at23 + j] = dr[j] / masses[at23];
            }
        }
        logger.info(String.format(" %d converged at iteration %d", nConverged, constraintIter));
    }
    time += System.nanoTime();
    logger.info(String.format(" Application of CCMA constraint: %10.4g sec", (time * 1.0E-9)));*/
  }

  private class MatrixWalker implements RealVectorPreservingVisitor {
    private final double[] constraintDelta; // DO NOT MODIFY.
    private double sum = 0.0;

    MatrixWalker(final double[] cDelta) {
      constraintDelta = cDelta;
    }

    /** {@inheritDoc} */
    @Override
    public double end() {
      return sum;
    }

    /** {@inheritDoc} */
    @Override
    public void start(int dimension, int start, int end) {
    }

    /** {@inheritDoc} */
    @Override
    public void visit(int index, double value) {
      sum += (value * constraintDelta[index]);
    }
  }
}
