// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.bonded;

import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.DoubleMath.scale;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.potential.parameters.BondType.units;

import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.switching.ConstantSwitch;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.parameters.BondType;
import java.util.logging.Logger;

/**
 * RestraintBond class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *     <p>TODO: RestraintBond should extend the Bond class.
 */
public class RestraintBond extends BondedTerm implements LambdaInterface {

  public static final double DEFAULT_RB_LAM_START = 0.75;
  public static final double DEFAULT_RB_LAM_END = 1.0;
  private static final Logger logger = Logger.getLogger(RestraintBond.class.getName());
  private final double restraintLambdaStart;
  private final double restraintLambdaStop;
  private final double restraintLambdaWindow;
  private final double rlwInv;
  private final UnivariateSwitchingFunction switchingFunction;
  public BondType bondType = null;
  private boolean lambdaTerm;
  private double lambda = 1.0;
  private double switchVal = 1.0;
  private double switchdUdL = 1.0;
  private double switchd2UdL2 = 1.0;
  private double dEdL = 0.0;
  private double d2EdL2 = 0.0;
  private double[][] dEdXdL = new double[2][3];
  private Crystal crystal;

  /**
   * Creates a distance restraint between two Atoms.
   *
   * @param a1 First Atom.
   * @param a2 Second Atom.
   * @param crystal Any Crystal used by the system.
   * @param lambdaTerm Whether lambda affects this restraint.
   * @param lamStart At what lambda does the restraint begin to take effect?
   * @param lamEnd At what lambda does the restraint hit full strength?
   * @param sf Switching function determining lambda dependence; null produces a ConstantSwitch.
   */
  public RestraintBond(
      Atom a1,
      Atom a2,
      Crystal crystal,
      boolean lambdaTerm,
      double lamStart,
      double lamEnd,
      UnivariateSwitchingFunction sf) {
    restraintLambdaStart = lamStart;
    restraintLambdaStop = lamEnd;
    assert lamEnd > lamStart;
    restraintLambdaWindow = lamEnd - lamStart;
    rlwInv = 1.0 / restraintLambdaWindow;

    atoms = new Atom[2];

    this.crystal = crystal;

    int i1 = a1.getIndex();
    int i2 = a2.getIndex();
    if (i1 < i2) {
      atoms[0] = a1;
      atoms[1] = a2;
    } else {
      atoms[0] = a2;
      atoms[1] = a1;
    }
    setID_Key(false);
    switchingFunction = (sf == null ? new ConstantSwitch() : sf);
    this.lambdaTerm = lambdaTerm;
    this.setLambda(1.0);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Evaluate this Bond energy.
   */
  @Override
  public double energy(
      boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {

    double[] a0 = new double[3];
    double[] a1 = new double[3];
    // The vector from Atom 1 to Atom 0.
    double[] v10 = new double[3];
    // Gradient on Atoms 0 & 1.
    double[] g0 = new double[3];
    double[] g1 = new double[3];

    atoms[0].getXYZ(a0);
    atoms[1].getXYZ(a1);

    sub(a0, a1, v10);

    if (crystal != null) {
      crystal.image(v10);
    }

    // value is the magnitude of the separation vector
    value = length(v10);
    double dv = value - bondType.distance; // bondType.distance = ideal bond length
    if (bondType.bondFunction.hasFlatBottom()) {
      if (dv > 0) {
        dv = Math.max(0, dv - bondType.flatBottomRadius);
      } else if (dv < 0) {
        dv = Math.min(0, dv + bondType.flatBottomRadius);
      } // Else, no adjustments needed.
    }

    double dv2 = dv * dv;
    double kx2 = units * bondType.forceConstant * dv2 * esvLambda;

    energy = switchVal * kx2;
    dEdL = switchdUdL * kx2;
    d2EdL2 = switchd2UdL2 * kx2;
    double deddt = 2.0 * units * bondType.forceConstant * dv * esvLambda;
    double de = 0.0;

    if (value > 0.0) {
      de = deddt / value;
    }

    scale(v10, switchVal * de, g0);
    scale(v10, -switchVal * de, g1);
    if (gradient) {
      grad.add(threadID, atoms[0].getIndex() - 1, g0[0], g0[1], g0[2]);
      grad.add(threadID, atoms[1].getIndex() - 1, g1[0], g1[1], g1[2]);
    }

    // Remove the factor of rL3
    scale(v10, switchdUdL * de, g0);
    scale(v10, -switchdUdL * de, g1);
    dEdXdL[0][0] = g0[0];
    dEdXdL[0][1] = g0[1];
    dEdXdL[0][2] = g0[2];
    dEdXdL[1][0] = g1[0];
    dEdXdL[1][1] = g1[1];
    dEdXdL[1][2] = g1[2];

    value = dv;
    if (esvTerm) {
      final double esvLambdaInv = (esvLambda != 0.0) ? 1 / esvLambda : 1.0;
      setEsvDeriv(energy * dedesvChain * esvLambdaInv);
    }
    return energy;
  }

  /**
   * Find the other Atom in <b>this</b> Bond. These two atoms are said to be 1-2.
   *
   * @param a The known Atom.
   * @return The other Atom that makes up <b>this</b> Bond, or Null if Atom a is not part of
   *     <b>this</b> Bond.
   */
  public Atom get1_2(Atom a) {
    if (a == atoms[0]) {
      return atoms[1];
    }
    if (a == atoms[1]) {
      return atoms[0];
    }
    return null; // Atom not found in bond
  }

  /**
   * Getter for the field <code>bondType</code>.
   *
   * @return a {@link ffx.potential.parameters.BondType} object.
   */
  public BondType getBondType() {
    return bondType;
  }

  /**
   * Set a reference to the force field parameters.
   *
   * @param bondType a {@link ffx.potential.parameters.BondType} object.
   */
  public void setBondType(BondType bondType) {
    this.bondType = bondType;
  }

  /** {@inheritDoc} */
  @Override
  public double getLambda() {
    return lambda;
  }

  /** {@inheritDoc} */
  @Override
  public void setLambda(double lambda) {
    this.lambda = lambda;
    if (lambdaTerm) {
      if (lambda < restraintLambdaStart) {
        switchVal = 0.0;
        switchdUdL = 0;
        switchd2UdL2 = 0;
      } else if (lambda > restraintLambdaStop) {
        switchVal = 1.0;
        switchdUdL = 0;
        switchd2UdL2 = 0;
      } else {
        double restraintLambda = (lambda - restraintLambdaStart) / restraintLambdaWindow;
        switchVal = switchingFunction.valueAt(restraintLambda);
        switchdUdL = rlwInv * switchingFunction.firstDerivative(restraintLambda);
        switchd2UdL2 = rlwInv * rlwInv * switchingFunction.secondDerivative(restraintLambda);
      }
    } else {
      switchVal = 1.0;
      switchdUdL = 0.0;
      switchd2UdL2 = 0.0;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    return d2EdL2;
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {

    return dEdL;
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] gradient) {
    int i1 = atoms[0].getIndex() - 1;
    int index = i1 * 3;
    gradient[index++] += dEdXdL[0][0];
    gradient[index++] += dEdXdL[0][1];
    gradient[index] += dEdXdL[0][2];
    int i2 = atoms[1].getIndex() - 1;
    index = i2 * 3;
    gradient[index++] += dEdXdL[1][0];
    gradient[index++] += dEdXdL[1][1];
    gradient[index] += dEdXdL[1][2];
  }

  @Override
  public boolean isLambdaScaled() {
    return lambdaTerm;
  }

  /** Log details for this Bond energy term. */
  public void log() {
    logger.info(
        String.format(
            " %s %6d-%s %6d-%s %6.4f  %6.4f  %10.4f",
            "Restraint-Bond",
            atoms[0].getIndex(),
            atoms[0].getAtomType().name,
            atoms[1].getIndex(),
            atoms[1].getAtomType().name,
            bondType.distance,
            value,
            energy));
    if (!(switchingFunction instanceof ConstantSwitch)) {
      logger.info(
          String.format(
              " Switching function (lambda dependence): %s", switchingFunction.toString()));
    }
  }

  @Override
  public String toString() {
    StringBuilder sb =
        new StringBuilder(
            String.format(
                " Distance restraint between atoms %s-%d %s-%d, "
                    + "current distance %10.4g, optimum %10.4g with a %10.4g Angstrom flat bottom, with force constant %10.4g.",
                atoms[0],
                atoms[0].getIndex(),
                atoms[1],
                atoms[1].getIndex(),
                value,
                bondType.distance,
                bondType.flatBottomRadius,
                bondType.forceConstant));
    if (!(switchingFunction instanceof ConstantSwitch)) {
      sb.append(
          String.format(
              "\n Switching function (lambda dependence): %s", switchingFunction.toString()));
    }
    return sb.toString();
  }
}
