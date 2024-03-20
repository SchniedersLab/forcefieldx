// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.numerics.multipole;

import ffx.numerics.multipole.GKSource.GK_MULTIPOLE_ORDER;

/**
 * The GeneralizedKirkwoodTensor class contains utilities for generated Generalized Kirkwood
 * interaction tensors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GKTensorQI extends CoulombTensorQI {

  /**
   * The GK tensor can be constructed for monopoles (GB), dipoles or quadrupoles.
   */
  protected final GK_MULTIPOLE_ORDER multipoleOrder;

  /**
   * The Kirkwood dielectric function for the given multipole order.
   */
  private final double c;

  private final GKSource gkSource;

  /**
   * @param multipoleOrder The multipole order.
   * @param order The tensor order.
   * @param gkSource Generate the source terms for the GK recurrence.
   * @param Eh Homogeneous dielectric constant.
   * @param Es Solvent dielectric constant.
   */
  public GKTensorQI(GK_MULTIPOLE_ORDER multipoleOrder, int order, GKSource gkSource, double Eh,
      double Es) {
    super(order);
    this.multipoleOrder = multipoleOrder;
    this.gkSource = gkSource;

    // Load the dielectric function
    c = GKSource.cn(multipoleOrder.getOrder(), Eh, Es);
  }

  /**
   * GK Permanent multipole energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return the GK permanent multipole energy.
   */
  @Override
  public double multipoleEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    return switch (multipoleOrder) {
      default -> {
        chargeIPotentialAtK(mI, 2);
        double eK = multipoleEnergy(mK);
        chargeKPotentialAtI(mK, 2);
        double eI = multipoleEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
      case DIPOLE -> {
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 2);
        double eK = multipoleEnergy(mK);
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 2);
        double eI = multipoleEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
      case QUADRUPOLE -> {
        quadrupoleIPotentialAtK(mI, 2);
        double eK = multipoleEnergy(mK);
        quadrupoleKPotentialAtI(mK, 2);
        double eI = multipoleEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
    };
  }

  /**
   * GK Permanent multipole energy and gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  @Override
  public double multipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {
    return switch (multipoleOrder) {
      default -> monopoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
      case DIPOLE -> dipoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
      case QUADRUPOLE -> quadrupoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
    };
  }

  /**
   * Permanent multipole energy and gradient using the GK monopole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double monopoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    chargeIPotentialAtK(mI, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Compute the potential due to a multipole component at site K.
    chargeKPotentialAtI(mK, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    double scale = c * 0.5;
    Gi[0] = scale * (Gi[0] - Gk[0]);
    Gi[1] = scale * (Gi[1] - Gk[1]);
    Gi[2] = scale * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = scale * Ti[0];
    Ti[1] = scale * Ti[1];
    Ti[2] = scale * Ti[2];
    Tk[0] = scale * Tk[0];
    Tk[1] = scale * Tk[1];
    Tk[2] = scale * Tk[2];

    return scale * (eK + eI);
  }

  /**
   * Permanent multipole energy and gradient using the GK dipole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double dipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I dipole due to site K multipole.
    multipoleKPotentialAtI(mK, 1);
    dipoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K dipole due to site I multipole.
    multipoleIPotentialAtK(mI, 1);
    dipoleTorque(mK, Tk);

    double scale = c * 0.5;
    Gi[0] = scale * (Gi[0] - Gk[0]);
    Gi[1] = scale * (Gi[1] - Gk[1]);
    Gi[2] = scale * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = scale * Ti[0];
    Ti[1] = scale * Ti[1];
    Ti[2] = scale * Ti[2];
    Tk[0] = scale * Tk[0];
    Tk[1] = scale * Tk[1];
    Tk[2] = scale * Tk[2];

    return scale * (eK + eI);
  }

  /**
   * Permanent multipole energy and gradient using the GK quadrupole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double quadrupoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    quadrupoleIPotentialAtK(mI, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I quadrupole due to site K multipole.
    multipoleKPotentialAtI(mK, 2);
    quadrupoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    quadrupoleKPotentialAtI(mK, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K quadrupole due to site I multipole.
    multipoleIPotentialAtK(mI, 2);
    quadrupoleTorque(mK, Tk);

    double scale = c * 0.5;
    Gi[0] = scale * (Gi[0] - Gk[0]);
    Gi[1] = scale * (Gi[1] - Gk[1]);
    Gi[2] = scale * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = scale * Ti[0];
    Ti[1] = scale * Ti[1];
    Ti[2] = scale * Ti[2];
    Tk[0] = scale * Tk[0];
    Tk[1] = scale * Tk[1];
    Tk[2] = scale * Tk[2];

    return scale * (eK + eI);
  }

  /**
   * GK Permanent multipole Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double multipoleEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK) {
    generateTensor();
    return multipoleEnergy(mI, mK);
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param scaleEnergy This is ignored, since masking/scaling is not applied to GK
   *     interactions.
   * @return a double.
   */
  @Override
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK,
      double scaleEnergy) {
    return polarizationEnergy(mI, mK);
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    return switch (multipoleOrder) {
      default -> {
        // Find the GK charge potential of site I at site K.
        chargeIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent charge I.
        double eK = polarizationEnergy(mK);
        // Find the GK charge potential of site K at site I.
        chargeKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent charge K.
        double eI = polarizationEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
      case DIPOLE -> {
        // Find the GK dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 1);
        // Energy of induced dipole K in the field of permanent dipole I.
        double eK = polarizationEnergy(mK);
        // Find the GK induced dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
        // Energy of permanent multipole K in the field of induced dipole I.
        eK += 0.5 * multipoleEnergy(mK);
        // Find the GK dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 1);
        // Energy of induced dipole I in the field of permanent dipole K.
        double eI = polarizationEnergy(mI);
        // Find the GK induced dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
        // Energy of permanent multipole I in the field of induced dipole K.
        eI += 0.5 * multipoleEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
      case QUADRUPOLE -> {
        // Find the GK quadrupole potential of site I at site K.
        quadrupoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent quadrupole I.
        double eK = polarizationEnergy(mK);
        // Find the GK quadrupole potential of site K at site I.
        quadrupoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent quadrupole K.
        double eI = polarizationEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
    };
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double polarizationEnergyBorn(PolarizableMultipole mI, PolarizableMultipole mK) {
    return switch (multipoleOrder) {
      default -> {
        // Find the GK charge potential of site I at site K.
        chargeIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent charge I.
        double eK = polarizationEnergyS(mK);
        // Find the GK charge potential of site K at site I.
        chargeKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent charge K.
        double eI = polarizationEnergyS(mI);
        yield c * 0.5 * (eK + eI);
      }
      case DIPOLE -> {
        // Find the GK dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 1);
        // Energy of induced dipole K in the field of permanent dipole I.
        double eK = polarizationEnergyS(mK);
        // Find the GK induced dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
        // Energy of permanent multipole K in the field of induced dipole I.
        eK += 0.5 * multipoleEnergy(mK);
        // Find the GK dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 1);
        // Energy of induced dipole I in the field of permanent dipole K.
        double eI = polarizationEnergyS(mI);
        // Find the GK induced dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
        // Energy of permanent multipole I in the field of induced dipole K.
        eI += 0.5 * multipoleEnergy(mI);
        yield c * 0.5 * (eK + eI);
      }
      case QUADRUPOLE -> {
        // Find the GK quadrupole potential of site I at site K.
        quadrupoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent quadrupole I.
        double eK = polarizationEnergyS(mK);
        // Find the GK quadrupole potential of site K at site I.
        quadrupoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent quadrupole K.
        double eI = polarizationEnergyS(mI);
        yield c * 0.5 * (eK + eI);
      }
    };
  }

  /**
   * Polarization Energy and Gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param inductionMask This is ignored, since masking/scaling is not applied to GK
   *     interactions (everything is intermolecular).
   * @param energyMask This is ignored, since masking/scaling is not applied to GK interactions
   *     (everything is intermolecular).
   * @param mutualMask This should be set to zero for direction polarization.
   * @param Gi an array of {@link double} objects.
   * @param Ti an array of {@link double} objects.
   * @param Tk an array of {@link double} objects.
   * @return a double.
   */
  @Override
  public double polarizationEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double inductionMask, double energyMask, double mutualMask, double[] Gi, double[] Ti,
      double[] Tk) {
    return switch (multipoleOrder) {
      default -> monopolePolarizationEnergyAndGradient(mI, mK, Gi);
      case DIPOLE -> dipolePolarizationEnergyAndGradient(mI, mK, mutualMask, Gi, Ti, Tk);
      case QUADRUPOLE -> quadrupolePolarizationEnergyAndGradient(mI, mK, Gi, Ti, Tk);
    };
  }

  /**
   * Monopole Polarization Energy and Gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi an array of {@link double} objects.
   * @return a double.
   */
  public double monopolePolarizationEnergyAndGradient(PolarizableMultipole mI,
      PolarizableMultipole mK, double[] Gi) {

    // Find the permanent multipole potential at site k.
    chargeIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    double eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = -(mK.sx * E200 + mK.sy * E110 + mK.sz * E101);
    Gi[1] = -(mK.sx * E110 + mK.sy * E020 + mK.sz * E011);
    Gi[2] = -(mK.sx * E101 + mK.sy * E011 + mK.sz * E002);

    // Find the permanent multipole potential and derivatives at site i.
    chargeKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    double eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] += (mI.sx * E200 + mI.sy * E110 + mI.sz * E101);
    Gi[1] += (mI.sx * E110 + mI.sy * E020 + mI.sz * E011);
    Gi[2] += (mI.sx * E101 + mI.sy * E011 + mI.sz * E002);

    double scale = c * 0.5;
    Gi[0] *= scale;
    Gi[1] *= scale;
    Gi[2] *= scale;

    // Total polarization energy.
    return scale * (eI + eK);
  }

  /**
   * Dipole Polarization Energy and Gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param mutualMask This should be set to zero for direction polarization.
   * @param Gi an array of {@link double} objects.
   * @param Ti an array of {@link double} objects.
   * @param Tk an array of {@link double} objects.
   * @return a double.
   */
  public double dipolePolarizationEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double mutualMask, double[] Gi, double[] Ti, double[] Tk) {

    // Find the permanent dipole potential and derivatives at site K.
    dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 2);
    // Energy of induced dipole k in the field of dipole I.
    double eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom K.
    Gi[0] = -(mK.sx * E200 + mK.sy * E110 + mK.sz * E101);
    Gi[1] = -(mK.sx * E110 + mK.sy * E020 + mK.sz * E011);
    Gi[2] = -(mK.sx * E101 + mK.sy * E011 + mK.sz * E002);
    // Find the potential at K due to the averaged induced dipole at site I.
    dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
    dipoleTorque(mI, Ti);

    // Find the GK induced dipole potential of site I at site K.
    dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
    // Energy of permanent multipole K in the field of induced dipole I.
    eK += 0.5 * multipoleEnergy(mK);
    // Find the GK induced dipole potential of site I at site K.
    dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 3);
    double[] G = new double[3];
    multipoleGradient(mK, G);
    Gi[0] -= G[0];
    Gi[1] -= G[1];
    Gi[2] -= G[2];
    multipoleTorque(mK, Tk);

    // Find the permanent multipole potential and derivatives at site i.
    dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 2);
    // Energy of induced dipole i in the field of multipole k.
    double eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] += (mI.sx * E200 + mI.sy * E110 + mI.sz * E101);
    Gi[1] += (mI.sx * E110 + mI.sy * E020 + mI.sz * E011);
    Gi[2] += (mI.sx * E101 + mI.sy * E011 + mI.sz * E002);
    // Find the potential at I due to the averaged induced dipole at k.
    dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
    dipoleTorque(mK, Tk);

    // Find the GK induced dipole potential of site K at site I.
    dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
    // Energy of permanent multipole I in the field of induced dipole K.
    eI += 0.5 * multipoleEnergy(mI);
    // Find the GK induced dipole potential of site K at site I.
    dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 3);
    G = new double[3];
    multipoleGradient(mI, G);
    Gi[0] += G[0];
    Gi[1] += G[1];
    Gi[2] += G[2];
    multipoleTorque(mI, Ti);

    // Get the induced-induced portion of the force (Ud . dC/dX . Up).
    // This contribution does not exist for direct polarization (mutualMask == 0.0).
    if (mutualMask != 0.0) {
      // Find the potential and its derivatives at k due to induced dipole i.
      dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
      Gi[0] -= mutualMask * (mK.px * E200 + mK.py * E110 + mK.pz * E101);
      Gi[1] -= mutualMask * (mK.px * E110 + mK.py * E020 + mK.pz * E011);
      Gi[2] -= mutualMask * (mK.px * E101 + mK.py * E011 + mK.pz * E002);

      // Find the potential and its derivatives at i due to induced dipole k.
      dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
      Gi[0] += mutualMask * (mI.px * E200 + mI.py * E110 + mI.pz * E101);
      Gi[1] += mutualMask * (mI.px * E110 + mI.py * E020 + mI.pz * E011);
      Gi[2] += mutualMask * (mI.px * E101 + mI.py * E011 + mI.pz * E002);
    }

    // Total polarization energy.
    double scale = c * 0.5;
    double energy = scale * (eI + eK);
    Gi[0] *= scale;
    Gi[1] *= scale;
    Gi[2] *= scale;
    Ti[0] *= scale;
    Ti[1] *= scale;
    Ti[2] *= scale;
    Tk[0] *= scale;
    Tk[1] *= scale;
    Tk[2] *= scale;

    return energy;
  }

  /**
   * Quadrupole Polarization Energy and Gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi an array of {@link double} objects.
   * @param Ti an array of {@link double} objects.
   * @param Tk an array of {@link double} objects.
   * @return a double.
   */
  public double quadrupolePolarizationEnergyAndGradient(PolarizableMultipole mI,
      PolarizableMultipole mK, double[] Gi, double[] Ti, double[] Tk) {

    // Find the permanent multipole potential and derivatives at site k.
    quadrupoleIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    double eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = -(mK.sx * E200 + mK.sy * E110 + mK.sz * E101);
    Gi[1] = -(mK.sx * E110 + mK.sy * E020 + mK.sz * E011);
    Gi[2] = -(mK.sx * E101 + mK.sy * E011 + mK.sz * E002);

    // Find the permanent multipole potential and derivatives at site i.
    quadrupoleKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    double eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] += (mI.sx * E200 + mI.sy * E110 + mI.sz * E101);
    Gi[1] += (mI.sx * E110 + mI.sy * E020 + mI.sz * E011);
    Gi[2] += (mI.sx * E101 + mI.sy * E011 + mI.sz * E002);

    double scale = c * 0.5;
    Gi[0] *= scale;
    Gi[1] *= scale;
    Gi[2] *= scale;

    // Find the potential and its derivatives at K due to the averaged induced dipole at site i.
    dipoleIPotentialAtK(scale * mI.sx, scale * mI.sy, scale * mI.sz, 2);
    quadrupoleTorque(mK, Tk);

    // Find the potential and its derivatives at I due to the averaged induced dipole at k.
    dipoleKPotentialAtI(scale * mK.sx, scale * mK.sy, scale * mK.sz, 2);
    quadrupoleTorque(mI, Ti);

    // Total polarization energy.
    return scale * (eI + eK);
  }

  /**
   * GK Direct Polarization Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return Partial derivative of the Polarization energy with respect to a Born grad.
   */
  public double polarizationEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK) {
    generateTensor();
    return 2.0 * polarizationEnergyBorn(mI, mK);
  }

  /**
   * GK Mutual Polarization Contribution to the Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return Mutual Polarization contribution to the partial derivative with respect to a Born grad.
   */
  public double mutualPolarizationEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK) {
    double db = 0.0;
    if (multipoleOrder == GK_MULTIPOLE_ORDER.DIPOLE) {
      // Find the potential and its derivatives at k due to induced dipole i.
      dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
      db = 0.5 * (mK.px * E100 + mK.py * E010 + mK.pz * E001);

      // Find the potential and its derivatives at i due to induced dipole k.
      dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
      db += 0.5 * (mI.px * E100 + mI.py * E010 + mI.pz * E001);
    }
    return c * db;
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  @Override
  protected void source(double[] work) {
    gkSource.source(work, multipoleOrder);
  }
}
