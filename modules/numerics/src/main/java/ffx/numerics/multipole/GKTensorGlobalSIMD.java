// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import jdk.incubator.vector.DoubleVector;

/**
 * The GeneralizedKirkwoodTensor class contains utilities for generated Generalized Kirkwood
 * interaction tensors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GKTensorGlobalSIMD extends CoulombTensorGlobalSIMD {

  /**
   * The GK tensor can be constructed for monopoles (GB), dipoles or quadrupoles.
   */
  protected final GKMultipoleOrder multipoleOrder;

  /**
   * The Kirkwood dielectric function for the given multipole order.
   */
  private final double c;

  private final GKSourceSIMD gkSource;

  private final DoubleVector zero = DoubleVector.zero(DoubleVector.SPECIES_PREFERRED);

  /**
   * Construct a new GKTensorGlobal object.
   *
   * @param multipoleOrder The multipole order.
   * @param order          The number of derivatives to complete.
   * @param gkSource       Generate the source terms for the GK recurrence.
   * @param Eh             Homogeneous dielectric constant.
   * @param Es             Solvent dielectric constant.
   */
  public GKTensorGlobalSIMD(GKMultipoleOrder multipoleOrder, int order, GKSourceSIMD gkSource, double Eh, double Es) {
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
   * @return a double.
   */
  @Override
  public DoubleVector multipoleEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    return switch (multipoleOrder) {
      default -> {
        chargeIPotentialAtK(mI, 2);
        DoubleVector eK = multipoleEnergy(mK);
        chargeKPotentialAtI(mK, 2);
        DoubleVector eI = multipoleEnergy(mI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case DIPOLE -> {
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 2);
        DoubleVector eK = multipoleEnergy(mK);
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 2);
        DoubleVector eI = multipoleEnergy(mI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case QUADRUPOLE -> {
        quadrupoleIPotentialAtK(mI, 2);
        DoubleVector eK = multipoleEnergy(mK);
        quadrupoleKPotentialAtI(mK, 2);
        DoubleVector eI = multipoleEnergy(mI);
        yield eK.add(eI).mul(c * 0.5);
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
  public DoubleVector multipoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                 DoubleVector[] Gi, DoubleVector[] Gk,
                                                 DoubleVector[] Ti, DoubleVector[] Tk) {
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
  protected DoubleVector monopoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                   DoubleVector[] Gi, DoubleVector[] Gk,
                                                   DoubleVector[] Ti, DoubleVector[] Tk) {

    // Compute the potential due to a multipole component at site I.
    chargeIPotentialAtK(mI, 3);
    DoubleVector eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Compute the potential due to a multipole component at site K.
    chargeKPotentialAtI(mK, 3);
    DoubleVector eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    double scale = c * 0.5;
    Gi[0] = Gi[0].sub(Gk[0]).mul(scale);
    Gi[1] = Gi[1].sub(Gk[1]).mul(scale);
    Gi[2] = Gi[2].sub(Gk[2]).mul(scale);
    Gk[0] = Gi[0].neg();
    Gk[1] = Gi[1].neg();
    Gk[2] = Gi[2].neg();

    Ti[0] = Ti[0].mul(scale);
    Ti[1] = Ti[1].mul(scale);
    Ti[2] = Ti[2].mul(scale);
    Tk[0] = Tk[0].mul(scale);
    Tk[1] = Tk[1].mul(scale);
    Tk[2] = Tk[2].mul(scale);

    return eK.add(eI).mul(scale);
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
  protected DoubleVector dipoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                 DoubleVector[] Gi, DoubleVector[] Gk,
                                                 DoubleVector[] Ti, DoubleVector[] Tk) {

    // Compute the potential due to a multipole component at site I.
    dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 3);
    DoubleVector eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I pole due to site K multipole.
    // Only torque on the site I dipole.
    multipoleKPotentialAtI(mK, 1);
    dipoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 3);
    DoubleVector eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K pole due to multipole on site I.
    // Only torque on the site K dipole.
    multipoleIPotentialAtK(mI, 1);
    dipoleTorque(mK, Tk);

    double scale = c * 0.5;
    Gi[0] = Gi[0].sub(Gk[0]).mul(scale);
    Gi[1] = Gi[1].sub(Gk[1]).mul(scale);
    Gi[2] = Gi[2].sub(Gk[2]).mul(scale);
    Gk[0] = Gi[0].neg();
    Gk[1] = Gi[1].neg();
    Gk[2] = Gi[2].neg();

    Ti[0] = Ti[0].mul(scale);
    Ti[1] = Ti[1].mul(scale);
    Ti[2] = Ti[2].mul(scale);
    Tk[0] = Tk[0].mul(scale);
    Tk[1] = Tk[1].mul(scale);
    Tk[2] = Tk[2].mul(scale);

    return eK.add(eI).mul(scale);
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
  protected DoubleVector quadrupoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                     DoubleVector[] Gi, DoubleVector[] Gk,
                                                     DoubleVector[] Ti, DoubleVector[] Tk) {

    // Compute the potential due to a multipole component at site I.
    quadrupoleIPotentialAtK(mI, 3);
    DoubleVector eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I quadrupole due to site K multipole.
    multipoleKPotentialAtI(mK, 2);
    quadrupoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    quadrupoleKPotentialAtI(mK, 3);
    DoubleVector eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K quadrupole due to site I multipole.
    multipoleIPotentialAtK(mI, 2);
    quadrupoleTorque(mK, Tk);

    double scale = c * 0.5;
    Gi[0] = Gi[0].sub(Gk[0]).mul(scale);
    Gi[1] = Gi[1].sub(Gk[1]).mul(scale);
    Gi[2] = Gi[2].sub(Gk[2]).mul(scale);
    Gk[0] = Gi[0].neg();
    Gk[1] = Gi[1].neg();
    Gk[2] = Gi[2].neg();

    Ti[0] = Ti[0].mul(scale);
    Ti[1] = Ti[1].mul(scale);
    Ti[2] = Ti[2].mul(scale);
    Tk[0] = Tk[0].mul(scale);
    Tk[1] = Tk[1].mul(scale);
    Tk[2] = Tk[2].mul(scale);

    return eK.add(eI).mul(scale);
  }

  /**
   * GK Permanent multipole Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public DoubleVector multipoleEnergyBornGrad(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    return multipoleEnergy(mI, mK);
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI          PolarizableMultipole at site I.
   * @param mK          PolarizableMultipole at site K.
   * @param scaleEnergy This is ignored, since masking/scaling is not applied to GK interactions
   *                    (everything is intermolecular).
   * @return a double.
   */
  public DoubleVector polarizationEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK, DoubleVector scaleEnergy) {
    return polarizationEnergy(mI, mK);
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public DoubleVector polarizationEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    return switch (multipoleOrder) {
      default -> {
        // Find the GK charge potential of site I at site K.
        chargeIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent charge I.
        DoubleVector eK = polarizationEnergy(mK);
        // Find the GK charge potential of site K at site I.
        chargeKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent charge K.
        DoubleVector eI = polarizationEnergy(mI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case DIPOLE -> {
        // Find the GK dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 1);
        // Energy of induced dipole K in the field of permanent dipole I.
        DoubleVector eK = polarizationEnergy(mK);
        // Find the GK induced dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
        // Energy of permanent multipole K in the field of induced dipole I.
        eK = multipoleEnergy(mK).mul(0.5).add(eK);
        // Find the GK dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 1);
        // Energy of induced dipole I in the field of permanent dipole K.
        DoubleVector eI = polarizationEnergy(mI);
        // Find the GK induced dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
        // Energy of permanent multipole I in the field of induced dipole K.
        eI = multipoleEnergy(mI).mul(0.5).add(eI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case QUADRUPOLE -> {
        // Find the GK quadrupole potential of site I at site K.
        quadrupoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent quadrupole I.
        DoubleVector eK = polarizationEnergy(mK);
        // Find the GK quadrupole potential of site K at site I.
        quadrupoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent quadrupole K.
        DoubleVector eI = polarizationEnergy(mI);
        yield eK.add(eI).mul(c * 0.5);
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
  public DoubleVector polarizationEnergyBorn(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    return switch (multipoleOrder) {
      default -> {
        // Find the GK charge potential of site I at site K.
        chargeIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent charge I.
        DoubleVector eK = polarizationEnergyS(mK);
        // Find the GK charge potential of site K at site I.
        chargeKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent charge K.
        DoubleVector eI = polarizationEnergyS(mI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case DIPOLE -> {
        // Find the GK dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 1);
        // Energy of induced dipole K in the field of permanent dipole I.
        DoubleVector eK = polarizationEnergyS(mK);
        // Find the GK induced dipole potential of site I at site K.
        dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
        // Energy of permanent multipole K in the field of induced dipole I.
        eK = multipoleEnergy(mK).mul(0.5).add(eK);
        // Find the GK dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 1);
        // Energy of induced dipole I in the field of permanent dipole K.
        DoubleVector eI = polarizationEnergyS(mI);
        // Find the GK induced dipole potential of site K at site I.
        dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
        // Energy of permanent multipole I in the field of induced dipole K.
        eI = multipoleEnergy(mI).mul(0.5).add(eI);
        yield eK.add(eI).mul(c * 0.5);
      }
      case QUADRUPOLE -> {
        // Find the GK quadrupole potential of site I at site K.
        quadrupoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent quadrupole I.
        DoubleVector eK = polarizationEnergyS(mK);
        // Find the GK quadrupole potential of site K at site I.
        quadrupoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent quadrupole K.
        DoubleVector eI = polarizationEnergyS(mI);
        yield eK.add(eI).mul(c * 0.5);
      }
    };
  }

  /**
   * Polarization Energy and Gradient.
   *
   * @param mI            PolarizableMultipole at site I.
   * @param mK            PolarizableMultipole at site K.
   * @param inductionMask This is ignored, since masking/scaling is not applied to GK
   *                      interactions (everything is intermolecular).
   * @param energyMask    This is ignored, since masking/scaling is not applied to GK interactions
   *                      (everything is intermolecular).
   * @param mutualMask    This should be set to zero for direction polarization.
   * @param Gi            an array of double values.
   * @param Ti            an array of double values.
   * @param Tk            an array of double values.
   * @return a double.
   */
  @Override
  public DoubleVector polarizationEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                    DoubleVector inductionMask, DoubleVector energyMask, DoubleVector mutualMask,
                                                    DoubleVector[] Gi, DoubleVector[] Ti, DoubleVector[] Tk) {
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
   * @param Gi an array of double values.
   * @return a double.
   */
  public DoubleVector monopolePolarizationEnergyAndGradient(PolarizableMultipoleSIMD mI,
                                                            PolarizableMultipoleSIMD mK, DoubleVector[] Gi) {
    // Find the permanent multipole potential at site k.
    chargeIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    DoubleVector eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = mK.sx.mul(E200).add(mK.sy.mul(E110)).add(mK.sz.mul(E101)).neg();
    Gi[1] = mK.sx.mul(E110).add(mK.sy.mul(E020)).add(mK.sz.mul(E011)).neg();
    Gi[2] = mK.sx.mul(E101).add(mK.sy.mul(E011)).add(mK.sz.mul(E002)).neg();

    // Find the permanent multipole potential and derivatives at site i.
    chargeKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    DoubleVector eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] = Gi[0].add(mI.sx.mul(E200)).add(mI.sy.mul(E110)).add(mI.sz.mul(E101));
    Gi[1] = Gi[1].add(mI.sx.mul(E110)).add(mI.sy.mul(E020)).add(mI.sz.mul(E011));
    Gi[2] = Gi[2].add(mI.sx.mul(E101)).add(mI.sy.mul(E011)).add(mI.sz.mul(E002));

    double scale = c * 0.5;
    Gi[0] = Gi[0].mul(scale);
    Gi[1] = Gi[1].mul(scale);
    Gi[2] = Gi[2].mul(scale);

    // Total polarization energy.
    return eI.add(eK).mul(scale);
  }

  /**
   * Dipole Polarization Energy and Gradient.
   *
   * @param mI         PolarizableMultipole at site I.
   * @param mK         PolarizableMultipole at site K.
   * @param mutualMask This should be set to zero for direction polarization.
   * @param Gi         an array of double values.
   * @param Ti         an array of double values.
   * @param Tk         an array of double values.
   * @return a double.
   */
  public DoubleVector dipolePolarizationEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                          DoubleVector mutualMask, DoubleVector[] Gi,
                                                          DoubleVector[] Ti, DoubleVector[] Tk) {

    // Find the permanent multipole potential and derivatives at site k.
    dipoleIPotentialAtK(mI.dx, mI.dy, mI.dz, 2);
    // Energy of induced dipole k in the field of multipole i.
    DoubleVector eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = mK.sx.mul(E200).add(mK.sy.mul(E110)).add(mK.sz.mul(E101)).neg();
    Gi[1] = mK.sx.mul(E110).add(mK.sy.mul(E020)).add(mK.sz.mul(E011)).neg();
    Gi[2] = mK.sx.mul(E101).add(mK.sy.mul(E011)).add(mK.sz.mul(E002)).neg();

    // Find the potential at K due to the averaged induced dipole at site i.
    dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
    dipoleTorque(mK, Tk);

    // Find the GK induced dipole potential of site I at site K.
    dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 3);
    // Energy of permanent multipole K in the field of induced dipole I.
    eK = eK.add(multipoleEnergy(mK).mul(0.5));

    DoubleVector[] G = new DoubleVector[3];
    multipoleGradient(mK, G);
    Gi[0] = Gi[0].sub(G[0]);
    Gi[1] = Gi[1].sub(G[1]);
    Gi[2] = Gi[2].sub(G[2]);
    multipoleTorque(mK, Tk);

    // Find the permanent multipole potential and derivatives at site i.
    dipoleKPotentialAtI(mK.dx, mK.dy, mK.dz, 2);
    // Energy of induced dipole i in the field of multipole k.
    DoubleVector eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] = Gi[0].add(mI.sx.mul(E200)).add(mI.sy.mul(E110)).add(mI.sz.mul(E101));
    Gi[1] = Gi[1].add(mI.sx.mul(E110)).add(mI.sy.mul(E020)).add(mI.sz.mul(E011));
    Gi[2] = Gi[2].add(mI.sx.mul(E101)).add(mI.sy.mul(E011)).add(mI.sz.mul(E002));

    // Find the potential at I due to the averaged induced dipole at k.
    dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
    dipoleTorque(mI, Ti);

    // Find the GK induced dipole potential of site K at site I.
    dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 3);
    // Energy of permanent multipole I in the field of induced dipole K.
    eI = eI.add(multipoleEnergy(mI).mul(0.5));

    multipoleGradient(mI, G);
    Gi[0] = Gi[0].add(G[0]);
    Gi[1] = Gi[1].add(G[1]);
    Gi[2] = Gi[2].add(G[2]);
    multipoleTorque(mI, Ti);

    // Get the induced-induced portion of the force (Ud . dC/dX . Up).
    // This contribution does not exist for direct polarization (mutualMask == 0.0).
    if (!mutualMask.eq(zero).allTrue()) {
      // Find the potential and its derivatives at k due to induced dipole i.
      dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
      Gi[0] = Gi[0].sub((mK.px.mul(E200).add(mK.py.mul(E110)).add(mK.pz.mul(E101)).mul(mutualMask)));
      Gi[1] = Gi[1].sub((mK.px.mul(E110).add(mK.py.mul(E020)).add(mK.pz.mul(E011)).mul(mutualMask)));
      Gi[2] = Gi[2].sub((mK.px.mul(E101).add(mK.py.mul(E011)).add(mK.pz.mul(E002)).mul(mutualMask)));

      // Find the potential and its derivatives at i due to induced dipole k.
      dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
      Gi[0] = Gi[0].sub((mI.px.mul(E200).add(mI.py.mul(E110)).add(mI.pz.mul(E101)).mul(mutualMask)));
      Gi[1] = Gi[1].sub((mI.px.mul(E110).add(mI.py.mul(E020)).add(mI.pz.mul(E011)).mul(mutualMask)));
      Gi[2] = Gi[2].sub((mI.px.mul(E101).add(mI.py.mul(E011)).add(mI.pz.mul(E002)).mul(mutualMask)));
    }

    // Total polarization energy.
    double scale = c * 0.5;
    DoubleVector energy = eI.add(eK).mul(scale);
    Gi[0] = Gi[0].mul(scale);
    Gi[1] = Gi[1].mul(scale);
    Gi[2] = Gi[2].mul(scale);
    Ti[0] = Ti[0].mul(scale);
    Ti[1] = Ti[1].mul(scale);
    Ti[2] = Ti[2].mul(scale);
    Tk[0] = Tk[0].mul(scale);
    Tk[1] = Tk[1].mul(scale);
    Tk[2] = Tk[2].mul(scale);
    return energy;
  }

  /**
   * Quadrupole Polarization Energy and Gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi an array of double values.
   * @param Ti an array of double values.
   * @param Tk an array of double values.
   * @return a double.
   */
  public DoubleVector quadrupolePolarizationEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                              DoubleVector[] Gi, DoubleVector[] Ti, DoubleVector[] Tk) {

    // Find the permanent multipole potential and derivatives at site k.
    quadrupoleIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    DoubleVector eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = (mK.sx.mul(E200).add(mK.sy.mul(E110)).add(mK.sz.mul(E101))).neg();
    Gi[1] = (mK.sx.mul(E110).add(mK.sy.mul(E020)).add(mK.sz.mul(E011))).neg();
    Gi[2] = (mK.sx.mul(E101).add(mK.sy.mul(E011)).add(mK.sz.mul(E002))).neg();

    // Find the permanent multipole potential and derivatives at site i.
    quadrupoleKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    DoubleVector eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] = Gi[0].add(mI.sx.mul(E200)).add(mI.sy.mul(E110)).add(mI.sz.mul(E101));
    Gi[1] = Gi[1].add(mI.sx.mul(E110)).add(mI.sy.mul(E020)).add(mI.sz.mul(E011));
    Gi[2] = Gi[2].add(mI.sx.mul(E101)).add(mI.sy.mul(E011)).add(mI.sz.mul(E002));

    double scale = c * 0.5;
    Gi[0] = Gi[0].mul(scale);
    Gi[1] = Gi[1].mul(scale);
    Gi[2] = Gi[2].mul(scale);

    // Find the potential and its derivatives at K due to the averaged induced dipole at site i.
    dipoleIPotentialAtK(mI.sx.mul(scale), mI.sy.mul(scale), mI.sz.mul(scale), 2);
    quadrupoleTorque(mK, Tk);

    // Find the potential and its derivatives at I due to the averaged induced dipole at k.
    dipoleKPotentialAtI(mK.sx.mul(scale), mK.sy.mul(scale), mK.sz.mul(scale), 2);
    quadrupoleTorque(mI, Ti);

    // Total polarization energy.
    return eI.add(eK).mul(scale);
  }

  /**
   * GK Direct Polarization Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return Partial derivative of the Polarization energy with respect to a Born grad.
   */
  public DoubleVector polarizationEnergyBornGrad(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    return polarizationEnergyBorn(mI, mK).mul(2.0);
  }

  /**
   * GK Mutual Polarization Contribution to the Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return Mutual Polarization contribution to the partial derivative with respect to a Born grad.
   */
  public DoubleVector mutualPolarizationEnergyBornGrad(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    DoubleVector db = zero;
    if (multipoleOrder == GKMultipoleOrder.DIPOLE) {
      // Find the potential and its derivatives at k due to induced dipole i.
      dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
      db = mK.px.mul(E100).add(mK.py.mul(E010)).add(mK.pz.mul(E001)).mul(0.5);


      // Find the potential and its derivatives at i due to induced dipole k.
      dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
      db = db.add(mI.px.mul(E100).add(mI.py.mul(E010)).add(mI.pz.mul(E001)).mul(0.5));

    }
    return db.mul(c);
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  @Override
  protected void source(DoubleVector[] work) {
    gkSource.source(work, multipoleOrder);
  }

}
