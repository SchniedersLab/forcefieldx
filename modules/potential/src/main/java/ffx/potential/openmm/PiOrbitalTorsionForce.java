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
package ffx.potential.openmm;

import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.openmm.CustomCompoundBondForce;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.parameters.PiOrbitalTorsionType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static java.lang.String.format;

/**
 * OpenMM Pi-Orbital Torsion Force.
 */
public class PiOrbitalTorsionForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(PiOrbitalTorsionForce.class.getName());

  /**
   * Create an Pi-Orbital Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the out-of-plane bends.
   */
  public PiOrbitalTorsionForce(OpenMMEnergy openMMEnergy) {
    super(6, openMMEnergy.getPiOrbitalTorsionEnergyString());
    PiOrbitalTorsion[] piOrbitalTorsions = openMMEnergy.getPiOrbitalTorsions();
    if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
      return;
    }

    addPerBondParameter("k");
    setName("PiOrbitalTorsion");

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
      int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
      int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
      int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
      PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
      double k = OpenMM_KJPerKcal * type.forceConstant * piOrbitalTorsion.piOrbitalTorsionType.piTorsUnit;
      particles.append(a1);
      particles.append(a2);
      particles.append(a3);
      particles.append(a4);
      particles.append(a5);
      particles.append(a6);
      parameters.append(k);
      addBond(particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("PI_ORBITAL_TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Pi-Orbital Torsions  \t%6d\t\t%1d", piOrbitalTorsions.length, forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Pi-Orbital Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the pi-orbital torsions.
   * @return An OpenMM Pi-Orbital Torsion Force, or null if there are no pi-orbital torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    PiOrbitalTorsion[] piOrbitalTorsions = openMMEnergy.getPiOrbitalTorsions();
    if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
      return null;
    }
    return new PiOrbitalTorsionForce(openMMEnergy);
  }

  /**
   * Update the Pi-Orbital Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the pi-orbital torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    PiOrbitalTorsion[] piOrbitalTorsions = openMMEnergy.getPiOrbitalTorsions();
    if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
      return;
    }

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
      int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
      int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
      int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
      PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
      double k = OpenMM_KJPerKcal * type.forceConstant * piOrbitalTorsion.piOrbitalTorsionType.piTorsUnit;
      particles.append(a1);
      particles.append(a2);
      particles.append(a3);
      particles.append(a4);
      particles.append(a5);
      particles.append(a6);
      parameters.append(k);
      setBondParameters(index++, particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }
}
