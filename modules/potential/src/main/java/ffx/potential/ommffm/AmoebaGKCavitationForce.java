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
package ffx.potential.ommffm;

import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.amoeba.GKCavitationForce;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.implicit.DispersionRegion;

import java.util.logging.Logger;

/**
 * AMOEBA GK Cavitation Force compatibility wrapper.
 */
public class AmoebaGKCavitationForce extends GKCavitationForce {

  private static final Logger logger = Logger.getLogger(AmoebaGKCavitationForce.class.getName());

  /**
   * Constructor.
   *
   * @param openMMEnergy OpenMM energy.
   */
  public AmoebaGKCavitationForce(OpenMMEnergy openMMEnergy) {
    logger.severe(" The AmoebaGKCavitationForce is not currently supported.");
  }

  /**
   * Convenience method to construct an AMOEBA Cavitation Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the cavitation information.
   * @return An AMOEBA Cavitation Force, or null if not supported.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      return null;
    }
    DispersionRegion dispersionRegion = gk.getDispersionRegion();
    if (dispersionRegion == null) {
      return null;
    }
    return null;
  }

  /**
   * Update the Cavitation force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy The OpenMM energy term.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    // Unsupported.
  }

}
