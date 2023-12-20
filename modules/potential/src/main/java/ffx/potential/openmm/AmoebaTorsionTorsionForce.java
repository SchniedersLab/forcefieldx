// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
import ffx.openmm.DoubleArray3D;
import ffx.openmm.Force;
import ffx.openmm.amoeba.TorsionTorsionForce;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.parameters.TorsionTorsionType;

import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static java.lang.String.format;

/**
 * OpenMM TorsionTorsion Force.
 */
public class AmoebaTorsionTorsionForce extends TorsionTorsionForce {

  private static final Logger logger = Logger.getLogger(AmoebaTorsionTorsionForce.class.getName());

  /**
   * Create an OpenMM TorsionTorsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the torsion-torsions.
   */
  public AmoebaTorsionTorsionForce(OpenMMEnergy openMMEnergy) {
    TorsionTorsion[] torsionTorsions = openMMEnergy.getTorsionTorsions();
    if (torsionTorsions == null || torsionTorsions.length < 1) {
      destroy();
      return;
    }

    // Load the torsion-torsions.
    int nTypes = 0;
    LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();

    for (TorsionTorsion torsionTorsion : torsionTorsions) {
      int ia = torsionTorsion.getAtom(0).getXyzIndex() - 1;
      int ib = torsionTorsion.getAtom(1).getXyzIndex() - 1;
      int ic = torsionTorsion.getAtom(2).getXyzIndex() - 1;
      int id = torsionTorsion.getAtom(3).getXyzIndex() - 1;
      int ie = torsionTorsion.getAtom(4).getXyzIndex() - 1;

      TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
      String key = torsionTorsionType.getKey();

      // Check if the TorTor parameters have already been added to the Hash.
      int gridIndex = 0;
      if (torTorTypes.containsKey(key)) {

        // If the TorTor has been added, get its (ordered) index in the Hash.
        int index = 0;
        for (String entry : torTorTypes.keySet()) {
          if (entry.equalsIgnoreCase(key)) {
            gridIndex = index;
            break;
          } else {
            index++;
          }
        }
      } else {
        // Add the new TorTor.
        torTorTypes.put(key, torsionTorsionType);
        gridIndex = nTypes;
        nTypes++;
      }

      Atom atom = torsionTorsion.getChiralAtom();
      int iChiral = -1;
      if (atom != null) {
        iChiral = atom.getXyzIndex() - 1;
      }
      addTorsionTorsion(ia, ib, ic, id, ie, iChiral, gridIndex);
    }

    // Load the Torsion-Torsion parameters.
    DoubleArray values = new DoubleArray(6);
    int gridIndex = 0;
    for (String key : torTorTypes.keySet()) {
      TorsionTorsionType torTorType = torTorTypes.get(key);
      int nx = torTorType.nx;
      int ny = torTorType.ny;
      double[] tx = torTorType.tx;
      double[] ty = torTorType.ty;
      double[] f = torTorType.energy;
      double[] dx = torTorType.dx;
      double[] dy = torTorType.dy;
      double[] dxy = torTorType.dxy;

      // Create the 3D grid.
      DoubleArray3D grid3D = new DoubleArray3D(nx, ny, 6);
      int xIndex = 0;
      int yIndex = 0;
      for (int j = 0; j < nx * ny; j++) {
        int addIndex = 0;
        values.set(addIndex++, tx[xIndex]);
        values.set(addIndex++, ty[yIndex]);
        values.set(addIndex++, OpenMM_KJPerKcal * f[j]);
        values.set(addIndex++, OpenMM_KJPerKcal * dx[j]);
        values.set(addIndex++, OpenMM_KJPerKcal * dy[j]);
        values.set(addIndex, OpenMM_KJPerKcal * dxy[j]);
        grid3D.set(yIndex, xIndex, values);
        xIndex++;
        if (xIndex == nx) {
          xIndex = 0;
          yIndex++;
        }
      }
      setTorsionTorsionGrid(gridIndex++, grid3D.getPointer());
      grid3D.destroy();
    }
    values.destroy();

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("TORSION_TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Torsion-Torsions  \t%6d\t\t%1d", torsionTorsions.length, forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Torsion-Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the torsion-torsions.
   * @return A Torsion-Torsion Force, or null if there are no torsion-torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    TorsionTorsion[] torsionTorsion = openMMEnergy.getTorsionTorsions();
    if (torsionTorsion == null || torsionTorsion.length < 1) {
      return null;
    }
    return new AmoebaTorsionTorsionForce(openMMEnergy);
  }
}
