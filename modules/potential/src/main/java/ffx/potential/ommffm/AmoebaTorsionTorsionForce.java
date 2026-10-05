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

import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.amoeba.DoubleArray3D;
import ffx.openmm.ffm.amoeba.TorsionTorsionForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.terms.TorsionTorsionPotentialEnergy;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static java.lang.String.format;

/**
 * OpenMM Torsion-Torsion Force backed by FFM {@link TorsionTorsionForce}.
 */
public class AmoebaTorsionTorsionForce extends TorsionTorsionForce {

  private static final Logger logger = Logger.getLogger(AmoebaTorsionTorsionForce.class.getName());

  /**
   * Create an OpenMM TorsionTorsion Force.
   *
   * @param torsionTorsionPotentialEnergy The TorsionTorsionPotentialEnergy instance.
   */
  public AmoebaTorsionTorsionForce(TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy) {
    TorsionTorsion[] torsionTorsions = torsionTorsionPotentialEnergy.getTorsionTorsionArray();

    // Load the torsion-torsions.
    Map<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();
    for (TorsionTorsion torsionTorsion : torsionTorsions) {
      int ia = torsionTorsion.getAtom(0).getArrayIndex();
      int ib = torsionTorsion.getAtom(1).getArrayIndex();
      int ic = torsionTorsion.getAtom(2).getArrayIndex();
      int id = torsionTorsion.getAtom(3).getArrayIndex();
      int ie = torsionTorsion.getAtom(4).getArrayIndex();

      TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
      String key = torsionTorsionType.getKey();

      // Check if the TorTor parameters have already been added to the map.
      int gridIndex;
      if (torTorTypes.containsKey(key)) {
        int index = 0;
        gridIndex = 0;
        for (String entry : torTorTypes.keySet()) {
          if (entry.equalsIgnoreCase(key)) {
            gridIndex = index;
            break;
          } else {
            index++;
          }
        }
      } else {
        torTorTypes.put(key, torsionTorsionType);
        gridIndex = torTorTypes.size() - 1;
      }

      Atom atom = torsionTorsion.getChiralAtom();
      int iChiral = -1;
      if (atom != null) {
        iChiral = atom.getArrayIndex();
      }
      addTorsionTorsion(ia, ib, ic, id, ie, iChiral, gridIndex);
    }

    // Load the Torsion-Torsion parameters.
    try (DoubleArray values = new DoubleArray(6)) {
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

        try (DoubleArray3D grid3D = new DoubleArray3D(nx, ny, 6)) {
          int xIndex = 0;
          int yIndex = 0;
          for (int j = 0; j < nx * ny; j++) {
            int addIndex = 0;
            values.set(addIndex++, tx[xIndex]);
            values.set(addIndex++, ty[yIndex]);
            values.set(addIndex++, KJ_PER_KCAL * f[j]);
            values.set(addIndex++, KJ_PER_KCAL * dx[j]);
            values.set(addIndex++, KJ_PER_KCAL * dy[j]);
            values.set(addIndex, KJ_PER_KCAL * dxy[j]);
            grid3D.set(yIndex, xIndex, values);
            xIndex++;
            if (xIndex == nx) {
              xIndex = 0;
              yIndex++;
            }
          }
          setTorsionTorsionGrid(gridIndex++, grid3D);
        }
      }
    }

    int forceGroup = torsionTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Torsion-Torsions:                  %10d", torsionTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a Dual Topology OpenMM TorsionTorsion Force.
   *
   * @param torsionTorsionPotentialEnergy The TorsionTorsionPotentialEnergy instance.
   * @param topology                      The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy      The OpenMMDualTopologyEnergy instance.
   */
  public AmoebaTorsionTorsionForce(TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy,
                                   int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    TorsionTorsion[] torsionTorsions = torsionTorsionPotentialEnergy.getTorsionTorsionArray();

    Map<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();
    for (TorsionTorsion torsionTorsion : torsionTorsions) {
      int ia = torsionTorsion.getAtom(0).getArrayIndex();
      int ib = torsionTorsion.getAtom(1).getArrayIndex();
      int ic = torsionTorsion.getAtom(2).getArrayIndex();
      int id = torsionTorsion.getAtom(3).getArrayIndex();
      int ie = torsionTorsion.getAtom(4).getArrayIndex();

      TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
      String key = torsionTorsionType.getKey();

      // Check if the TorTor parameters have already been added to the map.
      int gridIndex;
      if (torTorTypes.containsKey(key)) {
        int index = 0;
        gridIndex = 0;
        for (String entry : torTorTypes.keySet()) {
          if (entry.equalsIgnoreCase(key)) {
            gridIndex = index;
            break;
          } else {
            index++;
          }
        }
      } else {
        torTorTypes.put(key, torsionTorsionType);
        gridIndex = torTorTypes.size() - 1;
      }

      Atom atom = torsionTorsion.getChiralAtom();
      int iChiral = -1;
      if (atom != null) {
        iChiral = atom.getArrayIndex();
        iChiral = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, iChiral);
      }
      ia = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ia);
      ib = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ib);
      ic = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ic);
      id = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, id);
      ie = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ie);
      addTorsionTorsion(ia, ib, ic, id, ie, iChiral, gridIndex);
    }

    // Load the Torsion-Torsion parameters.
    try (DoubleArray values = new DoubleArray(6)) {
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

        try (DoubleArray3D grid3D = new DoubleArray3D(nx, ny, 6)) {
          int xIndex = 0;
          int yIndex = 0;
          for (int j = 0; j < nx * ny; j++) {
            int addIndex = 0;
            values.set(addIndex++, tx[xIndex]);
            values.set(addIndex++, ty[yIndex]);
            values.set(addIndex++, KJ_PER_KCAL * f[j]);
            values.set(addIndex++, KJ_PER_KCAL * dx[j]);
            values.set(addIndex++, KJ_PER_KCAL * dy[j]);
            values.set(addIndex, KJ_PER_KCAL * dxy[j]);
            grid3D.set(yIndex, xIndex, values);
            xIndex++;
            if (xIndex == nx) {
              xIndex = 0;
              yIndex++;
            }
          }
          setTorsionTorsionGrid(gridIndex++, grid3D);
        }
      }
    }

    int forceGroup = torsionTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Torsion-Torsions:                  %10d", torsionTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Torsion-Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the torsion-torsions.
   * @return A Torsion-Torsion Force, or null if there are no torsion-torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy =
        openMMEnergy.getTorsionTorsionPotentialEnergy();
    if (torsionTorsionPotentialEnergy == null) {
      return null;
    }
    return new AmoebaTorsionTorsionForce(torsionTorsionPotentialEnergy);
  }

  /**
   * Convenience method to construct a Dual Topology OpenMM Torsion-Torsion Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return A Torsion-Torsion Force, or null if there are no torsion-torsions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy = forceFieldEnergy.getTorsionTorsionPotentialEnergy();
    if (torsionTorsionPotentialEnergy == null) {
      return null;
    }
    return new AmoebaTorsionTorsionForce(torsionTorsionPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Convenience method to construct a Dual Topology OpenMM Torsion-Torsion Force.
   *
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return A Torsion-Torsion Force, or null if there are no torsion-torsions.
   */
  public static Force constructForce(OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(0);
    TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy = forceFieldEnergy.getTorsionTorsionPotentialEnergy();
    if (torsionTorsionPotentialEnergy == null) {
      return null;
    }
    return new AmoebaTorsionTorsionForce(torsionTorsionPotentialEnergy, 0, openMMDualTopologyEnergy);
  }
}
