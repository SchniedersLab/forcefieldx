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
package ffx.openmm;

import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_getLocalPosition;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_getOriginWeights;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_getXWeights;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalCoordinatesSite_getYWeights;

/**
 * This is a VirtualSite that computes the particle location based on a local coordinate system.
 * The site is placed at a fixed position relative to a local coordinate system defined by
 * other particles. This is useful for representing sites like lone pairs or pi-electron
 * centers that have well-defined positions relative to the molecular framework.
 * <p>
 * The local coordinate system is defined by an origin and two direction vectors (x and y axes).
 * Each of these is specified as a weighted average of particle positions. The z axis is taken
 * to be the cross product of the x and y axes.
 */
public class LocalCoordinatesSite extends VirtualSite {

  /**
   * Create a LocalCoordinatesSite using weighted averages to define the coordinate system.
   *
   * @param originWeights The weights for computing the origin as a weighted average of particle positions.
   * @param xWeights      The weights for computing the x direction as a weighted average of particle positions.
   * @param yWeights      The weights for computing the y direction as a weighted average of particle positions.
   * @param zWeights      The weights for computing the z direction as a weighted average of particle positions.
   * @param localPosition The position of the virtual site in the local coordinate system.
   */
  public LocalCoordinatesSite(PointerByReference originWeights, PointerByReference xWeights,
                              PointerByReference yWeights, PointerByReference zWeights, OpenMM_Vec3 localPosition) {
    super(OpenMM_LocalCoordinatesSite_create(originWeights, xWeights, yWeights, zWeights, localPosition));
  }

  /**
   * Create a LocalCoordinatesSite using specific particles and vectors to define the coordinate system.
   *
   * @param particle1     The index of the first particle.
   * @param particle2     The index of the second particle.
   * @param particle3     The index of the third particle.
   * @param originWeights The weights for the origin calculation.
   * @param xdir          The x direction vector.
   * @param ydir          The y direction vector.
   * @param localPosition The position of the virtual site in the local coordinate system.
   */
  public LocalCoordinatesSite(int particle1, int particle2, int particle3, OpenMM_Vec3 originWeights,
                              OpenMM_Vec3 xdir, OpenMM_Vec3 ydir, OpenMM_Vec3 localPosition) {
    super(OpenMM_LocalCoordinatesSite_create_2(particle1, particle2, particle3, originWeights, xdir, ydir, localPosition));
  }

  /**
   * Destroy the virtual site.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_LocalCoordinatesSite_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the position of the virtual site in the local coordinate system.
   *
   * @return The local position of the virtual site.
   */
  public OpenMM_Vec3 getLocalPosition() {
    return OpenMM_LocalCoordinatesSite_getLocalPosition(pointer);
  }

  /**
   * Get the weights used for computing the origin of the local coordinate system.
   *
   * @param weights The weights for the origin calculation (output).
   */
  public void getOriginWeights(PointerByReference weights) {
    OpenMM_LocalCoordinatesSite_getOriginWeights(pointer, weights);
  }

  /**
   * Get the weights used for computing the x direction of the local coordinate system.
   *
   * @param weights The weights for the x direction calculation (output).
   */
  public void getXWeights(PointerByReference weights) {
    OpenMM_LocalCoordinatesSite_getXWeights(pointer, weights);
  }

  /**
   * Get the weights used for computing the y direction of the local coordinate system.
   *
   * @param weights The weights for the y direction calculation (output).
   */
  public void getYWeights(PointerByReference weights) {
    OpenMM_LocalCoordinatesSite_getYWeights(pointer, weights);
  }
}