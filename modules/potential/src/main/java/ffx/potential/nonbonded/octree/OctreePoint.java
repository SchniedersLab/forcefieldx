//******************************************************************************
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
//******************************************************************************
package ffx.potential.nonbonded.octree;

import java.util.logging.Logger;

/**
 * OctreePoint: Object class for Octree method presented in the Fast Multipole Method (FMM) tutorial
 * from the Barba Group: https://github.com/barbagroup/FMM_tutorial
 */
public class OctreePoint {

  private static final Logger logger = Logger.getLogger(OctreePoint.class.getName());

  /** Coordinates of the point */
  private double x;

  private double y;
  private double z;

  /** Maximum value for random point coordinate Default = 1.0 */
  private double domain = 1.0;

  public OctreePoint(double coords[], double domain) {

    // Set max for random point coordinate generation
    setDomain(domain);

    // Assign input coordinates, if they are given
    // Otherwise, assign random coordinates between 0 and domain
    if (coords.length > 0) {
      if (coords.length == 3) {
        this.x = coords[0];
        this.y = coords[1];
        this.z = coords[2];
      } else {
        logger.warning("Coordinate array must have three points");
      }
    } else {
      this.x = this.domain * Math.random();
      this.y = this.domain * Math.random();
      this.z = this.domain * Math.random();
    }
  }

  public double distance(OctreePoint other) {
    return Math.sqrt(
        Math.pow((this.x - other.x), 2)
            + Math.pow((this.y - other.y), 2)
            + Math.pow((this.z - other.z), 2));
  }

  public double distance(OctreeCell other) {
    return Math.sqrt(
        Math.pow((this.x - other.getX()), 2)
            + Math.pow((this.y - other.getY()), 2)
            + Math.pow((this.z - other.getZ()), 2));
  }

  public double getX() {
    return this.x;
  }

  public double getY() {
    return this.y;
  }

  public double getZ() {
    return this.z;
  }

  public void setDomain(double domain) {
    this.domain = domain;
  }
}
