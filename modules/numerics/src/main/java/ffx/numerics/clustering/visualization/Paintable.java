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
package ffx.numerics.clustering.visualization;

import java.awt.Graphics2D;

/**
 * Implemented by visual components of the dendrogram.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Paintable {

  /**
   * Paints this visual component onto the provided Graphics2D context.
   *
   * @param g              the Graphics2D to draw into
   * @param xDisplayOffset x-axis pixel offset applied to model coordinates
   * @param yDisplayOffset y-axis pixel offset applied to model coordinates
   * @param xDisplayFactor scale factor to convert model X to pixels
   * @param yDisplayFactor scale factor to convert model Y to pixels
   * @param decorated      if true, draw additional decorations (e.g., distance labels)
   */
  void paint(Graphics2D g, int xDisplayOffset, int yDisplayOffset, double xDisplayFactor, double yDisplayFactor, boolean decorated);

}
