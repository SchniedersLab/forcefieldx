// ******************************************************************************
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
// ******************************************************************************
package ffx.ui;

import ffx.potential.bonded.MSNode;
import java.util.ArrayList;
import java.util.Enumeration;

/**
 * The Selection class will be used to make recursive multiscale selections, however its
 * implementation is not yet complete.
 *
 * @author Michael J. Schnieders
 */
public class Selection {

  MSNode m;
  private Class scale;
  private String criteria;
  private ArrayList<MSNode> selected;

  private Selection(MSNode m, Class scale, String criteria) {
    this.scale = scale;
    this.criteria = criteria;
    this.m = m;
    selected = new ArrayList<>();
  }

  private Selection(MSNode m, String scale, String criteria) {
    try {
      this.scale = Class.forName(scale);
    } catch (Exception e) {
      this.scale = null;
    }
    this.criteria = criteria;
    this.m = m;
    selected = new ArrayList<>();
  }

  /**
   * select
   *
   * @param m a {@link ffx.potential.bonded.MSNode} object.
   * @param scale a {@link java.lang.Class} object.
   * @param criteria a {@link java.lang.String} object.
   * @return a {@link ffx.ui.Selection} object.
   */
  public static Selection select(MSNode m, Class scale, String criteria) {
    Selection s = new Selection(m, scale, criteria);
    s.evaluate();
    return s;
  }

  /**
   * select
   *
   * @param m a {@link ffx.potential.bonded.MSNode} object.
   * @param scale a {@link java.lang.String} object.
   * @param criteria a {@link java.lang.String} object.
   * @return a {@link ffx.ui.Selection} object.
   */
  public static Selection select(MSNode m, String scale, String criteria) {
    Selection s = new Selection(m, scale, criteria);
    s.evaluate();
    return s;
  }

  /**
   * and
   *
   * @param scale a {@link java.lang.Class} object.
   * @param criteria a {@link java.lang.String} object.
   * @return a {@link ffx.ui.Selection} object.
   */
  public Selection and(Class scale, String criteria) {
    this.scale = scale;
    this.criteria = criteria;
    evaluate();
    return this;
  }

  /**
   * or
   *
   * @param scale a {@link java.lang.Class} object.
   * @param criteria a {@link java.lang.String} object.
   * @return a {@link ffx.ui.Selection} object.
   */
  public Selection or(Class scale, String criteria) {
    this.scale = scale;
    this.criteria = criteria;
    evaluate();
    return this;
  }

  private void evaluate() {
    Enumeration e = m.depthFirstEnumeration();
    while (e.hasMoreElements()) {
      MSNode n = (MSNode) e.nextElement();
      if (scale.isInstance(n)) {
        selected.add(n);
      }
    }
  }
}
