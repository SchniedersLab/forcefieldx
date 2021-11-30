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
package ffx.ui.behaviors;

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.util.Iterator;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.WakeupCriterion;
import org.jogamp.java3d.WakeupOnAWTEvent;
import org.jogamp.vecmath.Matrix4d;
import org.jogamp.vecmath.Vector3d;
import org.jogamp.vecmath.Vector3f;

/**
 * The MouseOrbit class implements a mouse orbit behavior.
 *
 * @author Michael J. Schnieders
 */
public class MouseOrbit extends MouseBehavior {

  double xAngle, yAngle;
  double xFactor = 0.01; // .03;
  double yFactor = 0.01; // .03;
  private TransformGroup tg_ghost;
  private MouseBehaviorCallback callback = null;

  /**
   * Constructor for MouseOrbit.
   *
   * @param flags a int.
   * @param VPTG a {@link org.jogamp.java3d.TransformGroup} object.
   */
  public MouseOrbit(int flags, TransformGroup VPTG) {
    super(flags, VPTG);
  }

  /**
   * Return the x-axis movement multipler.
   *
   * @return a double.
   */
  public double getXFactor() {
    return xFactor;
  }

  /**
   * Return the y-axis movement multipler.
   *
   * @return a double.
   */
  public double getYFactor() {
    return yFactor;
  }

  /** initialize */
  public void initialize() {
    super.initialize();
    xAngle = 0;
    yAngle = 0;
    if ((flags & INVERT_INPUT) == INVERT_INPUT) {
      invert = true;
      xFactor *= -1;
      yFactor *= -1;
    }
  }

  /** {@inheritDoc} */
  public void processStimulus(Iterator<WakeupCriterion> criteria) {
    while (criteria.hasNext()) {
      WakeupCriterion wakeup = criteria.next();
      if (wakeup instanceof WakeupOnAWTEvent) {
        AWTEvent[] event = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
        for (AWTEvent awtEvent : event) {
          processMouseEvent((MouseEvent) awtEvent);
          if (((buttonPress) && ((flags & MANUAL_WAKEUP) == 0))
              || ((wakeUp) && ((flags & MANUAL_WAKEUP) != 0))) {
            int id = awtEvent.getID();
            if ((id == MouseEvent.MOUSE_DRAGGED)) {
              x = ((MouseEvent) awtEvent).getX();
              y = ((MouseEvent) awtEvent).getY();
              int dx = x - xLast;
              int dy = y - yLast;
              if (!reset) {
                Transform3D tempT3D = new Transform3D();
                Transform3D orbitT3D = new Transform3D();
                tempT3D.rotX(-dy * yFactor);
                orbitT3D.mul(tempT3D);
                tempT3D.rotY(-dx * xFactor);
                orbitT3D.mul(tempT3D);
                Transform3D tg_ghost_T3D = new Transform3D();
                tg_ghost.getTransform(tg_ghost_T3D);
                Vector3f tg_ghost_vec3f = new Vector3f();
                tg_ghost_T3D.get(tg_ghost_vec3f);
                Matrix4d tg_ghost_mat4d = new Matrix4d();
                tg_ghost_T3D.get(tg_ghost_mat4d);
                Transform3D VPTG_ghost_T3D_inverted = new Transform3D();
                Transform3D VPTG_ghost_T3D_noninverted = new Transform3D();
                // (super.ViewerTG).getTransform(VPTG_ghost_T3D_inverted);
                ViewerTG.getTransform(VPTG_ghost_T3D_inverted);
                // (super.ViewerTG).getTransform(VPTG_ghost_T3D_noninverted);
                ViewerTG.getTransform(VPTG_ghost_T3D_noninverted);
                VPTG_ghost_T3D_inverted.setTranslation(new Vector3d(0.0, 0.0, 0.0));
                VPTG_ghost_T3D_noninverted.setTranslation(new Vector3d(0.0, 0.0, 0.0));
                VPTG_ghost_T3D_inverted.invert();
                tg_ghost_T3D.mul(VPTG_ghost_T3D_inverted, tg_ghost_T3D);
                tg_ghost_T3D.setTranslation(new Vector3d(0.0, 0.0, 0.0));
                if (invert) {
                  tg_ghost_T3D.mul(tg_ghost_T3D, orbitT3D);
                } else {
                  tg_ghost_T3D.mul(orbitT3D, tg_ghost_T3D);
                }
                tg_ghost_T3D.mul(VPTG_ghost_T3D_noninverted, tg_ghost_T3D);
                tg_ghost_T3D.setTranslation(tg_ghost_vec3f);
                tg_ghost.setTransform(tg_ghost_T3D);
                Transform3D VPTG_ghost_T3D = new Transform3D();
                // (super.ViewerTG).getTransform(VPTG_ghost_T3D);
                ViewerTG.getTransform(VPTG_ghost_T3D);
                Vector3f VPTG_ghost_vec3f = new Vector3f();
                VPTG_ghost_T3D.get(VPTG_ghost_vec3f);
                Vector3f temp_vec3f = new Vector3f();
                temp_vec3f.x = VPTG_ghost_vec3f.x - tg_ghost_vec3f.x;
                temp_vec3f.y = VPTG_ghost_vec3f.y - tg_ghost_vec3f.y;
                temp_vec3f.z = VPTG_ghost_vec3f.z - tg_ghost_vec3f.z;
                VPTG_ghost_T3D.setTranslation(temp_vec3f);
                VPTG_ghost_T3D.mul(VPTG_ghost_T3D_inverted, VPTG_ghost_T3D);
                if (invert) {
                  VPTG_ghost_T3D.mul(VPTG_ghost_T3D, orbitT3D);
                } else {
                  VPTG_ghost_T3D.mul(orbitT3D, VPTG_ghost_T3D);
                }
                VPTG_ghost_T3D.mul(VPTG_ghost_T3D_noninverted, VPTG_ghost_T3D);
                VPTG_ghost_T3D.get(temp_vec3f);
                temp_vec3f.x = temp_vec3f.x + tg_ghost_vec3f.x;
                temp_vec3f.y = temp_vec3f.y + tg_ghost_vec3f.y;
                temp_vec3f.z = temp_vec3f.z + tg_ghost_vec3f.z;
                VPTG_ghost_T3D.setTranslation(temp_vec3f);
                // (super.ViewerTG).setTransform(VPTG_ghost_T3D);
                ViewerTG.setTransform(VPTG_ghost_T3D);
                transformChanged(currXform);
                if (callback != null) {
                  callback.transformChanged(MouseBehaviorCallback.ORBIT, currXform);
                }
              } else {
                reset = false;
              }
              xLast = x;
              yLast = y;
            } else if (id == MouseEvent.MOUSE_PRESSED) {
              xLast = ((MouseEvent) awtEvent).getX();
              yLast = ((MouseEvent) awtEvent).getY();
            }
          }
        }
      }
    }
    wakeupOn(mouseCriterion);
  }

  /**
   * Set the x-axis amd y-axis movement multipler with factor.
   *
   * @param factor a double.
   */
  public void setFactor(double factor) {
    xFactor = yFactor = factor;
  }

  /**
   * Set the x-axis amd y-axis movement multipler with xFactor and yFactor respectively.
   *
   * @param xFactor a double.
   * @param yFactor a double.
   */
  public void setFactor(double xFactor, double yFactor) {
    this.xFactor = xFactor;
    this.yFactor = yFactor;
  }

  /**
   * setTransformGroups
   *
   * @param tg a {@link org.jogamp.java3d.TransformGroup} object.
   * @param VPTG a {@link org.jogamp.java3d.TransformGroup} object.
   */
  public void setTransformGroups(TransformGroup tg, TransformGroup VPTG) {
    super.ViewerTG = VPTG;
    tg_ghost = new TransformGroup();
    Transform3D tgT3D = new Transform3D();
    tg.getTransform(tgT3D);
    // Make a ghost TG since no transform on object is to occur
    tg_ghost.setTransform(tgT3D);
  }

  /**
   * The transformChanged method in the callback class will be called every time the transform is
   * updated
   *
   * @param c a {@link ffx.ui.behaviors.MouseBehaviorCallback} object.
   */
  public void setupCallback(MouseBehaviorCallback c) {
    callback = c;
  }

  /**
   * transformChanged
   *
   * @param transform a {@link org.jogamp.java3d.Transform3D} object.
   */
  public void transformChanged(Transform3D transform) {}
}
