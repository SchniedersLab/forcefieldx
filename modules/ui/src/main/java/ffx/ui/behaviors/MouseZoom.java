/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui.behaviors;

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.util.Enumeration;

import javax.media.j3d.Behavior;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.WakeupCriterion;
import javax.media.j3d.WakeupOnAWTEvent;
import javax.vecmath.Vector3d;

/**
 * The MouseZoom class implements a Mouse Zoom behavior.
 */
public class MouseZoom extends MouseBehavior {
	double z_factor = 0.0002;
	Vector3d translation = new Vector3d();
	private MouseBehaviorCallback callback = null;
	int mouseButton = MouseEvent.BUTTON2_DOWN_MASK;
	int doneID = 0;
	boolean first = true;

	public MouseZoom(int flags, TransformGroup VPTG) {
		super(flags, VPTG);
	}

	public MouseZoom(int flags, TransformGroup VPTG, Behavior behavior,
			int postID, int dID) {
		super(flags, VPTG, behavior, postID);
		doneID = dID;
	}

	/*
	 * Return the y-axis movement multipler.
	 */
	public double getFactor() {
		return z_factor;
	}

	public void initialize() {
		super.initialize();
		if ((flags & INVERT_INPUT) == INVERT_INPUT) {
			z_factor *= -1;
			invert = true;
		}
	}

	public void setMouseButton(int button) {
		mouseButton = button;
	}

	public void processStimulus(Enumeration criteria) {
		AWTEvent[] event;
		boolean done = false;
		while (criteria.hasMoreElements()) {
			WakeupCriterion wakeup = (WakeupCriterion) criteria.nextElement();
			if (wakeup instanceof WakeupOnAWTEvent) {
				event = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
				for (int i = 0; i < event.length; i++) {
					processMouseEvent((MouseEvent) event[i]);
					int id = event[i].getID();
					MouseEvent mevent = (MouseEvent) event[i];
					int mod = mevent.getModifiersEx();
					boolean middleButton = ((mod & mouseButton) == mouseButton);
					if (!middleButton) {
						middleButton = ((mod & MouseEvent.ALT_DOWN_MASK) == MouseEvent.ALT_DOWN_MASK);
					}
					if ((id == MouseEvent.MOUSE_DRAGGED) && middleButton) {
						y = ((MouseEvent) event[i]).getY();
						int dy = y - y_last;
						if (!reset) {
							transformGroup.getTransform(currXform);
							double z = (-1.0) * dy * z_factor;
							double scale = currXform.getScale() + z;
							if (scale > 0) {
								currXform.setScale(scale);
								transformGroup.setTransform(currXform);
								transformChanged(currXform);
							}
							if (callback != null) {
								callback.transformChanged(
										MouseBehaviorCallback.ZOOM, currXform);
							}
						} else {
							reset = false;
						}
						x_last = x;
						y_last = y;
					}
					if (id == MouseEvent.MOUSE_PRESSED) {
						x_last = ((MouseEvent) event[i]).getX();
						y_last = ((MouseEvent) event[i]).getY();
					} else if (id == MouseEvent.MOUSE_RELEASED) {
						done = true;
					}
				}
			}
		}
		if (!done) {
			wakeupOn(mouseCriterion);
		} else {
			reset = true;
			mouseButton = MouseEvent.BUTTON2_DOWN_MASK;
			postId(doneID);
			wakeupOn(postCriterion);
		}
	}

	/*
	 * Set the y-axis movement multipler with factor.
	 */
	public void setFactor(double factor) {
		z_factor = factor;
	}

	/*
	 * The transformChanged method in the callback class will be called every
	 * time the transform is updated
	 */
	public void setupCallback(MouseBehaviorCallback c) {
		callback = c;
	}

	public void transformChanged(Transform3D transform) {
	}
}
/*
 * Copyright (c) 1996-1998 Sun Microsystems, Inc. All Rights Reserved. Sun
 * grants you ("Licensee") a non-exclusive, royalty free, license to use, modify
 * and redistribute this software in source and binary code form, provided that
 * i) this copyright notice and license appear on all copies of the software;
 * and ii) Licensee does not utilize the software in a manner which is
 * disparaging to Sun. This software is provided "AS IS," without a warranty of
 * any kind. ALL EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES,
 * INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE OR NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL
 * NOT BE LIABLE FOR ANY DAMAGES SUFFXRED BY LICENSEE AS A RESULT OF USING,
 * MODIFYING OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL
 * SUN OR ITS LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR
 * DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES,
 * HOWEVER CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE
 * USE OF OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES. This software is not designed or intended for
 * use in on-line control of aircraft, air traffic, aircraft navigation or
 * aircraft communications; or in the design, construction, operation or
 * maintenance of any nuclear facility. Licensee represents and warrants that it
 * will not use or redistribute the Software for such purposes.
 */
