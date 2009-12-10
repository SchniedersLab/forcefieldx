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

import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.WakeupCriterion;
import javax.media.j3d.WakeupOnAWTEvent;

/**
 * The MouseSelection class implements a mouse selection behavior.
 */
public class MouseSelection extends MouseBehavior {
	double x_angle, y_angle;
	double x_factor = .03;
	double y_factor = .03;
	private MouseBehaviorCallback callback = null;

	public MouseSelection(int flags, TransformGroup VPTG) {
		super(flags, VPTG);
	}

	/*
	 * Return the x-axis movement multipler.
	 */
	public double getXFactor() {
		return x_factor;
	}

	/*
	 * Return the y-axis movement multipler.
	 */
	public double getYFactor() {
		return y_factor;
	}

	public void initialize() {
		super.initialize();
		x_angle = 0;
		y_angle = 0;
		if ((flags & INVERT_INPUT) == INVERT_INPUT) {
			invert = true;
			x_factor *= -1;
			y_factor *= -1;
		}
	}

	public void processStimulus(Enumeration criteria) {
		WakeupCriterion wakeup;
		AWTEvent[] event;
		int id;
		while (criteria.hasMoreElements()) {
			wakeup = (WakeupCriterion) criteria.nextElement();
			if (wakeup instanceof WakeupOnAWTEvent) {
				event = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
				for (int i = 0; i < event.length; i++) {
					processMouseEvent((MouseEvent) event[i]);
					if (((buttonPress) && ((flags & MANUAL_WAKEUP) == 0))
							|| ((wakeUp) && ((flags & MANUAL_WAKEUP) != 0))) {
						id = event[i].getID();
						if ((id == MouseEvent.MOUSE_DRAGGED)) {
							if (!reset) {
								transformChanged(currXform);
								if (callback != null) {
									callback.transformChanged(
											MouseBehaviorCallback.SELECTION,
											currXform);
								}
							} else {
								reset = false;
							}
							x_last = ((MouseEvent) event[i]).getX();
							y_last = ((MouseEvent) event[i]).getY();
						}
						if (id == MouseEvent.MOUSE_PRESSED) {
							x_last = ((MouseEvent) event[i]).getX();
							y_last = ((MouseEvent) event[i]).getY();
						}
					}
				}
			}
		}
		wakeupOn(mouseCriterion);
	}

	/*
	 * Set the x-axis amd y-axis movement multipler with factor.
	 */
	public void setFactor(double factor) {
		x_factor = y_factor = factor;
	}

	/*
	 * Set the x-axis amd y-axis movement multipler with xFactor and yFactor
	 * respectively.
	 */
	public void setFactor(double xFactor, double yFactor) {
		x_factor = xFactor;
		y_factor = yFactor;
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
