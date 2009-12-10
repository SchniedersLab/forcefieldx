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

import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;

import com.sun.j3d.utils.picking.PickResult;

/**
 * The PickOrbitBehavior class implements a mouse orbit behavior.
 */
public class PickOrbitBehavior extends PickMouseBehavior implements
		MouseBehaviorCallback {
	public MouseOrbit orbit;
	private PickingCallback callback = null;
	private TransformGroup currentTG;

	public PickOrbitBehavior(BranchGroup root, Canvas3D canvas, Bounds bounds,
			TransformGroup VPTG, int pickMode) {
		super(canvas, root, bounds);
		orbit = new MouseOrbit(MouseOrbit.MANUAL_WAKEUP, VPTG);
		orbit.setTransformGroup(currGrp);
		currGrp.addChild(orbit);
		orbit.setSchedulingBounds(bounds);
		setSchedulingBounds(bounds);
		pickCanvas.setMode(pickMode);
	}

	/*
	 * Return the pickMode component of this PickTranslateBehavior.
	 */
	public int getPickMode() {
		return pickCanvas.getMode();
	}

	/*
	 * Sets the pickMode component of this PickTranslateBehavior to the value of
	 * the passed pickMode. @param pickMode the pickMode to be copied.
	 */
	public void setPickMode(int pickMode) {
		pickCanvas.setMode(pickMode);
	}

	public void setTransformGroups(TransformGroup StarTG, TransformGroup VPTG) {
		orbit.setTransformGroups(StarTG, VPTG);
	}

	/*
	 * Register the class @param callback to be called each time the picked
	 * object moves
	 */
	public void setupCallback(PickingCallback c) {
		callback = c;
		if (callback == null) {
			orbit.setupCallback(null);
		} else {
			orbit.setupCallback(this);
		}
	}

	/*
	 * Callback method from MouseOrbit This is used when the Picking callback is
	 * enabled
	 */
	public void transformChanged(int type, Transform3D transform) {
		callback.transformChanged(PickingCallback.ORBIT, currentTG);
	}

	public void transformClicked(int type, Transform3D transform) {
		callback.transformClicked(PickingCallback.ORBIT, currentTG);
	}

	public void transformDoubleClicked(int type, Transform3D transform) {
		callback.transformDoubleClicked(PickingCallback.ORBIT, currentTG);
	}

	/*
	 * Update the scene to manipulate any nodes. This is not meant to be called
	 * by users. Behavior automatically calls this. You can call this only if
	 * you know what you are doing.
	 * @param xpos Current mouse X pos. @param ypos Current mouse Y pos.
	 */
	public void updateScene(int xpos, int ypos) {
		TransformGroup tg = null;
		if (mevent.isMetaDown() && !mevent.isAltDown()) {
			pickCanvas.setShapeLocation(xpos, ypos);
			PickResult r = pickCanvas.pickClosest();
			if (r != null) {
				tg = (TransformGroup) r.getNode(PickResult.TRANSFORM_GROUP);
				if ((tg != null)
						&& (tg
								.getCapability(TransformGroup.ALLOW_TRANSFORM_READ))
						&& (tg
								.getCapability(TransformGroup.ALLOW_TRANSFORM_WRITE))) {
					orbit.setTransformGroup(tg);
					orbit.wakeup();
					currentTG = tg;
					if (callback != null) {
						callback.transformClicked(PickingCallback.ORBIT,
								currentTG);
					}
				}
			} else if (callback != null) {
				callback.transformChanged(PickingCallback.NO_PICK, null);
			}
		}
	}
}
