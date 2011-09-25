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

import java.awt.event.MouseEvent;

import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;

import com.sun.j3d.utils.picking.PickResult;

/**
 * The PickTranslateBehavior class implements a translation behavior on a picked
 * scenegraph object.
 *
 * @author schnied
 * @version $Id: $
 */
public class PickTranslateBehavior extends PickMouseBehavior implements
		MouseBehaviorCallback {
	public MouseTranslate translate;
	private PickingCallback callback = null;
	private TransformGroup currentTG;

	/**
	 * <p>Constructor for PickTranslateBehavior.</p>
	 *
	 * @param root a {@link javax.media.j3d.BranchGroup} object.
	 * @param canvas a {@link javax.media.j3d.Canvas3D} object.
	 * @param bounds a {@link javax.media.j3d.Bounds} object.
	 * @param VPTG a {@link javax.media.j3d.TransformGroup} object.
	 * @param pickMode a int.
	 */
	public PickTranslateBehavior(BranchGroup root, Canvas3D canvas,
			Bounds bounds, TransformGroup VPTG, int pickMode) {
		super(canvas, root, bounds);
		translate = new MouseTranslate(MouseBehavior.MANUAL_WAKEUP, VPTG);
		translate.setTransformGroup(currGrp);
		translate.setFactor(0.1);
		currGrp.addChild(translate);
		translate.setSchedulingBounds(bounds);
		setSchedulingBounds(bounds);
		pickCanvas.setMode(pickMode);
	}

	/*
	 * Return the pickMode component of this PickTranslateBehavior.
	 */
	/**
	 * <p>getPickMode</p>
	 *
	 * @return a int.
	 */
	public int getPickMode() {
		return pickCanvas.getMode();
	}

	/*
	 * Sets the pickMode component of this PickTranslateBehavior to the value of
	 * the passed pickMode. @param pickMode the pickMode to be copied.
	 */
	/**
	 * <p>setPickMode</p>
	 *
	 * @param pickMode a int.
	 */
	public void setPickMode(int pickMode) {
		pickCanvas.setMode(pickMode);
	}

	/*
	 * Register the class @param callback to be called each time the picked
	 * object moves
	 */
	/**
	 * <p>setupCallback</p>
	 *
	 * @param callback a {@link ffx.ui.behaviors.PickingCallback} object.
	 */
	public void setupCallback(PickingCallback callback) {
		this.callback = callback;
		if (callback == null) {
			translate.setupCallback(null);
		} else {
			translate.setupCallback(this);
		}
	}

	/** {@inheritDoc} */
	public void transformChanged(int type, Transform3D transform) {
		callback.transformChanged(PickingCallback.TRANSLATE, currentTG);
	}

	/** {@inheritDoc} */
	public void transformClicked(int type, Transform3D transform) {
		callback.transformClicked(PickingCallback.TRANSLATE, currentTG);
	}

	/** {@inheritDoc} */
	public void transformDoubleClicked(int type, Transform3D transform) {
		callback.transformDoubleClicked(PickingCallback.TRANSLATE, currentTG);
	}

	/*
	 * Update the scene to manipulate any nodes. This is not meant to be called
	 * by users. Behavior automatically calls this. You can call this only if
	 * you know what you are doing.
	 * @param xpos Current mouse X pos. @param ypos Current mouse Y pos.
	 */
	/** {@inheritDoc} */
	public void updateScene(int xpos, int ypos) {
		if ((mevent.getModifiersEx() & MouseEvent.BUTTON3_DOWN_MASK) == MouseEvent.BUTTON3_DOWN_MASK) {
			pickCanvas.setShapeLocation(xpos, ypos);
			PickResult r = pickCanvas.pickClosest();
			if (r != null) {
				if (callback != null) {
					callback.transformChanged(PickingCallback.NO_PICK, null);
				}
			}
		}
	}
}
/*
 * Copyright (c) 1996-1998 Sun MicroFSystems, Inc. All Rights Reserved. Sun
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
