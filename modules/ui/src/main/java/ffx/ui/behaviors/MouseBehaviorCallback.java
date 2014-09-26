/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.ui.behaviors;

import javax.media.j3d.Transform3D;

/**
 * The MouseBehaviorCallback interface is implemented by classes that want to
 * receive callbacks when transforms are updated.
 *
 * @author Michael J. Schnieders
 *
 */
public interface MouseBehaviorCallback {

    /**
     * Constant <code>ROTATE=0</code>
     */
    public final static int ROTATE = 0;
    /**
     * Constant <code>TRANSLATE=1</code>
     */
    public final static int TRANSLATE = 1;
    /**
     * Constant <code>ZOOM=2</code>
     */
    public final static int ZOOM = 2;
    /**
     * Constant <code>SELECTION=4</code>
     */
    public final static int SELECTION = 4;
    /**
     * Constant <code>PROPERTIES=5</code>
     */
    public final static int PROPERTIES = 5;
    /**
     * Constant <code>ORBIT=6</code>
     */
    public final static int ORBIT = 6;

    /*
     * Classes implementing this interface that are registered with one of the
     * MouseBehaviors will be called every time the behavior updates the
     * Transform @param type will be one of ROTATE, TRANSLATE or ZOOM
     */
    /**
     * <p>
     * transformChanged</p>
     *
     * @param type a int.
     * @param transform a {@link javax.media.j3d.Transform3D} object.
     */
    public void transformChanged(int type, Transform3D transform);

    /**
     * <p>
     * transformClicked</p>
     *
     * @param type a int.
     * @param transform a {@link javax.media.j3d.Transform3D} object.
     */
    public void transformClicked(int type, Transform3D transform);

    /**
     * <p>
     * transformDoubleClicked</p>
     *
     * @param type a int.
     * @param transform a {@link javax.media.j3d.Transform3D} object.
     */
    public void transformDoubleClicked(int type, Transform3D transform);
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
