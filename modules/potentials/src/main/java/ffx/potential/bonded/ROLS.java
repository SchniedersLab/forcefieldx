/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.List;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.vecmath.Color3f;

import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;

/**
 * The ROLS Interace defines "Recursive Over Length Scales" (ROLS) Methods.
 */
public interface ROLS {

    public static final int MaxLengthScale = 5;
    public static int LengthScale = MaxLengthScale;

    public void drawLabel(Canvas3D graphics, J3DGraphics2D g2d, Node node);

    public double[] getCenter(boolean w);

    public ArrayList<ROLS> getList(Class c, ArrayList<ROLS> nodes);

    public long getMSCount(Class c, long count);

    public ROLS getMSNode(Class c);

    public double getMW();

    public void setColor(ColorModel colorModel, Color3f color, Material mat);

    public void setView(ViewModel viewModel, List<BranchGroup> newShapes);

    public void update();
}
