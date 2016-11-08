/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray;

/**
 * <p>
 * RealSpaceRefinementData class.</p>
 *
 * @author Timothy D. Fenn
 */
public class RealSpaceRefinementData {

    private final int[] origin;
    private final int[] extent;
    private final int[] ni;
    private double[] data;
    private double densityScore;
    private boolean periodic;

    /**
     * <p>
     * Constructor for RealSpaceRefinementData.</p>
     *
     */
    public RealSpaceRefinementData() {
        origin = new int[3];
        extent = new int[3];
        ni = new int[3];
        periodic = false;
    }

    /**
     * <p>
     * getDataIndex</p>
     *
     * @param x a int.
     * @param y a int.
     * @param z a int.
     * @return a double.
     */
    public double getDataIndex(int x, int y, int z) {
        int index = x + extent[0] * (y + extent[1] * z);
        return data[index];
    }

    /**
     * @return the origin
     */
    public int[] getOrigin() {
        return origin;
    }

    /**
     * @return the extent
     */
    public int[] getExtent() {
        return extent;
    }

    /**
     * @return the ni
     */
    public int[] getNi() {
        return ni;
    }

    /**
     * @return the data
     */
    public double[] getData() {
        return data;
    }

    /**
     * @param data the data to set
     */
    public void setData(double[] data) {
        this.data = data;
    }

    /**
     * @return the densityScore
     */
    public double getDensityScore() {
        return densityScore;
    }

    /**
     * @param densityScore the densityScore to set
     */
    public void setDensityScore(double densityScore) {
        this.densityScore = densityScore;
    }

    /**
     * @return the periodic
     */
    public boolean isPeriodic() {
        return periodic;
    }

    /**
     * @param periodic the periodic to set
     */
    public void setPeriodic(boolean periodic) {
        this.periodic = periodic;
    }

    public void setOrigin(int x, int y, int z) {
        origin[0] = x;
        origin[1] = y;
        origin[2] = z;
    }

    public void setExtent(int x, int y, int z) {
        extent[0] = x;
        extent[1] = y;
        extent[2] = z;
    }

    public void setNI(int x, int y, int z) {
        ni[0] = x;
        ni[1] = y;
        ni[2] = z;
    }
}
