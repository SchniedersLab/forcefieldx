/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.ui.commands;

import java.io.Serializable;

/**
 * The TinkerSystem class is a serializable wrapper that specifies a TINKER
 * system.
 *
 * @author Michael J. Schnieders
 *
 */
public class TinkerSystem implements Serializable {

    private static final long serialVersionUID = 1L;
    public boolean read = true;
    // System definition
    public int numatoms;
    public int numkeys;
    public String file;
    public String forcefield;
    public String[] keywords;
    public double[][] coordinates;
    public int[][] connectivity;
    public int[] types;
    public String[] name;
    public String[] story;
    public double[] charge;
    public double[] mass;
    public int[] atomic;

    /**
     * Constructor that allocates space for a TINKER system
     *
     * @param a The number of atoms
     * @param k The number of keywords
     */
    public TinkerSystem(int a, int k) {
        numatoms = a;
        numkeys = k;
        keywords = new String[k];
        coordinates = new double[3][a];
        connectivity = new int[4][a];
        types = new int[a];
        name = new String[a];
        story = new String[a];
        charge = new double[a];
        mass = new double[a];
        atomic = new int[a];
    }

    /**
     * <p>print</p>
     */
    public void print() {
        System.out.println(this.toString());
    }

    /**
     * <p>toString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toString() {
        return new String("Atoms: " + numatoms + " Keywords: " + numkeys);
    }
}
