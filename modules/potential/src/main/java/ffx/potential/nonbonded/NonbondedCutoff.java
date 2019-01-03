/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential.nonbonded;

/**
 * This class contains fields and methods for maintaining details of non-bonded
 * cutoffs.
 *
 * @author Michael J. Schnieders
 *
 */
public class NonbondedCutoff {

    /**
     * Non-bonded Cutoff constructor.
     *
     * @param off All vdW interactions are 0 at the distance <code>off</code>.
     * @param cut At the distance <code>cut</code>, a multiplicative switch
     * begins to be applied.
     * @param buff A buffer added to the cut-off distance <code>off</code> to
     * define neighbors included when collecting Verlet lists.
     */
    public NonbondedCutoff(double off, double cut, double buff) {
        this.cut = cut;
        this.cut2 = cut * cut;
        this.off = off;
        this.off2 = off * off;
        this.buff = buff;
    }

    /**
     * Returns a NonbondedCutoff that does not cut off anything.
     *
     * @return No-cutoff cutoff.
     */
    public static NonbondedCutoff noCutoffFactory() {
        return new NonbondedCutoff(Double.MAX_VALUE, Double.MAX_VALUE, 0.0);
    }

    /**
     * At the distance "cut", a multiplicative switch begins to be applied.
     */
    public final double cut;
    /**
     * The distance cut squared.
     */
    public final double cut2;
    /**
     * All vdW interactions are 0 at the distance "off".
     */
    public final double off;
    /**
     * The distance off squared.
     */
    public final double off2;
    /**
     * A buffer added to the cut-off distance <code>off</code> to define
     * neighbors included when collecting Verlet lists.
     */
    public final double buff;

}
