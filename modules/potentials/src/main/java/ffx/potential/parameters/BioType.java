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
package ffx.potential.parameters;

import java.util.Comparator;

/**
 * The BioType class maps PDB identifiers to atom types.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class BioType extends BaseType implements Comparator<String> {

    public final int index;
    public final String PDB;
    public final String residue;
    public final int atomType;

    /**
     * BioType Constructor.
     *
     * @param index
     *            int
     * @param PDB
     *            String
     * @param residue
     *            String
     * @param atomType
     *            int
     */
    public BioType(int index, String PDB, String residue, int atomType) {
        super(ForceField.ForceFieldType.BIOTYPE, Integer.toString(index));
        this.index = index;
        this.PDB = PDB;
        if (residue != null) {
            this.residue = residue.replace(',', ' ');
        } else {
            this.residue = null;
        }
        this.atomType = atomType;
    }

    /**
     * Nicely formatted biotype.
     *
     * @return String
     */
    @Override
    public String toString() {
        return String.format("biotype  %5d  %-4s  %-25s  %5d", index, PDB,
                residue, atomType);
    }

    @Override
    public int compare(String s1, String s2) {

        int t1 = Integer.parseInt(s1);
        int t2 = Integer.parseInt(s2);

        if (t1 < t2) {
            return -1;
        }
        if (t1 > t2) {
            return 1;
        }

        return 0;
    }

    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null || !(other instanceof BioType)) {
            return false;
        }
        BioType bioType = (BioType) other;
        if (bioType.index == this.index) {
            return true;
        }

        return false;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 17 * hash + index;
        return hash;
    }
}
