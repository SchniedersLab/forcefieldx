/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

/**
 * The BioType class maps PDB identifiers to atom types.
 */
public final class BioType extends BaseType {
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
		this.residue = residue;
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
}
