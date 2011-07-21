/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

package ffx.xray;

import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;

/**
 *
 * @author Tim Fenn
 */
public interface DataContainer {
    public Atom[] getAtomArray();

    public ArrayList<ArrayList<Residue>> getAltResidues();

    public ArrayList<ArrayList<Molecule>> getAltMolecules();

    public MolecularAssembly[] getMolecularAssembly();

    public RefinementModel getRefinementModel();

    public double getWeight();

    public String printOptimizationHeader();

    public String printOptimizationUpdate();

    public String printEnergyUpdate();
}
