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
package ffx.xray;

import java.util.ArrayList;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;

/**
 * <p>DataContainer interface.</p>
 *
 * @author Timothy D. Fenn
 *
 */
public interface DataContainer {

    /**
     * <p>getAtomArray</p>
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getAtomArray();

    /**
     * <p>getAltResidues</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ArrayList<Residue>> getAltResidues();

    /**
     * <p>getAltMolecules</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ArrayList<Molecule>> getAltMolecules();

    /**
     * <p>getMolecularAssembly</p>
     *
     * @return an array of {@link ffx.potential.bonded.MolecularAssembly}
     * objects.
     */
    public MolecularAssembly[] getMolecularAssembly();

    /**
     * <p>getRefinementModel</p>
     *
     * @return a {@link ffx.xray.RefinementModel} object.
     */
    public RefinementModel getRefinementModel();

    /**
     * <p>getWeight</p>
     *
     * @return the current data weight.
     */
    public double getWeight();

    /**
     * <p>setWeight</p>
     *
     * set the overall weight of the data.
     */
    public void setWeight(double weight);

    /**
     * <p>printOptimizationHeader</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printOptimizationHeader();

    /**
     * <p>printOptimizationUpdate</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printOptimizationUpdate();

    /**
     * <p>printEnergyUpdate</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printEnergyUpdate();
}
