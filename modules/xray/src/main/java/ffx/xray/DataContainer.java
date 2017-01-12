/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.util.ArrayList;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;

/**
 * <p>
 * DataContainer interface.</p>
 *
 * @author Timothy D. Fenn
 */
public interface DataContainer {

    /**
     * <p>
     * getAtomArray</p>
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getAtomArray();

    /**
     * <p>
     * getActiveAtomArray</p>
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getActiveAtomArray();

    /**
     * <p>
     * getAltResidues</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ArrayList<Residue>> getAltResidues();

    /**
     * <p>
     * getAltMolecules</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ArrayList<Molecule>> getAltMolecules();

    /**
     * <p>
     * getMolecularAssemblies</p>
     *
     * @return an array of {@link ffx.potential.MolecularAssembly} objects.
     */
    public MolecularAssembly[] getMolecularAssemblies();

    /**
     * <p>
     * getRefinementModel</p>
     *
     * @return a {@link ffx.xray.RefinementModel} object.
     */
    public RefinementModel getRefinementModel();

    /**
     * <p>
     * getWeight</p>
     *
     * @return the current data weight.
     */
    public double getWeight();

    /**
     * <p>
     * setWeight</p>
     *
     * @param weight set the overall weight of the data.
     */
    public void setWeight(double weight);

    /**
     * <p>
     * printOptimizationHeader</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printOptimizationHeader();

    /**
     * <p>
     * printOptimizationUpdate</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printOptimizationUpdate();

    /**
     * <p>
     * printEnergyUpdate</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String printEnergyUpdate();
}
