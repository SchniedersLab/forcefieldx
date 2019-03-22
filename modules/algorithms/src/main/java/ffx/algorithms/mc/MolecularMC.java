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
package ffx.algorithms.mc;

import java.util.logging.Logger;
import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.potential.AssemblyState;
import ffx.potential.MolecularAssembly;

/**
 * The MolecularMC class is a framework to take Monte Carlo steps on a molecular
 * system. It does not implement an MC algorithm, nor does it implement move
 * sets; it is used to evaluate a single MC step with movements defined by
 * implementations of MCMove.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MolecularMC extends BoltzmannMC {

    private static final Logger logger = Logger.getLogger(MolecularMC.class.getName());

    /**
     * The MolecularAssembly to operate on.
     */
    private final MolecularAssembly molecularAssembly;

    /**
     * The potential energy for the molecular assembly.
     */
    private final Potential potential;

    /**
     * Atomic coordinates.
     */
    private double[] x;

    /**
     * Initial state of the MolecularAssembly.
     */
    private AssemblyState initialState;

    /**
     * Constructs a DefaultMC instance with a molecular assembly and its
     * PotentialEnergy. Fancy footwork will be required if we ever need to use
     * multiple assemblies at once.
     *
     * @param molecularAssembly MolecularAssembly to operate on.
     */
    public MolecularMC(MolecularAssembly molecularAssembly) {
        this(molecularAssembly, molecularAssembly.getPotentialEnergy());
    }

    /**
     * Constructs a DefaultMC instance with a molecular assembly and a specific
     * Potential.
     *
     * @param molecularAssembly MolecularAssembly to operate on.
     * @param potential         a {@link ffx.numerics.Potential} object.
     */
    public MolecularMC(MolecularAssembly molecularAssembly, Potential potential) {
        this.molecularAssembly = molecularAssembly;
        this.potential = potential;
    }

    /**
     * Returns the associated MolecularAssembly.
     *
     * @return MolecularAssembly
     */
    public MolecularAssembly getMolecularAssembly() {
        return molecularAssembly;
    }

    /**
     * Returns the associated Potential.
     *
     * @return Potential.
     */
    public Potential getPotential() {
        return potential;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertStep() {
        initialState.revertState();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Calculates the energy at the current state; identical to
     * RotamerOptimization method of same name.
     */
    @Override
    protected double currentEnergy() {
        if (x == null) {
            int nVar = potential.getNumberOfVariables();
            x = new double[nVar * 3];
        }
        try {
            potential.getCoordinates(x);
            return potential.energy(x);
        } catch (ArithmeticException ex) {
            logger.warning(ex.getMessage());
            return 1e100;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Default Metropolis Monte Carlo implementation\nTemperature: ");
        sb.append(getTemperature());
        sb.append(format("\ne1: %10.6f   e2: %10.6f\nMolecular Assembly", getE1(), getE2()));
        sb.append(molecularAssembly.toString()).append("\nPotential: ").append(potential.toString());
        return sb.toString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void storeState() {
        initialState = new AssemblyState(molecularAssembly);
    }
}
