/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.algorithms.mc;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.ResidueState;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * The CoordShakeMove class implements a simplistic atomic coordinate shake. 
 * At present, simply adds a random number from a normal distribution to each
 * Cartesian coordinate; in the future, will use a move in polar coordinates.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 *
 */
public class CoordShakeMove implements MCMove {
    private Atom[] atoms;
    private double[][] originalCoords;
    private double sigma = 0.1;
    private NormalDistribution dist;
    
    
    public CoordShakeMove(MolecularAssembly assembly) {
        this(assembly.getAtomArray());
    }
    
    public CoordShakeMove(Atom[] atoms) {
        int nAtoms = atoms.length;
        this.atoms = new Atom[nAtoms];
        System.arraycopy(atoms, 0, this.atoms, 0, nAtoms);
        originalCoords = ResidueState.storeAtomicCoordinates(this.atoms);
        dist = new NormalDistribution(0, sigma);
    }
    
    public void setSigma(double sigma) {
        this.sigma = sigma;
        dist = new NormalDistribution(0, sigma);
    }
    
    public void setAtoms(Atom[] atoms) {
        int nAtoms = atoms.length;
        this.atoms = new Atom[nAtoms];
        System.arraycopy(atoms, 0, this.atoms, 0, nAtoms);
        originalCoords = ResidueState.storeAtomicCoordinates(this.atoms);
    }

    @Override
    public double move() {
        originalCoords = ResidueState.storeAtomicCoordinates(this.atoms);
        int nAtoms = atoms.length;

        // Perform the shake.
        // At some point, I really should change this to a move in polar coordinates.
        for (int i = 0; i < nAtoms; i++) {
            double[] xyz = new double[3];
            atoms[i].getXYZ(xyz);
            for (int j = 0; j < 3; j++) {
                xyz[j] += dist.sample();
            }
            atoms[i].setXYZ(xyz);
        }
        return 0.0;
    }

    @Override
    public double revertMove() {
        ResidueState.revertAtomicCoordinates(atoms, originalCoords);
        return 0.0;
    }

    @Override
    public double getEcorrection() {
        return 0.0;
    }
    
    @Override
    public String getDescription() {
        return this.toString();
    }
    
    @Override
    public String toString() {
        return String.format("Coordinate randomization: normal distribution with sigma %10.6f", sigma);
    }
}
