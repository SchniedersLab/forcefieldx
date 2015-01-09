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
package ffx.ui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.MultipoleType;

public final class Simulation {

    private static final double rotmat[][] = new double[3][3];
    private static final double oldQuad[][] = new double[3][3];

    /**
     * Rotate atomic multipoles into the global frame.
     *
     * @param numberOfMultipoles a int.
     * @param x an array of double.
     * @param multipoleReference an array of int.
     * @param multipoleIndex an array of int.
     * @param frameDefinition an array of
     * {@link ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition}
     * objects.
     * @param localDipole an array of double.
     * @param localQuadrupole an array of double.
     * @param dipole an array of double.
     * @param quadrupole an array of double.
     */
    public static void rotateMultipoles(int numberOfMultipoles, double x[][],
            int multipoleReference[][], int multipoleIndex[],
            MultipoleType.MultipoleFrameDefinition frameDefinition[],
            double localDipole[][], double localQuadrupole[][],
            double dipole[][], double quadrupole[][]) {
        for (int a = 0; a < numberOfMultipoles; a++) {
            double localOrigin[] = x[a];
            double zaxis[] = x[multipoleReference[a][0]];
            double xaxis[] = x[multipoleReference[a][1]];
            int index = multipoleIndex[a];
            MultipoleType.MultipoleFrameDefinition frame = frameDefinition[index];
            double oldDipole[] = localDipole[index];
            double oldQuadRef[] = localQuadrupole[index];
            oldQuad[0][0] = oldQuadRef[0];
            oldQuad[0][1] = oldQuadRef[1];
            oldQuad[0][2] = oldQuadRef[3];
            oldQuad[1][0] = oldQuadRef[1];
            oldQuad[1][1] = oldQuadRef[2];
            oldQuad[1][2] = oldQuadRef[4];
            oldQuad[2][0] = oldQuadRef[3];
            oldQuad[2][1] = oldQuadRef[4];
            oldQuad[2][2] = oldQuadRef[5];
            double newDipole[] = dipole[a];
            double newQuadrupole[] = quadrupole[a];
            newDipole[0] = 0.0;
            newDipole[1] = 0.0;
            newDipole[2] = 0.0;
            newQuadrupole[0] = 0.0;
            newQuadrupole[1] = 0.0;
            newQuadrupole[2] = 0.0;
            newQuadrupole[3] = 0.0;
            newQuadrupole[4] = 0.0;
            newQuadrupole[5] = 0.0;
            VectorMath.diff(zaxis, localOrigin, zaxis);
            VectorMath.norm(zaxis, zaxis);
            VectorMath.diff(xaxis, localOrigin, xaxis);
            // Separate differences between the Z-THEN-X definition
            // and BISECTOR methods for finding the Z elements of the
            // rotation matrix.
            switch (frame) {
                case ZTHENX:
                    rotmat[0][2] = zaxis[0];
                    rotmat[1][2] = zaxis[1];
                    rotmat[2][2] = zaxis[2];
                    break;
                case BISECTOR:
                    VectorMath.norm(xaxis, xaxis);
                    // Take the norm of the sum of unit vectors to find their
                    // bisector.
                    VectorMath.sum(zaxis, xaxis, zaxis);
                    VectorMath.norm(zaxis, zaxis);
                    rotmat[0][2] = zaxis[0];
                    rotmat[1][2] = zaxis[1];
                    rotmat[2][2] = zaxis[2];
                    break;
            }
            // Find the X elements.
            double dot = VectorMath.dot(zaxis, xaxis);
            VectorMath.scalar(zaxis, dot, zaxis);
            VectorMath.diff(xaxis, zaxis, xaxis);
            VectorMath.norm(xaxis, xaxis);
            rotmat[0][0] = xaxis[0];
            rotmat[1][0] = xaxis[1];
            rotmat[2][0] = xaxis[2];
            // Finally the Y elements.
            rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0]
                    * rotmat[2][2];
            rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0]
                    * rotmat[0][2];
            rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0]
                    * rotmat[1][2];
            // Do the rotation.
            for (int i = 0; i < 3; i++) {
                double[] rotmati = rotmat[i];
                for (int j = 0; j < 3; j++) {
                    double[] rotmatj = rotmat[j];
                    newDipole[i] += rotmati[j] * oldDipole[j];
                    if (j > i) {
                        continue;
                    }
                    int newIndex = i + j;
                    for (int k = 0; k < 3; k++) {
                        double oldQuadrupolek[] = oldQuad[k];
                        newQuadrupole[newIndex] += rotmati[k]
                                * (rotmatj[0] * oldQuadrupolek[0] + rotmatj[1]
                                * oldQuadrupolek[1] + rotmatj[2]
                                * oldQuadrupolek[2]);
                    }
                }
            }
        }
    }
    // Coordinates.
    private final int numberOfAtoms;
    private final double[][] x;
    // Local definition of permanent electrostatics.
    private final int localMultipoleIndex[];
    private final MultipoleType.MultipoleFrameDefinition frameDefinition[];
    private final int multipoleReference[][];
    private final double localDipole[][];
    private final double localQuadrupole[][];
    // Global instance of permanent electrostatics.
    private final int numberOfMultipoles;
    private final int iMultipole[];
    private final double charge[];
    private final double dipole[][];
    private final double quadrupole[][];
    // Logging
    private final Logger logger = Logger.getLogger("ffx");

    /**
     * <p>
     * Constructor for Simulation.</p>
     *
     * @param atoms a {@link java.util.List} object.
     */
    public Simulation(List<Atom> atoms) {
        numberOfAtoms = atoms.size();
        x = new double[numberOfAtoms][3];
        // Create a non-redundant Map of local multipole values.
        Map<MultipoleType, Integer> localMultipoleMap = new HashMap<MultipoleType, Integer>();
        int multipoleIndex = 0;
        // Count the number of multipole sites.
        List<Integer> multipoleSites = new ArrayList<Integer>();
        List<Atom> multipoleAtoms = new ArrayList<Atom>();
        for (int i = 0; i < numberOfAtoms; i++) {
            Atom atom = atoms.get(i);
            int atomIndex = atom.getXYZIndex() - 1;
            atom.getXYZ(x[atomIndex]);
            MultipoleType multipoleType = atom.getMultipoleType();
            if (multipoleType != null) {
                multipoleSites.add(i);
                multipoleAtoms.add(atom);
            }
        }
        numberOfMultipoles = multipoleSites.size();
        iMultipole = new int[numberOfMultipoles];
        localMultipoleIndex = new int[numberOfMultipoles];
        multipoleReference = new int[numberOfMultipoles][4];
        charge = new double[numberOfMultipoles];
        dipole = new double[numberOfMultipoles][3];
        quadrupole = new double[numberOfMultipoles][6];
        int index = 0;
        for (Integer i : multipoleSites) {
            iMultipole[index++] = i;
        }
        for (int i = 0; i < numberOfMultipoles; i++) {
            Atom atom = multipoleAtoms.get(i);
            MultipoleType multipoleType = atom.getMultipoleType();
            int atomIndex = atom.getXYZIndex() - 1;
            atom.getXYZ(x[atomIndex]);
            int referenceIndex = 0;
            for (Atom referenceAtom : atom.getMultipoleReferenceSites()) {
                multipoleReference[i][referenceIndex++] = referenceAtom
                        .getXYZIndex() - 1;
            }
            if (!localMultipoleMap.containsKey(multipoleType)) {
                localMultipoleIndex[i] = multipoleIndex;
                localMultipoleMap.put(multipoleType, multipoleIndex++);
            } else {
                localMultipoleIndex[i] = localMultipoleMap.get(multipoleType);
            }
            charge[atomIndex] = multipoleType.charge;
        }
        // Store the local multipole values into an array.
        int mapSize = localMultipoleMap.size();
        localDipole = new double[mapSize][3];
        localQuadrupole = new double[mapSize][6];
        frameDefinition = new MultipoleType.MultipoleFrameDefinition[mapSize];
        for (MultipoleType multipoleType : localMultipoleMap.keySet()) {
            index = localMultipoleMap.get(multipoleType);
            // Store frame definition
            frameDefinition[index] = multipoleType.frameDefinition;
            // Store Dipoles
            localDipole[index][0] = multipoleType.dipole[0];
            localDipole[index][1] = multipoleType.dipole[1];
            localDipole[index][2] = multipoleType.dipole[2];
            // Store non-redundant quadrupole values
            // xx, yx, yy, zx, zy, zz
            localQuadrupole[index][0] = multipoleType.quadrupole[0][0];
            localQuadrupole[index][1] = multipoleType.quadrupole[1][0];
            localQuadrupole[index][2] = multipoleType.quadrupole[1][1];
            localQuadrupole[index][3] = multipoleType.quadrupole[2][0];
            localQuadrupole[index][4] = multipoleType.quadrupole[2][1];
            localQuadrupole[index][5] = multipoleType.quadrupole[2][2];
        }
        long start = System.currentTimeMillis();
        rotateMultipoles(numberOfAtoms, x, multipoleReference,
                this.localMultipoleIndex, frameDefinition, localDipole,
                localQuadrupole, dipole, quadrupole);
        long done = System.currentTimeMillis();
        logger.info(String.format("Rotation time: %d (msec)", done - start));
    }
}
