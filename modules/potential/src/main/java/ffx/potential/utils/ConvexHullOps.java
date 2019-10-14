//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.utils;

import com.github.quickhull3d.QuickHull3D;

import ffx.numerics.math.VectorMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Constants;

import java.util.Arrays;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * This ConvexHullOps class uses the QuickHull3D package by John E. Lloyd to
 * construct and operate on 3D convex hulls: the minimal convex polyhedron
 * that contains all points in a set of points. This is especially useful
 * for max-dist operations, as the most distant points in a set are guaranteed
 * to be part of the convex polyhedron.
 *
 * The QuickHull3D package website is at quickhull3d.github.io/quickhull3d/
 * The algorithm it uses is described in Barber, Dobkin, and Huhdanpaa,
 * "The Quickhull Algorithm for Convex Hulls" (ACM Transactions on Mathematical
 * Software, Vol. 22, No. 4, December 1996), and the code is based on the C
 * package qhull.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0.0
 */
public class ConvexHullOps {
    private static final Logger logger = Logger.getLogger(ConvexHullOps.class.getName());

    /**
     * Constructs a convex hull from a MolecularAssembly.
     * @param ma MolecularAssembly to build a convex hull for.
     * @return   A QuickHull3D implementation of convex hulls.
     */
    public static QuickHull3D constructHull(MolecularAssembly ma) {
        Atom[] atoms = ma.getAtomArray();
        return constructHull(atoms);
    }

    /**
     * Constructs a convex hull from a set of atoms.
     * @param atoms Atoms to build a convex hull for.
     * @return      A QuickHull3D implementation of convex hulls.
     */
    public static QuickHull3D constructHull(Atom[] atoms) {
        int nAts = atoms.length;
        double[] xyz = new double[nAts * 3];
        for (int i = 0; i < nAts; i++) {
            Atom at = atoms[i];
            int i3 = 3*i;
            xyz[i3] = at.getX();
            xyz[i3 + 1] = at.getY();
            xyz[i3 + 2] = at.getZ();
        }
        return new QuickHull3D(xyz);
    }

    /**
     * Find the maximum pairwise distance between vertex points on a convex hull.
     * @param qh A QuickHull3D object.
     * @return   Maximum vertex-vertex distance.
     */
    public static double maxDist(QuickHull3D qh) {
        long time = -System.nanoTime();
        int nVerts = qh.getNumVertices();
        assert nVerts > 1;
        double[] vertPoints = new double[3*nVerts];
        qh.getVertices(vertPoints);
        double maxDist = IntStream.range(0, nVerts).
                parallel().
                mapToDouble((int i) -> {
                    double[] xyz = new double[3];
                    System.arraycopy(vertPoints, 3*i, xyz, 0, 3);
                    double mij = 0;
                    for (int j = i + 1; j < nVerts; j++) {
                        double[] xyzJ = new double[3];
                        System.arraycopy(vertPoints, 3*j, xyzJ, 0, 3);
                        double distIJ = VectorMath.dist2(xyz, xyzJ);
                        mij = Math.max(mij, distIJ);
                    }
                    return mij;
                }).
                max().getAsDouble();
        time += System.nanoTime();
        if (time > 1E9) {
            logger.warning(String.format(" Required %12.6g sec to find max distance on a convex hull." +
                    " It may be time to further optimize this!", Constants.NS2SEC * time));
        }
        return maxDist;
    }

    /**
     * UNTESTED: Identifies atoms forming the convex hull.
     * @param qh       A QuickHull3D.
     * @param allAtoms Atoms used in building the QuickHull3D.
     * @return         Atoms forming the convex hull.
     */
    public static Atom[] identifyHullAtoms(QuickHull3D qh, Atom[] allAtoms) {
        int[] indices = qh.getVertexPointIndices();
        return Arrays.stream(indices).
                mapToObj((int i) -> allAtoms[i]).
                toArray(Atom[]::new);
    }
}
