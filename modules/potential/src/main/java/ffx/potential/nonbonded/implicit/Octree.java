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
package ffx.potential.nonbonded.implicit;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.jogamp.vecmath.Point3d;

import ffx.potential.bonded.Atom;

/**
 * <p>Initial implementation of an Octree decomposition.</p>
 *
 * @author Stephen D. LuCore
 * @since 1.0
 */
public class Octree {

    private static final Logger logger = Logger.getLogger(Octree.class.getName());
    private static boolean warnedAboutSplit = false;
    private static boolean leaveStraddlersInParent = false;

    private static int maxAtomsPerVolume = 20;
    private static int maxTreeDepth = 10;

    private int depth;
    private Octree[] children;
    private List<Atom> contents = new ArrayList<>();
    private Point3d corner;
    private double edgeLength;

    /**
     * <p>Constructor for Octree.</p>
     *
     * @param depth      a int.
     * @param corner     a {@link org.jogamp.vecmath.Point3d} object.
     * @param edgeLength a double.
     */
    public Octree(int depth, Point3d corner, double edgeLength) {
        this.depth = depth;
        this.corner = corner;
        this.edgeLength = edgeLength;
    }

    /**
     * <p>debugPrintStats.</p>
     *
     * @param writePartitionFile a boolean.
     * @param partFile           a {@link java.io.File} object.
     */
    public void debugPrintStats(boolean writePartitionFile, File partFile) {
        List<Octree> nodes = new ArrayList<>();
        List<Octree> leaves = new ArrayList<>();
        debugFindNodes(nodes);
        int maxDepth = 0;
        int minContents = Integer.MAX_VALUE, maxContents = 0;
        double minEdgeLength = Double.MAX_VALUE, maxEdgeLength = 0;
        double minConc = Double.MAX_VALUE, maxConc = 0;
        for (Octree node : nodes) {
            if (node.getChildren() == null) {
                leaves.add(node);
                int nodeContents = node.getContents().size();
                double nodeEdge = node.getEdgeLength();
                double nodeVolume = Math.pow(nodeEdge, 3);
                double nodeConc = nodeContents / nodeVolume;
                maxDepth = (node.getDepth() > maxDepth) ? node.getDepth() : maxDepth;
                minContents = (nodeContents < minContents) ? nodeContents : minContents;
                maxContents = (nodeContents > maxContents) ? nodeContents : maxContents;
                minEdgeLength = (nodeEdge < minEdgeLength) ? nodeEdge : minEdgeLength;
                maxEdgeLength = (nodeEdge > maxEdgeLength) ? nodeEdge : maxEdgeLength;
                minConc = (nodeConc < minConc) ? nodeConc : minConc;
                maxConc = (nodeConc > maxConc) ? nodeConc : maxConc;
            }
        }
        logger.info(" Octree Leaf Stats:");
        logger.info(format("    Max Depth:        %10d", maxDepth));
        logger.info(format("    Min/Max Atoms:    %10d", minContents));
        logger.info(format("                      %10d", maxContents));
        logger.info(format("    Min/Max Volume:   %10.2g", Math.pow(minEdgeLength, 3)));
        logger.info(format("                      %10.2g", Math.pow(maxEdgeLength, 3)));
        logger.info(format("    Min/Max Conc:     %10.2g", minConc));
        logger.info(format("                      %10.2g", maxConc));
        logger.info(" ");
        if (writePartitionFile) {
            try {
                logger.info(format(" Writing atom partition file to %s", partFile.getName()));
                BufferedWriter bw = new BufferedWriter(new FileWriter(partFile));
                int i = 1;
                for (Octree leaf : leaves) {
                    Point3d p = leaf.getCorner();
                    bw.write(format("Leaf %d:  (%6.2f, %6.2f, %6.2f) x %6.2f\n",
                            i++, p.x, p.y, p.z, leaf.getEdgeLength()));
                    for (Atom atom : leaf.getContents()) {
                        bw.write(format("  %s\n", atom.toString()));
                    }
                }
                bw.flush();
                bw.close();
            } catch (IOException ex) {
                logger.warning("Exception writing atom partition file.");
            }
        }
    }

    private void debugFindNodes(List<Octree> nodes) {
        nodes.add(this);
        if (children != null) {
            for (Octree child : children) {
                child.debugFindNodes(nodes);
            }
        }
    }

    /**
     * <p>addAtoms.</p>
     *
     * @param atoms a {@link java.util.List} object.
     */
    public void addAtoms(List<Atom> atoms) {
        for (Atom atom : atoms) {
            addAtom(atom);
        }
    }

    /**
     * <p>Setter for the field <code>maxAtomsPerVolume</code>.</p>
     *
     * @param max a int.
     */
    public void setMaxAtomsPerVolume(int max) {
        maxAtomsPerVolume = max;
    }

    /**
     * <p>Setter for the field <code>maxTreeDepth</code>.</p>
     *
     * @param max a int.
     */
    public void setMaxTreeDepth(int max) {
        maxTreeDepth = max;
    }

    /**
     * <p>Setter for the field <code>leaveStraddlersInParent</code>.</p>
     *
     * @param set a boolean.
     */
    public void setLeaveStraddlersInParent(boolean set) {
        leaveStraddlersInParent = set;
    }

    private void split() {
        double sel = edgeLength / 2;
        double x = corner.x;    // for Java3d 1.3 compatibility, use corner.getX() in 1.5+
        double y = corner.y;
        double z = corner.z;
        children = new Octree[8];
        children[0] = new Octree(depth + 1,
                new Point3d(x, y, z), sel);
        children[1] = new Octree(depth + 1,
                new Point3d(x + sel, y, z), sel);
        children[2] = new Octree(depth + 1,
                new Point3d(x, y + sel, z), sel);
        children[3] = new Octree(depth + 1,
                new Point3d(x, y, z + sel), sel);
        children[4] = new Octree(depth + 1,
                new Point3d(x + sel, y + sel, z), sel);
        children[5] = new Octree(depth + 1,
                new Point3d(x + sel, y, z + sel), sel);
        children[6] = new Octree(depth + 1,
                new Point3d(x, y + sel, z + sel), sel);
        children[7] = new Octree(depth + 1,
                new Point3d(x + sel, y + sel, z + sel), sel);
    }

    private void addAtom(Atom atom) {
        if (children != null) {
            if (leaveStraddlersInParent) {
                int index = findAtomIndex(atom);
                if (index != -1) {
                    children[index].addAtom(atom);
                    return;
                }
            } else {
                boolean indices[] = partitionAtom(atom);
                for (int i = 0; i < 8; i++) {
                    if (indices[i]) {
                        children[i].addAtom(atom);
                    }
                }
                return;
            }
        }
        contents.add(atom);
        if (contents.size() > maxAtomsPerVolume && depth < maxTreeDepth) {
            if (children != null) {
                if (!warnedAboutSplit) {
                    logger.warning("Octree couldn't split due to straddlers; maxAtomsPerVolume may be violated.");
                    warnedAboutSplit = true;
                }
                return;
            }
            split();
            for (int i = 0; i < contents.size(); ) {    // NOTE: intentional lack of auto-increment
                if (leaveStraddlersInParent) {
                    int index = findAtomIndex(contents.get(i));
                    if (index != -1) {
                        children[index].addAtom(contents.remove(i));
                    } else {
                        i++;
                    }
                } else {
                    boolean[] indices = partitionAtom(contents.get(i));
                    Atom current = contents.remove(i);
                    for (int j = 0; j < 8; j++) {
                        if (indices[j]) {
                            children[j].addAtom(current);
                        }
                    }
                }
            }
        }
    }

    /**
     * Finds the index of the single child octree to which this atom belongs;
     * returns -1 for straddlers.
     *
     * @param atom query atom
     * @return index of children[] octree inside which this atom fits completely
     */
    private int findAtomIndex(Atom atom) {
        int index = -1;
        boolean b[] = partitionAtom(atom);
        for (int i = 0; i < 8; i++) {
            if (b[i]) {
                if (index != -1) {
                    return -1;
                }
                index = i;
            }
        }
        return index;
    }

    /**
     * Finds the indices of the children octrees to which this atom belongs
     * (wholly or partially).
     *
     * @param atom query atom
     * @return boolean array such that if (b[i] == true) then children[i]
     * touches this atom
     */
    private boolean[] partitionAtom(Atom atom) {
        boolean b[] = new boolean[8];
        for (int i = 0; i < 8; i++) {
            b[i] = true;
        }
        double midX = corner.x + (edgeLength / 2);  // for Java3d 1.3 compatibility, use corner.getX() in 1.5+
        double midY = corner.y + (edgeLength / 2);
        double midZ = corner.z + (edgeLength / 2);
        double vdwr = atom.getVDWR();
        double atomX = atom.getX();
        double atomY = atom.getY();
        double atomZ = atom.getZ();
        boolean leftSide = (atomX + vdwr < midX);
        boolean rightSide = (atomX - vdwr > midX);
        boolean topSide = (atomY + vdwr < midY);
        boolean botSide = (atomY - vdwr > midY);
        boolean frontSide = (atomZ + vdwr < midZ);
        boolean backSide = (atomZ - vdwr > midZ);
        if (leftSide) {
            b[1] = false;
            b[4] = false;
            b[5] = false;
            b[7] = false;
        } else if (rightSide) {
            b[0] = false;
            b[2] = false;
            b[3] = false;
            b[6] = false;
        }
        if (topSide) {
            b[2] = false;
            b[4] = false;
            b[6] = false;
            b[7] = false;
        } else if (botSide) {
            b[0] = false;
            b[1] = false;
            b[3] = false;
            b[5] = false;
        }
        if (frontSide) {
            b[3] = false;
            b[5] = false;
            b[6] = false;
            b[7] = false;
        } else if (backSide) {
            b[0] = false;
            b[1] = false;
            b[2] = false;
            b[4] = false;
        }
        return b;
    }

    /**
     * <p>Getter for the field <code>depth</code>.</p>
     *
     * @return a int.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * <p>Getter for the field <code>children</code>.</p>
     *
     * @return an array of {@link Octree} objects.
     */
    public Octree[] getChildren() {
        return children;
    }

    /**
     * <p>Getter for the field <code>contents</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<Atom> getContents() {
        return contents;
    }

    /**
     * <p>Getter for the field <code>corner</code>.</p>
     *
     * @return a {@link org.jogamp.vecmath.Point3d} object.
     */
    public Point3d getCorner() {
        return corner;
    }

    /**
     * <p>Getter for the field <code>edgeLength</code>.</p>
     *
     * @return a double.
     */
    private double getEdgeLength() {
        return edgeLength;
    }

}
