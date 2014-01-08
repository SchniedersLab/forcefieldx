/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.potential.nonbonded;

import ffx.potential.bonded.Atom;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point3d;
import java.util.logging.Logger;
import static java.lang.String.format;

/**
 *
 * @author Salvo
 */
public class Octree {
    
    private static final Logger logger = Logger.getLogger(Octree.class.getName());
    private static boolean warnedAboutSplit = false;
    
    private static int maxAtomsPerVolume = 20;
    private static int maxTreeDepth = 10;
    
    private int depth;
    private Octree children[];
    private List<Atom> contents = new ArrayList<>();
    private Point3d corner;
    private double edgeLength;
    
    public Octree(int depth, Point3d corner, double edgeLength) {
        this.depth = depth;
        this.corner = corner;
        this.edgeLength = edgeLength;
    }

    public void debugPrintStats(boolean writePartitionFile, File partFile) {
        List<Octree> nodes = new ArrayList<>();
        List<Octree> leaves = new ArrayList<>();
        debugFindNodes(nodes);
        int maxDepth = 0;
//        for (Octree node : nodes) {
//            maxDepth = (node.getDepth() > maxDepth) ? node.getDepth() : maxDepth;
//        }
        int minContents = Integer.MAX_VALUE, maxContents = 0;
        double minEdgeLength = Double.MAX_VALUE, maxEdgeLength = 0;
        double minVolume = Double.MAX_VALUE, maxVolume = 0;
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
                maxEdgeLength = (nodeEdge > maxEdgeLength) ? nodeEdge: maxEdgeLength;
                minConc = (nodeConc < minConc) ? nodeConc : minConc;
                maxConc = (nodeConc > maxConc) ? nodeConc : maxConc;
            }
        }
        logger.info(format(" Octree Leaf Stats:"));
        logger.info(format("    Max Depth:        %10d", maxDepth));
        logger.info(format("    Min/Max Atoms:    %10d", minContents));
        logger.info(format("                      %10d", maxContents));
        logger.info(format("    Min/Max Volume:   %10.2g", Math.pow(minEdgeLength,3)));
        logger.info(format("                      %10.2g", Math.pow(maxEdgeLength,3)));
        logger.info(format("    Min/Max Conc:     %10.2g", minConc));
        logger.info(format("                      %10.2g", maxConc));
        logger.info(format(" "));
        if (writePartitionFile) {
            try {
                logger.info(format(" Writing atom partition file to %s", partFile.getName()));
                BufferedWriter bw = new BufferedWriter(new FileWriter(partFile));
                int i = 1;
                for (Octree leaf : leaves) {
                    Point3d p = leaf.getCorner();
                    bw.write(format("Leaf %d:  (%6.2f, %6.2f, %6.2f) x %6.2f\n",
                            i++, p.getX(), p.getY(), p.getZ(), leaf.getEdgeLength()));
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
    
    public void addAtoms(List<Atom> atoms) {
        for (Atom atom : atoms) {
            addAtom(atom);
        }
    }
    
    public void setMaxAtomsPerVolume(int max) {
        maxAtomsPerVolume = max;
    }
    
    public void setMaxTreeDepth(int max) {
        maxTreeDepth = max;
    }
    
    private void split() {
        double sel = edgeLength / 2;
        double x = corner.getX();
        double y = corner.getY();
        double z = corner.getZ();
        children = new Octree[8];
        children[0] = new Octree(depth+1,
                new Point3d(x,y,z), sel);
        children[1] = new Octree(depth+1,
                new Point3d(x+sel, y, z), sel);
        children[2] = new Octree(depth+1,
                new Point3d(x, y+sel, z), sel);
        children[3] = new Octree(depth+1,
                new Point3d(x, y, z+sel), sel);
        children[4] = new Octree(depth+1,
                new Point3d(x+sel, y+sel, z), sel);
        children[5] = new Octree(depth+1,
                new Point3d(x+sel, y, z+sel), sel);
        children[6] = new Octree(depth+1,
                new Point3d(x, y+sel, z+sel), sel);
        children[7] = new Octree(depth+1,
                new Point3d(x+sel, y+sel, z+sel), sel);
    }
    
    private void addAtom(Atom atom) {
        if (children != null) {
            int index = findAtomIndex(atom);
            if (index != -1) {
                children[index].addAtom(atom);
                return;
            }
        }
        contents.add(atom);
        //logger.info(format("atom %d:    %d / %d    %d / %d", atom.getXYZIndex(), contents.size(), maxAtomsPerVolume, depth, maxTreeDepth));
        if (contents.size() > maxAtomsPerVolume && depth < maxTreeDepth) {
            if (children != null) {
                if (!warnedAboutSplit) {
                    logger.warning("Octree couldn't split: too many straddling atoms.");
                    warnedAboutSplit = true;
                }
                return;
            }
            split();
            for (int i = 0; i < contents.size(); ) {
                int index = findAtomIndex(contents.get(i));
                if (index != -1) {
                    children[index].addAtom(contents.remove(i));
                } else {
                    i++;
                }
            }
        }
    }
    
    private int findAtomIndex(Atom atom) {
        int index = -1;
        double midX = corner.getX() + (edgeLength / 2);
        double midY = corner.getY() + (edgeLength / 2);
        double midZ = corner.getZ() + (edgeLength / 2);
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
        if (leftSide && topSide && frontSide) {
            index = 0;
        } else if (rightSide && topSide && frontSide) {
            index = 1;
        } else if (leftSide && botSide && frontSide) {
            index = 2;
        } else if (leftSide && topSide && backSide) {
            index = 3;
        } else if (rightSide && botSide && frontSide) {
            index = 4;
        } else if (rightSide && topSide && backSide) {
            index = 5;
        } else if (leftSide && botSide && backSide) {
            index = 6;
        } else if (rightSide && botSide && backSide) {
            index = 7;
        }        
        return index;
    }
    
    public int getDepth() { return depth; }
    public Octree[] getChildren() { return children; }
    public List<Atom> getContents() { return contents; }
    public Point3d getCorner() { return corner; }
    public double getEdgeLength() { return edgeLength; }
    
}
