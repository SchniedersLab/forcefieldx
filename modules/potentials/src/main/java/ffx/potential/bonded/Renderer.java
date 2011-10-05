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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.logging.Logger;

import javax.media.j3d.Behavior;
import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.WakeupOnBehaviorPost;
import javax.swing.JLabel;

/**
 * The Renderer class attempts to maximize throughput of graphics operations on
 * MolecularAssembly instances.
 *
 * @author schnied
 * @version $Id: $
 */
public class Renderer extends Behavior {

    private static final Logger logger = Logger.getLogger(Renderer.class.getName());
    private static long frameNumber = 0;
    private static long frameDuration;
    private ArrayList<MSNode> nodesToUpdate = null;
    private ArrayList<MSNode> nodesCache = null;
    private boolean doTransform, doView, doColor;
    private boolean doTransformCache, doViewCache, doColorCache;
    private JLabel statusBar = null;
    private RendererCache.ViewModel viewModel, viewModelCache;
    private RendererCache.ColorModel colorModel, colorModelCache;
    private WakeupOnBehaviorPost postid;
    private boolean timer = false;
    private boolean gc = false;

    /**
     * Constructor
     *
     * @param b
     *            Bounds of this behavior
     * @param s
     *            JLabel for status
     */
    public Renderer(Bounds b, JLabel s) {
        setSchedulingBounds(b);
        statusBar = s;
        if (System.getProperty("ffe.timer", "false").equalsIgnoreCase("true")) {
            timer = true;
            if (System.getProperty("ffe.timer.gc", "false").equalsIgnoreCase(
                    "true")) {
                gc = true;
            }
        }
    }

    /**
     * This node arms UpdateBehavior with a graphics operation to carry out
     *
     * @param nodes
     *            Nodes where the operation will be performed
     * @param t
     *            True for a change in atomic position
     * @param v
     *            True for a change in rendering method
     * @param vtype
     *            The rendering method to use
     * @param c
     *            True for a change in rendering color
     * @param ctype
     *            The coloring method to use
     */
    public void arm(ArrayList<MSNode> nodes, boolean t, boolean v,
            RendererCache.ViewModel vtype, boolean c,
            RendererCache.ColorModel ctype) {
        // If the node isn't null, the last rendering
        // operation hasn't finished so one operation will be cached
        if (nodesToUpdate != null) {
            nodesCache = nodes;
            doTransformCache = t;
            doViewCache = v;
            viewModelCache = vtype;
            doColorCache = c;
            colorModelCache = ctype;
        } else {
            nodesToUpdate = nodes;
            doTransform = t;
            doView = v;
            viewModel = vtype;
            doColor = c;
            colorModel = ctype;
            postId(1);
        }
    }

    /**
     * <p>arm</p>
     *
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     * @param t a boolean.
     * @param v a boolean.
     * @param vtype a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     * @param c a boolean.
     * @param ctype a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void arm(MSNode node, boolean t, boolean v,
            RendererCache.ViewModel vtype, boolean c,
            RendererCache.ColorModel ctype) {
        ArrayList<MSNode> temp = new ArrayList<MSNode>();
        temp.add(node);
        arm(temp, t, v, vtype, c, ctype);
    }

    /**
     * {@inheritDoc}
     *
     * Initialize this behavior to respond to postID messages
     */
    @Override
    public void initialize() {
        postid = new WakeupOnBehaviorPost(this, 1);
        wakeupOn(postid);
    }

    /**
     * Check to see if a graphics operation is pending/executing
     *
     * @return Whether a node has been cued
     */
    public boolean isArmed() {
        if (nodesToUpdate != null) {
            return true;
        }
        return false;
    }

    /**
     * <p>isCacheFull</p>
     *
     * @return a boolean.
     */
    public boolean isCacheFull() {
        if (nodesCache != null) {
            return true;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     *
     * This method is called by the Java3D Behavior thread after the following
     * sequence of events: 1.) A graphics operation is loaded using the "arm"
     * method. 2.) The PostID call is processed by the Java3D Behavior Thread.
     */
    @Override
    public void processStimulus(Enumeration parm1) {
        // Do not perform two operations before the frame has be refreshed.
        if (getView().getFrameNumber() == frameNumber) {
            System.out.print(".");
            wakeupOn(postid);
            postId(1);
            return;
        }
        // Check that the requested View and Color Models are known.
        String viewString = null;
        String colorString = null;
        if (viewModel != null) {
            try {
                viewString = viewModel.toString();
            } catch (Exception e) {
                statusBar.setText("Unknown ViewModel: " + viewModel);
                return;
            }
        }
        if (colorModel != null) {
            try {
                colorString = colorModel.toString();
            } catch (Exception e) {
                statusBar.setText("Unknown ColorModel: " + colorModel);
                return;
            }
        }
        if (timer) {
            startTimer();
            if (viewString != null) {
                logger.info("Applying ViewModel Change: " + viewString);
            }
            if (colorString != null) {
                System.out.println("Applying ColorModel Change: " + colorString);
            }
        }
        // Perform the requested rendering operation
        ArrayList<ArrayList<BranchGroup>> newChildren = new ArrayList<ArrayList<BranchGroup>>();
        for (MSNode nodeToUpdate : nodesToUpdate) {
            if (nodeToUpdate == null) {
                continue;
            }
            if (doTransform) {
                nodeToUpdate.update();
            }
            if (doColor) {
                nodeToUpdate.setColor(colorModel, null, null);
                if (statusBar != null) {
                    statusBar.setText("  Color by \"" + colorString + "\" applied to " + nodeToUpdate.toString());
                }
            }
            if (doView) {
                ArrayList<BranchGroup> newShapes = new ArrayList<BranchGroup>();
                newChildren.add(newShapes);
                nodeToUpdate.setView(viewModel, newShapes);
                if (statusBar != null) {
                    statusBar.setText("  Style \"" + viewString + "\" applied to " + nodeToUpdate.toString());
                }
            }
        }
        // Wait for the parallel nodes to finish
        try {
            if (ROLSP.GO_PARALLEL && ROLSP.parallelNotDone > 0) {
                logger.info("Renderer waiting for " + ROLSP.parallelNotDone + " processes...");
            }
            while (ROLSP.GO_PARALLEL && ROLSP.parallelNotDone > 0) {
                synchronized (this) {
                    wait(10);
                }
            }
        } catch (Exception e) {
            System.out.println("Exception Waiting for Parallel MultiScale Methods to Finish");
        } finally {
            // If there are new children, they can not be added in parallel
            // because Java3D does not seem to be thread safe.
            // (There are ArrayList that are not synchronized).
            // Here we will add them one at a time. The cases are setView being
            // called on nodes below the
            // Scenegraph attachment points (MolecularAssemblies), setView being
            // called on the root node, setView
            // being called on a ParallelMSM node, or setView being called on
            // the MolecularAssembly.
            for (int i = 0; i < nodesToUpdate.size(); i++) {
                if (newChildren.isEmpty()) {
                    break;
                }
                MSNode nodeToUpdate = nodesToUpdate.get(i);
                if (nodeToUpdate == null) {
                    continue;
                }
                if (nodeToUpdate instanceof MolecularAssembly) {
                    MolecularAssembly ma = (MolecularAssembly) nodeToUpdate;
                    ma.sceneGraphChange(null);
                } else if (nodeToUpdate instanceof ROLSP) {
                    MolecularAssembly ma = (MolecularAssembly) nodeToUpdate.getChildAt(0);
                    ma.sceneGraphChange(null);
                } else if (nodeToUpdate instanceof MSRoot) {
                    for (Enumeration e = nodeToUpdate.children(); e.hasMoreElements();) {
                        MSNode updatedNode = (MSNode) e.nextElement();
                        if (updatedNode instanceof ROLSP) {
                            MolecularAssembly ma = (MolecularAssembly) updatedNode.getChildAt(0);
                            ma.sceneGraphChange(null);
                        } else {
                            MolecularAssembly ma = (MolecularAssembly) updatedNode;
                            ma.sceneGraphChange(null);
                        }
                    }
                } else {
                    ArrayList<BranchGroup> newShapes = newChildren.get(i);
                    if (!newShapes.isEmpty()) {
                        MolecularAssembly ma = (MolecularAssembly) nodeToUpdate.getMSNode(MolecularAssembly.class);
                        ma.sceneGraphChange(newShapes);
                    }
                }
            }
        }
        if (timer) {
            stopTimer();
        }
        nodesToUpdate = null;
        wakeupOn(postid);
        if (nodesCache != null) {
            nodesToUpdate = nodesCache;
            doTransform = doTransformCache;
            doView = doViewCache;
            viewModel = viewModelCache;
            doColor = doColorCache;
            colorModel = colorModelCache;
            nodesCache = null;
            postId(1);
        }
    }

    private void startTimer() {
        Runtime runtime = Runtime.getRuntime();
        frameDuration = getView().getLastFrameDuration();
        if (gc) {
            System.out.print("Running Finalization and GC for acccurate memory usage...");
            runtime.runFinalization();
            runtime.gc();
            logger.info(" Done\nProceeding with graphics operation...");
        }
    }

    private void stopTimer() {
        Runtime runtime = Runtime.getRuntime();
        frameNumber = getView().getFrameNumber();
        frameDuration = getView().getLastFrameDuration();
        logger.info("Frame Duration After Op: " + frameDuration / 1000);
        if (gc) {
            runtime.runFinalization();
            runtime.gc();
            logger.info(" Done");
        }
    }
}
