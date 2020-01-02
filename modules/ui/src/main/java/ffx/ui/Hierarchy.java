//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.ui;

import javax.swing.JLabel;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.DefaultTreeSelectionModel;
import javax.swing.tree.RowMapper;
import javax.swing.tree.TreePath;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.logging.Logger;

import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MSRoot;
import ffx.potential.bonded.ROLSP;
import ffx.potential.bonded.RendererCache;

/**
 * The Hierarchy Class creates and manages a JTree view of the data structure.
 * It is used for synchronization, handles the selection mechanism, and sets the
 * active system and nodes.
 *
 * @author Michael J. Schnieders
 */
public final class Hierarchy extends JTree implements TreeSelectionListener {

    private static final Logger logger = Logger.getLogger(Hierarchy.class.getName());
    private final MSRoot root;
    private final MainPanel mainPanel;
    private DefaultTreeModel hierarchyModel;
    private DefaultTreeSelectionModel treeSelectionModel;
    // Reference to the active FFXSystem
    private FFXSystem activeSystem = null;
    // Reference to the last Selected Node
    private MSNode activeNode = null;
    private final ArrayList<MSNode> activeNodes = new ArrayList<>();
    private JLabel status = null;
    private JLabel step = null;
    private JLabel energy = null;
    private ArrayList<TreePath> newPaths = new ArrayList<>();
    private ArrayList<TreePath> removedPaths = new ArrayList<>();
    private ArrayList<TreePath> previousPaths = new ArrayList<>();
    private TreePath nullPath = null;

    /**
     * The Default Constructor initializes JTree properties, then updates is
     * representation based on the Structure of the Tree that extends from the
     * Root argument.
     *
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    Hierarchy(MainPanel f) {
        mainPanel = f;
        root = mainPanel.getDataRoot();
        initTree();
    }

    /**
     * <p>
     * addSelection</p>
     *
     * @param f a {@link ffx.potential.bonded.MSNode} object.
     */
    private void addSelection(MSNode f) {
        if (f == null) {
            return;
        }
        synchronized (this) {
            TreePath path = new TreePath(f.getPath());
            addSelectionPath(path);
            f.setSelected(true);
        }
    }

    /**
     * <p>
     * addSelections</p>
     *
     * @param a a {@link java.util.ArrayList} object.
     */
    public void addSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                addSelection(f);
            }
        }
    }

    /**
     * <p>
     * addSystemNode</p>
     *
     * @param newSystem a {@link ffx.ui.FFXSystem} object.
     */
    void addSystemNode(FFXSystem newSystem) {
        addTreeNode(newSystem, root, root.getChildCount());
    }

    /**
     * <p>
     * addTreeNode</p>
     *
     * @param nodeToAdd a {@link ffx.potential.bonded.MSNode} object.
     * @param parent    a {@link ffx.potential.bonded.MSNode} object.
     * @param index     a int.
     */
    private void addTreeNode(MSNode nodeToAdd, MSNode parent, int index) {
        synchronized (this) {
            if (nodeToAdd == null || nodeToAdd.getParent() != null) {
                return;
            }
            int childCount = parent.getChildCount();
            if (index < 0 || index > childCount) {
                index = parent.getChildCount();
            }

            String name = nodeToAdd.getName();

            for (int i = 0; i < parent.getChildCount(); i++) {
                MSNode node = (MSNode) parent.getChildAt(i);
                if (node.getName().equals(name)) {
                    logger.info(" Parent already has a node with the name " + name);
                }
            }

            // Add a parallel node if the ffe.lang.parallel flag was set
            if (ROLSP.GO_PARALLEL) {
                ROLSP parallelNode = new ROLSP();
                parallelNode.add(nodeToAdd);
                hierarchyModel.insertNodeInto(parallelNode, parent, index);
            } else {
                hierarchyModel.insertNodeInto(nodeToAdd, parent, index);
            }
            if (nodeToAdd instanceof FFXSystem) {
                attach((FFXSystem) nodeToAdd);
                hierarchyModel.nodeStructureChanged(nodeToAdd);
            }
            onlySelection(nodeToAdd);
            if (!isRootVisible()) {
                setRootVisible(true);
            }
        }
    }

    private void attach(FFXSystem newModel) {
        if (newModel == null) {
            return;
        }
        newModel.finalize(true, newModel.getForceField());
        GraphicsCanvas graphics = mainPanel.getGraphics3D();
        if (graphics != null) {
            graphics.attachModel(newModel);
            if (newModel.getBondList().isEmpty()) {
                mainPanel.getGraphics3D().updateScene(newModel, false, true,
                        RendererCache.ViewModel.SPACEFILL, false, null);
            }
        }
    }

    /**
     * <p>
     * collapseAll</p>
     */
    private void collapseAll() {
        int row = getRowCount() - 1;
        while (row >= 0) {
            collapseRow(row);
            row--;
        }
    }

    /**
     * Returns the active FSystem.
     *
     * @return a {@link ffx.ui.FFXSystem} object.
     */
    public FFXSystem getActive() {
        return activeSystem;
    }

    /**
     * <p>
     * Getter for the field <code>activeNode</code>.</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getActiveNode() {
        return activeNode;
    }

    /**
     * <p>
     * Getter for the field <code>activeNodes</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    ArrayList<MSNode> getActiveNodes() {
        return activeNodes;
    }

    /**
     * <p>
     * getNonActiveSystems</p>
     *
     * @return an array of {@link ffx.ui.FFXSystem} objects.
     */
    FFXSystem[] getNonActiveSystems() {
        synchronized (this) {
            int childCount = root.getChildCount();
            if (childCount == 0) {
                return null;
            }
            FFXSystem[] systems = new FFXSystem[childCount - 1];
            int index = 0;
            for (Enumeration e = root.children(); e.hasMoreElements(); ) {
                FFXSystem system = (FFXSystem) e.nextElement();
                if (system != getActive()) {
                    systems[index++] = system;
                }
            }
            return systems;
        }
    }

    /**
     * <p>
     * getSystems</p>
     *
     * @return an array of {@link ffx.ui.FFXSystem} objects.
     */
    public FFXSystem[] getSystems() {
        synchronized (this) {
            int childCount = root.getChildCount();
            if (childCount == 0) {
                return null;
            }
            FFXSystem[] systems = new FFXSystem[childCount];
            int index = 0;
            for (Enumeration e = root.children(); e.hasMoreElements(); ) {
                systems[index++] = (FFXSystem) e.nextElement();
            }
            return systems;
        }
    }

    /**
     * <p>
     * groupSelection</p>
     *
     * @param f1 a {@link ffx.potential.bonded.MSNode} object.
     * @param f2 a {@link ffx.potential.bonded.MSNode} object.
     */
    public void groupSelection(MSNode f1, MSNode f2) {
        if (f1 == null || f2 == null) {
            return;
        }
        synchronized (this) {
            TreePath[] paths = new TreePath[2];
            paths[0] = new TreePath(f1.getPath());
            paths[1] = new TreePath(f2.getPath());
            RowMapper rm = treeSelectionModel.getRowMapper();
            int[] rows = rm.getRowsForPaths(paths);
            setSelectionInterval(rows[0], rows[1]);
        }
    }

    /**
     * Initialize the Tree representation based on the Root data node
     */
    private void initTree() {
        addTreeSelectionListener(this);
        setExpandsSelectedPaths(true);
        setScrollsOnExpand(true);
        //setLargeModel(true);
        setEditable(false);
        putClientProperty("JTree.lineStyle", "Angled");
        setShowsRootHandles(true);
        DefaultTreeCellRenderer tcr = new DefaultTreeCellRenderer();
        tcr.setBackgroundSelectionColor(Color.yellow);
        tcr.setBorderSelectionColor(Color.black);
        tcr.setTextSelectionColor(Color.black);
        setCellRenderer(tcr);
        hierarchyModel = new DefaultTreeModel(root);
        treeSelectionModel = new DefaultTreeSelectionModel();
        setModel(hierarchyModel);
        setSelectionModel(treeSelectionModel);
        setRootVisible(false);
    }

    /**
     * <p>
     * onlySelection</p>
     *
     * @param f a {@link ffx.potential.bonded.MSNode} object.
     */
    void onlySelection(MSNode f) {
        synchronized (this) {
            int num = activeNodes.size();
            TreePath[] paths = new TreePath[num];
            for (int i = 0; i < num; i++) {
                paths[i] = new TreePath((activeNodes.get(i)).getPath());
            }
            removeSelectionPaths(paths);
            collapseAll();
            addSelection(f);
        }
    }

    /**
     * <p>
     * removeSelection</p>
     *
     * @param f a {@link ffx.potential.bonded.MSNode} object.
     */
    private void removeSelection(MSNode f) {
        synchronized (this) {
            if (f == null) {
                return;
            }
            TreePath path = new TreePath(f.getPath());
            for (Enumeration e = getExpandedDescendants(path); e.hasMoreElements(); ) {
                TreePath treePath = new TreePath(e.nextElement());
                collapsePath(treePath);
            }
            removeSelectionPath(path);
            f.setSelected(false);
        }
    }

    /**
     * <p>
     * removeSelections</p>
     *
     * @param a a {@link java.util.ArrayList} object.
     */
    public void removeSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                removeSelection(f);
            }
        }
    }

    /**
     * <p>
     * removeTreeNode</p>
     *
     * @param nodeToRemove a {@link ffx.potential.bonded.MSNode} object.
     */
    void removeTreeNode(MSNode nodeToRemove) {
        synchronized (this) {
            if (nodeToRemove == null) {
                return;
            }

            hierarchyModel.removeNodeFromParent(nodeToRemove);

            /*
              The DefaultTreeModel and DefaultTreeSelectionModel classes retain
              references to removed nodes. To work around this, we create new
              instances of these classes whenever an FFXSystem is removed.
             */
            if (nodeToRemove instanceof FFXSystem) {
                hierarchyModel = new DefaultTreeModel(root);
                treeSelectionModel = new DefaultTreeSelectionModel();
                setModel(hierarchyModel);
                setSelectionModel(treeSelectionModel);
            }

            // Whenever a node is removed, clear and reset activeNodes and path instances.
            activeNodes.clear();
            previousPaths.clear();
            newPaths.clear();
            removedPaths.clear();

            if (getActive() == nodeToRemove && root.getChildCount() != 0) {
                FFXSystem m = (FFXSystem) root.getChildAt(0);
                setActive(m);
                onlySelection(activeSystem);
            } else {
                setActive(null);
            }

            if (root.getChildCount() <= 1) {
                setRootVisible(false);
            }
        }
    }

    /**
     * <p>
     * selectAll</p>
     */
    void selectAll() {
        synchronized (this) {
            if (activeSystem == null) {
                return;
            }
            onlySelection(root);
        }
    }

    /**
     * Sets the FFXSystem s to be active.
     *
     * @param ffxSystem a {@link ffx.ui.FFXSystem} object.
     */
    public void setActive(FFXSystem ffxSystem) {
        synchronized (this) {
            if (ffxSystem == activeSystem) {
                return;
            }
            if (ffxSystem != null) {
                // A closed (closing) system cannot be set active.
                if (ffxSystem.isClosing()) {
                    // Set the most recently opened system to be active.
                    for (int i = root.getChildCount() - 1; i >= 0; i--) {
                        FFXSystem child = (FFXSystem) root.getChildAt(i);
                        if (!child.isClosing()) {
                            setActive(child);
                            return;
                        }
                    }
                    setActive(null);
                    return;
                }
            }
            activeSystem = ffxSystem;
            updateStatus();
            if (mainPanel.getKeywordPanel() != null) {
                mainPanel.getKeywordPanel().loadActive(activeSystem);
            }
            if (mainPanel.getModelingPanel() != null) {
                mainPanel.getModelingPanel().loadActive(activeSystem);
            }
            if (mainPanel.getModelingShell() != null) {
                mainPanel.getModelingShell().sync();
            }
        }
    }

    /**
     * <p>
     * setActive</p>
     *
     * @param i a int.
     */
    public void setActive(int i) {
        synchronized (this) {
            if (i < root.getChildCount()) {
                setActive((FFXSystem) root.getChildAt(i));
            } else if (root.getChildCount() == 0) {
                setActive(null);
            }
        }
    }

    /**
     * <p>
     * setHighlighting</p>
     *
     * @param h a boolean.
     */
    void setHighlighting(boolean h) {
        synchronized (this) {
            if (RendererCache.highlightSelections != h) {
                RendererCache.highlightSelections = h;
                for (MSNode node : activeNodes) {
                    node.setSelected(h);
                }
                mainPanel.getGraphics3D().updateScene(activeNodes, false,
                        false, null, true, RendererCache.ColorModel.SELECT);
            }
        }
    }

    /**
     * <p>
     * Setter for the field <code>status</code>.</p>
     *
     * @param s a {@link javax.swing.JLabel} object.
     * @param t a {@link javax.swing.JLabel} object.
     * @param e a {@link javax.swing.JLabel} object.
     */
    void setStatus(JLabel s, JLabel t, JLabel e) {
        status = s;
        step = t;
        energy = e;
    }

    /**
     * <p>
     * toggleSelection</p>
     *
     * @param f a {@link ffx.potential.bonded.MSNode} object.
     */
    void toggleSelection(MSNode f) {
        synchronized (this) {
            if (f == null) {
                return;
            }
            TreePath path = new TreePath(f.getPath());
            if (isPathSelected(path)) {
                removeSelectionPath(path);
            } else {
                addSelectionPath(path);
            }
        }
    }

    /**
     * <p>
     * toggleSelections</p>
     *
     * @param a a {@link java.util.ArrayList} object.
     */
    public void toggleSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                toggleSelection(f);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return "Structural Hierarchy";
    }

    /**
     * <p>
     * updateStatus</p>
     */
    void updateStatus() {
        if (activeSystem == null) {
            status.setText("  ");
            step.setText("  ");
            energy.setText("  ");
            return;
        }
        if (activeSystem.getFile() != null) {
            status.setText("  " + activeSystem.toFFString());
        } else {
            status.setText("  " + activeSystem.toString());
        }

        if (activeSystem.getCycles() > 1) {
            step.setText("" + activeSystem.getCurrentCycle() + "/" + activeSystem.getCycles());
        } else {
            step.setText("");
        }

        if (activeSystem.getCycles() > 1) {
            energy.setText("");
        } else {
            energy.setText("");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void valueChanged(TreeSelectionEvent e) {
        synchronized (this) {

            // Determine the Active System
            MSNode lastNode = (MSNode) getLastSelectedPathComponent();
            if (lastNode != null) {
                activeNode = lastNode;
                FFXSystem s = (FFXSystem) activeNode.getMSNode(FFXSystem.class);
                if (s != null) {
                    setActive(s);
                }
            }
            TreePath[] paths = e.getPaths();
            if (paths == null) {
                return;
            }

            // Reuse the same ArrayLists
            ArrayList<TreePath> temp = previousPaths;
            previousPaths = newPaths;
            newPaths = temp;
            // Determine new and removed paths
            newPaths.clear();
            removedPaths.clear();
            for (int i = 0; i < paths.length; i++) {
                if (e.isAddedPath(i)) {
                    newPaths.add(paths[i]);
                } else {
                    removedPaths.add(paths[i]);
                }
            }
            // Create a non-redundant set of new/removed paths
            TreePath pathi, pathj;
            for (int i = 0; i < newPaths.size(); i++) {
                pathi = newPaths.get(i);
                if (pathi == nullPath) {
                    continue;
                }
                for (int j = i + 1; j < newPaths.size(); j++) {
                    pathj = newPaths.get(j);
                    if (pathi == nullPath || pathj == nullPath) {
                        continue;
                    }
                    if (pathi.isDescendant(pathj)) {
                        newPaths.set(j, nullPath);
                    } else if (pathj.isDescendant(pathi)) {
                        newPaths.set(i, nullPath);
                    }
                }
            }
            boolean check = true;
            while (check) {
                check = newPaths.remove(nullPath);
            }
            for (int i = 0; i < removedPaths.size(); i++) {
                pathi = removedPaths.get(i);
                if (pathi == nullPath) {
                    continue;
                }
                for (int j = i + 1; j < removedPaths.size(); j++) {
                    pathj = removedPaths.get(j);
                    if (pathi == nullPath || pathj == nullPath) {
                        continue;
                    }
                    if (pathi.isDescendant(pathj)) {
                        removedPaths.set(j, nullPath);
                    } else if (pathj.isDescendant(pathi)) {
                        removedPaths.set(i, nullPath);
                    }
                }
            }
            check = true;
            while (check) {
                check = removedPaths.remove(nullPath);
            }
            // Remove the RemovedPaths from the Existing List
            for (int i = 0; i < previousPaths.size(); i++) {
                pathi = previousPaths.get(i);
                for (TreePath removedPath : removedPaths) {
                    if (removedPath.isDescendant(pathi)) {
                        previousPaths.set(i, nullPath);
                        break;
                    }
                }
            }
            check = true;
            while (check) {
                check = previousPaths.remove(nullPath);
            }
            // Combine new Paths and Existing Paths non-redundantly
            for (int i = 0; i < newPaths.size(); i++) {
                pathi = newPaths.get(i);
                if (pathi == nullPath) {
                    continue;
                }
                for (int j = 0; j < previousPaths.size(); j++) {
                    pathj = previousPaths.get(j);
                    if (pathj == nullPath) {
                        continue;
                    }
                    if (pathi == nullPath) {
                        continue;
                    }
                    if (pathi.isDescendant(pathj)) {
                        previousPaths.set(j, nullPath);
                    } else if (pathj.isDescendant(pathi)) {
                        newPaths.set(i, nullPath);
                    }
                }
            }
            check = true;
            while (check) {
                check = newPaths.remove(nullPath);
            }
            check = true;
            while (check) {
                check = previousPaths.remove(nullPath);
            }
            newPaths.addAll(previousPaths);
            activeNodes.clear();
            for (TreePath newPath : newPaths) {
                pathi = newPath;
                activeNodes.add((MSNode) pathi.getLastPathComponent());
            }
            if (activeNode != null) {
                TreePath activePath = new TreePath(activeNode);
                expandPath(activePath.getParentPath());
                makeVisible(activePath);
                scrollPathToVisible(activePath);
            }
            // We now have a non-redundant set of Active Paths; and a
            // non-redundant set of removed paths
            ArrayList<MSNode> picks = new ArrayList<>();
            // Clear highlight of de-selected nodes
            for (TreePath r : removedPaths) {
                boolean change = true;
                for (TreePath n : newPaths) {
                    if (n.isDescendant(r)) {
                        change = false;
                    }
                }
                if (change) {
                    MSNode f = (MSNode) r.getLastPathComponent();
                    f.setSelected(false);
                    picks.add(f);
                }
            }
            for (TreePath n : newPaths) {
                boolean change = true;
                for (TreePath p : previousPaths) {
                    if (p.isDescendant(n)) {
                        change = false;
                    }
                }
                if (change) {
                    MSNode f = (MSNode) n.getLastPathComponent();
                    f.setSelected(true);
                    picks.add(f);
                }
            }
            if (RendererCache.highlightSelections) {
                mainPanel.getGraphics3D().updateScene(picks, false, false,
                        null, true, RendererCache.ColorModel.SELECT);
            } else if (RendererCache.labelAtoms || RendererCache.labelResidues) {
                mainPanel.getGraphics3D().setLabelsUpdated();
            }
        }
    }
}
