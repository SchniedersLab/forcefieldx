/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Biophysics Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
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
package ffx.ui;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Enumeration;

import javax.swing.JLabel;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.DefaultTreeSelectionModel;
import javax.swing.tree.RowMapper;
import javax.swing.tree.TreePath;

import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MSRoot;
import ffx.potential.bonded.ROLSP;
import ffx.potential.bonded.RendererCache;

/**
 * The Hierarchy Class creates and manages a JTree view of the data structure.
 * It is used for synchronization, handles the selection mechanism, and sets the
 * active system and nodes.
 */
public final class Hierarchy extends JTree implements TreeSelectionListener {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private MSRoot root;
    private MainPanel mainPanel;
    private DefaultTreeModel treeModel;
    private DefaultTreeSelectionModel treeSelectionModel;
    private FFXSystem activeSystem = null; // Reference to the active FSystem
    private MSNode activeNode = null; // Reference to the last Selected Node
    private ArrayList<MSNode> activeNodes = new ArrayList<MSNode>();
    private JLabel status = null;
    private JLabel step = null;
    private JLabel energy = null;
    ArrayList<TreePath> newPaths = new ArrayList<TreePath>();
    ArrayList<TreePath> removedPaths = new ArrayList<TreePath>();
    ArrayList<TreePath> prePaths = new ArrayList<TreePath>();
    ArrayList<MSNode> picks = new ArrayList<MSNode>();
    TreePath nullPath = null;

    /**
     * The Default Constructor initializes JTree properties, then updates is
     * representation based on the Structure of the Tree that extends from the
     * Root argument.
     */
    public Hierarchy(MainPanel f) {
        mainPanel = f;
        root = mainPanel.getDataRoot();
        initTree();
    }

    public void addSelection(MSNode f) {
        if (f == null) {
            return;
        }
        synchronized (this) {
            TreePath path = new TreePath(f.getPath());
            addSelectionPath(path);
            f.setSelected(true);
        }
    }

    public void addSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                addSelection(f);
            }
        }
    }

    public void addSystemNode(FFXSystem newSystem) {
        addTreeNode(newSystem, root, root.getChildCount());
    }

    public void addTreeNode(MSNode nodeToAdd, MSNode parent, int index) {
        if (nodeToAdd == null || nodeToAdd.getParent() != null) {
            return;
        }
        synchronized (this) {
            int childCount = parent.getChildCount();
            if (index < 0 || index > childCount) {
                index = parent.getChildCount();
            }
            // Add a parallel node if the ffe.lang.parallel flag was set
            if (ROLSP.GO_PARALLEL) {
                ROLSP parallelNode = new ROLSP();
                parallelNode.add(nodeToAdd);
                treeModel.insertNodeInto(parallelNode, parent, index);
            } else {
                treeModel.insertNodeInto(nodeToAdd, parent, index);
            }
            if (nodeToAdd instanceof FFXSystem) {
                attach((FFXSystem) nodeToAdd);
                treeModel.nodeStructureChanged(nodeToAdd);
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
        newModel.finalize(true);
        GraphicsCanvas graphics = mainPanel.getGraphics3D();
        if (graphics != null) {
            graphics.attachModel(newModel);
            if (newModel.getBondList().size() == 0) {
                mainPanel.getGraphics3D().updateScene(newModel, false, true,
                        RendererCache.ViewModel.SPACEFILL, false, null);
            }
        }
    }

    public void collapseAll() {
        int row = getRowCount() - 1;
        while (row >= 0) {
            collapseRow(row);
            row--;
        }
    }

    /**
     * Returns the active FSystem.
     */
    public FFXSystem getActive() {
        return activeSystem;
    }

    public MSNode getActiveNode() {
        return activeNode;
    }

    public ArrayList<MSNode> getActiveNodes() {
        return activeNodes;
    }

    public FFXSystem[] getNonActiveSystems() {
        synchronized (this) {
            int childCount = root.getChildCount();
            if (childCount == 0) {
                return null;
            }
            FFXSystem[] systems = new FFXSystem[childCount - 1];
            int index = 0;
            for (Enumeration e = root.children(); e.hasMoreElements();) {
                FFXSystem system = (FFXSystem) e.nextElement();
                if (system != getActive()) {
                    systems[index++] = system;
                }
            }
            return systems;
        }
    }

    public FFXSystem[] getSystems() {
        synchronized (this) {
            int childCount = root.getChildCount();
            if (childCount == 0) {
                return null;
            }
            FFXSystem[] systems = new FFXSystem[childCount];
            int index = 0;
            for (Enumeration e = root.children(); e.hasMoreElements();) {
                systems[index++] = (FFXSystem) e.nextElement();
            }
            return systems;
        }
    }

    public void groupSelection(MSNode f1, MSNode f2) {
        if (f1 == null || f2 == null) {
            return;
        }
        synchronized (this) {
            TreePath paths[] = new TreePath[2];
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
    public void initTree() {
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
        treeModel = new DefaultTreeModel(root);
        treeSelectionModel = new DefaultTreeSelectionModel();
        setModel(treeModel);
        setSelectionModel(treeSelectionModel);
        setRootVisible(false);
    }

    public void onlySelection(MSNode f) {
        synchronized (this) {
            if (activeNodes != null) {
                int num = activeNodes.size();
                TreePath paths[] = new TreePath[num];
                for (int i = 0; i < num; i++) {
                    paths[i] = new TreePath((activeNodes.get(i)).getPath());
                }
                removeSelectionPaths(paths);
            }
            collapseAll();
            addSelection(f);
        }
    }

    public void removeSelection(MSNode f) {
        if (f == null) {
            return;
        }
        synchronized (this) {
            TreePath path = new TreePath(f.getPath());
            for (Enumeration e = getExpandedDescendants(path); e.hasMoreElements();) {
                TreePath treePath = new TreePath((DefaultMutableTreeNode) e.nextElement());
                collapsePath(treePath);
            }
            removeSelectionPath(path);
            f.setSelected(false);
        }
    }

    public void removeSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                removeSelection(f);
            }
        }
    }

    public void removeTreeNode(MSNode nodeToRemove) {
        if (nodeToRemove == null) {
            return;
        }
        synchronized (this) {
            if (root.getChildCount() <= 1) {
                setRootVisible(false);
            }
            treeModel.removeNodeFromParent(nodeToRemove);
            if (getActive() == nodeToRemove && root.getChildCount() != 0) {
                FFXSystem m = (FFXSystem) root.getChildAt(0);
                setActive(m);
                onlySelection(activeSystem);
            } else {
                setActive(null);
            }
        }
    }

    public void selectAll() {
        if (activeSystem == null) {
            return;
        }
        synchronized (this) {
            onlySelection(root);
        }
    }

    /**
     * Sets the FSystem s to be active.
     */
    public void setActive(FFXSystem s) {
        if (s == activeSystem) {
            return;
        }
        synchronized (this) {
            activeSystem = s;
            updateStatus();
            if (mainPanel.getKeywordPanel() != null) {
                mainPanel.getKeywordPanel().loadActive(activeSystem);
            }
            /*
            if (mainPanel.getModelingPanel() != null) {
                mainPanel.getModelingPanel().loadActive(activeSystem);
            } */
            if (mainPanel.getModelingShell() != null) {
                mainPanel.getModelingShell().sync();
            }
        }
    }

    public void setActive(int i) {
        synchronized (this) {
            if (i < root.getChildCount()) {
                setActive((FFXSystem) root.getChildAt(i));
            } else if (root.getChildCount() == 0) {
                setActive(null);
            }
        }
    }

    public void setHighlighting(boolean h) {
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

    public void setStatus(JLabel s, JLabel t, JLabel e) {
        status = s;
        step = t;
        energy = e;
    }

    public void toggleSelection(MSNode f) {
        if (f == null) {
            return;
        }
        synchronized (this) {
            TreePath path = new TreePath(f.getPath());
            if (isPathSelected(path)) {
                removeSelectionPath(path);
            } else {
                addSelectionPath(path);
            }
        }
    }

    public void toggleSelections(ArrayList<MSNode> a) {
        synchronized (this) {
            for (MSNode f : a) {
                toggleSelection(f);
            }
        }
    }

    @Override
    public String toString() {
        return "Structural Hierarchy";
    }

    public void updateStatus() {
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
            ArrayList<TreePath> temp = prePaths;
            prePaths = newPaths;
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
            for (int i = 0; i < prePaths.size(); i++) {
                pathi = prePaths.get(i);
                for (int j = 0; j < removedPaths.size(); j++) {
                    pathj = removedPaths.get(j);
                    if (pathj.isDescendant(pathi)) {
                        prePaths.set(i, nullPath);
                        break;
                    }
                }
            }
            check = true;
            while (check) {
                check = prePaths.remove(nullPath);
            }
            // Combine new Paths and Existing Paths non-redundantly
            for (int i = 0; i < newPaths.size(); i++) {
                pathi = newPaths.get(i);
                if (pathi == nullPath) {
                    continue;
                }
                for (int j = 0; j < prePaths.size(); j++) {
                    pathj = prePaths.get(j);
                    if (pathj == nullPath) {
                        continue;
                    }
                    if (pathi == nullPath || pathj == nullPath) {
                        continue;
                    }
                    if (pathi.isDescendant(pathj)) {
                        prePaths.set(j, nullPath);
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
                check = prePaths.remove(nullPath);
            }
            newPaths.addAll(prePaths);
            activeNodes.clear();
            for (int i = 0; i < newPaths.size(); i++) {
                pathi = newPaths.get(i);
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
            picks = new ArrayList<MSNode>();
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
                for (TreePath p : prePaths) {
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
