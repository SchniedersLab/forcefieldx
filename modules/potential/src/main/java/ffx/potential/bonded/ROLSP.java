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
package ffx.potential.bonded;

import javax.media.j3d.BranchGroup;

import java.util.Enumeration;
import java.util.List;
import java.util.logging.Logger;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The ROLSP class is used for Proof-Of-Concept Parallel Recusive Over Length
 * Scales (ROLS) Methods (currently only on shared memory systems). Simply
 * Simply inserting a ParallelMSM node into the Hierarchy causes a seperate
 * thread of execution to be created for all operations on nodes below the ROLSP
 * node. This is very preliminary code, but a useful concept for parallelizing
 * ROLS in ffe.lang.
 *
 * @author Michael J. Schnieders
 *
 */
public class ROLSP extends MSNode implements ROLS, Runnable {

    private static final Logger logger = Logger.getLogger(ROLSP.class.getName());

    public enum PARALLELMETHOD {

        SETVIEW, NONE;
    }
    private static final long serialVersionUID = 1L;
    /**
     * Constant <code>GO_PARALLEL=false</code>
     */
    public static boolean GO_PARALLEL = false;
    /**
     * Constant <code>parallelNotDone=0</code>
     */
    public static int parallelNotDone = 0;

    static {
        try {
            boolean b = Boolean.parseBoolean(System.getProperty(
                    "ffx.lang.parallel", "false"));
            GO_PARALLEL = b;
        } catch (Exception e) {
            GO_PARALLEL = false;
        }
    }
    private PARALLELMETHOD parallelMethod = PARALLELMETHOD.NONE;
    private Thread thread = null;
    private long startTime = 0;
    private long threadTime = 0;
    private RendererCache.ViewModel viewModel = null;
    private List<BranchGroup> newShapes = null;

    /**
     * <p>
     * Constructor for ROLSP.</p>
     */
    public ROLSP() {
        super("Parallel Node");
    }

    /**
     * {@inheritDoc}
     *
     * Overidden equals method.
     */
    @Override
    public boolean equals(Object object) {
        if (!(object instanceof ROLSP)) {
            return false;
        }
        if (this == object) {
            return true;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        MSNode child = (MSNode) getChildAt(0);
        if (child == null) {
            return hash(SEED, "none".hashCode());
        }
        return hash(SEED, child.hashCode());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        switch (parallelMethod) {
            case SETVIEW:
                setView(viewModel, newShapes);
                break;
            default:
                return;
        }
        threadTime = System.currentTimeMillis() - startTime;
        logger.info("Start Time: " + startTime + " Total Time: " + threadTime);
        parallelNotDone--;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel viewModel,
            List<BranchGroup> newShapes) {
        // Set Up the Parallel setView Method
        if (parallelMethod == PARALLELMETHOD.NONE) {
            startTime = System.currentTimeMillis();
            this.viewModel = viewModel;
            this.newShapes = newShapes;
            parallelMethod = PARALLELMETHOD.SETVIEW;
            thread = new Thread(this);
            thread.setName(getParent().toString() + ": Parallel setView MSM");
            thread.setPriority(Thread.MAX_PRIORITY);
            parallelNotDone++;
            thread.start();
        } else if (parallelMethod == PARALLELMETHOD.SETVIEW) {
            // setView has been called from within the 'run' method of the
            // "setView" thread
            for (Enumeration e = children(); e.hasMoreElements();) {
                MSNode node = (MSNode) e.nextElement();
                node.setView(viewModel, newShapes);
            }
            parallelMethod = PARALLELMETHOD.NONE;
            thread = null;
            viewModel = null;
            newShapes = null;
        } else {
            logger.info("Parallel setView method called by: " + parallelMethod);
            return;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (threadTime != 0) {
            return "Parallel Node " + threadTime + " (msec)";
        }
        return "Parallel Node";
    }
}
