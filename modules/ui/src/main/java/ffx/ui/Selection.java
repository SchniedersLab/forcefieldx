/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
package ffx.ui;

import java.util.ArrayList;
import java.util.Enumeration;

import ffx.potential.bonded.MSNode;

/**
 * The Selection class will be used to make recursive multiscale selections,
 * however its implementation is not yet complete.
 *
 * @author Michael J. Schnieders
 *
 */
public class Selection {

    /**
     * <p>select</p>
     *
     * @param m a {@link ffx.potential.bonded.MSNode} object.
     * @param scale a {@link java.lang.Class} object.
     * @param criteria a {@link java.lang.String} object.
     * @return a {@link ffx.ui.Selection} object.
     */
    public static Selection select(MSNode m, Class scale, String criteria) {
        Selection s = new Selection(m, scale, criteria);
        s.evaluate();
        return s;
    }

    /**
     * <p>select</p>
     *
     * @param m a {@link ffx.potential.bonded.MSNode} object.
     * @param scale a {@link java.lang.String} object.
     * @param criteria a {@link java.lang.String} object.
     * @return a {@link ffx.ui.Selection} object.
     */
    public static Selection select(MSNode m, String scale, String criteria) {
        Selection s = new Selection(m, scale, criteria);
        s.evaluate();
        return s;
    }
    Class scale = null;
    String criteria = null;
    MSNode m = null;
    private ArrayList<MSNode> selected = null;

    private Selection(MSNode m, Class scale, String criteria) {
        this.scale = scale;
        this.criteria = criteria;
        this.m = m;
        selected = new ArrayList<MSNode>();
    }

    private Selection(MSNode m, String scale, String criteria) {
        try {
            this.scale = Class.forName(scale);
        } catch (Exception e) {
            this.scale = null;
        }
        this.criteria = criteria;
        this.m = m;
        selected = new ArrayList<MSNode>();
    }

    /**
     * <p>and</p>
     *
     * @param scale a {@link java.lang.Class} object.
     * @param criteria a {@link java.lang.String} object.
     * @return a {@link ffx.ui.Selection} object.
     */
    public Selection and(Class scale, String criteria) {
        this.scale = scale;
        this.criteria = criteria;
        evaluate();
        return this;
    }

    private void evaluate() {
        Enumeration e = m.depthFirstEnumeration();
        while (e.hasMoreElements()) {
            MSNode n = (MSNode) e.nextElement();
            if (scale.isInstance(n)) {
                selected.add(n);
            }
        }
    }

    /**
     * <p>or</p>
     *
     * @param scale a {@link java.lang.Class} object.
     * @param criteria a {@link java.lang.String} object.
     * @return a {@link ffx.ui.Selection} object.
     */
    public Selection or(Class scale, String criteria) {
        this.scale = scale;
        this.criteria = criteria;
        evaluate();
        return this;
    }
}
