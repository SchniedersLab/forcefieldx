/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
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

import java.util.ArrayList;
import java.util.Enumeration;

import ffx.potential.bonded.MSNode;


/**
 * The Selection class will be used to make recursive multiscale selections,
 * however its implementation is not yet complete.
 */
public class Selection {
	public static Selection select(MSNode m, Class scale, String criteria) {
		Selection s = new Selection(m, scale, criteria);
		s.evaluate();
		return s;
	}

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

	public Selection or(Class scale, String criteria) {
		this.scale = scale;
		this.criteria = criteria;
		evaluate();
		return this;
	}
}
