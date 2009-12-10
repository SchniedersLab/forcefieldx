/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.potential.parameters;

import java.util.Vector;
import java.util.logging.Logger;

/**
 * The Keyword class holds a single Force Field X keyword entry.
 */
public class Keyword {
	private static final Logger logger = Logger.getLogger(Keyword.class.getName());

        private String keyword = null;
	private Vector<String> data = null;

	public Keyword(String k) {
		keyword = k;
		data = new Vector<String>();
	}

	public Keyword(String k, String entry) {
		this(k);
		data.add(entry);
	}

	public Keyword(String k, String entry[]) {
		this(k);
		for (String s : entry) {
			data.add(s);
		}
	}

	public void append(String entry) {
		data.add(entry);
	}

	public void append(String entry[]) {
		for (String s : entry) {
			data.add(s);
		}
	}

	public void clear() {
		data.clear();
	}

	public Vector<String> getEntries() {
		return data;
	}

	public String getEntry(int i) {
		return data.get(i);
	}

	public String getKeyword() {
		return keyword;
	}

	public void print() {
		logger.info(this.toString());
	}

        @Override
	public String toString() {
		StringBuffer sb = new StringBuffer(keyword + " ");
		for (String s : data) {
			sb.append(s);
		}
		return sb.toString();
	}
}
