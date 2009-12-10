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
package ffx.ui.commands;

import java.io.Serializable;

/**
 * The FFXMessage class is used to pass simple messages between a TinkerServer
 * and its FFXClient(s).
 */
public class FFXMessage implements Serializable {
	private static final long serialVersionUID = 1L;
	public static int SYSTEM = 0;
	public static int UPDATE = 1;
	public static int CLOSING = 2;
	public static int OK = 3;
	private int message = 0;
	private int step = -1;
	private double time = -1.0;
	private int type = 0;

	public FFXMessage(int m) {
		message = m;
	}

	public int getMessage() {
		return message;
	}

	public int getStep() {
		return step;
	}

	public double getTime() {
		return time;
	}

	public int getType() {
		return type;
	}

	public void print() {
		System.out.println(toString());
	}

	public void setMessage(int m) {
		message = m;
	}

	public void setStep(int s) {
		step = s;
	}

	public void setTime(double t) {
		time = t;
	}

	public void setType(int t) {
		type = t;
	}

	public String toString() {
		if (message == 0) {
			return new String("SYSTEM").intern();
		} else if (message == 1) {
			return new String("UPDATE").intern();
		} else if (message == 2) {
			return new String("CLOSING").intern();
		} else {
			return new String("OK").intern();
		}
	}
}
