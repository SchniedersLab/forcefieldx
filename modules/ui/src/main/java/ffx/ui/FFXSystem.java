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

import java.io.File;
import java.util.Hashtable;

import ffx.ui.commands.TinkerUpdate;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.utilities.Keyword;

/**
 * The FFXSystem class contains extensions to the generic
 * ffe.lang.MolecularAssembly class specific to Force Field X interacting
 * with TINKER.
 */
public class FFXSystem extends MolecularAssembly {
	private static final long serialVersionUID = 50L;
	public static final int MultiScaleLevel = 4;
	// Log file being used for modeling commands
	private File logFile;
	// Key file for this system
	private File keyFile;
	private Hashtable<String, Keyword> keywords = new Hashtable<String, Keyword>();
	// Command Description if this System is the result of a TINKER commad
	private String commandDescription = null;
	// Archive
	private Trajectory trajectory = null;
	// Simulation type
	private int simulation = 0;
	// Simulation data
	private double time, temperature, energy;
	private int step;
	// Flag to indicate this System is being closed
	private boolean closing = false;

	/**
	 * FFXSystem Constructor
	 * 
	 * @param name
	 *            String
	 */
	public FFXSystem(String name, String description, File file) {
		super(name);
		setFile(file);
		commandDescription = description;
	}

	public void addKeyword(Keyword k) {
		if (keywords.containsKey(k.getKeyword())) {
			return;
		}
		keywords.put(k.getKeyword(), k);
	}

	public boolean destroy() {
		setClosing(true);
		return super.destroy();
	}

	public double getEnergy() {
		return energy;
	}

	public String getEnergyString() {
		return String.format("Energy: %9.3f kcal/mole", energy);
	}

	public File getKeyFile() {
		return keyFile;
	}

	public Keyword getKeyword(String k) {
		return keywords.get(k);
	}

	public Hashtable<String, Keyword> getKeywords() {
		return keywords;
	}

	public File getLogFile() {
		if (logFile == null) {
			if (getFile() == null) {
				return null;
			}
			String fileName = getFile().getName();
			int dot = fileName.lastIndexOf(".");
			fileName = fileName.subSequence(0, dot) + ".log";
			logFile = new File(fileName);
		}
		return logFile;
	}

	public String getStepString() {
		return String.format("Step: %12d", step);
	}

	public double getTemperature() {
		return temperature;
	}

	public double getTime() {
		return time;
	}

	public String getTimeString() {
		return String.format("Time: %9.3f picoseconds", this.time);
	}

	public Trajectory getTrajectory() {
		return trajectory;
	}

	public boolean isClosing() {
		return closing;
	}

	public boolean isOptimization() {
		if (simulation == TinkerUpdate.OPTIMIZATION) {
			return true;
		}
		return false;
	}

	public boolean isSimulation() {
		if (simulation == TinkerUpdate.SIMULATION) {
			return true;
		}
		return false;
	}

	public boolean isStale() {
		for (Atom a : getAtomList()) {
			if (a.isStale()) {
				return true;
			}
		}
		return false;
	}

	public void removeKeyword(Keyword kd) {
		if (keywords.containsKey(kd.getKeyword())) {
			keywords.remove(kd.getKeyword());
		}
	}

	public void setClosing(boolean b) {
		closing = b;
	}

	public void setCommandDescription(String command) {
		commandDescription = command;
	}

	public void setEnergy(double e) {
		energy = e;
	}

	public void setKeyFile(File f) {
		keyFile = f;
	}

	public void setKeywords(Hashtable<String, Keyword> k) {
		keywords = k;
	}

	public void setLogFile(File f) {
		logFile = f;
	}

	public void setSimulation(int type) {
		simulation = type;
	}

	public void setStep(int s) {
		step = s;
	}

	public void setTemperature(double t) {
		temperature = t;
	}

	public void setTime(double t) {
		time = t;
	}

	public void setTrajectory(Trajectory t) {
		trajectory = t;
	}

	public String toFFString() {
		StringBuffer sb = new StringBuffer(toString());
		if (forceField != null) {
			String ff = forceField.toString("forcefield");
			if (ff != null) {
				ff = ff.substring(10).trim();
				sb.append(" (");
				sb.append(ff);
				sb.append(")");
			}
		}
		return sb.toString();
	}

	public String toFileString() {
		if (getFile() == null) {
			return toFFString();
		}
		StringBuffer sb = new StringBuffer(getFile().getAbsolutePath());
		if (forceField != null) {
			String ff = forceField.toString("forcefield");
			if (ff != null) {
				ff = ff.substring(10).trim();
				sb.append(" (");
				sb.append(ff);
				sb.append(")");
			}
		}
		return sb.toString();
	}

	public String toString() {
		if (getFile() != null) {
			if (commandDescription != null) {
				return getFile().getName() + " (" + commandDescription + ")";
			}
			return getFile().getName();
		}
		if (getName() != null) {
			if (commandDescription != null) {
				return getName() + commandDescription;
			}
			return getName();
		}
		return "FFX System";
	}
}
