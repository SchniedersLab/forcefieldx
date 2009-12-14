/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 * @author Michael J. Schnieders
 * @version 0.1
 */
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Atom;

/**
 * The DYNFilter class parses TINKER Restart (*.DYN) files.
 */
public class DYNFilter {
	MolecularAssembly fsystem;
	File file;

	public DYNFilter(MolecularAssembly s, File f) {
		fsystem = s;
		file = f;
	}

	public boolean read() {
		if (!file.exists() || !file.canRead()) {
			return false;
		}
		try {
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			br.readLine();
			String data = br.readLine().trim();
			String tokens[] = data.split(" +");
			if (tokens.length == 0) {
				return false;
			}
			int numatoms = Integer.parseInt(tokens[0]);
			if (numatoms != fsystem.getAtomList().size()) {
				return false;
			}
			br.readLine();
			data = br.readLine().trim();
			tokens = data.split(" +");
			if (tokens.length != 3) {
				return false;
			}
			double d[] = new double[3];
			d[0] = Double.parseDouble(tokens[0]);
			d[1] = Double.parseDouble(tokens[1]);
			d[2] = Double.parseDouble(tokens[2]);
			fsystem.setBox(d);
			data = br.readLine().trim();
			tokens = data.split(" +");
			if (tokens.length != 3) {
				return false;
			}
			d[0] = Double.parseDouble(tokens[0]);
			d[1] = Double.parseDouble(tokens[1]);
			d[2] = Double.parseDouble(tokens[2]);
			fsystem.setAngle(d);
			// Positions
			br.readLine();
			double x[][] = new double[numatoms][3];
			for (int i = 0; i < numatoms; i++) {
				data = br.readLine().trim();
				tokens = data.split(" +");
				if (tokens.length != 3) {
					return false;
				}
				x[i][0] = Double.parseDouble(tokens[0]);
				x[i][1] = Double.parseDouble(tokens[1]);
				x[i][2] = Double.parseDouble(tokens[2]);
			}
			List<Atom> atoms = fsystem.getAtomList();
			for (Atom a : atoms) {
				int j = a.getXYZIndex() - 1;
				a.moveTo(x[j][0], x[j][1], x[j][2]);
			}
			// Velocities
			br.readLine();
			for (int i = 0; i < numatoms; i++) {
				data = br.readLine().trim();
				tokens = data.split(" +");
				if (tokens.length != 3) {
					return false;
				}
				x[i][0] = Double.parseDouble(tokens[0]);
				x[i][1] = Double.parseDouble(tokens[1]);
				x[i][2] = Double.parseDouble(tokens[2]);
			}
			for (Atom a : atoms) {
				int j = a.getXYZIndex() - 1;
				// a.setVeclocity(x[j][0], x[j][1], x[j][2]);
			}
			// Accelerations
			br.readLine();
			for (int i = 0; i < numatoms; i++) {
				data = br.readLine().trim();
				tokens = data.split(" +");
				if (tokens.length != 3) {
					return false;
				}
				x[i][0] = Double.parseDouble(tokens[0]);
				x[i][1] = Double.parseDouble(tokens[1]);
				x[i][2] = Double.parseDouble(tokens[2]);
			}
			for (Atom a : atoms) {
				int j = a.getXYZIndex() - 1;
				// a.setAcceleration(x[j][0], x[j][1], x[j][2]);
			}
			br.close();
			fr.close();
		} catch (Exception e) {
			return false;
		}
		return true;
	}
}
