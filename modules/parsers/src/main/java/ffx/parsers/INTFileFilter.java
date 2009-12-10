/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 * @author Michael J. Schnieders
 * @version 0.1
 */
package ffx.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.filechooser.FileFilter;

/**
 * The INTFileFilter class is used to choose TINKER Internal Coordinate (*.INT)
 * files
 */
public final class INTFileFilter extends FileFilter {
	/**
	 * Default Constructor
	 */
	public INTFileFilter() {
	}

	/**
	 * This method determines whether or not the parm File parameter is a Tinker
	 * *.xyz or not, returning true if it is. (Also returns true for any
	 * directory)
	 */
	public boolean accept(File parm) {
		if (parm.isDirectory()) {
			return true;
		}
		String filename = parm.getName().toLowerCase();
		int dot = filename.lastIndexOf(".");
		return filename.regionMatches(false, dot + 1, "int", 0, 3);
	}

	public boolean acceptDeep(File parm) {
		try {
			if (parm == null || parm.isDirectory() || !parm.canRead()) {
				return false;
			}
			FileReader fr = new FileReader(parm);
			BufferedReader br = new BufferedReader(fr);
			if (!br.ready()) {
				// Empty File?
				return false;
			}
			// If the first token is not an integer this file is not
			// an Internal Coordinates File.
			String rawdata = br.readLine();
			String header[] = rawdata.trim().split(" +");
			if (header == null || header.length == 0) {
				return false;
			}
			try {
				Integer.parseInt(header[0]);
			} catch (Exception e) {
				return false;
			}
			// If the the first Atom line does not begin with an integer and
			// contain
			// three tokens, it is not an internal coordinate file.
			String firstAtom = br.readLine();
			if (firstAtom == null) {
				return false;
			}
			br.close();
			fr.close();
			String data[] = firstAtom.trim().split(" +");
			if (data == null || data.length != 3) {
				return false;
			}
			try {
				Integer.parseInt(data[0]);
			} catch (Exception e) {
				return false;
			}
			return true;
		} catch (Exception e) {
			return true;
		}
	}

	/**
	 * Provides a description of this FileFilter
	 */
	public String getDescription() {
		return new String("Tinker Internal Coordinate Files: *.int");
	}
}
