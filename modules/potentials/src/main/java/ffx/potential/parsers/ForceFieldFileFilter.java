/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 * @author Michael J. Schnieders
 * @version 0.1
 */
package ffx.potential.parsers;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * The ForceFieldFileFilter class is used to choose TINKER Parameter (*.PRM)
 * files
 */
public final class ForceFieldFileFilter extends FileFilter {
	/**
	 * Default Constructor
	 */
	public ForceFieldFileFilter() {
	}

	/**
	 * This method determines whether or not the parm File parameter is a Tinker
	 * *.key or not, returning true if it is. (Also returns true for any
	 * directory)
	 */
	public boolean accept(File parm) {
		if (parm.isDirectory()) {
			return true;
		}
		String filename = parm.getName();
		int dot = filename.lastIndexOf(".");
		if (filename.regionMatches(false, dot + 1, "prm", 0, 3)) {
			return true;
		}
		return false;
	}

	/**
	 * Provides a description of this FileFilter
	 */
	public String getDescription() {
		return new String("Tinker Parameter Files: *.prm");
	}
}
