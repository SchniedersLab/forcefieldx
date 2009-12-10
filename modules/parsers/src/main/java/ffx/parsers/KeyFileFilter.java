/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 * @author Michael J. Schnieders
 * @version 0.1
 */
package ffx.parsers;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * The KeyFileFilter class is used to choose TINKER Keyword (*.KEY) files
 */
public final class KeyFileFilter extends FileFilter {
	/**
	 * Default Constructor
	 */
	public KeyFileFilter() {
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
		if (filename.regionMatches(false, dot + 1, "key", 0, 3)
				|| filename.regionMatches(false, dot + 1, "prm", 0, 3)) {
			return true;
		}
		return false;
	}

	/**
	 * Provides a description of this FileFilter
	 */
	public String getDescription() {
		return new String("Tinker Key Files: *.key");
	}
}
