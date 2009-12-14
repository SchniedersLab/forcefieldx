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
 * The PDBXFileFilter class is used to choose a PDBML File (*.XML).
 */
public final class PDBMLFileFilter extends FileFilter {
	/**
	 * Default Constructor
	 */
	public PDBMLFileFilter() {
	}

	/**
	 * This method determines whether or not the File is PDBML or not, returning
	 * true if it is. (Also returns true for any directory)
	 */
	public boolean accept(File parm) {
		if (parm.isDirectory()) {
			return true;
		}
		String filename = parm.getName().toUpperCase();
		int dot = filename.lastIndexOf(".");
		if (filename.regionMatches(false, dot + 1, "XML", 0, 3)) {
			return true;
		}
		return false;
	}

	/**
	 * Provides a description of this FileFilter
	 */
	public String getDescription() {
		return new String("PDBML Files: *.xml");
	}
}
