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
 * The PDBFileFilter class is used to choose Protein Databank (*.PDB) files
 */
public final class PDBFileFilter extends FileFilter {
	/**
	 * Default Constructor
	 */
	public PDBFileFilter() {
	}

	/**
	 * This method determines whether or not the parm File parameter is a
	 * Protein Database *.pdb file or not, returning true if it is. (Also
	 * returns true for any directory)
	 */
	public boolean accept(File parm) {
		if (parm.isDirectory()) {
			return true;
		}
		String filename = parm.getName();
		int dot = filename.lastIndexOf(".");
		return filename.regionMatches(false, dot + 1, "pdb", 0, 3);
	}

	/**
	 * Provides a description of this FileFilter
	 */
	public String getDescription() {
		return new String("Protein Data Bank Files: *.pdb");
	}
}
