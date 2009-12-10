/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 * @author Michael J. Schnieders
 * @version 0.1
 */
package ffx.parsers;

import java.util.ArrayList;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;

/**
 * The MergeFilter class allows Force Field X to treat merging of Systems
 * just like opening a file from a hard disk or socket.
 */
public class MergeFilter extends SystemFilter {
	public MergeFilter(MolecularAssembly f, ArrayList<Atom> a, ArrayList<Bond> b) {
		super(f);
		atomList = a;
		bondList = b;
	}

	/**
	 * 
	 */
	public boolean readFile() {
		return true;
	}

	/**
	 * 
	 */
	public boolean writeFile() {
		return false;
	}
}
