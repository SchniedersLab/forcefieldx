/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.bonded;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * The AtomIterator class provides an Iterator over an MSNode's Atoms, so as to
 * fulfill Biojava API. Intended to be used with classes which implement the
 * Biojava Group interface.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 *
 */
public class AtomIterator implements Iterator<org.biojava.nbio.structure.Atom> {
        
        private final Atom[] atoms;
        private final int nAtoms;
        private int count;
        
        public AtomIterator(MSNode node) {
            List<Atom> atList = node.getAtomList();
            nAtoms = atList.size();
            atoms = new Atom[nAtoms];
            atList.toArray(atoms);
            count = 0;
        }
        
        @Override
        public boolean hasNext() {
            return count < nAtoms;
        }

        @Override
        public org.biojava.nbio.structure.Atom next() throws NoSuchElementException {
            if (count >= nAtoms) {
                throw new NoSuchElementException("End of Residue atoms reached.");
            }
            return atoms[count++];
        }
        
        @Override
        public void remove() throws IllegalArgumentException {
            if (count == 0) {
                throw new IllegalStateException("No atoms returned yet; cannot remove!");
            }
            if (atoms[count-1] == null) {
                throw new IllegalStateException("Last atom already removed!");
            }
            atoms[count-1].removeFromParent();
            atoms[count-1].setGroup(null, true);
            atoms[count-1] = null;
        }
    }

