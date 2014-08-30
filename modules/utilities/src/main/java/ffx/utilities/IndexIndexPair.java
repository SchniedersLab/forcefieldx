/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.utilities;

/**
 * <p>
 * IndexIndexPair class.</p>
 *
 * @author Jacob M. Litman
 *
 * @since 1.0
 */
public class IndexIndexPair implements Comparable {

    private final int sortedIndex;
    private final int referenceIndex;

    /**
     * Pass in an int to be sorted on, then an int referring to some original
     * list. Enables sorting based on one integer value while remembering where
     * it was in the original list.
     *
     * @param sortedIndex Value to be sorted upon
     * @param referenceIndex Original index
     */
    public IndexIndexPair(int sortedIndex, int referenceIndex) {
        this.sortedIndex = sortedIndex;
        this.referenceIndex = referenceIndex;
    }

    public int getSortedIndex() {
        return sortedIndex;
    }

    public int getReferenceIndex() {
        return referenceIndex;
    }

    @Override
    public int compareTo(Object o) {
        if (o == null) {
            return 0;
        }
        if (!(o instanceof IndexIndexPair)) {
            return 0;
        }
        IndexIndexPair other = (IndexIndexPair) o;
        if (sortedIndex < other.sortedIndex) {
            return -1;
        } else if (sortedIndex > other.sortedIndex) {
            return 1;
        } else {
            return 0;
        }
    }
}
