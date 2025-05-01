// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.utilities;

import java.util.ArrayList;
import java.util.List;

/**
 * ObjectPair record.
 * <p>
 * Implements a comparable pair, where a non-comparable T can be compared via a comparable S.
 *
 * @param <T> The class of the non-comparable object.
 * @param <S> The class of the indexing, comparable object.
 * @param val The non-comparable object.
 * @param key The indexing, comparable object.
 * @author Jacob Litman
 * @since 1.0
 */
public record ObjectPair<T, S extends Comparable<S>>(T val, S key) implements Comparable<ObjectPair<T, S>> {

  /**
   * sortAndReturn.
   *
   * @param theList a {@link List} object.
   * @param <U>     Some object.
   * @param <V>     Some indexing, comparable object
   * @return the sorted {@link List} of objects.
   */
  public static <U, V extends Comparable<V>> List<U> sortAndReturn(List<ObjectPair<U, V>> theList) {
    theList.sort(null);
    List<U> retList = new ArrayList<>(theList.size());
    for (ObjectPair<U, V> e : theList) {
      retList.add(e.val());
    }
    return retList;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compareTo(ObjectPair<T, S> o) {
    return key.compareTo(o.key());
  }

}
