// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
 * Implements a comparable pair, where a non-comparable T can be compared via a comparable S.
 *
 * @param <T> Some object
 * @param <S> Some indexing, comparable object
 * @author Jacob Litman
 */
public class ObjectPair<T, S extends Comparable<S>> implements Comparable<ObjectPair<T, S>> {

  private final T val;
  private final S key;

  /**
   * Constructor for ObjectPair.
   *
   * @param val a T object.
   * @param key a S object.
   */
  public ObjectPair(T val, S key) {
    this.val = val;
    this.key = key;
  }

  /**
   * sortAndReturn.
   *
   * @param theList a {@link java.util.List} object.
   * @param <U> a U object.
   * @param <V> a V object.
   * @return a {@link java.util.List} object.
   */
  public static <U, V extends Comparable<V>> List<U> sortAndReturn(List<ObjectPair<U, V>> theList) {
    theList.sort(null);
    List<U> retList = new ArrayList<>(theList.size());
    for (ObjectPair<U, V> e : theList) {
      retList.add(e.getVal());
    }
    return retList;
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(ObjectPair<T, S> o) {
    return key.compareTo(o.getKey());
  }

  /**
   * Getter for the field <code>key</code>.
   *
   * @return a S object.
   */
  public S getKey() {
    return key;
  }

  /**
   * Getter for the field <code>val</code>.
   *
   * @return a T object.
   */
  public T getVal() {
    return val;
  }
}
