// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.numerics.quickhull;

/**
 * Maintains a double-linked list of vertices for use by QuickHull3D.
 *
 * @author John E. Lloyd, Fall 2004
 * @author Michael J. Schnieders
 * @since 1.0
 */
class VertexList {

  private Vertex head;

  private Vertex tail;

  /**
   * Clears this list.
   */
  public void clear() {
    head = tail = null;
  }

  /**
   * Adds a vertex to the end of this list.
   *
   * @param vtx Vertex to add
   */
  public void add(Vertex vtx) {
    if (head == null) {
      head = vtx;
    } else {
      tail.next = vtx;
    }
    vtx.prev = tail;
    vtx.next = null;
    tail = vtx;
  }

  /**
   * Adds a chain of vertices to the end of this list.
   * The provided vertex is assumed to be the head of the chain.
   *
   * @param vtx head of the vertex chain to append
   */
  public void addAll(Vertex vtx) {
    if (head == null) {
      head = vtx;
    } else {
      tail.next = vtx;
    }
    vtx.prev = tail;
    while (vtx.next != null) {
      vtx = vtx.next;
    }
    tail = vtx;
  }

  /**
   * Deletes a single vertex from this list.
   *
   * @param vtx the vertex to remove
   */
  public void delete(Vertex vtx) {
    if (vtx.prev == null) {
      head = vtx.next;
    } else {
      vtx.prev.next = vtx.next;
    }
    if (vtx.next == null) {
      tail = vtx.prev;
    } else {
      vtx.next.prev = vtx.prev;
    }
  }

  /**
   * Deletes a contiguous chain of vertices from this list.
   *
   * @param vtx1 the first vertex in the chain to remove
   * @param vtx2 the last vertex in the chain to remove
   */
  public void delete(Vertex vtx1, Vertex vtx2) {
    if (vtx1.prev == null) {
      head = vtx2.next;
    } else {
      vtx1.prev.next = vtx2.next;
    }
    if (vtx2.next == null) {
      tail = vtx1.prev;
    } else {
      vtx2.next.prev = vtx1.prev;
    }
  }

  /**
   * Inserts a vertex into this list before another specified vertex.
   *
   * @param vtx  the vertex to insert
   * @param next the vertex before which to insert vtx
   */
  public void insertBefore(Vertex vtx, Vertex next) {
    vtx.prev = next.prev;
    if (next.prev == null) {
      head = vtx;
    } else {
      next.prev.next = vtx;
    }
    vtx.next = next;
    next.prev = vtx;
  }

  /**
   * Returns the first vertex in this list (head), or null if empty.
   *
   * @return the first Vertex, or null if the list is empty
   */
  public Vertex first() {
    return head;
  }

  /**
   * Returns true if this list is empty.
   *
   * @return true if there are no vertices in the list
   */
  public boolean isEmpty() {
    return head == null;
  }
}
