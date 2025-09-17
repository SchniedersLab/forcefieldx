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
package ffx.numerics.quickhull;

/**
 * Represents the half-edges that surround each face in a counter-clockwise
 * direction.
 *
 * @author John E. Lloyd, Fall 2004
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class HalfEdge {

  /**
   * The vertex associated with the head of this half-edge.
   */
  protected Vertex vertex;

  /**
   * Triangular face associated with this half-edge.
   */
  protected Face face;

  /**
   * Next half-edge in the triangle.
   */
  protected HalfEdge next;

  /**
   * Previous half-edge in the triangle.
   */
  protected HalfEdge prev;

  /**
   * Half-edge associated with the opposite triangle adjacent to this edge.
   */
  protected HalfEdge opposite;

  /**
   * Constructs a HalfEdge with head vertex <code>v</code> and left-hand
   * triangular face <code>f</code>.
   *
   * @param v head vertex
   * @param f left-hand triangular face
   */
  public HalfEdge(Vertex v, Face f) {
    vertex = v;
    face = f;
  }

  /**
   * Constructs an uninitialized HalfEdge.
   * Fields may be set later via mutators.
   */
  public HalfEdge() {
  }

  /**
   * Sets the value of the next edge adjacent (counter-clockwise) to this one
   * within the triangle.
   *
   * @param edge next adjacent edge
   */
  public void setNext(HalfEdge edge) {
    next = edge;
  }

  /**
   * Gets the value of the next edge adjacent (counter-clockwise) to this one
   * within the triangle.
   *
   * @return next adjacent edge
   */
  public HalfEdge getNext() {
    return next;
  }

  /**
   * Sets the value of the previous edge adjacent (clockwise) to this one
   * within the triangle.
   *
   * @param edge previous adjacent edge
   */
  public void setPrev(HalfEdge edge) {
    prev = edge;
  }

  /**
   * Gets the value of the previous edge adjacent (clockwise) to this one
   * within the triangle.
   *
   * @return previous adjacent edge
   */
  public HalfEdge getPrev() {
    return prev;
  }

  /**
   * Returns the triangular face located to the left of this half-edge.
   *
   * @return left-hand triangular face
   */
  public Face getFace() {
    return face;
  }

  /**
   * Returns the half-edge opposite to this half-edge.
   *
   * @return opposite half-edge
   */
  public HalfEdge getOpposite() {
    return opposite;
  }

  /**
   * Sets the half-edge opposite to this half-edge.
   *
   * @param edge opposite half-edge
   */
  public void setOpposite(HalfEdge edge) {
    opposite = edge;
    edge.opposite = this;
  }

  /**
   * Returns the head vertex associated with this half-edge.
   *
   * @return head vertex
   */
  public Vertex head() {
    return vertex;
  }

  /**
   * Returns the tail vertex associated with this half-edge.
   *
   * @return tail vertex
   */
  public Vertex tail() {
    return prev != null ? prev.vertex : null;
  }

  /**
   * Returns the opposite triangular face associated with this half-edge.
   *
   * @return opposite triangular face
   */
  public Face oppositeFace() {
    return opposite != null ? opposite.face : null;
  }

  /**
   * Produces a string identifying this half-edge by the point index values of
   * its tail and head vertices.
   *
   * @return identifying string
   */
  public String getVertexString() {
    if (tail() != null) {
      return tail().index + "-" + head().index;
    } else {
      return "?-" + head().index;
    }
  }

  /**
   * Returns the length of this half-edge.
   *
   * @return half-edge length
   */
  public double length() {
    if (tail() != null) {
      return head().pnt.distance(tail().pnt);
    } else {
      return -1;
    }
  }

  /**
   * Returns the length squared of this half-edge.
   *
   * @return half-edge length squared
   */
  public double lengthSquared() {
    if (tail() != null) {
      return head().pnt.distanceSquared(tail().pnt);
    } else {
      return -1;
    }
  }
}
