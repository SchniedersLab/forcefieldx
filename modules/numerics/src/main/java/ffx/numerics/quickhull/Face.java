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
 * Basic triangular face used to form the hull.
 * <p>
 * The information stored for each face consists of a planar normal, a planar
 * offset, and a doubly-linked list of three <a href=HalfEdge>HalfEdges</a>
 * which surround the face in a counter-clockwise direction.
 *
 * @author John E. Lloyd, Fall 2004
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Face {

  protected static final int DELETED = 3;

  protected static final int NON_CONVEX = 2;

  protected static final int VISIBLE = 1;

  protected double area;

  protected HalfEdge he0;

  protected int mark;

  protected Face next;

  protected int numVerts;

  protected Vertex outside;

  protected double planeOffset;

  private final Point3d centroid;

  private final Vector3d normal;

  public Face() {
    normal = new Vector3d();
    centroid = new Point3d();
    mark = VISIBLE;
  }

  /**
   * Creates a Face by linking the specified indices of a vertex array into a closed
   * counter-clockwise half-edge loop and computing its normal and centroid.
   *
   * @param vtxArray array of vertices
   * @param indices  indices into vtxArray specifying the face vertices in CCW order
   * @return the newly created Face
   */
  public static Face create(Vertex[] vtxArray, int[] indices) {
    Face face = new Face();
    HalfEdge hePrev = null;
    for (int index : indices) {
      HalfEdge he = new HalfEdge(vtxArray[index], face);
      if (hePrev != null) {
        he.setPrev(hePrev);
        hePrev.setNext(he);
      } else {
        face.he0 = he;
      }
      hePrev = he;
    }
    face.he0.setPrev(hePrev);
    hePrev.setNext(face.he0);

    // compute the normal and offset
    face.computeNormalAndCentroid();
    return face;
  }

  /**
   * Convenience method to create a triangular Face from three vertices.
   *
   * @param v0 first vertex
   * @param v1 second vertex
   * @param v2 third vertex
   * @return the newly created triangular Face
   */
  public static Face createTriangle(Vertex v0, Vertex v1, Vertex v2) {
    return createTriangle(v0, v1, v2, 0);
  }

  /**
   * Constructs a triangular Face from vertices v0, v1, and v2 and computes its
   * normal and centroid. If the computed area is below minArea, the normal is
   * adjusted for robustness against near-collinearity.
   *
   * @param v0      first vertex
   * @param v1      second vertex
   * @param v2      third vertex
   * @param minArea minimum area threshold used to stabilize the normal
   * @return the newly created triangular Face
   */
  public static Face createTriangle(Vertex v0, Vertex v1, Vertex v2, double minArea) {
    Face face = new Face();
    HalfEdge he0 = new HalfEdge(v0, face);
    HalfEdge he1 = new HalfEdge(v1, face);
    HalfEdge he2 = new HalfEdge(v2, face);

    he0.prev = he2;
    he0.next = he1;
    he1.prev = he0;
    he1.next = he2;
    he2.prev = he1;
    he2.next = he0;

    face.he0 = he0;

    // compute the normal and offset
    face.computeNormalAndCentroid(minArea);
    return face;
  }

  /**
   * Computes the centroid (arithmetic mean of vertices) of this face into the
   * provided point.
   *
   * @param centroid output parameter to receive the centroid
   */
  public void computeCentroid(Point3d centroid) {
    centroid.setZero();
    HalfEdge he = he0;
    do {
      centroid.add(he.head().pnt);
      he = he.next;
    } while (he != he0);
    centroid.scale(1 / (double) numVerts);
  }

  /**
   * Computes a unit-length normal vector for this face using the vertex winding
   * and writes it into the provided vector. Also updates the face area and number
   * of vertices encountered.
   *
   * @param normal output parameter to receive the unit normal
   */
  public void computeNormal(Vector3d normal) {
    HalfEdge he1 = he0.next;
    HalfEdge he2 = he1.next;

    Point3d p0 = he0.head().pnt;
    Point3d p2 = he1.head().pnt;

    double d2x = p2.x - p0.x;
    double d2y = p2.y - p0.y;
    double d2z = p2.z - p0.z;

    normal.setZero();

    numVerts = 2;

    while (he2 != he0) {
      double d1x = d2x;
      double d1y = d2y;
      double d1z = d2z;

      p2 = he2.head().pnt;
      d2x = p2.x - p0.x;
      d2y = p2.y - p0.y;
      d2z = p2.z - p0.z;

      normal.x += d1y * d2z - d1z * d2y;
      normal.y += d1z * d2x - d1x * d2z;
      normal.z += d1x * d2y - d1y * d2x;

      he1 = he2;
      he2 = he2.next;
      numVerts++;
    }
    area = normal.norm();
    normal.scale(1 / area);
  }

  /**
   * Computes a unit-length normal for this face, and if the preliminary area is
   * below the specified minArea, adjusts the normal to be more orthogonal to the
   * longest edge to improve robustness.
   *
   * @param normal  output parameter to receive the unit normal
   * @param minArea minimum area threshold used to stabilize the normal computation
   */
  public void computeNormal(Vector3d normal, double minArea) {
    computeNormal(normal);

    if (area < minArea) {
      // make the normal more robust by removing
      // components parallel to the longest edge

      HalfEdge hedgeMax = null;
      double lenSqrMax = 0;
      HalfEdge hedge = he0;
      do {
        double lenSqr = hedge.lengthSquared();
        if (lenSqr > lenSqrMax) {
          hedgeMax = hedge;
          lenSqrMax = lenSqr;
        }
        hedge = hedge.next;
      } while (hedge != he0);

      Point3d p2 = hedgeMax.head().pnt;
      Point3d p1 = hedgeMax.tail().pnt;
      double lenMax = Math.sqrt(lenSqrMax);
      double ux = (p2.x - p1.x) / lenMax;
      double uy = (p2.y - p1.y) / lenMax;
      double uz = (p2.z - p1.z) / lenMax;
      double dot = normal.x * ux + normal.y * uy + normal.z * uz;
      normal.x -= dot * ux;
      normal.y -= dot * uy;
      normal.z -= dot * uz;

      normal.normalize();
    }
  }

  /**
   * Computes the distance from a point p to the plane of this face.
   *
   * @param p the point
   * @return distance from the point to the plane
   */
  public double distanceToPlane(Point3d p) {
    return normal.x * p.x + normal.y * p.y + normal.z * p.z - planeOffset;
  }

  /**
   * Finds the half-edge within this face which has tail <code>vt</code> and
   * head <code>vh</code>.
   *
   * @param vt tail point
   * @param vh head point
   * @return the half-edge, or null if none is found.
   */
  public HalfEdge findEdge(Vertex vt, Vertex vh) {
    HalfEdge he = he0;
    do {
      if (he.head() == vh && he.tail() == vt) {
        return he;
      }
      he = he.next;
    } while (he != he0);
    return null;
  }

  /**
   * Returns the centroid previously computed for this face.
   *
   * @return the centroid point
   */
  public Point3d getCentroid() {
    return centroid;
  }

  /**
   * Gets the i-th half-edge associated with the face.
   *
   * @param i the half-edge index, in the range 0-2.
   * @return the half-edge
   */
  public HalfEdge getEdge(int i) {
    HalfEdge he = he0;
    while (i > 0) {
      he = he.next;
      i--;
    }
    while (i < 0) {
      he = he.prev;
      i++;
    }
    return he;
  }

  /**
   * Returns the first half-edge of this face (the start of the circular list).
   *
   * @return the first HalfEdge in the face
   */
  public HalfEdge getFirstEdge() {
    return he0;
  }

  /**
   * Returns the normal of the plane associated with this face.
   *
   * @return the planar normal
   */
  public Vector3d getNormal() {
    return normal;
  }

  /**
   * Writes the vertex indices of this face into the provided array in CCW order.
   * The array must be large enough to hold numVertices() entries.
   *
   * @param idxs output array receiving the vertex indices
   */
  public void getVertexIndices(int[] idxs) {
    HalfEdge he = he0;
    int i = 0;
    do {
      idxs[i++] = he.head().index;
      he = he.next;
    } while (he != he0);
  }

  /**
   * Returns a space-separated string of the vertex indices defining this face
   * in counter-clockwise order.
   *
   * @return string listing the vertex indices in CCW order
   */
  public String getVertexString() {
    StringBuilder s = null;
    HalfEdge he = he0;
    do {
      if (s == null) {
        s = new StringBuilder("" + he.head().index);
      } else {
        s.append(" ").append(he.head().index);
      }
      he = he.next;
    } while (he != he0);
    return s.toString();
  }

  /**
   * Merges this face with the adjacent face across the specified half-edge if possible.
   * Updates connectivity, recomputes normals/centroids, and records any faces that
   * become redundant in the provided discarded array.
   *
   * @param hedgeAdj the half-edge along which to merge with the adjacent face
   * @param discarded an output array to collect faces that are deleted by the merge
   * @return the number of faces recorded in discarded (0, 1, or 2)
   */
  public int mergeAdjacentFace(HalfEdge hedgeAdj, Face[] discarded) {
    Face oppFace = hedgeAdj.oppositeFace();
    int numDiscarded = 0;

    discarded[numDiscarded++] = oppFace;
    oppFace.mark = DELETED;

    HalfEdge hedgeOpp = hedgeAdj.getOpposite();

    HalfEdge hedgeAdjPrev = hedgeAdj.prev;
    HalfEdge hedgeAdjNext = hedgeAdj.next;
    HalfEdge hedgeOppPrev = hedgeOpp.prev;
    HalfEdge hedgeOppNext = hedgeOpp.next;

    while (hedgeAdjPrev.oppositeFace() == oppFace) {
      hedgeAdjPrev = hedgeAdjPrev.prev;
      hedgeOppNext = hedgeOppNext.next;
    }

    while (hedgeAdjNext.oppositeFace() == oppFace) {
      hedgeOppPrev = hedgeOppPrev.prev;
      hedgeAdjNext = hedgeAdjNext.next;
    }

    HalfEdge hedge;

    for (hedge = hedgeOppNext; hedge != hedgeOppPrev.next; hedge = hedge.next) {
      hedge.face = this;
    }

    if (hedgeAdj == he0) {
      he0 = hedgeAdjNext;
    }

    // handle the half edges at the head
    Face discardedFace;

    discardedFace = connectHalfEdges(hedgeOppPrev, hedgeAdjNext);
    if (discardedFace != null) {
      discarded[numDiscarded++] = discardedFace;
    }

    // handle the half edges at the tail
    discardedFace = connectHalfEdges(hedgeAdjPrev, hedgeOppNext);
    if (discardedFace != null) {
      discarded[numDiscarded++] = discardedFace;
    }

    computeNormalAndCentroid();
    checkConsistency();

    return numDiscarded;
  }

  /**
   * Returns the number of vertices (half-edges) bounding this face.
   *
   * @return the vertex count of this face
   */
  public int numVertices() {
    return numVerts;
  }

  /**
   * Triangulates this (convex polygonal) face into a fan of triangles sharing the first
   * vertex. Newly created faces are added to newFaces, and appropriate opposite links are
   * established. Normals/centroids are updated using the given minArea for stability.
   *
   * @param newFaces list that collects newly created triangular faces
   * @param minArea  minimum area threshold used when computing normals
   */
  public void triangulate(FaceList newFaces, double minArea) {
    HalfEdge hedge;

    if (numVertices() < 4) {
      return;
    }

    Vertex v0 = he0.head();

    hedge = he0.next;
    HalfEdge oppPrev = hedge.opposite;
    Face face0 = null;

    for (hedge = hedge.next; hedge != he0.prev; hedge = hedge.next) {
      Face face = createTriangle(v0, hedge.prev.head(), hedge.head(), minArea);
      face.he0.next.setOpposite(oppPrev);
      face.he0.prev.setOpposite(hedge.opposite);
      oppPrev = face.he0;
      newFaces.add(face);
      if (face0 == null) {
        face0 = face;
      }
    }
    hedge = new HalfEdge(he0.prev.prev.head(), this);
    hedge.setOpposite(oppPrev);

    hedge.prev = he0;
    hedge.prev.next = hedge;

    hedge.next = he0.prev;
    hedge.next.prev = hedge;

    computeNormalAndCentroid(minArea);
    checkConsistency();

    for (Face face = face0; face != null; face = face.next) {
      face.checkConsistency();
    }

  }

  /**
   * Computes the squared area of the triangle defined by hedge0 (tail->head)
   * and the point at the head of hedge1. Useful for robust comparisons without
   * taking a square root.
   *
   * @param hedge0 the reference half-edge whose tail and head form two triangle vertices
   * @param hedge1 a half-edge providing the third vertex via its head
   * @return the squared area of the resulting triangle
   */
  public double areaSquared(HalfEdge hedge0, HalfEdge hedge1) {
    Point3d p0 = hedge0.tail().pnt;
    Point3d p1 = hedge0.head().pnt;
    Point3d p2 = hedge1.head().pnt;

    double dx1 = p1.x - p0.x;
    double dy1 = p1.y - p0.y;
    double dz1 = p1.z - p0.z;

    double dx2 = p2.x - p0.x;
    double dy2 = p2.y - p0.y;
    double dz2 = p2.z - p0.z;

    double x = dy1 * dz2 - dz1 * dy2;
    double y = dz1 * dx2 - dx1 * dz2;
    double z = dx1 * dy2 - dy1 * dx2;

    return x * x + y * y + z * z;
  }

  private void computeNormalAndCentroid() {
    computeNormal(normal);
    computeCentroid(centroid);
    planeOffset = normal.dot(centroid);
    int numv = 0;
    HalfEdge he = he0;
    do {
      numv++;
      he = he.next;
    } while (he != he0);
    if (numv != numVerts) {
      throw new InternalErrorException("face " + getVertexString() + " numVerts=" + numVerts + " should be " + numv);
    }
  }

  private void computeNormalAndCentroid(double minArea) {
    computeNormal(normal, minArea);
    computeCentroid(centroid);
    planeOffset = normal.dot(centroid);
  }

  private Face connectHalfEdges(HalfEdge hedgePrev, HalfEdge hedge) {
    Face discardedFace = null;

    if (hedgePrev.oppositeFace() == hedge.oppositeFace()) { // then there is
      // a redundant
      // edge that we
      // can get rid
      // off

      Face oppFace = hedge.oppositeFace();
      HalfEdge hedgeOpp;

      if (hedgePrev == he0) {
        he0 = hedge;
      }
      if (oppFace.numVertices() == 3) { // then we can get rid of the
        // opposite face altogether
        hedgeOpp = hedge.getOpposite().prev.getOpposite();

        oppFace.mark = DELETED;
        discardedFace = oppFace;
      } else {
        hedgeOpp = hedge.getOpposite().next;

        if (oppFace.he0 == hedgeOpp.prev) {
          oppFace.he0 = hedgeOpp;
        }
        hedgeOpp.prev = hedgeOpp.prev.prev;
        hedgeOpp.prev.next = hedgeOpp;
      }
      hedge.prev = hedgePrev.prev;
      hedge.prev.next = hedge;

      hedge.opposite = hedgeOpp;
      hedgeOpp.opposite = hedge;

      // oppFace was modified, so need to recompute
      oppFace.computeNormalAndCentroid();
    } else {
      hedgePrev.next = hedge;
      hedge.prev = hedgePrev;
    }
    return discardedFace;
  }

  /**
   * Performs sanity checks on this face's topology and geometry, verifying
   * opposite half-edges, face markings, and vertex counts.
   *
   * @throws InternalErrorException if inconsistencies are detected
   */
  void checkConsistency() {
    // do a sanity check on the face
    HalfEdge hedge = he0;
    double maxd = 0;
    int numv = 0;

    if (numVerts < 3) {
      throw new InternalErrorException("degenerate face: " + getVertexString());
    }
    do {
      HalfEdge hedgeOpp = hedge.getOpposite();
      if (hedgeOpp == null) {
        throw new InternalErrorException("face " + getVertexString() + ": " + "unreflected half edge " + hedge.getVertexString());
      } else if (hedgeOpp.getOpposite() != hedge) {
        throw new InternalErrorException("face " + getVertexString() + ": " + "opposite half edge " + hedgeOpp.getVertexString() + " has opposite "
            + hedgeOpp.getOpposite().getVertexString());
      }
      if (hedgeOpp.head() != hedge.tail() || hedge.head() != hedgeOpp.tail()) {
        throw new InternalErrorException("face " + getVertexString() + ": " + "half edge " + hedge.getVertexString() + " reflected by " + hedgeOpp.getVertexString());
      }
      Face oppFace = hedgeOpp.face;
      if (oppFace == null) {
        throw new InternalErrorException("face " + getVertexString() + ": " + "no face on half edge " + hedgeOpp.getVertexString());
      } else if (oppFace.mark == DELETED) {
        throw new InternalErrorException("face " + getVertexString() + ": " + "opposite face " + oppFace.getVertexString() + " not on hull");
      }
      double d = Math.abs(distanceToPlane(hedge.head().pnt));
      if (d > maxd) {
        maxd = d;
      }
      numv++;
      hedge = hedge.next;
    } while (hedge != he0);

    if (numv != numVerts) {
      throw new InternalErrorException("face " + getVertexString() + " numVerts=" + numVerts + " should be " + numv);
    }

  }
}
