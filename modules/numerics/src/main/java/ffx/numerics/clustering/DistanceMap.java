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
package ffx.numerics.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.logging.Logger;

/**
 * Container for linkages
 * with the minimal methods needed in the package
 * Created by Alexandre Masselot on 7/18/14.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DistanceMap {

  public static final Logger logger = Logger.getLogger(DistanceMap.class.getName());

  private final Map<String, Item> pairHash;
  private final PriorityQueue<Item> data;

  private class Item implements Comparable<Item> {
    final ClusterPair pair;
    final String hash;
    boolean removed = false;

    Item(ClusterPair p) {
      pair = p;
      hash = hashCodePair(p);
    }

    @Override
    public int compareTo(Item o) {
      return pair.compareTo(o.pair);
    }

    @Override
    public String toString() {
      return hash;
    }
  }

  public DistanceMap() {
    data = new PriorityQueue<>();
    pairHash = new HashMap<>();
  }

  /**
   * Returns a snapshot list of all cluster pairs currently stored.
   *
   * @return list of ClusterPair entries in this map
   */
  public List<ClusterPair> list() {
    List<ClusterPair> l = new ArrayList<>(data.size());
    for (Item clusterPair : data) {
      l.add(clusterPair.pair);
    }
    return l;
  }

  /**
   * Finds the ClusterPair for the two provided clusters.
   *
   * @param c1 the first cluster
   * @param c2 the second cluster
   * @return the matching ClusterPair, or null if absent
   */
  public ClusterPair findByCodePair(Cluster c1, Cluster c2) {
    String inCode = hashCodePair(c1, c2);
    Item item = pairHash.get(inCode);
    return item == null ? null : item.pair;
  }

  /**
   * Removes and returns the minimal-distance pair (according to priority queue ordering).
   *
   * @return the next ClusterPair, or null if none
   */
  public ClusterPair removeFirst() {
    Item poll = data.poll();
    while (poll != null && poll.removed) {
      poll = data.poll();
    }
    if (poll == null) {
      return null;
    }
    ClusterPair link = poll.pair;
    pairHash.remove(poll.hash);
    return link;
  }

  /**
   * Marks the given ClusterPair as removed (lazy removal) and drops it from the hash index.
   *
   * @param link the ClusterPair to remove
   * @return true if the pair was present and marked removed; false otherwise
   */
  public boolean remove(ClusterPair link) {
    Item remove = pairHash.remove(hashCodePair(link));
    if (remove == null) {
      return false;
    }
    remove.removed = true;
    // data.remove(remove);  // bottleneck
    return true;
  }

  /**
   * Adds a new ClusterPair if no equivalent pair already exists.
   *
   * @param link the pair to add
   * @return true if added; false if a duplicate existed
   */
  public boolean add(ClusterPair link) {
    Item e = new Item(link);
    Item existingItem = pairHash.get(e.hash);
    if (existingItem != null) {
      logger.info("hashCode = " + existingItem.hash +
          " adding redundant link:" + link + " (exist:" + existingItem + ")");
      return false;
    } else {
      pairHash.put(e.hash, e);
      data.add(e);
      return true;
    }
  }

  /**
   * Peek at the minimal linkage distance currently in the map.
   *
   * @return the smallest linkage distance, or null if empty
   */
  public Double minDist() {
    Item peek = data.peek();
    if (peek != null)
      return peek.pair.getLinkageDistance();
    else
      return null;
  }

  /**
   * Compute some kind of unique ID for a given cluster pair.
   *
   * @return The ID
   */
  private String hashCodePair(ClusterPair link) {
    return hashCodePair(link.getlCluster(), link.getrCluster());
  }

  private String hashCodePair(Cluster lCluster, Cluster rCluster) {
    return hashCodePairNames(lCluster.getName(), rCluster.getName());
  }

  private String hashCodePairNames(String lName, String rName) {
    if (lName.compareTo(rName) < 0) {
      return lName + "~~~" + rName;//getlCluster().hashCode() + 31 * (getrCluster().hashCode());
    } else {
      return rName + "~~~" + lName;//return getrCluster().hashCode() + 31 * (getlCluster().hashCode());
    }
  }

  @Override
  public String toString() {
    return data.toString();
  }
}
