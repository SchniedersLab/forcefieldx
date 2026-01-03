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
package ffx.numerics.clustering;

import ffx.utilities.FFXTest;
import org.junit.Test;

import static org.junit.Assert.*;

public class DistanceMapTest extends FFXTest {
    DistanceMap map = new DistanceMap();
    ClusterPair ab = new ClusterPair(new Cluster("a"), new Cluster("b"), 1.0);
    ClusterPair bc = new ClusterPair(new Cluster("b"), new Cluster("c"), 2.0);
    ClusterPair ca = new ClusterPair(new Cluster("c"), new Cluster("a"), 3.0);

    @Test
    public void testMapWorksWithSameDistance() throws Exception {
        this.map.add(ab);
        this.map.add(ab);  //add the same link twice. This seems to be an error case
        assertEquals(1,this.map.list().size());
        ClusterPair remove = this.map.removeFirst();
        assertNotNull(remove);
        assertEquals(0,this.map.list().size());  //still exists in the map(even though removeFirst will return null now)
        ClusterPair remove2 = this.map.removeFirst();
        assertNull(remove2);
    }
    @Test
    public void testMapRemovalFront() throws Exception {
        this.map.add(ca);
        this.map.add(bc);
        this.map.add(ab);

        ClusterPair removeFirst = this.map.removeFirst();
        assertEquals(ab, removeFirst);
    }
    @Test
    public void testMapRemovalByObjectPollLoop() throws Exception {
        this.map.add(ca);
        this.map.add(bc);
        this.map.add(ab);

        assertTrue(this.map.remove(ab)); //Doesn't actually remove from prioQueue
        ClusterPair removeFirst = this.map.removeFirst();
        assertEquals(bc, removeFirst);  //removeFirst should now skip the ab
    }
    @Test
    public void testMapRemovalByObjectPollLoopHandlesAllEmpty() throws Exception {
        this.map.add(ca);
        this.map.add(bc);
        this.map.add(ab);

        assertTrue(this.map.remove(ab)); //Doesn't actually remove from prioQueue
        assertTrue(this.map.remove(bc)); //Doesn't actually remove from prioQueue
        assertTrue(this.map.remove(ca)); //Doesn't actually remove from prioQueue
        assertFalse(this.map.remove(ab)); //Doesn't actually remove from prioQueue
        ClusterPair removeFirst = this.map.removeFirst();
        assertNull(removeFirst);  //removeFirst should now skip the ab
    }
}
