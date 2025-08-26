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

/*
 * #%L
 * A Robust 3D Convex Hull Algorithm in Java
 * %%
 * Copyright (C) 2004 - 2014 John E. Lloyd
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

import ffx.utilities.FFXTest;
import org.junit.Assert;
import org.junit.Test;

public class ConvexHullTest extends FFXTest {

    @Test
    public void testConvexHull() throws Exception {
        Point3d[] points = new Point3d[]{
            new Point3d(0.3215426810286406, 0.1678336189760208, -0.2203710966001927),
            new Point3d(0.2229772524190855, -0.4213242506806965, -0.1966818060695024),
            new Point3d(0.3688830163971363, -0.1831502133823468, -0.2056387967482571),
            new Point3d(-0.1712592515826777, -0.3542439228428937, 0.2223876390814666),
            new Point3d(-0.3309556113844324, -0.370961861099081, 0.2439994981922204),
            new Point3d(-0.1004397059794885, -0.09014152417903909, -0.008600084584765189),
            new Point3d(0.458374538420117, -0.09914027349943322, -0.2505798421339875),
            new Point3d(-0.4954086979808367, -0.3339869997780649, -0.3195065691317492),
            new Point3d(0.053091190339151, 0.3036317017894533, 0.1380056861210668),
            new Point3d(0.4615616439483703, 0.4665423151725366, 0.1766835406205464),
            new Point3d(-0.4797380864431505, 0.0419809916447671, -0.4254776681079321),
            new Point3d(-0.003168473023146823, -0.2525299883005488, -0.27151530400991),
            new Point3d(-0.3577162826971303, -0.1375644040643837, -0.04494194644032229),
            new Point3d(-0.3392973838740004, 0.4288679723896719, -0.01599531622230571),
            new Point3d(0.1667164640191164, 0.003605551555385444, -0.4014989499947977),
            new Point3d(0.00714666676441833, 0.1140243407469469, 0.407090128778564),
            new Point3d(-0.03621271768232132, 0.3728502838619522, 0.4947140370446388),
            new Point3d(-0.3411871756810576, -0.3328629143842151, -0.4270033635450559),
            new Point3d(0.3544683273457627, -0.450828987127942, -0.0827870439577727),
            new Point3d(-0.4018510635028137, 0.08917494033386464, -0.2367824197158054),
            new Point3d(0.3978697768392692, -0.002667689232777493, 0.1641431727112673),
            new Point3d(-0.245701439441835, 0.495905311308713, -0.3194406286994373),
            new Point3d(0.161352035739787, -0.1563404972258401, 0.3852604361113724),
            new Point3d(0.07214279572678994, -0.4960366976410492, 0.1112227161519441),
            new Point3d(0.3201855824516951, 0.359077846965825, 0.02136723140381946),
            new Point3d(0.1190541238701475, -0.05734495917087884, 0.2032677509852384),
            new Point3d(0.3210853052521919, 0.4807189479290684, 0.4433501688235907),
            new Point3d(0.3862800354941562, 0.2085496142586224, 0.09336129957191763),
            new Point3d(0.1233572616459404, 0.265491605052251, 0.117400122450106),
            new Point3d(0.1438531872293476, -0.2594872752758556, -0.2026374435076839),
            new Point3d(0.2724846394476338, -0.3506708492996831, 0.2750346518820475),
            new Point3d(-0.4926118841325975, -0.3279366743079728, 0.3683135596740186),
            new Point3d(0.2459906458351674, 0.3647787136629026, -0.1641662355178652),
            new Point3d(-0.141922976953837, -0.2994764654892278, -0.3009570467294725),
            new Point3d(-0.1850859398814719, 0.2606059478228967, 0.004159106876849283),
            new Point3d(-0.09789466634196664, -0.3156603563722785, -0.303610991503681),
            new Point3d(0.2100642609503719, -0.4499717643018549, 0.3245569875692548),
            new Point3d(-0.1707163766685095, -0.2301452446078371, -0.05112823569320907),
            new Point3d(-0.312260808713977, -0.1674135249735914, 0.2808831662692904),
            new Point3d(-0.1966306233747216, 0.2291105671125563, -0.3387042454804333)
        };
        Point3d[] expected = {
            new Point3d(0.3215426810286406, 0.1678336189760208, -0.2203710966001927),
            new Point3d(0.2229772524190855, -0.4213242506806965, -0.1966818060695024),
            new Point3d(0.458374538420117, -0.09914027349943322, -0.2505798421339875),
            new Point3d(-0.4954086979808367, -0.3339869997780649, -0.3195065691317492),
            new Point3d(0.4615616439483703, 0.4665423151725366, 0.1766835406205464),
            new Point3d(-0.4797380864431505, 0.0419809916447671, -0.4254776681079321),
            new Point3d(-0.3392973838740004, 0.4288679723896719, -0.01599531622230571),
            new Point3d(0.1667164640191164, 0.003605551555385444, -0.4014989499947977),
            new Point3d(-0.03621271768232132, 0.3728502838619522, 0.4947140370446388),
            new Point3d(-0.3411871756810576, -0.3328629143842151, -0.4270033635450559),
            new Point3d(0.3544683273457627, -0.450828987127942, -0.0827870439577727),
            new Point3d(0.3978697768392692, -0.002667689232777493, 0.1641431727112673),
            new Point3d(-0.245701439441835, 0.495905311308713, -0.3194406286994373),
            new Point3d(0.161352035739787, -0.1563404972258401, 0.3852604361113724),
            new Point3d(0.07214279572678994, -0.4960366976410492, 0.1112227161519441),
            new Point3d(0.3210853052521919, 0.4807189479290684, 0.4433501688235907),
            new Point3d(0.2724846394476338, -0.3506708492996831, 0.2750346518820475),
            new Point3d(-0.4926118841325975, -0.3279366743079728, 0.3683135596740186),
            new Point3d(0.2459906458351674, 0.3647787136629026, -0.1641662355178652),
            new Point3d(0.2100642609503719, -0.4499717643018549, 0.3245569875692548)
        };

        QuickHull3D covexHull = new QuickHull3D();
        covexHull.build(points);

        logger.info("Vertices:");

        Point3d[] vertices = covexHull.getVertices();

        Assert.assertEquals(expected.length, vertices.length);
        for (int i = 0; i < vertices.length; i++) {
            Point3d pnt = vertices[i];
            logger.info(pnt.x + " " + pnt.y + " " + pnt.z);
            boolean found = false;
            for (Point3d point3d : expected) {
                double diff = Math.abs(pnt.x - point3d.x) + Math.abs(pnt.y - point3d.y) + Math.abs(pnt.y - point3d.y);
                if (diff < 0.0000001) {
                    found = true;
                }
            }
            Assert.assertTrue(found);
        }

    }
}
