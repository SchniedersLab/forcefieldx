/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential.grooy;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import ffx.potential.groovy.Frac2Cart;

import groovy.lang.Binding;

/**
 * Test the Cart2Frac script.
 */
public class Frac2CartTest {

    Binding binding;
    Frac2Cart frac2Cart;

    @Before
    public void before() {
        binding = new Binding();
        frac2Cart = new Frac2Cart();
        frac2Cart.setBinding(binding);
    }

    @Test
    public void testFrac2CartHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        frac2Cart.run();

        // Pull out the biotype results to check.
        Assert.assertNull(frac2Cart.cartCoordinates);
        Assert.assertNull(frac2Cart.fracCoordinates);
    }

    @Test
    public void testFrac2Cart() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"src/main/java/ffx/potential/structures/acetanilide.frac.xyz"};
        binding.setVariable("args", args);

        // Evaluate the script.
        frac2Cart.writeFiles = false;
        frac2Cart.run();

        // Pull out the Cart2Frac results to check.
        double cartCoordinates[][] = frac2Cart.cartCoordinates;
        Assert.assertNotNull(cartCoordinates);
        Assert.assertEquals(19, cartCoordinates.length);
        Assert.assertEquals(7.98011035, cartCoordinates[0][0], 1.0e-6);
        Assert.assertEquals(0.70504091, cartCoordinates[0][1], 1.0e-6);
        Assert.assertEquals(0.99860734, cartCoordinates[0][2], 1.0e-6);

        double fracCoordinates[][] = frac2Cart.fracCoordinates;
        Assert.assertNotNull(fracCoordinates);
        Assert.assertEquals(19, fracCoordinates.length);
        Assert.assertEquals(0.4063192642, fracCoordinates[0][0], 1.0e-6);
        Assert.assertEquals(0.0743478761, fracCoordinates[0][1], 1.0e-6);
        Assert.assertEquals(0.1251544479, fracCoordinates[0][2], 1.0e-6);
    }
}
