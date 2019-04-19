//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.xray;

import java.util.logging.Logger;

import ffx.algorithms.misc.PJDependentTest;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;


import ffx.realspace.cli.RealSpaceOptions;
import ffx.realspace.groovy.Alchemical;

import groovy.lang.Binding;

/**
 * Test the Energy script.
 */
public class AlchemicalTest extends PJDependentTest {

    private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

    Binding binding;
    Alchemical alchemical;

    @Before
    public void before() {
        binding = new Binding();
        alchemical = new Alchemical();
        alchemical.setBinding(binding);
    }

    @After
    public void after() {
        alchemical.destroyPotentials();
    }

    @Test
    public void testAlchemicalHelp() {
        // Set-up the input arguments for the Alchemical script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        alchemical.run();
    }

    @Test
    public void testAlchemical() {

        // Set-up the input arguments for the Alchemical script.
        String[] args = {"-N", "-n", "10", "-r", "0.01", "src/main/java/ffx/xray/structures/5zck.pdb",
                "src/main/java/ffx/xray/structures/5zck_ffx_2fofc.map"};
        binding.setVariable("args", args);

        alchemical.run();
    }
}
