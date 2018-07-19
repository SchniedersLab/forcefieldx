package ffx.xray;
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

import edu.rit.pj.Comm;

import ffx.realspace.cli.RealSpaceOptions;
import ffx.realspace.groovy.Alchemical;
import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;

import groovy.lang.Binding;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Test the Energy script.
 */
public class AlchemicalTest {

    private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

    Binding binding;
    Alchemical alchemical;

    @Before
    public void before() {
        binding = new Binding();
        alchemical = new Alchemical();
        alchemical.setBinding(binding);
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

        // Initialize Parallel Java if needed
        try {
            Comm.world();
        } catch (IllegalStateException ise) {
            try {
                String args[] = new String[0];
                Comm.init(args);
            } catch (Exception e) {
                String message = String.format(" Exception starting up the Parallel Java communication layer.");
                logger.log(Level.WARNING, message, e.toString());
                message = String.format(" Skipping alchemical test.");
                logger.log(Level.WARNING, message, e.toString());
                return;
            }
        }
        // Set-up the input arguments for the Alchemical script.
        String[] args = {"-N", "/Users/sohali/Desktop/alchemicalTestCases/5zck/5zck.pdb",
                "/Users/sohali/Desktop/alchemicalTestCases/5zck/5zck_ffx_2fofc.map"};
        binding.setVariable("args", args);

        alchemical.run();
    }
}
