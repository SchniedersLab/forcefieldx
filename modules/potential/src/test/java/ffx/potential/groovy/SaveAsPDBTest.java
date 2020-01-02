/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.potential.groovy;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import ffx.potential.groovy.SaveAsPDB;
import ffx.utilities.DirectoryUtils;

import groovy.lang.Binding;

/**
 * Test the SaveAsPDB script.
 */
public class SaveAsPDBTest {

    Binding binding;
    SaveAsPDB saveAsPDB;

    @Before
    public void before() {
        binding = new Binding();
        saveAsPDB = new SaveAsPDB();
        saveAsPDB.setBinding(binding);
    }

    @After
    public void after() {
        saveAsPDB.destroyPotentials();
        System.gc();
    }

    @Test
    public void testSaveAsPDBHelp() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        saveAsPDB.run();
    }

    @Test
    public void testSaveAsPDB() {
        // Set-up the input arguments for the Biotype script.
        String[] args = {"src/main/java/ffx/potential/structures/peptide-amber99sb.xyz"};
        binding.setVariable("args", args);

        Path path = null;
        try {
            path = Files.createTempDirectory("SaveAsPDB");
            saveAsPDB.setBaseDir(path.toFile());
        } catch (IOException e) {
            Assert.fail(" Could not create a temporary directory.");
        }

        // Evaluate the script.
        saveAsPDB.run();

        // Delete all created space group directories.
        try {
            DirectoryUtils.deleteDirectoryTree(path);
        } catch (IOException e) {
            System.out.println(e.toString());
            Assert.fail(" Exception deleting files created by SaveAsPDB.");
        }
    }
}
