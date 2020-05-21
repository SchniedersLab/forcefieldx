// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.potential.groovy;

import static java.lang.String.format;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import ffx.potential.utils.PotentialTest;
import ffx.utilities.DirectoryUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.MapConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.junit.After;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

@RunWith(Parameterized.class)
public class SolvatorTest extends PotentialTest {

  /**
   * Set of default options, such as "-d", "0.5" to specify a half-femtosecond timestep. Can be
   * over-ridden in test constructors.
   */
  private static final Map<String, String> DEFAULT_OPTIONS;

  static {
    String[] opts = {
        "-b", "2.5", // Standard 2.5 A bumpcheck.
        "-p", "7.0", // Use a shorter padding than ordinary to reduce compute cost.
    };

    int nOpts = opts.length;
    Map<String, String> optMap = new HashMap<>(nOpts);
    for (int i = 0; i < nOpts; i += 2) {
      optMap.put(opts[i], opts[i + 1]);
    }
    DEFAULT_OPTIONS = Collections.unmodifiableMap(optMap);
  }

  private final String info;
  private final SolvatorTestMode mode;
  private final File solvatedTestFile;
  private final Map<String, String> opts;
  /** Configuration containing the properties to be used by the solvator. */
  private final Configuration algorithmConfig;

  private final List<String> flags;
  private final File tempDir;
  private final List<File> copiedFiles = new ArrayList<>();
  Solvator solvator;

  public SolvatorTest(
      String info,
      SolvatorTestMode mode,
      String soluteFile,
      String solvatedTestFileName,
      String solventFileName,
      String ionFileName,
      String[] options,
      String[] properties,
      String[] flagArray)
      throws IOException {
    this.info = info;
    this.mode = mode;
    this.solvatedTestFile = new File("src/main/java/" + solvatedTestFileName);

    if (mode == SolvatorTestMode.HELP) {
      opts = Collections.emptyMap();
      algorithmConfig = null;
      flags = Collections.emptyList();
      tempDir = null;
    } else {
      File temp = null;
      try {
        temp = Files.createTempDirectory("Cart2Frac").toFile();
      } catch (java.io.IOException e) {
        fail(" Could not create a temporary directory.");
      }
      tempDir = temp;
      logger.fine(format(" Running test %s in directory %s", info, tempDir));
      String[] copiedExtensions = new String[] {"key", "properties", "ions"};
      String[] filesUsed;
      if (ionFileName == null) {
        filesUsed = new String[] {soluteFile, solventFileName};
      } else {
        filesUsed = new String[] {soluteFile, solventFileName, ionFileName};
      }
      for (String fname : filesUsed) {
        File srcFile = new File("src/main/java/" + fname);
        File tempFile = new File(tempDir.getAbsolutePath() + "/" + FilenameUtils.getName(fname));

        FileUtils.copyFile(srcFile, tempFile);
        copiedFiles.add(tempFile);
        logger.info(format(" Copied file %s to %s", srcFile, tempFile));

        for (String ext : copiedExtensions) {
          srcFile =
              new File(
                  format("%s.%s", FilenameUtils.removeExtension(srcFile.getPath()), ext));
          if (srcFile.exists()) {
            logger.fine(" Copying extension " + ext);
            tempFile =
                new File(
                    format("%s.%s", FilenameUtils.removeExtension(tempFile.getPath()), ext));
            logger.fine(format(" Copied file %s to %s", srcFile, tempFile));
            FileUtils.copyFile(srcFile, tempFile);
          }
        }
      }

      int nOpts = options.length;
      int nProps = properties.length;
      int nFlags = flagArray.length;

      if (nOpts > 0) {
        assertEquals(
            format("Unmatched option key %s for test %s", options[nOpts - 1], info),
            0,
            options.length % 2);
      }
      if (nProps > 0) {
        assertEquals(
            format("Unmatched property key %s for test %s", properties[nProps - 1], info),
            0,
            properties.length % 2);
      }

      Pattern validOption = Pattern.compile("^--?[^D]");
      Pattern validProperty = Pattern.compile("^[^-]");

      Map<String, String> groovyOpts = new HashMap<>(DEFAULT_OPTIONS);
      for (int i = 0; i < nOpts; i += 2) {
        String opti = options[i];
        assertTrue(
            format(" Option %s for test %s does not look like a Groovy option!", opti, info),
            validOption.matcher(opti).find());
        groovyOpts.put(opti, options[i + 1]);
      }
      groovyOpts.put("--sFi", copiedFiles.get(1).getPath());
      if (ionFileName != null) {
        groovyOpts.put("--iFi", copiedFiles.get(2).getPath());
      }
      this.opts = Collections.unmodifiableMap(groovyOpts);

      Map<String, String> addedProps = new HashMap<>();
      for (int i = 0; i < nProps; i += 2) {
        String propi = properties[i];
        assertTrue(
            format(" Property %s for test %s does not look like a property!", propi, info),
            validProperty.matcher(propi).find());
        addedProps.put(propi, properties[i + 1]);
      }
      algorithmConfig = new MapConfiguration(addedProps);

      Map<String, Boolean> addedFlags = new HashMap<>();
      for (int i = 0; i < nFlags; i += 2) {
        String flagi = flagArray[i];
        assertTrue(
            format(" Flag %s for test %s does not look like a flag!", flagi, info),
            validOption.matcher(flagi).find());
        String vali = flagArray[i + 1];
        assertTrue(
            format(
                " Value %s for flag %s in test %s is not a true/false value!", vali, flagi, info),
            vali.equalsIgnoreCase("TRUE") || vali.equalsIgnoreCase("FALSE"));
        addedFlags.put(flagi, Boolean.parseBoolean(vali));
      }
      this.flags =
          Collections.unmodifiableList(
              addedFlags.entrySet().stream()
                  .filter(Map.Entry::getValue)
                  .map(Map.Entry::getKey)
                  .collect(Collectors.toList()));
    }
  }

  @Parameterized.Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "Solvator Help Message Test",
                SolvatorTestMode.HELP,
                "",
                null,
                null,
                null,
                new String[] {},
                new String[] {},
                new String[] {"-h"}
            },
            {
                "Aspartate Pure Water Solvation",
                SolvatorTestMode.SOLVATE,
                "ffx/potential/structures/capAsp.xyz",
                "ffx/potential/structures/capAsp.pureWater.pdb",
                "ffx/potential/structures/watertiny.xyz",
                null,
                new String[] {"-b", "2.5", "-p", "8.0", "-s", "42"},
                new String[] {},
                new String[] {}
            },
            {
                "Aspartate Neutralizing NaCl (+200 mM) Solvation",
                SolvatorTestMode.SOLVATE,
                "ffx/potential/structures/capAsp.xyz",
                "ffx/potential/structures/capAsp.neutNaCl.pdb",
                "ffx/potential/structures/watertiny.xyz",
                "ffx/potential/structures/nacl.pdb",
                new String[] {"-b", "2.5", "-p", "8.0", "-s", "42"},
                new String[] {},
                new String[] {}
            },
            {
                "Aspartate Charged NaCl (200 mM) Solvation",
                SolvatorTestMode.SOLVATE,
                "ffx/potential/structures/capAsp.xyz",
                "ffx/potential/structures/capAsp.chargedNaCl.pdb",
                "ffx/potential/structures/watertiny.xyz",
                "ffx/potential/structures/naclCharged.pdb",
                new String[] {"-b", "2.5", "-p", "8.0", "-s", "42"},
                new String[] {},
                new String[] {}
            }
        });
  }

  @After
  public void after() {
    if (mode != SolvatorTestMode.HELP) {
      // Clean up the temporary directory if it exists.
      // Delete temporary directory.
      try {
        DirectoryUtils.deleteDirectoryTree(tempDir.toPath());
      } catch (IOException e) {
        System.out.println(e.toString());
        fail(" Exception deleting files created by Cart2Frac.");
      }
    }
  }

  @Test
  public void testSolvator() {
    switch (mode) {
      case HELP:
        testHelp();
        break;
      case SOLVATE:
        testSolvate();
        break;
    }
  }

  @Override
  public String toString() {
    return info;
  }

  private void testHelp() {
    String[] args = {"-h"};
    binding.setVariable("args", args);
    solvator = new Solvator(binding).run();
    potentialScript = solvator;
  }

  private void testSolvate() {
    // Prepare the Binding with input arguments.
    List<String> argList = new ArrayList<>(flags);
    opts.forEach(
        (String k, String v) -> {
          argList.add(k);
          argList.add(v);
        });
    argList.add(copiedFiles.get(0).getPath());
    String[] args = argList.toArray(new String[0]);
    binding.setVariable("args", args);
    solvator = new Solvator(binding);
    solvator.setProperties(algorithmConfig);

    solvator.run();
    potentialScript = solvator;

    File written = solvator.getWrittenFile();
    try (BufferedReader expectedReader = new BufferedReader(new FileReader(solvatedTestFile));
        BufferedReader writtenReader = new BufferedReader(new FileReader(written))) {
      boolean same = IOUtils.contentEqualsIgnoreEOL(expectedReader, writtenReader);
      assertTrue(" File written by test did not match the expected file!", same);
    } catch (IOException ex) {
      fail(format(" Exception %s in attempting to compare expected file %s to written file %s",
          ex.toString(), solvatedTestFile, written.toString()));
    }
  }

  private enum SolvatorTestMode {
    HELP,
    SOLVATE;
  }
}
