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
package ffx.algorithms.groovy;

import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static org.apache.commons.io.FileUtils.copyFile;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import ffx.algorithms.misc.AlgorithmsTest;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import ffx.crystal.CrystalPotential;
import ffx.potential.bonded.LambdaInterface;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.MapConfiguration;
import org.apache.commons.io.FilenameUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

/**
 * Tests the functionality of the OST algorithm in Force Field X, both by running a very simple and
 * quick OST run, and setting up a number of systems, comparing output energies and gradients from
 * the starting algorithmConfig.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class ThermodynamicsTest extends AlgorithmsTest {

  // Enables behavior where it prints out instead of testing information.
  private static final boolean debugMode = false;
  /**
   * Set of default options, such as "-d", "0.5" to specify a half-femtosecond timestep. Can be
   * over-ridden in test constructors.
   */
  private static final Map<String, String> DEFAULT_OPTIONS;
  /**
   * Set of default properties, such as "ost-alwaystemper", "true" to use standard tempering scheme.
   * Can be over-ridden in test constructors.
   */
  private static final Map<String, String> DEFAULT_PROPERTIES;
  /**
   * Set of default boolean options, such as "-y", false to specify "do not add the -y --synchronous
   * flag". Can be over-ridden in test constructors. Maps a flag to whether it should be included by
   * default.
   */
  private static final Map<String, Boolean> DEFAULT_FLAGS;
  /**
   * Default number of atomic gradients to test; pick the first N (3) if ffx.ci is false. Otherwise,
   * use the whole array to test.
   */
  private static final int DEFAULT_GRADIENT_EVALS = 3;
  // Currently, everything is just 1.0E-2 kcal/mol per unit.
  private static final double DEFAULT_PE_TOL = 1.0E-2; // kcal/mol.
  private static final double DEFAULT_DUDL_TOL = 1.0E-2; // per lambda (unit-less).
  private static final double DEFAULT_D2UDL2_TOL = 1.0E-2; // per lambda^2 (unit-less).
  private static final double DEFAULT_DUDX_TOL = 1.0E-2; // kcal/mol/A.
  private static final double DEFAULT_D2UDXDL_TOL = 1.0E-2; // kcal/mol/A.

  static {
    String[] opts = {
        "--bM",
        "0.05",
        "--dw",
        "OFF",
        "--lf",
        "1.0E-18",
        "--lm",
        "1.0E-18",
        "--np",
        "1",
        "--sf",
        "1.0",
        "--tp",
        "2.0", // Set low to encourage fast tempering for toy cases.
        "-b",
        "ADIABATIC", // Used for Langevin dynamics.
        "-C",
        "5000", // By default, evaluate 0 dynamics steps, and trigger adding counts manually.
        "-d",
        "0.5", // Most of the tests are in vacuum, thus short timestep.
        "-F",
        "XYZ",
        "-i",
        "STOCHASTIC", // Stochastic because most toy systems work best with Langevin dynamics.
        "-k",
        "5.0", // Infrequent because we should not need to restart.
        "-l",
        "1.0", // Set walker to start at the usually-physical topology. Likely commonly set by test.
        "-n",
        "0", // Do not run any steps by default, just evaluate energy.
        "-p",
        "0", // Use NVT, not NPT dynamics.
        "-Q",
        "0", // Do not run any steps by default, just evaluate energy.
        "-r",
        "2.0", // Very infrequent printouts, since it should work.
        "-t",
        "298.15",
        "-w",
        "10.0"
    }; // Infrequent since we should not need the trajectory.

    int nOpts = opts.length;
    Map<String, String> optMap = new HashMap<>(nOpts);
    for (int i = 0; i < nOpts; i += 2) {
      optMap.put(opts[i], opts[i + 1]);
    }
    DEFAULT_OPTIONS = Collections.unmodifiableMap(optMap);

    String[] props = {
        "ost-alwaystemper",
        "true",
        "ost-temperOffset",
        "0.5", // Small to temper fast.
        "print-on-failure",
        "false",
        "disable-neighbor-updates",
        "false",
        "vdw-lambda-alpha",
        "0.25",
        "vdw-lambda-exponent",
        "3",
        "permanent-lambda-alpha",
        "1.0",
        "permanent-lambda-exponent",
        "3"
    };

    int nProps = props.length;
    Map<String, String> propMap = new HashMap<>(nProps);
    for (int i = 0; i < nProps; i += 2) {
      propMap.put(props[i], props[i + 1]);
    }
    DEFAULT_PROPERTIES = Collections.unmodifiableMap(propMap);

    DEFAULT_FLAGS = Map.of("--rn", false, "-y", false, "-h", false);
  }

  private final String info;
  private final ThermoTestMode mode;
  private final boolean doTest;
  private final File[] copiedFiles;
  private final Map<String, String> opts;
  private final List<String> flags;
  /**
   * Free energy associated with this test, if calculated. Not used alongside testing gradients, and
   * values are meaningless in that case.
   */
  private final double freeEnergy;
  private final double feTol;
  private final int numGradAtoms;
  /** Potential energy of the underlying Potential (0), and OST after one bias drop (1). */
  private final double[] pe;
  private final double[] dudl;
  private final double[] d2udl2;
  /**
   * Gradient of the underlying Potential ([0][0..n][0-2]), and OST after one bias drop
   * ([1][0..n][0-2]). Axes: before/after bias, atoms, X/Y/Z.
   */
  private final double[][][] dudx;
  private final double[][][] d2udxdl;
  // TODO: Make these settable if needed.
  private final double peTol = DEFAULT_PE_TOL;
  private final double dudlTol = DEFAULT_DUDL_TOL;
  private final double d2udl2Tol = DEFAULT_D2UDL2_TOL;
  private final double dudxTol = DEFAULT_DUDX_TOL;
  private final double d2udxdlTol = DEFAULT_D2UDXDL_TOL;
  /** Configuration containing the properties to be used by OST. */
  Configuration algorithmConfig;

  public ThermodynamicsTest(
      String info,
      String[] filenames,
      ThermoTestMode mode,
      boolean ciOnly,
      double freeEnergy,
      double feTol,
      int[] gradAtomIndices,
      double[] pe,
      double[] dudl,
      double[] d2udl2,
      double[][][] dudx,
      double[][][] d2udxdl,
      String[] options,
      String[] properties,
      String[] flags)
      throws IOException {
    this.info = info;
    this.mode = mode;
    doTest = (ffxCI || !ciOnly);

    File tempDir;
    if (doTest) {
      int nFiles = filenames.length;
      tempDir = registerTemporaryDirectory().toFile();
      String tempDirName = tempDir.getAbsolutePath() + File.separator;
      logger.fine(format(" Running test %s in directory %s", info, tempDirName));
      copiedFiles = new File[nFiles];

      String[] copiedExtensions = new String[] {"dyn", "key", "properties", "his", "lam", "prm"};
      for (int i = 0; i < nFiles; i++) {
        File srcFile = new File("src/main/java/" + filenames[i]);
        File tempFile = new File(tempDirName + FilenameUtils.getName(filenames[i]));
        copyFile(srcFile, tempFile);
        logger.fine(format(" Copied file %s to %s", srcFile, tempFile));
        copiedFiles[i] = tempFile;

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
            copyFile(srcFile, tempFile);
          }
        }
      }

      switch (mode) {
        case HELP:
          assertTrue(
              format(
                  "Help tests must have no file arguments, found %d for %s", nFiles, info),
              nFiles == 0);
          break;
        default:
          assertTrue(
              format(
                  "Must have 1, 2, or 4 distinct filenames, found %d for test %s", nFiles, info),
              nFiles == 1 || nFiles == 2 || nFiles == 4);
          for (int i = 0; i < nFiles; i++) {
            for (int j = i + 1; j < nFiles; j++) {
              assertNotEquals(
                  format(
                      " Filenames %d and %d matched in test %s: files %s and %s",
                      i, j, info, filenames[i], filenames[j]),
                  filenames[i],
                  filenames[j]);
            }
          }
          break;
      }
      int nOpts = options.length;
      int nProps = properties.length;
      int nFlags = flags.length;

      if (nOpts > 0) {
        assertTrue(
            format("Unmatched option key %s for test %s", options[nOpts - 1], info),
            options.length % 2 == 0);
      }
      if (nProps > 0) {
        assertTrue(
            format("Unmatched property key %s for test %s", properties[nProps - 1], info),
            properties.length % 2 == 0);
      }
      if (nFlags > 0) {
        assertTrue(
            format("Unmatched flag key %s for test %s", flags[nFlags - 1], info),
            flags.length % 2 == 0);
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
      this.opts = Collections.unmodifiableMap(groovyOpts);

      Map<String, String> addedProps = new HashMap<>(DEFAULT_PROPERTIES);
      for (int i = 0; i < nProps; i += 2) {
        String propi = properties[i];
        assertTrue(
            format(" Property %s for test %s does not look like a property!", propi, info),
            validProperty.matcher(propi).find());
        addedProps.put(propi, properties[i + 1]);
      }
      algorithmConfig = new MapConfiguration(addedProps);

      Map<String, Boolean> addedFlags = new HashMap<>(DEFAULT_FLAGS);
      for (int i = 0; i < nFlags; i += 2) {
        String flagi = flags[i];
        assertTrue(
            format(" Flag %s for test %s does not look like a flag!", flagi, info),
            validOption.matcher(flagi).find());
        String vali = flags[i + 1];
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

      // Only meaningful for free energy evaluations.
      this.freeEnergy = freeEnergy;
      this.feTol = feTol;
      if (mode == ThermoTestMode.FREE) {
        assertTrue(
            format(" Free energy tolerance for test %s was %10.4g <= 0!", info, feTol),
            feTol > 0);
      }

      // Only meaningful for gradient evaluations.
      this.dudx = new double[2][][];
      this.d2udxdl = new double[2][][];
      boolean isGradient = (mode == ThermoTestMode.GRAD);
      if (isGradient) {
        numGradAtoms =
            ffxCI
                ? gradAtomIndices.length
                : Math.min(DEFAULT_GRADIENT_EVALS, gradAtomIndices.length);

        if (debugMode) {
          this.pe = this.dudl = this.d2udl2 = null;
        } else {
          assertNotNull(gradAtomIndices);
          assertNotNull(pe);
          assertNotNull(dudl);
          assertNotNull(d2udl2);

          this.pe = copyOf(pe, 2);
          this.dudl = copyOf(dudl, 2);
          this.d2udl2 = copyOf(d2udl2, 2);
          for (int i = 0; i < 2; i++) {
            assertNotNull(dudx[i]);
            assertNotNull(d2udxdl[i]);
            this.dudx[i] = new double[numGradAtoms][3];
            this.d2udxdl[i] = new double[numGradAtoms][3];

            for (int j = 0; j < numGradAtoms; j++) {
              assertNotNull(dudx[i][j]);
              assertNotNull(d2udxdl[i][j]);
              arraycopy(dudx[i][j], 0, this.dudx[i][j], 0, 3);
              arraycopy(d2udxdl[i][j], 0, this.d2udxdl[i][j], 0, 3);
            }
          }
        }
      } else {
        // this.gradAtomIndices = null;
        numGradAtoms = 0;
        this.pe = null;
        this.dudl = null;
        this.d2udl2 = null;
        for (int i = 0; i < 2; i++) {
          this.dudx[i] = null;
          this.d2udxdl[i] = null;
        }
      }
    } else {
      logger.fine(" Skipping test " + info);
      tempDir = null;
      copiedFiles = new File[0];
      opts = Collections.emptyMap();
      this.flags = Collections.emptyList();
      this.freeEnergy = 0;
      this.feTol = 0;
      numGradAtoms = 0;
      this.pe = new double[0];
      this.dudl = new double[0];
      this.d2udl2 = new double[0];
      this.dudx = new double[0][0][0];
      this.d2udxdl = new double[0][0][0];
    }
  }

  @Parameterized.Parameters
  public static Collection<Object[]> data() {
    /*
     Test info, filenames, mode, dG, tol(dG), grad atoms, PE, dU/dL, d2U/dL2,
     dU/dX, d2U/dXdL, Groovy options, properties, Groovy flags
    */
    return Arrays.asList(
        new Object[][] {
            {
                "Thermodynamics Help Message Test",
                new String[] {},
                ThermoTestMode.HELP,
                false,
                0,
                0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {},
                new String[] {},
                new String[] {"-h", "true"}
            },
            {
                "Acetamide Implicit Solvation Free Energy: -10.5 kcal/mol",
                new String[] {"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.FREE,
                false,
                -9.2,
                1.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-C", "10", "--ac", "1-9", "-d", "1.0", "-n", "20000", "-w", "5", "--bM", "0.25",
                    "--tp", "2.0"
                },
                new String[] {"randomseed", "42"},
                new String[] {}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the
                // ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 0.0",
                new String[] {
                    "ffx/algorithms/structures/4icb_ca_a.xyz",
                    "ffx/algorithms/structures/4icb_ca_b.xyz",
                    "ffx/algorithms/structures/4icb_mg_a.xyz",
                    "ffx/algorithms/structures/4icb_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                new int[] {1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-8422.43552052, 0},
                new double[] {0, 0},
                new double[] {-187.708505969645246, Double.NaN},
                new double[][][] {
                    {
                        {1.51602997670091, -0.16798107175016064, 1.1909011229485826},
                        {-1.2287108493416157, 0.3317143880292477, -1.168707082997385},
                        {1.330781894105006, -0.8783264761133562, 0.1096093788802257},
                        {-1.0136367827497708, 0.737243265172173, -0.5962871564399377},
                        {-0.17816595179992722, -0.35265175625985923, -0.6111655937237352},
                        {-0.25917771149474067, 0.07021711090701199, 0.15836021218566376},
                        {-0.5324396799274522, 0.7187373577290392, -0.6184469870748088},
                        {0.2993938879707754, -0.3907323336942321, 0.4139579476752431},
                        {0.054758010534996515, 0.3356578298638757, 0.9215488358561414},
                        {-0.882274050563856, 0.5899806670077385, -0.9015065344963205},
                        {-0.09630989055561034, 0.9691000112562431, -0.2005707212562129},
                        {0.167130858362869, -0.62725014243849, 1.175516097449706}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {-2.9163424369832507E-16, 3.2314026482573344E-17, -2.2909015893365526E-16},
                        {2.3636350519364876E-16, -6.38109246937693E-17, 2.248207565107174E-16},
                        {-2.5599861294253615E-16, 1.6896109016191808E-16, -2.1085234990896332E-17},
                        {1.9499033730541748E-16, -1.4182132633554818E-16, 1.1470601278862118E-16},
                        {3.427326200965616E-17, 6.783858486066127E-17, 1.1756813416574204E-16},
                        {4.985725680683067E-17, -1.3507459844962662E-17, -3.0463289923312223E-17},
                        {1.0242386084510416E-16, -1.3826139915463767E-16, 1.1896883446563482E-16},
                        {-5.75935248168813E-17, 7.516403394171084E-17, -7.963187723760774E-17},
                        {-1.0533638011269165E-17, -6.456951304273904E-17, -1.7727564883699187E-16},
                        {1.6972047349009049E-16, -1.1349285189852184E-16, 1.7342017006092682E-16},
                        {1.8526851397737774E-17, -1.8642292909391106E-16, 3.8583201849920577E-17},
                        {-3.215047342492132E-17, 1.2066227166417432E-16, -2.2613058666622065E-16}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "0.0",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "329-345.857-972.1008-1022.1204.1208-1213",
                    "--uaB",
                    "329-345.857-972.1008-1022.1204.1208-1213"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the
                // ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 0.5",
                new String[] {
                    "ffx/algorithms/structures/4icb_ca_a.xyz",
                    "ffx/algorithms/structures/4icb_ca_b.xyz",
                    "ffx/algorithms/structures/4icb_mg_a.xyz",
                    "ffx/algorithms/structures/4icb_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                false,
                0,
                0,
                new int[] {1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-8441.45436853, 0},
                new double[] {-59.7520122, 0},
                new double[] {0, Double.NaN},
                new double[][][] {
                    {
                        {1.516028943632357, -0.16798187601812842, 1.1909000012817135},
                        {-1.2287100421217505, 0.33171692634032457, -1.1687028725533608},
                        {1.3307953311709682, -0.8783244163319, 0.10962280681203751},
                        {-1.0136504919823237, 0.7372413665734416, -0.5963003206269726},
                        {-0.17816209958414242, -0.3526516261736795, -0.6111638455124022},
                        {-0.25917449280899385, 0.07021855608468686, 0.1583629255853478},
                        {-0.5324365378704783, 0.718739459636424, -0.6184448359875749},
                        {0.2993955300471849, -0.3907325347437327, 0.4139595896177193},
                        {0.05475347093067029, 0.335655947922745, 0.9215447577884727},
                        {-0.88227515311464, 0.5899832887545906, -0.901509658156534},
                        {-0.09630970969676556, 0.9690984223562502, -0.20056743541628919},
                        {0.16711991894934464, -0.6272519054249652, 1.1755119757400072}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {-3.2454805314330315E-6, -2.5266823584502163E-6, -3.5238203937026924E-6},
                        {2.535955989202421E-6, 7.974339453653556E-6, 1.3227500003987203E-5},
                        {4.221378770274953E-5, 6.470994293295007E-6, 4.218509194231501E-5},
                        {-4.306882429183645E-5, -5.964623852605655E-6, -4.135651328152079E-5},
                        {1.2102092812327214E-5, 4.0867779382836034E-7, 5.492167883147658E-6},
                        {1.0111799477741101E-5, 4.540159575405767E-6, 8.524396519327126E-6},
                        {9.871063099353705E-6, 6.603336803578941E-6, 6.757839859261594E-6},
                        {5.158735174148887E-6, -6.316156557772956E-7, 5.1583144315969776E-6},
                        {-1.4261587612196536E-5, -5.912292426302201E-6, -1.2811627438935602E-5},
                        {-3.463765445133049E-6, 8.236460654842404E-6, -9.81326797955262E-6},
                        {5.681848183058946E-7, -4.991676597398964E-6, 1.0322770563675476E-5},
                        {-3.436718115423787E-5, -5.538585353903613E-6, -1.2948732912576588E-5}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "0.5",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "329-345.857-972.1008-1022.1204.1208-1213",
                    "--uaB",
                    "329-345.857-972.1008-1022.1204.1208-1213"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the
                // ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 1.0",
                new String[] {
                    "ffx/algorithms/structures/4icb_ca_a.xyz",
                    "ffx/algorithms/structures/4icb_ca_b.xyz",
                    "ffx/algorithms/structures/4icb_mg_a.xyz",
                    "ffx/algorithms/structures/4icb_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                new int[] {1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-8460.47321653, 0},
                new double[] {0, 0},
                new double[] {187.70850596942546, Double.NaN},
                new double[][][] {
                    {
                        {1.5160279105638157, -0.16798268028609797, 1.1908988796148456},
                        {-1.2287092349018947, 0.3317194646514121, -1.1686986621093438},
                        {1.3308087682369107, -0.8783223565504423, 0.10963623474384576},
                        {-1.0136642012148873, 0.7372394679746987, -0.5963134848140097},
                        {-0.1781582473683545, -0.35265149608749957, -0.6111620973010643},
                        {-0.2591712741232506, 0.07022000126236438, 0.15836563898503053},
                        {-0.5324333958135048, 0.718741561543808, -0.6184426849003377},
                        {0.29939717212359174, -0.3907327357932422, 0.4139612315602017},
                        {0.05474893132634229, 0.3356540659816183, 0.921540679720797},
                        {-0.8822762556654258, 0.5899859105014462, -0.9015127818167441},
                        {-0.09630952883792121, 0.9690968334562342, -0.20056414957637436},
                        {0.1671089795358247, -0.6272536684114378, 1.1755078540303066}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {2.9163384624158903E-16, -3.231433591192038E-17, 2.290897273901192E-16},
                        {-2.3636319462860973E-16, 6.381190126870034E-17, -2.248191366091609E-16},
                        {2.560037826405337E-16, -1.6896029769367465E-16, 2.109040117467388E-17},
                        {-1.9499561171519864E-16, 1.4182059587979724E-16, -1.147110775007824E-16},
                        {-3.4271779930734056E-17, -6.783853481206671E-17, -1.1756746156916025E-16},
                        {-4.9856018468543316E-17, 1.3508015854151776E-17, 3.0464333860802826E-17},
                        {-1.0242265198852129E-16, 1.38262207830166E-16, -1.1896800686893987E-16},
                        {5.759415657973307E-17, -7.51641112923199E-17, 7.963250894893407E-17},
                        {1.0531891470506986E-17, 6.456878899573994E-17, 1.7727407986513553E-16},
                        {-1.6972089767901735E-16, 1.1349386057403977E-16, -1.734213718396476E-16},
                        {-1.8526781815166296E-17, 1.8642231778983455E-16, -3.858193767513045E-17},
                        {3.214626465908195E-17, -1.2066294994525673E-16, 2.2612900090378923E-16}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "1.0",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "329-345.857-972.1008-1022.1204.1208-1213",
                    "--uaB",
                    "329-345.857-972.1008-1022.1204.1208-1213"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the
                // ions, and some water.
                "Carp Parvalbumin Ca/Mg Simultaneous Bookending: L = 0.0",
                new String[] {
                    "ffx/algorithms/structures/5cpv_ca_a.xyz",
                    "ffx/algorithms/structures/5cpv_ca_b.xyz",
                    "ffx/algorithms/structures/5cpv_mg_a.xyz",
                    "ffx/algorithms/structures/5cpv_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                new int[] {1, 303, 1401, 1402, 1482, 1488, 1489, 1602, 1603, 1604, 1605, 1606},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-10378.104156403408, 0},
                new double[] {0, 0},
                new double[] {-197.40065865800716, Double.NaN},
                new double[][][] {
                    {
                        {-1.1610229302444828, 0.8877842196621852, -0.38190472079527416},
                        {0.17141739386080523, -0.6936874602017804, -0.2646850671524603},
                        {1.7845002387380973, 0.06252966517014968, 0.6441113717050584},
                        {-0.8091327019583163, -0.03278919597876184, -0.6934281103867814},
                        {-0.009204859689755152, -0.03031752596955206, -0.2231311840527026},
                        {0.44438169596232324, -0.03623403581690088, -0.09120393985608466},
                        {0.6143253257220098, -0.559457100728987, 0.7795004807627701},
                        {-0.4100603133810665, 0.27711514086384526, -0.3596263242692279},
                        {0.2060111442204393, -1.6995161698485726, 0.316094377538672},
                        {0.17044049656085825, 1.3626089754787625, 0.18294188325894245},
                        {-0.1369658985076594, 0.28860743562734736, -0.23182335542775423},
                        {-0.3317851447377951, 0.10693619585058839, 0.07599282491286408}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {2.2334257856503093E-16, -1.707804485712817E-16, 7.346589191878714E-17},
                        {-3.2975061696419736E-17, 1.334426237792466E-16, 5.09166906746972E-17},
                        {-3.4327908122000154E-16, -1.202864843760528E-17, -1.2390581692419626E-16},
                        {1.5565048660891405E-16, 6.307561537504449E-18, 1.33392733415393E-16},
                        {1.7707118825002895E-18, 5.832093621377665E-18, 4.292309196176272E-17},
                        {-8.548440453491081E-17, 6.970235281623787E-18, 1.754463462531147E-17},
                        {-1.181759625502523E-16, 1.0762112290669118E-16, -1.4995022305853868E-16},
                        {7.888210075095321E-17, -5.330782752663655E-17, 6.918026206876034E-17},
                        {-3.962976007169519E-17, 3.269309449444003E-16, -6.080614905213747E-17},
                        {-3.278713882575124E-17, -2.62121095313077E-16, -3.519199394795007E-17},
                        {2.6347728499845737E-17, -5.55185665906526E-17, 4.459517949566215E-17},
                        {6.382453595443383E-17, -2.0571002605587047E-17, -1.4618517021800703E-17}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "0.0",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "810-829.1339-1447.1476-1490.1603.1605-1610",
                    "--uaB",
                    "810-829.1339-1447.1476-1490.1603.1605-1610"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                "Carp Parvalbumin Ca/Mg Simultaneous Bookending: L = 0.5",
                new String[] {
                    "ffx/algorithms/structures/5cpv_ca_a.xyz",
                    "ffx/algorithms/structures/5cpv_ca_b.xyz",
                    "ffx/algorithms/structures/5cpv_mg_a.xyz",
                    "ffx/algorithms/structures/5cpv_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                new int[] {1, 303, 1401, 1402, 1482, 1488, 1489, 1602, 1603, 1604, 1605, 1606},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-10398.105024790417, 0},
                new double[] {-62.8370972544, 0},
                new double[] {0, Double.NaN},
                new double[][][] {
                    {
                        {-1.1610220532282816, 0.8877825087765396, -0.38190561321870153},
                        {0.17141666802985345, -0.6936868435232721, -0.26468548230330735},
                        {1.7844987915596682, 0.0625323990316684, 0.6441061559020982},
                        {-0.8091313886502285, -0.032792749245311725, -0.6934234392602314},
                        {-0.009206939878777853, -0.030315101956861712, -0.2231346110988761},
                        {0.4443810630713887, -0.036232529267353986, -0.09120578149088898},
                        {0.6143248825840462, -0.5594553644984039, 0.7794994337695262},
                        {-0.4100603484233871, 0.2771150828568485, -0.3596269108370689},
                        {0.2060121473201555, -1.6995185434801163, 0.31609655101286327},
                        {0.1704401628733358, 1.3626093888666202, 0.18294114185436605},
                        {-0.13696614734124157, 0.2886080953858223, -0.23182359983306977},
                        {-0.33178545555019134, 0.10693683664021286, 0.07599211616273438}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {2.755227653583603E-6, -5.374905775568095E-6, -2.803630881231811E-6},
                        {-2.2802651855613476E-6, 1.9373526711774502E-6, -1.3042348507852353E-6},
                        {-4.54644512082325E-6, 8.588679262899745E-6, -1.638592826225249E-5},
                        {4.125879040195457E-6, -1.1162916089269004E-5, 1.46747768532407E-5},
                        {-6.535106551908143E-6, 7.6152604604473595E-6, -1.0766383081950437E-5},
                        {-1.988285510456933E-6, 4.732964988818367E-6, -5.78566637177147E-6},
                        {-1.3921589704368742E-6, 5.4545292442753635E-6, -3.289226282676694E-6},
                        {-1.1008869904571839E-7, -1.8223435560571488E-7, -1.842757221481861E-6},
                        {3.151330700390531E-6, -7.456983421860741E-6, 6.8281705551953564E-6},
                        {-1.0483102688141344E-6, 1.2986962580896488E-6, -2.329191170602485E-6},
                        {-7.817337568383209E-7, 2.0726923768421557E-6, -7.678219411388909E-7},
                        {-9.764459392158642E-7, 2.0130999769385483E-6, -2.226604200572524E-6}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "0.5",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "810-829.1339-1447.1476-1490.1603.1605-1610",
                    "--uaB",
                    "810-829.1339-1447.1476-1490.1603.1605-1610"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                "Carp Parvalbumin Ca/Mg Simultaneous Bookending: L = 1.0",
                new String[] {
                    "ffx/algorithms/structures/5cpv_ca_a.xyz",
                    "ffx/algorithms/structures/5cpv_ca_b.xyz",
                    "ffx/algorithms/structures/5cpv_mg_a.xyz",
                    "ffx/algorithms/structures/5cpv_mg_b.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                new int[] {1, 303, 1401, 1402, 1482, 1488, 1489, 1602, 1603, 1604, 1605, 1606},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-10418.105893177431, 0},
                new double[] {0, 0},
                new double[] {197.40065865800716, Double.NaN},
                new double[][][] {
                    {
                        {-1.1610211762120812, 0.887780797890894, -0.3819065056421289},
                        {0.17141594219890166, -0.6936862268447639, -0.2646858974541544},
                        {1.78449734438124, 0.06253513289318713, 0.644100940099138},
                        {-0.8091300753421411, -0.03279630251186161, -0.6934187681336814},
                        {-0.009209020067800555, -0.030312677944171362, -0.22313803814504962},
                        {0.44438043018045414, -0.03623102271780709, -0.0912076231256933},
                        {0.614324439446083, -0.5594536282678209, 0.7794983867762828},
                        {-0.4100603834657077, 0.27711502484985173, -0.35962749740490985},
                        {0.20601315041987167, -1.699520917111661, 0.31609872448705456},
                        {0.17043982918581357, 1.362609802254478, 0.18294040044978965},
                        {-0.13696639617482465, 0.2886087551442973, -0.2318238442383853},
                        {-0.3317857663625876, 0.10693747742983745, 0.07599140741260468}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new double[][][] {
                    {
                        {-2.2334224114695824E-16, 1.7077979033516612E-16, -7.346623526454586E-17},
                        {3.297478244447361E-17, -1.3344238652197185E-16, -5.091685039740076E-17},
                        {3.4327852444105517E-16, 1.202970024746214E-17, 1.2390381022673677E-16},
                        {-1.5564998133445806E-16, -6.308928600450258E-18, -1.3339093627354278E-16},
                        {-1.7715122022323975E-18, -5.831161020942922E-18, -4.292441046342064E-17},
                        {8.548416104016219E-17, -6.969655660581429E-18, -1.7545343165091794E-17},
                        {1.181757920599496E-16, -1.0762045491951309E-16, 1.499498202444948E-16},
                        {-7.888211423293022E-17, 5.3307805209364545E-17, -6.918048774143343E-17},
                        {3.963014599840059E-17, -3.2693185816148815E-16, 6.080698526185848E-17},
                        {3.2787010444769745E-17, 2.621212543574982E-16, 3.5191708704298895E-17},
                        {-2.6347824234619856E-17, 5.551882042226125E-17, -4.459527352673071E-17},
                        {-6.382465553457344E-17, 2.057124913923137E-17, 1.4618244341429974E-17}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[12][3]
                },
                new String[] {
                    "-l",
                    "1.0",
                    "--sf",
                    "TRIG",
                    "--uaA",
                    "810-829.1339-1447.1476-1490.1603.1605-1610",
                    "--uaB",
                    "810-829.1339-1447.1476-1490.1603.1605-1610"
                },
                new String[] {"disable-neighbor-updates", "true"},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.0",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {2.61926128937, 0},
                new double[] {17.7107919457, 0},
                new double[] {-403.664911547, Double.NaN},
                new double[][][] {
                    {
                        {-5.677212457329941, 1.4246439517642553, -2.2981580040498444},
                        {-16.262337142221504, -4.01152481099988, 2.8309493281407083},
                        {20.234575856650853, -2.575942513845196, 0.3415783090214405},
                        {-4.107688955265742E-5, -3.0491094574426836E-4, 1.5088773834574196E-4},
                        {3.855038222953342E-5, 3.0678248144100465E-4, -1.5278355755630972E-4},
                        {2.5265073231240012E-6, -1.8715356967362725E-6, 1.895819210567746E-6},
                        {1.7049737429005947, 5.162823373080822, -0.8743696331123039}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {34.717862394027094, -6.918833553692483, 13.770559137316184},
                        {44.6501751603733, 28.43111292615495, -17.0699023179925},
                        {-125.31576661843137, 11.778366126287878, -3.1016651209654698},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {45.94772906403098, -33.29064549875034, 6.401008301641785}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.0", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.1",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {2.850550906, 0},
                new double[] {-8.71629068942, 0},
                new double[] {-141.329823284, Double.NaN},
                new double[][][] {
                    {
                        {-2.9145419263022183, 0.8087080881505835, -1.1826969161418945},
                        {-10.75486342640376, -1.8502906018217953, 1.441002156956333},
                        {10.24415407077068, -1.5066777027112754, 0.13037812748160396},
                        {-3.940188541217106E-5, -3.046447158662098E-4, 1.504274409501606E-4},
                        {4.348535004930246E-5, 3.0537557756647634E-4, -1.5274632518595476E-4},
                        {5.712456879867664E-6, -9.857865557163607E-7, 1.145744081896436E-6},
                        {3.425241486013779, 2.548260471307343, -0.38868219515588814}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {21.219764941164165, -5.301043079657732, 8.702685150224386},
                        {59.404183336652935, 15.559560209446712, -10.86515770297192},
                        {-76.64822167834306, 9.426371612748719, -1.3109489861942616},
                        {8.074643396228061E-5, 1.2233389408248271E-5, -2.0734472804670857E-5},
                        {2.2173087875276693E-4, -6.414406402810961E-5, 2.6740754942568507E-6},
                        {1.7176239889004114E-4, 3.3106310006263015E-5, -3.056405648898283E-5},
                        {-3.9762008391856365, -19.684869938173087, 3.4734701633955916}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.1", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.25",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                false,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {0.974157909812, 0},
                new double[] {-10.6503844408, 0},
                new double[] {85.6111397243, Double.NaN},
                new double[][][] {
                    {
                        {-0.8087917281519832, 0.230920162100479, -0.31677660566831417},
                        {-3.2828492735876607, -0.4439185686570825, 0.3694313047851273},
                        {2.7179567720851567, -0.4453048785572774, 0.02913967164983769},
                        {1.1191654242300788E-4, -2.8436948793288743E-4, 1.1810143809370248E-4},
                        {3.948395026959898E-4, 1.992169489195181E-4, -1.434611972003E-4},
                        {4.178187790317917E-4, 1.7909105729778537E-5, -3.118731586794815E-5},
                        {1.3727596548303362, 0.6583705285471645, -0.08173782369167623}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {8.203573000276744, -2.407260166945433, 3.3325694163231305},
                        {34.07610421632959, 4.853369454929932, -4.037664708947347},
                        {-28.691594190062066, 4.581997199828557, -0.30050567027187414},
                        {0.0030722784633713677, 4.013278190369727E-4, -6.302447025815397E-4},
                        {0.006927249436665227, -0.0021139142799777464, 2.0753922987915073E-4},
                        {0.008777938573427544, 2.1069797116961348E-4, -5.578096868382651E-4},
                        {-13.606860493017741, -7.026604599323282, 1.006581478055631}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.25", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.4",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {0.110037916292, 0},
                new double[] {-2.26607999096, 0},
                new double[] {30.1632932351, Double.NaN},
                new double[][][] {
                    {
                        {-0.1135918356399739, 0.031116766115897736, -0.041752199556192734},
                        {-0.4844299250151845, -0.052731371080211295, 0.044681832632104526},
                        {0.34363344315228866, -0.058702201304027674, 0.003020221799192616},
                        {0.0015984934497727645, -9.5805104515793E-5, -1.7187445022315994E-4},
                        {0.0036936320812564253, -8.190880189616524E-4, -2.9551973622635126E-5},
                        {0.004965647007413021, 1.9626489810683822E-5, -2.4565648661058724E-4},
                        {0.24413054496442746, 0.08121207290200799, -0.005502771964648001}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {1.9790078473632116, -0.5466360053457613, 0.7442354166761601},
                        {7.5148029466680395, 1.0124425518627698, -0.8308708232627937},
                        {-6.43514862353287, 1.0832749443479637, -0.07362559817071525},
                        {0.020724339018201313, 0.002599652594099935, -0.003961656523433774},
                        {0.04599278248257384, -0.014257367173055743, 0.0016722631878542601},
                        {0.0653051230131746, -5.218785659658921E-4, -0.002695798950450919},
                        {-3.1906844150123317, -1.5369018977200506, 0.16524619704337926}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.4", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.5",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-0.00563228503, 0},
                new double[] {-0.384474276341, 0},
                new double[] {9.92344157784, Double.NaN},
                new double[][][] {
                    {
                        {-0.009836012990445158, 0.0029589899860514987, -0.0034155479511113907},
                        {-0.10184729530732066, 5.327593924820649E-4, 0.0057353060076526016},
                        {0.016303644915009106, -0.003689130235597261, -0.0025312483345677947},
                        {0.0038107513266330815, 0.003205959229414844, -8.159303166222604E-4},
                        {0.015696774651107862, -0.0037059448135031083, 2.7367574864321055E-4},
                        {0.013110199299505135, -0.0030608519096883596, -6.447454233376755E-4},
                        {0.06276193810551062, 0.00375821835084032, 0.0013984902693433081}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {0.34392939856659827, -0.09158200189901379, 0.12910527928175775},
                        {1.1760470096366977, 0.23616310164517534, -0.05057458142688201},
                        {-1.0121014533251116, 0.15248505767669873, -0.07620425570370314},
                        {0.011511611202389766, 0.10091080897132744, -0.010683387778191795},
                        {0.26561965255269504, -0.057037692991738595, 0.004170051640670491},
                        {0.06748923520446536, -0.09680681673784088, -0.004103281313994225},
                        {-0.8524954538377344, -0.24413245666460826, 0.008290175300342928}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.5", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.6",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {-0.010482021204, 0},
                new double[] {0.180367758296, 0},
                new double[] {3.43489697899, Double.NaN},
                new double[][][] {
                    {
                        {-0.008649580249889975, 0.0022485600112194097, -8.650726609048702E-4},
                        {-0.11266057761820013, 0.024694362034725272, 0.0260361864833212},
                        {0.015218945213069365, -0.00914553962283122, -0.02282021358689983},
                        {2.5066123005239126E-4, 0.031224802161732783, -0.0034016265001422597},
                        {0.08190367312390062, -0.018650512458739874, 3.4326425602877264E-4},
                        {0.0064502702618650975, -0.028665281917196323, -2.4419049163833308E-5},
                        {0.017486608039202628, -0.0017063902089100456, 7.31881057760824E-4}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {-0.2058541574059063, 0.04582469500939303, -0.032210244958406097},
                        {-1.0401985420569162, 0.37454108615584125, 0.47159019722227746},
                        {0.621062316619004, -0.20684943525738805, -0.39772097441817644},
                        {-0.10747378808409644, 0.5344440256235907, -0.051334182756052026},
                        {1.2489110819012346, -0.2984439482946815, -0.007849154127628813},
                        {-0.2968358770072303, -0.4711447951314603, 0.027030911889720764},
                        {-0.21961103396608955, 0.02162837189470477, -0.009506552851734766}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.6", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.75",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                false,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {0.0829830017977, 0},
                new double[] {1.27996801166, 0},
                new double[] {11.9168667081, Double.NaN},
                new double[][][] {
                    {
                        {-0.08520546853249629, 0.02002391005107648, -0.015337962582326878},
                        {-0.5072588908421117, 0.17043258439271164, 0.2050240947794576},
                        {0.25349423271871935, -0.0805529898238775, -0.17425804293531943},
                        {-0.0506141349576535, 0.22437919802386463, -0.02673413275920512},
                        {0.5612744057443229, -0.15045901450002216, -0.007975329992730817},
                        {-0.17340305926115046, -0.1836639508243553, 0.019217340035074935},
                        {0.0017129151303696705, -1.5973731939765025E-4, 6.403345504972649E-5}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {-0.9272841243095468, 0.22051093159987895, -0.19112897890018363},
                        {-4.804523425931252, 1.8204016841668706, 2.205227241175209},
                        {2.994591363051441, -0.8476847561428067, -1.8543652609150763},
                        {-0.6630357716342526, 2.3143378077667545, -0.3143349864556007},
                        {5.981286181791292, -1.7558797902974297, -0.1325990985326137},
                        {-2.5478012735172113, -1.7548767300847088, 0.28852239893802845},
                        {-0.033232949450470756, 0.003190852991440346, -0.0013213153097627725}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.75", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 0.9",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {0.397347160918, 0},
                new double[] {2.75351951766, 0},
                new double[] {5.37105253479, Double.NaN},
                new double[][][] {
                    {
                        {-0.31777515190383454, 0.06654081254942476, -0.061155191586619},
                        {-1.6398611619313017, 0.5966978877811928, 0.7834547554508479},
                        {1.0488049692859205, -0.2533323908962118, -0.6671198866431937},
                        {-0.2740673658300599, 0.8081103620907575, -0.10766129298323851},
                        {2.047913761610311, -0.6075094312234967, -0.06048821698815812},
                        {-0.8650389649167101, -0.6105053956998422, 0.11296922972554538},
                        {2.3913685674263586E-5, -1.84460182457579E-6, 6.030248161425271E-7}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {-2.2518864370570797, 0.34060365423349903, -0.39403631649939214},
                        {-10.009715065439687, 3.6614921088174124, 5.8097209404351},
                        {7.999725481943763, -1.2317995270963276, -5.0518709728060704},
                        {-2.770170337863771, 5.712033354860753, -0.7489455032193499},
                        {13.829283560251287, -4.313426777210048, -0.6802886194434065},
                        {-6.796164346239962, -4.168990151630222, 1.065450769329658},
                        {-0.001072855594543585, 8.733802493437212E-5, -3.0297796539917875E-5}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "0.9", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Softcoring Test: L = 1.0",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.GRAD,
                true,
                0,
                0,
                intRange(1, 8),
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[] {0.678991455919, 0},
                new double[] {2.62955540214, 0},
                new double[] {-9.66072732215, Double.NaN},
                new double[][][] {
                    {
                        {-0.5967303838923714, 0.09294823835932634, -0.10065639850802355},
                        {-2.7182918130174105, 0.9717379671351778, 1.5205922317233873},
                        {2.059991765566265, -0.34808314525235373, -1.325262659176951},
                        {-0.6963970726796157, 1.5240570186788684, -0.18878016013974244},
                        {3.628422351985993, -1.0902413269297175, -0.16292088956721837},
                        {-1.6769948479628611, -1.150418751991301, 0.2570278756685478},
                        {0.0, 0.0, 0.0}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new double[][][] {
                    {
                        {-3.3452690489458763, 0.1345828814666643, -0.35806274829197704},
                        {-10.93663217405725, 3.541368037656955, 8.988383737482236},
                        {12.27657065565445, -0.4907488150661682, -8.232590002516137},
                        {-5.952373992671734, 8.67661156047595, -0.8110392583955152},
                        {17.201456551782357, -5.049385958157128, -1.4130599622896023},
                        {-9.243751991761942, -6.81242770637627, 1.826368234010993},
                        {0.0, 0.0, 0.0}
                    },
                    // Fill in once I can get the bias-deposition actually working.
                    new double[7][3]
                },
                new String[] {"-l", "1.0", "--ac", "1-3", "--ac2", "1"},
                new String[] {},
                new String[] {}
            },
            {
                "Water-Sodium to Water Dimer Free Energy Test",
                new String[] {
                    "ffx/algorithms/structures/water-dimer.xyz",
                    "ffx/algorithms/structures/water-na.xyz"
                },
                ThermoTestMode.FREE,
                true,
                16.0,
                16.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-d", "1.0", "-l", "0.5", "--ac", "1-3", "--ac2", "1", "-n", "20000", "-Q",
                    "5000",
                    "-k", "20.0", "-w", "20.0", "-r", "5.0", "-C", "10"
                },
                new String[] {
                    "disable-neighbor-updates",
                    "true",
                    "lambda-bin-width",
                    "0.025",
                    "flambda-bin-width",
                    "5.0",
                    "randomseed",
                    "2019"
                },
                new String[] {}
            },
            {
                "Dual Well Gradient Test: L = 0.5",
                new String[] {
                    "ffx/algorithms/structures/dualWell-cis.xyz",
                    "ffx/algorithms/structures/dualWell-trans.xyz"
                },
                ThermoTestMode.GRAD,
                false,
                0,
                0,
                new int[] {0, 1, 2, 3},
                new double[] {1.0, 1.0},
                new double[] {2.0, 2.0},
                new double[] {0, Double.NaN},
                new double[][][] {
                    {
                        {8.561159028211611E-6, 5.743324479222663E-6, -1.4984967045837012E-10},
                        {-7.001815021252073E-6, -4.427074139684514E-6, 1.2219984010349385E-10},
                        {-2.1489889641599235E-6, 7.92540413855834E-9, 8.443778889936074E-12},
                        {5.896449572003855E-7, -1.3241757436767075E-6, 1.9206051464940198E-11}
                    },
                    new double[4][3]
                },
                new double[][][] {
                    {
                        {-4.504521213498369E-18, -1.6250327093073752E-17, -8.651645616845026E-13},
                        {-4.504521213498369E-18, -1.6250327093073752E-17, -8.651646266038085E-13},
                        {4.504097697024742E-18, 1.6249944935630753E-17, 8.651645891779005E-13},
                        {4.504203576143149E-18, 1.6249903576600125E-17, 8.651645991104025E-13}
                    },
                    new double[4][3]
                },
                new String[] {"-l", "0.5"},
                new String[] {
                    "pj.nt",
                    "1",
                    "lambda-bin-width",
                    "0.02",
                    "flambda-bin-width",
                    "0.20",
                    "disable-neighbor-updates",
                    "true",
                    "ost-temperOffset",
                    "6.0",
                    "randomseed",
                    "2020"
                },
                new String[0]
            },
            {
                "Dual Well 1-Step MC-OST Test",
                new String[] {
                    "ffx/algorithms/structures/dualWell-cis.xyz",
                    "ffx/algorithms/structures/dualWell-trans.xyz"
                },
                ThermoTestMode.FREE,
                true,
                0,
                1.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-l",
                    "0.5",
                    "--bM",
                    "0.1",
                    "-b",
                    "ADIABATIC",
                    "-i",
                    "VERLET",
                    "-d",
                    "2.0",
                    "-k",
                    "1000.0",
                    "-w",
                    "1000.0",
                    "-r",
                    "1.0",
                    "-Q",
                    "1000",
                    "-n",
                    "250000",
                    "-t",
                    "298.15",
                    "--tp",
                    "6.0",
                    "--mcMD",
                    "10",
                    "--mcL",
                    "0.10"
                },
                new String[] {
                    "pj.nt",
                    "1",
                    "lambda-bin-width",
                    "0.02",
                    "flambda-bin-width",
                    "0.20",
                    "disable-neighbor-updates",
                    "true",
                    "ost-temperOffset",
                    "6.0",
                    "randomseed",
                    "445"
                },
                new String[] {"--mc", "true"}
            },
            {
                "Dual Well 2-Step MC-OST Test",
                new String[] {
                    "ffx/algorithms/structures/dualWell-cis.xyz",
                    "ffx/algorithms/structures/dualWell-trans.xyz"
                },
                ThermoTestMode.FREE,
                true,
                0,
                1.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-l",
                    "0.5",
                    "--bM",
                    "0.1",
                    "-b",
                    "ADIABATIC",
                    "-i",
                    "VERLET",
                    "-d",
                    "2.0",
                    "-k",
                    "1000.0",
                    "-w",
                    "1000.0",
                    "-r",
                    "1.0",
                    "-Q",
                    "1000",
                    "-n",
                    "250000",
                    "-t",
                    "298.15",
                    "--tp",
                    "6.0",
                    "--mcMD",
                    "10",
                    "--mcL",
                    "0.10"
                },
                new String[] {
                    "pj.nt",
                    "1",
                    "lambda-bin-width",
                    "0.02",
                    "flambda-bin-width",
                    "0.20",
                    "disable-neighbor-updates",
                    "true",
                    "ost-temperOffset",
                    "6.0",
                    "randomseed",
                    "445"
                },
                new String[] {"--mc", "true", "--ts", "true"}
            },
            {
                "Short 1-Step MC-OST Test",
                new String[] {
                    "ffx/algorithms/structures/dualWell-cis.xyz",
                    "ffx/algorithms/structures/dualWell-trans.xyz"
                },
                ThermoTestMode.FREE,
                false,
                0,
                2.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-l",
                    "0.5",
                    "--bM",
                    "0.1",
                    "-b",
                    "ADIABATIC",
                    "-i",
                    "VERLET",
                    "-d",
                    "2.0",
                    "-k",
                    "1000.0",
                    "-w",
                    "1000.0",
                    "-r",
                    "0.01",
                    "-Q",
                    "100",
                    "-n",
                    "200",
                    "-t",
                    "100.0",
                    "--tp",
                    "6.0",
                    "--mcMD",
                    "10",
                    "--mcL",
                    "0.10"
                },
                new String[] {
                    "pj.nt",
                    "1",
                    "lambda-bin-width",
                    "0.02",
                    "flambda-bin-width",
                    "0.20",
                    "disable-neighbor-updates",
                    "true",
                    "ost-temperOffset",
                    "6.0",
                    "randomseed",
                    "445"
                },
                new String[] {"--mc", "true"}
            },
            {
                "Short 2-Step MC-OST Test",
                new String[] {
                    "ffx/algorithms/structures/dualWell-cis.xyz",
                    "ffx/algorithms/structures/dualWell-trans.xyz"
                },
                ThermoTestMode.FREE,
                true,
                0,
                2.0,
                null,
                null,
                null,
                null,
                null,
                null,
                new String[] {
                    "-l",
                    "0.5",
                    "--bM",
                    "0.1",
                    "-b",
                    "ADIABATIC",
                    "-i",
                    "VERLET",
                    "-d",
                    "2.0",
                    "-k",
                    "1000.0",
                    "-w",
                    "1000.0",
                    "-r",
                    "0.01",
                    "-Q",
                    "100",
                    "-n",
                    "200",
                    "-t",
                    "100.0",
                    "--tp",
                    "6.0",
                    "--mcMD",
                    "10",
                    "--mcL",
                    "0.10"
                },
                new String[] {
                    "pj.nt",
                    "1",
                    "lambda-bin-width",
                    "0.02",
                    "flambda-bin-width",
                    "0.20",
                    "disable-neighbor-updates",
                    "true",
                    "ost-temperOffset",
                    "6.0",
                    "randomseed",
                    "445"
                },
                new String[] {"--mc", "true", "--ts", "true"}
            }
        });
  }

  /**
   * Checks if two double values are approximately equal to within a tolerance.
   *
   * @param v1 One value to compare.
   * @param v2 Second value to compare.
   * @param absTol Tolerance for inequality (absolute, not relative).
   * @return True if v1 approximately equal to v2.
   */
  private static boolean approxEquals(double v1, double v2, double absTol) {
    double diff = v1 - v2;
    return Math.abs(diff) < absTol;
  }

  private static int[] intRange(int low, int high) {
    return IntStream.range(low, high).toArray();
  }

  @Test
  public void testThermodynamics() {
    if (doTest) {
      logger.info(format(" Thermodynamics test: %s\n", info));
      switch (mode) {
        case HELP:
          testHelp();
          break;
        case FREE:
          testFreeEnergy();
          break;
        case GRAD:
          testStaticGradients();
          break;
        default:
          throw new IllegalStateException(
              format(" Thermodynamics test mode %s not recognized!", mode));
      }
    } else {
      logger.info(format(" Skipping test %s: use ffx.ci true to enable!", info));
    }
  }

  private void testHelp() {
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the Thermodynamics script.
    Thermodynamics thermodynamics = new Thermodynamics(binding).run();
    algorithmsScript = thermodynamics;
  }

  private String[] assembleArgs() {
    List<String> argList = new ArrayList<>(flags);
    opts.forEach(
        (String k, String v) -> {
          argList.add(k);
          argList.add(v);
        });
    argList.addAll(Arrays.stream(copiedFiles).map(File::getPath).collect(Collectors.toList()));
    return argList.toArray(new String[0]);
  }

  /**
   * Test a free energy evaluation; very expensive, so likely only done for one or a few tests.
   * Currently not implemented.
   */
  private void testFreeEnergy() {
    binding.setVariable("args", assembleArgs());

    // Construct and evaluate the Thermodynamics script.
    Thermodynamics thermodynamics = new Thermodynamics(binding);
    thermodynamics.setProperties(algorithmConfig);
    algorithmsScript = thermodynamics;

    thermodynamics.run();
    OrthogonalSpaceTempering orthogonalSpaceTempering = thermodynamics.getOST();
    Histogram histogram = orthogonalSpaceTempering.getHistogram();
    double delG = histogram.updateFLambda(true, false);
    assertEquals(
        format(" Test %s: not within tolerance %12.5g", info, feTol),
        freeEnergy,
        delG,
        feTol);
  }

  /** Tests gradients & energies for a static structure, before and after dropping a bias. */
  private void testStaticGradients() {

    binding.setVariable("args", assembleArgs());

    // Construct and evaluate the Thermodynamics script.
    Thermodynamics thermodynamics = new Thermodynamics(binding);
    thermodynamics.setProperties(algorithmConfig);
    algorithmsScript = thermodynamics;
    thermodynamics.run();

    OrthogonalSpaceTempering orthogonalSpaceTempering = thermodynamics.getOST();
    orthogonalSpaceTempering.setPropagateLambda(false);
    CrystalPotential under = thermodynamics.getPotential();
    int nVars = orthogonalSpaceTempering.getNumberOfVariables();

    double[] x = new double[nVars];
    x = orthogonalSpaceTempering.getCoordinates(x);
    double[] gOSTPre = new double[nVars];
    double[] gUnderPre = new double[nVars];

    logger.info(" Testing the OST potential before bias added.");
    EnergyResult ostPre = testGradientSet("Unbiased OST", orthogonalSpaceTempering, x, gOSTPre, 0);
    logger.info(" Testing the underlying CrystalPotential before bias added.");
    EnergyResult underPre = testGradientSet("Unbiased potential", under, x, gUnderPre, 0);

    // Assert that, before biases, OST and underlying potential are equal.
    ostPre.assertResultsEqual(underPre);
  }

  /**
   * Generates and tests an EnergyResult against tabulated values.
   *
   * @param description Description (such as unbiased OST)
   * @param potential A CrystalPotential (either OST or underlying).
   * @param x Coordinates.
   * @param g Array to add gradients to.
   * @param tableIndex 0 for unbiased potential, 1 for a biased potential.
   * @return The generated EnergyResult.
   */
  private EnergyResult testGradientSet(
      String description, CrystalPotential potential, double[] x, double[] g, int tableIndex) {
    assertTrue(
        format(" Potential %s is not a lambda interface!", potential),
        potential instanceof LambdaInterface);

    EnergyResult er = new EnergyResult(description, potential, x, g);

    checkThGradScalar(er.energy, pe, tableIndex, peTol, "potential energy");
    checkThGradScalar(er.firstLam, dudl, tableIndex, dudlTol, "dU/dL");
    checkThGradArray(er.gradient, dudx, tableIndex, dudxTol, "dU/dX gradient");
    if (er.hasSecondDerivatives) {
      checkThGradScalar(er.secondLam, d2udl2, tableIndex, d2udl2Tol, "d2U/dL2");
      checkThGradArray(er.lamGradient, d2udxdl, tableIndex, d2udxdlTol, "d2U/dXdL gradient");
    }

    return er;
  }

  /**
   * Checks a scalar value against its expected value. Mostly a convenience formatting method.
   *
   * @param actual Value from the test.
   * @param expected Array of expected values (see tableIndex).
   * @param tableIndex 0 for unbiased potential, 1 for biased OST potential.
   * @param tol Tolerance for this test.
   * @param description Scalar to be tested.
   */
  private void checkThGradScalar(
      double actual, double[] expected, int tableIndex, double tol, String description) {
    if (debugMode) {
      logger.info(format(" %s is %20.12g", description, actual));
    } else {
      assertEquals(
          format(
              " Expected %s %12.6g, received %12.6g from test %s",
              description, expected[tableIndex], actual, info),
          expected[tableIndex],
          actual,
          tol);
    }
  }

  /**
   * Checks an array value (generally 1-D flat) against its expected value (generally 3-D; indices
   * tableIndex, then atom number, then X/Y/Z).
   *
   * @param actual Array from the test, flat.
   * @param expected Array of expected values, indexed by tableIndex, atoms, and XYZ.
   * @param tableIndex 0 for unbiased potential, 1 for biased OST potential.
   * @param tol Tolerance for this test.
   * @param description Array to be tested.
   */
  private void checkThGradArray(
      double[] actual, double[][][] expected, int tableIndex, double tol, String description) {
    double[] actualSlice = new double[3];
    for (int i = 0; i < numGradAtoms; i++) {
      int i3 = i * 3;
      arraycopy(actual, i3, actualSlice, 0, 3);
      if (debugMode) {
        logger.info(
            format(" %s at atom %d is %s", description, i, Arrays.toString(actualSlice)));
      } else {
        double[] exp = expected[tableIndex][i];
        assertArrayEquals(
            format(
                " Discrepancy found on array of %s for test %s on atom %d."
                    + "\n Expected: %s\n Found: %s",
                description, info, i, Arrays.toString(exp), Arrays.toString(actualSlice)),
            exp,
            actualSlice,
            tol);
      }
    }
  }

  private enum ThermoTestMode {
    HELP,
    FREE,
    GRAD
  }

  /** Contains the result of an energy evaluation: potential energy and several derivatives. */
  private class EnergyResult {

    final double energy;
    final double firstLam;
    final double secondLam;
    final int nVars;
    final boolean hasSecondDerivatives;
    private final String description;

    private final double[] gradient;
    private final double[] lamGradient;

    public EnergyResult(String description, CrystalPotential potential, double[] x, double[] g) {
      this.description = description;

      LambdaInterface linter = (LambdaInterface) potential;
      energy = potential.energyAndGradient(x, g, false);
      firstLam = linter.getdEdL();
      nVars = g.length;
      gradient = copyOf(g, nVars);

      hasSecondDerivatives = !(potential instanceof OrthogonalSpaceTempering);
      if (hasSecondDerivatives) {
        secondLam = linter.getd2EdL2();
        lamGradient = new double[g.length];
        linter.getdEdXdL(lamGradient);
      } else {
        secondLam = 0;
        lamGradient = null;
      }
    }

    /**
     * Asserts that two energy results are equivalent to each other. Always tests U, dU/dX, and
     * dU/dL. Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
     *
     * <p>All values must be identical to tolerance.
     *
     * @param other Another EnergyResult
     */
    public void assertResultsEqual(EnergyResult other) {
      assertEquals(
          format(
              " Test %s: potential energy for %s did not " + "match %s, %12.6g != %12.6g.",
              info, this.toString(), other.toString(), this.energy, other.energy),
          this.energy,
          other.energy,
          peTol);
      assertEquals(
          format(
              " Test %s: dU/dL for %s did not " + "match %s, %12.6g != %12.6g.",
              info, this.toString(), other.toString(), this.firstLam, other.firstLam),
          this.firstLam,
          other.firstLam,
          dudlTol);
      assertArrayEquals(
          format(
              " Test %s: dU/dX for %s did not match %s.", info, this.toString(), other.toString()),
          this.gradient,
          other.gradient,
          dudxTol);

      if (hasSecondDerivatives && other.hasSecondDerivatives) {
        assertEquals(
            format(
                " Test %s: d2U/dL2 for %s did not " + "match %s, %12.6g != %12.6g.",
                info, toString(), other.toString(), this.secondLam, other.secondLam),
            this.secondLam,
            other.secondLam,
            d2udl2Tol);
        assertArrayEquals(
            format(
                " Test %s: d2U/dXdL for %s did not match %s.",
                info, this.toString(), other.toString()),
            this.lamGradient,
            other.lamGradient,
            d2udxdlTol);
      }
    }

    /**
     * Asserts that two energy results are not equivalent to each other. Always tests U, dU/dX, and
     * dU/dL. Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
     *
     * <p>Prints all equalities, fails if no inequalities are found. At least one value must not be
     * identical to tolerance.
     *
     * @param other Another EnergyResult
     */
    public void assertResultsInequal(EnergyResult other) {
      int diffsFound = 0;
      StringBuilder sb = new StringBuilder(" Equalities found between ");
      sb = sb.append(toString()).append(" and ").append(other.toString());
      sb = sb.append(" for test ").append(info);
      boolean equalFound = false;

      if (approxEquals(energy, other.energy, peTol)) {
        equalFound = true;
        sb.append(
            format(
                "\n Potential energy %12.6g == other potential energy %12.6g",
                energy, other.energy));
      } else {
        ++diffsFound;
      }

      if (approxEquals(firstLam, other.firstLam, dudlTol)) {
        equalFound = true;
        sb.append(
            format(
                "\n Lambda derivative %12.6g == other lambda derivative %12.6g",
                firstLam, other.firstLam));
      } else {
        ++diffsFound;
      }

      for (int i = 0; i < nVars; i++) {
        if (approxEquals(gradient[i], other.gradient[i], dudxTol)) {
          equalFound = true;
          sb.append(
              format(
                  "\n Gradient %d %12.6g == other gradient %12.6g",
                  i, gradient[i], other.gradient[i]));
        } else {
          ++diffsFound;
        }
      }

      if (hasSecondDerivatives && other.hasSecondDerivatives) {
        if (approxEquals(secondLam, other.secondLam, d2udl2Tol)) {
          equalFound = true;
          sb.append(
              format(
                  "\n Lambda derivative %12.6g == other lambda derivative %12.6g",
                  secondLam, other.secondLam));
        } else {
          ++diffsFound;
        }

        for (int i = 0; i < nVars; i++) {
          if (approxEquals(lamGradient[i], other.lamGradient[i], d2udxdlTol)) {
            equalFound = true;
            sb.append(
                format(
                    "\n Lambda gradient %d %12.6g == other lambda gradient %12.6g",
                    i, lamGradient[i], other.lamGradient[i]));
          } else {
            ++diffsFound;
          }
        }
      }

      if (equalFound) {
        logger.info(sb.toString());
      }
      assertTrue(
          format(
              " No inequalities found between %s and %s for test %s!",
              this.toString(), other.toString(), info),
          diffsFound > 0);
    }

    @Override
    public String toString() {
      return description;
    }
  }
}
