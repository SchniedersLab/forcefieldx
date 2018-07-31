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
package ffx.algorithms.groovy;

import ffx.algorithms.AbstractOSRW;
import ffx.algorithms.PJDependentTest;
import ffx.algorithms.TransitionTemperedOSRW;
import ffx.algorithms.groovy.NewThermodynamics;

import ffx.crystal.CrystalPotential;
import ffx.potential.bonded.LambdaInterface;
import groovy.lang.Binding;
import static org.junit.Assert.*;

import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.MapConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;
import java.util.stream.IntStream;
import java.util.stream.Collectors;

/**
 * Tests the functionality of the Transition-Tempered OSRW algorithm in Force Field X,
 * both by running a very simple and quick TT-OSRW run, and setting up a number of systems,
 * comparing output energies and gradients from the starting algorithmConfig.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class ThermodynamicsTest extends PJDependentTest {

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        /**
         * Test info, filenames, mode, dG, tol(dG), grad atoms, PE, dU/dL, d2U/dL2,
         * dU/dX, d2U/dXdL, Groovy options, properties, Groovy flags
         */
        return Arrays.asList(new Object[][]{
            {
                "Thermodynamics Help Message Test", new String[]{}, ThermoTestMode.HELP, 0, 0, null, null,
                null, null, null, null, new String[]{}, new String[]{}, new String[]{"-h", "true"}
            },
            {
                "Acetamide Implicit Solvation Free Energy: -10.8 kcal/mol", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.FREE, -10.8, 1.0, null, null, null, null, null, null,
                new String[]{"-C", "10", "-d", "1.0", "-n", "50000", "-w", "5", "--bM", "0.25", "--tp", "4"},
                new String[]{}, new String[]{}
            },
            {
                "Acetamide Implicit Solvation Gradients: L = 0.9", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.GRAD, 0, 0, intRange(1, 9),
                new double[]{-14.3232441, -14.1547628}, new double[]{-47.0065247, 2.88020118},
                new double[]{-626.753662, Double.NaN},
                new double[][][]{
                    {
                        {
                            0.3454325973545622,-0.3556977754400722,0.22632332784499146
                        },
                        {
                            -2.243711749882479,-0.19356781608568252,-0.7607652107647813
                        },
                        {
                            -0.4691314984840982,0.41778512359780806,-0.29489746509549386
                        },
                        {
                            -0.13799993351680084,0.10993253720542051,-0.08060781131332917
                        },
                        {
                            -0.12916866666434035,-0.1292595807344552,-0.19923738828440363
                        },
                        {
                            -0.14207609439404922,0.1277438883781784,-0.09272512156557061
                        },
                        {
                            -0.23004200820343945,-0.024067408991374717,0.12377960757109271
                        },
                        {
                            0.7198849556438509,-0.058542226119262075,0.28168889775425177
                        }
                    },
                    {
                        {
                            0.3245837863444,-0.3346459581204537,0.21247087522245325
                        },
                        {
                            -2.1088804470365896,-0.18290359797842182,-0.7158477527631418
                        },
                        {
                            -0.4416189323792037,0.3919514976823785,-0.27720438558513366
                        },
                        {
                            -0.1294804178454077,0.10345020607616696,-0.07562848810923806
                        },
                        {
                            -0.12170071648317575,-0.12226365328373233,-0.18708155978049676
                        },
                        {
                            -0.1338286311435258,0.12067970333834456,-0.08762615563248287
                        },
                        {
                            -0.21577447110516357,-0.022826463023161717,0.11687334540079644
                        },
                        {
                            0.6773203131299521,-0.05386903789837927,0.26508674941152166
                        }
                    }
                }, new double[][][]{
                    {
                        {
                            6.925790792595476,-6.993227694785508,4.601662357689225
                        },
                        {
                            -44.78976740441211,-3.542559023701438,-14.921182646958457
                        },
                        {
                            -9.139431352563745,8.581702257132498,-5.877484673144383
                        },
                        {
                            -2.8301078255996117,2.1533731217409113,-1.6540871699322346
                        },
                        {
                            -2.4807870616248175,-2.3239852814755535,-4.03805881724825
                        },
                        {
                            -2.7397344153056427,2.3466598494381614,-1.6938314273137203
                        },
                        {
                            -4.7395497527455746,-0.412231285350944,2.294199267621269
                        },
                        {
                            14.139598131960081,-1.5523918335908964,5.515086978981156
                        }
                    }, new double[8][3]
                },
                new String[]{"-l", "0.9"}, new String[]{}, new String[]{}
            },
            {
                "Acetamide Implicit Solvation Gradients: L = 1.0", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.GRAD, 0, 0, intRange(1, 9),
                new double[]{-22.8540579, -22.4300811}, new double[]{-130.57368, -4.28073063},
                new double[]{-1044.58944, Double.NaN},
                new double[][][]{
                    {
                        {
                            1.602335371,-1.624839098,1.06143983
                        },
                        {
                            -10.37222509,-0.836476676,-3.468683543
                        },
                        {
                            -2.12776904,1.975205163,-1.361552091
                        },
                        {
                            -0.651612094,0.500729882,-0.380794001
                        },
                        {
                            -0.579385578,-0.551019873,-0.932070285
                        },
                        {
                            -0.639287155,0.553619194,-0.400124158
                        },
                        {
                            -1.090182519,-0.098879753,0.540134289
                        },
                        {
                            3.285960172,-0.340272596,1.282575053
                        }
                    },
                    {
                        {
                            1.670608498,-1.693777005,1.106802141
                        },
                        {
                            -10.81375409,-0.871398548,-3.61577372
                        },
                        {
                            -2.217863814,2.05980195,-1.419491216
                        },
                        {
                            -0.679510758,0.521957424,-0.397099678
                        },
                        {
                            -0.603840704,-0.573929277,-0.9718767
                        },
                        {
                            -0.666294935,0.57675212,-0.416821626
                        },
                        {
                            -1.136904098,-0.102943451,0.562750069
                        },
                        {
                            3.425345638,-0.355575779,1.33694173
                        }
                    }
                }, new double[][][]{
                    {
                        {
                            19.23830776,-19.42563249,12.78239544
                        },
                        {
                            -124.4160206,-9.840441733,-41.44772957
                        },
                        {
                            -25.38730931,23.83806183,-16.32634631
                        },
                        {
                            -7.861410627,5.981592005,-4.594686583
                        },
                        {
                            -6.891075171,-6.455514671,-11.21683005
                        },
                        {
                            -7.610373376,6.518499582,-4.705087298
                        },
                        {
                            -13.16541598,-1.145086904,6.372775743
                        },
                        {
                            39.27666148,-4.312199538,15.31968605
                        }
                    }, new double[8][3]
                },
                new String[]{"-l", "1.0"}, new String[]{}, new String[]{}
            },
            /*{
                "Acetamide Crystallization Gradients: L = 1.0", new String[]{"ffx/algorithms/structures/acetamide.xtal.xyz"},
                ThermoTestMode.GRAD, 0, 0, intRange(1, 9),
                new double[]{-31.0594281, -29.4146512}, new double[]{-618.267325, -621.952472},
                new double[]{-2210.57780, Double.NaN},
                new double[][][]{
                    {
                        {
                            0.0035512833043158665,-0.003887887320341954,0.003235545226927705
                        },
                        {
                            0.0021285361170626516,-4.49299151844218E-4,0.0021704548936494206
                        },
                        {
                            0.006734914213447851,-0.01339944324718978,6.90125313571599E-4
                        },
                        {
                            0.015842367597945994,-0.0016054178146653886,0.003997206107637208
                        },
                        {
                            -0.0015637067839852914,0.007092319275139369,-0.0027402606404022256
                        },
                        {
                            -0.0121741042114627,-0.006474092419587185,-0.005616416838037305
                        },
                        {
                            -0.01403124580288001,0.005646687947437945,-0.014774772293374261
                        },
                        {
                            0.001034108856084348,0.008938210799334989,-0.0032032891059294855
                        }
                    },
                    {

                        {
                            0.01230167056545174,0.013766350856164476,0.005324850693954207
                        },
                        {
                            -0.01592928349640763,0.013108020828172438,-0.023387576357417714
                        },
                        {
                            0.016673516267039817,0.03472511805283299,-0.029486485771014235
                        },
                        {
                            -0.03579194456319919,0.002462865399243487,-0.008218460449024072
                        },
                        {
                            -0.0060197970719967665,-0.0016748080087555465,-0.006106837158237614
                        },
                        {
                            -0.018223478181111117,-0.010099259335422165,-0.0078850771316659
                        },
                        {
                            -0.016771705343051856,-0.0033391238062166328,-0.01475968020657157
                        },
                        {
                            0.029666060795302252,-0.008765333776555077,0.032179277056848456
                        }
                    }
                }, new double[][][]{
                    {
                        {
                            5.249020006689424,10.590096943813828,1.2532938108026226
                        },
                        {
                            -10.832189890529506,8.132513646368627,-15.331277732618949
                        },
                        {
                            5.961784257198348,28.868069210936916,-18.10178573702526
                        },
                        {
                            -30.97343345811145,2.440406275220606,-7.327707475959026
                        },
                        {
                            -2.6730367897283713,-5.259061701045617,-2.0194794777431255
                        },
                        {
                            -3.6287862521994168,-2.1745978893208817,-1.3608818574839396
                        },
                        {
                            -1.6438960384957468,-5.390242084572063,0.009053161100865514
                        },
                        {
                            17.175204259475937,-10.61967395213345,21.224637508480413
                        }
                    }, new double[8][3]
                },
                new String[]{"-l", "1.0", "--la1", "1-9"}, new String[]{}, new String[]{}
            }*/
        });
    }

    private static int[] intRange(int low, int high) {
        return IntStream.range(low, high).toArray();
    }

    // Enables behavior where it prints out instead of testing information.
    private static final boolean debugMode = false;
    /**
     * Set of default options, such as "-d", "0.5" to specify a half-femtosecond timestep.
     * Can be over-ridden in test constructors.
     */
    private static final Map<String, String> DEFAULT_OPTIONS;
    /**
     * Set of default properties, such as "ttosrw-alwaystemper", "true" to use standard tempering scheme.
     * Can be over-ridden in test constructors.
     */
    private static final Map<String, String> DEFAULT_PROPERTIES;
    /**
     * Set of default boolean options, such as "-y", false to specify "do not add the -y --synchronous flag".
     * Can be over-ridden in test constructors.
     * Maps a flag to whether it should be included by default.
     */
    private static final Map<String, Boolean> DEFAULT_FLAGS;

    /**
     * Default number of atomic gradients to test; pick the first N (3) if ffx.ci is false. Otherwise, use
     * the whole array to test.
     */
    private static final int DEFAULT_GRADIENT_EVALS = 3;

    static {
        String[] opts = {"--bM", "0.05",
                "--dw", "OFF",
                "--lf", "1.0E-18",
                "--lm", "1.0E-18",
                "--np", "1",
                "--sf", "1.0",
                "--tp", "2.0",      // Set low to encourage fast tempering for toy cases.
                "-b", "ADIABATIC",  // Used for Langevin dynamics.
                "-C", "5000",       // By default, evaluate 0 dynamics steps, and trigger adding counts manually.
                "-d", "0.5",        // Most of the tests are in vacuum, thus short timestep.
                "-F", "XYZ",
                "-i", "STOCHASTIC", // Stochastic because most toy systems work best with Langevin dynamics.
                "-k", "5.0",        // Infrequent because we should not need to restart.
                "-l", "1.0",        // Set walker to start at the usually-physical topology. Likely commonly set by test.
                "-n", "0",          // Do not run any steps by default, just evaluate energy.
                "-p", "0",          // Use NVT, not NPT dynamics.
                "-Q", "0",          // Do not run any steps by default, just evaluate energy.
                "-r", "2.0",        // Very infrequent printouts, since it should work.
                "-t", "298.15",
                "-w", "10.0"};      // Infrequent since we should not need the trajectory.

        int nOpts = opts.length;
        Map<String, String> optMap = new HashMap<>(nOpts);
        for (int i = 0; i < nOpts; i += 2) {
            optMap.put(opts[i], opts[i+1]);
        }
        DEFAULT_OPTIONS = Collections.unmodifiableMap(optMap);

        String[] props = {"ttosrw-alwaystemper", "true",
                "ttosrw-temperOffset", "0.5", // Small to temper fast.
                "print-on-failure", "false",
                "disable-neighbor-updates", "false",
                "vdw-lambda-alpha", "0.25",
                "vdw-lambda-exponent", "3",
                "permanent-lambda-alpha", "1.0",
                "permanent-lambda-exponent", "3"
        };

        int nProps = props.length;
        Map<String, String> propMap = new HashMap<>(nProps);
        for (int i = 0; i < nProps; i += 2) {
            propMap.put(props[i], props[i+1]);
        }
        DEFAULT_PROPERTIES = Collections.unmodifiableMap(propMap);

        Map<String, Boolean> flags = new HashMap<>();
        flags.put("--rn", false);
        flags.put("-y", false);
        flags.put("-h", false);
        DEFAULT_FLAGS = Collections.unmodifiableMap(flags);
    }

    Binding binding;
    NewThermodynamics thermo;

    private final String info;
    private final ThermoTestMode mode;
    private final String[] filenames;
    private final File tempDir;
    private final File[] copiedFiles;
    private final Map<String, String> opts;
    private final List<String> flags;
    /**
     * Configuration containing the properties to be used by TT-OSRW.
     */
    Configuration algorithmConfig;
    private final boolean ffxCI;

    /**
     * Free energy associated with this test, if calculated.
     * Not used alongside testing gradients, and values are meaningless in that case.
     */
    private final double freeEnergy;
    private final double feTol;

    /**
     * Set of atoms for which tabulated gradients (dU/dL, d2U/dXdL) should be evaluated.
     * Not used alongside testing free energies, and values are meaningless (and often null) in that case.
     */
    private final int[] gradientAtoms;
    private final int numGradAtoms;

    /**
     * Potential energy of the underlying Potential (0), and OSRW after one bias drop (1).
     */
    private final double[] pe;
    private final double[] dudl;
    private final double[] d2udl2;
    /**
     * Gradient of the underlying Potential ([0][0..n][0-2]), and OSRW after one bias drop ([1][0..n][0-2]).
     * Axes: before/after bias, atoms, X/Y/Z.
     */
    private final double[][][] dudx;
    private final double[][][] d2udxdl;

    // Currently, everything is just 1.0E-5 kcal/mol per unit.
    private static final double DEFAULT_PE_TOL = 1.0E-5;      // Units: none.
    private static final double DEFAULT_DUDL_TOL = 1.0E-5;    // Units: per lambda (quasi-unitless).
    private static final double DEFAULT_D2UDL2_TOL = 1.0E-5;  // Units: per lambda^2 (quasi-unitless).
    private static final double DEFAULT_DUDX_TOL = 1.0E-5;    // Units: per Angstrom.
    private static final double DEFAULT_D2UDXDL_TOL = 1.0E-5; // Units: per Angstrom per lambda.

    // TODO: Make these settable if needed.
    private final double peTol = DEFAULT_PE_TOL;
    private final double dudlTol = DEFAULT_DUDL_TOL;
    private final double d2udl2Tol = DEFAULT_D2UDL2_TOL;
    private final double dudxTol = DEFAULT_DUDX_TOL;
    private final double d2udxdlTol = DEFAULT_D2UDXDL_TOL;

    public ThermodynamicsTest(String info, String[] filenames, ThermoTestMode mode,
                              double freeEnergy, double feTol,
                              int[] gradientAtoms,
                              double[] pe, double[] dudl, double[] d2udl2,
                              double[][][] dudx, double[][][] d2udxdl,
                              String[] options, String[] properties, String[] flags) throws IOException {
        this.info = info;
        this.mode = mode;
        int nFiles = filenames.length;

        String tempDirName = String.format("temp-%016x/", new Random().nextLong());
        tempDir = new File(tempDirName);
        tempDir.mkdir();
        logger.fine(String.format(" Running test %s in directory %s", info, tempDir));
        this.filenames = Arrays.copyOf(filenames, nFiles);
        copiedFiles = new File[nFiles];

        String[] copiedExtensions = new String[]{"dyn", "key", "properties", "his", "lam"};
        for (int i = 0; i < nFiles; i++) {
            File srcFile = new File("src/main/java/" + filenames[i]);
            File tempFile = new File(tempDirName + FilenameUtils.getName(filenames[i]));
            FileUtils.copyFile(srcFile, tempFile);
            logger.info(String.format(" Copied file %s to %s", srcFile, tempFile));
            copiedFiles[i] = tempFile;

            for (String ext : copiedExtensions) {
                srcFile = new File(String.format("%s.%s", FilenameUtils.removeExtension(srcFile.getPath()), ext));
                if (srcFile.exists()) {
                    logger.fine(" Copying extension " + ext);
                    tempFile = new File(String.format("%s.%s", FilenameUtils.removeExtension(tempFile.getPath()), ext));
                    logger.info(String.format(" Copied file %s to %s", srcFile, tempFile));
                    FileUtils.copyFile(srcFile, tempFile);
                }
            }
        }
        ffxCI = Boolean.parseBoolean(System.getProperty("ffx.ci", "false"));

        switch (mode) {
            case HELP:
                assertTrue(String.format("Help tests must have no file arguments, found %d for %s", nFiles, info), nFiles == 0);
                break;
            default:
                assertTrue(String.format("Must have 1, 2, or 4 distinct filenames, found %d for test %s", nFiles, info),
                        nFiles == 1 || nFiles == 2 || nFiles == 4);
                for (int i = 0; i < nFiles; i++) {
                    for (int j = i+1; j < nFiles; j++) {
                        assertNotEquals(String.format(" Filenames %d and %d matched in test %s: files %s and %s", i, j, info, filenames[i], filenames[j]),
                                filenames[i], filenames[j]);
                    }
                }
                break;
        }
        int nOpts = options.length;
        int nProps = properties.length;
        int nFlags = flags.length;

        if (nOpts > 0) {
            assertTrue(String.format("Unmatched option key %s for test %s", options[nOpts - 1], info),
                    options.length % 2 == 0);
        }
        if (nProps > 0) {
            assertTrue(String.format("Unmatched property key %s for test %s", properties[nProps - 1], info),
                    properties.length % 2 == 0);
        }
        if (nFlags > 0) {
            assertTrue(String.format("Unmatched flag key %s for test %s", flags[nFlags - 1], info),
                    flags.length % 2 == 0);
        }

        Pattern validOption = Pattern.compile("^--?[^D]");
        Pattern validProperty = Pattern.compile("^[^-]");

        Map<String, String> groovyOpts = new HashMap<>(DEFAULT_OPTIONS);
        for (int i = 0; i < nOpts; i += 2) {
            String opti = options[i];
            assertTrue(String.format(" Option %s for test %s does not look like a Groovy option!", opti, info),
                    validOption.matcher(opti).find());
            groovyOpts.put(opti, options[i+1]);
        }
        this.opts = Collections.unmodifiableMap(groovyOpts);

        Map<String, String> addedProps = new HashMap<>(DEFAULT_PROPERTIES);
        for (int i = 0; i < nProps; i += 2) {
            String propi = properties[i];
            assertTrue(String.format(" Property %s for test %s does not look like a property!", propi, info),
                    validProperty.matcher(propi).find());
            addedProps.put(propi, properties[i+1]);
        }
        algorithmConfig = new MapConfiguration(addedProps);


        Map<String, Boolean> addedFlags = new HashMap<>(DEFAULT_FLAGS);
        for (int i = 0; i < nFlags; i+= 2) {
            String flagi = flags[i];
            assertTrue(String.format(" Flag %s for test %s does not look like a flag!", flagi, info),
                    validOption.matcher(flagi).find());
            String vali = flags[i+1];
            assertTrue(String.format(" Value %s for flag %s in test %s is not a true/false value!", vali, flagi, info),
                    vali.equalsIgnoreCase("TRUE") || vali.equalsIgnoreCase("FALSE"));
            addedFlags.put(flagi, Boolean.parseBoolean(vali));
        }
        this.flags = Collections.unmodifiableList(
                addedFlags.entrySet().stream().
                        filter(Map.Entry::getValue).
                        map(Map.Entry::getKey).
                        collect(Collectors.toList()));
        
        // Only meaningful for free energy evaluations.
        this.freeEnergy = freeEnergy;
        this.feTol = feTol;
        if (mode == ThermoTestMode.FREE) {
            assertTrue(String.format(" Free energy tolerance for test %s was %10.4g <= 0!", info, feTol), feTol > 0);
        }

        // Only meaningful for gradient evaluations.
        this.dudx = new double[2][][];
        this.d2udxdl = new double[2][][];
        boolean isGradient = (mode == ThermoTestMode.GRAD);
        if (isGradient) {
            numGradAtoms = ffxCI ? gradientAtoms.length : Math.min(DEFAULT_GRADIENT_EVALS, gradientAtoms.length);
            this.gradientAtoms = Arrays.copyOf(gradientAtoms, numGradAtoms);

            if (debugMode) {
                this.pe = this.dudl = this.d2udl2 = null;
            } else {
                assertNotNull(gradientAtoms);
                assertNotNull(pe);
                assertNotNull(dudl);
                assertNotNull(d2udl2);

                this.pe = Arrays.copyOf(pe, 2);
                this.dudl = Arrays.copyOf(dudl, 2);
                this.d2udl2 = Arrays.copyOf(d2udl2, 2);
                for (int i = 0; i < 2; i++) {
                    assertNotNull(dudx[i]);
                    assertNotNull(d2udxdl[i]);
                    this.dudx[i] = new double[numGradAtoms][3];
                    this.d2udxdl[i] = new double[numGradAtoms][3];

                    for (int j = 0; j < numGradAtoms; j++) {
                        assertNotNull(dudx[i][j]);
                        assertNotNull(d2udxdl[i][j]);
                        System.arraycopy(dudx[i][j], 0, this.dudx[i][j], 0, 3);
                        System.arraycopy(d2udxdl[i][j], 0, this.d2udxdl[i][j], 0, 3);
                    }
                }
            }
        } else {
            this.gradientAtoms = null;
            numGradAtoms = 0;
            this.pe = null;
            this.dudl = null;
            this.d2udl2 = null;
            for (int i = 0; i < 2; i++) {
                this.dudx[i] = null;
                this.d2udxdl[i] = null;
            }
        }
    }

    @Before
    public void before() {
        binding = new Binding();
        thermo = new NewThermodynamics();
        thermo.setBinding(binding);
    }

    @After
    public void after() throws IOException {
        // Clean up the temporary directory if it exists.
        if (tempDir != null && tempDir.exists()) {
            if (tempDir.isDirectory()) {
                for (File file : tempDir.listFiles()) {
                    if (file.isDirectory()) {
                        for (File file2 : file.listFiles()) {
                            // Should never have to go more than 2-deep into the temporary directory.
                            file2.delete();
                        }
                    }
                    file.delete();
                }
            } else {
                logger.warning(String.format(" Expected %s (temporary directory) to be a directory!", tempDir));
            }
            tempDir.delete();
        }

        if (thermo.getOSRW() == null) {
            assert mode == ThermoTestMode.HELP;
        } else {
            thermo.destroyPotentials();
        }
        System.gc();
    }

    @Test
    public void testThermodynamics() {
        logger.info(String.format(" Thermodynamics test: %s\n", info));
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
                throw new IllegalStateException(String.format(" Thermodynamics test mode %s not recognized!", mode));
        }
    }

    private void testHelp() {
        String[] args = {"-h"};
        binding.setVariable("args", args);
        thermo.run();
    }

    private String[] assembleArgs() {
        List<String> argList = new ArrayList<>(flags);
        opts.forEach((String k, String v) -> {
            argList.add(k);
            argList.add(v);
        });
        argList.addAll(Arrays.stream(copiedFiles).map(File::getPath).collect(Collectors.toList()));
        return argList.toArray(new String[argList.size()]);
    }

    private void assembleThermo() {
        String[] args = assembleArgs();
        binding.setVariable("args", args);
        thermo = new NewThermodynamics();
        thermo.setBinding(binding);
        thermo.setProperties(algorithmConfig);
    }

    /**
     * Test a free energy evaluation; very expensive, so likely only done for one or a few tests.
     * Currently not implemented.
     */
    private void testFreeEnergy() {
        assembleThermo();
        thermo.run();
        AbstractOSRW osrw = thermo.getOSRW();
        double delG = osrw.lastFreeEnergy();
        assertEquals(String.format(" Test %s: not within tolerance %12.5g", info, feTol), freeEnergy, delG, feTol);
    }

    /**
     * Tests gradients & energies for a static structure, before and after dropping a bias.
     */
    private void testStaticGradients() {
        assembleThermo();
        thermo.run();
        AbstractOSRW osrw = thermo.getOSRW();
        osrw.setPropagateLambda(false);
        CrystalPotential under = thermo.getPotential();
        int nVars = osrw.getNumberOfVariables();

        double[] x = new double[nVars];
        x = osrw.getCoordinates(x);
        double[] gOSRWPre = new double[nVars];
        double[] gUnderPre = new double[nVars];
        double[] gOSRWPost = new double[nVars];
        double[] gUnderPost = new double[nVars];

        logger.info(" Testing the OSRW potential before bias added.");
        EnergyResult osrwPre = testGradientSet("Unbiased OSRW", osrw, x, gOSRWPre, 0);
        logger.info(" Testing the underlying CrystalPotential before bias added.");
        EnergyResult underPre = testGradientSet("Unbiased potential",under, x, gUnderPre, 0);

        // Assert that, before biases, OSRW and underlying potential are equal.
        osrwPre.assertResultsEqual(underPre);

        double currentdUdL = osrwPre.firstLam;
        logger.info(String.format(" Adding an OSRW bias at lambda %8.4g, dU/dL %14.8g", osrw.getLambda(), currentdUdL));
        osrw.addBias(currentdUdL, x, gOSRWPre);

        // Wait for the bias to be received by the OSRW object.
        boolean biasReceived = false;
        for (int i = 0; i < 200; i++) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                logger.warning(String.format(" Interrupted at %d 100-msec wait for bias to be added!", i));
            }
            if (osrw.getCountsReceived() > 0) {
                logger.fine(String.format(" Required %d 100-msec waits for bias to be added.", i));
                biasReceived = true;
                break;
            }
        }
        try {
            Thread.sleep(500);
        } catch (InterruptedException ex) {
            logger.warning(" Interrupted at final 500-msec wait for bias to be added!");
        }
        assertTrue(" No bias was received by the OSRW over 20 seconds!", biasReceived);

        logger.info(" Testing the OSRW potential after bias added.");
        EnergyResult osrwPost = testGradientSet("Biased OSRW", osrw, x, gOSRWPost, 1);
        logger.info(" Testing the underlying CrystalPotential after bias added.");
        EnergyResult underPost = testGradientSet("Biased potential", under, x, gUnderPost, 0);

        // Assert that bias only affects the OSRW potential, and does affect the OSRW potential.
        underPre.assertResultsEqual(underPost);
        osrwPre.assertResultsInequal(osrwPost);
    }

    /**
     * Generates and tests an EnergyResult against tabulated values.
     *
     * @param description Description (such as unbiased OSRW)
     * @param potential A CrystalPotential (either OSRW or underlying).
     * @param x Coordinates.
     * @param g Array to add gradients to.
     * @param tableIndex 0 for unbiased potential, 1 for a biased potential.
     * @return The generated EnergyResult.
     */
    private EnergyResult testGradientSet(String description, CrystalPotential potential, double[] x, double[] g, int tableIndex) {
        assertTrue(String.format(" Potential %s is not a lambda interface!", potential), potential instanceof LambdaInterface);

        EnergyResult er = new EnergyResult(description, potential, x, g);

        checkThGradScalar(er.energy, pe, tableIndex, peTol, "potential energy");
        checkThGradScalar(er.firstLam, dudl, tableIndex, dudlTol, "dU/dL");
        checkThGradArray(er.gradient, dudx, tableIndex, dudxTol, "dU/dX gradient");
        if (er.hasSecondDerivatives) {
            checkThGradScalar(er.secondLam, d2udl2, tableIndex, d2udl2Tol, "d2U/dL2");
            checkThGradArray(er.lamGradient, d2udxdl, tableIndex, d2udxdlTol, "d2U/dXdL");
        }

        return er;
    }

    /**
     * Checks a scalar value against its expected value. Mostly a convenience formatting method.
     *
     * @param actual Value from the test.
     * @param expected Array of expected values (see tableIndex).
     * @param tableIndex 0 for unbiased potential, 1 for biased OSRW potential.
     * @param tol Tolerance for this test.
     * @param description Scalar to be tested.
     */
    private void checkThGradScalar(double actual, double[] expected, int tableIndex, double tol, String description) {
        if (debugMode) {
            logger.info(String.format(" %s is %15.9g", description, actual));
        } else {
            assertEquals(String.format(" Expected %s %12.6g, received %12.6g from test %s",
                    description, expected[tableIndex], actual, info), expected[tableIndex], actual, tol);
        }
    }

    /**
     * Checks an array value (generally 1-D flat) against its expected value (generally 3-D; indices
     * tableIndex, then atom number, then X/Y/Z).
     *
     * @param actual Array from the test, flat.
     * @param expected Array of expected values, indexed by tableIndex, atoms, and XYZ.
     * @param tableIndex 0 for unbiased potential, 1 for biased OSRW potential.
     * @param tol Tolerance for this test.
     * @param description Array to be tested.
     */
    private void checkThGradArray(double[] actual, double[][][] expected, int tableIndex, double tol, String description) {
        double[] actualSlice = new double[3];
        for (int i = 0; i < numGradAtoms; i++) {
            int i3 = i*3;
            System.arraycopy(actual, i3, actualSlice, 0, 3);
            if (debugMode) {
                logger.info(String.format(" %s at atom %d is %s", description, i, Arrays.toString(actualSlice)));
            } else {
                double[] exp = expected[tableIndex][i];
                assertArrayEquals(String.format(" Discrepancy found on array of %s for test %s on atom %d." +
                                "\n Expected: %s\n Found: %s", description, info, i, Arrays.toString(exp), Arrays.toString(actualSlice)),
                        exp, actualSlice, tol);
            }
        }
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

    /**
     * Contains the result of an energy evaluation: potential energy and several derivatives.
     */
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
            gradient = Arrays.copyOf(g, nVars);

            hasSecondDerivatives = !(potential instanceof TransitionTemperedOSRW);
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
         * Asserts that two energy results are equivalent to each other. Always tests U, dU/dX, and dU/dL.
         * Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
         *
         * All values must be identical to tolerance.
         *
         * @param other Another EnergyResult
         */
        public void assertResultsEqual(EnergyResult other) {
            assertEquals(String.format(" Test %s: potential energy for %s did not " +
                            "match %s, %12.6g != %12.6g.", info, this.toString(), other.toString(), this.energy, other.energy),
                    this.energy, other.energy, peTol);
            assertEquals(String.format(" Test %s: dU/dL for %s did not " +
                            "match %s, %12.6g != %12.6g.", info, this.toString(), other.toString(), this.firstLam, other.firstLam),
                    this.firstLam, other.firstLam, dudlTol);
            assertArrayEquals(String.format(" Test %s: dU/dX for %s did not match %s.",
                    info, this.toString(), other.toString()),
                    this.gradient, other.gradient, dudxTol);

            if (hasSecondDerivatives && other.hasSecondDerivatives) {
                assertEquals(String.format(" Test %s: d2U/dL2 for %s did not " +
                                "match %s, %12.6g != %12.6g.", info, toString(), other.toString(), this.secondLam, other.secondLam),
                        this.secondLam, other.secondLam, d2udl2Tol);
                assertArrayEquals(String.format(" Test %s: d2U/dXdL for %s did not match %s.",
                        info, this.toString(), other.toString()),
                        this.lamGradient, other.lamGradient, d2udxdlTol);
            }
        }

        /**
         * Asserts that two energy results are not equivalent to each other. Always tests U, dU/dX, and dU/dL.
         * Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
         *
         * Prints all equalities, fails if no inequalities are found. At least one value must not be identical to tolerance.
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
                sb.append(String.format("\n Potential energy %12.6g == other potential energy %12.6g", energy, other.energy));
            } else {
                ++diffsFound;
            }

            if (approxEquals(firstLam, other.firstLam, dudlTol)) {
                equalFound = true;
                sb.append(String.format("\n Lambda derivative %12.6g == other lambda derivative %12.6g", firstLam, other.firstLam));
            } else {
                ++diffsFound;
            }

            for (int i = 0; i < nVars; i++) {
                if (approxEquals(gradient[i], other.gradient[i], dudxTol)) {
                    equalFound = true;
                    sb.append(String.format("\n Gradient %d %12.6g == other gradient %12.6g", i, gradient[i], other.gradient[i]));
                } else {
                    ++diffsFound;
                }
            }

            if (hasSecondDerivatives && other.hasSecondDerivatives) {
                if (approxEquals(secondLam, other.secondLam, d2udl2Tol)) {
                    equalFound = true;
                    sb.append(String.format("\n Lambda derivative %12.6g == other lambda derivative %12.6g", secondLam, other.secondLam));
                } else {
                    ++diffsFound;
                }

                for (int i = 0; i < nVars; i++) {
                    if (approxEquals(lamGradient[i], other.lamGradient[i], d2udxdlTol)) {
                        equalFound = true;
                        sb.append(String.format("\n Lambda gradient %d %12.6g == other lambda gradient %12.6g", i, lamGradient[i], other.lamGradient[i]));
                    } else {
                        ++diffsFound;
                    }
                }
            }

            if (equalFound) {
                logger.info(sb.toString());
            }
            assertTrue(String.format(" No inequalities found between %s and %s for test %s!",
                    this.toString(), other.toString(), info), diffsFound > 0);
        }

        @Override
        public String toString() {
            return description;
        }
    }

    private enum ThermoTestMode {
        HELP, FREE, GRAD
    }
}
