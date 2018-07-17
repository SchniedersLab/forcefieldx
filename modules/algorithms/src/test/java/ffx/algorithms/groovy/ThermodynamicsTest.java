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

import ffx.algorithms.groovy.NewThermodynamics;
import static ffx.potential.utils.PotentialsFunctions.logger;

import groovy.lang.Binding;
import static org.junit.Assert.*;

import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.MapConfiguration;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
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
public class ThermodynamicsTest {

/*


    public ThermodynamicsTest(String info, String[] filenames,
                              double pe, double dudl, double d2udl2,
                              double[] dudx, double[] d2udxdl,
                              double freeEnergy, double feTol,
                              String[] options, String[] properties, String[] flags) {
 */
    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        // Test info, file names, expected initial energy, dU/dL, d2U/dL2, dU/dX gradient, d2udxdl gradient,
        // expected free energy, tolerance for free energy, Groovy options, properties, Groovy boolean flags.
        return Arrays.asList(new Object[][]{
                {
                        "Thermodynamics Help Test", new String[]{}, 0, 0, 0, new double[]{}, new double[]{}, 0, 0, new String[]{}, new String[]{}, new String[]{"-h", "true"}
                }
        });
    }

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

    static {
        String[] opts = {"--bM", "0.05",
                "--dw", "OFF",
                "--lf", "1.0E-18",
                "--lm", "1.0E-18",
                "--np", "1",
                "--sf", "1.0",
                "--tp", "2.0",      // Set low to encourage fast tempering for toy cases.
                "-b", "ADIABATIC",  // Used for Langevin dynamics.
                "-C", "2",          // By default, evaluate 0 dynamics steps. Evaluate energy three times, should be different on the third step.
                "-d", "0.5",        // Most of the tests are in vacuum, thus short timestep.
                "-F", "XYZ",
                "-i", "STOCHASTIC", // Stochastic because most toy systems work best with Langevin dynamics.
                "-k", "1.0",
                "-l", "1.0",        // Set walker to start at the usually-physical topology. Likely commonly set by test.
                "-n", "0",          // Do not run any steps by default, just evaluate energy.
                "-p", "0",          // Use NVT, not NPT dynamics.
                "-Q", "0",          // Do not run any steps by default, just evaluate energy.
                "-r", "0.25",
                "-t", "298.15",
                "-w", "10.0"};

        int nOpts = opts.length;
        Map<String, String> optMap = new HashMap<>(nOpts);
        for (int i = 0; i < nOpts; i += 2) {
            optMap.put(opts[i], opts[i+1]);
        }
        DEFAULT_OPTIONS = Collections.unmodifiableMap(optMap);

        String[] props = {"ttosrw-alwaystemper", "true",
                "ttosrw-temperOffset", "0.5", // Small to temper fast.
                "print-on-failure", "false",
                "disable-neighbor-updates", "false"};

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
    private final String[] filenames;
    private final Map<String, String> opts;
    /**
     * Configuration containing the properties to be used by TT-OSRW.
     */
    Configuration algorithmConfig;
    private final List<String> flags;

    private final double pe;
    private final double dudl;
    private final double d2udl2;
    private final boolean runFreeEnergy;
    private final double freeEnergy;
    private final double feTol;

    private final int nVars;
    
    /**
     * Set of atoms for which tabulated gradients (dU/dL, d2U/dXdL) should be evaluated.
     */
    // private final int[] gradientAtoms;
    // private final boolean tabulatedGradients;
    private final double[] dudx;
    private final double[] d2udxdl;

    private final boolean onlyHelp;

    public ThermodynamicsTest(String info, String[] filenames,
                              double pe, double dudl, double d2udl2,
                              double[] dudx, double[] d2udxdl,
                              double freeEnergy, double feTol,
                              String[] options, String[] properties, String[] flags) {
        this.info = info;
        int nFiles = filenames.length;
        this.filenames = Arrays.copyOf(filenames, nFiles);

        assertTrue(String.format("Must have 0, 1, 2, or 4 distinct filenames, found %d for test %s", nFiles, info),
                nFiles < 3 || nFiles == 4);
        for (int i = 0; i < nFiles; i++) {
            for (int j = i+1; j < nFiles; j++) {
                assertNotEquals(String.format(" Filenames %d and %d matched in test %s: files %s and %s", i, j, info, filenames[i], filenames[j]),
                        filenames[i], filenames[j]);
            }
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

        this.pe = pe;
        this.dudl = dudl;
        this.d2udl2 = d2udl2;
        this.nVars = dudx.length; // Later, assert vs. the CrystalPotential.

        this.dudx = Arrays.copyOf(dudx, nVars);
        this.d2udxdl = Arrays.copyOf(d2udxdl, nVars);

        this.freeEnergy = freeEnergy;
        this.feTol = feTol;
        int nSteps = Integer.parseInt(opts.get("-n"));
        this.runFreeEnergy = (nSteps > 0);

        this.onlyHelp = addedFlags.getOrDefault("-h", false);
        assertTrue(String.format("Test %s lacks files, and is not a test of the help message!", info), onlyHelp || nFiles > 0);
    }

    @Before
    public void before() {
        binding = new Binding();
        thermo = new NewThermodynamics();
        thermo.setBinding(binding);
    }

    @Test
    public void testThermodynamicsHelp() {
        if (onlyHelp) {
            String[] args = {"-h"};
            binding.setVariable("args", args);
            thermo.run();
        }
    }
}
