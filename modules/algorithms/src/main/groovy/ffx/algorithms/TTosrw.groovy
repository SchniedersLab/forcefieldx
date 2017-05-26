
package ffx.algorithms;

// Java Imports
import java.util.regex.Pattern;
import java.util.stream.Collectors;

// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;
import groovy.util.CliBuilder;

// Parallel Java Imports
import edu.rit.pj.Comm;

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

// FFX Imports
import ffx.algorithms.Barostat;
import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.AlgorithmUtils;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.TransitionTemperedOSRW;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.RotamerOptimization;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.SymOp;

import ffx.numerics.Potential;
import ffx.numerics.PowerSwitch;
import ffx.numerics.SquaredTrigSwitch;
import ffx.numerics.UnivariateSwitchingFunction;

import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.OctTopologyEnergy;
import ffx.potential.QuadTopologyEnergy;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;

import ffx.potential.nonbonded.MultiplicativeSwitch;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;

/**
 * The TTosrw script uses the Transition-Tempered Orthogonal Space Random Walk
 * algorithm to estimate a free energy.
 * <br>
 * Usage:
 * <br>
 * ffxc TTosrw [options] &lt;filename [file2...]&gt;
 */
class TTosrw extends Script {

    /**
     * Options for the TTosrw Script.
     * <br>
     * Usage:
     * <br>
     * ffxc TTosrw [options] &lt;filename [file2...]&gt;
     */
    class Options {

        /*private final Closure parseThermostat = {str ->
        try {
        return Thermostats.valueOf(str.toUpperCase());
        } catch (Exception e) {
        logger.warning(String.format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
        return Thermostats.BERENDSEN;
        }
        };

        private final Closure parseIntegrator = {str ->
        try {
        return Integrators.valueOf(str.toUpperCase());
        } catch (Exception e) {
        logger.warning(String.format(" Could not parse %s as an integrator; defaulting to Beeman.", str));
        return Integrators.BEEMAN;
        }
        }*/

        private static parseThermo(String str) {
            try {
                return Thermostats.valueOf(str.toUpperCase());
            } catch (Exception e) {
                //logger.warning(String.format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
                System.stderr.println(String.format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
                return Thermostats.BERENDSEN;
            }
        }

        private static parseIntegrator(String str) {
            try {
                return Integrators.valueOf(str.toUpperCase());
            } catch (Exception e) {
                //logger.warning(String.format(" Could not parse %s as an integrator; defaulting to Beeman.", str));
                System.stderr.println(String.format(" Could not parse %s as an integrator; defaulting to Beeman.", str));
                return Integrators.BEEMAN;
            }
        }

        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -mc or --monteCarlo to use the experimental MC-TTOSRW
         */
        @Option(shortName='mc', longName='monteCarlo', defaultValue='false', description='Use experimental MC-TTOSRW sampling') boolean mc;
        /**
         * -a or --async sets asynchronous walker communication (recommended)
         */
        @Option(shortName='a', defaultValue='false', description='Walker communication is asynchronous.') boolean async;
        /**
         * -o or --optimize saves low-energy snapshots discovered (only for single topology simulations).
         */
        @Option(shortName='o', defaultValue='false', description='Optimize and save low-energy snapshots.') boolean optimize;
        /**
         * -n or --steps sets the number of molecular dynamics steps (default is 1 nsec).
         */
        @Option(shortName='n', defaultValue='1000000', description='Number of molecular dynamics steps') int steps;
        /**
         * -q or --equilibrate sets the number of equilibration steps prior to
         * production ttOSRW counts begin.
         */
        @Option(shortName='q', longName='equilibrate', defaultValue='1000', description='Equilibration steps prior to OSRW counts.') int nEquil;
        /**
         * -d or --dt sets the timestep; 1 femtosecond by default. 2.0 is
         * suggested for the RESPA integrator; otherwise, increase only if you
         * know what you are doing.
         */
        @Option(shortName='d', defaultValue='1.0', description='Time discretization step (fsec)') double dt;
        /**
         * -r or --report sets the thermodynamics reporting frequency in picoseconds (0.1 psec default).
         */
        @Option(shortName='r', longName='report', defaultValue='0.25', description='Interfal to report thermodynamics (psec).') double report;
        /**
         * -w or --write sets snapshot save frequency in picoseconds (1.0 psec default).
         */
        @Option(shortName='w', longName='write', defaultValue='10.0', description='Interval to write out coordinates (psec).') double write;
        /**
         * -t or --temperature sets the simulation temperature (Kelvin).
         */
        @Option(shortName='t', longName='temperature', defaultValue='298.15', description='Temperature (Kelvin)') double temp;
        /**
         * -b or --thermostat sets the desired thermostat: current choices are Adiabatic, Berendsen, or Bussi.
         */
        //@Option(shortName='b', longName='thermostat', defaultValue='Berendsen', description='Thermostat: [Adiabatic / Berendsen / Bussi].') String tstat;
        @Option(shortName='b', longName='thermostat', convert = {s -> return parseThermo(s);}, defaultValue='Berendsen', description='Thermostat: [Adiabatic / Berendsen / Bussi].') Thermostats tstat;
        /**
         * -i or --integrator sets the desired integrator: current choices are Beeman, RESPA, or Stochastic (AKA Langevin dynamics).
         */
        //@Option(shortName='i', longName='integrator', defaultValue='Beeman', description='Integrator: [Beeman / Respa / Stochastic]') String integrator;
        @Option(shortName='i', longName='integrator', convert = {s -> return parseIntegrator(s);}, defaultValue='Beeman', description='Integrator: [Beeman / Respa / Stochastic]') Integrators integrator;
        /**
         * -s1 or --start1 defines the first softcored atom for the first topology.
         */
        @Option(shortName='s1', longName='start1', defaultValue='0', description='Starting ligand atom for 1st topology') int s1;
        /**
         * -s2 or --start2 defines the first softcored atom for the second topology.
         */
        @Option(shortName='s2', longName='start2', defaultValue='0', description='Starting ligand atom for 2nd topology') int s2;
        /**
         * -f1 or --final1 defines the last softcored atom for the first topology.
         */
        @Option(shortName='f1', longName='final1', defaultValue='-1', description='Final ligand atom for the 1st topology') int f1;
        /**
         * -f2 or --final2 defines the last softcored atom for the second topology.
         */
        @Option(shortName='f2', longName='final2', defaultValue='-1', description='Final ligand atom for the 2nd topology') int f2;
        /**
         * --la1 or -ligAtoms1 allows for multiple ranges and/or singletons of ligand atoms in the first topology, separated by periods.
         */
        @Option(shortName='la1', longName='ligAtoms1', description='Period-separated ranges of 1st toplogy ligand atoms (e.g. 40-50.72-83)') String ligAt1;
        /**
         * --la2 or -ligAtoms2 allows for multiple ranges and/or singletons of ligand atoms in the second topology, separated by periods.
         */
        @Option(shortName='la2', longName='ligAtoms2', description='Period-separated ranges of 2nd toplogy ligand atoms (e.g. 40-50.72-83)') String ligAt2;
        /**
         * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
         */
        @Option(shortName='es1', longName='noElecStart1', defaultValue='1', description='Starting no-electrostatics atom for 1st topology') int es1;
        /**
         * -es2 or --noElecStart2 defines the first atom of the second topology to have no electrostatics.
         */
        @Option(shortName='es2', longName='noElecStart2', defaultValue='1', description='Starting no-electrostatics atom for 2nd topology') int es2;
        /**
         * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
         */
        @Option(shortName='ef1', longName='noElecFinal1', defaultValue='-1', description='Final no-electrostatics atom for 1st topology') int ef1;
        /**
         * -ef2 or --noElecFinal2 defines the last atom of the second topology to have no electrostatics.
         */
        @Option(shortName='ef2', longName='noElecFinal2', defaultValue='-1', description='Final no-electrostatics atom for 2nd topology') int ef2;
        /**
         * -l or --lambda sets the lambda value to minimize at.
         */
        @Option(shortName='l', longName='lambda', defaultValue='-1',
            description='Initial lambda value (> 1.0 distributes lambda across walkers)') double lambda;
        /**
         * -c or --count sets the number of time steps between OSRW counts.
         */
        @Option(shortName='c', longName='count', defaultValue='10',
            description='Time steps between OSRW counts.') int countFreq;
        /**
         * -b or --bias sets the initial Gaussian bias magnitude in kcal/mol.
         */
        @Option(shortName='g', longName='bias', defaultValue='0.002',
            description='Gaussian bias magnitude (kcal/mol)') double biasMag;
        /**
         * -m or --mass sets the lambda particle mass.
         */
        @Option(shortName='m', longName='mass', defaultValue='1.0e-18',
            description='Lambda particle mass') double lamMass;
        /**
         * -x or --friction sets the friction on the lambda particle.
         */
        @Option(shortName='x', longName='friction', defaultValue='1.0e-18',
            description='Lambda particle friction') double lamFric;
        /**
         * -p or --npt Specify use of a MC Barostat at the given pressure (default 1.0 atm).
         */
        @Option(shortName='p', longName='npt', defaultValue='0',
            description='Specify use of a MC Barostat at the given pressure (default of 0 = disabled)') double pressure;
        /**
         * -sym or --symOp to apply a random Cartesian symmetry operator with the specified translation range -X .. X (no default).
         */
        @Option(shortName='rsym', longName='randomSymOp', defaultValue='-1.0',
            description='Apply a random Cartesian symmetry operator with a random translation in the range -X .. X.') double symScalar
        /**
         * -ruc or --unitCell random unit cell axes will be used achieve the specified density (g/cc) (no default density).
         */
        @Option(shortName='ruc', longName='randomUnitCell', defaultValue='-1.0',
            description='Apply random unit cell axes to achieve the specified density (g/cc).') double ucDensity
        /**
         * -ld or --minDensity sets a tin box constraint on the barostat, preventing over-expansion of the box (particularly in vapor phase), permitting an analytic correction.
         */
        @Option(shortName='ld', longName='minDensity', defaultValue='0.75',
            description='Minimum density allowed by the barostat') double minDensity;
        /**
         * -hd or --maxDensity sets a maximum density on the barostat, preventing under-expansion of the box.
         */
        @Option(shortName='hd', longName='maxDensity', defaultValue='1.6',
            description='Maximum density allowed by the barostat') double maxDensity;
        /**
         * -sm or --maxSideMove sets the width of proposed crystal side length moves (rectangularly distributed) in Angstroms.
         */
        @Option(shortName='sm', longName='maxSideMove', defaultValue='0.25',
            description='Maximum side move allowed by the barostat in Angstroms') double maxSideMove;
        /**
         * -am or --maxAngleMove sets the width of proposed crystal angle moves (rectangularly distributed) in degrees.
         */
        @Option(shortName='am', longName='maxAngleMove', defaultValue='0.5',
            description='Maximum angle move allowed by the barostat in degrees') double maxAngleMove;
        /**
         * -mi or --meanInterval sets the mean number of MD steps (Poisson distribution) between barostat move proposals.
         */
        @Option(shortName='mi', longName='meanInterval', defaultValue='10',
            description='Mean number of MD steps between barostat move proposals.') int meanInterval;

        /**
         * -W or --traversals sets writing lambda traversal snapshots.
         */
        @Option(shortName='W', longName='traversals', defaultValue='false',
            description='Write out lambda-traversal snapshots.') boolean traversals;
        /**
         * -rt or --reset resets the OSRW histogram once, at lambda 0.99.
         */
        @Option(shortName='rt', longName='reset', defaultValue='false',
            description='Reset OSRW histogram once, when lambda reaches 0.99.') boolean reset;
        /**
         * -tp or --temperingParam sets the Dama et al tempering rate parameter,
         * in multiples of kBT.
         */
        @Option(shortName='tp', longName='temperingParam', defaultValue='4.0',
            description='Dama et al tempering rate parameter in multiples of kBT') double temperParam;
        /**
         * -rn or --resetNumSteps, ignores steps detected in .lam lambda-restart
         * files and thus resets the histogram; use -rn false to continue from
         * the end of any prior simulation.
         */
        @Option(shortName='rn', longName='resetNumSteps', defaultValue='true',
            description='Ignore prior steps logged in .lam files') String resetStepsString;
        /**
         * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, 4, or 8.
         */
        @Option(shortName='np', longName='nParallel', defaultValue='1',
            description='Number of topologies to evaluate in parallel') int nPar;
        /**
         * -uaA or --unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaA', longName='unsharedA',
            description='Unshared atoms in the A dual topology (period-separated hyphenated ranges)') String unsharedA;
        /**
         * -uaB or --unsharedB sets atoms unique to the B dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaB', longName='unsharedB',
            description='Unshared atoms in the B dual topology (period-separated hyphenated ranges)') String unsharedB;
        /**
         * -dw or --distributeWalkers allows walkers to start from multiple
         * conformations; AUTO picks up per-walker conformations as
         * filename.pdb_(walker number), and specifying a residue starts a
         * rotamer optimization to generate side-chain configurations to start
         * from.
         */
        @Option(shortName='dw', longName='distributeWalkers',
            description='AUTO: Pick up per-walker configurations as [filename.pdb]_[num], or specify a residue to distribute on.') String distWalksString;
        /**
         * -le or --lambdaExponent sets the power of lambda used by dual
         * topologies. Now deprecated by numeric argument to -sf.
         */

        @Option(shortName='le', longName='lambdaExponent', defaultValue='1.0',
            description='DEPRECATED: Exponent to apply to dual topology lambda.') double lamExp;
        /**
         * -sf or --switchingFunction sets the switching function to be used by
         * dual topologies; TRIG produces the function sin^2(pi/2*lambda)*E1(lambda)
         * + cos^2(pi/2*lambda)*E2(1-lambda), MULT uses a 5'th-order polynomial
         * switching function with zero first and second derivatives at the end
         * (same function as used for van der Waals switch), and a number uses
         * the original function, of l^beta*E1(lambda) + (1-lambda)^beta*E2(1-lambda).
         *
         * All of these are generalizations of Udt = f(l)*E1(l) +
         * f(1-l)*E2(1-lambda), where f(l) is a continuous switching function
         * such that f(0) = 0, f(1) = 1, and 0 <= f(l) <= 1 for lambda 0-1.
         * The trigonometric switch can be restated thusly, since
         * cos^2(pi/2*lambda) is identical to sin^2(pi/2*(1-lambda)), f(1-l).
         */
        @Option(shortName='sf', longName='switchingFunction', defaultValue='1.0',
            description='Switching function to use for dual topology: options are TRIG, MULT, or a number (original behavior with specified lambda exponent)') String lambdaFunction;


        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames;
    }

    private static final Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");
    private static double maxdUdL = 1000.0;

    private int threadsAvail = edu.rit.pj.ParallelTeam.getDefaultThreadCount();
    private int threadsPer = threadsAvail;
    private AlgorithmFunctions aFuncts;
    def ranges1 = []; // Groovy mechanism for creating an untyped ArrayList.
    def ranges2 = [];
    def rangesA = [];
    def rangesB = [];
    def topologies = []; // MolecularAssembly
    def properties = []; // CompositeConfiguration
    def energies = [];   // ForceFieldEnergy

    boolean autoDist = false;
    def distResidues = [];
    RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
    def final static nochainMatcher = ~/^\ ?([0-9]+)$/;
    def final static chainMatcher = ~/^([a-zA-Z])([0-9]+)$/;

    int size = 1;
    int rank = 0;

    /**
     * Opens a file and processes it for ttOSRW
     *
     * @param options TTosrw Options
     * @param toOpen Filename to open
     * @param structFile Structure file to use for coordinates
     * @param topNum Number of the topology to open
     */
    private void openFile(Options options, String toOpen, File structFile, int topNum) {
        if (autoDist) {
            String openName = String.format("%s_%d", toOpen, rank+1);
            File testFile = new File(openName);
            if (testFile.exists()) {
                toOpen = openName;
            } else {
                logger.warning(String.format(" File %s does not exist; using default %s", openName, toOpen));
            }
        }

        MolecularAssembly[] opened = aFuncts.openAll(toOpen, threadsPer);
        MolecularAssembly mola = aFuncts.getActiveAssembly();
        processFile(options, mola, structFile, topNum);
    }

    /**
     * Split off from openFile to ensure molecular assemblies obtained from the
     * UI are properly treated.
     *
     * @param options TTosrw Options
     * @param mola Molecular assembly to process
     * @param structFile Structure file to use for coordinates
     * @param topNum Number of the topology to open
     */
    private void processFile(Options options, MolecularAssembly mola, File structFile, int topNum) {
        if (size > 1) {
            mola.setFile(structFile);
        }
        ForceFieldEnergy energy = mola.getPotentialEnergy();
        Atom[] atoms = mola.getAtomArray();
        int remainder = (topNum % 2) + 1;
        switch(remainder) {
        case 1:
            /**
             * Improve this logic once @Option annotations are more finished.
             */
            if (options.s1 > 0) {
                for (int i = options.s1; i <= options.f1; i++) {
                    Atom ai = atoms[i-1];
                    ai.setApplyLambda(true);
                    ai.print();
                }
            }
            if (ranges1) {
                for (range in ranges1) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        // Don't need to worry about negative numbers; rangeregex just won't match.
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            Atom ai = atoms[i-1];
                            ai.setApplyLambda(true);
                            ai.print();
                        }
                    } else {
                        logger.warning(" Could not recognize ${range} as a valid range; skipping");
                    }
                }
            }

            // Apply the no electrostatics atom selection
            int noElecStart = options.es1;
            noElecStart = (noElecStart < 1) ? 1 : noElecStart;

            int noElecStop = options.ef1;
            noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop;

            for (int i = noElecStart; i <= noElecStop; i++) {
                Atom ai = atoms[i - 1];
                ai.setElectrostatics(false);
                ai.print();
            }
            break;
        case 2:
            /**
             * Improve this logic once @Option annotations are more finished.
             */
            if (options.s2 > 0) {
                for (int i = options.s2; i <= options.f2; i++) {
                    Atom ai = atoms[i-1];
                    ai.setApplyLambda(true);
                    ai.print();
                }
            }
            if (ranges2) {
                for (range in ranges2) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        // Don't need to worry about negative numbers; rangeregex just won't match.
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            Atom ai = atoms[i-1];
                            ai.setApplyLambda(true);
                            ai.print();
                        }
                    } else {
                        logger.warning(" Could not recognize ${range} as a valid range; skipping");
                    }
                }
            }

            // Apply the no electrostatics atom selection
            int noElecStart2 = options.es2;
            noElecStart2 = (noElecStart2 < 1) ? 1 : noElecStart2;

            int noElecStop2 = options.ef2;
            noElecStop2 = (noElecStop2 > atoms.length) ? atoms.length : noElecStop2;

            for (int i = noElecStart2; i <= noElecStop2; i++) {
                Atom ai = atoms[i - 1];
                ai.setElectrostatics(false);
                ai.print();
            }
            break;
        }

        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        energy.getCrystal().setSpecialPositionCutoff(0.0);
        // Save a reference to the topology.
        properties[topNum] = mola.getProperties();
        topologies[topNum] = mola;
        energies[topNum] = energy;
    }

    /**
     * Distribute side-chain conformations of mola.
     *
     * @param mola To distribute
     * @param pot Potential to use
     */
    private void optStructure(MolecularAssembly mola, Potential pot) {
        if (!distResidues) {
            throw new IllegalArgumentException(" Programming error: Must have list of residues to split on!");
        }

        LambdaInterface linter = (pot instanceof LambdaInterface) ? (LambdaInterface) pot : null;
        double initialLambda = linter ? linter.getLambda() : -1.0;
        linter?.setLambda(0.5);
        // Safe navigation operator ?. operates only if LHS is non-null.

        def resList = [];
        Polymer[] polymers = mola.getChains();

        for (String ts : distResidues) {
            Character chainID = 'A';
            def m = chainMatcher.matcher(ts);
            int resNum = -1;
            if (m.find()) {
                chainID = m.group(1);
                resNum = Integer.parseInt(m.group(2));
            } else {
                m = nochainMatcher.matcher(ts);
                if (m.find()) {
                    resNum = Integer.parseInt(m.group(1));
                } else {
                    logger.warning(String.format(" Could not parse %s as a valid residue!", ts));
                    continue;
                }
            }

            logger.info(String.format(" Looking for chain %c residue %d", chainID, resNum));

            for (Polymer p : mola.getChains()) {
                if (p.getChainID() == chainID) {
                    for (Residue r : p.getResidues()) {
                        if (r.getResidueNumber() == resNum && r.getRotamers(rLib) != null) {
                            resList.add(r);
                        }
                    }
                }
            }
        }

        if (!resList) {
            throw new IllegalArgumentException(" No valid entries for distWalkers!");
        }

        AlgorithmListener alist = aFuncts.getDefaultListener();
        RotamerOptimization ropt = new RotamerOptimization(mola, pot, alist);

        ropt.setThreeBodyEnergy(false);
        ropt.setVerboseEnergies(true);
        if (System.getProperty("ro-ensembleNumber") == null && System.getProperty("ro-ensembleEnergy") == null) {
            logger.info(String.format(" Setting ensemble to default of number of walkers %d", size));
            ropt.setEnsemble(size);
        }
        ropt.setPrintFiles(false);
        def addedResList = ropt.setResiduesIgnoreNull(resList);

        rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        rLib.setUseOrigCoordsRotamer(false);
        RotamerLibrary.measureRotamers(resList, false);

        String oldLazyMat = System.getProperty("ro-lazyMatrix");
        System.setProperty("ro-lazyMatrix", "true");

        ropt.optimize(RotamerOptimization.Algorithm.GLOBAL_DEE);
        ropt.setCoordinatesToEnsemble(rank);

        // One final energy call to ensure the coordinates are properly set at the
        // end of rotamer optimization.
        double[] xyz = new double[pot.getNumberOfVariables()];
        pot.getCoordinates(xyz);
        logger.info(" Final Optimized Energy:");
        pot.energy(xyz, true);

        linter?.setLambda(initialLambda);

        if (oldLazyMat) {
            System.setProperty("ro-lazyMatrix", oldLazyMat);
        } else {
            System.clearProperty("ro-lazyMatrix");
        }
    }

    def run() {
        def cli = new CliBuilder(usage:' ffxc TTosrw [options] <filename> [file2...]', header:' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }

        try {
            aFuncts = getAlgorithmUtils();
        } catch (MissingMethodException ex) {
            aFuncts = new AlgorithmUtils();
        }

        if (options.distWalksString) {
            if (options.distWalksString.equalsIgnoreCase("AUTO")) {
                autoDist = true;
            } else {
                distResidues = options.distWalksString.split("\\.");
            }
        }

        List<String> arguments = options.filenames;
        // Check nArgs; should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1;
        nArgs = (nArgs < 1) ? 1 : nArgs;

        int numParallel = options.nPar;
        if (threadsAvail % numParallel != 0) {
            logger.warning(String.format(" Number of threads available %d not evenly divisible by np %d; reverting to sequential", threadsAvail, numParallel));
            numParallel = 1;
        } else if (nArgs % numParallel != 0) {
            logger.warning(String.format(" Number of topologies %d not evenly divisible by np %d; reverting to sequential", arguments.size(), numParallel));
            numParallel = 1;
        } else {
            threadsPer = threadsAvail / numParallel;
        }

        if (options.ligAt1) {
            ranges1 = options.ligAt1.tokenize(".");
        }
        if (options.ligAt2) {
            ranges2 = options.ligAt2.tokenize(".");
        }

        MolecularAssembly fromUI = null;
        if (!arguments || arguments.isEmpty()) {
            fromUI = aFuncts.getActiveAssembly();
            if (fromUI == null) {
                return cli.usage();
            }
            arguments = new ArrayList<>();
            arguments.add(fromUI.getFile().getName());
            // Do not open the files yet!
        }

        CrystalPotential potential;

        UnivariateSwitchingFunction sf;
        if (options.lambdaFunction) {
            String lf = options.lambdaFunction.toUpperCase();
            switch (lf) {
            case ~/^-?[0-9]*\.?[0-9]+/:
                double exp = Double.parseDouble(lf);
                sf = new PowerSwitch(1.0, exp);
                break;
            case "TRIG":
                sf = new SquaredTrigSwitch(false);
                break;
            case "MULT":
                sf = new MultiplicativeSwitch(0.0, 1.0);
                break;
            default:
                try {
                    double beta = Double.parseDouble(lf);
                    sf = new PowerSwitch(1.0, beta);
                } catch (NumberFormatException ex) {
                    logger.warning(String.format("Argument to option -sf %s could not be properly parsed; using default linear switch", options.lambdaFunction));
                    sf = new PowerSwitch(1.0, 1.0);
                }
            }
        } else {
            sf = new PowerSwitch(1.0, options.lamExp);
        }

        String filename = arguments[0];
        File structureFile = new File(FilenameUtils.normalize(filename));
        structureFile = new File(structureFile.getAbsolutePath());
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
        File histogramRestart = new File(baseFilename + ".his");
        File lambdaOneFile = null;
        File lambdaZeroFile = null;
        File lambdaRestart = null;
        File dyn = null;

        Comm world = Comm.world();
        size = world.size();
        rank = 0;

        // For a multi-process job, try to get the restart files from rank sub-directories.
        String withRankName = baseFilename;
        if (size > 1) {
            rank = world.rank();
            File rankDirectory = new File(structureFile.getParent() + File.separator
                + Integer.toString(rank));
            if (!rankDirectory.exists()) {
                rankDirectory.mkdir();
            }
            withRankName = rankDirectory.getPath() + File.separator + baseFilename;
            /*lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
            dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
            lambdaOneFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam1");
            lambdaZeroFile = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam0");*/
            structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
        }/* else {
        // For a single process job, try to get the restart files from the current directory.
        lambdaRestart = new File(baseFilename + ".lam");
        dyn = new File(baseFilename + ".dyn");
        lambdaOneFile = new File(baseFilename + ".lam1");
        lambdaZeroFile = new File(baseFilename + ".lam0");
        }*/
        lambdaRestart = new File(withRankName + ".lam");
        dyn = new File(withRankName + ".dyn");
        lambdaOneFile = new File(withRankName + ".lam1");
        lambdaZeroFile = new File(withRankName + ".lam0");

        if (!dyn.exists()) {
            dyn = null;
        }

        // Turn on computation of lambda derivatives if softcore atoms exist or only a single topology.
        boolean lambdaTerm = options.ligAt1 || options.ligAt2 || (options.s1 > 0) || (options.s2 > 0) || nArgs == 1;
        if (lambdaTerm) {
            System.setProperty("lambdaterm","true");
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec","false");
        }

        double lambda = options.lambda;

        // Apply the command line lambda value if a lambda restart file does not exist.
        //if (!lambdaRestart.exists()) {
        if (lambda < 0.0 || lambda > 1.0) {
            if (size > 1) {
                //dL = 1.0 / (size - 1.0);
                //lambda = rank * dL;
                dL = 1.0 / (size + 1.0);
                lambda = dL * (rank + 1);
                if (lambda > 1.0) {
                    lambda = 1.0;
                }
                if (lambda < 0.0) {
                    lambda = 0.0;
                }
                logger.info(String.format(" Setting lambda to %5.3f.", lambda));
            } else {
                lambda = 0.5;
                logger.info(String.format(" Setting lambda to %5.3f", lambda));
            }
        }
        //}

        if (fromUI != null) {
            processFile(options, fromUI, structureFile, 0);
        } else {
            for (int i = 0; i < arguments.size(); i++) {
                openFile(options, arguments.get(i), structureFile, i);
            }
        }

        List<Integer> uniqueA;
        List<Integer> uniqueB;
        if (nArgs >= 4) {
            uniqueA = new ArrayList<>();
            uniqueB = new ArrayList<>();

            if (options.unsharedA) {
                def ra = [] as Set;
                String[] toksA = options.unsharedA.tokenize(".");
                Atom[] atA1 = topologies[0].getAtomArray();
                for (range in toksA) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        logger.info(String.format("Range %s for A, start %d end %d", range, rangeStart, rangeEnd));
                        logger.info(String.format(" First atom in range: %s", atA1[rangeStart-1]));
                        if (rangeEnd > rangeStart) {
                            logger.info(String.format(" Last atom in range: %s", atA1[rangeEnd-1]));
                        }
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            ra.add(i-1);
                        }
                    }
                }
                int counter = 0;
                def raAdj = [] as Set; // Indexed by common variables in dtA.
                for (int i = 0; i < atA1.length; i++) {
                    Atom ai = atA1[i];
                    if (i in ra) {
                        if (ai.applyLambda()) {
                            logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                        } else {
                            logger.fine(String.format(" Unshared A: %d variables %d-%d", i, counter, counter+2));
                            for (int j = 0; j < 3; j++) {
                                raAdj.add(new Integer(counter + j));
                            }
                        }
                    }
                    if (! ai.applyLambda()) {
                        counter += 3;
                    }
                }
                if (raAdj) {
                    uniqueA.addAll(raAdj);
                }
            }
            if (options.unsharedB) {
                def rb = [] as Set;
                String[] toksB = options.unsharedB.tokenize(".");
                Atom[] atB1 = topologies[2].getAtomArray();
                for (range in toksB) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        logger.info(String.format("Range %s for B, start %d end %d", range, rangeStart, rangeEnd));
                        logger.info(String.format(" First atom in range: %s", atB1[rangeStart-1]));
                        if (rangeEnd > rangeStart) {
                            logger.info(String.format(" Last atom in range: %s", atB1[rangeEnd-1]));
                        }
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            rb.add(i-1);
                        }
                    }
                }
                int counter = 0;
                def rbAdj = [] as Set; // Indexed by common variables in dtA.
                for (int i = 0; i < atB1.length; i++) {
                    Atom bi = atB1[i];
                    if (i in rb) {
                        if (bi.applyLambda()) {
                            logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                        } else {
                            logger.fine(String.format(" Unshared B: %d variables %d-%d", i, counter, counter+2));
                            for (int j = 0; j < 3; j++) {
                                rbAdj.add(counter + j);
                            }
                        }
                    }
                    if (! bi.applyLambda()) {
                        counter += 3;
                    }
                }
                if (rbAdj) {
                    uniqueB.addAll(rbAdj);
                }
            }
        }

        TransitionTemperedOSRW osrw = null;
        def dualTopologies = []; // Used for distResidues on quad/oct topologies
        StringBuilder sb = new StringBuilder("\n Running Transition-Tempered Orthogonal Space Random Walk for ");
        switch (nArgs) {
        case 1:
            if (options.symScalar > 0.0) {
                SymOp symOp = SymOp.randomSymOpFactory(options.symScalar);
                logger.info(String.format("\n Applying random Cartesian SymOp:\n%s", symOp.toString()));
                Crystal crystal = topologies[0].getCrystal();
                Atom[] atoms = topologies[0].getAtomArray();
                double[] xyz = new double[3];
                for (int i=0; i<atoms.length; i++) {
                    atoms[i].getXYZ(xyz);
                    crystal.applyCartesianSymOp(xyz, xyz, symOp);
                    atoms[i].setXYZ(xyz);
                }
            }

            if (options.ucDensity > 0.0) {
                logger.info(String.format("\n Applying random unit cell axes with target density of %6.3f\n",
                        options.ucDensity));
                Crystal crystal = topologies[0].getCrystal();
                if (!crystal.aperiodic()) {
                    double mass = topologies[0].getMass();
                    crystal.randomParameters(options.ucDensity, mass);
                    energies[0].setCrystal(crystal);
                }
            }
            potential = topologies[0].getPotentialEnergy();
            break;
        case 2:
            sb.append("dual topology ");
            DualTopologyEnergy dte = new DualTopologyEnergy(topologies[0], topologies[1], sf);
            if (numParallel == 2) {
                dte.setParallel(true);
            }
            potential = dte;
            break;
        case 4:
            sb.append("quad topology ");

            DualTopologyEnergy dta = new DualTopologyEnergy(topologies[0], topologies[1], sf);
            DualTopologyEnergy dtb = new DualTopologyEnergy(topologies[3], topologies[2], sf);
            QuadTopologyEnergy qte = new QuadTopologyEnergy(dta, dtb, uniqueA, uniqueB);
            if (numParallel >= 2) {
                qte.setParallel(true);
                if (numParallel == 4) {
                    dta.setParallel(true);
                    dtb.setParallel(true);
                }
            }
            potential = qte;
            dualTopologies[0] = dta;
            dualTopologies[1] = dtb;
            break;
        case 8:
            sb.append("oct-topology ");

            DualTopologyEnergy dtga = new DualTopologyEnergy(topologies[0], topologies[1], sf);
            DualTopologyEnergy dtgb = new DualTopologyEnergy(topologies[3], topologies[2], sf);
            QuadTopologyEnergy qtg = new QuadTopologyEnergy(dtga, dtgb, uniqueA, uniqueB);

            DualTopologyEnergy dtda = new DualTopologyEnergy(topologies[4], topologies[5], sf);
            DualTopologyEnergy dtdb = new DualTopologyEnergy(topologies[7], topologies[6], sf);
            QuadTopologyEnergy qtd = new QuadTopologyEnergy(dtda, dtdb, uniqueA, uniqueB);

            OctTopologyEnergy ote = new OctTopologyEnergy(qtg, qtd, true);
            if (numParallel >= 2) {
                ote.setParallel(true);
                if (numParallel >= 4) {
                    qtg.setParallel(true);
                    qtd.setParallel(true);
                    if (numParallel == 8) {
                        dtga.setParallel(true);
                        dtgb.setParallel(true);
                        dtda.setParallel(true);
                        dtdb.setParallel(true);
                    }
                }
            }
            potential = ote;
            dualTopologies[0] = dtga;
            dualTopologies[1] = dtgb;
            dualTopologies[2] = dtda;
            dualTopologies[3] = dtdb;
            break;
        default:
            logger.severe(" Must have 1, 2, 4, or 8 topologies!");
            break;
        }
        sb.append(topologies.stream().map{t -> t.getFile().getName()}.collect(Collectors.joining(",", "[", "]")));
        sb.append("\n");
        logger.info(sb.toString());

        logger.info(" Starting energy (before .dyn restart loaded):");
        boolean updatesDisabled = topologies[0].getForceField().getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (updatesDisabled) {
            logger.info(" This ensures neighbor list is properly constructed from the source file, before coordinates updated by .dyn restart");
        }
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);
        LambdaInterface linter = (LambdaInterface) potential;
        linter.setLambda(lambda);

        potential.energy(x, true);

        if (distResidues) {
            logger.info(" Distributing walker conformations.");
            switch (nArgs) {
            case 1:
                optStructure(topologies[0], energies[0]);
                break;
            case 2:
                if (potential.getNumSharedVariables() == potential.getNumberOfVariables()) {
                    logger.info(" Generating starting structures based on dual-topology:");
                    optStructure(topologies[0], potential);
                } else {
                    logger.info(" Generating separate starting structures for each topology of the dual toplogy:");
                    optStructure(topologies[0], energies[0]);
                    optStructure(topologies[1], energies[1]);
                }
                break;
            case 4:
                optStructure(topologies[0], dualTopologies[0]);
                optStructure(topologies[3], dualTopologies[1]);
                break;
            case 8:
                optStructure(topologies[0], dualTopologies[0]);
                optStructure(topologies[3], dualTopologies[1]);

                // More elegant would be to copy coordinates from 0-4 and 3-7.
                // More elegant is currently low priority, and thus will be accomplished as t approaches infinity.
                optStructure(topologies[4], dualTopologies[2]);
                optStructure(topologies[7], dualTopologies[3]);
                break;
            default:
                logger.severe(" First: must have 1, 2, 4, or 8 topologies. Second, how did this script not fail earlier?");
                break;
            }
        }

        boolean resetNumSteps = true;
        if (options.resetStepsString) {
            if (options.nEquil > 0) {
                logger.info(" Ignoring resetNumSteps input due to equilibration");
            } else if (options.resetStepsString.equalsIgnoreCase("false")) {
                resetNumSteps = false;
            }
        }

        osrw = new TransitionTemperedOSRW(potential, potential, lambdaRestart, histogramRestart,
            topologies[0].getProperties(), options.temp, options.dt, options.report,
            options.write, options.async, resetNumSteps, aFuncts.getDefaultListener());

        osrw.setResetStatistics(options.reset);
        if (options.traversals) {
            if (nArgs == 1) {
                osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[0]);
            } else if (nArgs == 2) {
                osrw.setTraversalOutput(lambdaOneFile, topologies[0], lambdaZeroFile, topologies[1]);
            }
        }
        osrw.setDeltaT(options.temperParam);

        if (options.optimize) {
            osrw.setOptimization(true, topologies[0]);
        }

        if (!lambdaRestart.exists()) {
            logger.info(String.format(" Setting lambda to %5.3f", lambda));
            osrw.setLambda(lambda);
        }

        // Apply the command line OSRW values if a histogram restart file does not exist.
        if (!histogramRestart.exists()) {
            osrw.setThetaFrication(options.lamFric);
            osrw.setThetaMass(options.lamMass);
            osrw.setCountInterval(options.countFreq);
            osrw.setBiasMagnitude(options.biasMag);
        }

        if (options.pressure) {
            Barostat barostat = new Barostat(topologies[0], osrw);
            barostat.setPressure(options.pressure);
            barostat.setMaxDensity(options.maxDensity);
            barostat.setMinDensity(options.minDensity);
            double dens = barostat.density();
            if (dens < options.minDensity || dens > options.maxDensity) {
                barostat.setDensity(1.0);
            }
            barostat.setMaxSideMove(options.maxSideMove);
            barostat.setMaxAngleMove(options.maxAngleMove);
            barostat.setMeanBarostatInterval(options.meanInterval);
            potential = barostat;
        } else {
            potential = osrw;
        }

        if (options.mc) {
            MonteCarloOSRW mcOSRW = new MonteCarloOSRW(osrw.getPotentialEnergy(), osrw, topologies[0],
                topologies[0].getProperties(), null, Thermostats.ADIABATIC, Integrators.VELOCITYVERLET);
            mcOSRW.sample();
        } else {
            // Create the MolecularDynamics instance.
            MolecularDynamics molDyn = new MolecularDynamics(topologies[0], potential,
                topologies[0].getProperties(), null, options.tstat, options.integrator);
            for (int i = 1; i < topologies.size(); i++) {
                molDyn.addAssembly(topologies.get(i), properties.get(i));
            }

            boolean initVelocities = true;
            double restartInterval = 0.1;
            String fileType = "XYZ";
            int nSteps = options.steps;
            // Start sampling.
            if (options.nEquil > 0) {
                logger.info(" Beginning equilibration");
                osrw.setPropagateLambda(false);
                molDyn.dynamic(options.nEquil, options.dt, options.report, options.write, options.temp, initVelocities, dyn);
                logger.info(" Beginning Transition-Tempered OSRW sampling");
                osrw.setPropagateLambda(true);
                molDyn.dynamic(nSteps, options.dt, options.report, options.write, options.temp, false,
                    fileType, restartInterval, dyn);
            } else {
                logger.info(" Beginning Transition-Tempered OSRW sampling without equilibration");
                boolean resetSteps = true;
                if (options.resetStepsString) {
                    if (options.resetStepsString.equalsIgnoreCase("false")) {
                        resetSteps = false;
                    }
                }
                if (!resetSteps) {
                    int nEnergyCount = osrw.getEnergyCount();
                    if (nEnergyCount > 0) {
                        nSteps -= nEnergyCount;
                        logger.info(String.format(" Lambda file: %12d steps picked up, now sampling %12d steps", nEnergyCount, nSteps));
                        initVelocities = false;
                    }
                }
                if (nSteps > 0) {
                    molDyn.dynamic(nSteps, options.dt, options.report, options.write, options.temp, initVelocities,
                        fileType, restartInterval, dyn);
                } else {
                    logger.info(" No steps remaining for this process!");
                }
            }
        }
    }
}

