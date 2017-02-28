
package ffx.algorithms;

// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;
import groovy.util.CliBuilder;

// Java Imports
import java.util.regex.Pattern;
import java.util.stream.Collectors;

// FFX Imports
import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.AlgorithmUtils;
import ffx.algorithms.Minimize;

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
import ffx.potential.nonbonded.MultiplicativeSwitch;

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;

/**
 * The Minimize script uses a limited-memory BFGS algorithm to minimize the 
 * energy of a molecular system or supersystem.
 * <br>
 * Usage:
 * <br>
 * ffxc Minimizer [options] &lt;filename [file2...]&gt;
 */
class Minimizer extends Script {

    /**
     * Options for the Minimize Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Minimize [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
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
        @Option(shortName='l', longName='lambda', defaultValue='-1', description='Lambda value') double lambda;
        /**
         * -e or --eps sets the RMS gradient convergence criteria (in kcal/mol/A^2)
         */
        @Option(shortName='e', longName='eps', defaultValue='1.0', description='RMS gradient convergence criteria') double eps;
        /**
         * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, or 4.
         */
        @Option(shortName='np', longName='nParallel', defaultValue='1', description='Number of topologies to evaluate in parallel') int nPar;
        /**
         * -uaA or --unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaA', longName='unsharedA', description='Unshared atoms in the A dual topology (period-separated hyphenated ranges)') String unsharedA;
        //@Option(shortName='uaA', longName='unsharedA', numberOfArgumentsString='*', valueSeparator='.', description='Unshared atoms in the A dual topology (period-separated hyphenated ranges)') String[] unsharedA;
        /**
         * -uaB or --unsharedB sets atoms unique to the B dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaB', longName='unsharedB', description='Unshared atoms in the B dual topology (period-separated hyphenated ranges)') String unsharedB;
        //@Option(shortName='uaB', longName='unsharedB', numberOfArgumentsString='*', valueSeparator='.', description='Unshared atoms in the B dual topology (period-separated hyphenated ranges)') String[] unsharedB;
        // I'm not sure it's worthwhile to clutter up the options with something that can be set by property.
        //@Option(shortName='p', longName='polarization', defaultValue='MUTUAL', description='Polarization model: [none/direct/mutual]') String polarizationModel;
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

    private void openFile(Options options, String toOpen, int topNum) {
        MolecularAssembly[] opened = aFuncts.openAll(toOpen, threadsPer);
        MolecularAssembly mola = aFuncts.getActiveAssembly();
        processFile(options, mola, topNum);
    }
    
    private void processFile(Options options, MolecularAssembly mola, int topNum) {
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
    
    def run() {

        def cli = new CliBuilder(usage:' ffxc Minimize [options] <filename> [file2...]', header:' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }
        
        try {
            // getAlgorithmUtils is a magic variable/closure passed in from ModelingShell
            aFuncts = getAlgorithmUtils();
        } catch (MissingMethodException ex) {
            // This is the fallback, which does everything necessary without magic names
            aFuncts = new AlgorithmUtils();
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
        
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = aFuncts.getActiveAssembly();
            if (mola == null) {
                return cli.usage();
            }
            arguments = new ArrayList<>();
            arguments.add(mola.getFile().getName());
            
            processFile(options, mola, 0);
        } else {
            for (int i = 0; i < nArgs; i++) {
                openFile(options, arguments.get(i), i);
            }
        }

        boolean lambdaTerm = false;
        if ((options.lambda >= 0.0 && options.lambda <= 1.0) || nArgs > 1) {
            lambdaTerm = true;
            System.setProperty("lambdaterm", 'true');
            // Check that options.lambda is within the limit 0..1
            if (options.lambda < 0.0 || options.lambda > 1.0) {
                options.lambda = 0.0;
            }
        }
        
        Potential potential;
        
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
        
        List<Integer> uniqueA;
        List<Integer> uniqueB;
        if (nArgs >= 4) {
            uniqueA = new ArrayList<>();
            uniqueB = new ArrayList<>();
            
            if (options.unsharedA) {
                //rangesA = options.uaA.tokenize(".");
                def ra = [] as Set;
                String[] toksA = options.unsharedA.split(".");
                for (range in toksA) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        logger.info(String.format("Range %s A rangeStart %d groupCount %d", range, rangeStart, m.groupCount()));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            ra.add(i-1);
                        }
                    }
                }
                Atom[] atA1 = topologies[0].getAtomArray();
                int counter = 0;
                def raAdj = [] as Set; // Indexed by common variables in dtA.
                for (int i = 0; i < atA1.length; i++) {
                    Atom ai = atA1[i];
                    if (i in ra) {
                        if (ai.applyLambda()) {
                            logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                        } else {
                            logger.info(String.format(" Unshared A: %d variables %d-%d", i, counter, counter+2));
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
                String[] toksB = options.unsharedB.split(".");
                for (range in toksB) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        for (int i = rangeStart; i <= rangeEnd; i++) {
                            rb.add(i-1);
                        }
                    }
                }
                Atom[] atB1 = topologies[2].getAtomArray();
                int counter = 0;
                def rbAdj = [] as Set; // Indexed by common variables in dtA.
                for (int i = 0; i < atB1.length; i++) {
                    Atom bi = atB1[i];
                    if (i in rb) {
                        if (bi.applyLambda()) {
                            logger.warning(String.format(" Ranges defined in uaA should not overlap with ligand atoms; they are assumed to not be shared."));
                        } else {
                            logger.info(String.format(" Unshared B: %d variables %d-%d", i, counter, counter+2));
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
        
        switch (nArgs) {
            case 1:
                logger.info("\n Running minimize on " + topologies[0].getFile().getName());
                potential = energies[0];
                break;
            case 2:
                logger.info(String.format("\n Running minimize on dual topology [%s,%s]", topologies[0].getFile().getName(), topologies[1].getFile().getName()));
                DualTopologyEnergy dte = new DualTopologyEnergy(topologies[0], topologies[1], sf);
                if (numParallel == 2) {
                    dte.setParallel(true);
                }
                potential = dte;
                break;
            case 4:
                StringBuilder sb = new StringBuilder("\n Running minimize on quad topology ");
                sb.append(topologies.stream().map{t -> t.getFile().getName()}.collect(Collectors.joining(",", "[", "]")));
                logger.info(sb.toString());
                
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
                break;
            case 8:
                logger.warning("Oct-topology minimization not presently functional");
                sb = new StringBuilder("\n Running minimize on oct-topology ");
                sb.append(topologies.stream().map{t -> t.getFile().getName()}.collect(Collectors.joining(",", "[", "]")));
                logger.info(sb.toString());
                
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
                break;
            default:
                logger.severe(" Must have 1, 2, 4, or 8 topologies!");
                break;
        }
        logger.info(" RMS gradient convergence criteria: " + options.eps);
        
        LambdaInterface linter = (LambdaInterface) potential;
        
        double[] x = new double[potential.getNumberOfVariables()];
        potential.getCoordinates(x);
        potential.energy(x, true);
        
        AlgorithmListener theListener = aFuncts.getDefaultListener();
        Minimize minimize = new Minimize(topologies[0], potential, theListener);
        minimize.minimize(options.eps);
        
        // This does not work, and seems to have poorly-scaled output
        //potential.getCoordinates(x);
        //potential.energy(x, true);
        
        // This does work... for singular topologies.
        //energy();
        
        for (mola in topologies) {
            String filename = mola.getFile().getName();
            String ext = FilenameUtils.getExtension(filename);
            filename = FilenameUtils.removeExtension(filename);
            
            if (ext.toUpperCase().contains("XYZ")) {
                aFuncts.saveAsXYZ(mola, new File(filename + ".xyz"));
            } else {
                aFuncts.saveAsPDB(mola, new File(filename + ".pdb"));
            }
        }
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
