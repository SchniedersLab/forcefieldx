
package ffx.potential;

// Java Imports
import java.util.regex.Pattern;
import java.util.stream.Collectors;

// Groovy Imports
import groovy.cli.Option;
import groovy.cli.Unparsed;
import groovy.util.CliBuilder;

// FFX Imports
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
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.nonbonded.MultiplicativeSwitch;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;

/**
 * The TestLambdaGradient script tests numeric gradients w.r.t. lambda against
 * numeric, finite-difference gradients
 * <br>
 * Usage:
 * <br>
 * ffxc TestLambdaGradient [options] &lt;filename [file2...]&gt;
 */
class TestLambdaGradient extends Script {

    /**
     * Options for the TestLambdaGradient Script.
     * <br>
     * Usage:
     * <br>
     * ffxc TestLambdaGradient [options] &lt;filename&gt;
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
         * -as or --activeStart starts an active set of atoms for single-topology lambda gradients.
         */
        @Option(shortName='as', longName='activeStart', defaultValue='1', description='Starting active atom (single-topology only).') int actStart;
        /**
         * -af or --activeFinal ends an active set of atoms for single-topology lambda gradients.
         */
        @Option(shortName='af', longName='activeFinal', defaultValue='-1', description='Starting active atom (single-topology only).') int actFinal;
        /**
         * -l or --lambda sets the lambda value to minimize at.
         */
        @Option(shortName='l', longName='lambda', defaultValue='-1', description='Lambda value') double initialLambda;
        /**
         * -np or --nParallel sets the number of topologies to evaluate in parallel; currently 1, 2, or 4.
         */
        @Option(shortName='np', longName='nParallel', defaultValue='1', description='Number of topologies to evaluate in parallel') int nPar;
        /**
         * -uaA or -unsharedA sets atoms unique to the A dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaA', longName='unsharedA', description='Unshared atoms in the A dual topology (period-separated hyphenated ranges)') String unsharedA;
        /**
         * -uaB or -unsharedB sets atoms unique to the B dual-topology, as period-separated hyphenated ranges or singletons.
         */
        @Option(shortName='uaB', longName='unsharedB', description='Unshared atoms in the B dual topology (period-separated hyphenated ranges)') String unsharedB;
        /**
         * -v or --verbose is a flag to print out energy at each step.
         */
        @Option(shortName='v', longName='verbose', defaultValue='false', description='Print out the energy for each step') boolean print;
        /**
         * -qi or --quasi-internal sets use of the quasi-internal multipole tensor formulation; not presently recommended for production use.
         */
        @Option(shortName='qi', longName='quasi-internal', defaultValue='false', description='Use quasi-internal multipole tensors.') boolean qi;
        /**
         * -d or --dx sets the finite-difference step size (in Angstroms).
         */
        @Option(shortName='d', longName='dx', defaultValue='1.0e-5', description='Finite-difference step size (in Angstroms)') double step;
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
         * -sdX or --skipdX disables the very long atomic gradient tests and 
         * only performs dU/dL, dU2/dL2, and dU2/dLdX tests. Useful if the 
         * script is being used to generate a dU/dL profile.
         */
        @Option(shortName='sdX', longName='skipdX', defaultValue='false', 
            description='Skip calculating per-atom dUdX values and only test lambda gradients') boolean skipAtomGradients;
        
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames;
    }
    
    
    // Following variables are largely intended to be shared by helper methods such as openFile.
    private static final Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");
    private int threadsAvail = edu.rit.pj.ParallelTeam.getDefaultThreadCount();
    private int threadsPer = threadsAvail;
    private PotentialsFunctions pFuncts;
    def ranges1 = []; // Groovy mechanism for creating an untyped ArrayList.
    def ranges2 = [];
    def rangesA = [];
    def rangesB = [];
    def topologies = []; // MolecularAssembly
    def properties = []; // CompositeConfiguration
    def energies = [];   // ForceFieldEnergy

    private void openFile(Options options, String toOpen, int topNum) {
        MolecularAssembly[] opened = pFuncts.openAll(toOpen, threadsPer);
        MolecularAssembly mola = pFuncts.getActiveAssembly();
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
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
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
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
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
     * Compares MSNodes based on:
     * Residue vs. Molecule (must match)
     * Residue name
     * Residue number
     * Chain ID
     */
    /*boolean fuzzyCompareNodes(MSNode nodeA, MSNode nodeB) {
        if (nodeA instanceof Residue) {
            if (nodeB instanceof Residue) {
                Residue resA = (Residue) nodeA;
                Residue resB = (Residue) nodeB;
                if (!(resA.getName().equals(resB.getName()))) {
                    return false;
                }
                if (resA.getResidueNumber() != resB.getResidueNumber()) {
                    return false;
                }
                if (!(resA.getChainID().equals(resB.getChainID()))) {
                    return false;
                }
                return true;
            } else {
                return false;
            }
        } else if (nodeA instanceof Molecule) {
            if (nodeB instanceof Molecule) {
                Molecule molecA = (Molecule) nodeA;
                Molecule molecB = (Molecule) nodeB;
                if (!(molecA.getName().equals(molecB.getName()))) {
                    return false;
                }
                if (molecA.getResidueNumber() != resB.getResidueNumber()) {
                    return false;
                }
                if (!(molecA.getChainID().equals(molecB.getChainID()))) {
                    return false;
                }
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }*/
    
    /**
     * Script run method.
     */
    def run() {

        def cli = new CliBuilder(usage:' ffxc TestLambdaGradient [options] <filename> [file2...]', header:' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }
        
        try {
            pFuncts = getPotentialsUtils();
        } catch (MissingMethodException ex) {
            pFuncts = new PotentialsUtils();
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

        if (options.qi) {
            System.setProperty("pme-qi","true");
            System.setProperty("no-ligand-condensed-scf","false");
            System.setProperty("ligand-vapor-elec","false");
            System.setProperty("polarization","NONE");
        }
        boolean debug = (System.getProperty("debug") != null);

        // Turn on computation of lambda derivatives
        System.setProperty("lambdaterm","true");

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Condensed phase polarization is evaluated over the entire range.
            System.setProperty("polarization-lambda-start","0.0");
            // Polarization energy is not scaled individually by lambda, but
            // along with the overall potential energy of a topology.
            System.setProperty("polarization-lambda-exponent","0.0");
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec","false");
            // Condensed phase polarization, without the ligand present, is unecessary.
            System.setProperty("no-ligand-condensed-scf","false");
        }
        
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = aFuncts.getActiveAssembly();
            if (mola == null) {
                return cli.usage();
            }
            arguments = new ArrayList<>();
            arguments.add(mola.getFile().getName());
            
            processFile(mola);
        } else {
            logger.info(String.format(" Initializing %d topologies...", nArgs));
            for (int i = 0; i < nArgs; i++) {
                openFile(options, arguments.get(i), i);
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
                String[] toksA = options.unsharedA.tokenize(".");
                for (range in toksA) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        logger.info(String.format("Range %s for A, start %d end %d", range, rangeStart, rangeEnd));
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
                for (range in toksB) {
                    def m = rangeregex.matcher(range);
                    if (m.find()) {
                        int rangeStart = Integer.parseInt(m.group(1));
                        int rangeEnd = (m.groupCount() > 1) ? Integer.parseInt(m.group(2)) : rangeStart;
                        if (rangeStart > rangeEnd) {
                            logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                        }
                        logger.info(String.format("Range %s for A, start %d end %d", range, rangeStart, rangeEnd));
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
            
            /*if (options.unsharedRadius > 0.0) {
                Atom[] atomsA = topologies[0].getAtomArray();
                Atom[] atomsB = topologies[2].getAtomArray();
                NeighborList nlistA = energies[0].getVdwNode().getNeighborList();
                NeighborList nlistB = energies[2].getVdwNode().getNeighborList();
                Set<Integer> nearbyA = nlistA.getNeighborIndices(uniqueA, options.unsharedRadius);
                Set<Integer> nearbyB = nlistB.getNeighborIndices(uniqueB, options.unsharedRadius);
                
                Set<MSNode> nodesA = nearbyA.stream().map({i -> return atomsA[i];}).
                    map({at -> return at.getMoleculeOrResidue();}).
                    filter({node -> return node != null;}).
                    collect(Collectors.toSet());
                
                Set<MSNode> nodesB = nearbyB.stream().map({i -> return atomsB[i];}).
                    map({at -> return at.getMoleculeOrResidue();}).
                    filter({node -> return node != null;}).
                    collect(Collectors.toSet());
            }*/
        }
        
        StringBuilder sb = new StringBuilder("\n Testing lambda derivatives for ");
        switch (nArgs) {
            case 1:
                potential = energies[0];
                if (options.actFinal > 0) {
                    // Apply active atom selection
                    int nAtoms1 = (energies[0].getNumberOfVariables()) / 3;
                    if (options.actFinal > options.actStart && options.actStart > 0 && options.actFinal <= nAtoms1) {
                        // Make all atoms inactive.
                        for (int i = 0; i <= nAtoms1; i++) {
                            Atom ai = atoms[i - 1];
                            ai.setActive(false);
                        }
                        // Make requested atoms active.
                        for (int i = options.actStart; i <= options.actFinal; i++) {
                            Atom ai = atoms[i - 1];
                            ai.setActive(true);
                        }
                    } 
               }
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
                break;
            default:
                logger.severe(" Must have 1, 2, 4, or 8 topologies!");
                break;
        }
        sb.append(topologies.stream().map{t -> t.getFile().getName()}.collect(Collectors.joining(",", "[", "]")));
        logger.info(sb.toString());
        
        LambdaInterface linter = (LambdaInterface) potential;
        
        // End boilerplate open-topologies code.

        // Reset the number of variables for the case of dual topology.
        int n = potential.getNumberOfVariables();
        double[] x = new double[n];
        double[] gradient = new double[n];
        double[] lambdaGrad = new double[n];
        double[][] lambdaGradFD = new double[2][n];

        // Number of independent atoms.
        assert(n % 3 == 0);
        int nAtoms = n / 3;

        // Compute the Lambda = 0.0 energy.
        double lambda = 0.0;
        linter.setLambda(lambda);
        potential.getCoordinates(x);
        double e0 = potential.energyAndGradient(x,gradient);

        // Compute the Lambda = 1.0 energy.
        lambda = 1.0;
        linter.setLambda(lambda);
        double e1 = potential.energyAndGradient(x,gradient);

        logger.info(String.format(" E(0):      %20.8f.", e0));
        logger.info(String.format(" E(1):      %20.8f.", e1));
        logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1-e0));

        // Finite-difference step size.
        double width = 2.0 * options.step;

        // Error tolerence
        double errTol = 1.0e-3;
        // Upper bound for typical gradient sizes (expected gradient)
        double expGrad = 1000.0;

        // Test Lambda gradient in the neighborhood of the lambda variable.
        for (int j=0; j<3; j++) {
            lambda = options.initialLambda - 0.01 + 0.01 * j;

            if (lambda - options.step < 0.0) {
                continue;
            }
            if (lambda + options.step > 1.0) {
                continue;
            }

            logger.info(String.format(" Current lambda value %6.4f", lambda));
            linter.setLambda(lambda);

            // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
            double e = potential.energyAndGradient(x,gradient);

            // Analytic dEdL, d2E/dL2 and dE/dL/dX
            double dEdL = linter.getdEdL();
            double d2EdL2 = linter.getd2EdL2();
            for (int i = 0; i < n; i++) {
                lambdaGrad[i] = 0.0;
            }
            potential.getdEdXdL(lambdaGrad);

            // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
            linter.setLambda(lambda + options.step);
            double lp = potential.energyAndGradient(x,lambdaGradFD[0]);
            double dedlp = linter.getdEdL();
            linter.setLambda(lambda - options.step);
            double lm = potential.energyAndGradient(x,lambdaGradFD[1]);
            double dedlm = linter.getdEdL();

            double dEdLFD = (lp - lm) / width;
            double d2EdL2FD = (dedlp - dedlm) / width;

            if (debug) {
                logger.info(String.format(" db> Lambda FD Test  lower,center,upper: %g %g %g", lambda - options.step , lambda, lambda + options.step));
                logger.info(String.format(" db> dE/dL   Numeric lp,lm,width,val: %+9.6g %+9.6g %4.2g (%+9.6g)", lp, lm, width, dEdLFD));
                logger.info(String.format(" db> d2E/dL2 Numeric lp,lm,width,val: %+9.6g %+9.6g %4.2g [%+9.6g]", dedlp, dedlm, width, d2EdL2FD));
                logger.info(String.format(" db> Analytic vers   l,dEdL,d2EdL2: %4.2f (%+9.6g) [%+9.6g]", lambda, dEdL, d2EdL2));
            }

            double err = Math.abs(dEdLFD - dEdL);
            if (err < errTol) {
                logger.info(String.format(" dE/dL passed:   %10.6f", err));
            } else {
                logger.info(String.format(" dE/dL failed: %10.6f", err));
            }
            logger.info(String.format(" Numeric:   %15.8f", dEdLFD));
            logger.info(String.format(" Analytic:  %15.8f", dEdL));

            err = Math.abs(d2EdL2FD - d2EdL2);
            if (err < errTol) {
                logger.info(String.format(" d2E/dL2 passed: %10.6f", err));
            } else {
                logger.info(String.format(" d2E/dL2 failed: %10.6f", err));
            }
            logger.info(String.format(" Numeric:   %15.8f", d2EdL2FD));
            logger.info(String.format(" Analytic:  %15.8f", d2EdL2));

            //boolean passed = true;
            int ndEdXdLFails = 0;

            double rmsError = 0;
            for (int i = 0; i < nAtoms; i++) {
                int ii = i * 3;
                double dX = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dXa = lambdaGrad[ii];
                double eX = dX - dXa;
                ii++;
                double dY = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dYa = lambdaGrad[ii];
                double eY = dY - dYa;
                ii++;
                double dZ = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width;
                double dZa = lambdaGrad[ii];
                double eZ = dZ - dZa;

                //double error = Math.sqrt(eX * eX + eY * eY + eZ * eZ);
                double error = eX * eX + eY * eY + eZ * eZ;
                rmsError += error;
                error = Math.sqrt(error);
                if (error < errTol) {
                    logger.fine(String.format(" dE/dX/dL for Atom %d passed: %10.6f", i + 1, error));
                } else {
                    logger.info(String.format(" dE/dX/dL for Atom %d failed: %10.6f", i + 1, error));
                    logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa,dYa,dZa));
                    logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ));
                    ndEdXdLFails++;
                }
            }
            rmsError = Math.sqrt(rmsError / nAtoms);
            if (ndEdXdLFails == 0) {
                logger.info(String.format(" dE/dX/dL passed for all atoms: RMS error %15.8f", rmsError));
            } else {
                logger.info(String.format(" dE/dX/dL failed for %d of %d atoms: RMS error %15.8f", ndEdXdLFails, nAtoms, rmsError));
            }

            logger.info("");
        }

        linter.setLambda(options.initialLambda);
        potential.getCoordinates(x);
        potential.energyAndGradient(x,gradient, options.print);

        if (!options.skipAtomGradients) {
            logger.info(String.format(" Checking Cartesian coordinate gradient"));
            double[] numeric = new double[3];
            double avLen = 0.0;
            int nFailures = 0;
            double avGrad = 0.0;
            for (int i=0; i<nAtoms; i++) {
                int i3 = i*3;
                int i0 = i3 + 0;
                int i1 = i3 + 1;
                int i2 = i3 + 2;

                // Find numeric dX
                double orig = x[i0];
                x[i0] = x[i0] + options.step;
                double e = potential.energyAndGradient(x,lambdaGradFD[0], options.print);
                x[i0] = orig - options.step;
                e -= potential.energyAndGradient(x,lambdaGradFD[1], options.print);
                x[i0] = orig;
                numeric[0] = e / width;

                // Find numeric dY
                orig = x[i1];
                x[i1] = x[i1] + options.step;
                e = potential.energyAndGradient(x,lambdaGradFD[0], options.print);
                x[i1] = orig - options.step;
                e -= potential.energyAndGradient(x,lambdaGradFD[1], options.print);
                x[i1] = orig;
                numeric[1] = e / width;

                // Find numeric dZ
                orig = x[i2];
                x[i2] = x[i2] + options.step;
                e = potential.energyAndGradient(x,lambdaGradFD[0], options.print);
                x[i2] = orig - options.step;
                e -= potential.energyAndGradient(x,lambdaGradFD[1], options.print);
                x[i2] = orig;
                numeric[2] = e / width;

                double dx = gradient[i0] - numeric[0];
                double dy = gradient[i1] - numeric[1];
                double dz = gradient[i2] - numeric[2];
                double len = dx * dx + dy * dy + dz * dz;
                avLen += len;
                len = Math.sqrt(len);

                double grad2 = gradient[i0] * gradient[i0] + gradient[i1] * gradient[i1] + gradient[i2] * gradient[i2];
                avGrad += grad2;

                if (len > errTol) {
                    logger.info(String.format(" Atom %d failed: %10.6f.",i+1,len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
                        + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
                    ++nFailures;
                    //return;
                } else {
                    logger.info(String.format(" Atom %d passed: %10.6f.",i+1,len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1], gradient[i2])
                        + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]));
                }

                if (grad2 > expGrad) {
                    logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i+1, grad2));
                }
                logger.info("\n");
            }

            avLen = avLen / nAtoms;
            avLen = Math.sqrt(avLen);
            if (avLen > errTol) {
                logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol));
            } else {
                logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol));
            }
            logger.info(String.format(" Number of atoms failing gradient test: %d", nFailures));

            avGrad = avGrad / nAtoms;
            avGrad = Math.sqrt(avGrad);
            if (avGrad > expGrad) {
                logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad));
            } else {
                logger.info(String.format(" RMS gradient: %10.6f", avGrad));
            }
        } else {
            logger.info(" Skipping atomic dU/dX gradients.");
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
