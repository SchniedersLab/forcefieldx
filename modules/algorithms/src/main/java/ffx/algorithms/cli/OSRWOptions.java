package ffx.algorithms.cli;

import groovy.cli.Option;
import picocli.CommandLine;

import java.util.logging.Logger;

public class OSRWOptions {

    private static final Logger logger = Logger.getLogger(OSRWOptions.class.getName());

    /**
     * -a or --async sets asynchronous walker communication (recommended)
     */
    @CommandLine.Option(names = {"-a", "--asynchronous"},
            description = "Walker communication is asynchronous")
    private boolean async = false;

    /**
     * -o or --optimize saves low-energy snapshots discovered (only for single topology simulations).
     */
    @CommandLine.Option(names = {"-o", "--optimize"},
            description = "Optimize and save low-energy snapshots.")
    private boolean optimize = false;

    /**
     * -q or --equilibrate sets the number of equilibration steps prior to
     * production ttOSRW counts begin.
     */
    @CommandLine.Option(names = {"-q", "--equilibrate"}, paramLabel = "1000", description = "Number of equilibration steps before evaluation of thermodynamics.")
    private int nEquil = 1000;

    /**
     * -k or --checkpoint sets the restart save frequency in picoseconds (1.0 psec default).
     */
    @CommandLine.Option(names = {"-k", "--checkpoint"}, paramLabel = "1.0", description = "Interval to write out histogram and .dyn restart files.")
    private double checkpoint = 1.0;

    /**
     * -c or --count sets the number of time steps between OSRW counts.
     */
    @CommandLine.Option(names = {"-c", "--count"}, paramLabel = "10", description = "Time steps between MD-OSRW counts.")
    private int countFreq = 10;

    /**
     * -b or --bias sets the initial Gaussian bias magnitude in kcal/mol.
     */
    @CommandLine.Option(names = {"-g", "--biasMag"}, paramLabel = "0.05", description = "OSRW Gaussian bias magnitude (kcal/mol).")
    private double biasMag = 0.05;

    /**
     * --tp or --temperingParam sets the Dama et al tempering rate parameter,
     * in multiples of kBT.
     */
    @CommandLine.Option(names = {"--tp", "--temperingParam"}, paramLabel = "8.0", description = "Dama et al tempering rate parameter in multiples of kBT")
    private double temperParam = 8.0;

    /**
     * -rn or --resetNumSteps, ignores steps detected in .lam lambda-restart
     * files and thus resets the histogram; use -rn false to continue from
     * the end of any prior simulation.
     */
    @CommandLine.Option(names = {"--rn", "--resetNumSteps"},
            description = "Ignore prior steps logged in .lam files")
    private boolean resetNumSteps = false;

    /**
     * -dw or --distributeWalkers allows walkers to start from multiple
     * conformations; AUTO picks up per-walker conformations as
     * filename.pdb_(walker number), and specifying a residue starts a
     * rotamer optimization to generate side-chain configurations to start
     * from.
     */
    @CommandLine.Option(names = {"--dw", "--distributeWalkers"}, paramLabel = "OFF", description = "AUTO: Pick up per-walker configurations as [filename.pdb]_[num], or specify a residue to distribute on.")
    private String distributeWalkersString;
}
