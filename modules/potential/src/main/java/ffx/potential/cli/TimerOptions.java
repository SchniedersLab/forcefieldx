package ffx.potential.cli;

import picocli.CommandLine.Option;

public class TimerOptions {
    /**
     * -n or --iterations to set the number of iterations
     */
    @Option(names = {"-n", "--iterations"}, paramLabel = "5", description = "Number of iterations.")
    private int iterations = 5;
    /**
     * -c or --threads to set the number of SMP threads (the default of 0 specifies use of all CPU cores)
     */
    @Option(names = {"-c", "--threads"}, paramLabel = "0", description = "Number of SMP threads (0 specifies use of all CPU cores)")
    private int threads = 0;
    /**
     * -g or --gradient to ignore computation of the atomic coordinates gradient
     */
    @Option(names = {"-g", "--gradient"}, description = "Ignore computation of the atomic coordinates gradient")
    private boolean gradient = false;
    /**
     * -q or --quiet to suppress printing of the energy for each iteration
     */
    @Option(names = {"-q", "--quiet"}, description = "Suppress printing of the energy for each iteration")
    private boolean quiet = false;
}
