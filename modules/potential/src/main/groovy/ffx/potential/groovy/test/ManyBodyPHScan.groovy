package ffx.potential.groovy.test

import edu.rit.pj.Comm
import ffx.potential.cli.PotentialScript
import ffx.utilities.FFXScript
import picocli.CommandLine

import static java.lang.String.format
import static java.lang.String.format

@CommandLine.Command(description = " Run a pH Scan with ManyBody", name = "ffxc test.ManyBodyPHScan")
class ManyBodyPHScan extends PotentialScript {

    /**
     * --spH --startpH Lower end of the pH range to be evaluated (exclusive).
     */
    @CommandLine.Option(names = ['--spH', "--startpH"], paramLabel = "0.0", defaultValue = "0.0",
            description = 'Lower end of the pH range to be evaluated.')
    double start = 0.0

    /**
     * --epH --endpH Upper end of the pH range to be evaluated (inclusive).
     */
    @CommandLine.Option(names = ['--epH', "--endpH"], paramLabel = "14.0", defaultValue = "14.0",
            description = 'Upper end of the pH range to be evaluated.')
    double end = 14.0

    /**
     * --ns --nSteps Number of steps for the given pH range.
     */
    @CommandLine.Option(names = ['--ns', "--nSteps"], paramLabel = "2.0", defaultValue = "2.0",
            description = 'Number of steps for a given pH range.')
    double nSteps = 2.0

    /**
     * The final argument(s) should be one or more filenames.
     */
    @CommandLine.Unmatched
    List<String> unmatched = null

    /**
     * ManyBodyPHScan Constructor.
     */
    ManyBodyPHScan() {
        this(new Binding())
    }

    /**
     * ManyBodyPHScan Constructor.
     * @param binding The Groovy Binding to use.
     */
    ManyBodyPHScan(Binding binding) {
        super(binding)
    }

    @Override
    ManyBodyPHScan run() {

        if (!init()) {
            return this
        }

        // Set a flag to avoid double use of MPI in downstream commands.
        System.setProperty("pj.use.mpi", "false")

        Class<? extends FFXScript> script = getScript(unmatched.get(0))
        Comm world = Comm.world()
        int numProc = world.size()
        int rank = world.rank()

        if (numProc > 1) {
            logger.info(format(" Number of processes:  %d", numProc))
            logger.info(format(" Rank of this process: %d", rank))
        }

        // Remove ManyBodyPHScan command.
        unmatched.remove(0)

        double[] pHValues = new double[nSteps]
        double stepSize = (end - start) / nSteps
        for(int i=1; i < nSteps+1; i++){
            pHValues[i-1] = start + (stepSize*i)
        }

        int pHIndex = unmatched.indexOf("0.0")
        for (int i = 0; i < nSteps; i++) {
            List<String> commandArgs = new ArrayList<>()
            for (String arg : unmatched) {
                unmatched.set(pHIndex, pHValues[i].toString())
                commandArgs.add(arg)
            }

            // Create a Binding for command line arguments.
            Binding binding = new Binding()
            binding.setVariable("args", commandArgs)

            Script groovyScript = script.getDeclaredConstructor().newInstance()
            groovyScript.setBinding(binding)

            try {
                groovyScript.run()
            } catch (Exception e) {
                logger.info(format(" Exception for pH value: %s", pHValues[i].toString()))
            }
        }
        // Clear the pj.use.mpi flag.
        System.clearProperty("pj.use.mpi")

        return this

    }
}
