package ffx.algorithms.cli;

import ffx.algorithms.dynamics.MolecularDynamics;
import picocli.CommandLine;

import java.util.logging.Logger;

public class RepExOptions {

    private static final Logger logger = Logger.getLogger(RepExOptions.class.getName());


    /**
     * The ArgGroup keeps the RepExOptions together when printing help.
     */
    @CommandLine.ArgGroup(heading = "%n Replica Exchange Options%n", validate = false)
    public RepExOptions.RepExOptionGroup group = new RepExOptions.RepExOptionGroup();

    public boolean getRepEx() { return group.repEx; }

    public void setRepEx(boolean repEx) { this.group.repEx = repEx; }

    public int getReplicaSteps() { return group.replicaSteps; }

    public void setReplicaSteps(int numSteps) { this.group.replicaSteps = numSteps; }

    public double getExponent() { return group.exponent; }

    public void setExponent(double exponent) { this.group.exponent = exponent; }




    private static class RepExOptionGroup {

        /**
         * -r or --repEx to execute temperature replica exchange.
         */
        @CommandLine.Option(
                names = {"-x", "--repEx"},
                paramLabel = "false",
                description = "Execute temperature replica exchange")
        boolean repEx = false;

        /**
         * --rs or --replicaSteps number of steps for temperature replica exchange
         */
        @CommandLine.Option(
                names = {"--rs", "--replicaSteps"},
                paramLabel = "100",
                defaultValue = "100",
                description = "Number of steps for replica exchange.")
        private int replicaSteps = 100;

        /**
         * -e or --exponent exponent to set the exponential temperature ladder for temperature replica exchange
         */
        @CommandLine.Option(
                names = {"-e", "--exponent"},
                paramLabel = "0.05",
                defaultValue = "0.05",
                description = "Exponent to set the exponential temperature ladder.")
        private double exponent = 0.05;



    }
}

