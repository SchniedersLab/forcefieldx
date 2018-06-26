/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms.cli;

import picocli.CommandLine.Option;

/**
 * Barostat options with Pico CLI format.
 * @author hbernabe
 */
public class BarostatOptions {

    /**
     * -p or --npt Specify use of a MC Barostat at the given pressure (default
     * 0.0 atm).
     */
    @Option(names = {"-p", "--npt"}, paramLabel = "0",
            description = "Specify use of a MC Barostat at the given pressure; the default 0 disables NPT (atm).")
    double pressure = 0;
    /**
     * --ld or --minDensity sets a tin box constraint on the barostat, preventing over-expansion of the box (particularly in vapor phase), permitting an analytic correction.
     */
    @Option(names = {"--ld", "--minDensity"}, paramLabel = "0.75",
            description = "Minimum density allowed by the barostat")
    double minDensity = 0.75;
    /**
     * --hd or --maxDensity sets a maximum density on the barostat, preventing under-expansion of the box.
     */
    @Option(names = {"--hd", "--maxDensity"}, paramLabel = "1.6",
            description = "Maximum density allowed by the barostat")
    double maxDensity = 1.6;
    /**
     * --sm or --maxSideMove sets the width of proposed crystal side length moves (rectangularly distributed) in Angstroms.
     */
    @Option(names = {"--sm", "--maxSideMove"}, paramLabel = "0.25",
            description = "Maximum side move allowed by the barostat in Angstroms")
    double maxSideMove = 0.25;
    /**
     * --am or --maxAngleMove sets the width of proposed crystal angle moves (rectangularly distributed) in degrees.
     */
    @Option(names = {"--am", "--maxAngleMove"}, paramLabel = "0.5",
            description = "Maximum angle move allowed by the barostat in degrees")
    double maxAngleMove = 0.5;
    /**
     * --mi or --meanInterval sets the mean number of MD steps (Poisson distribution) between barostat move proposals.
     */
    @Option(names = {"--mi", "--meanInterval"}, paramLabel = "10",
            description = "Mean number of MD steps between barostat move proposals.")
    int meanInterval = 10;
}
