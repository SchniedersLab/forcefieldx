/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.algorithms.cli;

import ffx.algorithms.Barostat;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import picocli.CommandLine.Option;

import java.util.logging.Logger;

/**
 * Represents command line options for scripts that use a barostat/NPT.
 *
 * @author Michael J. Schnieders
 * @author Hernan V. Bernabe
 * @since 1.0
 */
public class BarostatOptions {

    private static final Logger logger = Logger.getLogger(BarostatOptions.class.getName());

    /**
     * -p or --npt Specify use of a MC Barostat at the given pressure (default 0 = constant volume).
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

    /**
     * If pressure has been set > 0, creates a Barostat around a CrystalPotential, else
     * returns the original, unmodified CrystalPotential.
     *
     * @param assembly Primary assembly of the CrystalPotential.
     * @param cpot A CrystalPotential.
     * @return Either a Barostat (NPT enabled) or cpot.
     */
    public CrystalPotential checkNPT(MolecularAssembly assembly, CrystalPotential cpot) {
        if (pressure > 0) {
            return createBarostat(assembly, cpot);
        } else {
            return cpot;
        }
    }

    /**
     * Creates a Barostat around a CrystalPotential.
     *
     * @param assembly Primary assembly of the CrystalPotential.
     * @param cpot A CrystalPotential.
     * @return An NPT potential.
     * @thros IllegalArgumentException If this BarostatOptions has pressure <= 0.
     */
    public Barostat createBarostat(MolecularAssembly assembly, CrystalPotential cpot) throws IllegalArgumentException {
        if (pressure > 0) {
            Barostat barostat = new Barostat(assembly, cpot);
            barostat.setPressure(pressure);
            barostat.setMaxDensity(maxDensity);
            barostat.setMinDensity(minDensity);
            double dens = barostat.density();
            if (dens < minDensity) {
                logger.info(String.format(" Barostat: initial density %9.4g < minimum density %9.4g, resetting to minimum density", dens, minDensity));
                barostat.setDensity(minDensity);
            } else if (dens > maxDensity) {
                logger.info(String.format(" Barostat: initial density %9.4g > maximum density %9.4g, resetting to maximum density", dens, minDensity));
                barostat.setDensity(maxDensity);
            }
            barostat.setMaxSideMove(maxSideMove);
            barostat.setMaxAngleMove(maxAngleMove);
            barostat.setMeanBarostatInterval(meanInterval);
            return barostat;
        } else {
            throw new IllegalArgumentException(" Pressure is <= 0; cannot create a Barostat!");
        }
    }
}
