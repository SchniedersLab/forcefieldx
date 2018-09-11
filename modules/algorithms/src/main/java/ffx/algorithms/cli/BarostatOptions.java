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
package ffx.algorithms.cli;

import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.Barostat;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;

import picocli.CommandLine.Option;

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
     * Sets a tin box constraint on the barostat, preventing over-expansion of the box (particularly in vapor phase), permitting an analytic correction.
     */
    double minDensity = 0.75;

    /**
     * Sets a maximum density on the barostat, preventing under-expansion of the box.
     */
    double maxDensity = 1.6;

    /**
     * Sets the width of proposed crystal side length moves (rectangularly distributed) in Angstroms.
     */
    double maxSideMove = 0.25;

    /**
     * Sets the width of proposed crystal angle moves (rectangularly distributed) in degrees.
     */
    double maxAngleMove = 0.5;

    /**
     * Sets the mean number of MD steps (Poisson distribution) between barostat move proposals.
     */
    int meanInterval = 10;

    /**
     * If pressure has been set &gt; 0, creates a Barostat around a CrystalPotential, else
     * returns the original, unmodified CrystalPotential.
     *
     * @param assembly         Primary assembly of the CrystalPotential.
     * @param crystalPotential A CrystalPotential.
     * @return Either a Barostat (NPT enabled) or cpot.
     */
    public CrystalPotential checkNPT(MolecularAssembly assembly, CrystalPotential crystalPotential) {
        if (pressure > 0) {
            return createBarostat(assembly, crystalPotential);
        } else {
            return crystalPotential;
        }
    }

    /**
     * Creates a Barostat around a CrystalPotential.
     *
     * @param assembly         Primary assembly of the CrystalPotential.
     * @param crystalPotential A CrystalPotential.
     * @return An NPT potential.
     * @throws IllegalArgumentException If this BarostatOptions has pressure &lt;= 0.
     */
    public Barostat createBarostat(MolecularAssembly assembly,
                                   CrystalPotential crystalPotential) throws IllegalArgumentException {
        if (pressure > 0) {

            CompositeConfiguration properties = assembly.getProperties();

            minDensity = properties.getDouble("minDensity", 0.75);
            maxDensity = properties.getDouble("maxDensity", 1.6);
            maxSideMove = properties.getDouble("maxSideMove", 0.25);
            maxAngleMove = properties.getDouble("maxAngleMove", 0.5);
            meanInterval = properties.getInt("meanInterval", 10);

            Barostat barostat = new Barostat(assembly, crystalPotential);
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
