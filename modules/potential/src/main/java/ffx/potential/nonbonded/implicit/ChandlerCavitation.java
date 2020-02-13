//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.nonbonded.implicit;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.pow;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;

/**
 * The ChandlerCavitation class smoothly switches between a volume based dependence for small solutes to a
 * surface area dependence for large solutes.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ChandlerCavitation {

    private static final Logger logger = Logger.getLogger(ChandlerCavitation.class.getName());

    /**
     * Surface area (Ang^2).
     */
    private double surfaceArea;
    /**
     * Surface area energy (kcal/mol).
     */
    private double surfaceAreaEnergy;
    /**
     * Volume (Ang^3).
     */
    private double volume;
    /**
     * Volume energy (kcal/mol).
     */
    private double volumeEnergy;
    /**
     * Cavitation energy, which is a function of volume and/or surface area.
     */
    private double cavitationEnergy;
    /**
     * Effective radius probe.
     * <p>
     * In cavitation volume scaling regime, approximate solvent excluded volume
     * and effective radius are computed as follow.
     * <p>
     * 1) GaussVol vdW volume is computed from defined radii.
     * 2) Effective radius is computed as Reff = cbrt(3.0 * volume / (4.0 * PI)) + effectiveRadiusProbe.
     * 3) Solvent Excluded Volume = 4/3 * Pi * Reff^3
     */
    private double effectiveRadius;
    /**
     * Solvent pressure in kcal/mol/Ang^3.
     * Original value from Schnieders thesis work: 0.0327
     * Value based on testing with Schnieders thesis test set, Sept 2019: 0.11337
     */
    private double solventPressure = GeneralizedKirkwood.DEFAULT_SOLVENT_PRESSURE;
    /**
     * Surface tension in kcal/mol/Ang^2.
     */
    private double surfaceTension = GeneralizedKirkwood.DEFAULT_CAVDISP_SURFACE_TENSION;
    /**
     * Radius where volume dependence crosses over to surface area dependence (approximately at 1 nm).
     * Originally 3.0*surfaceTension/solventPressure
     * Reset to 7.339 to match Tinker
     */
    private double crossOver = GeneralizedKirkwood.DEFAULT_CROSSOVER;
    private double switchRange = 3.5;
    private double saSwitchRangeOff = 3.9;
    /**
     * Begin turning off the Volume term.
     */
    private double volumeOff = crossOver - switchRange;
    /**
     * Volume term is zero at the cut-off.
     */
    private double volumeCut = crossOver + switchRange;
    /**
     * Begin turning off the SA term.
     */
    private double surfaceAreaOff = crossOver + saSwitchRangeOff;
    /**
     * SA term is zero at the cut-off.
     */
    private double surfaceAreaCut = crossOver - switchRange;
    /**
     * Volume multiplicative switch.
     */
    private MultiplicativeSwitch volumeSwitch = new MultiplicativeSwitch(volumeCut, volumeOff);
    /**
     * Surface area multiplicative switch.
     */
    private MultiplicativeSwitch surfaceAreaSwitch = new MultiplicativeSwitch(surfaceAreaCut, surfaceAreaOff);

    private int nAtoms;
    private Atom[] atoms;
    private final GaussVol gaussVol;
    private final ConnollyRegion connollyRegion;

    public ChandlerCavitation(Atom[] atoms, GaussVol gaussVol) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.gaussVol = gaussVol;
        this.connollyRegion = null;
    }

    public ChandlerCavitation(Atom[] atoms, ConnollyRegion connollyRegion) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.connollyRegion = connollyRegion;
        this.gaussVol = null;
    }

    public GaussVol getGaussVol() {
        return gaussVol;
    }

    /**
     * Compute molecular volume and surface area.
     *
     * @param positions Atomic positions to use.
     * @param gradient  Atomic coordinate gradient.
     * @return The cavitation energy.
     */
    public double energyAndGradient(double[][] positions, AtomicDoubleArray3D gradient) {
        if (gaussVol != null) {
            return energyAndGradientGausVol(positions, gradient);
        } else {
            return energyAndGradientConnolly(positions, gradient);
        }
    }

    /**
     * Compute the cavitation energy.
     */
    public double energyAndGradientConnolly(double[][] positions, AtomicDoubleArray3D gradient) {

        connollyRegion.init(atoms, true);
        connollyRegion.runVolume();

        // Calculate a purely surface area based cavitation energy.
        surfaceArea = connollyRegion.getSurfaceArea();
        surfaceAreaEnergy = surfaceArea * surfaceTension;

        // Calculate a purely volume based cavitation energy.
        volume = connollyRegion.getVolume();
        volumeEnergy = volume * solventPressure;

        // effectiveRadius = 0.5 * sqrt(surfaceArea / PI);
        // double reff = effectiveRadius;
        // double reff2 = reff * reff;
        // double reff3 = reff2 * reff;
        // double reff4 = reff3 * reff;
        // double reff5 = reff4 * reff;
        // double dReffdvdW = reff / (2.0 * surfaceArea);

        // Use Volume to find an effective cavity radius.
        effectiveRadius = cbrt(3.0 * volume / (4.0 * PI));
        double reff = effectiveRadius;
        double reff2 = reff * reff;
        double reff3 = reff2 * reff;
        double reff4 = reff3 * reff;
        double reff5 = reff4 * reff;
        double vdWVolPI23 = pow(volume / PI, 2.0 / 3.0);
        double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * vdWVolPI23);

        // Find the cavitation energy using a combination of volume and surface area dependence.
        if (reff < volumeOff) {
            // Find cavity energy from only the molecular volume.
            cavitationEnergy = volumeEnergy;
            double[][] volumeGradient = connollyRegion.getVolumeGradient();
            for (int i = 0; i < nAtoms; i++) {
                double dx = solventPressure * volumeGradient[0][i];
                double dy = solventPressure * volumeGradient[1][i];
                double dz = solventPressure * volumeGradient[2][i];
                gradient.add(0, i, dx, dy, dz);
            }
        } else if (reff <= volumeCut) {
            // Include a tapered molecular volume.
            double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
            cavitationEnergy = taper * volumeEnergy;
            double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            double factor = dtaper * volumeEnergy + taper * solventPressure;
            double[][] volumeGradient = connollyRegion.getVolumeGradient();
            for (int i = 0; i < nAtoms; i++) {
                double dx = factor * volumeGradient[0][i];
                double dy = factor * volumeGradient[1][i];
                double dz = factor * volumeGradient[2][i];
                gradient.add(0, i, dx, dy, dz);
            }

        }

        // TODO: We need to include the surface area contribution to the gradient.
        if (reff > surfaceAreaOff) {
            // Find cavity energy from only SA.
            cavitationEnergy = surfaceAreaEnergy;
            // addSurfaceAreaGradient(0.0, surfaceTension, gradient);
        } else if (reff >= surfaceAreaCut) {
            // Include a tapered surface area term.
            double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
            cavitationEnergy += taperSA * surfaceAreaEnergy;
            // double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            // addSurfaceAreaGradient(dtaperSA * surfaceAreaEnergy, taperSA * surfaceTension, gradient);
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("\n Volume:              %8.3f (Ang^3)", volume));
            logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
            logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
            logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
            logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
            logger.fine(format(" Effective Radius:    %8.3f (Ang)", reff));
        }

        return cavitationEnergy;
    }

    /**
     * Compute molecular volume and surface area.
     *
     * @param positions Atomic positions to use.
     * @param gradient  Atomic coordinate gradient.
     * @return The cavitation energy.
     */
    public double energyAndGradientGausVol(double[][] positions, AtomicDoubleArray3D gradient) {

        gaussVol.computeVolumeAndSA(positions);

        // Save the total molecular volume.
        volume = gaussVol.getVolume();
        double[] volumeGradient = gaussVol.getVolumeGradient();

        // Calculate a purely volume based cavitation energy.
        volumeEnergy = volume * solventPressure;

        // Calculate the surface area.
        surfaceArea = gaussVol.getSurfaceArea();
        double[] surfaceAreaGradient = gaussVol.getSurfaceAreaGradient();

        // Calculate a purely surface area based cavitation energy.
        surfaceAreaEnergy = surfaceArea * surfaceTension;

        // Use SA to find an effective cavity radius.
//        effectiveRadius = 0.5 * sqrt(surfaceArea / PI) + effectiveRadiusProbe;
//        double reff = effectiveRadius;
//        double reff2 = reff * reff;
//        double reff3 = reff2 * reff;
//        double reff4 = reff3 * reff;
//        double reff5 = reff4 * reff;
//        double dreff = reff / (2.0 * surfaceArea);

        // Use Volume to find an effective cavity radius.
        double vdwVolume = volume;

        effectiveRadius = cbrt(3.0 * vdwVolume / (4.0 * PI));
        double reff = effectiveRadius;
        double reff2 = reff * reff;
        double reff3 = reff2 * reff;
        double reff4 = reff3 * reff;
        double reff5 = reff4 * reff;
        double vdWVolPI23 = pow(vdwVolume / PI, 2.0 / 3.0);
        double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * vdWVolPI23);

        // double sevSASAChainRule = 4.0 * pow(2.0, 1.0/3.0) * reff / (pow(3.0, 2.0/3.0) * volPI23);
        // logger.info(format(" Volume Reff: %16.8f, %16.8f", reff, dreff));

        // Find the cavitation energy using a combination of volume and surface area dependence.
        if (reff < volumeOff) {
            // Find cavity energy from only the molecular volume.
            cavitationEnergy = volumeEnergy;
            addVolumeGradient(solventPressure, gradient, volumeGradient);
        } else if (reff <= volumeCut) {
            // Include a tapered molecular volume.
            double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
            double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            cavitationEnergy = taper * volumeEnergy;
            double factor = dtaper * volumeEnergy + taper * solventPressure;
            addVolumeGradient(factor, gradient, volumeGradient);
        }

        if (reff > surfaceAreaOff) {
            // Find cavity energy from only SA.
            cavitationEnergy = surfaceAreaEnergy;
            addSurfaceAreaGradient(0.0, surfaceTension, gradient, volumeGradient, surfaceAreaGradient);
        } else if (reff >= surfaceAreaCut) {
            // Include a tapered surface area term.
            double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
            double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            cavitationEnergy += taperSA * surfaceAreaEnergy;
            addSurfaceAreaGradient(dtaperSA * surfaceAreaEnergy, taperSA * surfaceTension, gradient,
                    volumeGradient, surfaceAreaGradient);
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("\n Volume:              %8.3f (Ang^3)", volume));
            logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
            logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
            logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
            logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
            logger.fine(format(" Effective Radius:    %8.3f (Ang)", reff));
        }

        return cavitationEnergy;
    }

    /**
     * Collect the volume base cavitation energy and gradient contributions.
     *
     * @param factor   dTaper/dReff * dReff/dVolvdW * Vsev + Tapered volume (A^3) * dVolSEV/dVolvdW
     * @param gradient Array to accumulate derivatives.
     */
    private void addVolumeGradient(double factor, AtomicDoubleArray3D gradient, double[] volumeGradient) {
        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            double gx = factor * volumeGradient[index++];
            double gy = factor * volumeGradient[index++];
            double gz = factor * volumeGradient[index];
            gradient.add(0, i, gx, gy, gz);
        }
    }

    /**
     * Collect the surface area based cavitation energy and gradient contributions.
     *
     * @param volFactor Factor to multiply vdW volume gradient by.
     * @param saFactor  Factor to multiply vdw surface area by.
     * @param gradient  Array to accumulate derivatives.
     */
    private void addSurfaceAreaGradient(double volFactor, double saFactor, AtomicDoubleArray3D gradient,
                                        double[] volumeGradient, double[] surfaceAreaGradient) {
        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            double gx = volFactor * volumeGradient[index] + saFactor * surfaceAreaGradient[index++];
            double gy = volFactor * volumeGradient[index] + saFactor * surfaceAreaGradient[index++];
            double gz = volFactor * volumeGradient[index] + saFactor * surfaceAreaGradient[index];
            gradient.add(0, i, gx, gy, gz);
        }
    }

    public double getEnergy() {
        return cavitationEnergy;
    }

    public double getSolventPressure() {
        return solventPressure;
    }

    public double getSurfaceTension() {
        return surfaceTension;
    }

    public void setSolventPressure(double solventPressure) {
        this.solventPressure = solventPressure;
    }

    public void setSurfaceTension(double surfaceTension) {
        this.surfaceTension = surfaceTension;
    }

    public void setCrossOver(double crossOver) {
        this.crossOver = crossOver;
        volumeOff = crossOver - switchRange;
        volumeCut = crossOver + switchRange;
        surfaceAreaOff = crossOver + saSwitchRangeOff;
        surfaceAreaCut = crossOver - switchRange;
        volumeSwitch = new MultiplicativeSwitch(volumeCut, volumeOff);
        surfaceAreaSwitch = new MultiplicativeSwitch(surfaceAreaCut, surfaceAreaOff);
    }

    public double getCrossOver() {
        return crossOver;
    }

    /**
     * Return Volume based cavitation energy.
     *
     * @return Volume based cavitation energy.
     */
    public double getVolumeEnergy() {
        return volumeEnergy;
    }

    /**
     * Return Surface Area based cavitation energy.
     *
     * @return Surface Area based cavitation energy.
     */
    public double getSurfaceAreaEnergy() {
        return surfaceAreaEnergy;
    }

    /**
     * Return Volume (A^3).
     *
     * @return Volume (A^3).
     */
    public double getVolume() {
        return volume;
    }

    /**
     * Return Surface Area (A^2).
     *
     * @return Surface Area (A^2).
     */
    public double getSurfaceArea() {
        return surfaceArea;
    }

    public double getEffectiveRadius() {
        return effectiveRadius;
    }

}
