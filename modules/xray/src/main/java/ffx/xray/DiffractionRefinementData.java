//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.xray;

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.ReflectionList;
import ffx.numerics.math.ComplexNumber;

/**
 * <p>
 * DiffractionRefinementData class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class DiffractionRefinementData {

    private static final Logger logger = Logger.getLogger(DiffractionRefinementData.class.getName());
    /**
     * Number of reflections in the data set.
     */
    public final int n;
    /**
     * Number of scale parameters.
     */
    final int nScale;
    /**
     * 2D array of F/sigF data.
     */
    final double[][] fSigF;
    /**
     * Array of R free flags;
     */
    public final int[] freeR;
    /**
     * Calculated atomic structure factors.
     */
    public final double[][] fc;
    /**
     * Calculated bulk solvent structure factors.
     */
    public final double[][] fs;
    /**
     * Scaled sum of Fc and Fs
     */
    final double[][] fcTot;
    /**
     * Figure of merit and phase.
     */
    public final double[][] fomPhi;
    /**
     * 2mFo - DFc coefficients.
     */
    final double[][] foFc2;
    /**
     * mFo - DFc coefficients.
     */
    public final double[][] foFc1;
    /**
     * Derivatives with respect to Fc.
     */
    final double[][] dFc;
    /**
     * Derivatives with respect to Fs.
     */
    final double[][] dFs;
    /**
     * Log Likelihoods.
     */
    double llkR, llkF;
    /**
     * Reciprocal space reference for structure factor calculations and computing derivatives.
     */
    CrystalReciprocalSpace crystalReciprocalSpaceFc;
    /**
     * Reciprocal space reference for bulk solvent structure factor calculations and computing derivatives.
     */
    CrystalReciprocalSpace crystalReciprocalSpaceFs;
    /**
     * Number of resolution bins.
     */
    final int nBins;
    /**
     * Spine scaling coefficients.
     */
    public double[] spline;
    /**
     * SigmaA coefficient - s.
     */
    public double[] sigmaA;
    /**
     * SigmaA coefficient - w.
     */
    public double[] sigmaW;
    /**
     * Scaled, E-like Fc.
     */
    double[] esqFc;
    /**
     * Scaled, E-like Fo.
     */
    double[] esqFo;
    /**
     * Bulk solvent parameters.
     */
    double solventA, solventB;
    /**
     * bulk solvent scale and Ueq
     */
    double bulkSolventK, bulkSolventUeq;
    /**
     * Overall model scale.
     */
    double modelScaleK;
    /**
     * model anisotropic B
     */
    double[] modelAnisoB = new double[6];
    /**
     * Duplicated settings - these are also in DiffractionData, but duplicated
     * here until settings are put in their own class.
     */
    public int rFreeFlag;
    public final double fSigFCutoff;

    /**
     * allocate data given a {@link ffx.crystal.ReflectionList}
     *
     * @param properties     configuration properties
     * @param reflectionList {@link ffx.crystal.ReflectionList} to use to allocate data
     */
    public DiffractionRefinementData(CompositeConfiguration properties,
                                     ReflectionList reflectionList) {

        int rflag = -1;
        if (properties != null) {
            rflag = properties.getInt("rfreeflag", -1);
            fSigFCutoff = properties.getDouble("fsigfcutoff", -1.0);
        } else {
            fSigFCutoff = -1.0;
        }

        this.n = reflectionList.hkllist.size();
        this.nScale = reflectionList.crystal.scaleN;
        this.rFreeFlag = rflag;
        this.nBins = reflectionList.nbins;
        fSigF = new double[n][2];
        freeR = new int[n];
        fc = new double[n][2];
        fs = new double[n][2];
        fcTot = new double[n][2];
        fomPhi = new double[n][2];
        foFc2 = new double[n][2];
        foFc1 = new double[n][2];
        dFc = new double[n][2];
        dFs = new double[n][2];

        for (int i = 0; i < n; i++) {
            fSigF[i][0] = fSigF[i][1] = Double.NaN;
            fcTot[i][0] = fcTot[i][1] = Double.NaN;
        }

        spline = new double[nBins * 2];
        sigmaA = new double[nBins];
        sigmaW = new double[nBins];
        esqFc = new double[nBins];
        esqFo = new double[nBins];
        for (int i = 0; i < nBins; i++) {
            spline[i] = spline[i + nBins] = sigmaA[i] = esqFc[i] = esqFo[i] = 1.0;
            sigmaW[i] = 0.05;
        }

        // Initial guess for scaling/bulk solvent correction.
        bulkSolventK = 0.33;
        bulkSolventUeq = 50.0 / (8.0 * PI * PI);
        modelScaleK = 0.0;
    }

    /**
     * <p>
     * setCrystalReciprocalSpace_fc</p>
     *
     * @param crystalReciprocalSpace a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    void setCrystalReciprocalSpaceFc(CrystalReciprocalSpace crystalReciprocalSpace) {
        this.crystalReciprocalSpaceFc = crystalReciprocalSpace;
    }

    /**
     * <p>
     * setCrystalReciprocalSpace_fs</p>
     *
     * @param crystalReciprocalSpace a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    void setCrystalReciprocalSpaceFs(CrystalReciprocalSpace crystalReciprocalSpace) {
        this.crystalReciprocalSpaceFs = crystalReciprocalSpace;
    }

    /**
     * Mark 5% of reflections for cross validation (R free flags).
     */
    public void generateRFree() {
        Random generator = new Random();
        int free;
        int nonfree;
        int nfree = 0;

        if (rFreeFlag == 0) {
            free = 0;
            nonfree = 1;
        } else {
            free = 1;
            nonfree = 0;
        }

        for (int i = 0; i < n; i++) {
            if (isNaN(fSigF[i][0])) {
                freeR[i] = nonfree;
                continue;
            }

            int randomi = generator.nextInt(100);
            if (randomi < 5) {
                freeR[i] = free;
                nfree++;
            } else {
                freeR[i] = nonfree;
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            logger.info(format(" Assigned %d of %d reflections to the R free set (5%%).\n", nfree, n));
        }
    }

    /**
     * Generate average F/sigF from anomalous F/sigF.
     *
     * @param anomalousFsigF an array of double.
     */
    public void generateFsigFfromAnomalousFsigF(double[][] anomalousFsigF) {
        for (int i = 0; i < n; i++) {
            if (isNaN(anomalousFsigF[i][0]) && isNaN(anomalousFsigF[i][2])) {
                // Do Nothing.
            } else if (isNaN(anomalousFsigF[i][0])) {
                fSigF[i][0] = anomalousFsigF[i][2];
                fSigF[i][1] = anomalousFsigF[i][3];
            } else if (isNaN(anomalousFsigF[i][2])) {
                fSigF[i][0] = anomalousFsigF[i][0];
                fSigF[i][1] = anomalousFsigF[i][1];
            } else {
                fSigF[i][0] = (anomalousFsigF[i][0] + anomalousFsigF[i][2]) / 2.0;
                fSigF[i][1] = sqrt(anomalousFsigF[i][1] * anomalousFsigF[i][1]
                        + anomalousFsigF[i][3] * anomalousFsigF[i][3]);
            }
        }
    }

    /**
     * Generate amplitudes from intensities.
     * <p>
     * This does NOT use French and Wilson caling, but just a simple square root.
     */
    public void intensitiesToAmplitudes() {
        double tmp;

        for (int i = 0; i < n; i++) {
            if (fSigF[i][0] > 0.0) {
                tmp = fSigF[i][0];
                fSigF[i][0] = sqrt(tmp);
                if (fSigF[i][1] < tmp) {
                    fSigF[i][1] = fSigF[i][0]
                            - sqrt(tmp - fSigF[i][1]);
                } else {
                    fSigF[i][1] = fSigF[i][0];
                }
            } else if (!isNaN(fSigF[i][0])) {
                fSigF[i][0] = 0.0;
                fSigF[i][1] = 0.0;
            }
        }

        if (logger.isLoggable(Level.WARNING)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\n Internally converting intensities to amplitudes.\n");
            sb.append(" This does not use French & Wilson scaling.\n");
            logger.info(sb.toString());
        }
    }

    /**
     * Set amplitude (F).
     *
     * @param i reflection to set
     * @param f value of F desired
     */
    public void setF(int i, double f) {
        fSigF[i][0] = f;
    }

    /**
     * <p>
     * getF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getF(int i) {
        return fSigF[i][0];
    }

    /**
     * Set amplitude sigma (sigF).
     *
     * @param i    reflection to set
     * @param sigF value of sigF desired
     */
    public void setSigF(int i, double sigF) {
        fSigF[i][1] = sigF;
    }

    /**
     * <p>
     * getSigF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getSigF(int i) {
        return fSigF[i][1];
    }

    /**
     * Set amplitude and sigF.
     *
     * @param i    reflection to set
     * @param f    value of F desired
     * @param sigf value of sigF desired
     */
    public void setFSigF(int i, double f, double sigf) {
        fSigF[i][0] = f;
        fSigF[i][1] = sigf;
    }

    /**
     * <p>
     * getFSigF</p>
     *
     * @param i a int.
     * @return an array of double.
     */
    public double[] getFSigF(int i) {
        return fSigF[i];
    }

    /**
     * Set FreeR value flag.
     *
     * @param i If FreeR value is i, it is marked for cross validation.
     */
    public void setFreeRFlag(int i) {
        rFreeFlag = i;
    }

    /**
     * Set FreeR value flag of a reflection.
     *
     * @param i reflection to set
     * @param f FreeR value to set reflection to
     */
    public void setFreeR(int i, int f) {
        freeR[i] = f;
    }

    /**
     * <p>
     * getFreeR</p>
     *
     * @param i a int.
     * @return a int.
     */
    public int getFreeR(int i) {
        return freeR[i];
    }

    /**
     * <p>
     * isFreeR</p>
     *
     * @param i a int.
     * @param f a int.
     * @return a boolean.
     */
    public boolean isFreeR(int i, int f) {
        return (freeR[i] == f);
    }

    /**
     * <p>
     * isFreeR</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    boolean isFreeR(int i) {
        return (freeR[i] == rFreeFlag);
    }

    /**
     * set complex Fc
     *
     * @param i reflection to set
     * @param c {@link ComplexNumber} to set reflection to
     */
    public void setFc(int i, ComplexNumber c) {
        fc[i][0] = c.re();
        fc[i][1] = c.im();
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @return newly allocated {@link ComplexNumber}
     */
    ComplexNumber getFc(int i) {
        return new ComplexNumber(fc[i][0], fc[i][1]);
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @param c {@link ComplexNumber} to fill
     */
    void getFcIP(int i, ComplexNumber c) {
        c.re(fc[i][0]);
        c.im(fc[i][1]);
    }

    /**
     * get the amplitude of a complex Fc
     *
     * @param i reflection to get
     * @return amplitude of Fc
     */
    public double fcF(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
        return c.abs();
    }

    /**
     * get the phase of a complex Fc
     *
     * @param i reflection to get
     * @return phase of Fc
     */
    public double fcPhi(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
        return c.phase();
    }

    /**
     * <p>
     * setFs</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void setFs(int i, ComplexNumber c) {
        fs[i][0] = c.re();
        fs[i][1] = c.im();
    }

    /**
     * <p>
     * getFs</p>
     *
     * @param i a int.
     * @return a {@link ComplexNumber} object.
     */
    public ComplexNumber getFs(int i) {
        return new ComplexNumber(fs[i][0], fs[i][1]);
    }

    /**
     * <p>
     * getFsIP</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    void getFsIP(int i, ComplexNumber c) {
        c.re(fs[i][0]);
        c.im(fs[i][1]);
    }

    /**
     * <p>
     * fsF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fsF(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);
        return c.abs();
    }

    /**
     * <p>
     * fsPhi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fsPhi(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);
        return c.phase();
    }

    /**
     * <p>
     * setFcTot</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void setFcTot(int i, ComplexNumber c) {
        fcTot[i][0] = c.re();
        fcTot[i][1] = c.im();
    }

    /**
     * <p>
     * getFcTot</p>
     *
     * @param i a int.
     * @return a {@link ComplexNumber} object.
     */
    public ComplexNumber getFcTot(int i) {
        return new ComplexNumber(fcTot[i][0], fcTot[i][1]);
    }

    /**
     * <p>
     * getFcTotIP</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    void getFcTotIP(int i, ComplexNumber c) {
        c.re(fcTot[i][0]);
        c.im(fcTot[i][1]);
    }

    /**
     * <p>
     * fcTotF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fcTotF(int i) {
        ComplexNumber c = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
        return c.abs();
    }

    /**
     * <p>
     * fcTotPhi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fcTotPhi(int i) {
        ComplexNumber c = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
        return c.phase();
    }

    /**
     * <p>
     * setFoFc2</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void setFoFc2(int i, ComplexNumber c) {
        foFc2[i][0] = c.re();
        foFc2[i][1] = c.im();
    }

    /**
     * <p>
     * getFoFc2</p>
     *
     * @param i a int.
     * @return a {@link ComplexNumber} object.
     */
    public ComplexNumber getFoFc2(int i) {
        return new ComplexNumber(foFc2[i][0], foFc2[i][1]);
    }

    /**
     * <p>
     * getFoFc2IP</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void getFoFc2IP(int i, ComplexNumber c) {
        c.re(foFc2[i][0]);
        c.im(foFc2[i][1]);
    }

    /**
     * <p>
     * FoFc2F</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double FoFc2F(int i) {
        ComplexNumber c = new ComplexNumber(foFc2[i][0], foFc2[i][1]);
        return c.abs();
    }

    /**
     * <p>
     * FoFc2Phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double FoFc2Phi(int i) {
        ComplexNumber c = new ComplexNumber(foFc2[i][0], foFc2[i][1]);
        return c.phase();
    }

    /**
     * <p>
     * setFoFc1</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void setFoFc1(int i, ComplexNumber c) {
        foFc1[i][0] = c.re();
        foFc1[i][1] = c.im();
    }

    /**
     * <p>
     * getFoFc1</p>
     *
     * @param i a int.
     * @return a {@link ComplexNumber} object.
     */
    public ComplexNumber getFoFc1(int i) {
        return new ComplexNumber(foFc1[i][0], foFc1[i][1]);
    }

    /**
     * <p>
     * getFoFc1IP</p>
     *
     * @param i a int.
     * @param c a {@link ComplexNumber} object.
     */
    public void getFoFc1IP(int i, ComplexNumber c) {
        c.re(foFc1[i][0]);
        c.im(foFc1[i][1]);
    }

    /**
     * <p>
     * foFc1F</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double foFc1F(int i) {
        ComplexNumber c = new ComplexNumber(foFc1[i][0], foFc1[i][1]);
        return c.abs();
    }

    /**
     * <p>
     * foFc1Phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double foFc1Phi(int i) {
        ComplexNumber c = new ComplexNumber(foFc1[i][0], foFc1[i][1]);
        return c.phase();
    }

    /**
     * return the current likelihood
     *
     * @return the work likelihood (non-Rfree based)
     */
    public double likelihoodWork() {
        return llkR;
    }

    /**
     * return the current likelihood
     *
     * @return the free likelihood (Rfree based)
     */
    public double likelihoodFree() {
        return llkF;
    }
}
