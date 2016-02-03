/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.xray;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.CCP4MapWriter;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.parsers.PDBFilter;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.MTZWriter.MTZType;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.xray.CrystalReciprocalSpace.SolventModel.POLYNOMIAL;

import static java.util.Arrays.fill;
import static java.util.Arrays.fill;
import static java.util.Arrays.fill;
import static java.util.Arrays.fill;
import static java.util.Arrays.fill;
import static java.util.Arrays.fill;
import static java.util.Arrays.fill;

/**
 * <p>
 * DiffractionData class.</p>
 *
 * @author Timothy D. Fenn
 *
 */
public class DiffractionData implements DataContainer {

    private static final Logger logger = Logger.getLogger(DiffractionData.class.getName());
    protected final MolecularAssembly assembly[];
    protected final String modelName;
    protected final DiffractionFile[] dataFiles;
    protected final int n;
    protected final Crystal crystal[];
    protected final Resolution resolution[];
    protected final ReflectionList[] reflectionList;
    protected final DiffractionRefinementData[] refinementData;
    protected final CrystalReciprocalSpace crs_fc[];
    protected final CrystalReciprocalSpace crs_fs[];
    public final SolventModel solventModel;
    protected final RefinementModel refinementModel;
    protected ScaleBulkMinimize[] scaleBulkMinimize;
    protected SigmaAMinimize[] sigmaAMinimize;
    protected SplineMinimize[] splineMinimize;
    protected CrystalStats[] crystalStats;
    protected ParallelTeam parallelTeam;
    protected boolean scaled[];
    // settings
    public int rFreeFlag;
    public final double fsigfCutoff;
    public final boolean use_3g;
    public final double aRadBuff;
    public final double xrayScaleTol;
    public final double sigmaATol;
    public double xWeight;
    public final double bsimweight;
    public final double bnonzeroweight;
    public final double bmass;
    public final boolean residuebfactor;
    public final int nresiduebfactor;
    public final boolean addanisou;
    public final boolean refinemolocc;
    public final double occmass;
    public final boolean lambdaTerm;
    /**
     * if true, grid search bulk solvent params
     */
    public boolean gridsearch;
    /**
     * fit a scaling spline between Fo and Fc?
     */
    public boolean splinefit;

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object, used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties,
                POLYNOMIAL, new DiffractionFile(assembly));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object, used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties,
                POLYNOMIAL, datafile);
    }

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object, used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired - see
     * {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, SolventModel solventmodel) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                new DiffractionFile(assembly));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object, used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired - see
     * {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, SolventModel solventmodel,
            DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                datafile);
    }

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object array
     * (typically containing alternate conformer assemblies), used as the atomic
     * model for comparison against the data
     * @param properties system properties file
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties) {
        this(assembly, properties, POLYNOMIAL,
                new DiffractionFile(assembly[0]));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object array
     * (typically containing alternate conformer assemblies), used as the atomic
     * model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(assembly, properties, POLYNOMIAL, datafile);
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular assembly} object array
     * (typically containing alternate conformer assemblies), used as the atomic
     * model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired - see
     * {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, SolventModel solventmodel,
            DiffractionFile... datafile) {

        this.assembly = assembly;
        this.solventModel = solventmodel;
        this.modelName = assembly[0].getFile().getName();
        this.dataFiles = datafile;
        this.n = datafile.length;

        int rflag = properties.getInt("rfreeflag", -1);
        fsigfCutoff = properties.getDouble("fsigfcutoff", -1.0);
        gridsearch = properties.getBoolean("gridsearch", false);
        splinefit = properties.getBoolean("splinefit", true);
        use_3g = properties.getBoolean("use_3g", true);
        aRadBuff = properties.getDouble("aradbuff", 0.75);
        double sampling = properties.getDouble("sampling", 0.6);
        xrayScaleTol = properties.getDouble("xrayscaletol", 1e-4);
        sigmaATol = properties.getDouble("sigmaatol", 0.05);
        xWeight = properties.getDouble("xweight", 1.0);
        bsimweight = properties.getDouble("bsimweight", 1.0);
        bnonzeroweight = properties.getDouble("bnonzeroweight", 1.0);
        bmass = properties.getDouble("bmass", 5.0);
        residuebfactor = properties.getBoolean("residuebfactor", false);
        nresiduebfactor = properties.getInt("nresiduebfactor", 1);
        addanisou = properties.getBoolean("addanisou", false);
        refinemolocc = properties.getBoolean("refinemolocc", false);
        occmass = properties.getDouble("occmass", 10.0);
        lambdaTerm = properties.getBoolean("lambdaterm", false);

        crystal = new Crystal[n];
        resolution = new Resolution[n];
        reflectionList = new ReflectionList[n];
        refinementData = new DiffractionRefinementData[n];
        scaleBulkMinimize = new ScaleBulkMinimize[n];
        sigmaAMinimize = new SigmaAMinimize[n];
        splineMinimize = new SplineMinimize[n];
        crystalStats = new CrystalStats[n];

        // read in Fo/sigFo/FreeR
        File tmp;
        Crystal crystalinit = Crystal.checkProperties(properties);
        Resolution resolutioninit = Resolution.checkProperties(properties);
        if (crystalinit == null || resolutioninit == null) {
            for (int i = 0; i < n; i++) {
                tmp = new File(datafile[i].filename);
                reflectionList[i] = datafile[i].diffractionfilter.getReflectionList(tmp, properties);

                if (reflectionList[i] == null) {
                    logger.info(" Crystal information from the PDB or property file will be used.");
                    crystalinit = assembly[i].getCrystal().getUnitCell();
                    double res = datafile[i].diffractionfilter.getResolution(tmp, crystalinit);
                    if (res < 0.0) {
                        logger.severe("MTZ/CIF/CNS file does not contain full crystal information!");
                    } else {
                        resolutioninit = new Resolution(res, sampling);
                        reflectionList[i] = new ReflectionList(crystalinit, resolutioninit, properties);
                    }
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                reflectionList[i] = new ReflectionList(crystalinit, resolutioninit,
                        properties);
            }
        }

        for (int i = 0; i < n; i++) {
            crystal[i] = reflectionList[i].crystal;
            resolution[i] = reflectionList[i].resolution;
            refinementData[i] = new DiffractionRefinementData(properties, reflectionList[i]);
            tmp = new File(datafile[i].filename);
            datafile[i].diffractionfilter.readFile(tmp, reflectionList[i],
                    refinementData[i], properties);
        }

        if (!crystal[0].equals(assembly[0].getCrystal().getUnitCell())) {
            logger.severe("PDB and reflection file crystal information do not match! (check CRYST1 record?)");
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Refinement Settings\n\n");
            sb.append("  Target Function\n");
            sb.append("  X-ray refinement weight (xweight): ").append(xWeight).append("\n");
            sb.append("  Use cctbx 3 Gaussians (use_3g): ").append(use_3g).append("\n");
            sb.append("  Atomic form factor radius buffer (aradbuff): ").append(aRadBuff).append("\n");
            sb.append("  Reciprocal space sampling rate (sampling): ").append(sampling).append("\n");
            sb.append("  Resolution dependent spline scale (splinefit): ").append(splinefit).append("\n");
            sb.append("  Solvent grid search (gridsearch): ").append(gridsearch).append("\n");
            sb.append("  X-ray scale fit tolerance (xrayscaletol): ").append(xrayScaleTol).append("\n");
            sb.append("  Sigma A fit tolerance (sigmaatol): ").append(sigmaATol).append("\n\n");
            sb.append("  Reflections\n");
            sb.append("  F/sigF cutoff (fsigfcutoff): ").append(fsigfCutoff).append("\n");
            sb.append("  R Free flag (rfreeflag) (if -1, value will be updated when data is read in): ").
                    append(rflag).append("\n");
            sb.append("  Number of bins (nbins): ").append(reflectionList[0].nbins).append("\n\n");
            sb.append("  B-Factors\n");
            sb.append("  Similarity weight (bsimweight): ").append(bsimweight).append("\n");
            sb.append("  Non-zero weight (bnonzeroweight): ").append(bnonzeroweight).append("\n");
            sb.append("  Lagrangian mass (bmass): ").append(bmass).append("\n");
            sb.append("  Refined by residue (residuebfactor): ").append(residuebfactor).append("\n");
            sb.append("    (if true, num. residues per B (nresiduebfactor): ").append(nresiduebfactor).append(")\n");
            sb.append("  Add ANISOU for refinement (addanisou): ").append(addanisou).append("\n\n");
            sb.append("  Occupancies\n");
            sb.append("  Refine on molecules (HETATMs - refinemolocc): ").append(refinemolocc).append("\n");
            sb.append("  Lagrangian mass (occmass): ").append(occmass).append("\n");

            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementModel = new RefinementModel(assembly, refinemolocc);

        // initialize atomic form factors
        for (Atom a : refinementModel.usedAtoms) {
            a.setFormFactorIndex(-1);
            XRayFormFactor atomff
                    = new XRayFormFactor(a, use_3g, 2.0);
            a.setFormFactorIndex(atomff.ffIndex);

            if (a.getOccupancy() == 0.0) {
                a.setFormFactorWidth(1.0);
                continue;
            }

            double arad = a.getVDWType().radius * 0.5;
            double xyz[] = new double[3];
            xyz[0] = a.getX() + arad;
            xyz[1] = a.getY();
            xyz[2] = a.getZ();
            double rho = atomff.rho(0.0, 1.0, xyz);
            while (rho > 0.001) {
                arad += 0.1;
                xyz[0] = a.getX() + arad;
                rho = atomff.rho(0.0, 1.0, xyz);
            }
            arad += aRadBuff;
            a.setFormFactorWidth(arad);
        }

        // set up FFT and run it
        crs_fc = new CrystalReciprocalSpace[n];
        crs_fs = new CrystalReciprocalSpace[n];
        parallelTeam = new ParallelTeam();
        for (int i = 0; i < n; i++) {
            crs_fc[i] = new CrystalReciprocalSpace(reflectionList[i],
                    refinementModel.usedAtoms, parallelTeam, parallelTeam,
                    false, dataFiles[i].neutron);
            refinementData[i].setCrystalReciprocalSpace_fc(crs_fc[i]);
            crs_fc[i].setUse3G(use_3g);
            crs_fc[i].setWeight(dataFiles[i].weight);
            crs_fc[i].lambdaTerm = lambdaTerm;
            crs_fs[i] = new CrystalReciprocalSpace(reflectionList[i],
                    refinementModel.usedAtoms, parallelTeam, parallelTeam,
                    true, dataFiles[i].neutron, solventmodel);
            refinementData[i].setCrystalReciprocalSpace_fs(crs_fs[i]);
            crs_fs[i].setUse3G(use_3g);
            crs_fs[i].setWeight(dataFiles[i].weight);
            crs_fs[i].lambdaTerm = lambdaTerm;

            crystalStats[i] = new CrystalStats(reflectionList[i],
                    refinementData[i]);
        }

        scaled = new boolean[n];
        fill(scaled, false);
    }

    /**
     * read in a different assembly to average in structure factors
     *
     * @param assembly the
     * {@link ffx.potential.MolecularAssembly molecular assembly} object array
     * (typically containing alternate conformer assemblies), used as the atomic
     * model to average in with previous assembly
     * @param index the current data index (for cumulative average purposes)
     */
    public void AverageFc(MolecularAssembly assembly[], int index) {
        RefinementModel tmprefinementmodel = new RefinementModel(assembly, refinemolocc);

        // initialize atomic form factors
        for (Atom a : tmprefinementmodel.usedAtoms) {
            a.setFormFactorIndex(-1);
            XRayFormFactor atomff
                    = new XRayFormFactor(a, use_3g, 2.0);
            a.setFormFactorIndex(atomff.ffIndex);

            if (a.getOccupancy() == 0.0) {
                a.setFormFactorWidth(1.0);
                continue;
            }

            double arad = a.getVDWType().radius * 0.5;
            double xyz[] = new double[3];
            xyz[0] = a.getX() + arad;
            xyz[1] = a.getY();
            xyz[2] = a.getZ();
            double rho = atomff.rho(0.0, 1.0, xyz);
            while (rho > 0.001) {
                arad += 0.1;
                xyz[0] = a.getX() + arad;
                rho = atomff.rho(0.0, 1.0, xyz);
            }
            arad += aRadBuff;
            a.setFormFactorWidth(arad);
        }

        // set up FFT and run it
        for (int i = 0; i < n; i++) {
            crs_fc[i] = new CrystalReciprocalSpace(reflectionList[i],
                    tmprefinementmodel.usedAtoms, parallelTeam, parallelTeam,
                    false, dataFiles[i].neutron);
            refinementData[i].setCrystalReciprocalSpace_fc(crs_fc[i]);
            crs_fs[i] = new CrystalReciprocalSpace(reflectionList[i],
                    tmprefinementmodel.usedAtoms, parallelTeam, parallelTeam,
                    true, dataFiles[i].neutron, solventModel);
            refinementData[i].setCrystalReciprocalSpace_fs(crs_fs[i]);
        }

        int nhkl = refinementData[0].n;
        double fc[][] = new double[nhkl][2];
        double fs[][] = new double[nhkl][2];

        // run FFTs
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeDensity(fc);
            if (solventModel != SolventModel.NONE) {
                crs_fs[i].computeDensity(fs);
            }

            // average in with current data
            for (int j = 0; j < refinementData[i].n; j++) {
                refinementData[i].fc[j][0] += (fc[j][0] - refinementData[i].fc[j][0]) / index;
                refinementData[i].fc[j][1] += (fc[j][1] - refinementData[i].fc[j][1]) / index;

                refinementData[i].fs[j][0] += (fs[j][0] - refinementData[i].fs[j][0]) / index;
                refinementData[i].fs[j][1] += (fs[j][1] - refinementData[i].fs[j][1]) / index;
            }
            //refinementdata[i].setCrystalReciprocalSpace_fc(null);
            //refinementdata[i].setCrystalReciprocalSpace_fs(null);
            //crs_fc[i] = null;
            //crs_fs[i] = null;
        }

        tmprefinementmodel = null;
        fc = null;
        fs = null;
    }

    /**
     * move the atomic coordinates for the FFT calculation - this is independent
     * of the {@link Atom} object coordinates as it uses a linearized array
     *
     * @param x array of coordinates to move atoms to
     * @see CrystalReciprocalSpace#setCoordinates(double[])
     */
    public void setFFTCoordinates(double x[]) {
        for (int i = 0; i < n; i++) {
            crs_fc[i].setCoordinates(x);
            crs_fs[i].setCoordinates(x);
        }
    }

    /**
     * determine the total atomic gradients due to the diffraction data -
     * performs an inverse FFT
     *
     * @param refinementMode the
     * {@link RefinementMinimize.RefinementMode refinement mode} requested
     * @see CrystalReciprocalSpace#computeAtomicGradients(double[][], int[],
     * int, ffx.xray.RefinementMinimize.RefinementMode, boolean)
     * computeAtomicGradients
     */
    public void computeAtomicGradients(RefinementMode refinementMode) {
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeAtomicGradients(refinementData[i].dfc,
                    refinementData[i].freer, refinementData[i].rfreeflag,
                    refinementMode);
            crs_fs[i].computeAtomicGradients(refinementData[i].dfs,
                    refinementData[i].freer, refinementData[i].rfreeflag,
                    refinementMode);
        }
    }

    /**
     * parallelized call to compute atomic density on a grid, followed by FFT to
     * compute structure factors.
     *
     * @see CrystalReciprocalSpace#computeDensity(double[][], boolean)
     */
    public void computeAtomicDensity() {
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeDensity(refinementData[i].fc);
            if (solventModel != SolventModel.NONE) {
                crs_fs[i].computeDensity(refinementData[i].fs);
            }
        }
    }

    /**
     * compute total crystallographic likelihood target
     *
     * @return the total -log likelihood
     * @see SigmaAMinimize#calculateLikelihood()
     */
    public double computeLikelihood() {
        double e = 0.0;
        for (int i = 0; i < n; i++) {
            e += dataFiles[i].weight * sigmaAMinimize[i].calculateLikelihood();
        }
        return e;
    }

    /**
     * {@inheritDoc}

 return the usedAtoms for the model associated with this data
     */
    @Override
    public Atom[] getAtomArray() {
        return refinementModel.usedAtoms;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return refinementModel.altResidues;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return refinementModel.altMolecules;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public MolecularAssembly[] getMolecularAssemblies() {
        return assembly;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public RefinementModel getRefinementModel() {
        return refinementModel;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getWeight() {
        return xWeight;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setWeight(double weight) {
        this.xWeight = weight;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printOptimizationHeader() {
        return "R  Rfree";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printOptimizationUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%6.2f %6.2f ",
                    crystalStats[i].getR(),
                    crystalStats[i].getRFree()));
        }

        return sb.toString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printEnergyUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("     dataset %d (weight: %5.1f): R: %6.2f Rfree: %6.2f chemical energy: %8.2f likelihood: %8.2f free likelihood: %8.2f\n",
                    i + 1,
                    dataFiles[i].weight,
                    crystalStats[i].getR(),
                    crystalStats[i].getRFree(),
                    assembly[0].getPotentialEnergy().getTotalEnergy(),
                    dataFiles[i].weight * sigmaAMinimize[i].calculateLikelihood(),
                    dataFiles[i].weight * refinementData[i].llkf));
        }

        return sb.toString();
    }

    /**
     * print scale and R statistics for all datasets associated with the model
     */
    public void printScaleAndR() {
        for (int i = 0; i < n; i++) {
            if (!scaled[i]) {
                scaleBulkFit(i);
            }
            crystalStats[i].printScaleStats();
            crystalStats[i].printRStats();
        }
    }

    /**
     * print all statistics for all datasets associated with the model
     */
    public void printStats() {
        int nat = 0;
        int nnonh = 0;
        for (Atom a : refinementModel.usedAtomList) {
            if (a.getOccupancy() == 0.0) {
                continue;
            }
            nat++;
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            nnonh++;
        }

        for (int i = 0; i < n; i++) {
            if (!scaled[i]) {
                scaleBulkFit(i);
            }

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(" Statistics for Data Set %d of %d\n\n"
                    + "  Weight:     %6.2f\n  Neutron data: %4s\n"
                    + "  Model:        %s\n  Data file:    %s\n",
                    i + 1, n, dataFiles[i].weight, dataFiles[i].neutron,
                    modelName, dataFiles[i].filename));
            logger.info(sb.toString());

            crystalStats[i].printScaleStats();
            crystalStats[i].printDPIStats(nnonh, nat);
            crystalStats[i].printHKLStats();
            crystalStats[i].printSNStats();
            crystalStats[i].printRStats();
        }
    }

    /**
     * scale model and fit bulk solvent to all data
     */
    public void scaleBulkFit() {
        for (int i = 0; i < n; i++) {
            scaleBulkFit(i);
        }
    }

    /**
     * scale model and fit bulk solvent to dataset i of n
     *
     * @param i a int.
     */
    public void scaleBulkFit(int i) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Scaling Data Set %d of %d\n\n  Weight: %6.2f\n  Neutron data: %s\n  Model: %s\n  Data file: %s",
                i + 1, n, dataFiles[i].weight, dataFiles[i].neutron,
                modelName, dataFiles[i].filename));
        logger.info(sb.toString());

        // reset some values
        refinementData[i].solvent_k = 0.33;
        refinementData[i].solvent_ueq = 50.0 / (8.0 * Math.PI * Math.PI);
        refinementData[i].model_k = 0.0;
        for (int j = 0; j < 6; j++) {
            refinementData[i].model_b[j] = 0.0;
        }

        // run FFTs
        crs_fc[i].computeDensity(refinementData[i].fc);
        if (solventModel != SolventModel.NONE) {
            crs_fs[i].computeDensity(refinementData[i].fs);
        }

        // initialize minimizers
        scaleBulkMinimize[i] = new ScaleBulkMinimize(reflectionList[i],
                refinementData[i], crs_fs[i], parallelTeam);
        splineMinimize[i] = new SplineMinimize(reflectionList[i],
                refinementData[i], refinementData[i].spline,
                SplineEnergy.Type.FOFC);

        // minimize
        if (solventModel != SolventModel.NONE && gridsearch) {
            scaleBulkMinimize[i].minimize(6, xrayScaleTol);
            scaleBulkMinimize[i].GridOptimize();
        }
        scaleBulkMinimize[i].minimize(6, xrayScaleTol);

        // sigmaA / LLK calculation
        sigmaAMinimize[i] = new SigmaAMinimize(reflectionList[i], refinementData[i], parallelTeam);
        sigmaAMinimize[i].minimize(7, sigmaATol);

        if (splinefit) {
            splineMinimize[i].minimize(7, 1e-5);
        }

        scaled[i] = true;
    }

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        for (int i = 0; i < n; i++) {
            crs_fc[i].setLambda(lambda);
            crs_fs[i].setLambda(lambda);
        }
    }

    /**
     * set the bulk solvent parameters for a given bulk solvent model
     *
     * @param a typically the width of the atom
     * @param b typically the rate with which the atom transitions to bulk
     * @see CrystalReciprocalSpace#setSolventAB(double, double)
     */
    public void setSolventAB(double a, double b) {
        for (int i = 0; i < n; i++) {
            if (solventModel != SolventModel.NONE) {
                refinementData[i].solvent_a = a;
                refinementData[i].solvent_b = b;
                crs_fs[i].setSolventAB(a, b);
            }
        }
    }

    /**
     * perform 10 Fc calculations for the purposes of timings
     */
    public void timings() {
        logger.info("performing 10 Fc calculations for timing...");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 10; j++) {
                crs_fc[i].computeDensity(refinementData[i].fc, true);
                crs_fs[i].computeDensity(refinementData[i].fs, true);
                crs_fc[i].computeAtomicGradients(refinementData[i].dfc,
                        refinementData[i].freer, refinementData[i].rfreeflag,
                        RefinementMode.COORDINATES, true);
                crs_fs[i].computeAtomicGradients(refinementData[i].dfs,
                        refinementData[i].freer, refinementData[i].rfreeflag,
                        RefinementMode.COORDINATES, true);
            }
        }
    }

    /**
     * write current model to PDB file
     *
     * @param filename output PDB filename
     */
    public void writeModel(String filename) {
        StringBuilder remark = new StringBuilder();
        File file = new File(filename);
        PDBFilter pdbFilter = new PDBFilter(file, Arrays.asList(assembly), null, null);

        Date now = new Date();
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
        remark.append("REMARK FFX output ISO-8601 date: " + sdf.format(now) + "\n");
        remark.append("REMARK\n");
        remark.append("REMARK   3\n");
        remark.append("REMARK   3 REFINEMENT\n");
        remark.append("REMARK   3   PROGRAM     : FORCE FIELD X\n");
        remark.append("REMARK   3\n");
        for (int i = 0; i < n; i++) {
            remark.append("REMARK   3  DATA SET " + (i + 1) + "\n");
            if (dataFiles[i].isNeutron()) {
                remark.append("REMARK   3   DATA SET TYPE   : NEUTRON\n");
            } else {
                remark.append("REMARK   3   DATA SET TYPE   : X-RAY\n");
            }
            remark.append("REMARK   3   DATA SET WEIGHT : " + dataFiles[i].getWeight() + "\n");
            remark.append("REMARK   3\n");
            remark.append(crystalStats[i].getPDBHeaderString());
        }
        for (int i = 0; i < assembly.length; i++) {
            remark.append("REMARK   3  CHEMICAL SYSTEM " + (i + 1) + "\n");
            remark.append(assembly[i].getPotentialEnergy().getPDBHeaderString());
        }
        pdbFilter.writeFileWithHeader(file, remark);
    }

    /**
     * write current datasets to MTZ files
     *
     * @param filename output filename, or filename root for multiple datasets
     * @see MTZWriter#write()
     */
    public void writeData(String filename) {
        if (n == 1) {
            writeData(filename, 0);
        } else {
            for (int i = 0; i < n; i++) {
                writeData("" + FilenameUtils.removeExtension(filename) + "_" + i + ".mtz", i);
            }
        }
    }

    /**
     * write dataset i to MTZ file
     *
     * @param filename output filename
     * @param i dataset to write out
     * @see MTZWriter#write()
     */
    public void writeData(String filename, int i) {
        MTZWriter mtzwriter;
        if (scaled[i]) {
            mtzwriter = new MTZWriter(reflectionList[i], refinementData[i], filename);
        } else {
            mtzwriter = new MTZWriter(reflectionList[i], refinementData[i], filename, MTZType.DATAONLY);
        }
        mtzwriter.write();
    }

    /**
     * write 2Fo-Fc and Fo-Fc maps for all datasets
     *
     * @param filename output root filename for Fo-Fc and 2Fo-Fc maps
     */
    public void writeMaps(String filename) {
        if (n == 1) {
            writeMaps(filename, 0);
        } else {
            for (int i = 0; i < n; i++) {
                writeMaps("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
            }
        }
    }

    /**
     * write 2Fo-Fc and Fo-Fc maps for a datasets
     *
     * @param filename output root filename for Fo-Fc and 2Fo-Fc maps
     * @param i a int.
     */
    public void writeMaps(String filename, int i) {
        if (!scaled[i]) {
            scaleBulkFit(i);
        }

        // Fo-Fc
        crs_fc[i].computeAtomicGradients(refinementData[i].fofc1,
                refinementData[i].freer, refinementData[i].rfreeflag,
                RefinementMode.COORDINATES);
        double[] densityGrid = crs_fc[i].densityGrid;
        int extx = (int) crs_fc[i].getXDim();
        int exty = (int) crs_fc[i].getYDim();
        int extz = (int) crs_fc[i].getZDim();

        CCP4MapWriter mapwriter = new CCP4MapWriter(extx, exty, extz,
                crystal[i], FilenameUtils.removeExtension(filename) + "_fofc.map");
        mapwriter.write(densityGrid);

        // 2Fo-Fc
        crs_fc[i].computeAtomicGradients(refinementData[i].fofc2,
                refinementData[i].freer, refinementData[i].rfreeflag,
                RefinementMode.COORDINATES);
        densityGrid = crs_fc[i].densityGrid;
        extx = (int) crs_fc[i].getXDim();
        exty = (int) crs_fc[i].getYDim();
        extz = (int) crs_fc[i].getZDim();

        mapwriter = new CCP4MapWriter(extx, exty, extz,
                crystal[i], FilenameUtils.removeExtension(filename) + "_2fofc.map");
        mapwriter.write(densityGrid);
    }

    /**
     * write bulk solvent mask for all datasets to a CNS map file
     *
     * @param filename output filename, or output root filename for multiple
     * datasets
     */
    public void writeSolventMaskCNS(String filename) {
        if (n == 1) {
            writeSolventMaskCNS(filename, 0);
        } else {
            for (int i = 0; i < n; i++) {
                writeSolventMaskCNS("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
            }
        }
    }

    /**
     * write bulk solvent mask for dataset i to a CNS map file
     *
     * @param filename output filename
     * @param i dataset to write out
     */
    public void writeSolventMaskCNS(String filename, int i) {
        if (solventModel != SolventModel.NONE) {
            try {
                PrintWriter cnsfile = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
                cnsfile.println(" ANOMalous=FALSE");
                cnsfile.println(" DECLare NAME=FS DOMAin=RECIprocal TYPE=COMP END");
                for (HKL ih : reflectionList[i].hkllist) {
                    int j = ih.index();
                    cnsfile.printf(" INDE %d %d %d FS= %.4f %.4f\n",
                            ih.h(), ih.k(), ih.l(),
                            refinementData[i].fsF(j),
                            Math.toDegrees(refinementData[i].fsPhi(j)));
                }
                cnsfile.close();
            } catch (Exception e) {
                String message = "Fatal exception evaluating structure factors.\n";
                logger.log(Level.SEVERE, message, e);
                System.exit(-1);
            }
        }
    }

    /**
     * write bulk solvent mask for all datasets to a CCP4 map file
     *
     * @param filename output filename, or output root filename for multiple
     * datasets
     * @see CCP4MapWriter#write(double[], boolean)
     */
    public void writeSolventMask(String filename) {
        if (n == 1) {
            writeSolventMask(filename, 0);
        } else {
            for (int i = 0; i < n; i++) {
                writeSolventMask("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
            }
        }
    }

    /**
     * write bulk solvent mask for dataset i to a CCP4 map file
     *
     * @param filename output filename
     * @param i dataset to write out
     * @see CCP4MapWriter#write(double[], boolean)
     */
    public void writeSolventMask(String filename, int i) {
        if (solventModel != SolventModel.NONE) {
            CCP4MapWriter mapwriter = new CCP4MapWriter((int) crs_fs[i].getXDim(),
                    (int) crs_fs[i].getYDim(), (int) crs_fs[i].getZDim(),
                    crystal[i], filename);
            mapwriter.write(crs_fs[i].solventGrid);
        }
    }
}
