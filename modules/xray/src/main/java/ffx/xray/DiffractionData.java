/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.xray;

import ffx.crystal.CCP4MapWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.parsers.PDBFilter;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 *
 * @author Tim Fenn
 */
public class DiffractionData implements DataContainer {

    private static final Logger logger = Logger.getLogger(DiffractionData.class.getName());
    protected final MolecularAssembly assembly[];
    protected final String modelname;
    protected final DiffractionFile dataname[];
    protected final int n;
    protected final Crystal crystal[];
    protected final Resolution resolution[];
    protected final ReflectionList reflectionlist[];
    protected final DiffractionRefinementData refinementdata[];
    protected final CrystalReciprocalSpace crs_fc[];
    protected final CrystalReciprocalSpace crs_fs[];
    public final int solventmodel;
    protected final RefinementModel refinementmodel;
    protected ScaleBulkMinimize scalebulkminimize[];
    protected SigmaAMinimize sigmaaminimize[];
    protected SplineMinimize splineminimize[];
    protected CrystalStats crystalstats[];
    protected boolean scaled[];
    // settings
    public int rfreeflag;
    public final double fsigfcutoff;
    public final boolean use_3g;
    public final double aradbuff;
    public final double xrayscaletol;
    public final double sigmaatol;
    public final double xweight;
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
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties,
                SolventModel.POLYNOMIAL, new DiffractionFile(assembly));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties,
                SolventModel.POLYNOMIAL, datafile);
    }

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired -
     * see {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                new DiffractionFile(assembly));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired -
     * see {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel,
            DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                datafile);
    }

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object array (typically containing alternate conformer assemblies), used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties) {
        this(assembly, properties, SolventModel.POLYNOMIAL,
                new DiffractionFile(assembly[0]));
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object array (typically containing alternate conformer assemblies), used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(assembly, properties, SolventModel.POLYNOMIAL, datafile);
    }

    /**
     * construct a diffraction data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object array (typically containing alternate conformer assemblies), used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param solventmodel the type of solvent model desired -
     * see {@link CrystalReciprocalSpace.SolventModel bulk solvent model} selections
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, int solventmodel,
            DiffractionFile... datafile) {

        this.assembly = assembly;
        this.solventmodel = solventmodel;
        this.modelname = assembly[0].getFile().getName();
        this.dataname = datafile;
        this.n = datafile.length;

        int rflag = properties.getInt("rfreeflag", -1);
        fsigfcutoff = properties.getDouble("fsigfcutoff", -1.0);
        gridsearch = properties.getBoolean("gridsearch", false);
        splinefit = properties.getBoolean("splinefit", true);
        use_3g = properties.getBoolean("use_3g", true);
        aradbuff = properties.getDouble("aradbuff", 0.5);
        xrayscaletol = properties.getDouble("xrayscaletol", 1e-4);
        sigmaatol = properties.getDouble("sigmaatol", 1.0);
        xweight = properties.getDouble("xweight", 1.0);
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
        reflectionlist = new ReflectionList[n];
        refinementdata = new DiffractionRefinementData[n];
        scalebulkminimize = new ScaleBulkMinimize[n];
        sigmaaminimize = new SigmaAMinimize[n];
        splineminimize = new SplineMinimize[n];
        crystalstats = new CrystalStats[n];

        // read in Fo/sigFo/FreeR
        File tmp;
        Crystal crystalinit = Crystal.checkProperties(properties);
        Resolution resolutioninit = Resolution.checkProperties(properties);
        if (crystalinit == null || resolutioninit == null) {
            for (int i = 0; i < n; i++) {
                tmp = new File(datafile[i].filename);
                reflectionlist[i] = datafile[i].diffractionfilter.getReflectionList(tmp, properties);

                if (reflectionlist[i] == null) {
                    logger.info("Using crystal information from molecular assembly to generate crystal information");
                    crystalinit = assembly[i].getCrystal().getUnitCell();
                    double res = datafile[i].diffractionfilter.getResolution(tmp, crystalinit);
                    if (res < 0.0) {
                        logger.severe("MTZ/CIF/CNS file does not contain full crystal information!");
                    } else {
                        resolutioninit = new Resolution(res);
                        reflectionlist[i] = new ReflectionList(crystalinit, resolutioninit, properties);
                    }
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                reflectionlist[i] = new ReflectionList(crystalinit, resolutioninit,
                        properties);
            }
        }

        for (int i = 0; i < n; i++) {
            crystal[i] = reflectionlist[i].crystal;
            resolution[i] = reflectionlist[i].resolution;
            refinementdata[i] = new DiffractionRefinementData(properties, reflectionlist[i]);
            tmp = new File(datafile[i].filename);
            datafile[i].diffractionfilter.readFile(tmp, reflectionlist[i],
                    refinementdata[i], properties);
        }

        if (!crystal[0].equals(assembly[0].getCrystal().getUnitCell())) {
            logger.severe("PDB and reflection file crystal information do not match! (check CRYST1 record?)");
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  using cctbx 3 Gaussians (use_3g): " + use_3g + "\n");
            sb.append("  atomic form factor radius buffer (aradbuff): " + aradbuff + "\n");
            sb.append("  resolution dependent spline scale (splinefit): " + splinefit + "\n");
            sb.append("  F/sigF cutoff (fsigfcutoff): " + fsigfcutoff + "\n");
            sb.append("  R Free flag (rfreeflag) (if -1, value will be updated when data is read in): " + rflag + "\n");
            sb.append("  n bins (nbins): " + reflectionlist[0].nbins + "\n");
            sb.append("  solvent grid search (gridsearch): " + gridsearch + "\n");
            sb.append("  X-ray scale fit tolerance (xrayscaletol): " + xrayscaletol + "\n");
            sb.append("  sigma A fit tolerance (sigmaatol): " + sigmaatol + "\n");
            sb.append("  X-ray refinement weight (xweight): " + xweight + "\n");
            sb.append("  B similarity weight (bsimweight): " + bsimweight + "\n");
            sb.append("  B non-zero weight (bnonzeroweight): " + bnonzeroweight + "\n");
            sb.append("  B Lagrangian mass (bmass): " + bmass + "\n");
            sb.append("  B factors refined by residue (residuebfactor): " + residuebfactor + "\n");
            sb.append("    (if true, num. residues per B (nresiduebfactor): " + nresiduebfactor + ")\n");
            sb.append("  add ANISOU for refinement (addanisou): " + addanisou + "\n");
            sb.append("  refine occupancies on molecules (HETATMs - refinemolocc): " + refinemolocc + "\n");
            sb.append("  Occupancy Lagrangian mass (occmass): " + occmass + "\n");
            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementmodel = new RefinementModel(assembly, refinemolocc);

        // initialize atomic form factors
        for (Atom a : refinementmodel.atomarray) {
            a.setFormFactorIndex(-1);
            XRayFormFactor atomff =
                    new XRayFormFactor(a, use_3g, 2.0);
            a.setFormFactorIndex(atomff.ffindex);

            if (a.getOccupancy() == 0.0) {
                a.setFormFactorWidth(1.0);
                continue;
            }

            double arad = 2.4;
            // double arad = 2.0;
            double xyz[] = new double[3];
            xyz[0] = a.getX() + arad;
            xyz[1] = a.getY();
            xyz[2] = a.getZ();
            while (true) {
                double rho = atomff.rho(0.0, 1.0, xyz);
                if (rho > 0.1) {
                    arad += 0.5;
                } else if (rho > 0.001) {
                    arad += 0.1;
                } else {
                    arad += aradbuff;
                    a.setFormFactorWidth(arad);
                    break;
                }
                xyz[0] = a.getX() + arad;
            }
        }

        // set up FFT and run it
        crs_fc = new CrystalReciprocalSpace[n];
        crs_fs = new CrystalReciprocalSpace[n];
        ParallelTeam parallelTeam = new ParallelTeam();
        for (int i = 0; i < n; i++) {
            crs_fc[i] = new CrystalReciprocalSpace(reflectionlist[i],
                    refinementmodel.atomarray, parallelTeam, parallelTeam,
                    false, dataname[i].neutron);
            refinementdata[i].setCrystalReciprocalSpace_fc(crs_fc[i]);
            crs_fc[i].setUse3G(use_3g);
            crs_fc[i].setWeight(dataname[i].weight);
            crs_fc[i].lambdaTerm = lambdaTerm;
            crs_fs[i] = new CrystalReciprocalSpace(reflectionlist[i],
                    refinementmodel.atomarray, parallelTeam, parallelTeam,
                    true, dataname[i].neutron, solventmodel);
            refinementdata[i].setCrystalReciprocalSpace_fs(crs_fs[i]);
            crs_fs[i].setUse3G(use_3g);
            crs_fs[i].setWeight(dataname[i].weight);
            crs_fs[i].lambdaTerm = lambdaTerm;

            crystalstats[i] = new CrystalStats(reflectionlist[i],
                    refinementdata[i]);
        }

        scaled = new boolean[n];
        for (int i = 0; i < n; i++) {
            scaled[i] = false;
        }
    }

    /**
     * move the atomic coordinates for the FFT calculation - this is independent
     * of the {@link Atom} object coordinates as it uses a linearized array
     *
     * @param x array of coordinates to move atoms to
     *
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
     *
     * @see CrystalReciprocalSpace#computeAtomicGradients(double[][], int[], int, ffx.xray.RefinementMinimize.RefinementMode, boolean) computeAtomicGradients
     */
    public void computeAtomicGradients(RefinementMode refinementMode) {
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeAtomicGradients(refinementdata[i].dfc,
                    refinementdata[i].freer, refinementdata[i].rfreeflag,
                    refinementMode);
            crs_fs[i].computeAtomicGradients(refinementdata[i].dfs,
                    refinementdata[i].freer, refinementdata[i].rfreeflag,
                    refinementMode);
        }
    }

    /**
     * parallelized call to compute atomic density on a grid, followed by FFT
     * to compute structure factors.
     *
     * @see CrystalReciprocalSpace#computeDensity(double[][], boolean)
     */
    public void computeAtomicDensity() {
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeDensity(refinementdata[i].fc);
            crs_fs[i].computeDensity(refinementdata[i].fs);
        }
    }

    /**
     * compute total crystallographic likelihood target
     *
     * @return the total -log likelihood
     *
     * @see SigmaAMinimize#calculateLikelihood()
     */
    public double computeLikelihood() {
        double e = 0.0;
        for (int i = 0; i < n; i++) {
            e += dataname[i].weight * sigmaaminimize[i].calculateLikelihood();
        }
        return e;
    }

    /**
     * return the atomarray for the model associated with this data
     *
     * @return an {@link ffx.potential.bonded.Atom} array
     */
    @Override
    public Atom[] getAtomArray() {
        return refinementmodel.atomarray;
    }

    @Override
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return refinementmodel.altresidues;
    }

    @Override
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return refinementmodel.altmolecules;
    }

    @Override
    public MolecularAssembly[] getMolecularAssembly() {
        return assembly;
    }

    @Override
    public RefinementModel getRefinementModel() {
        return refinementmodel;
    }

    @Override
    public double getWeight() {
        return xweight;
    }

    @Override
    public String printOptimizationHeader() {
        return "R  Rfree";
    }

    @Override
    public String printOptimizationUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%6.2f %6.2f ",
                    crystalstats[i].getR(),
                    crystalstats[i].getRFree()));
        }

        return sb.toString();
    }

    @Override
    public String printEnergyUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("     dataset %d (weight: %5.1f): R: %6.2f Rfree: %6.2f chemical energy: %8.2f likelihood: %8.2f free likelihood: %8.2f\n",
                    i + 1,
                    dataname[i].weight,
                    crystalstats[i].getR(),
                    crystalstats[i].getRFree(),
                    assembly[0].getPotentialEnergy().getTotal(),
                    dataname[i].weight * sigmaaminimize[i].calculateLikelihood(),
                    dataname[i].weight * refinementdata[i].llkf));
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
            crystalstats[i].printScaleStats();
            crystalstats[i].printRStats();
        }
    }

    /**
     * print all statistics for all datasets associated with the model
     */
    public void printStats() {
        int nat = 0;
        int nnonh = 0;
        for (Atom a : refinementmodel.atomlist) {
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
            sb.append(String.format("statistics for data set %d of %d\nweight: %6.2f is neutron: %s\nmodel: %s data file: %s\n",
                    i + 1, n, dataname[i].weight, dataname[i].neutron,
                    modelname, dataname[i].filename));
            logger.info(sb.toString());

            crystalstats[i].printScaleStats();
            crystalstats[i].printDPIStats(nnonh, nat);
            crystalstats[i].printHKLStats();
            crystalstats[i].printSNStats();
            crystalstats[i].printRStats();
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
     */
    public void scaleBulkFit(int i) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("scaling data set %d of %d\nweight: %6.2f is neutron: %s\nmodel: %s data file: %s\n",
                i + 1, n, dataname[i].weight, dataname[i].neutron,
                modelname, dataname[i].filename));
        logger.info(sb.toString());

        // reset some values
        refinementdata[i].solvent_k = 0.33;
        refinementdata[i].solvent_ueq = 50.0 / (8.0 * Math.PI * Math.PI);
        refinementdata[i].model_k = 0.0;
        for (int j = 0; j < 6; j++) {
            refinementdata[i].model_b[j] = 0.0;
        }

        // run FFTs
        crs_fc[i].computeDensity(refinementdata[i].fc);
        if (solventmodel != SolventModel.NONE) {
            crs_fs[i].computeDensity(refinementdata[i].fs);
        }

        // initialize minimizers
        scalebulkminimize[i] = new ScaleBulkMinimize(reflectionlist[i],
                refinementdata[i], crs_fs[i]);
        splineminimize[i] = new SplineMinimize(reflectionlist[i],
                refinementdata[i], refinementdata[i].spline,
                SplineEnergy.Type.FOFC);

        // minimize
        if (solventmodel != SolventModel.NONE && gridsearch) {
            scalebulkminimize[i].minimize(7, 1e-2);
            scalebulkminimize[i].GridOptimize();
        }
        scalebulkminimize[i].minimize(7, xrayscaletol);

        // sigmaA / LLK calculation
        sigmaaminimize[i] = new SigmaAMinimize(reflectionlist[i], refinementdata[i]);
        sigmaaminimize[i].minimize(7, sigmaatol);

        if (splinefit) {
            splineminimize[i].minimize(7, 1e-5);
        }

        scaled[i] = true;
    }

    /**
     * Set the current value of the state variable.
     * 
     * @param lambda
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
     *
     * @see CrystalReciprocalSpace#setSolventAB(double, double)
     */
    public void setSolventAB(double a, double b) {
        for (int i = 0; i < n; i++) {
            if (solventmodel != SolventModel.NONE) {
                refinementdata[i].solvent_a = a;
                refinementdata[i].solvent_b = b;
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
                crs_fc[i].computeDensity(refinementdata[i].fc, true);
                crs_fs[i].computeDensity(refinementdata[i].fs, true);
                crs_fc[i].computeAtomicGradients(refinementdata[i].dfc,
                        refinementdata[i].freer, refinementdata[i].rfreeflag,
                        RefinementMode.COORDINATES, true);
                crs_fs[i].computeAtomicGradients(refinementdata[i].dfs,
                        refinementdata[i].freer, refinementdata[i].rfreeflag,
                        RefinementMode.COORDINATES, true);
            }
        }
    }

    /**
     * write current model to PDB file
     * 
     * @param filename output PDB filename
     * 
     * @see PDBFilter#writeFile()
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
            if (dataname[i].isNeutron()) {
                remark.append("REMARK   3   DATA SET TYPE   : NEUTRON\n");
            } else {
                remark.append("REMARK   3   DATA SET TYPE   : X-RAY\n");
            }
            remark.append("REMARK   3   DATA SET WEIGHT : " + dataname[i].getWeight() + "\n");
            remark.append("REMARK   3\n");
            remark.append(crystalstats[i].getPDBHeaderString());
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
     *
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
     *
     * @see MTZWriter#write()
     */
    public void writeData(String filename, int i) {
        MTZWriter mtzwriter;
        if (scaled[i]) {
            mtzwriter = new MTZWriter(reflectionlist[i], refinementdata[i], filename);
        } else {
            mtzwriter = new MTZWriter(reflectionlist[i], refinementdata[i], filename, true);
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
     */
    public void writeMaps(String filename, int i) {
        if (!scaled[i]) {
            scaleBulkFit(i);
        }

        // Fo-Fc
        crs_fc[i].computeAtomicGradients(refinementdata[i].fofc1,
                refinementdata[i].freer, refinementdata[i].rfreeflag,
                RefinementMode.COORDINATES);
        double[] densityGrid = crs_fc[i].densityGrid;
        int extx = (int) crs_fc[i].getXDim();
        int exty = (int) crs_fc[i].getYDim();
        int extz = (int) crs_fc[i].getZDim();

        CCP4MapWriter mapwriter = new CCP4MapWriter(extx, exty, extz,
                crystal[i], FilenameUtils.removeExtension(filename) + "_fofc.map");
        mapwriter.write(densityGrid);

        // 2Fo-Fc
        crs_fc[i].computeAtomicGradients(refinementdata[i].fofc2,
                refinementdata[i].freer, refinementdata[i].rfreeflag,
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
        if (solventmodel != SolventModel.NONE) {
            try {
                PrintWriter cnsfile = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
                cnsfile.println(" ANOMalous=FALSE");
                cnsfile.println(" DECLare NAME=FS DOMAin=RECIprocal TYPE=COMP END");
                for (HKL ih : reflectionlist[i].hkllist) {
                    int j = ih.index();
                    cnsfile.printf(" INDE %d %d %d FS= %.4f %.4f\n",
                            ih.h(), ih.k(), ih.l(),
                            refinementdata[i].fs_f(j),
                            Math.toDegrees(refinementdata[i].fs_phi(j)));
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
     *
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
     *
     * @see CCP4MapWriter#write(double[], boolean)
     */
    public void writeSolventMask(String filename, int i) {
        if (solventmodel != SolventModel.NONE) {
            CCP4MapWriter mapwriter = new CCP4MapWriter((int) crs_fs[i].getXDim(),
                    (int) crs_fs[i].getYDim(), (int) crs_fs[i].getZDim(),
                    crystal[i], filename);
            mapwriter.write(crs_fs[i].solventGrid);
        }
    }
}
