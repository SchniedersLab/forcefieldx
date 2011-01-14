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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.List;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 *
 * @author Tim Fenn
 */
public class DiffractionData {

    private static final Logger logger = Logger.getLogger(DiffractionData.class.getName());
    protected final String modelname;
    protected final DiffractionFile dataname[];
    protected final int n;
    protected final Crystal crystal[];
    protected final Resolution resolution[];
    protected final ReflectionList reflectionlist[];
    protected final RefinementData refinementdata[];
    protected final CrystalReciprocalSpace crs_fc[];
    protected final CrystalReciprocalSpace crs_fs[];
    public final int solventmodel;
    protected List<Atom> atomlist;
    protected final Atom atomarray[];
    protected List<Integer> xindex[];
    protected ArrayList<ArrayList<Residue>> altresidues;
    protected ArrayList<ArrayList<Molecule>> altmolecules;
    protected ScaleBulkMinimize scalebulkminimize[];
    protected SigmaAMinimize sigmaaminimize[];
    protected SplineMinimize splineminimize[];
    protected CrystalStats crystalstats[];
    protected boolean scaled[];
    protected XRayEnergy xrayenergy = null;

    /**
     * construct a diffraction data assembly, assumes an X-ray data set with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly molecular assembly, used as the atomic model for
     * comparison against the data
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
     * @param assembly molecular assembly, used as the atomic model for
     * comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link DiffractionFile} to be refined against
     */
    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties,
                SolventModel.POLYNOMIAL, datafile);
    }

    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                new DiffractionFile(assembly));
    }

    public DiffractionData(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel,
            DiffractionFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel,
                datafile);
    }

    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties) {
        this(assembly, properties, SolventModel.POLYNOMIAL,
                new DiffractionFile(assembly[0]));
    }

    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, DiffractionFile... datafile) {
        this(assembly, properties, SolventModel.POLYNOMIAL, datafile);
    }

    public DiffractionData(MolecularAssembly assembly[],
            CompositeConfiguration properties, int solventmodel,
            DiffractionFile... datafile) {
        List<Atom> alist;

        this.solventmodel = solventmodel;
        this.modelname = assembly[0].getFile().getName();
        this.dataname = datafile;
        this.n = datafile.length;

        crystal = new Crystal[n];
        resolution = new Resolution[n];
        reflectionlist = new ReflectionList[n];
        refinementdata = new RefinementData[n];
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
                    logger.severe("MTZ/CIF/CNS file does not contain full crystal information!");
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
            refinementdata[i] = new RefinementData(properties, reflectionlist[i]);
            tmp = new File(datafile[i].filename);
            datafile[i].diffractionfilter.readFile(tmp, reflectionlist[i], refinementdata[i]);
        }

        // FIXME: assembly crystal can have replicates (and when PDB is written, too)
        if (!crystal[0].equals(assembly[0].getCrystal())) {
            // logger.severe("PDB and reflection file crystal information do not match! (check CRYST1 record?)");
        }

        // build alternate conformer list for occupancy refinement (if necessary)
        altresidues = new ArrayList<ArrayList<Residue>>();
        altmolecules = new ArrayList<ArrayList<Molecule>>();
        ArrayList<MSNode> nlist0 = assembly[0].getNodeList();
        ArrayList<Residue> rtmp = null;
        ArrayList<Molecule> mtmp = null;
        Residue r0 = null;
        Molecule m0 = null;
        boolean altconf;
        for (int i = 0; i < nlist0.size(); i++) {
            altconf = false;
            MSNode node = nlist0.get(i);
            if (node instanceof Residue) {
                r0 = (Residue) node;
                m0 = null;
                alist = r0.getAtomList();
                for (Atom a : alist) {
                    if (!a.getAltLoc().equals(' ')
                            || a.getOccupancy() < 1.0) {
                        rtmp = new ArrayList<Residue>();
                        rtmp.add(r0);
                        altresidues.add(rtmp);
                        altconf = true;
                        break;
                    }
                }
            } else if (node instanceof Molecule
                    && refinementdata[0].refinemolocc) {
                r0 = null;
                m0 = (Molecule) node;
                alist = m0.getAtomList();
                for (Atom a : alist) {
                    if (!a.getAltLoc().equals(' ')
                            || a.getOccupancy() < 1.0) {
                        rtmp = new ArrayList<Residue>();
                        rtmp.add(r0);
                        altresidues.add(rtmp);
                        altconf = true;
                        break;
                    }
                }
            }
            if (altconf) {
                for (int j = 1; j < assembly.length; j++) {
                    ArrayList<MSNode> nlist = assembly[j].getNodeList();
                    Residue r;
                    Molecule m;
                    if (r0 != null
                            && nlist.size() > i) {
                        r = (Residue) nlist.get(i);
                        alist = r.getAtomList();
                        for (Atom a : alist) {
                            if (!a.getAltLoc().equals(' ')
                                    || a.getOccupancy() < 1.0) {
                                if (rtmp != null) {
                                    rtmp.add(r);
                                }
                                break;
                            }
                        }
                    } else if (m0 != null
                            && nlist.size() > i) {
                        m = (Molecule) nlist.get(i);
                        alist = m.getAtomList();
                        for (Atom a : alist) {
                            if (!a.getAltLoc().equals(' ')
                                    || a.getOccupancy() < 1.0) {
                                if (mtmp != null) {
                                    mtmp.add(m);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        for (ArrayList<Residue> list : altresidues) {
            if (list.size() == 1) {
                Residue r = list.get(0);
                logger.info("residue: " + r.toString() + ": single conformer, non-unity occupancy: occupancy will be refined independently!");
            }
        }

        for (ArrayList<Molecule> list : altmolecules) {
            if (list.size() == 1) {
                Molecule m = list.get(0);
                logger.info("molecule: " + m.toString() + ": single conformer, non-unity occupancy: occupancy will be refined independently!");
            }
        }

        xindex = new List[assembly.length];
        for (int i = 0; i < assembly.length; i++) {
            xindex[i] = new ArrayList<Integer>();
        }
        atomlist = new ArrayList<Atom>();
        int index = 0;
        alist = assembly[0].getAtomList();
        for (Atom a : alist) {
            for (int i = 0; i < assembly.length; i++) {
                xindex[i].add(index);
            }
            index++;
            atomlist.add(a);
        }

        for (int i = 1; i < assembly.length; i++) {
            alist = assembly[i].getAtomList();
            int subindex = 0;
            for (Atom a : alist) {
                Character c = a.getAltLoc();
                if (c == null) {
                    continue;
                }
                if (!c.equals(' ')
                        && !c.equals('A')) {
                    xindex[i].set(subindex, index++);
                    atomlist.add(a);
                }
                subindex++;
            }
        }
        atomarray = atomlist.toArray(new Atom[atomlist.size()]);

        // initialize atomic form factors
        for (Atom a : atomarray) {
            XRayFormFactor atomff =
                    new XRayFormFactor(a, refinementdata[0].use_3g, 2.0);
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
                double rho = atomff.rho(0.0, xyz);
                if (rho > 0.1) {
                    arad += 0.5;
                } else if (rho > 0.001) {
                    arad += 0.1;
                } else {
                    arad += 0.2;
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
            crs_fc[i] = new CrystalReciprocalSpace(reflectionlist[i], atomarray,
                    parallelTeam, parallelTeam, false, dataname[i].neutron);
            refinementdata[i].setCrystalReciprocalSpace_fc(crs_fc[i]);
            crs_fc[i].setUse3G(refinementdata[i].use_3g);
            crs_fc[i].setWeight(dataname[i].weight);
            crs_fs[i] = new CrystalReciprocalSpace(reflectionlist[i], atomarray,
                    parallelTeam, parallelTeam, true, dataname[i].neutron,
                    solventmodel);
            refinementdata[i].setCrystalReciprocalSpace_fs(crs_fs[i]);
            crs_fs[i].setUse3G(refinementdata[i].use_3g);
            crs_fs[i].setWeight(dataname[i].weight);

            crystalstats[i] = new CrystalStats(reflectionlist[i],
                    refinementdata[i]);
        }

        scaled = new boolean[n];
        for (int i = 0; i < n; i++) {
            scaled[i] = false;
        }
    }

    public void setFFTCoordinates(double x[]) {
        for (int i = 0; i < n; i++) {
            crs_fc[i].setCoordinates(x);
            crs_fs[i].setCoordinates(x);
        }
    }

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

    public void computeAtomicDensity() {
        for (int i = 0; i < n; i++) {
            crs_fc[i].computeDensity(refinementdata[i].fc);
            crs_fs[i].computeDensity(refinementdata[i].fs);
        }
    }

    public double computeLikelihood() {
        double e = 0.0;
        for (int i = 0; i < n; i++) {
            e += dataname[i].weight * sigmaaminimize[i].calculateLikelihood();
        }
        return e;
    }

    public XRayEnergy getXRayEnergy() {
        return xrayenergy;
    }

    public void setXRayEnergy(XRayEnergy xrayenergy) {
        this.xrayenergy = xrayenergy;
    }

    public Atom[] getAtomArray() {
        return atomarray;
    }

    public void printscaleandr() {
        for (int i = 0; i < n; i++) {
            if (!scaled[i]) {
                scalebulkfit(i);
            }
            crystalstats[i].print_scalestats();
            crystalstats[i].print_rstats();
        }
    }

    public void printstats() {
        int nat = 0;
        int nnonh = 0;
        for (Atom a : atomlist) {
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
                scalebulkfit(i);
            }

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("statistics for data set %d of %d\nweight: %6.2f is neutron: %s\nmodel: %s data file: %s\n",
                    i + 1, n, dataname[i].weight, dataname[i].neutron,
                    modelname, dataname[i].filename));
            logger.info(sb.toString());

            crystalstats[i].print_scalestats();
            crystalstats[i].print_dpistats(nnonh, nat);
            crystalstats[i].print_hklstats();
            crystalstats[i].print_snstats();
            crystalstats[i].print_rstats();
        }
    }

    public void scalebulkfit() {
        for (int i = 0; i < n; i++) {
            scalebulkfit(i);
        }
    }

    public void scalebulkfit(int i) {
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
        if (solventmodel != SolventModel.NONE
                && refinementdata[i].gridsearch) {
            scalebulkminimize[i].minimize(7, 1e-2);
            scalebulkminimize[i].GridOptimize();
        }
        scalebulkminimize[i].minimize(7, refinementdata[i].xrayscaletol);

        // sigmaA / LLK calculation
        sigmaaminimize[i] = new SigmaAMinimize(reflectionlist[i], refinementdata[i]);
        sigmaaminimize[i].minimize(7, refinementdata[i].sigmaatol);

        if (refinementdata[i].splinefit) {
            splineminimize[i].minimize(7, 1e-5);
        }

        scaled[i] = true;
    }

    public void setSolventAB(double a, double b) {
        for (int i = 0; i < n; i++) {
            if (solventmodel != SolventModel.NONE) {
                refinementdata[i].solvent_a = a;
                refinementdata[i].solvent_b = b;
                crs_fs[i].setSolventAB(a, b);
            }
        }
    }

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

    public void writedata(String filename) {
        writedata(filename, 0);
    }

    public void writedata(String filename, int i) {
        MTZWriter mtzwriter;
        if (scaled[i]) {
            mtzwriter = new MTZWriter(reflectionlist[i], refinementdata[i], filename);
        } else {
            mtzwriter = new MTZWriter(reflectionlist[i], refinementdata[i], filename, true);
        }
        mtzwriter.write();
    }

    public void writeSolventMaskCNS(String filename) {
        writeSolventMaskCNS(filename, 0);
    }

    public void writeSolventMaskCNS(String filename, int i) {
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

    public void writeSolventMask(String filename) {
        writeSolventMask(filename, 0);
    }

    public void writeSolventMask(String filename, int i) {
        if (solventmodel != SolventModel.NONE) {
            CCP4MapWriter mapwriter = new CCP4MapWriter((int) crs_fs[i].getXDim(),
                    (int) crs_fs[i].getYDim(), (int) crs_fs[i].getZDim(),
                    crystal[i], filename);
            mapwriter.write(crs_fs[i].solventGrid);
        }
    }
}
