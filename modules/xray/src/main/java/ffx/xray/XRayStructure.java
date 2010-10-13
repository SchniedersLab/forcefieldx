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

import static org.apache.commons.io.FilenameUtils.removeExtension;

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
public class XRayStructure {

    private static final Logger logger = Logger.getLogger(XRayStructure.class.getName());
    protected final Crystal crystal;
    protected final Resolution resolution;
    protected final ReflectionList reflectionlist;
    protected final RefinementData refinementdata;
    protected final CrystalReciprocalSpace crs_fc;
    protected final CrystalReciprocalSpace crs_fs;
    public final int solventmodel;
    protected List<Atom> atomlist;
    protected final Atom atomarray[];
    protected List<Integer> xindex[];
    protected ArrayList<ArrayList<Residue>> altresidues;
    protected ArrayList<ArrayList<Molecule>> altmolecules;
    protected ScaleBulkMinimize scalebulkminimize;
    protected SigmaAMinimize sigmaaminimize;
    protected SplineMinimize splineminimize;
    protected CrystalStats crystalstats;
    protected boolean scaled = false;
    protected XRayEnergy xrayenergy = null;

    public XRayStructure(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties, SolventModel.POLYNOMIAL);
    }

    public XRayStructure(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel) {
        this(new MolecularAssembly[]{assembly}, properties, solventmodel);
    }

    public XRayStructure(MolecularAssembly assembly[],
            CompositeConfiguration properties) {
        this(assembly, properties, SolventModel.POLYNOMIAL);
    }

    public XRayStructure(MolecularAssembly assembly[],
            CompositeConfiguration properties, int solventmodel) {
        // String name = assembly.getName();
        String name = assembly[0].getFile().getPath();
        name = removeExtension(name);
        this.solventmodel = solventmodel;

        // load the structure
        File tmp = new File(name + ".mtz");
        File mtzfile = null, ciffile = null, cnsfile = null;
        if (tmp.exists()) {
            logger.info("data file: " + tmp.getName());
            mtzfile = tmp;
        } else {
            tmp = new File(name + ".cif");
            if (tmp.exists()) {
                logger.info("data file: " + tmp.getName());
                ciffile = tmp;
            } else {
                tmp = new File(name + ".ent");
                if (tmp.exists()) {
                    logger.info("data file: " + tmp.getName());
                    ciffile = tmp;
                } else {
                    tmp = new File(name + ".cns");
                    if (tmp.exists()) {
                        logger.info("data file: " + tmp.getName());
                        cnsfile = tmp;
                    } else {
                        tmp = new File(name + ".hkl");
                        if (tmp.exists()) {
                            logger.info("data file: " + tmp.getName());
                            cnsfile = tmp;
                        } else {
                            logger.severe("no input data found!");
                        }
                    }
                }
            }
        }

        // read in Fo/sigFo/FreeR
        MTZFilter mtzfilter = new MTZFilter();
        CIFFilter ciffilter = new CIFFilter();
        CNSFilter cnsfilter = new CNSFilter();
        Crystal crystalinit = Crystal.checkProperties(properties);
        Resolution resolutioninit = Resolution.checkProperties(properties);
        if (crystalinit == null || resolutioninit == null) {
            if (mtzfile != null) {
                reflectionlist = mtzfilter.getReflectionList(mtzfile, properties);
            } else if (ciffile != null) {
                reflectionlist = ciffilter.getReflectionList(ciffile, properties);
            } else if (cnsfile != null) {
                reflectionlist = cnsfilter.getReflectionList(cnsfile, properties);
            } else {
                reflectionlist = null;
            }

            if (reflectionlist == null) {
                logger.severe("MTZ/CIF/CNS file does not contain full crystal information!");
            }
        } else {
            reflectionlist = new ReflectionList(crystalinit, resolutioninit);
        }

        crystal = reflectionlist.crystal;
        resolution = reflectionlist.resolution;
        refinementdata = new RefinementData(properties, reflectionlist);
        if (mtzfile != null) {
            mtzfilter.readFile(mtzfile, reflectionlist, refinementdata);
        } else if (ciffile != null) {
            ciffilter.readFile(ciffile, reflectionlist, refinementdata);
        } else {
            cnsfilter.readFile(cnsfile, reflectionlist, refinementdata);
        }

        // FIXME: assembly crystal can have replicates (and when PDB is written, too)
        if (!crystal.equals(assembly[0].getCrystal())) {
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
        Character c0 = ' ';
        double occ;
        for (int i = 0; i < nlist0.size(); i++) {
            MSNode node = nlist0.get(i);
            if (node instanceof Residue) {
                r0 = (Residue) node;
                m0 = null;
                c0 = r0.getAtomList().get(0).getAltLoc();
                occ = r0.getAtomList().get(0).getOccupancy();
                if (!c0.equals(' ')
                        || occ < 1.0) {
                    rtmp = new ArrayList<Residue>();
                    rtmp.add(r0);
                    altresidues.add(rtmp);
                }
            } else if (node instanceof Molecule
                    && refinementdata.refinemolocc) {
                r0 = null;
                m0 = (Molecule) node;
                c0 = m0.getAtomList().get(0).getAltLoc();
                occ = m0.getAtomList().get(0).getOccupancy();
                if (!c0.equals(' ')
                        || occ < 1.0) {
                    mtmp = new ArrayList<Molecule>();
                    mtmp.add(m0);
                    altmolecules.add(mtmp);
                }
            } else {
                c0 = ' ';
            }
            if (!c0.equals(' ')) {
                for (int j = 1; j < assembly.length; j++) {
                    ArrayList<MSNode> nlist = assembly[j].getNodeList();
                    Residue r;
                    Molecule m;
                    Character c = 'A';
                    if (r0 != null
                            && nlist.size() > i) {
                        r = (Residue) nlist.get(i);
                        c = r.getAtomList().get(0).getAltLoc();
                        if (!c.equals(' ')
                                && !c.equals('A')) {
                            if (rtmp != null) {
                                rtmp.add(r);
                            }
                        }
                    } else if (m0 != null
                            && nlist.size() > i) {
                        m = (Molecule) nlist.get(i);
                        c = m.getAtomList().get(0).getAltLoc();
                        if (!c.equals(' ')
                                && !c.equals('A')) {
                            if (mtmp != null) {
                                mtmp.add(m);
                            }
                        }
                    }
                }
            }
        }

        for (ArrayList<Residue> list : altresidues) {
            if (list.size() == 1){
                Residue r = list.get(0);
                logger.info("residue: " + r.toString() + ": single conformer, non-unity occupancy: occupancy will be refined independently!");
            }
        }

        for (ArrayList<Molecule> list : altmolecules) {
            if (list.size() == 1){
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
        List<Atom> alist = assembly[0].getAtomList();
        for (Atom a : alist) {
            for (int i = 0; i < assembly.length; i++) {
                xindex[i].add(index);
            }
            index++;
            atomlist.add(a);
        }

        for (int i = 1; i < assembly.length; i++) {
            alist = assembly[i].getAtomList();
            for (Atom a : alist) {
                Character c = a.getAltLoc();
                if (c == null) {
                    continue;
                }
                if (!c.equals(' ')
                        && !c.equals('A')) {
                    xindex[i].add(index++);
                    atomlist.add(a);
                }
            }
        }
        atomarray = atomlist.toArray(new Atom[atomlist.size()]);

        // initialize atomic form factors
        for (Atom a : atomarray) {
            XRayFormFactor atomff =
                    new XRayFormFactor(a, refinementdata.use_3g, 2.0);
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
        ParallelTeam parallelTeam = new ParallelTeam();
        crs_fc = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, false);
        refinementdata.setCrystalReciprocalSpace_fc(crs_fc);
        crs_fc.setUse3G(refinementdata.use_3g);
        crs_fs = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, true, solventmodel);
        refinementdata.setCrystalReciprocalSpace_fs(crs_fs);
        crs_fs.setUse3G(refinementdata.use_3g);

        crystalstats = new CrystalStats(reflectionlist, refinementdata);
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
        if (!scaled) {
            scalebulkfit();
        }
        crystalstats.print_scalestats();
        crystalstats.print_rstats();
    }

    public void printstats() {
        if (!scaled) {
            scalebulkfit();
        }

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
        crystalstats.print_scalestats();
        crystalstats.print_dpistats(nnonh, nat);
        crystalstats.print_hklstats();
        crystalstats.print_snstats();
        crystalstats.print_rstats();
    }

    public void scalebulkfit() {
        // reset some values
        refinementdata.solvent_k = 0.33;
        refinementdata.solvent_ueq = 50.0 / (8.0 * Math.PI * Math.PI);
        refinementdata.model_k = 0.0;
        for (int i = 0; i < 6; i++) {
            refinementdata.model_b[i] = 0.0;
        }

        // run FFTs
        crs_fc.computeDensity(refinementdata.fc);
        if (solventmodel != SolventModel.NONE) {
            crs_fs.computeDensity(refinementdata.fs);
        }

        // initialize minimizers
        scalebulkminimize = new ScaleBulkMinimize(reflectionlist, refinementdata, crs_fs);
        splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.spline, SplineEnergy.Type.FOFC);

        // minimize
        if (solventmodel != SolventModel.NONE
                && refinementdata.gridsearch) {
            scalebulkminimize.minimize(7, 1e-2);
            scalebulkminimize.GridOptimize();
        }
        scalebulkminimize.minimize(7, refinementdata.xrayscaletol);

        // sigmaA / LLK calculation
        sigmaaminimize = new SigmaAMinimize(reflectionlist, refinementdata);
        sigmaaminimize.minimize(7, refinementdata.sigmaatol);

        if (refinementdata.splinefit) {
            splineminimize.minimize(7, 1e-5);
        }

        scaled = true;
    }

    public void setSolventAB(double a, double b) {
        if (solventmodel != SolventModel.NONE) {
            refinementdata.solvent_a = a;
            refinementdata.solvent_b = b;
            crs_fs.setSolventAB(a, b);
        }
    }

    public void timings() {
        logger.info("performing 10 Fc calculations for timing...");
        for (int i = 0; i < 10; i++) {
            crs_fc.computeDensity(refinementdata.fc, true);
            crs_fs.computeDensity(refinementdata.fs, true);
            crs_fc.computeAtomicGradients(refinementdata.dfc,
                    refinementdata.freer, refinementdata.rfreeflag,
                    RefinementMode.COORDINATES, true);
            crs_fs.computeAtomicGradients(refinementdata.dfs,
                    refinementdata.freer, refinementdata.rfreeflag,
                    RefinementMode.COORDINATES, true);
        }
    }

    public void writedata(String filename) {
        MTZWriter mtzwriter;
        if (scaled) {
            mtzwriter = new MTZWriter(reflectionlist, refinementdata, filename);
        } else {
            mtzwriter = new MTZWriter(reflectionlist, refinementdata, filename, true);
        }
        mtzwriter.write();
    }

    public void writeSolventMaskCNS(String filename) {
        try {
            PrintWriter cnsfile = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
            cnsfile.println(" ANOMalous=FALSE");
            cnsfile.println(" DECLare NAME=FS DOMAin=RECIprocal TYPE=COMP END");
            for (HKL ih : reflectionlist.hkllist) {
                int i = ih.index();
                cnsfile.printf(" INDE %d %d %d FS= %.4f %.4f\n",
                        ih.h(), ih.k(), ih.l(),
                        refinementdata.fs_f(i),
                        Math.toDegrees(refinementdata.fs_phi(i)));
            }
            cnsfile.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    public void writeSolventMask(String filename) {
        if (solventmodel != SolventModel.NONE) {
            CCP4MapWriter mapwriter = new CCP4MapWriter((int) crs_fs.getXDim(),
                    (int) crs_fs.getYDim(), (int) crs_fs.getZDim(),
                    crystal, filename);
            mapwriter.write(crs_fs.solventGrid);
        }
    }
}
