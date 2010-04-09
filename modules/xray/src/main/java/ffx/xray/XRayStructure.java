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

import java.io.File;
import java.util.List;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

/**
 *
 * @author Tim Fenn
 */
public class XRayStructure {

    private static final Logger logger = Logger.getLogger(XRayStructure.class.getName());
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final CrystalReciprocalSpace crs_fc;
    private final CrystalReciprocalSpace crs_fs;
    private ScaleBulkMinimize scalebulkminimize;
    private SigmaAMinimize sigmaaminimize;
    private SplineMinimize splineminimize;
    private CrystalStats crystalstats;
    private boolean scaled = false;

    public XRayStructure(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(assembly, properties, SolventModel.POLYNOMIAL);
    }

    public XRayStructure(MolecularAssembly assembly,
            CompositeConfiguration properties, int solventmodel) {
        // String name = assembly.getName();
        String name = assembly.getFile().getPath();
        name = name.substring(0, name.lastIndexOf(".pdb"));

        // load the structure
        File tmp = new File(name + ".mtz");
        File mtzfile = null, ciffile = null;
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
                    logger.severe("no input data found!");
                }
            }
        }

        // read in Fo/sigFo/FreeR
        MTZFilter mtzfilter = new MTZFilter();
        CIFFilter ciffilter = new CIFFilter();
        Crystal crystalinit = Crystal.checkProperties(properties);
        Resolution resolutioninit = Resolution.checkProperties(properties);
        if (crystalinit == null || resolutioninit == null) {
            if (mtzfile != null) {
                reflectionlist = mtzfilter.getReflectionList(mtzfile);
            } else {
                reflectionlist = ciffilter.getReflectionList(ciffile);
            }
            if (reflectionlist == null) {
                logger.severe("MTZ/CIF file does not contain full crystal information!");
            }
        } else {
            reflectionlist = new ReflectionList(crystalinit, resolutioninit);
        }

        crystal = reflectionlist.crystal;
        resolution = reflectionlist.resolution;
        refinementdata = new RefinementData(properties, reflectionlist);
        if (mtzfile != null) {
            mtzfilter.readFile(mtzfile, reflectionlist, refinementdata);
        } else {
            ciffilter.readFile(ciffile, reflectionlist, refinementdata);
        }

        List<Atom> atomlist = assembly.getAtomList();
        Atom atomarray[] = atomlist.toArray(new Atom[atomlist.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        crs_fc = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, false);
        refinementdata.setCrystalReciprocalSpaceFc(crs_fc);
        crs_fs = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, true, solventmodel);
        refinementdata.setCrystalReciprocalSpaceFs(crs_fs);

        crystalstats = new CrystalStats(reflectionlist, refinementdata);
    }

    public void setSolventBinaryrad(double rad) {
        crs_fs.setSolventBinaryrad(rad);
        refinementdata.solvent_binaryrad = rad;
    }

    public void setSolventA(double a) {
        crs_fs.setSolventA(a);
        refinementdata.solvent_a = a;
    }

    public void setSolventsd(double sd) {
        crs_fs.setSolventsd(sd);
        refinementdata.solvent_sd = sd;
    }

    public void timings() {
        logger.info("performing 10 Fc calculations for timing...");
        for (int i = 0; i < 10; i++) {
            crs_fc.computeDensity(refinementdata.fc);
        }

        logger.info("performing 10 Fs calculations for timing...");
        for (int i = 0; i < 10; i++) {
            crs_fs.computeDensity(refinementdata.fs);
        }
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
        crs_fs.computeDensity(refinementdata.fs);

        // initialize minimizers
        scalebulkminimize = new ScaleBulkMinimize(reflectionlist, refinementdata, crs_fs);
        splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.spline, SplineEnergy.Type.FOFC);

        // minimize
        if (refinementdata.solvent_n > 1
                && refinementdata.gridsearch) {
            scalebulkminimize.minimize(7, 1e-2);
            scalebulkminimize.GridOptimize();
        }
        scalebulkminimize.minimize(7, 1e-4);

        // sigmaA / LLK calculation
        sigmaaminimize = new SigmaAMinimize(reflectionlist, refinementdata);
        sigmaaminimize.minimize(7, 1.0);

        splineminimize.minimize(7, 1e-5);

        scaled = true;
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
        crystalstats.print_scalestats();
        crystalstats.print_hklstats();
        crystalstats.print_snstats();
        crystalstats.print_rstats();
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
}
