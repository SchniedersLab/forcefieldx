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
    final RefinementData refinementdata;
    final CrystalReciprocalSpace crs_fc;
    final CrystalReciprocalSpace crs_fs;
    List<Atom> atomlist;
    ScaleBulkMinimize scalebulkminimize;
    SigmaAMinimize sigmaaminimize;
    SplineMinimize splineminimize;
    CrystalStats crystalstats;
    private boolean scaled = false;

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
                reflectionlist = mtzfilter.getReflectionList(mtzfile, properties);
            } else {
                reflectionlist = ciffilter.getReflectionList(ciffile, properties);
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

        atomlist = assembly[0].getAtomList();
        for (int i = 1; i < assembly.length; i++) {
            List<Atom> tmplist = assembly[i].getAtomList();
            int tmpsize = tmplist.size();
            for (int j = 0; j < tmpsize; j++) {
                Atom tmpatom = tmplist.get(j);
                Character c = tmpatom.getAltLoc();
                if (c == null) {
                    continue;
                }
                if (!c.equals(' ')
                        && !c.equals('A')) {
                    atomlist.add(tmpatom);
                }
            }
        }
        Atom atomarray[] = atomlist.toArray(new Atom[atomlist.size()]);

        // initialize atomic form factors
        for (int i = 0; i < atomarray.length; i++) {
            FormFactor atomff =
                    new FormFactor(atomarray[i], refinementdata.use_3g, 2.0);

            double arad = 2.4;
            double xyz[] = new double[3];
            xyz[0] = atomarray[i].getX() + arad;
            xyz[1] = atomarray[i].getY();
            xyz[2] = atomarray[i].getZ();
            while (true) {
                double rho = atomff.rho(xyz);
                if (rho > 0.1) {
                    arad += 0.5;
                } else if (rho > 0.001) {
                    arad += 0.1;
                } else {
                    arad += 0.2;
                    atomarray[i].setFormFactorWidth(arad);
                    break;
                }
                xyz[0] = atomarray[i].getX() + arad;
            }
        }

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        crs_fc = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, false);
        refinementdata.setCrystalReciprocalSpaceFc(crs_fc);
        if (refinementdata.bulksolvent) {
            crs_fs = new CrystalReciprocalSpace(reflectionlist, atomarray,
                    parallelTeam, parallelTeam, true, solventmodel);
            refinementdata.setCrystalReciprocalSpaceFs(crs_fs);
        } else {
            crs_fs = null;
        }
        crs_fc.setUse3G(refinementdata.use_3g);

        crystalstats = new CrystalStats(reflectionlist, refinementdata);
    }

    public void setSolventBinaryrad(double rad) {
        if (refinementdata.bulksolvent) {
            crs_fs.setSolventBinaryrad(rad);
            refinementdata.solvent_binaryrad = rad;
        }
    }

    public void setSolventA(double a) {
        if (refinementdata.bulksolvent) {
            crs_fs.setSolventA(a);
            refinementdata.solvent_a = a;
        }
    }

    public void setSolventsd(double sd) {
        if (refinementdata.bulksolvent) {
            crs_fs.setSolventsd(sd);
            refinementdata.solvent_sd = sd;
        }
    }

    public void timings() {
        logger.info("performing 10 Fc calculations for timing...");
        for (int i = 0; i < 10; i++) {
            crs_fc.computeDensity(refinementdata.fc);
        }

        if (refinementdata.bulksolvent) {
            logger.info("performing 10 Fs calculations for timing...");
            for (int i = 0; i < 10; i++) {
                crs_fs.computeDensity(refinementdata.fs);
            }
        }
    }

    public void writeSolventMask(String filename) {
        // CNSMapWriter mapwriter = new CNSMapWriter((int) crs_fs.getXDim(),
        if (refinementdata.bulksolvent) {
            CCP4MapWriter mapwriter = new CCP4MapWriter((int) crs_fs.getXDim(),
                    (int) crs_fs.getYDim(), (int) crs_fs.getZDim(),
                    crystal, filename);
            mapwriter.write(crs_fs.df_map);
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
        if (refinementdata.bulksolvent) {
            crs_fs.computeDensity(refinementdata.fs);
        }

        // initialize minimizers
        scalebulkminimize = new ScaleBulkMinimize(reflectionlist, refinementdata, crs_fs);
        splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.spline, SplineEnergy.Type.FOFC);

        // minimize
        if (refinementdata.bulksolvent
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
