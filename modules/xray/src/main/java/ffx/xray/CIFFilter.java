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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;

/**
 * CIF file reader
 *
 * @author Timothy D. Fenn
 *
 */
public class CIFFilter implements DiffractionFileFilter {

    private static final Logger logger = Logger.getLogger(CIFFilter.class.getName());
    private double cell[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    private double reshigh = -1.0;
    private String sgname = null;
    private int sgnum = -1;
    private int h = -1, k = -1, l = -1, fo = -1, sigfo = -1, io = -1, sigio = -1, rfree = -1;
    private int nall, nobs;

    private static enum Header {

        d_resolution_high, number_all, number_obs, length_a, length_b, length_c,
        angle_alpha, angle_beta, angle_gamma, Int_Tables_number,
        space_group_name_H_M, crystal_id, wavelength_id, scale_group_code,
        status, index_h, index_k, index_l, F_meas, F_meas_au, F_meas_sigma,
        F_meas_sigma_au, intensity_meas, intensity_sigma, NOVALUE;

        public static Header toHeader(String str) {
            try {
                return valueOf(str.replace('-', '_'));
            } catch (Exception ex) {
                return NOVALUE;
            }
        }
    }

    /**
     * null constructor
     */
    public CIFFilter() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ReflectionList getReflectionList(File cifFile) {
        return getReflectionList(cifFile, null);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ReflectionList getReflectionList(File cifFile, CompositeConfiguration properties) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));

            String str;
            while ((str = br.readLine()) != null) {
                // reached reflections, break
                if (str.startsWith("_refln.")) {
                    break;
                }

                String strarray[] = str.split("\\s+");

                if (strarray[0].startsWith("_reflns")
                        || strarray[0].startsWith("_cell")
                        || strarray[0].startsWith("_symmetry")) {
                    String cifarray[] = strarray[0].split("\\.+");

                    switch (Header.toHeader(cifarray[1])) {
                        case d_resolution_high:
                            reshigh = Double.parseDouble(strarray[1]);
                            break;
                        case number_all:
                            nall = Integer.parseInt(strarray[1]);
                            break;
                        case number_obs:
                            nobs = Integer.parseInt(strarray[1]);
                            break;
                        case length_a:
                            cell[0] = Double.parseDouble(strarray[1]);
                            break;
                        case length_b:
                            cell[1] = Double.parseDouble(strarray[1]);
                            break;
                        case length_c:
                            cell[2] = Double.parseDouble(strarray[1]);
                            break;
                        case angle_alpha:
                            cell[3] = Double.parseDouble(strarray[1]);
                            break;
                        case angle_beta:
                            cell[4] = Double.parseDouble(strarray[1]);
                            break;
                        case angle_gamma:
                            cell[5] = Double.parseDouble(strarray[1]);
                            break;
                        case Int_Tables_number:
                            sgnum = Integer.parseInt(strarray[1]);
                            break;
                        case space_group_name_H_M:
                            String sgnamearray[] = str.split("'+");
                            if (sgnamearray.length > 1) {
                                sgname = sgnamearray[1];
                            } else if (strarray.length > 1) {
                                sgname = strarray[1];
                            }
                            break;
                    }
                }
            }

            br.close();
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return null;
        }

        if (sgnum < 0 && sgname != null) {
            sgnum = SpaceGroup.spaceGroupNumber(SpaceGroup.pdb2ShortName(sgname));
        }

        if (sgnum < 0 || reshigh < 0 || cell[0] < 0) {
            logger.info(" The CIF header contains insufficient information to generate the reflection list.");
            return null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\nOpening %s\n", cifFile.getName()));
            sb.append(String.format("setting up Reflection List based on CIF:\n"));
            sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                    sgnum, SpaceGroup.spaceGroupNames[sgnum - 1]));
            sb.append(String.format("  resolution: %8.3f\n", 0.999999 * reshigh));
            sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]));
            sb.append(String.format("\n  CIF # HKL (observed): %d\n", nobs));
            sb.append(String.format("  CIF # HKL (all):      %d\n", nall));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(cell[0], cell[1], cell[2],
                cell[3], cell[4], cell[5], SpaceGroup.spaceGroupNames[sgnum - 1]);
        double sampling = 1.0 / 1.5;
        if (properties != null) {
            sampling = properties.getDouble("sampling", 1.0 / 1.5);
        }
        Resolution resolution = new Resolution(0.999999 * reshigh, sampling);

        ReflectionList reflectionlist = new ReflectionList(crystal, resolution,
                properties);
        return reflectionlist;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getResolution(File cifFile, Crystal crystal) {
        double res = Double.POSITIVE_INFINITY;

        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));

            String str;
            int ncol = 0;
            boolean inhkl = false;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                if (strarray[0].startsWith("_refln.")) {
                    inhkl = true;
                    br.mark(0);
                    String cifarray[] = strarray[0].split("\\.+");
                    switch (Header.toHeader(cifarray[1])) {
                        case index_h:
                            h = ncol;
                            break;
                        case index_k:
                            k = ncol;
                            break;
                        case index_l:
                            l = ncol;
                            break;
                    }

                    ncol++;
                } else if (inhkl) {
                    if (h < 0 || k < 0 || l < 0) {
                        String message = "Fatal error in CIF file - no H K L indexes?\n";
                        logger.log(Level.SEVERE, message);
                        return -1.0;
                    }
                    break;
                }
            }

            // go back to where the reflections start
            br.reset();
            HKL hkl = new HKL();
            while ((str = br.readLine()) != null) {
                // reached end, break
                if (str.trim().startsWith("#END")) {
                    break;
                } else if (str.trim().startsWith("data")) {
                    break;
                } else if (str.trim().startsWith("#")) {
                    continue;
                }

                String strarray[] = str.trim().split("\\s+");
                // some files split data on to multiple lines
                while (strarray.length < ncol) {
                    str = str + " " + br.readLine();
                    strarray = str.trim().split("\\s+");
                }

                int ih = Integer.parseInt(strarray[h]);
                int ik = Integer.parseInt(strarray[k]);
                int il = Integer.parseInt(strarray[l]);

                hkl.h(ih);
                hkl.k(ik);
                hkl.l(il);
                res = Math.min(res, Crystal.res(crystal, hkl));
            }
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return -1.0;
        }

        return res;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readFile(File cifFile, ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, CompositeConfiguration properties) {
        int nread, nnan, nres, nignore, ncifignore, nfriedel, ncut;
        boolean transpose = false;
        boolean intensitiesToAmplitudes = false;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Opening %s\n", cifFile.getName()));
        if (refinementdata.rfreeflag < 0) {
            refinementdata.setFreeRFlag(1);
            sb.append(String.format(" Setting R free flag to CIF default: %d\n", refinementdata.rfreeflag));
        }

        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));

            String str;
            int ncol = 0;
            boolean inhkl = false;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                if (strarray[0].startsWith("_refln.")) {
                    inhkl = true;
                    br.mark(0);
                    String cifarray[] = strarray[0].split("\\.+");
                    switch (Header.toHeader(cifarray[1])) {
                        case index_h:
                            h = ncol;
                            break;
                        case index_k:
                            k = ncol;
                            break;
                        case index_l:
                            l = ncol;
                            break;
                        case F_meas:
                        case F_meas_au:
                            fo = ncol;
                            break;
                        case F_meas_sigma:
                        case F_meas_sigma_au:
                            sigfo = ncol;
                            break;
                        case intensity_meas:
                            io = ncol;
                            break;
                        case intensity_sigma:
                            sigio = ncol;
                            break;
                        case status:
                            rfree = ncol;
                            break;
                    }

                    ncol++;
                } else if (inhkl) {
                    if (h < 0 || k < 0 || l < 0) {
                        String message = "Fatal error in CIF file - no H K L indexes?\n";
                        logger.log(Level.SEVERE, message);
                        return false;
                    }
                    break;
                }
            }

            if (fo < 0 && sigfo < 0 && io > 0 && sigio > 0) {
                intensitiesToAmplitudes = true;
            }

            if (fo < 0 && io < 0) {
                logger.severe("Reflection data (I/F) not found in CIF file!");
            }

            // go back to where the reflections start
            br.reset();

            // check if HKLs need to be transposed or not
            HKL mate = new HKL();
            int nposignore = 0;
            int ntransignore = 0;
            while ((str = br.readLine()) != null) {
                // reached end, break
                if (str.trim().startsWith("#END")) {
                    break;
                } else if (str.trim().startsWith("data")) {
                    break;
                } else if (str.trim().startsWith("#")) {
                    continue;
                }

                String strarray[] = str.trim().split("\\s+");
                // some files split data on to multiple lines
                while (strarray.length < ncol) {
                    str = str + " " + br.readLine();
                    strarray = str.trim().split("\\s+");
                }

                if (rfree > 0) {
                    // ignored cases
                    if (strarray[rfree].charAt(0) == 'x'
                            || strarray[rfree].charAt(0) == '<'
                            || strarray[rfree].charAt(0) == '-'
                            || strarray[rfree].charAt(0) == 'h'
                            || strarray[rfree].charAt(0) == 'l') {
                        continue;
                    }
                }
                int ih = Integer.parseInt(strarray[h]);
                int ik = Integer.parseInt(strarray[k]);
                int il = Integer.parseInt(strarray[l]);
                boolean friedel = reflectionlist.findSymHKL(ih, ik, il, mate, false);
                HKL hklpos = reflectionlist.getHKL(mate);
                if (hklpos == null) {
                    nposignore++;
                }

                friedel = reflectionlist.findSymHKL(ih, ik, il, mate, true);
                HKL hkltrans = reflectionlist.getHKL(mate);
                if (hkltrans == null) {
                    ntransignore++;
                }
            }
            if (nposignore > ntransignore) {
                transpose = true;
            }

            // reopen to start at beginning
            br = new BufferedReader(new FileReader(cifFile));
            inhkl = false;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                if (strarray[0].startsWith("_refln.")) {
                    br.mark(0);
                    inhkl = true;
                } else if (inhkl) {
                    break;
                }
            }

            // go back to where the reflections start
            br.reset();

            // read in data
            double anofsigf[][] = new double[refinementdata.n][4];
            for (int i = 0; i < refinementdata.n; i++) {
                anofsigf[i][0] = anofsigf[i][1] = anofsigf[i][2] = anofsigf[i][3] = Double.NaN;
            }
            nread = nnan = nres = nignore = ncifignore = nfriedel = ncut = 0;
            while ((str = br.readLine()) != null) {
                // reached end, break
                if (str.trim().startsWith("#END")) {
                    break;
                } else if (str.trim().startsWith("data")) {
                    break;
                } else if (str.trim().startsWith("#")) {
                    continue;
                }

                String strarray[] = str.trim().split("\\s+");
                // some files split data on to multiple lines
                while (strarray.length < ncol) {
                    str = str + " " + br.readLine();
                    strarray = str.trim().split("\\s+");
                }

                int ih = Integer.parseInt(strarray[h]);
                int ik = Integer.parseInt(strarray[k]);
                int il = Integer.parseInt(strarray[l]);
                boolean friedel = reflectionlist.findSymHKL(ih, ik, il, mate, transpose);
                HKL hkl = reflectionlist.getHKL(mate);

                if (hkl != null) {
                    boolean isnull = false;

                    if (rfree > 0) {
                        if (strarray[rfree].charAt(0) == 'o') {
                            refinementdata.setFreeR(hkl.index(), 0);
                        } else if (strarray[rfree].charAt(0) == 'f') {
                            refinementdata.setFreeR(hkl.index(), 1);
                        } else if (strarray[rfree].charAt(0) == 'x') {
                            isnull = true;
                            nnan++;
                        } else if (strarray[rfree].charAt(0) == '<'
                                || strarray[rfree].charAt(0) == '-'
                                || strarray[rfree].charAt(0) == 'h'
                                || strarray[rfree].charAt(0) == 'l') {
                            isnull = true;
                            ncifignore++;
                        } else {
                            refinementdata.setFreeR(hkl.index(),
                                    Integer.parseInt(strarray[rfree]));
                        }
                    }

                    if (!intensitiesToAmplitudes && !isnull) {
                        if (strarray[fo].charAt(0) == '?'
                                || strarray[sigfo].charAt(0) == '?') {
                            isnull = true;
                            nnan++;
                            continue;
                        }

                        if (refinementdata.fsigfcutoff > 0.0) {
                            double f1 = Double.parseDouble(strarray[fo]);
                            double sigf1 = Double.parseDouble(strarray[sigfo]);
                            if ((f1 / sigf1) < refinementdata.fsigfcutoff) {
                                ncut++;
                                continue;
                            }
                        }

                        if (friedel) {
                            anofsigf[hkl.index()][2] = Double.parseDouble(strarray[fo]);
                            anofsigf[hkl.index()][3] = Double.parseDouble(strarray[sigfo]);
                            nfriedel++;
                        } else {
                            anofsigf[hkl.index()][0] = Double.parseDouble(strarray[fo]);
                            anofsigf[hkl.index()][1] = Double.parseDouble(strarray[sigfo]);
                        }
                    }

                    if (intensitiesToAmplitudes && !isnull) {
                        if (strarray[io].charAt(0) == '?'
                                || strarray[sigio].charAt(0) == '?') {
                            isnull = true;
                            nnan++;
                            continue;
                        }

                        if (friedel) {
                            anofsigf[hkl.index()][2] = Double.parseDouble(strarray[io]);
                            anofsigf[hkl.index()][3] = Double.parseDouble(strarray[sigio]);
                            nfriedel++;
                        } else {
                            anofsigf[hkl.index()][0] = Double.parseDouble(strarray[io]);
                            anofsigf[hkl.index()][1] = Double.parseDouble(strarray[sigio]);
                        }
                    }

                    nread++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (!reflectionlist.resolution.inInverseResSqRange(Crystal.invressq(reflectionlist.crystal, tmp))) {
                        nres++;
                    } else {
                        nignore++;
                    }
                }
            }

            br.close();

            // set up fsigf from F+ and F-
            refinementdata.generate_fsigf_from_anofsigf(anofsigf);

            if (intensitiesToAmplitudes) {
                refinementdata.intensities_to_amplitudes();
            }
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return false;
        }

        sb.append(String.format(" HKL data is %s\n",
                transpose ? "transposed" : "not transposed"));
        sb.append(String.format(" HKL read in:                             %d\n",
                nread));
        sb.append(String.format(" HKL read as friedel mates:               %d\n",
                nfriedel));
        sb.append(String.format(" HKL with NaN (ignored):                  %d\n",
                nnan));
        sb.append(String.format(" HKL NOT read in (status <, -, h or l):   %d\n",
                ncifignore));
        sb.append(String.format(" HKL NOT read in (too high resolution):   %d\n",
                nres));
        sb.append(String.format(" HKL NOT read in (not in internal list?): %d\n",
                nignore));
        sb.append(String.format(" HKL NOT read in (F/sigF cutoff):         %d\n",
                ncut));
        sb.append(String.format(" HKL in internal list:                    %d\n",
                reflectionlist.hkllist.size()));
        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }

        if (rfree < 0) {
            refinementdata.generateRFree();
        }

        return true;
    }
}
