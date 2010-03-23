/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;

/**
 *
 * @author fenn
 *
 * CIF file reader
 */
public class CIFFilter {

    private static final Logger logger = Logger.getLogger(MTZFilter.class.getName());
    private double cell[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    private double reshigh = -1.0;
    private int sgnum = -1;
    private int h = -1, k = -1, l = -1, fo = -1, sigfo = -1, rfree = -1;
    private int nall, nobs;

    private static enum Header {

        d_resolution_high, number_all, number_obs, length_a, length_b, length_c,
        angle_alpha, angle_beta, angle_gamma, Int_Tables_number, crystal_id,
        wavelength_id, scale_group_code, status, index_h, index_k, index_l,
        F_meas_au, F_meas_sigma_au, NOVALUE;

        public static Header toHeader(String str) {
            try {
                return valueOf(str);
            } catch (Exception ex) {
                return NOVALUE;
            }
        }
    }

    // null constructor
    public CIFFilter() {
    }

    public ReflectionList getReflectionList(File cifFile) {
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
                    }
                }
            }

            br.close();
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return null;
        }

        if (sgnum < 0 || reshigh < 0 || cell[0] < 0) {
            logger.info("insufficient information in CIF header to generate Reflection List");
            return null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append(String.format("\nOpening %s\n", cifFile.getName()));
            sb.append(String.format("setting up Reflection List based on CIF:\n"));
            sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                    sgnum, SpaceGroup.spaceGroupNames[sgnum - 1]));
            sb.append(String.format("  resolution: %8.3f\n", 0.9999 * reshigh));
            sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]));
            sb.append(String.format("\n  CIF # HKL (observed): %d\n", nobs));
            sb.append(String.format("  CIF # HKL (all):      %d\n", nall));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(cell[0], cell[1], cell[2],
                cell[3], cell[4], cell[5], SpaceGroup.spaceGroupNames[sgnum - 1]);
        Resolution resolution = new Resolution(0.9999 * reshigh);

        ReflectionList reflectionlist = new ReflectionList(crystal, resolution);
        return reflectionlist;
    }

    public boolean readFile(File cifFile, ReflectionList reflectionlist,
            RefinementData refinementdata) {
        int nread, nnan, nres, nignore;

        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));

            String str;
            int ncol = 0;
            boolean inhkl = false;
            int hklline;
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
                        case F_meas_au:
                            fo = ncol;
                            break;
                        case F_meas_sigma_au:
                            sigfo = ncol;
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

            // go back to where the reflections start
            br.reset();

            // read in data
            nread = nnan = nres = nignore = 0;
            while ((str = br.readLine()) != null) {
                // reached end, break
                if (str.startsWith("#END")) {
                    break;
                }

                String strarray[] = str.split("\\s+");

                int ih = Integer.parseInt(strarray[h]);
                int ik = Integer.parseInt(strarray[k]);
                int il = Integer.parseInt(strarray[l]);
                HKL hkl = reflectionlist.getHKL(ih, ik, il);
                if (hkl != null) {
                    boolean isnull = false;

                    if (rfree > 0) {
                        if (strarray[rfree].charAt(0) == 'o') {
                            refinementdata.freer(hkl.index(), 0);
                        } else if (strarray[rfree].charAt(0) == 'f') {
                            refinementdata.freer(hkl.index(), 1);
                        } else if (strarray[rfree].charAt(0) == 'x') {
                            isnull = true;
                            nnan++;
                        } else {
                            refinementdata.freer(hkl.index(),
                                    Integer.parseInt(strarray[rfree]));
                        }
                    }

                    if (fo > 0 && sigfo > 0 && !isnull) {
                        refinementdata.fsigf(hkl.index(),
                                Double.parseDouble(strarray[fo]),
                                Double.parseDouble(strarray[sigfo]));
                    } else {
                        refinementdata.fsigf(hkl.index(),
                                Double.NaN,
                                Double.NaN);
                    }

                    nread++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (Crystal.invressq(reflectionlist.crystal, tmp)
                            > reflectionlist.resolution.invressq_limit()) {
                        nres++;
                    } else {
                        nignore++;
                    }
                }
            }

            br.close();
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return false;
        }

        StringBuffer sb = new StringBuffer();
        sb.append(String.format("\nOpening %s\n", cifFile.getName()));
        sb.append(String.format("# HKL read in:                             %d\n",
                nread));
        sb.append(String.format("# HKL with NaN (ignored):                  %d\n",
                nnan));
        sb.append(String.format("# HKL NOT read in (too high resolution):   %d\n",
                nres));
        sb.append(String.format("# HKL NOT read in (not in internal list?): %d\n",
                nignore));
        sb.append(String.format("# HKL in internal list:                    %d\n",
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
