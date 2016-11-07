/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;

import static org.apache.commons.math3.util.FastMath.min;

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
    private final double cell[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    private double resHigh = -1.0;
    private String spacegroupName = null;
    private int spacegroupNum = -1;
    private int h = -1, k = -1, l = -1;
    private int fo = -1, sigFo = -1;
    private int io = -1, sigIo = -1;
    private int rFree = -1;
    private int nAll, nObs;

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
            String string;
            while ((string = br.readLine()) != null) {

                // Reached reflections, break.
                if (string.startsWith("_refln.")) {
                    break;
                }
                String strArray[] = string.split("\\s+");
                if (strArray[0].startsWith("_reflns")
                        || strArray[0].startsWith("_cell")
                        || strArray[0].startsWith("_symmetry")) {
                    String cifArray[] = strArray[0].split("\\.+");
                    switch (Header.toHeader(cifArray[1])) {
                        case d_resolution_high:
                            resHigh = Double.parseDouble(strArray[1]);
                            break;
                        case number_all:
                            nAll = Integer.parseInt(strArray[1]);
                            break;
                        case number_obs:
                            nObs = Integer.parseInt(strArray[1]);
                            break;
                        case length_a:
                            cell[0] = Double.parseDouble(strArray[1]);
                            break;
                        case length_b:
                            cell[1] = Double.parseDouble(strArray[1]);
                            break;
                        case length_c:
                            cell[2] = Double.parseDouble(strArray[1]);
                            break;
                        case angle_alpha:
                            cell[3] = Double.parseDouble(strArray[1]);
                            break;
                        case angle_beta:
                            cell[4] = Double.parseDouble(strArray[1]);
                            break;
                        case angle_gamma:
                            cell[5] = Double.parseDouble(strArray[1]);
                            break;
                        case Int_Tables_number:
                            spacegroupNum = Integer.parseInt(strArray[1]);
                            break;
                        case space_group_name_H_M:
                            String spacegroupNameArray[] = string.split("'+");
                            if (spacegroupNameArray.length > 1) {
                                spacegroupName = spacegroupNameArray[1];
                            } else if (strArray.length > 1) {
                                spacegroupName = strArray[1];
                            }
                            break;
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            String message = " CIF IO Exception.";
            logger.log(Level.WARNING, message, e);
            return null;
        }

        if (spacegroupNum < 0 && spacegroupName != null) {
            spacegroupNum = SpaceGroup.spaceGroupNumber(SpaceGroup.pdb2ShortName(spacegroupName));
        }

        if (spacegroupNum < 0 || resHigh < 0 || cell[0] < 0) {
            logger.info(" The CIF header contains insufficient information to generate the reflection list.");
            return null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\nOpening %s\n", cifFile.getName()));
            sb.append(String.format("setting up Reflection List based on CIF:\n"));
            sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                    spacegroupNum, SpaceGroup.spaceGroupNames[spacegroupNum - 1]));
            sb.append(String.format("  resolution: %8.3f\n", 0.999999 * resHigh));
            sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]));
            sb.append(String.format("\n  CIF # HKL (observed): %d\n", nObs));
            sb.append(String.format("  CIF # HKL (all):      %d\n", nAll));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(cell[0], cell[1], cell[2],
                cell[3], cell[4], cell[5], SpaceGroup.spaceGroupNames[spacegroupNum - 1]);
        double sampling = 1.0 / 1.5;
        if (properties != null) {
            sampling = properties.getDouble("sampling", 1.0 / 1.5);
        }
        Resolution resolution = new Resolution(0.999999 * resHigh, sampling);
        ReflectionList reflectionlist = new ReflectionList(crystal,
                resolution, properties);
        return reflectionlist;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getResolution(File cifFile, Crystal crystal) {

        double resolution = Double.POSITIVE_INFINITY;

        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));
            String string;
            int nCol = 0;
            boolean inHKL = false;
            while ((string = br.readLine()) != null) {
                String strArray[] = string.split("\\s+");
                if (strArray[0].startsWith("_refln.")) {
                    inHKL = true;
                    br.mark(0);
                    String cifArray[] = strArray[0].split("\\.+");
                    switch (Header.toHeader(cifArray[1])) {
                        case index_h:
                            h = nCol;
                            break;
                        case index_k:
                            k = nCol;
                            break;
                        case index_l:
                            l = nCol;
                            break;
                    }
                    nCol++;
                } else if (inHKL) {
                    if (h < 0 || k < 0 || l < 0) {
                        String message = " Fatal error in CIF file - no H K L indexes?\n";
                        logger.log(Level.SEVERE, message);
                        return -1.0;
                    }
                    break;
                }
            }

            // Go back to where the reflections start.
            br.reset();
            HKL hkl = new HKL();
            while ((string = br.readLine()) != null) {

                // Reached end, break.
                if (string.trim().startsWith("#END")) {
                    break;
                } else if (string.trim().startsWith("data")) {
                    break;
                } else if (string.trim().startsWith("#")) {
                    continue;
                }

                // Some files split data on to multiple lines.
                String strArray[] = string.trim().split("\\s+");
                while (strArray.length < nCol) {
                    string = string + " " + br.readLine();
                    strArray = string.trim().split("\\s+");
                }

                int ih = Integer.parseInt(strArray[h]);
                int ik = Integer.parseInt(strArray[k]);
                int il = Integer.parseInt(strArray[l]);

                hkl.h(ih);
                hkl.k(ik);
                hkl.l(il);
                resolution = min(resolution, Crystal.res(crystal, hkl));
            }
        } catch (IOException e) {
            String message = " CIF IO Exception.";
            logger.log(Level.WARNING, message, e);
            return -1.0;
        }
        return resolution;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readFile(File cifFile, ReflectionList reflectionList,
            DiffractionRefinementData refinementData, CompositeConfiguration properties) {

        int nRead, nNAN, nRes;
        int nIgnore, nCIFIgnore;
        int nFriedel, nCut;
        boolean transpose = false;
        boolean intensitiesToAmplitudes = false;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Opening %s\n", cifFile.getName()));
        if (refinementData.rfreeflag < 0) {
            refinementData.setFreeRFlag(1);
            sb.append(format(" Setting R free flag to CIF default: %d\n", refinementData.rfreeflag));
        }

        try {
            BufferedReader br = new BufferedReader(new FileReader(cifFile));

            String string;
            int nCol = 0;
            boolean inHKL = false;
            while ((string = br.readLine()) != null) {
                String stringArray[] = string.split("\\s+");
                if (stringArray[0].startsWith("_refln.")) {
                    inHKL = true;
                    br.mark(0);
                    String cifArray[] = stringArray[0].split("\\.+");
                    switch (Header.toHeader(cifArray[1])) {
                        case index_h:
                            h = nCol;
                            break;
                        case index_k:
                            k = nCol;
                            break;
                        case index_l:
                            l = nCol;
                            break;
                        case F_meas:
                        case F_meas_au:
                            fo = nCol;
                            break;
                        case F_meas_sigma:
                        case F_meas_sigma_au:
                            sigFo = nCol;
                            break;
                        case intensity_meas:
                            io = nCol;
                            break;
                        case intensity_sigma:
                            sigIo = nCol;
                            break;
                        case status:
                            rFree = nCol;
                            break;
                    }
                    nCol++;
                } else if (inHKL) {
                    if (h < 0 || k < 0 || l < 0) {
                        String message = " Fatal error in CIF file - no H K L indexes?\n";
                        logger.log(Level.SEVERE, message);
                        return false;
                    }
                    break;
                }
            }

            if (fo < 0 && sigFo < 0 && io > 0 && sigIo > 0) {
                intensitiesToAmplitudes = true;
            }

            if (fo < 0 && io < 0) {
                logger.severe(" Reflection data (I/F) not found in CIF file!");
            }

            // Go back to where the reflections start.
            br.reset();

            // Check if HKLs need to be transposed or not.
            HKL mate = new HKL();
            int nPosIgnore = 0;
            int nTransIgnore = 0;
            while ((string = br.readLine()) != null) {

                if (string.trim().startsWith("#END")) {
                    // Reached end, break.
                    break;
                } else if (string.trim().startsWith("data")) {
                    break;
                } else if (string.trim().startsWith("#")) {
                    continue;
                }

                // Some files split data on to multiple lines.
                String strArray[] = string.trim().split("\\s+");
                while (strArray.length < nCol) {
                    string = string + " " + br.readLine();
                    strArray = string.trim().split("\\s+");
                }

                if (rFree > 0) {
                    // Ignored cases.
                    if (strArray[rFree].charAt(0) == 'x'
                            || strArray[rFree].charAt(0) == '<'
                            || strArray[rFree].charAt(0) == '-'
                            || strArray[rFree].charAt(0) == 'h'
                            || strArray[rFree].charAt(0) == 'l') {
                        continue;
                    }
                }

                int ih = Integer.parseInt(strArray[h]);
                int ik = Integer.parseInt(strArray[k]);
                int il = Integer.parseInt(strArray[l]);

                boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, false);
                HKL hklPos = reflectionList.getHKL(mate);
                if (hklPos == null) {
                    nPosIgnore++;
                }

                friedel = reflectionList.findSymHKL(ih, ik, il, mate, true);
                HKL hklTrans = reflectionList.getHKL(mate);
                if (hklTrans == null) {
                    nTransIgnore++;
                }
            }
            if (nPosIgnore > nTransIgnore) {
                transpose = true;
            }

            // Re-open to start at beginning.
            br = new BufferedReader(new FileReader(cifFile));
            inHKL = false;
            while ((string = br.readLine()) != null) {
                String strArray[] = string.split("\\s+");
                if (strArray[0].startsWith("_refln.")) {
                    br.mark(0);
                    inHKL = true;
                } else if (inHKL) {
                    break;
                }
            }

            // Go back to where the reflections start.
            br.reset();

            // Read in data.
            double anofSigF[][] = new double[refinementData.n][4];
            for (int i = 0; i < refinementData.n; i++) {
                anofSigF[i][0] = anofSigF[i][1] = anofSigF[i][2] = anofSigF[i][3] = Double.NaN;
            }
            nRead = nNAN = nRes = nIgnore = nCIFIgnore = nFriedel = nCut = 0;
            while ((string = br.readLine()) != null) {

                // Reached end, break.
                if (string.trim().startsWith("#END")) {
                    break;
                } else if (string.trim().startsWith("data")) {
                    break;
                } else if (string.trim().startsWith("#")) {
                    continue;
                }

                // Some files split data on to multiple lines.
                String strArray[] = string.trim().split("\\s+");
                while (strArray.length < nCol) {
                    string = string + " " + br.readLine();
                    strArray = string.trim().split("\\s+");
                }

                int ih = Integer.parseInt(strArray[h]);
                int ik = Integer.parseInt(strArray[k]);
                int il = Integer.parseInt(strArray[l]);

                boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, transpose);
                HKL hkl = reflectionList.getHKL(mate);
                if (hkl != null) {
                    boolean isnull = false;
                    if (rFree > 0) {
                        if (strArray[rFree].charAt(0) == 'o') {
                            refinementData.setFreeR(hkl.index(), 0);
                        } else if (strArray[rFree].charAt(0) == 'f') {
                            refinementData.setFreeR(hkl.index(), 1);
                        } else if (strArray[rFree].charAt(0) == 'x') {
                            isnull = true;
                            nNAN++;
                        } else if (strArray[rFree].charAt(0) == '<'
                                || strArray[rFree].charAt(0) == '-'
                                || strArray[rFree].charAt(0) == 'h'
                                || strArray[rFree].charAt(0) == 'l') {
                            isnull = true;
                            nCIFIgnore++;
                        } else {
                            refinementData.setFreeR(hkl.index(),
                                    Integer.parseInt(strArray[rFree]));
                        }
                    }

                    if (!intensitiesToAmplitudes && !isnull) {
                        if (strArray[fo].charAt(0) == '?'
                                || strArray[sigFo].charAt(0) == '?') {
                            isnull = true;
                            nNAN++;
                            continue;
                        }

                        if (refinementData.fsigfcutoff > 0.0) {
                            double f1 = Double.parseDouble(strArray[fo]);
                            double sigf1 = Double.parseDouble(strArray[sigFo]);
                            if ((f1 / sigf1) < refinementData.fsigfcutoff) {
                                nCut++;
                                continue;
                            }
                        }

                        if (friedel) {
                            anofSigF[hkl.index()][2] = Double.parseDouble(strArray[fo]);
                            anofSigF[hkl.index()][3] = Double.parseDouble(strArray[sigFo]);
                            nFriedel++;
                        } else {
                            anofSigF[hkl.index()][0] = Double.parseDouble(strArray[fo]);
                            anofSigF[hkl.index()][1] = Double.parseDouble(strArray[sigFo]);
                        }
                    }

                    if (intensitiesToAmplitudes && !isnull) {
                        if (strArray[io].charAt(0) == '?'
                                || strArray[sigIo].charAt(0) == '?') {
                            isnull = true;
                            nNAN++;
                            continue;
                        }

                        if (friedel) {
                            anofSigF[hkl.index()][2] = Double.parseDouble(strArray[io]);
                            anofSigF[hkl.index()][3] = Double.parseDouble(strArray[sigIo]);
                            nFriedel++;
                        } else {
                            anofSigF[hkl.index()][0] = Double.parseDouble(strArray[io]);
                            anofSigF[hkl.index()][1] = Double.parseDouble(strArray[sigIo]);
                        }
                    }

                    nRead++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (!reflectionList.resolution.inInverseResSqRange(
                            Crystal.invressq(reflectionList.crystal, tmp))) {
                        nRes++;
                    } else {
                        nIgnore++;
                    }
                }
            }
            br.close();

            // Set up fsigf from F+ and F-.
            refinementData.generate_fsigf_from_anofsigf(anofSigF);
            if (intensitiesToAmplitudes) {
                refinementData.intensities_to_amplitudes();
            }
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return false;
        }

        sb.append(String.format(" HKL data is %s\n",
                transpose ? "transposed" : "not transposed"));
        sb.append(String.format(" HKL read in:                             %d\n",
                nRead));
        sb.append(String.format(" HKL read as friedel mates:               %d\n",
                nFriedel));
        sb.append(String.format(" HKL with NaN (ignored):                  %d\n",
                nNAN));
        sb.append(String.format(" HKL NOT read in (status <, -, h or l):   %d\n",
                nCIFIgnore));
        sb.append(String.format(" HKL NOT read in (too high resolution):   %d\n",
                nRes));
        sb.append(String.format(" HKL NOT read in (not in internal list?): %d\n",
                nIgnore));
        sb.append(String.format(" HKL NOT read in (F/sigF cutoff):         %d\n",
                nCut));
        sb.append(String.format(" HKL in internal list:                    %d\n",
                reflectionList.hkllist.size()));

        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }

        if (rFree < 0) {
            refinementData.generateRFree();
        }
        return true;
    }
}
