/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.xray;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 *
 * @author fenn
 */
public class CNSFilter {

    private static final Logger logger = Logger.getLogger(CNSFilter.class.getName());
    private double cell[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    private double reshigh = -1.0;
    private String sgname = null;
    private int sgnum = -1;

    // null constructor
    public CNSFilter() {
    }

    public ReflectionList getReflectionList(File cnsFile) {
        return getReflectionList(cnsFile, null);
    }

    public ReflectionList getReflectionList(File cnsFile, CompositeConfiguration properties) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(cnsFile));

            String str;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                if (strarray[0].equalsIgnoreCase("{")) {
                    if (strarray[1].toLowerCase().startsWith("sg=")) {
                        sgname = strarray[1].substring(3);
                        cell[0] = Double.parseDouble(strarray[2].substring(2));
                        cell[1] = Double.parseDouble(strarray[3].substring(2));
                        cell[2] = Double.parseDouble(strarray[4].substring(2));
                        cell[3] = Double.parseDouble(strarray[5].substring(6));
                        cell[4] = Double.parseDouble(strarray[6].substring(5));
                        cell[5] = Double.parseDouble(strarray[7].substring(6));
                    }
                } else if (strarray[0].equalsIgnoreCase("CRYST1")) {
                    cell[0] = Double.parseDouble(strarray[1]);
                    cell[1] = Double.parseDouble(strarray[2]);
                    cell[2] = Double.parseDouble(strarray[3]);
                    cell[3] = Double.parseDouble(strarray[4]);
                    cell[4] = Double.parseDouble(strarray[5]);
                    cell[5] = Double.parseDouble(strarray[6]);
                    sgname = SpaceGroup.pdb2ShortName(str.substring(55, 65));
                } else if (strarray[0].toLowerCase().startsWith("inde")) {
                    break;
                }
            }

            br.close();
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return null;
        }

        Resolution resolution = null;
        if (properties != null) {
            resolution = Resolution.checkProperties(properties);
            reshigh = resolution.resolution;
        }

        if (sgname != null) {
            sgnum = SpaceGroup.spaceGroupNumber(sgname);
        }

        if (sgnum < 0 || cell[0] < 0 || resolution == null) {
            logger.info("insufficient information to generate Reflection List");
            return null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\nOpening %s\n", cnsFile.getName()));
            sb.append(String.format("setting up Reflection List based on CNS:\n"));
            sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                    sgnum, SpaceGroup.spaceGroupNames[sgnum - 1]));
            sb.append(String.format("  resolution: %8.3f\n", 0.9999 * reshigh));
            sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(cell[0], cell[1], cell[2],
                cell[3], cell[4], cell[5], SpaceGroup.spaceGroupNames[sgnum - 1]);

        ReflectionList reflectionlist = new ReflectionList(crystal, resolution,
                properties);
        return reflectionlist;
    }

    public boolean readFile(File cnsFile, ReflectionList reflectionlist,
            RefinementData refinementdata) {
        int nread, nres, nignore, nfriedel, ncut;
        boolean transpose = false;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\nOpening %s\n", cnsFile.getName()));
        if (refinementdata.rfreeflag < 0) {
            refinementdata.set_freerflag(1);
            sb.append(String.format("Setting R free flag to CNS default: %d\n", refinementdata.rfreeflag));
        }

        try {
            BufferedReader br = new BufferedReader(new FileReader(cnsFile));

            String str;
            boolean hashkl, hasfo, hassigfo, hasfree;
            int ih, ik, il, free;
            double fo, sigfo;

            hashkl = hasfo = hassigfo = hasfree = false;
            ih = ik = il = free = -1;
            fo = sigfo = -1.0;

            // check if HKLs need to be transposed or not
            HKL mate = new HKL();
            int nposignore = 0;
            int ntransignore = 0;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                for (int i = 0; i < strarray.length; i++) {
                    if (strarray[i].toLowerCase().startsWith("inde")) {
                        if (i < strarray.length - 3) {
                            ih = Integer.parseInt(strarray[i + 1]);
                            ik = Integer.parseInt(strarray[i + 2]);
                            il = Integer.parseInt(strarray[i + 3]);
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
                    }
                }
            }
            if (nposignore > ntransignore) {
                transpose = true;
            }

            // reopen to start at beginning
            br = new BufferedReader(new FileReader(cnsFile));

            // read in data
            nread = nres = nignore = nfriedel = ncut = 0;
            while ((str = br.readLine()) != null) {
                String strarray[] = str.split("\\s+");

                for (int i = 0; i < strarray.length; i++) {
                    if (strarray[i].toLowerCase().startsWith("inde")) {
                        if (hashkl && hasfo && hassigfo && hasfree) {
                            boolean friedel = reflectionlist.findSymHKL(ih, ik, il, mate, transpose);
                            HKL hkl = reflectionlist.getHKL(mate);
                            if (hkl != null) {
                                if (refinementdata.fsigfcutoff > 0.0
                                        && (fo / sigfo) < refinementdata.fsigfcutoff) {
                                    ncut++;
                                } else if (friedel) {
                                    refinementdata.set_ano_fsigfminus(hkl.index(), fo, sigfo);
                                    nfriedel++;
                                } else if (!friedel) {
                                    refinementdata.set_ano_fsigfplus(hkl.index(), fo, sigfo);
                                }
                                refinementdata.set_freer(hkl.index(), free);
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
                        hashkl = false;
                        hasfo = false;
                        hassigfo = false;
                        hasfree = false;
                        if (i < strarray.length - 3) {
                            ih = Integer.parseInt(strarray[i + 1]);
                            ik = Integer.parseInt(strarray[i + 2]);
                            il = Integer.parseInt(strarray[i + 3]);
                            hashkl = true;
                        }
                    }
                    if (strarray[i].toLowerCase().startsWith("fobs=")) {
                        fo = Double.parseDouble(strarray[i + 1]);
                        hasfo = true;
                    }
                    if (strarray[i].toLowerCase().startsWith("sigma=")) {
                        sigfo = Double.parseDouble(strarray[i + 1]);
                        hassigfo = true;
                    }
                    if (strarray[i].toLowerCase().startsWith("test=")) {
                        free = Integer.parseInt(strarray[i + 1]);
                        hasfree = true;
                    }
                }
            }

            br.close();
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return false;
        }

        // set up fsigf from F+ and F-
        refinementdata.generate_fsigf_from_anofsigf();

        sb.append(String.format("# HKL read in:                             %d\n",
                nread));
        sb.append(String.format("# HKL read as friedel mates:               %d\n",
                nfriedel));
        sb.append(String.format("# HKL NOT read in (too high resolution):   %d\n",
                nres));
        sb.append(String.format("# HKL NOT read in (not in internal list?): %d\n",
                nignore));
        sb.append(String.format("# HKL NOT read in (F/sigF cutoff):         %d\n",
                ncut));
        sb.append(String.format("# HKL in internal list:                    %d\n",
                reflectionlist.hkllist.size()));
        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }

        return true;
    }
}
