//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.crystal;

import java.util.HashMap;

import static java.lang.Integer.parseInt;

/**
 * Information about the 230 space groups.
 */
public class SpaceGroupInfo {

  private SpaceGroupInfo() {
    // Prevent instantiation.
  }

  /**
   * Names of the 230 three-dimensional space groups.
   *
   * @since 1.0
   */
  public static final String[] spaceGroupNames = {
      "P1", "P-1", "P2", "P21", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m", "P21/m", "C2/m", "P2/c",
      "P21/c", "C2/c", "P222", "P2221", "P21212", "P212121", "C2221", "C222", "F222", "I222",
      "I212121", "Pmm2", "Pmc21", "Pcc2", "Pma2", "Pca21", "Pnc2", "Pmn21", "Pba2", "Pna21", "Pnn2",
      "Cmm2", "Cmc21", "Ccc2", "Amm2", "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2",
      "Ima2", "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", "Pbcm",
      "Pnnm", "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm", "Cccm", "Cmma", "Ccca", "Fmmm",
      "Fddd", "Immm", "Ibam", "Ibca", "Imma", "P4", "P41", "P42", "P43", "I4", "I41", "P-4", "I-4",
      "P4/m", "P42/m", "P4/n", "P42/n", "I4/m", "I41/a", "P422", "P4212", "P4122", "P41212", "P4222",
      "P42212", "P4322", "P43212", "I422", "I4122", "P4mm", "P4bm", "P42cm", "P42nm", "P4cc", "P4nc",
      "P42mc", "P42bc", "I4mm", "I4cm", "I41md", "I41cd", "P-42m", "P-42c", "P-421m", "P-421c",
      "P-4m2", "P-4c2", "P-4b2", "P-4n2", "I-4m2", "I-4c2", "I-42m", "I-42d", "P4/mmm", "P4/mcc",
      "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", "P42/mmc", "P42/mcm", "P42/nbc",
      "P42/nnm", "P42/mbc", "P42/mnm", "P42/nmc", "P42/ncm", "I4/mmm", "I4/mcm", "I41/amd",
      "I41/acd", "P3", "P31", "P32", "H3", "P-3", "H-3", "P312", "P321", "P3112", "P3121", "P3212",
      "P3221", "H32", "P3m1", "P31m", "P3c1", "P31c", "H3m", "H3c", "P-31m", "P-31c", "P-3m1",
      "P-3c1", "H-3m", "H-3c", "P6", "P61", "P65", "P62", "P64", "P63", "P-6", "P6/m", "P63/m",
      "P622", "P6122", "P6522", "P6222", "P6422", "P6322", "P6mm", "P6cc", "P63cm", "P63mc", "P-6m2",
      "P-6c2", "P-62m", "P-62c", "P6/mmm", "P6/mcc", "P63/mcm", "P63/mmc", "P23", "F23", "I23",
      "P213", "I213", "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3", "P432", "P4232",
      "F432", "F4132", "I432", "P4332", "P4132", "I4132", "P-43m", "F-43m", "I-43m", "P-43n",
      "F-43c", "I-43d", "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c", "Fd-3m", "Fd-3c",
      "Im-3m", "Ia-3d"
  };

  /**
   * Space group frequency ranking for the 1,231,510 CSD structures for which the space group is fully
   * defined. Statistics for enantiomorphous space groups are as reported in the CSD.
   * <p>
   * 965,978 (78%) of structures adopt centrosymmetric space groups.
   * <p>
   * 265,597 (22%) adopt non-centrosymmetric space groups.
   * <p>
   * 199,748 (16%) structures adopt Sohncke space groups.
   * <p>
   * As of January 2023.
   */
  public static final double[] csdPercent = {
      0.9904, 25.1399, 0.0178, 5.1417, 0.8494, 0.0039, 0.4518, 0.0352, 1.0386, 0.0151, 0.4617,
      0.5087, 0.6512, 34.1340, 8.2625, 0.0033, 0.0091, 0.4014, 6.9461, 0.1727, 0.0075, 0.0037,
      0.0238, 0.0068, 0.0010, 0.0164, 0.0032, 0.0017, 0.7583, 0.0144, 0.0574, 0.0170, 1.3646, 0.0288,
      0.0011, 0.1350, 0.01220, 0.0050, 0.0059, 0.0191, 0.1041, 0.0079, 0.3311, 0.0074, 0.0578, 0.0142,
      0.0040, 0.0076, 0.0019, 0.0102, 0.0071, 0.1067, 0.0148, 0.0469, 0.0296, 0.3523, 0.0961, 0.0696,
      0.0281, 0.8229, 3.2168, 1.0023, 0.0958, 0.1223, 0.0220, 0.0131, 0.0077, 0.0461, 0.0137, 0.1117,
      0.0149, 0.0416, 0.0281, 0.0286, 0.0060, 0.0891, 0.0107, 0.0786, 0.0277, 0.0242, 0.0222, 0.1344,
      0.0049, 0.0102, 0.0846, 0.1279, 0.0744, 0.3553, 0.0015, 0.0097, 0.0090, 0.1879, 0.0007, 0.0210,
      0.0089, 0.1697, 0.0100, 0.0208, 0.0007, 0.0007, 0.0007, 0.0028, 0.0032, 0.0119, 0.0002, 0.0101,
      0.0024, 0.0045, 0.0050, 0.0395, 0.0007, 0.0030, 0.0292, 0.1147, 0.0006, 0.0042, 0.0083, 0.0169,
      0.0037, 0.0101, 0.0188, 0.0625, 0.0190, 0.0095, 0.0058, 0.0242, 0.0076, 0.0198, 0.0215, 0.0499,
      0.0138, 0.0022, 0.0044, 0.0057, 0.0101, 0.0195, 0.0125, 0.0138, 0.0292, 0.0136, 0.0257, 0.0530,
      0.0235, 0.0720, 0.0723, 0.1351, 0.1128, 0.8122, 0.0015, 0.0097, 0.0025, 0.0895, 0.0021, 0.0720,
      0.0468, 0.0006, 0.0023, 0.0097, 0.0354, 0.0299, 0.0985, 0.0059, 0.0463, 0.0140, 0.0707, 0.0668,
      0.1809, 0.0030, 0.0657, 0.0633, 0.0076, 0.0056, 0.0727, 0.0037, 0.0050, 0.1086, 0.0045, 0.0279,
      0.0231, 0.0076, 0.0047, 0.0180, 0.0003, 0.0013, 0.0032, 0.0142, 0.0041, 0.0015, 0.0044, 0.0188,
      0.0140, 0.0123, 0.0072, 0.0437, 0.0024, 0.0115, 0.0180, 0.0588, 0.0127, 0.0042, 0.0054, 0.0060,
      0.0124, 0.0156, 0.0922, 0.0140, 0.0042, 0.0006, 0.0070, 0.0060, 0.0071, 0.0058, 0.0061, 0.0031,
      0.0101, 0.0052, 0.0401, 0.0157, 0.0099, 0.0296, 0.0324, 0.0153, 0.0097, 0.0032, 0.0913, 0.0086,
      0.0209, 0.0150, 0.0210, 0.0122
  };

  /**
   * PDB space group ranking (as of Feb. 2017).
   */
  public static final HashMap<String, String> rank = new HashMap<>();

  /**
   * PDB space group names.
   */
  public static final String[] pdbSpaceGroupNames = {
      "P 1", "P -1", "P 1 2 1", "P 1 21 1", "C 1 2 1", "P 1 m 1", "P 1 c 1", "C 1 m 1", "C 1 c 1",
      "P 1 2/m 1", "P 1 21/m 1", "C 1 2/m 1", "P 1 2/c 1", "P 1 21/c 1", "C 1 2/c 1", "P 2 2 2",
      "P 2 2 21", "P 21 21 2", "P 21 21 21", "C 2 2 21", "C 2 2 2", "F 2 2 2", "I 2 2 2",
      "I 21 21 21", "P m m 2", "P m c 21", "P c c 2", "P m a 2", "P c a 21", "P n c 2", "P m n 21",
      "P b a 2", "P n a 21", "P n n 2", "C m m 2", "C m c 21", "C c c 2", "A m m 2", "A b m 2",
      "A m a 2", "A b a 2", "F m m 2", "F d d 2", "I m m 2", "I b a 2", "I m a 2", "P 2/m 2/m 2/m",
      "P 2/n 2/n 2/n", "P 2/c 2/c 2/m", "P 2/b 2/a 2/n", "P 21/m 2/m 2/a", "P 2/n 21/n 2/a",
      "P 2/m 2/n 21/a", "P 21/c 2/c 2/a", "P 21/b 21/a 2/m", "P 21/c 21/c 2/n", "P 2/b 21/c 21/m",
      "P 21/n 21/n 2/m", "P 21/m 21/m 2/n", "P 21/b 2/c 21/n", "P 21/b 21/c 21/a",
      "P 21/n 21/m 21/a", "C 2/m 2/c 21/m", "C 2/m 2/c 21/a", "C 2/m 2/m 2/m", "C 2/c 2/c 2/m",
      "C 2/m 2/m 2/a", "C 2/c 2/c 2/a", "F 2/m 2/m 2/m", "F 2/d 2/d 2/d", "I 2/m 2/m 2/m",
      "I 2/b 2/a 2/m", "I 21/b 21/c 21/a", "I 21/m 21/m 21/a", "P 4", "P 41", "P 42", "P 43", "I 4",
      "I 41", "P -4", "I -4", "P 4/m", "P 42/m", "P 4/n", "P 42/n", "I 4/m", "I 41/a", "P 4 2 2",
      "P 4 21 2", "P 41 2 2", "P 41 21 2", "P 42 2 2", "P 42 21 2", "P 43 2 2", "P 43 21 2",
      "I 4 2 2", "I 41 2 2", "P 4 m m", "P 4 b m", "P 42 c m", "P 42 n m", "P 4 c c", "P 4 n c",
      "P 42 m c", "P 42 b c", "I 4 m m", "I 4 c m", "I 41 m d", "I 41 c d", "P -4 2 m", "P -4 2 c",
      "P -4 21 m", "P -4 21 c", "P -4 m 2", "P -4 c 2", "P -4 b 2", "P -4 n 2", "I -4 m 2",
      "I -4 c 2", "I -4 2 m", "I -4 2 d", "P 4/m 2/m 2/m", "P 4/m 2/c 2/c", "P 4/n 2/b 2/m",
      "P 4/n 2/n 2/c", "P 4/m 21/b 2/m", "P 4/m 21/n 2/c", "P 4/n 21/m 2/m", "P 4/n 2/c 2/c",
      "P 42/m 2/m 2/c", "P 42/m 2/c 2/m", "P 42/n 2/b 2/c", "P 42/n 2/n 2/m", "P 42/m 21/b 2/c",
      "P 42/m 21/n 2/m", "P 42/n 21/m 2/c", "P 42/n 21/c 2/m", "I 4/m 2/m 2/m", "I 4/m 2/c 2/m",
      "I 41/a 2/m 2/d", "I 41/a 2/c 2/d", "P 3", "P 31", "P 32", "H 3", "P -3", "H -3", "P 3 1 2",
      "P 3 2 1", "P 31 1 2", "P 31 2 1", "P 32 1 2", "P 32 2 1", "H 3 2", "P 3 m 1", "P 3 1 m",
      "P 3 c 1", "P 3 1 c", "H 3 m", "H 3 c", "P -3 1 2/m", "P -3 1 2/c", "P -3 2/m 1", "P -3 2/c 1",
      "H -3 2/m", "H -3 2/c", "P 6", "P 61", "P 65", "P 62", "P 64", "P 63", "P -6", "P 6/m",
      "P 63/m", "P 6 2 2", "P 61 2 2", "P 65 2 2", "P 62 2 2", "P 64 2 2", "P 63 2 2", "P 6 m m",
      "P 6 c c", "P 63 c m", "P 63 m c", "P -6 m 2", "P -6 c 2", "P -6 2 m", "P -6 2 c",
      "P 6/m 2/m 2/m", "P 6/m 2/c 2/c", "P 63/m 2/c 2/m", "P 63/m 2/m 2/c", "P 2 3", "F 2 3",
      "I 2 3", "P 21 3", "I 21 3", "P 2/m -3", "P 2/n -3", "F 2/m -3", "F 2/d -3", "I 2/m -3",
      "P 21/a -3", "I 21/a -3", "P 4 3 2", "P 42 3 2", "F 4 3 2", "F 41 3 2", "I 4 3 2", "P 43 3 2",
      "P 41 3 2", "I 41 3 2", "P -4 3 m", "F -4 3 m", "I -4 3 m", "P -4 3 n", "F -4 3 c", "I -4 3 d",
      "P 4/m -3 2/m", "P 4/n -3 2/n", "P 42/m -3 2/n", "P 42/n -3 2/m", "F 4/m -3 2/m",
      "F 4/m -3 2/c", "F 41/d -3 2/m", "F 41/d -3 2/c", "I 4/m -3 2/m", "I 41/a -3 2/d"};

  static {
    rank.put("P 21 21 21", "1");
    rank.put("P 1 21 1", "2");
    rank.put("C 1 2 1", "3");
    rank.put("P 21 21 2", "4");
    rank.put("P 1", "5");
    rank.put("C 2 2 21", "6");
    rank.put("P 43 21 2", "7");
    rank.put("P 41 21 2", "8");
    rank.put("P 32 2 1", "9");
    rank.put("P 31 2 1", "10");
    rank.put("I 2 2 2", "11");
    rank.put("P 61 2 2", "12");
    rank.put("H 3", "13");
    rank.put("H 3 2", "14");
    rank.put("P 65 2 2", "15");
    rank.put("P 61", "16");
    rank.put("P 65", "17");
    rank.put("P 63", "18");
    rank.put("P 41", "19");
    rank.put("P 31", "20");
    rank.put("P 32", "21");
    rank.put("I 4 2 2", "22");
    rank.put("P 43", "23");
    rank.put("I 41 2 2", "24");
    rank.put("P 42 21 2", "25");
    rank.put("P 63 2 2", "26");
    rank.put("I 4", "27");
    rank.put("P 21 3", "28");
    rank.put("I 2 3", "29");
    rank.put("P 2 21 21", "30");
    rank.put("P 4 21 2", "31");
    rank.put("P 3 2 1", "32");
    rank.put("P 62 2 2", "33");
    rank.put("P 64 2 2", "34");
    rank.put("I 41", "35");
    rank.put("P 43 2 2", "36");
    rank.put("I 21 3", "37");
    rank.put("P 6", "38");
    rank.put("I 21 21 21", "39");
    rank.put("P 41 2 2", "40");
    rank.put("C 2 2 2", "41");
    rank.put("P 64", "42");
    rank.put("F 4 3 2", "43");
    rank.put("P 62", "44");
    rank.put("P 1 2 1", "45");
    rank.put("P 3", "46");
    rank.put("I 1 2 1", "47");
    rank.put("P 41 3 2", "48");
    rank.put("P 21 2 21", "49");
    rank.put("F 2 2 2", "50");
    rank.put("P 2 2 21", "51");
    rank.put("P 32 1 2", "52");
    rank.put("I 4 3 2", "53");
    rank.put("F 2 3", "54");
    rank.put("P 4", "55");
    rank.put("P 31 1 2", "56");
    rank.put("P 43 3 2", "57");
    rank.put("P 6 2 2", "58");
    rank.put("P 42", "59");
    rank.put("I 41 3 2", "60");
    rank.put("F 41 3 2", "61");
    rank.put("P 2 3", "62");
    rank.put("P 42 2 2", "63");
    rank.put("P 4 3 2", "64");
    rank.put("P 4 2 2", "65");
    rank.put("B 2", "66");
    rank.put("P 42 3 2", "67");
    rank.put("P -1", "68");
    rank.put("P 3 1 2", "69");
    rank.put("P 2 2 2", "70");
    rank.put("P 1 1 21", "71");
    rank.put("I 21", "72");
    rank.put("R 3 2", "73");
    rank.put("R 3", "74");
    rank.put("C 1 21 1", "75");
    rank.put("P 1 21/c 1", "76");
    rank.put("I 1 21 1", "77");
    rank.put("P 1 21/n 1", "78");
    rank.put("C 1 2/c 1", "79");
    rank.put("P 21 21 2 A", "80");
    rank.put("I 41/a", "81");
    rank.put("F 4 2 2", "82");
    rank.put("A 2", "83");
    rank.put("P 1 1 2", "84");
    rank.put("P b c a", "85");
    rank.put("C 4 21 2", "86");
    rank.put("H -3", "87");
    rank.put("P n n a", "88");
    rank.put("P 2 21 2", "89");
    rank.put("P -3", "90");
    rank.put("I -4 2 d", "91");
    rank.put("I -4 c 2", "92");
    rank.put("B 2 21 2", "93");
    rank.put("B 1 1 2", "94");
    rank.put("A 1", "95");
    rank.put("P 21 2 2", "96");
  }

  /**
   * PDB space group ranking (as of Feb. 2017).
   *
   * <p>Note that a ranking of 97 or higher indicates the space group has not been observed.
   *
   * @param sg a SpaceGroup instance.
   * @return the PDB space group ranking.
   */
  public static int getPDBRank(SpaceGroup sg) {
    String r = rank.getOrDefault(sg.pdbName, "97");
    return parseInt(r);
  }

  /**
   * Return the given space group representation in the CCDC.
   *
   * @param spaceGroup The space group (from 1 to 230).
   * @return Return the percentage.
   */
  public static double getCCDCPercent(int spaceGroup) {
    return csdPercent[spaceGroup - 1];
  }

  /**
   * Returns the space group name for the given PDB name.
   *
   * @param pdbName PDB space group name.
   * @return A short space group name.
   * @since 1.0
   */
  public static String pdb2ShortName(String pdbName) {
    if (pdbName == null) {
      return null;
    }
    String n = pdbName.trim();
    int num = pdbSpaceGroupNames.length;
    for (int i = 0; i < num; i++) {
      if (pdbSpaceGroupNames[i].equalsIgnoreCase(n)) {
        return spaceGroupNames[i];
      }
    }
    return pdbName;
  }

  /**
   * Check if the value of x is between lower and upper (inclusive).
   *
   * @param x     Input integer.
   * @param lower Lower limit for comparison.
   * @param upper Upper limit for comparison.
   * @return True if x is within the limits (inclusive).
   */
  private static boolean isBetween(int x, int lower, int upper) {
    return lower <= x && x <= upper;
  }

  /**
   * Sohncke groups respect chiral molecules (i.e., non-enantiogenic) and include space group numbers:
   * 1, 3-5, 16-24, 75-80, 89-98, 143-146, 149-155, 168-173, 177-182, 195-199 and 207-214.
   *
   * @param number Space group number.
   * @return true if the space group is a Sohncke Group.
   */
  public static boolean isSohnckeGroup(int number) {
    if (number == 1) {
      return true;
    } else if (isBetween(number, 3, 5)) {
      return true;
    } else if (isBetween(number, 16, 24)) {
      return true;
    } else if (isBetween(number, 75, 80)) {
      return true;
    } else if (isBetween(number, 89, 98)) {
      return true;
    } else if (isBetween(number, 143, 146)) {
      return true;
    } else if (isBetween(number, 149, 155)) {
      return true;
    } else if (isBetween(number, 168, 173)) {
      return true;
    } else if (isBetween(number, 177, 182)) {
      return true;
    } else if (isBetween(number, 195, 199)) {
      return true;
    } else {
      return isBetween(number, 207, 214);
    }
  }
}
