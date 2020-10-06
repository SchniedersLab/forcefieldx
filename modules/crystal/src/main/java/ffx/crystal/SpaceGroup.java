// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
// ******************************************************************************
package ffx.crystal;

import static ffx.crystal.CrystalSystem.MONOCLINIC;
import static ffx.crystal.CrystalSystem.TRIGONAL;
import static ffx.crystal.LatticeSystem.MONOCLINIC_LATTICE;
import static ffx.crystal.LatticeSystem.RHOMBOHEDRAL_LATTICE;
import static ffx.crystal.LaueSystem.L121;
import static ffx.crystal.LaueSystem.L32U;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.random;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * The Spacegroup class defines the symmetry of a crystal. There are 230 distinct space groups in
 * three dimensions.
 *
 * @author Michael J. Schnieders
 * @see <ul>
 *     <li><a href="http://it.iucr.org/Ab/" target="_blank"> International Tables for
 *     Crystallography Volume A: Space-group symmetry </a>
 *     <li><a href="http://www.ccp4.ac.uk/html/symlib.html" target="_blank"> CCP4 Symlib </a>
 *     </ul>
 * @since 1.0
 */
public class SpaceGroup {

  /**
   * Names of the 230 three dimensional space groups.
   *
   * @since 1.0
   */
  public static final String[] spaceGroupNames = {
      "P1", "P-1", "P2", "P21", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m", "P21/m", "C2/m", "P2/c",
      "P21/c", "C2/c", "P222", "P2221", "P21212", "P212121", "C2221", "C222", "F222", "I222",
      "I212121", "Pmm2", "Pmc21", "Pcc2", "Pma2", "Pca21", "Pnc2", "Pmn21", "Pba2", "Pna21", "Pnn2",
      "Cmm2", "Cmc21", "Ccc2", "Amm2", "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2",
      "Ima2", "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", "Pbcm",
      "Pnnm",
      "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm", "Cccm", "Cmma", "Ccca", "Fmmm", "Fddd",
      "Immm", "Ibam", "Ibca", "Imma", "P4", "P41", "P42", "P43", "I4", "I41", "P-4", "I-4", "P4/m",
      "P42/m", "P4/n", "P42/n", "I4/m", "I41/a", "P422", "P4212", "P4122", "P41212", "P4222",
      "P42212", "P4322", "P43212", "I422", "I4122", "P4mm", "P4bm", "P42cm", "P42nm", "P4cc", "P4nc",
      "P42mc", "P42bc", "I4mm", "I4cm", "I41md", "I41cd", "P-42m", "P-42c", "P-421m", "P-421c",
      "P-4m2", "P-4c2", "P-4b2", "P-4n2", "I-4m2", "I-4c2", "I-42m", "I-42d", "P4/mmm", "P4/mcc",
      "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", "P42/mmc", "P42/mcm", "P42/nbc",
      "P42/nnm", "P42/mbc", "P42/mnm", "P42/nmc", "P42/ncm", "I4/mmm", "I4/mcm", "I41/amd",
      "I41/acd", "P3", "P31", "P32", "H3", "P-3", "H-3", "P312", "P321", "P3112", "P3121", "P3212",
      "P3221",
      "H32", "P3m1", "P31m", "P3c1", "P31c", "H3m", "H3c", "P-31m", "P-31c", "P-3m1", "P-3c1",
      "H-3m", "H-3c", "P6", "P61", "P65", "P62", "P64", "P63", "P-6", "P6/m", "P63/m", "P622",
      "P6122",
      "P6522", "P6222", "P6422", "P6322", "P6mm", "P6cc", "P63cm", "P63mc", "P-6m2", "P-6c2",
      "P-62m", "P-62c", "P6/mmm", "P6/mcc", "P63/mcm", "P63/mmc", "P23", "F23", "I23", "P213",
      "I213", "Pm-3",
      "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3", "P432", "P4232", "F432", "F4132", "I432",
      "P4332", "P4132", "I4132", "P-43m", "F-43m", "I-43m", "P-43n", "F-43c", "I-43d", "Pm-3m",
      "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c", "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d"
  };

  /**
   * Space group frequency ranking for the 807,190 CSD structures for which the space group is fully
   * defined. Statistics for enantiomorphous space groups are as reported in the CSD. 631,404 (78.2
   * %) of structures adopt centrosymmetric space groups, 175,787 (21.8 %) adopt non-centrosymmetric
   * space groups, and 132,700 (16.4 %) structures adopt Sohncke space groups (as of January 2016.
   */
  public static final double[] csdPercent = {
      0.9492, 24.5313, 0.0176, 5.1773, 0.8456, 0.0026,
      0.4270, 0.0363, 1.0468, 0.0136, 0.4984, 0.5072,
      0.6482, 34.5694, 8.3542, 0.0043, 0.0098, 0.4080,
      7.2397, 0.1749, 0.0072, 0.0033, 0.0235, 0.0074,
      0.0014, 0.0162, 0.0024, 0.0016, 0.7394, 0.0131,
      0.0618, 0.0187, 1.3807, 0.0311, 0.0007, 0.1415,
      0.0130, 0.0027, 0.0063, 0.0180, 0.1077, 0.0085,
      0.3408, 0.0087, 0.0580, 0.0141, 0.0037, 0.0072,
      0.0022, 0.0098, 0.0064, 0.1056, 0.0142, 0.0473,
      0.0275, 0.3528, 0.0991, 0.0704, 0.0306, 0.8507,
      3.3389, 1.0820, 0.0974, 0.1250, 0.0138, 0.0110,
      0.0072, 0.0474, 0.0074, 0.1083, 0.0108, 0.0416,
      0.0285, 0.0212, 0.0056, 0.0919, 0.0104, 0.0767,
      0.0280, 0.0256, 0.0229, 0.1368, 0.0047, 0.0123,
      0.0880, 0.1313, 0.0659, 0.3626, 0.0010, 0.0072,
      0.0084, 0.1950, 0.0011, 0.0173, 0.0073, 0.1690,
      0.0059, 0.0103, 0.0004, 0.0005, 0.0007, 0.0030,
      0.0028, 0.0118, 0.0002, 0.0102, 0.0016, 0.0036,
      0.0046, 0.0380, 0.0007, 0.0035, 0.0261, 0.1243,
      0.0005, 0.0033, 0.0069, 0.0190, 0.0035, 0.0104,
      0.0183, 0.0644, 0.0166, 0.0092, 0.0026, 0.0217,
      0.0079, 0.0111, 0.0217, 0.0474, 0.0050, 0.0020,
      0.0042, 0.0059, 0.0099, 0.0172, 0.0119, 0.0120,
      0.0196, 0.0089, 0.0187, 0.0510, 0.0232, 0.0726,
      0.0712, 0.1264, 0.1152, 0.6463, 0.0019, 0.0098,
      0.0031, 0.0898, 0.0021, 0.0671, 0.0455, 0.0010,
      0.0019, 0.0092, 0.0342, 0.0296, 0.0979, 0.0031,
      0.0422, 0.0107, 0.0665, 0.0446, 0.1603, 0.0030,
      0.0612, 0.0577, 0.0073, 0.0056, 0.0660, 0.0027,
      0.0045, 0.1172, 0.0010, 0.0276, 0.0211, 0.0067,
      0.0043, 0.0128, 0.0005, 0.0011, 0.0033, 0.0133,
      0.0038, 0.0019, 0.0033, 0.0156, 0.0059, 0.0116,
      0.0063, 0.0230, 0.0019, 0.0095, 0.0176, 0.0618,
      0.0072, 0.0024, 0.0038, 0.0048, 0.0109, 0.0126,
      0.0907, 0.0123, 0.0011, 0.0006, 0.0038, 0.0046,
      0.0030, 0.0042, 0.0057, 0.0030, 0.0063, 0.0062,
      0.0334, 0.0146, 0.0108, 0.0300, 0.0159, 0.0128,
      0.0076, 0.0031, 0.0659, 0.0059, 0.0157, 0.0144,
      0.0161, 0.0102
  };

  /** PDB space group ranking (as of Feb. 2017). */
  private static final HashMap<String, String> rank = new HashMap<>();
  /** PDB space group names. */
  private static final String[] pdbSpaceGroupNames = {
      "P 1",
      "P -1",
      "P 1 2 1",
      "P 1 21 1",
      "C 1 2 1",
      "P 1 m 1",
      "P 1 c 1",
      "C 1 m 1",
      "C 1 c 1",
      "P 1 2/m 1",
      "P 1 21/m 1",
      "C 1 2/m 1",
      "P 1 2/c 1",
      "P 1 21/c 1",
      "C 1 2/c 1",
      "P 2 2 2",
      "P 2 2 21",
      "P 21 21 2",
      "P 21 21 21",
      "C 2 2 21",
      "C 2 2 2",
      "F 2 2 2",
      "I 2 2 2",
      "I 21 21 21",
      "P m m 2",
      "P m c 21",
      "P c c 2",
      "P m a 2",
      "P c a 21",
      "P n c 2",
      "P m n 21",
      "P b a 2",
      "P n a 21",
      "P n n 2",
      "C m m 2",
      "C m c 21",
      "C c c 2",
      "A m m 2",
      "A b m 2",
      "A m a 2",
      "A b a 2",
      "F m m 2",
      "F d d 2",
      "I m m 2",
      "I b a 2",
      "I m a 2",
      "P 2/m 2/m 2/m",
      "P 2/n 2/n 2/n",
      "P 2/c 2/c 2/m",
      "P 2/b 2/a 2/n",
      "P 21/m 2/m 2/a",
      "P 2/n 21/n 2/a",
      "P 2/m 2/n 21/a",
      "P 21/c 2/c 2/a",
      "P 21/b 21/a 2/m",
      "P 21/c 21/c 2/n",
      "P 2/b 21/c 21/m",
      "P 21/n 21/n 2/m",
      "P 21/m 21/m 2/n",
      "P 21/b 2/c 21/n",
      "P 21/b 21/c 21/a",
      "P 21/n 21/m 21/a",
      "C 2/m 2/c 21/m",
      "C 2/m 2/c 21/a",
      "C 2/m 2/m 2/m",
      "C 2/c 2/c 2/m",
      "C 2/m 2/m 2/a",
      "C 2/c 2/c 2/a",
      "F 2/m 2/m 2/m",
      "F 2/d 2/d 2/d",
      "I 2/m 2/m 2/m",
      "I 2/b 2/a 2/m",
      "I 21/b 21/c 21/a",
      "I 21/m 21/m 21/a",
      "P 4",
      "P 41",
      "P 42",
      "P 43",
      "I 4",
      "I 41",
      "P -4",
      "I -4",
      "P 4/m",
      "P 42/m",
      "P 4/n",
      "P 42/n",
      "I 4/m",
      "I 41/a",
      "P 4 2 2",
      "P 4 21 2",
      "P 41 2 2",
      "P 41 21 2",
      "P 42 2 2",
      "P 42 21 2",
      "P 43 2 2",
      "P 43 21 2",
      "I 4 2 2",
      "I 41 2 2",
      "P 4 m m",
      "P 4 b m",
      "P 42 c m",
      "P 42 n m",
      "P 4 c c",
      "P 4 n c",
      "P 42 m c",
      "P 42 b c",
      "I 4 m m",
      "I 4 c m",
      "I 41 m d",
      "I 41 c d",
      "P -4 2 m",
      "P -4 2 c",
      "P -4 21 m",
      "P -4 21 c",
      "P -4 m 2",
      "P -4 c 2",
      "P -4 b 2",
      "P -4 n 2",
      "I -4 m 2",
      "I -4 c 2",
      "I -4 2 m",
      "I -4 2 d",
      "P 4/m 2/m 2/m",
      "P 4/m 2/c 2/c",
      "P 4/n 2/b 2/m",
      "P 4/n 2/n 2/c",
      "P 4/m 21/b 2/m",
      "P 4/m 21/n 2/c",
      "P 4/n 21/m 2/m",
      "P 4/n 2/c 2/c",
      "P 42/m 2/m 2/c",
      "P 42/m 2/c 2/m",
      "P 42/n 2/b 2/c",
      "P 42/n 2/n 2/m",
      "P 42/m 21/b 2/c",
      "P 42/m 21/n 2/m",
      "P 42/n 21/m 2/c",
      "P 42/n 21/c 2/m",
      "I 4/m 2/m 2/m",
      "I 4/m 2/c 2/m",
      "I 41/a 2/m 2/d",
      "I 41/a 2/c 2/d",
      "P 3",
      "P 31",
      "P 32",
      "H 3",
      "P -3",
      "H -3",
      "P 3 1 2",
      "P 3 2 1",
      "P 31 1 2",
      "P 31 2 1",
      "P 32 1 2",
      "P 32 2 1",
      "H 3 2",
      "P 3 m 1",
      "P 3 1 m",
      "P 3 c 1",
      "P 3 1 c",
      "H 3 m",
      "H 3 c",
      "P -3 1 2/m",
      "P -3 1 2/c",
      "P -3 2/m 1",
      "P -3 2/c 1",
      "H -3 2/m",
      "H -3 2/c",
      "P 6",
      "P 61",
      "P 65",
      "P 62",
      "P 64",
      "P 63",
      "P -6",
      "P 6/m",
      "P 63/m",
      "P 6 2 2",
      "P 61 2 2",
      "P 65 2 2",
      "P 62 2 2",
      "P 64 2 2",
      "P 63 2 2",
      "P 6 m m",
      "P 6 c c",
      "P 63 c m",
      "P 63 m c",
      "P -6 m 2",
      "P -6 c 2",
      "P -6 2 m",
      "P -6 2 c",
      "P 6/m 2/m 2/m",
      "P 6/m 2/c 2/c",
      "P 63/m 2/c 2/m",
      "P 63/m 2/m 2/c",
      "P 2 3",
      "F 2 3",
      "I 2 3",
      "P 21 3",
      "I 21 3",
      "P 2/m -3",
      "P 2/n -3",
      "F 2/m -3",
      "F 2/d -3",
      "I 2/m -3",
      "P 21/a -3",
      "I 21/a -3",
      "P 4 3 2",
      "P 42 3 2",
      "F 4 3 2",
      "F 41 3 2",
      "I 4 3 2",
      "P 43 3 2",
      "P 41 3 2",
      "I 41 3 2",
      "P -4 3 m",
      "F -4 3 m",
      "I -4 3 m",
      "P -4 3 n",
      "F -4 3 c",
      "I -4 3 d",
      "P 4/m -3 2/m",
      "P 4/n -3 2/n",
      "P 42/m -3 2/n",
      "P 42/n -3 2/m",
      "F 4/m -3 2/m",
      "F 4/m -3 2/c",
      "F 41/d -3 2/m",
      "F 41/d -3 2/c",
      "I 4/m -3 2/m",
      "I 41/a -3 2/d"
  };

  private static final String[] alternativeSpaceGroupNames = {"P21/a"};

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

  /** Space group number. */
  public final int number;
  /** Number of primitive symmetry equivalents. */
  public final int numPrimitiveSymEquiv;
  /** Space group name. */
  public final String shortName;
  /**
   * Point group name. There are 32 distinct points groups, or crystal classes in three dimensions.
   */
  public final String pointGroupName;
  /** Crystal system. */
  public final CrystalSystem crystalSystem;
  /**
   * Tolerance for checking if the crystal system restrictions are satisfied.
   * <p>
   * Set this to 0.0 for strict checking of lattice parameters.
   * <p>
   * For an acetamide crystal minimization, 1.0e-15 was too small a tolerance for equivalent lattice
   * parameters to equate as equal.
   */
  private static final double tolerance = 1.0e-15;
  /** Laue group */
  public final LaueSystem laueSystem;
  /** Space group name under the PDB convention. */
  public final String pdbName;
  /** A List of SymOp instances. */
  public final List<SymOp> symOps;
  /** True for a Sohncke group (non-enantiogenic). */
  public final boolean respectsChirality;
  /** realspace ASU limit */
  public final ASULimit[] limit;

  private final double[] asulim;
  /** Number of symmetry equivalents. */
  private final int numSymEquiv;
  /** Lattice system. */
  public final LatticeSystem latticeSystem;

  /**
   * Immutable SpaceGroup instances are made available only through the factory method so this
   * constructor is private.
   *
   * @param number Space group number.
   * @param numSymEquiv Number of symmetry equivalents.
   * @param numPrimitiveSymEquiv Number of primitive symmetry equivalents.
   * @param shortName Short PDB name.
   * @param pointGroupName Point group name.
   * @param pdbName PDB space group name.
   * @param crystalSystem Crystal system.
   * @param laueSystem Laue System.
   * @param symOps Symmetry operators.
   * @param asulim Assymetric unit limit.
   * @param limit ASULimit instance.
   * @since 1.0
   */
  protected SpaceGroup(
      int number,
      int numSymEquiv,
      int numPrimitiveSymEquiv,
      String shortName,
      String pointGroupName,
      String pdbName,
      CrystalSystem crystalSystem,
      LatticeSystem latticeSystem,
      LaueSystem laueSystem,
      ASULimit[] limit,
      double[] asulim,
      SymOp... symOps) {
    this.number = number;
    this.numSymEquiv = numSymEquiv;
    this.numPrimitiveSymEquiv = numPrimitiveSymEquiv;
    this.shortName = shortName;
    this.pointGroupName = pointGroupName;
    this.crystalSystem = crystalSystem;
    this.latticeSystem = latticeSystem;
    this.laueSystem = laueSystem;
    this.limit = limit;
    this.asulim = asulim;
    this.pdbName = pdbName;
    this.respectsChirality = sohnckeGroup(number);
    this.symOps = new ArrayList<>(Arrays.asList(symOps));

    // ToDo: Crystal systems are subdivided into crystal classes. This info needs to be added to
    // each space group.
  }

  /**
   * Check the given point is within the real space ASU
   *
   * @param asulim a {@link ASULimit} object.
   * @param lim a double.
   * @param x a double.
   * @return True if the point is within the real space ASU, false otherwise
   */
  public static boolean checkASULimit(ASULimit asulim, double lim, double x) {
    switch (asulim) {
      case LT:
        if (x >= 0.0 && x < lim) {
          return true;
        } else {
          return false;
        }
      case LTE:
        if (x >= 0.0 && x <= lim) {
          return true;
        } else {
          return false;
        }
      default:
        return false;
    }
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
    return Integer.parseInt(r);
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
   * Returns a SpaceGroup instance corresponding to the number parameter. If number is not between
   * 1-230 inclusive then null is returned.
   *
   * @param number All 230 3D spacegroups are available.
   * @return The space group corresponding to the given number.
   * @since 1.0
   */
  public static SpaceGroup spaceGroupFactory(int number) {
    if (number > 0 && number <= 50) {
      return SpaceGroupDefinitions.getSpaceGroup1(number);
    } else if (number > 50 && number <= 100) {
      return SpaceGroupDefinitions.getSpaceGroup2(number);
    } else if (number > 100 && number <= 150) {
      return SpaceGroupDefinitions.getSpaceGroup3(number);
    } else if (number > 150 && number <= 200) {
      return SpaceGroupDefinitions.getSpaceGroup4(number);
    } else if (number > 200 && number <= 215) {
      return SpaceGroupDefinitions.getSpaceGroup5(number);
    } else if (number > 215 && number <= 230) {
      return SpaceGroupDefinitions.getSpaceGroup6(number);
    }
    return null;
  }

  /**
   * Return a SpaceGroup based on its name.
   *
   * @param name Available SpaceGroup names are given in the "spaceGroupName" array.
   * @return The space group corresponding to the given number.
   * @since 1.0
   */
  public static SpaceGroup spaceGroupFactory(String name) {
    SpaceGroup spaceGroup = spaceGroupFactory(spaceGroupNumber(name));
    if (spaceGroup == null) {
      spaceGroup = getAlternativeSpaceGroup(name);
    }
    return spaceGroup;
  }

  /**
   * Returns the space group number for a given space group name.
   *
   * @param name The space group name.
   * @return The space group number.
   * @since 1.0
   */
  public static int spaceGroupNumber(String name) {
    if (name == null) {
      return -1;
    }
    String n = name.trim();
    int num = spaceGroupNames.length;
    for (int i = 0; i < num; i++) {
      String s1 = spaceGroupNames[i];
      String s2 = pdbSpaceGroupNames[i];
      if (s1.equalsIgnoreCase(n) || s2.equalsIgnoreCase(n)) {
        return i + 1;
      }
    }
    return -1;
  }

  /**
   * Check the given HKL is valid given the crystal/Laue system
   *
   * @param laueSystem a {@link LaueSystem} object.
   * @param h a int.
   * @param k a int.
   * @param l a int.
   * @return True if the reflection is valid, false otherwise
   */
  static boolean checkLaueRestrictions(LaueSystem laueSystem, int h, int k, int l) {
    switch (laueSystem) {
      case L111:
        return (l > 0 || (l == 0 && (h > 0 || (h == 0 && k >= 0))));
      case L112:
        return (l >= 0 && (h > 0 || (h == 0 && k >= 0)));
      case L121:
        return (k >= 0 && (l > 0 || (l == 0 && h >= 0)));
      case L211:
        return (h >= 0 && (k > 0 || (k == 0 && l >= 0)));
      case L21U:
        return (h + k >= 0 && (l > 0 || (l == 0 && h - k >= 0)));
      case L21V:
        return (l + h >= 0 && (k > 0 || (k == 0 && l - h >= 0)));
      case L21W:
        return (k + l >= 0 && (h > 0 || (h == 0 && k - l >= 0)));
      case L21X:
        return (h - k >= 0 && (l > 0 || (l == 0 && h + k >= 0)));
      case L21Y:
        return (l - h >= 0 && (k > 0 || (k == 0 && l + h >= 0)));
      case L21Z:
        return (k - l >= 0 && (h > 0 || (h == 0 && k + l >= 0)));
      case L222:
        return (h >= 0 && k >= 0 && l >= 0);
      case L22U:
        return (h <= k && h >= -k && l >= 0);
      case L22V:
        return (l <= h && l >= -h && k >= 0);
      case L22W:
        return (k <= l && k >= -l && h >= 0);
      case L114:
        return (l >= 0 && ((h >= 0 && k > 0) || (h == 0 && k == 0)));
      case L141:
        return (k >= 0 && ((l >= 0 && h > 0) || (l == 0 && h == 0)));
      case L411:
        return (h >= 0 && ((k >= 0 && l > 0) || (k == 0 && l == 0)));
      case L224:
        return (h >= k && k >= 0 && l >= 0);
      case L242:
        return (l >= h && h >= 0 && k >= 0);
      case L422:
        return (k >= l && l >= 0 && h >= 0);
      case L113:
        return (h >= 0 && k > 0) || (h == 0 && k == 0 && l >= 0);
      case L131:
        return (l >= 0 && h > 0) || (l == 0 && h == 0 && k >= 0);
      case L311:
        return (k >= 0 && l > 0) || (k == 0 && l == 0 && h >= 0);
      case L11T:
        return (h <= 0 && k > 0) || (h == 0 && k == 0 && l >= 0);
      case L1T1:
        return (l <= 0 && h > 0) || (l == 0 && h == 0 && k >= 0);
      case LT11:
        return (k <= 0 && l > 0) || (k == 0 && l == 0 && h >= 0);
      case L31A:
        return (k - l >= 0 && l - h > 0) || (h == l && k == l && h + k + l >= 0);
      case L31B:
        return (k - l >= 0 && l + h > 0) || (-h == l && k == l && -h + k + l >= 0);
      case L31C:
        return (-k - l >= 0 && l - h > 0) || (h == l && -k == l && h - k + l >= 0);
      case L31D:
        return (k + l >= 0 && -l - h > 0) || (h == -l && k == -l && h + k - l >= 0);
      case L223:
        return (h >= k && k >= 0 && (k > 0 || l >= 0));
      case L232:
        return (l >= h && h >= 0 && (h > 0 || k >= 0));
      case L322:
        return (k >= l && l >= 0 && (l > 0 || h >= 0));
      case L32A:
        return (h >= k && k + l >= h + h && (k + l > h + h || h + k + l >= 0));
      case L32B:
        return (-h >= k && k + l >= -h - h && (k + l > -h - h || -h + k + l >= 0));
      case L32C:
        return (h >= -k && -k + l >= h + h && (-k + l > h + h || h - k + l >= 0));
      case L32D:
        return (h >= k && k - l >= h + h && (k - l > h + h || h + k - l >= 0));
      case L32U:
        return (h >= k && k >= 0 && (h > k || l >= 0));
      case L32V:
        return (k >= l && l >= 0 && (k > l || h >= 0));
      case L32W:
        return (l >= h && h >= 0 && (l > h || k >= 0));
      case L32X:
        return (-h >= k && k >= 0 && (-h > k || l >= 0));
      case L32Y:
        return (-k >= l && l >= 0 && (-k > l || h >= 0));
      case L32Z:
        return (-l >= h && h >= 0 && (-l > h || k >= 0));
      case LM3B:
        return (h >= 0 && ((l >= h && k > h) || (l == h && k == h)));
      case LM3M:
        return (k >= l && l >= h && h >= 0);
      default:
        assert (2 != 2);
        return false;
    }
  }

  /**
   * Check that the lattice parameters satisfy the restrictions of the lattice systems.
   *
   * @param latticeSystem a {@link LatticeSystem} object.
   * @param a the a-axis length.
   * @param b the b-axis length.
   * @param c the c-axis length.
   * @param alpha the alpha angle.
   * @param beta the beta angle.
   * @param gamma the gamma angle.
   * @return True if the restrictions are satisfied, false otherwise.
   */
  static boolean checkLatticeSystemRestrictions(LatticeSystem latticeSystem,
      double a, double b, double c, double alpha, double beta, double gamma) {
    switch (latticeSystem) {
      case TRICLINIC_LATTICE:
        // No restrictions.
        return true;
      case MONOCLINIC_LATTICE:
        // alpha = gamma = 90
        return check(alpha, 90.0) && check(gamma, 90.0);
      case ORTHORHOMBIC_LATTICE:
        // alpha = beta = gamma = 90
        return check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      case TETRAGONAL_LATTICE:
        // a = b, alpha = beta = gamma = 90
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      case RHOMBOHEDRAL_LATTICE:
        // a = b = c, alpha = beta = gamma.
        return check(a, b) && check(b, c) && check(alpha, beta) && check(beta, gamma);
      case HEXAGONAL_LATTICE:
        // a = b, alpha = beta = 90, gamma = 120
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 120.0);
      case CUBIC_LATTICE:
        // a = b = c; alpha = beta = gamma = 90
        return check(a, b) && check(b, c) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma,
            90.0);
      default:
        assert (2 != 2);
        return false;
    }
  }

  /**
   * If the two passed values are the same, within the tolerance, return true.
   *
   * @param x1 First value.
   * @param x2 Second value.
   * @return Return true if the two values are the same within specified tolerance.
   */
  public static boolean check(double x1, double x2) {
    return abs(x1 - x2) < tolerance;
  }

  private static SpaceGroup getAlternativeSpaceGroup(String name) {
    SpaceGroup spaceGroup = null;
    if (name.equalsIgnoreCase("P21/a")) {
      spaceGroup =
          new SpaceGroup(
              14,
              4,
              4,
              "P21/a",
              "PG2/m",
              "P 1 21/a 1",
              MONOCLINIC,
              MONOCLINIC_LATTICE,
              L121,
              new ASULimit[] {ASULimit.LT, ASULimit.LT, ASULimit.LT},
              new double[] {-1.0, -1.0, -1.0},
              new SymOp(SymOp.Rot_X_Y_Z, SymOp.Tr_0_0_0),
              new SymOp(SymOp.Rot_mX_mY_Z, SymOp.Tr_12_0_12),
              new SymOp(SymOp.Rot_mX_mY_mZ, SymOp.Tr_0_0_0),
              new SymOp(SymOp.Rot_X_Y_mZ, SymOp.Tr_12_0_12));
    } else if (name.equalsIgnoreCase("R3c")) {
      spaceGroup =
          new SpaceGroup(
              161,
              6,
              6,
              "R3c",
              "PG3m",
              "R 3 c",
              TRIGONAL,
              RHOMBOHEDRAL_LATTICE,
              L32U,
              new ASULimit[] {ASULimit.LT, ASULimit.LT, ASULimit.LT},
              new double[] {-1.0, -1.0, -1.0},
              new SymOp(SymOp.Rot_X_Y_Z, SymOp.Tr_0_0_0),
              new SymOp(SymOp.Rot_Z_X_Y, SymOp.Tr_0_0_0),
              new SymOp(SymOp.Rot_Y_Z_X, SymOp.Tr_0_0_0),
              new SymOp(SymOp.Rot_Y_X_Z, SymOp.Tr_12_12_12),
              new SymOp(SymOp.Rot_Z_Y_X, SymOp.Tr_12_12_12),
              new SymOp(SymOp.Rot_X_Z_Y, SymOp.Tr_12_12_12));
    }

    return spaceGroup;
  }

  private static boolean isBetween(int x, int lower, int upper) {
    return lower <= x && x <= upper;
  }

  /**
   * Reset lattice parameters for the given lattice systems.
   *
   * @param latticeSystem a {@link LatticeSystem} object.
   * @return New unit cell parameters.
   */
  static double[] resetUnitCellParams(LatticeSystem latticeSystem) {
    double alpha = 60.0 + random() * 60.0;
    double beta = 60.0 + random() * 60.0;
    double gamma = 60.0 + random() * 60.0;
    double[] params = {0.1 + random(), 0.1 + random(), 0.1 + random(), alpha, beta, gamma};
    switch (latticeSystem) {
      case TRICLINIC_LATTICE:
        break;
      case MONOCLINIC_LATTICE:
        // alpha = gamma = 90
        params[3] = 90.0;
        params[5] = 90.0;
        break;
      case ORTHORHOMBIC_LATTICE:
        // alpha = beta = gamma = 90
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
      case TETRAGONAL_LATTICE:
        // a = b, alpha = beta = gamma = 90
        double ab = 0.5 * (params[0] + params[1]);
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
      case RHOMBOHEDRAL_LATTICE:
        // a = b = c, alpha = beta = gamma.
        double abc = (params[0] + params[1] + params[2]) / 3.0;
        double angles = (params[3] + params[4] + params[5]) / 3.0;
        params[0] = abc;
        params[1] = abc;
        params[2] = abc;
        params[3] = angles;
        params[4] = angles;
        params[5] = angles;
        break;
      case HEXAGONAL_LATTICE:
        // a = b, alpha = beta = 90, gamma = 120
        ab = 0.5 * (params[0] + params[1]);
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 120.0;
        break;
      case CUBIC_LATTICE:
      default:
        // a = b = c, alpha = beta = gamma = 90
        abc = (params[0] + params[1] + params[2]) / 3.0;
        params[0] = abc;
        params[1] = abc;
        params[2] = abc;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
    }
    return params;
  }

  /**
   * Sohncke groups respect chiral molecules (i.e. non-enantiogenic) and include space group numbers:
   * 1, 3-5, 16-24, 75-80, 89-98, 143-146, 149-155, 168-173, 177-182, 195-199 and 207-214.
   *
   * @param number Space group number.
   * @return true if the space group is a Sohncke Group.
   */
  private static boolean sohnckeGroup(int number) {
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
    } else if (isBetween(number, 207, 214)) {
      return true;
    }

    return false;
  }

  /**
   * Return the number of symmetry operators.
   *
   * @return the number of symmetry operators.
   * @since 1.0
   */
  public int getNumberOfSymOps() {
    return symOps.size();
  }

  /**
   * Return the ith symmetry operator.
   *
   * @param i the symmetry operator number.
   * @return the SymOp
   * @since 1.0
   */
  public SymOp getSymOp(int i) {
    return symOps.get(i);
  }

  public double[] randomUnitCellParams() {
    return resetUnitCellParams(latticeSystem);
  }

  public boolean respectsChirality() {
    return respectsChirality;
  }

}
