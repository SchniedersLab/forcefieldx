// ******************************************************************************
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
// ******************************************************************************
package ffx.crystal;

/**
 * Enumeration of the different Laue systems. Some are only used for nonstandard cells.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public enum LaueSystem {
  /**
   * Laue System 111.
   */
  L111,
  /**
   * Laue System 112.
   */
  L112,
  /**
   * Laue System 121.
   */
  L121,
  /**
   * Laue System 211.
   */
  L211,
  /**
   * Laue System 21U.
   */
  L21U,
  /**
   * Laue System 21V.
   */
  L21V,
  /**
   * Laue System 21W.
   */
  L21W,
  /**
   * Laue System 21X.
   */
  L21X,
  /**
   * Laue System 21Y.
   */
  L21Y,
  /**
   * Laue System 21Z.
   */
  L21Z,
  /**
   * Laue System 222.
   */
  L222,
  /**
   * Laue System 22U.
   */
  L22U,
  /**
   * Laue System 22V.
   */
  L22V,
  /**
   * Laue System 22W.
   */
  L22W,
  /**
   * Laue System 114.
   */
  L114,
  /**
   * Laue System 141.
   */
  L141,
  /**
   * Laue System 411.
   */
  L411,
  /**
   * Laue System 224.
   */
  L224,
  /**
   * Laue System 242.
   */
  L242,
  /**
   * Laue System 422.
   */
  L422,
  /**
   * Laue System 113.
   */
  L113,
  /**
   * Laue System 131.
   */
  L131,
  /**
   * Laue System 311.
   */
  L311,
  /**
   * Laue System 11T.
   */
  L11T,
  /**
   * Laue System 1T1.
   */
  L1T1,
  /**
   * Laue System T11.
   */
  LT11,
  /**
   * Laue System 31A.
   */
  L31A,
  /**
   * Laue System 31B.
   */
  L31B,
  /**
   * Laue System 31C.
   */
  L31C,
  /**
   * Laue System 31D.
   */
  L31D,
  /**
   * Laue System 223.
   */
  L223,
  /**
   * Laue System 232.
   */
  L232,
  /**
   * Laue System 322.
   */
  L322,
  /**
   * Laue System 32A.
   */
  L32A,
  /**
   * Laue System 32B.
   */
  L32B,
  /**
   * Laue System 32C.
   */
  L32C,
  /**
   * Laue System 32D.
   */
  L32D,
  /**
   * Laue System 32U.
   */
  L32U,
  /**
   * Laue System 32V.
   */
  L32V,
  /**
   * Laue System 32W.
   */
  L32W,
  /**
   * Laue System 32X.
   */
  L32X,
  /**
   * Laue System 32Y.
   */
  L32Y,
  /**
   * Laue System 32Z.
   */
  L32Z,
  /**
   * Laue System M3B.
   */
  LM3B,
  /**
   * Laue System M3M.
   */
  LM3M;

  /**
   * Check the given HKL is valid given the Laue system.
   *
   * @param h an int.
   * @param k an int.
   * @param l an int.
   * @return True if the reflection is valid, false otherwise.
   */
  public boolean checkRestrictions(int h, int k, int l) {
    switch (this) {
      case L111 -> {
        return (l > 0 || (l == 0 && (h > 0 || (h == 0 && k >= 0))));
      }
      case L112 -> {
        return (l >= 0 && (h > 0 || (h == 0 && k >= 0)));
      }
      case L121 -> {
        return (k >= 0 && (l > 0 || (l == 0 && h >= 0)));
      }
      case L211 -> {
        return (h >= 0 && (k > 0 || (k == 0 && l >= 0)));
      }
      case L21U -> {
        return (h + k >= 0 && (l > 0 || (l == 0 && h - k >= 0)));
      }
      case L21V -> {
        return (l + h >= 0 && (k > 0 || (k == 0 && l - h >= 0)));
      }
      case L21W -> {
        return (k + l >= 0 && (h > 0 || (h == 0 && k - l >= 0)));
      }
      case L21X -> {
        return (h - k >= 0 && (l > 0 || (l == 0 && h + k >= 0)));
      }
      case L21Y -> {
        return (l - h >= 0 && (k > 0 || (k == 0 && l + h >= 0)));
      }
      case L21Z -> {
        return (k - l >= 0 && (h > 0 || (h == 0 && k + l >= 0)));
      }
      case L222 -> {
        return (h >= 0 && k >= 0 && l >= 0);
      }
      case L22U -> {
        return (h <= k && h >= -k && l >= 0);
      }
      case L22V -> {
        return (l <= h && l >= -h && k >= 0);
      }
      case L22W -> {
        return (k <= l && k >= -l && h >= 0);
      }
      case L114 -> {
        return (l >= 0 && ((h >= 0 && k > 0) || (h == 0 && k == 0)));
      }
      case L141 -> {
        return (k >= 0 && ((l >= 0 && h > 0) || (l == 0 && h == 0)));
      }
      case L411 -> {
        return (h >= 0 && ((k >= 0 && l > 0) || (k == 0 && l == 0)));
      }
      case L224 -> {
        return (h >= k && k >= 0 && l >= 0);
      }
      case L242 -> {
        return (l >= h && h >= 0 && k >= 0);
      }
      case L422 -> {
        return (k >= l && l >= 0 && h >= 0);
      }
      case L113 -> {
        return (h >= 0 && k > 0) || (h == 0 && k == 0 && l >= 0);
      }
      case L131 -> {
        return (l >= 0 && h > 0) || (l == 0 && h == 0 && k >= 0);
      }
      case L311 -> {
        return (k >= 0 && l > 0) || (k == 0 && l == 0 && h >= 0);
      }
      case L11T -> {
        return (h <= 0 && k > 0) || (h == 0 && k == 0 && l >= 0);
      }
      case L1T1 -> {
        return (l <= 0 && h > 0) || (l == 0 && h == 0 && k >= 0);
      }
      case LT11 -> {
        return (k <= 0 && l > 0) || (k == 0 && l == 0 && h >= 0);
      }
      case L31A -> {
        return (k - l >= 0 && l - h > 0) || (h == l && k == l && h + k + l >= 0);
      }
      case L31B -> {
        return (k - l >= 0 && l + h > 0) || (-h == l && k == l && -h + k + l >= 0);
      }
      case L31C -> {
        return (-k - l >= 0 && l - h > 0) || (h == l && -k == l && h - k + l >= 0);
      }
      case L31D -> {
        return (k + l >= 0 && -l - h > 0) || (h == -l && k == -l && h + k - l >= 0);
      }
      case L223 -> {
        return (h >= k && k >= 0 && (k > 0 || l >= 0));
      }
      case L232 -> {
        return (l >= h && h >= 0 && (h > 0 || k >= 0));
      }
      case L322 -> {
        return (k >= l && l >= 0 && (l > 0 || h >= 0));
      }
      case L32A -> {
        return (h >= k && k + l >= h + h && (k + l > h + h || h + k + l >= 0));
      }
      case L32B -> {
        return (-h >= k && k + l >= -h - h && (k + l > -h - h || -h + k + l >= 0));
      }
      case L32C -> {
        return (h >= -k && -k + l >= h + h && (-k + l > h + h || h - k + l >= 0));
      }
      case L32D -> {
        return (h >= k && k - l >= h + h && (k - l > h + h || h + k - l >= 0));
      }
      case L32U -> {
        return (h >= k && k >= 0 && (h > k || l >= 0));
      }
      case L32V -> {
        return (k >= l && l >= 0 && (k > l || h >= 0));
      }
      case L32W -> {
        return (l >= h && h >= 0 && (l > h || k >= 0));
      }
      case L32X -> {
        return (-h >= k && k >= 0 && (-h > k || l >= 0));
      }
      case L32Y -> {
        return (-k >= l && l >= 0 && (-k > l || h >= 0));
      }
      case L32Z -> {
        return (-l >= h && h >= 0 && (-l > h || k >= 0));
      }
      case LM3B -> {
        return (h >= 0 && ((l >= h && k > h) || (l == h && k == h)));
      }
      case LM3M -> {
        return (k >= l && l >= h && h >= 0);
      }
      default -> {
        return false;
      }
    }
  }
}
