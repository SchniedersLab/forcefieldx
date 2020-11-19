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

/**
 * Enumeration of the different Laue systems. Some are only used for nonstandard cells.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public enum LaueSystem {
  L111, L112, L121, L211, L21U, L21V, L21W, L21X, L21Y, L21Z, L222, L22U, L22V, L22W, L114, L141,
  L411, L224, L242, L422, L113, L131, L311, L11T, L1T1, LT11, L31A, L31B, L31C, L31D, L223, L232,
  L322, L32A, L32B, L32C, L32D, L32U, L32V, L32W, L32X, L32Y, L32Z, LM3B, LM3M;

  /**
   * Check the given HKL is valid given the Laue system.
   *
   * @param h a int.
   * @param k a int.
   * @param l a int.
   * @return True if the reflection is valid, false otherwise.
   */
  public boolean checkRestrictions(int h, int k, int l) {
    switch (this) {
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
}
