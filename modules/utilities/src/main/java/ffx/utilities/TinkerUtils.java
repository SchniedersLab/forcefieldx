//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.utilities;

import static java.lang.String.format;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A collection of Utility methods for compatibility with Tinker.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TinkerUtils {

  private static final Logger logger = Logger.getLogger(TinkerUtils.class.getName());

  /**
   * Tinker ranges begin with a negative number.
   *
   * <p>Here is the description from the Tinker Guide: "Several keywords take a list of integer
   * values (atom numbers, for example) as modifiers. For these keywords the integers can simply be
   * listed explicitly and separated by spaces, commas or tabs. If a range of numbers is desired, it
   * can be specified by listing the negative of the first number of the range, followed by a
   * separator and the last number of the range. For example, the keyword line ACTIVE 4 -9 17 23
   * could be used to add atoms 4, 9 through 17, and 23 to the set of active atoms during a TINKER
   * calculation."
   */
  private static final Pattern atomSelectionStartPattern = Pattern.compile("-(\\d+)");

  /** Tinker ranges end with a positive number. Positive numbers may also be a Singleton. */
  private static final Pattern atomSingletonPattern = Pattern.compile("(\\d+)");

  /**
   * Parse a Tinker selection list. No checking is done to make sure range i
   *
   * @param tokens A list of tokens to parse.
   * @param offset A constant offset that is added to each selected integer (i.e. an offset of -1
   *     will index atoms from 0 to n-1).
   * @param maxEntries The maximum number of entries to parse (e.g. "4 -9 17 23" is 3 entries). A
   *     value less than 1 includes all entries.
   * @return A list of selected integers.
   */
  public static List<Integer> parseTinkerAtomList(List<String> tokens, int offset, int maxEntries) {
    // Add selected atom to this list.
    List<Integer> list = new ArrayList<>();

    // No entries parsed to start.
    int nParsed = 0;

    // Loop over tokens.
    int n = tokens.size();
    for (int i = 0; i < n; i++) {
      // Check for Tinker range that begins with a negative number.
      Matcher m = atomSelectionStartPattern.matcher(tokens.get(i));
      if (m.matches()) {
        int start = Integer.parseInt(m.group(1));
        ++i;
        if (i == n) {
          logger.info(
              format(
                  " Attempted to parse -%d as a Tinker-style range, but it was ignored because it was the last token provided.",
                  start));
          continue;
        }
        int end = Integer.parseInt(tokens.get(i));
        if (end < start) {
          logger.info(
              format(
                  " Attempted to parse -%d to %d as a Tinker-style range, which is ignored due to being invalid.",
                  start, end));
          continue;
        }
        for (int j = start; j <= end; j++) {
          list.add(j + offset);
        }
        nParsed++;
      } else {
        // Check for a singleton entry.
        m = atomSingletonPattern.matcher(tokens.get(i));
        if (m.matches()) {
          int start = Integer.parseInt(m.group(1));
          list.add(start + offset);
          nParsed++;
        } else {
          logger.info(
              format(
                  " Attempted to parse %s as a Tinker-style range, but it was not recognized and ignored.",
                  tokens.get(i)));
        }
      }

      // Enforce the maximum number of entries to parse.
      if (maxEntries > 0 && nParsed >= maxEntries) {
        break;
      }
    }
    return list;
  }

  /**
   * Get the previous file based on the TINKER scheme.
   *
   * @param file Root file.
   * @return New File created previously to the passed file.
   */
  public static File previousVersion(File file) {
    if (file == null) {
      return null;
    }
    String fileName = file.getAbsolutePath();
    int dot = file.getAbsolutePath().lastIndexOf(".");
    int under = file.getAbsolutePath().lastIndexOf("_");
    File newFile = file;
    if (under > dot) {
      String name = fileName.substring(0, under);
      newFile = new File(name);
    }
    File baseFile = newFile;
    File previousFile = null;
    int i = 1;
    while (newFile.exists()) {
      i = i + 1;
      previousFile = newFile;
      newFile = tinkerFile(baseFile, i);
    }
    return previousFile;
  }

  /**
   * This follows the TINKER file versioning scheme.
   *
   * @param file File to find a version for.
   * @return File Versioned File.
   */
  public static File version(File file) {
    if (file == null) {
      return null;
    }
    if (!file.exists()) {
      return file;
    }
    String fileName = file.getAbsolutePath();
    int dot = file.getAbsolutePath().lastIndexOf(".");
    int under = file.getAbsolutePath().lastIndexOf("_");
    File newFile = file;
    if (under > dot) {
      String name = fileName.substring(0, under);
      newFile = new File(name);
    }
    File oldFile = newFile;
    int i = 1;
    while (newFile.exists()) {
      i = i + 1;
      newFile = tinkerFile(oldFile, i);
    }
    return newFile;
  }

  /**
   * Append the integer i to filename of File.
   *
   * @param file Root file.
   * @param i Integer to append.
   * @return New File with integer appended.
   */
  private static File tinkerFile(File file, int i) {
    return new File(file.getAbsolutePath() + '_' + i);
  }
}
