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
package ffx.potential.parsers;

import static ffx.potential.bonded.Bond.logNoBondType;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.OptionalDouble;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.jogamp.vecmath.Vector3d;

/**
 * The XYZFilter class parses TINKER Cartesian coordinate (*.XYZ) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class XYZFilter extends SystemFilter {

  private static final Logger logger = Logger.getLogger(XYZFilter.class.getName());
  private BufferedReader bufferedReader = null;
  private int snapShot;
  private String remarkLine;

  /**
   * Constructor for XYZFilter.
   *
   * @param files a {@link java.util.List} object.
   * @param system a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   */
  public XYZFilter(
      List<File> files,
      MolecularAssembly system,
      ForceField forceField,
      CompositeConfiguration properties) {
    super(files, system, forceField, properties);
    this.fileType = FileType.XYZ;
  }

  /**
   * Constructor for XYZFilter.
   *
   * @param file a {@link java.io.File} object.
   * @param system a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   */
  public XYZFilter(
      File file,
      MolecularAssembly system,
      ForceField forceField,
      CompositeConfiguration properties) {
    super(file, system, forceField, properties);
    this.fileType = FileType.XYZ;
  }

  /**
   * readOnto
   *
   * @param newFile a {@link java.io.File} object.
   * @param oldSystem a {@link ffx.potential.MolecularAssembly} object.
   * @return a boolean.
   */
  public static boolean readOnto(File newFile, MolecularAssembly oldSystem) {
    if (newFile == null || !newFile.exists() || oldSystem == null) {
      return false;
    }
    try (BufferedReader br = new BufferedReader(new FileReader(newFile))) {
      String data = br.readLine();
      if (data == null) {
        return false;
      }
      String[] tokens = data.trim().split(" +");
      int num_atoms = parseInt(tokens[0]);
      if (num_atoms != oldSystem.getAtomList().size()) {
        return false;
      }

      br.mark(10000);
      data = br.readLine();
      if (!readPBC(data, oldSystem)) {
        br.reset();
      }

      double[][] d = new double[num_atoms][3];
      for (int i = 0; i < num_atoms; i++) {
        if (!br.ready()) {
          return false;
        }
        data = br.readLine();
        if (data == null) {
          logger.warning(format(" Check atom %d.", (i + 1)));
          return false;
        }
        tokens = data.trim().split(" +");
        if (tokens.length < 6) {
          logger.warning(format(" Check atom %d.", (i + 1)));
          return false;
        }
        d[i][0] = parseDouble(tokens[2]);
        d[i][1] = parseDouble(tokens[3]);
        d[i][2] = parseDouble(tokens[4]);
      }
      List<Atom> atoms = oldSystem.getAtomList();
      for (Atom a : atoms) {
        int index = a.getIndex() - 1;
        a.setXYZ(d[index]);
      }
      oldSystem.center();
      oldSystem.setFile(newFile);
      return true;
    } catch (Exception e) {
      return false;
    }
  }

  private static boolean firstTokenIsInteger(String data) {
    if (data == null) {
      return false;
    }

    // Check for a blank line.
    data = data.trim();
    if (data.equals("")) {
      return false;
    }

    // Check if the first token in an integer.
    try {
      String[] tokens = data.split(" +");
      parseInt(tokens[0]);
      return true;
    } catch (NumberFormatException e) {
      return false;
    }
  }

  /**
   * Attempt to parse the String as unit cell parameters.
   *
   * @param data The String to parse.
   * @return false if the first token in the String is an integer and true otherwise.
   */
  private static boolean readPBC(String data, MolecularAssembly activeMolecularAssembly) {
    if (firstTokenIsInteger(data)) {
      return false;
    }

    String[] tokens = data.trim().split(" +");
    if (tokens.length == 6) {
      CompositeConfiguration config = activeMolecularAssembly.getProperties();
      double a = parseDouble(tokens[0]);
      double b = parseDouble(tokens[1]);
      double c = parseDouble(tokens[2]);
      double alpha = parseDouble(tokens[3]);
      double beta = parseDouble(tokens[4]);
      double gamma = parseDouble(tokens[5]);
      config.setProperty("a-axis", a);
      config.setProperty("b-axis", b);
      config.setProperty("c-axis", c);
      config.setProperty("alpha", alpha);
      config.setProperty("beta", beta);
      config.setProperty("gamma", gamma);

      Crystal crystal = activeMolecularAssembly.getCrystal();
      if (crystal != null) {
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
      }
    }
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void closeReader() {
    if (bufferedReader != null) {
      try {
        bufferedReader.close();
      } catch (IOException ex) {
        logger.warning(format(" Exception in closing XYZ filter: %s", ex.toString()));
      }
    }
  }

  @Override
  public int countNumModels() {
    File xyzFile = activeMolecularAssembly.getFile();
    int nAtoms = activeMolecularAssembly.getAtomArray().length;
    Pattern crystInfoPattern =
        Pattern.compile(
            "^ *(?:[0-9]+\\.[0-9]+ +){3}(?:-?[0-9]+\\.[0-9]+ +){2}(?:-?[0-9]+\\.[0-9]+) *$");

    try (BufferedReader br = new BufferedReader(new FileReader(xyzFile))) {
      String line = br.readLine();
      int nSnaps = 0;
      // For each header line, will read either X or X+1 lines, where X is the number of atoms.
      while (line != null) {
        assert parseInt(line.trim().split("\\s+")[0]) == nAtoms;
        // Read either the crystal information *or* the first line of the snapshot.
        line = br.readLine();
        Matcher m = crystInfoPattern.matcher(line);
        if (m.matches()) {
          // If this is crystal information, move onto the first line of the snapshot.
          br.readLine();
        }
        // Read lines 2-X of the XYZ.
        for (int i = 1; i < nAtoms; i++) {
          br.readLine();
        }

        ++nSnaps;
        line = br.readLine();
      }
      return nSnaps;
    } catch (Exception ex) {
      logger.log(
          Level.WARNING, String.format(" Exception reading trajectory file %s: %s", xyzFile, ex));
      return 1;
    }
  }

  /** {@inheritDoc} */
  @Override
  public OptionalDouble getLastReadLambda() {
    String[] toks = remarkLine.split("\\s+");
    int nToks = toks.length;
    for (int i = 0; i < (nToks - 1); i++) {
      if (toks[i].equals("Lambda:")) {
        return OptionalDouble.of(Double.parseDouble(toks[i + 1]));
      }
    }
    return OptionalDouble.empty();
  }

  @Override
  public String[] getRemarkLines() {
    return new String[] {remarkLine};
  }

  @Override
  public int getSnapshot() {
    return snapShot;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Parse the XYZ File
   */
  @Override
  public boolean readFile() {
    File xyzFile = activeMolecularAssembly.getFile();

    if (forceField == null) {
      logger.warning(format(" No force field is associated with %s.", xyzFile.toString()));
      return false;
    }
    try (BufferedReader br = new BufferedReader(new FileReader(xyzFile))) {
      String data = br.readLine();
      // Read blank lines at the top of the file
      while (data != null && data.trim().equals("")) {
        data = br.readLine();
      }
      if (data == null) {
        return false;
      }
      String[] tokens = data.trim().split(" +", 2);
      int numberOfAtoms = parseInt(tokens[0]);
      if (numberOfAtoms < 1) {
        return false;
      }
      if (tokens.length == 2) {
        getActiveMolecularSystem().setName(tokens[1]);
      }
      logger.info(format(" Opening %s with %d atoms\n", xyzFile.getName(), numberOfAtoms));
      remarkLine = data.trim();

      // The header line is reasonable. Check for periodic box dimensions.
      br.mark(10000);
      data = br.readLine();
      if (!readPBC(data, activeMolecularAssembly)) {
        br.reset();
      }

      // Prepare to parse atom lines.
      HashMap<Integer, Integer> labelHash = new HashMap<>();
      int[] label = new int[numberOfAtoms];
      int[][] bonds = new int[numberOfAtoms][8];
      double[][] d = new double[numberOfAtoms][3];
      boolean renumber = false;
      atomList = new ArrayList<>();
      // Loop over the expected number of atoms.
      for (int i = 0; i < numberOfAtoms; i++) {
        if (!br.ready()) {
          return false;
        }
        data = br.readLine();
        if (data == null) {
          logger.warning(
              format(
                  " Check atom %d in %s.", (i + 1), activeMolecularAssembly.getFile().getName()));
          return false;
        }
        tokens = data.trim().split(" +");
        if (tokens.length < 6) {
          logger.warning(
              format(
                  " Check atom %d in %s.", (i + 1), activeMolecularAssembly.getFile().getName()));
          return false;
        }
        // Valid number of tokens, so try to parse this line.
        label[i] = parseInt(tokens[0]);
        // Check for valid atom numbering, or flag for re-numbering.
        if (label[i] != i + 1) {
          renumber = true;
        }
        String atomName = tokens[1];
        d[i][0] = parseDouble(tokens[2]);
        d[i][1] = parseDouble(tokens[3]);
        d[i][2] = parseDouble(tokens[4]);
        int type = parseInt(tokens[5]);
        AtomType atomType = forceField.getAtomType(Integer.toString(type));
        if (atomType == null) {
          StringBuilder message = new StringBuilder("Check atom type ");
          message.append(type).append(" for Atom ").append(i + 1);
          message.append(" in ").append(activeMolecularAssembly.getFile().getName());
          logger.warning(message.toString());
          return false;
        }
        Atom a = new Atom(i + 1, atomName, atomType, d[i]);
        atomList.add(a);
        // Bond Data
        int numberOfBonds = tokens.length - 6;
        for (int b = 0; b < 8; b++) {
          if (b < numberOfBonds) {
            int bond = parseInt(tokens[6 + b]);
            bonds[i][b] = bond;
          } else {
            bonds[i][b] = 0;
          }
        }
      }
      // Check if this is an archive.
      if (br.ready()) {
        // Read past blank lines between archive files
        data = br.readLine().trim();
        while (data.equals("") && br.ready()) {
          data = br.readLine().trim();
        }
        tokens = data.split(" +", 2);
        if (tokens.length > 0) {
          try {
            int archiveNumberOfAtoms = parseInt(tokens[0]);
            if (archiveNumberOfAtoms == numberOfAtoms) {
              setType(FileType.ARC);
            }
          } catch (NumberFormatException e) {
            //
          }
        }
      }
      // Try to renumber
      if (renumber) {
        for (int i = 0; i < numberOfAtoms; i++) {
          if (labelHash.containsKey(label[i])) {
            logger.warning(format(" Two atoms have the same index: %d.", label[i]));
            return false;
          }
          labelHash.put(label[i], i + 1);
        }
        for (int i = 0; i < numberOfAtoms; i++) {
          int j = -1;
          while (j < 3 && bonds[i][++j] > 0) {
            bonds[i][j] = labelHash.get(bonds[i][j]);
          }
        }
      }
      bondList = new ArrayList<>();
      int[] c = new int[2];
      for (int a1 = 1; a1 <= numberOfAtoms; a1++) {
        int j = -1;
        while (j < 7 && bonds[a1 - 1][++j] > 0) {
          int a2 = bonds[a1 - 1][j];
          if (a1 < a2) {
            if (a2 > numberOfAtoms) {
              logger.warning(
                  format(
                      " Check the bond between %d and %d in %s.",
                      a1, a2, activeMolecularAssembly.getFile().getName()));
              return false;
            }
            // Check for bidirectional connection
            boolean bidirectional = false;
            int k = -1;
            while (k < 7 && bonds[a2 - 1][++k] > 0) {
              int a3 = bonds[a2 - 1][k];
              if (a3 == a1) {
                bidirectional = true;
                break;
              }
            }
            if (!bidirectional) {
              logger.warning(
                  format(
                      " Check the bond between %d and %d in %s.",
                      a1, a2, activeMolecularAssembly.getFile().getName()));
              return false;
            }
            Atom atom1 = atomList.get(a1 - 1);
            Atom atom2 = atomList.get(a2 - 1);
            if (atom1 == null || atom2 == null) {
              logger.warning(
                  format(
                      " Check the bond between %d and %d in %s.",
                      a1, a2, activeMolecularAssembly.getFile().getName()));
              return false;
            }
            Bond bond = new Bond(atom1, atom2);
            c[0] = atom1.getAtomType().atomClass;
            c[1] = atom2.getAtomType().atomClass;
            String key = BondType.sortKey(c);
            BondType bondType = forceField.getBondType(key);
            if (bondType == null) {
              logNoBondType(atom1, atom2, key);
            } else {
              bond.setBondType(bondType);
            }
            bondList.add(bond);
          }
        }
      }
      return true;
    } catch (IOException e) {
      logger.severe(e.toString());
    }
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public boolean readNext() {
    return readNext(false);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Reads the next snap-shot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  @Override
  public boolean readNext(boolean resetPosition) {
    return readNext(resetPosition, true);
  }

  /**
   * Reads the next snap-shot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  public boolean readNext(boolean resetPosition, boolean print) {
    try {
      String data;
      Atom[] atoms = activeMolecularAssembly.getAtomArray();
      int nSystem = atoms.length;

      if (bufferedReader == null && !resetPosition) {
        bufferedReader = new BufferedReader(new FileReader(currentFile));
        // Read past the first N + 1 lines that begin with an integer.
        for (int i = 0; i < nSystem + 1; i++) {
          data = bufferedReader.readLine();
          while (!firstTokenIsInteger(data)) {
            data = bufferedReader.readLine();
          }
        }
        snapShot = 1;
      } else if (resetPosition){
        //Reset the reader to the beginning of the file. Do not skip reading the first entry if resetPostion is true.
        bufferedReader = new BufferedReader(new FileReader(currentFile));
        snapShot = 0;
      }

      snapShot++;

      data = bufferedReader.readLine();
      // Read past blank lines
      while (data != null && data.trim().equals("")) {
        data = bufferedReader.readLine();
      }
      if (data == null) {
        return false;
      }

      if (print) {
        logger.info(format("\n Attempting to read snapshot %d.", snapShot));
      }
      try {
        int nArchive = parseInt(data.trim().split(" +")[0]);
        if (nArchive != nSystem) {
          String message =
              format("Number of atoms mismatch (Archive: %d, System: %d).", nArchive, nSystem);
          if (dieOnMissingAtom) {
            logger.severe(message);
          }
          logger.warning(message);
          return false;
        }
      } catch (NumberFormatException e) {
        logger.warning(e.toString());
        return false;
      }

      remarkLine = data;

      // The header line is reasonable. Check for periodic box dimensions.
      bufferedReader.mark(10000);
      data = bufferedReader.readLine();
      if (!readPBC(data, activeMolecularAssembly)) {
        bufferedReader.reset();
      }

      for (int i = 0; i < nSystem; i++) {
        data = bufferedReader.readLine();
        // Read past blank lines
        while (data != null && data.trim().equals("")) {
          data = bufferedReader.readLine();
        }
        String[] tokens = data.trim().split(" +");
        if (tokens.length < 6) {
          String message = format("Check atom %d in %s.", (i + 1), currentFile.getName());
          logger.warning(message);
          return false;
        }
        double x = parseDouble(tokens[2]);
        double y = parseDouble(tokens[3]);
        double z = parseDouble(tokens[4]);
        int xyzIndex = atoms[i].getIndex();
        if (xyzIndex != i + 1) {
          String message =
              format(
                  "Archive atom index %d being read onto system atom index %d.", i + 1, xyzIndex);
          logger.warning(message);
        }
        atoms[i].moveTo(x, y, z);
      }
      return true;
    } catch (FileNotFoundException e) {
      String message = format("Exception opening file %s.", currentFile);
      logger.log(Level.WARNING, message, e);
    } catch (IOException e) {
      String message = format("Exception reading from file %s.", currentFile);
      logger.log(Level.WARNING, message, e);
    }
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public boolean writeFile(File saveFile, boolean append, String[] extraLines) {
    if (saveFile == null) {
      return false;
    }

    File newFile = saveFile;
    if (!append) {
      newFile = version(saveFile);
    }
    activeMolecularAssembly.setFile(newFile);
    activeMolecularAssembly.setName(newFile.getName());

    try (FileWriter fw = new FileWriter(newFile, append && newFile.exists());
        BufferedWriter bw = new BufferedWriter(fw)) {
      // XYZ File First Line
      int numberOfAtoms = activeMolecularAssembly.getAtomList().size();
      StringBuilder sb =
          new StringBuilder(format("%7d  %s", numberOfAtoms, activeMolecularAssembly.toString()));
      if (extraLines != null) {
        for (String line : extraLines) {
          line = line.replaceAll("\n", " ");
          sb.append(" ").append(line);
        }
      }
      String output = sb.append("\n").toString();
      bw.write(output);

      Crystal crystal = activeMolecularAssembly.getCrystal();
      if (!crystal.aperiodic()) {
        Crystal uc = crystal.getUnitCell();
        String params =
            format(
                "%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f\n",
                uc.a, uc.b, uc.c, uc.alpha, uc.beta, uc.gamma);
        bw.write(params);
      }

      Atom a2;
      StringBuilder line;
      StringBuilder[] lines = new StringBuilder[numberOfAtoms];
      // XYZ File Atom Lines
      List<Atom> atoms = activeMolecularAssembly.getAtomList();
      Vector3d offset = activeMolecularAssembly.getOffset();
      for (Atom a : atoms) {
        if (vdwH) {
          line =
              new StringBuilder(
                  format(
                      "%7d %3s%14.8f%14.8f%14.8f%6d",
                      a.getIndex(),
                      a.getAtomType().name,
                      a.getRedX() - offset.x,
                      a.getRedY() - offset.y,
                      a.getRedZ() - offset.z,
                      a.getType()));
        } else {
          line =
              new StringBuilder(
                  format(
                      "%7d %3s%14.8f%14.8f%14.8f%6d",
                      a.getIndex(),
                      a.getAtomType().name,
                      a.getX() - offset.x,
                      a.getY() - offset.y,
                      a.getZ() - offset.z,
                      a.getType()));
        }
        for (Bond b : a.getBonds()) {
          a2 = b.get1_2(a);
          line.append(format("%8d", a2.getIndex()));
        }
        lines[a.getIndex() - 1] = line.append("\n");
      }
      try {
        for (int i = 0; i < numberOfAtoms; i++) {
          bw.write(lines[i].toString());
        }
      } catch (IOException e) {
        String message =
            format(
                " Their was an unexpected error writing to %s.",
                getActiveMolecularSystem().toString());
        logger.log(Level.WARNING, message, e);
        return false;
      }
    } catch (IOException e) {
      String message =
          format(
              " Their was an unexpected error writing to %s.",
              getActiveMolecularSystem().toString());
      logger.log(Level.WARNING, message, e);
      return false;
    }
    return true;
  }

  /**
   * writeFileAsP1
   *
   * @param saveFile a {@link java.io.File} object.
   * @param append a boolean.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @return a boolean.
   */
  public boolean writeFileAsP1(File saveFile, boolean append, Crystal crystal) {
    return writeFileAsP1(saveFile, append, crystal, null);
  }

  /**
   * writeFileAsP1
   *
   * @param saveFile a {@link java.io.File} object.
   * @param append a boolean.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param extraLines Additional lines to print in the header.
   * @return a boolean.
   */
  public boolean writeFileAsP1(
      File saveFile, boolean append, Crystal crystal, String[] extraLines) {
    if (saveFile == null) {
      return false;
    }

    File newFile = saveFile;
    if (!append) {
      newFile = version(saveFile);
    }
    activeMolecularAssembly.setFile(newFile);
    activeMolecularAssembly.setName(newFile.getName());

    try (FileWriter fw = new FileWriter(newFile, append && newFile.exists());
        BufferedWriter bw = new BufferedWriter(fw)) {
      int nSymm = crystal.spaceGroup.symOps.size();
      // XYZ File First Line
      int numberOfAtoms = activeMolecularAssembly.getAtomList().size() * nSymm;
      StringBuilder sb =
          new StringBuilder(format("%7d  %s", numberOfAtoms, activeMolecularAssembly.toString()));
      if (extraLines != null) {
        for (String line : extraLines) {
          line = line.replaceAll("\n", " ");
          sb.append(" ").append(line);
        }
      }
      String output = sb.append("\n").toString();
      bw.write(output);

      if (!crystal.aperiodic()) {
        Crystal uc = crystal.getUnitCell();
        String params =
            format(
                "%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f\n",
                uc.a, uc.b, uc.c, uc.alpha, uc.beta, uc.gamma);
        bw.write(params);
      }

      Atom a2;
      StringBuilder line;
      StringBuilder[] lines = new StringBuilder[numberOfAtoms];
      // XYZ File Atom Lines
      Atom[] atoms = activeMolecularAssembly.getAtomArray();
      double[] xyz = new double[3];
      for (int iSym = 0; iSym < nSymm; iSym++) {
        SymOp symOp = crystal.spaceGroup.getSymOp(iSym);
        int indexOffset = iSym * atoms.length;
        for (Atom a : atoms) {
          int index = a.getIndex() + indexOffset;
          String id = a.getAtomType().name;
          if (vdwH) {
            a.getRedXYZ(xyz);
          } else {
            xyz[0] = a.getX();
            xyz[1] = a.getY();
            xyz[2] = a.getZ();
          }
          crystal.applySymOp(xyz, xyz, symOp);
          int type = a.getType();
          line =
              new StringBuilder(
                  format("%7d %3s%14.8f%14.8f%14.8f%6d", index, id, xyz[0], xyz[1], xyz[2], type));
          for (Bond b : a.getBonds()) {
            a2 = b.get1_2(a);
            line.append(format("%8d", a2.getIndex() + indexOffset));
          }
          lines[index - 1] = line.append("\n");
        }
      }
      try {
        for (int i = 0; i < numberOfAtoms; i++) {
          bw.write(lines[i].toString());
        }
      } catch (IOException e) {
        String message =
            format(
                " Their was an unexpected error writing to %s.",
                getActiveMolecularSystem().toString());
        logger.log(Level.WARNING, message, e);
        return false;
      }
    } catch (IOException e) {
      String message =
          format(
              " Their was an unexpected error writing to %s.",
              getActiveMolecularSystem().toString());
      logger.log(Level.WARNING, message, e);
      return false;
    }
    return true;
  }
}
