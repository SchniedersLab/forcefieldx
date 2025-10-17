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
package ffx.potential.commands;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Residue;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.utils.GetProteinFeatures;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Create a Feature Map for a given protein structure.
 *
 * Usage:
 *   ffxc FeatureMap [options] &lt;pdb/xyz&gt; &lt;variants.csv&gt; [ddgFile]
 */
@Command(name = "FeatureMap", description = " Create a Feature Map for a given protein structure")
public class FeatureMap extends PotentialCommand {

  @Option(names = {"-d", "--delimiter"}, paramLabel = ",",
      description = "Delimiter of input variant list file")
  private String delimiter = ",";

  @Option(names = {"--iP", "--includePolarity"}, defaultValue = "false",
      description = "Include polarity change in feature map.")
  private boolean includePolarity = false;

  @Option(names = {"--iA", "--includeAcidity"}, defaultValue = "false",
      description = "Include acidity change in feature map.")
  private boolean includeAcidity = false;

  @Option(names = {"--iAng", "--includeAngles"}, defaultValue = "false",
      description = "Include ramachandran angles in feature map.")
  private boolean includeAngles = false;

  @Option(names = {"--iS", "--includeStructure"}, defaultValue = "false",
      description = "Include secondary structure annotations in feature map.")
  private boolean includeStructure = false;

  @Option(names = {"--iPPI", "--includePPI"}, defaultValue = "false",
      description = "Mark residue as being part of an interaction interface.")
  private boolean includePPI = false;

  @Option(names = {"--rR", "--reRun"}, defaultValue = "false",
      description = "Reading in a CSV a second time. Specifically for gathering ppi information")
  private boolean rerun = false;

  @Option(names = {"--ddG", "--ddgFile"}, paramLabel = "file",
      description = "Read in the folding free energy difference file")
  private String ddgFile;

  @Option(names = {"--mI", "--multiple isomers"}, defaultValue = "false",
      description = "Set this flag if the variant list contains variants from multiple isomers. Isomer should be in name of pdb file")
  private boolean multipleIsoforms = false;

  /** The final argument(s) should be one or more filenames. */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in PDB or XYZ format, variant list file, and a free energy file")
  private List<String> filenames = null;

  private List<Residue> residues;

  public FeatureMap() { super(); }
  public FeatureMap(FFXBinding binding) { super(binding); }
  public FeatureMap(String[] args) { super(args); }

  @Override
  public FeatureMap run() {
    // Init the context and bind variables.
    if (!init()) {
      return null;
    }

    // Enable GK surface area for confidence calculation as in Groovy script.
    System.setProperty("gkterm", "true");
    System.setProperty("cavmodel", "cav");
    System.setProperty("surface-tension", "1.0");

    // Load the MolecularAssembly.
    String structureFile = filenames != null && !filenames.isEmpty() ? filenames.get(0) : null;
    activeAssembly = getActiveAssembly(structureFile);
    if (activeAssembly == null) {
      logger.info(helpString());
      return null;
    }

    ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy();
    int nVars = forceFieldEnergy.getNumberOfVariables();
    double[] x = new double[nVars];
    forceFieldEnergy.getCoordinates(x);
    forceFieldEnergy.energy(x);

    residues = activeAssembly.getResidueList();
    GetProteinFeatures getProteinFeatures = new GetProteinFeatures(activeAssembly);

    // Handle multiple isoforms by parsing the PDB basename.
    String fileIsoform = null;
    if (multipleIsoforms) {
      String baseName = new File(structureFile).getName();
      baseName = baseName.replaceFirst("\\.pdb$", "");
      if (baseName.toUpperCase().contains("ENS")) {
        String[] geneSplit = baseName.split("_");
        if (geneSplit.length > 1) {
          fileIsoform = geneSplit[1];
        }
      } else {
        String[] geneSplit = baseName.split("_");
        int index = -1;
        for (int i = 0; i < geneSplit.length; i++) {
          if (geneSplit[i].equalsIgnoreCase("NP")) {
            index = i;
          }
        }
        if (index >= 0 && index + 1 < geneSplit.length) {
          fileIsoform = (geneSplit[index] + '_' + geneSplit[index + 1]).toUpperCase();
        }
      }
    }

    // Build per-residue feature list (surface area, normalized SA, confidence, optional others).
    List<String[]> featureList = new ArrayList<>();
    for (Residue residue : residues) {
      double residueSurfaceArea = forceFieldEnergy.getGK().getSurfaceAreaRegion().getResidueSurfaceArea(residue);
      featureList.add(getProteinFeatures.saveFeatures(residue, residueSurfaceArea, includeAngles, includeStructure, includePPI));
    }

    // Optional ddG input parsing for ddgun output.
    List<String> ddgunLines = new ArrayList<>();
    List<String> npChanges = new ArrayList<>();
    List<Double[]> ddGun = new ArrayList<>();
    List<String[]> polarityAndAcidityChange = new ArrayList<>();
    if (ddgFile != null) {
      try (BufferedReader txtReader = new BufferedReader(new FileReader(new File(ddgFile)))) {
        String line;
        while ((line = txtReader.readLine()) != null) {
          if (line.contains(".pdb")) {
            ddgunLines.add(line);
          }
        }
      } catch (IOException e) {
        logger.warning(" Error reading ddG file: " + e);
      }
      ddGun = getProteinFeatures.getDDGunValues(ddgunLines);
      npChanges = getProteinFeatures.ddgunToNPChange(ddgunLines);
      if (includePolarity || includeAcidity) {
        polarityAndAcidityChange = getProteinFeatures.getPolarityAndAcidityChange(npChanges, includePolarity, includeAcidity);
      }
    }

    // Write updated CSV alongside the input variants CSV (filenames[1]).
    BufferedReader br = null;
    BufferedWriter bw = null;
    try {
      String variantsCSV = filenames.get(1);
      File inputCSVFile = new File(variantsCSV);
      String inputCSVPath = inputCSVFile.getParent();
      if (inputCSVPath == null) inputCSVPath = ".";
      String newCSVFileName = "update_" + inputCSVFile.getName();
      File updatedFile = new File(inputCSVPath, newCSVFileName);

      br = new BufferedReader(new InputStreamReader(new FileInputStream(inputCSVFile)));
      bw = new BufferedWriter(new FileWriter(updatedFile, true));

      int i = 0;
      int npIndex = -1;
      int isoformIndex = -1;
      for (String line = br.readLine(); line != null; line = br.readLine(), i++) {
        StringBuilder newCSVLine = new StringBuilder();
        if (i == 0) {
          if (!rerun) {
            if (updatedFile.length() == 0) {
              newCSVLine.append(line).append(delimiter)
                  .append("Surface Area").append(delimiter)
                  .append("Normalized SA").append(delimiter)
                  .append("Confidence Score");
              if (ddgFile != null) {
                newCSVLine.append(delimiter).append("ddG").append(delimiter).append("|ddG|");
              }
              if (includeAcidity) {
                newCSVLine.append(delimiter).append("Acidity Change");
              }
              if (includePolarity) {
                newCSVLine.append(delimiter).append("Polarity Change");
              }
              if (includeAngles) {
                newCSVLine.append(delimiter).append("Phi").append(delimiter).append("Psi").append(delimiter).append("Omega");
              }
              if (includeStructure) {
                newCSVLine.append(delimiter).append("Secondary Structure Annotation");
              }
              if (includePPI) {
                newCSVLine.append(delimiter).append("Interface Residue");
              }
              bw.write(newCSVLine.toString());
            } else {
              bw.write(line);
              bw.newLine();
            }
          }
        } else {
          String[] splits = line.split(delimiter, -1);
          if (i == 1) {
            for (int j = 0; j < splits.length; j++) {
              if (splits[j].contains("p.")) {
                npIndex = j;
              }
              if (multipleIsoforms) {
                String sj = splits[j];
                String sju = sj.toUpperCase();
                if (sju.contains("NP") || sju.contains("XP") || sju.contains("ENS")) {
                  isoformIndex = j;
                }
              }
            }
          }
          int length = splits.length;
          String npChange = npIndex >= 0 && npIndex < splits.length ? splits[npIndex] : "";

          String[] feat = featureList.isEmpty() ? new String[0] : new String[featureList.get(0).length];
          Arrays.fill(feat, null);
          String[] ddG = new String[]{"null", "null"};
          String[] pA = new String[]{"null", "null"};

          if (npChange.contains("del")) {
            // leave feat as nulls and ddG/pA as "null" strings
          } else if (!featureList.isEmpty()) {
            int position = -1;
            int pIndex = npChange.indexOf('p');
            if (pIndex >= 0) {
              String proteinChange = npChange.substring(pIndex);
              String[] parts = proteinChange.split("p\\.");
              if (parts.length > 1) {
                String sub = parts[1];
                if (sub.length() >= 6) {
                  try {
                    position = Integer.parseInt(sub.substring(3, sub.length() - 3));
                  } catch (Exception ignored) {
                  }
                }
              }
              if (position > 0 && position <= residues.size()) {
                if (ddgFile != null) {
                  int idx = npChanges.indexOf(proteinChange);
                  if (idx != -1) {
                    Double[] d = ddGun.get(idx);
                    if (d != null && d.length >= 2) {
                      ddG = new String[]{String.valueOf(d[0]), String.valueOf(d[1])};
                    }
                    if (includeAcidity || includePolarity) {
                      String[] pa = polarityAndAcidityChange.get(idx);
                      if (pa != null && pa.length >= 2) {
                        pA = new String[]{pa[0], pa[1]};
                      }
                    }
                  }
                }
                feat = featureList.get(position - 1);
              }
            }
          }

          String isoform = null;
          if (multipleIsoforms && isoformIndex >= 0 && isoformIndex < splits.length) {
            String isoStr = splits[isoformIndex];
            isoform = isoStr.contains(":p.") ? isoStr.split(":")[0] : isoStr;
          }

          if (splits.length == length) {
            if (multipleIsoforms) {
              if (fileIsoform != null) {
                String isofU = (isoform == null ? "" : isoform).toUpperCase();
                if (!isofU.contains(fileIsoform.toUpperCase())) {
                  continue;
                }
              }
            }
            if (rerun) {
              newCSVLine.append(line);
            } else {
              if (ddgFile != null) {
                newCSVLine.append(line).append(delimiter)
                    .append(safe(feat, 0)).append(delimiter)
                    .append(safe(feat, 1)).append(delimiter)
                    .append(safe(feat, 2)).append(delimiter)
                    .append(ddG[0]).append(delimiter).append(ddG[1]);
              } else {
                newCSVLine.append(line).append(delimiter)
                    .append(safe(feat, 0)).append(delimiter)
                    .append(safe(feat, 1)).append(delimiter)
                    .append(safe(feat, 2));
              }
            }

            if (includeAcidity && !rerun) {
              newCSVLine.append(delimiter).append(pA[0]);
            }
            if (includePolarity && !rerun) {
              newCSVLine.append(delimiter).append(pA[1]);
            }
            if (includeAngles && !rerun) {
              newCSVLine.append(delimiter).append(safe(feat, 3))
                  .append(delimiter).append(safe(feat, 4))
                  .append(delimiter).append(safe(feat, 5));
              if (includeStructure) {
                newCSVLine.append(delimiter).append(safe(feat, 6));
              }
              if (includePPI) {
                newCSVLine.append(delimiter).append(safe(feat, 7));
              }
            } else if (includeStructure && !rerun) {
              newCSVLine.append(delimiter).append(safe(feat, 3));
              if (includePPI) {
                newCSVLine.append(delimiter).append(safe(feat, 4));
              }
            } else if (includePPI) {
              newCSVLine.append(delimiter).append(safe(feat, 3));
            }
            bw.newLine();
            bw.write(newCSVLine.toString());
          }
        }
      }
    } catch (Exception e) {
      logger.warning(" Error updating CSV: " + e);
    } finally {
      try {
        if (br != null) br.close();
      } catch (IOException ignored) {}
      try {
        if (bw != null) bw.close();
      } catch (IOException ignored) {}
    }

    logger.info(" Wrote variants with feature data to update_" + (filenames != null && filenames.size() > 1 ? new File(filenames.get(1)).getName() : "variants.csv"));
    return this;
  }

  private static String safe(String[] arr, int idx) {
    if (arr == null || idx < 0 || idx >= arr.length) return "";
    return arr[idx] == null ? "" : arr[idx];
  }
}
