//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

package ffx.potential.groovy

import ffx.numerics.math.DoubleMath
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Atom
import ffx.potential.bonded.NucleicAcidUtils
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import ffx.crystal.Crystal
import picocli.CommandLine.Option
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.toDegrees

@Command(description = "Calculates nucleic acid torsions as well as information regarding the sugar pucker.", name = "NucleicAcidAnalysis")
class NucleicAcidAnalysis extends PotentialScript {

  /**
   * --xo or --chiOnly Print all torsions and information.
   */
  @Option(names = ["--at", "--allTorsions"], paramLabel = "false", defaultValue = "false",
      description = "Print all torsions and information.")
  private boolean allTorsions = false

  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ or ARC format.')
  private List<String> filenames = null

  private List<Residue> residues

  NucleicAcidAnalysis() {
    this(new groovy.lang.Binding())
  }

  NucleicAcidAnalysis(Binding binding) {
    super(binding)
  }

  @Override
  NucleicAcidAnalysis run() {
    if (!init()) {
      return null
    }

    activeAssembly = getActiveAssembly(filenames[0])
    if (activeAssembly == null) {
      logger.info(helpString())
      return null
    }

    ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x)

    // Create NucleicAcidUtils instance once
    def naUtils = new NucleicAcidUtils()

    // Filter for nucleic acid residues
    residues = activeAssembly.getResidueList().findAll { residue ->
      def name = residue.name.toUpperCase()
      name in ["DA", "DC", "DG", "DT", "DU", "A", "C", "G", "T", "U",
               "DAD", "DCY", "DGU", "DTY", "URA"]
    }

    SystemFilter systemFilter = potentialFunctions.getFilter()
    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
      int numModels = systemFilter.countNumModels()

      if (numModels == 1 && allTorsions) {
        // Single frame - output all information
        printHeader()
        analyzeAndPrintResidues()
      } else {
        // Multiple frames (ARC file) - output only pseudorotation and sugar pucker per frame
        int frameNumber = 1
        println("\nFrame $frameNumber:")
        println("Residue    Name     χ          P         Sugar pucker")
        println("-----------------------------------------------------")

        for (Residue residue : residues) {
          double chi = getDihedral(residue, "O4'", "C1'", "N9", "C4") ?: getDihedral(residue, "O4'", "C1'", "N1", "C2")
          double P = calculateP(residue)
          String sugarPucker = determineSugarPucker(P)

          println(format("%-10s %-8s %-10s %-10s %-14s",
              residue.getResidueNumber(),
              residue.name,
              formatValue(chi),
              formatValue(P),
              sugarPucker
          ))
        }
        frameNumber++
        while (systemFilter.readNext()) {
          Crystal crystal = activeAssembly.getCrystal()
          forceFieldEnergy.setCrystal(crystal)
          forceFieldEnergy.getCoordinates(x)
          forceFieldEnergy.energy(x)

          // Filter residues again in case the model changed
          residues = activeAssembly.getResidueList().findAll { residue ->
            def name = residue.name.toUpperCase()
            name in ["DA", "DC", "DG", "DT", "DU", "A", "C", "G", "T", "U",
                     "DAD", "DCY", "DGU", "DTY", "URA"]
          }

          println("\nFrame $frameNumber:")
          println("Residue    Name     χ          P         Sugar pucker")
          println("-----------------------------------------------------")

          for (Residue residue : residues) {
            double chi = getDihedral(residue, "O4'", "C1'", "N9", "C4") ?: getDihedral(residue, "O4'", "C1'", "N1", "C2")
            double P = calculateP(residue)
            String sugarPucker = determineSugarPucker(P)

            println(format("%-10s %-8s %-10s %-10s %-14s",
                residue.getResidueNumber(),
                residue.name,
                formatValue(chi),
                formatValue(P),
                sugarPucker
            ))
          }
          frameNumber++
        }
      }
    }

    return this
  }

  private void printHeader() {
    println("Residue    Name      V0         V1         V2         V3         V4         P          νmax       χ          γ          α          β          δ          ε          ζ          Sugar pucker                  Stage")
    println("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
  }

  private void analyzeAndPrintResidues() {
    for (Residue residue : residues) {
      def v0 = getDihedral(residue, "C2'", "C1'", "O4'", "C4'")
      def v1 = getDihedral(residue, "O4'", "C1'", "C2'", "C3'")
      def v2 = getDihedral(residue, "C1'", "C2'", "C3'", "C4'")
      def v3 = getDihedral(residue, "O4'", "C4'", "C3'", "C2'")
      def v4 = getDihedral(residue, "C1'", "O4'", "C4'", "C3'")
      def chi = getDihedral(residue, "O4'", "C1'", "N9", "C4") ?: getDihedral(residue, "O4'", "C1'", "N1", "C2")
      def gamma = getDihedral(residue, "O5'", "C5'", "C4'", "C3'")
      def alpha = getDihedral(residue, "O3'", "P", "O5'", "C5'")
      def beta = getDihedral(residue, "P", "O5'", "C5'", "C4'")
      def delta = getDihedral(residue, "C5'", "C4'", "C3'", "O3'")
      def epsilon = getDihedral(residue, "C4'", "C3'", "O3'", "P")
      def zeta = getDihedral(residue, "C3'", "O3'", "P", "O5'")

      Double P = (v0 != null && v1 != null && v3 != null && v4 != null && v2 != null) ? calculateP(v0, v1, v2, v3, v4) : null
      double nuMax = (v2 != null && P != null) ? Math.abs(v2 / Math.cos(Math.toRadians(P))) : null
      String type = determineType(residue, P)
      String stage = determineStage(delta, chi, P)

      println(format("%-10s %-8s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-14s %-18s",
          residue.getResidueNumber(),
          residue.name,
          formatValue(v0), formatValue(v1), formatValue(v2),
          formatValue(v3), formatValue(v4),
          formatValue(P), formatValue(nuMax),
          formatValue(chi), formatValue(gamma),
          formatValue(alpha), formatValue(beta),
          formatValue(delta), formatValue(epsilon),
          formatValue(zeta),
          type,
          stage
      ))
    }
  }

  private static Double getDihedral(Residue residue, String atom1, String atom2, String atom3, String atom4) {
    try {
      Atom a1 = residue.getAtomByName(atom1, true)
      Atom a2 = residue.getAtomByName(atom2, true)
      Atom a3 = residue.getAtomByName(atom3, true)
      Atom a4 = residue.getAtomByName(atom4, true)

      if (a1 != null && a2 != null && a3 != null && a4 != null) {
        return toDegrees(DoubleMath.dihedralAngle(a1.getXYZ(null), a2.getXYZ(null), a3.getXYZ(null), a4.getXYZ(null)))
      }
    } catch (Exception ignored) {
      logger.warning("Could not calculate dihedral for ${residue} atoms $atom1-$atom2-$atom3-$atom4")
    }
    return null
  }

  private static String formatValue(Double value) {
    return value != null ? format("%.2f", value) : "N/A"
  }

  private static double calculateP(Residue residue) {
    def v0 = getDihedral(residue, "C2'", "C1'", "O4'", "C4'")
    def v1 = getDihedral(residue, "O4'", "C1'", "C2'", "C3'")
    def v2 = getDihedral(residue, "C1'", "C2'", "C3'", "C4'")
    def v3 = getDihedral(residue, "O4'", "C4'", "C3'", "C2'")
    def v4 = getDihedral(residue, "C1'", "O4'", "C4'", "C3'")

    double P = (v0 != null && v1 != null && v3 != null && v4 != null && v2 != null) ?
        calculateP(v0, v1, v2, v3, v4) : null

    return P
  }

  private static double calculateP(double v0, double v1, double v2, double v3, double v4) {
    try {
      double sin36 = Math.sin(Math.toRadians(36))
      double sin72 = Math.sin(Math.toRadians(72))
      double denominator = 2 * v2 * (sin36 + sin72)
      if (Math.abs(denominator) < 1e-10) {
        return Double.NaN
      }
      double p = ((v4 - v0) - (v3 - v1)) / denominator

      double P;
      if (v2 < 0) {
        P = Math.toDegrees(Math.atan(p)) + 180.0
      } else if (p < 0) {
        P = Math.toDegrees(Math.atan(p)) + 360.0
      } else {
        P = Math.toDegrees(Math.atan(p))
      }
      return P
    } catch (Exception ignored) {
      return Double.NaN
    }
  }

  private static String determineType(Residue residue, Double P) {
    if (P == null) return "Unknown"

    String base = switch (residue.name) {
      case "DAD" -> "Ade"
      case "DGU" -> "Gua"
      case "DCY" -> "Cyt"
      case "DTY" -> "Thy"
      case "URA" -> "Ura"
      default -> "Unknown"
    }

    String sugarPucker = determineSugarPucker(P)
    return base + ", " + sugarPucker
  }

  private static String determineSugarPucker(Double P) {
    if (P == null) return "Unknown"

    if (P >= 0 && P < 36) {
      return "C3'-endo"
    } else if (P >= 36 && P < 72) {
      return "C4'-endo"
    } else if (P >= 72 && P < 108) {
      return "O4'-endo"
    } else if (P >= 108 && P < 144) {
      return "C1'-exo"
    } else if (P >= 144 && P < 180) {
      return "C2'-endo"
    } else if (P >= 180 && P < 216) {
      return "C3'-exo"
    } else if (P >= 216 && P < 252) {
      return "C4'-exo"
    } else if (P >= 252 && P < 288) {
      return "O4'-exo"
    } else if (P >= 288 && P < 324) {
      return "C1'-endo"
    } else if (P >= 324 || P < 0) {
      return "C2'-endo"
    }
    return "Unknown"
  }

  private static String determineStage(Double delta, Double chi, Double P) {
    if (delta == null || chi == null || P == null) return "Unknown"

    if (delta > 120.0 && chi > -160) {
      return "Stage 1: Slide and Roll Change"
    } else if (delta < 120.0 && P >= 0 && P < 180.0) {
      return "Stage 2: Backbone and Sugar Pucker Change"
    } else if (delta < 120.0 && P >= 180.0 && P < 360.0) {
      return "Stage 3: Roll and Inclination to Finish Transition"
    }
    return "Unknown"
  }
}