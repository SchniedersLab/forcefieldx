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
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.math3.util.FastMath.toDegrees
import static org.apache.commons.math3.util.FastMath.abs
import static org.apache.commons.math3.util.FastMath.atan
import static org.apache.commons.math3.util.FastMath.sin
import static org.apache.commons.math3.util.FastMath.cos

@Command(description = "Nucleic Acid Analysis", name = "NucleicAcidAnalysis")
class NucleicAcidAnalysis extends PotentialScript {

    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private List<String> filenames = null

    private List<Residue> residues

    NucleicAcidAnalysis() {
        this(new Binding())
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

        residues = activeAssembly.getResidueList()
        println("Residue    Name      V0         V1         V2         V3         V4         P          νmax       χ          γ          α          β          δ          ε          ζ          Sugar pucker                  Stage")
        println("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

        for (Residue residue : residues) {
            def v0 = getDihedral(residue, "C2'", "C1'", "O4'", "C4'")
            def v1 = getDihedral(residue, "O4'", "C1'", "C2'", "C3'")
            def v2 = getDihedral(residue, "C1'", "C2'", "C3'", "C4'")
            def v3 = getDihedral(residue, "O4'", "C4'", "C3'", "C2'")
            def v4 = getDihedral(residue, "C1'", "O4'", "C4'", "C3'")
            def chi = getDihedral(residue, "O4'", "C1'", "N9", "C4") ?: getDihedral(residue, "O4'", "C1'", "N1", "C2") // The glycosidic angle is influenced by the δ angle and plays a role in the overall structural transition.
            // In B-form DNA, χ is typically around –160° to –120° (or 200°–240°).
            // In A-form DNA, χ is typically around –60° to –90° (or 270°–300°).

            def gamma = getDihedral(residue, "O5'", "C5'", "C4'", "C3'")

            // Calculate backbone torsion angles(specific type of dihedral angle)
            def alpha = getDihedral(residue, "O3'", "P", "O5'", "C5'") // Rotation around the P–O5' bond.
            def beta = getDihedral(residue, "P", "O5'", "C5'", "C4'") // Rotation around the O5'–C5' bond.
            def delta = getDihedral(residue, "C5'", "C4'", "C3'", "O3'") // Rotation around the C4'–C3' bond. A drop in this angle triggers the rotation of the glycosidic bond and subsequent sugar pucker conversion.
            // In B-form DNA, δ is typically around 140°–150°.
            // In A-form DNA, δ is typically around 80°–90°.

            def epsilon = getDihedral(residue, "C4'", "C3'", "O3'", "P") // Rotation around the C3'–O3' bond.
            def zeta = getDihedral(residue, "C3'", "O3'", "P", "O5'") // Rotation around the O3'–P bond.

            // Calculate pseudorotational phase angle
            Double P = (v0 != null && v1 != null && v3 != null && v4 != null && v2 != null) ? calculateP(v0, v1, v2, v3, v4) : null
            // C2′-endo (B-form): Around 162°.
            // C3′-endo (A-form): Around 18°.


            // Calculate νmax (puckering amplitude)
            double nuMax = (v2 != null && P != null) ? Math.abs(v2 / Math.cos(Math.toRadians(P))) : null

            // Determine the type
            String type = determineType(residue, P)

            // Determine the current stage of the B-to-A transition
            String stage = determineStage(delta, chi, P)

            // Align the format string with the headers
            println(String.format("%-10s %-8s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-14s %-18s",
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

        return this
    }

    private Double getDihedral(Residue residue, String atom1, String atom2, String atom3, String atom4) {
        Atom a1 = residue.getAtomByName(atom1, true)
        Atom a2 = residue.getAtomByName(atom2, true)
        Atom a3 = residue.getAtomByName(atom3, true)
        Atom a4 = residue.getAtomByName(atom4, true)

        if (a1 != null && a2 != null && a3 != null && a4 != null) {
            return toDegrees(DoubleMath.dihedralAngle(a1.getXYZ(null), a2.getXYZ(null), a3.getXYZ(null), a4.getXYZ(null)))
        }
        return null
    }

    private String formatValue(Double value) {
        return value != null ? String.format("%.2f", value) : "N/A"
    }

    private static double calculateP(double v0, double v1, double v2, double v3, double v4) {
        // Calculate p
        double sin36 = Math.sin(Math.toRadians(36))
        double sin72 = Math.sin(Math.toRadians(72))
        double denominator = 2 * v2 * (sin36 + sin72)
        double p = ((v4 - v0) - (v3 - v1)) / denominator

        // Calculate P
        double P;
        if (v2 < 0) {
            P = Math.toDegrees(Math.atan(p)) + 180.0
        } else if (p < 0) {
            P = Math.toDegrees(Math.atan(p)) + 360.0
        } else {
            P = Math.toDegrees(Math.atan(p))
        }
        return P
    }

    private String determineType(Residue residue, Double P) {
        if (P == null) return "Unknown"

        // Determine the base
        String base = switch (residue.name) {
            case "DAD" -> "Ade"
            case "DGU" -> "Gua"
            case "DCY" -> "Cyt"
            case "DTY" -> "Thy"
            case "URA" -> "Ura"
            default -> "Unknown"
        }
        if (!["DAD", "DGU", "DCY", "DTY"].contains(residue.name)) {
            println("DEBUG: Unknown residue name '${residue.name}'")
        }

        // Determine sugar pucker conformation
        String sugarPucker = "Unknown"
        if (P >= 0 && P < 36) {
            sugarPucker = "C3'-endo"
        } else if (P >= 36 && P < 72) {
            sugarPucker = "C4'-endo"
        } else if (P >= 72 && P < 108) {
            sugarPucker = "O4'-endo"
        } else if (P >= 108 && P < 144) {
            sugarPucker = "C1'-exo"
        } else if (P >= 144 && P < 180) {
            sugarPucker = "C2'-endo"
        } else if (P >= 180 && P < 216) {
            sugarPucker = "C3'-exo"
        } else if (P >= 216 && P < 252) {
            sugarPucker = "C4'-exo"
        } else if (P >= 252 && P < 288) {
            sugarPucker = "O4'-exo"
        } else if (P >= 288 && P < 324) {
            sugarPucker = "C1'-endo"
        } else if (P >= 324 || P < 0) {
            sugarPucker = "C2'-endo"
        }

        return base + ", " + sugarPucker
    }

    private String determineStage(Double delta, Double chi, Double P) {
        /*
        * Uses the δ angle, χ angle, and pseudorotation angle (P) to classify the current state into one of the three stages of the B-to-A transition.
        */

        if (delta == null || chi == null || P == null) return "Unknown"

        // Stage 1: Slide and Roll change (initial base pair movement)
        // The value 120.0 was chosen as a midpoint to distinguish between B-form-like (higher δ) and A-form-like (lower δ) states.
        // The value 200.0(-160) was chosen to distinguish between B-form-like (higher χ) and A-form-like (lower χ) states.
        if (delta > 120.0 && chi > -160) {
            return "Stage 1: Slide and Roll Change"
        }
        // Stage 2: Backbone and sugar pucker change
        else if (delta < 120.0 && P >= 0 && P < 180.0) {
            // This range includes C3′-endo and other pucker states closer to A-form
            return "Stage 2: Backbone and Sugar Pucker Change"
        }
        // Stage 3: Roll and inclination to finish the transition
        else if (delta < 120.0 && P >= 180.0 && P < 360.0) {
            // This range includes C2′-endo and other pucker states closer to B-form.
            return "Stage 3: Roll and Inclination to Finish Transition"
        }

        // The thresholds (120.0 for δ, 200.0 for χ, and 180.0 for P) are approximate midpoints between typical B-form and A-form values.

        // Stage 1: Slide and Roll Change:
        // The δ angle is still high (B-form-like), and the χ angle is high (B-form-like).
        // This stage involves initial base pair movement and groove width changes.

        // Stage 2: Backbone and Sugar Pucker Change:
        // The δ angle drops (backbone deformation), and the pseudorotation angle shifts toward A-form (C3′-endo).
        // This stage involves backbone deformation and sugar pucker conversion.

        // Stage 3: Roll and Inclination to Finish Transition:
        // The δ angle remains low (A-form-like), and the pseudorotation angle is in the A-form range.
        // This stage involves final adjustments to Roll and Inclination to complete the transition.

        return "Unknown"
    }
}