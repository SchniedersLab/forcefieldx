/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential.bonded;

/**
 * This exception is thrown when a nucleic acid Rotamer must be distorted to an
 * excessive degree to meet the prior residue's 3' end.
 *
 * @author Jacob M. Litman
 */
public class NACorrectionTooLargeException extends RuntimeException {

    public final double correctionThreshold;
    public final double correction;
    public final Residue residue;
    public final Rotamer rotamer;

    public NACorrectionTooLargeException(double correctionThreshold, double correction, Residue residue, Rotamer rotamer) {
        this.correctionThreshold = correctionThreshold;
        this.correction = correction;
        this.residue = residue;
        this.rotamer = rotamer;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(super.toString());
        if (residue == null) {
            sb.append("\n NACorrectionTooLargeException was mistakenly thrown on a null residue");
        } else {
            sb.append("\n Nucleic acid Residue " + residue.toString());
            Residue prevResidue = residue.getPreviousResidue();
            if (prevResidue != null) {
                sb.append(" had to be shifted too far to meet up with prior residue ");
                sb.append(prevResidue.toString());
                sb.append("\n Correction threshold was " + correctionThreshold + " Angstroms");
            } else {
                sb.append(" had NACorrectionTooLargeException mistakenly "
                        + "thrown, despite the absence of a prior Residue");
            }
        }
        return sb.toString();
    }

    public double getCorrection() {
        return correction;
    }
}
