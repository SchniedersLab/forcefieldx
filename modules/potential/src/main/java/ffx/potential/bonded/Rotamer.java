/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.bonded;

import static org.apache.commons.math3.util.FastMath.max;

import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import java.util.Arrays;

/**
 * The Rotamer Class usually represents one immutable amino acid Rotamer. It is
 * additionally being extended to represent one nucleic acid Rotamer.
 *
 * @author Ava M. Lynn
 * @author Jacob M. Litman
 */
public class Rotamer {

    /**
     * Torsions chi 1-4 are used for amino acids and nucleic acids.
     */
    public final double chi1;
    public final double chi2;
    public final double chi3;
    public final double chi4;
    /**
     * Torsions chi 5-7 are only currently used for nucleic acids.
     */
    public final double chi5;
    public final double chi6;
    public final double chi7;
    public final double angles[];
    public final double sigmas[];
    public final AminoAcid3 name;
    public final NucleicAcid3 nucleicName;
    public final int length;
    public final ResidueState originalState;
    public final boolean isState;
    
    /**
     * Constructs a Rotamer from a Residue, an array of torsions, and optionally
     * an array of torsion bin widths. Intended to be agnostic to AA vs. NA vs.
     * other, and not require explicitly passed-in sigmas.
     * 
     * @param res Residue to construct rotamer for
     * @param chis Torsion angles
     * @param sigmas Torsion angle bin widths (optional)
     */
    public Rotamer(Residue res, double[] chis, double[] sigmas) {
        int nChi = chis.length;
        angles = new double[nChi];
        Arrays.fill(angles, 0);
        System.arraycopy(chis, 0, angles, 0, nChi);
        
        // Hooray, 
        double[] tempVals = new double[7];
        Arrays.fill(tempVals, 0);
        System.arraycopy(chis, 0, tempVals, 0, nChi);
        
        chi1 = tempVals[0];
        chi2 = tempVals[1];
        chi3 = tempVals[2];
        chi4 = tempVals[3];
        chi5 = tempVals[4];
        chi6 = tempVals[5];
        chi7 = tempVals[6];
        
        this.sigmas = new double[nChi];
        if (sigmas != null) {
            System.arraycopy(sigmas, 0, this.sigmas, 0, nChi);
        } else {
            Arrays.fill(sigmas, 0);
        }
        
        switch(res.getResidueType()) {
            case AA:
                this.name = res.getAminoAcid3();
                this.nucleicName = null;
                break;
            case NA:
                this.name = null;
                this.nucleicName = res.getNucleicAcid3();
                break;
            case UNK:
            default:
                this.name = null;
                this.nucleicName = null;
                break;
        }
        
        length = nChi;
        originalState = null;
        isState = false;
    }

    public Rotamer(AminoAcid3 name, double... values) {
        length = values.length / 2;
        angles = new double[max(length, 4)];
        sigmas = new double[max(length, 4)];
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        this.name = name;
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = chi6 = chi7 = 0;
        nucleicName = null;
        originalState = null;
        isState = false;
    }

    public Rotamer(NucleicAcid3 name, double... values) {
        length = values.length / 2;
        angles = new double[max(length, 7)];
        sigmas = new double[max(length, 7)];
        nucleicName = name;
        this.name = null;
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = angles[4];
        chi6 = angles[5];
        chi7 = angles[6];
        originalState = null;
        isState = false;
    }
    
    /**
     * Constructor for unknown residue types.
     * @param values 
     */
    public Rotamer(double... values) {
        length = values.length / 2;
        angles = new double[max(length, 7)];
        sigmas = new double[max(length, 7)];
        nucleicName = null;
        this.name = null;
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = angles[4];
        chi6 = angles[5];
        chi7 = angles[6];
        originalState = null;
        isState = false;
    }

    public Rotamer(AminoAcid3 name, ResidueState residueState, double... values) {
        length = values.length / 2;
        angles = new double[max(length, 4)];
        sigmas = new double[max(length, 4)];
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        this.name = name;
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = chi6 = chi7 = 0;
        nucleicName = null;
        originalState = residueState;
        isState = true;
    }

    public Rotamer(NucleicAcid3 name, ResidueState residueState, double... values) {
        length = values.length / 2;
        angles = new double[max(length, 7)];
        sigmas = new double[max(length, 7)];
        nucleicName = name;
        this.name = null;
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = angles[4];
        chi6 = angles[5];
        chi7 = angles[6];
        originalState = residueState;
        isState = true;
    }

    /**
     * Constructor for unknown residue types.
     * @param residueState
     * @param values 
     */
    public Rotamer(ResidueState residueState, double... values) {
        length = values.length / 2;
        angles = new double[max(length, 7)];
        sigmas = new double[max(length, 7)];
        nucleicName = null;
        this.name = null;
        for (int i = 0; i < length; i++) {
            int ii = 2 * i;
            angles[i] = values[ii];
            sigmas[i] = values[ii + 1];
        }
        chi1 = angles[0];
        chi2 = angles[1];
        chi3 = angles[2];
        chi4 = angles[3];
        chi5 = angles[4];
        chi6 = angles[5];
        chi7 = angles[6];
        originalState = residueState;
        isState = true;
    }

    @Override
    public String toString() {
        StringBuilder sb;
        if (name != null) {
            sb = new StringBuilder(name.toString());
        } else {
            sb = new StringBuilder(nucleicName.toString());
        }
        int n = angles.length;
        for (int i = 0; i < n; i++) {
            sb.append(String.format(" %6.1f %4.1f", angles[i], sigmas[i]));
        }
        return sb.toString();
    }
}
