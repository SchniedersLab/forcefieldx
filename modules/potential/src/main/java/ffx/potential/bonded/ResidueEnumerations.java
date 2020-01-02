//******************************************************************************
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
//******************************************************************************
package ffx.potential.bonded;

import java.util.Arrays;
import java.util.List;

/**
 * <p>ResidueEnumerations class.</p>
 *
 * @author Ava M. Lynn
 * @since 1.0
 */
public class ResidueEnumerations {

    /**
     * <p>getAminoAcid.</p>
     *
     * @param residueName a {@link java.lang.String} object.
     * @return a {@link ffx.potential.bonded.ResidueEnumerations.AminoAcid3} object.
     */
    public static AminoAcid3 getAminoAcid(String residueName) {
        for (AminoAcid3 aminoAcid : aminoAcidList) {
            if (aminoAcid.toString().equalsIgnoreCase(residueName)) {
                return aminoAcid;
            }
        }
        return AminoAcid3.UNK;
    }

    /**
     * <p>getAminoAcidNumber.</p>
     *
     * @param residueName a {@link java.lang.String} object.
     * @return a int.
     */
    public static int getAminoAcidNumber(String residueName) {
        int aminoAcidNumber = -1;
        for (AminoAcid3 aminoAcid : aminoAcidList) {
            aminoAcidNumber++;
            if (aminoAcid.toString().equalsIgnoreCase(residueName)) {
                break;
            }
        }
        return aminoAcidNumber;
    }

    public enum AminoAcid1 {

        G, A, V, L, I, S, T, C, X, c,
        P, F, Y, y, W, H, U, Z, D, d,
        N, E, e, Q, M, K, k, R, O, B,
        J, t, f, a, o, n, m, x
    }

    public enum AminoAcid3 {

        GLY, ALA, VAL, LEU, ILE, SER, THR, CYS(true, false), CYX, CYD(true, true),
        PRO, PHE, TYR(true, false), TYD(true, true), TRP, HIS(true, false), HID(true, true), HIE(true, true), ASP(true, false), ASH(true, true),
        ASN, GLU(true, false), GLH(true, true), GLN, MET, LYS(true, false), LYD(true, true), ARG, ORN, AIB,
        PCA, H2N, FOR, ACE, COH, NH2, NME, UNK;

        public final boolean isTitratable;
        public final boolean nonstandardProtonation;

        AminoAcid3() {
            isTitratable = false;
            nonstandardProtonation = false;
        }

        AminoAcid3(boolean titr, boolean nonstd) {
            isTitratable = titr;
            nonstandardProtonation = nonstd;
            if (nonstandardProtonation && !isTitratable) {
                throw new IllegalArgumentException(String.format(" Amino acid class %s cannot be both nonstandard and non-titratable!", this));
            }
        }
    }

    public enum NucleicAcid1 {

        A, G, C, U, D, B, I, T, O, W, H, X
    }

    /**
     * Since enumeration values must start with a letter, an 'M' is added to
     * modified bases whose IUPAC name starts with an integer.
     */
    public enum NucleicAcid3 {

        ADE, GUA, CYT, URI, DAD, DGU, DCY, DTY, THY, MP1, DP2, TP3, UNK, M2MG,
        H2U, M2G, OMC, OMG, PSU, M5MC, M7MG, M5MU, M1MA, YYG
    }

    /**
     * Constant <code>aminoAcidList</code>
     */
    public static final List<AminoAcid3> aminoAcidList = Arrays.asList(AminoAcid3.values());

    /**
     * Constant <code>nucleicAcidList</code>
     */
    public static final List<NucleicAcid3> nucleicAcidList = Arrays.asList(NucleicAcid3.values());

    /**
     * Constant <code>aminoAcidHeavyAtoms</code>
     */
    public static final int[] aminoAcidHeavyAtoms = {
            4, 5, 7, 8, 8, 6, 7, 6, 6, 6,
            7, 11, 12, 12, 14, 10, 10, 10, 8, 8,
            8, 9, 9, 9, 8, 9, 9, 11, 8, 6,
            8, 0, 0, 0, 0, 0, 0, 0
    };
}
