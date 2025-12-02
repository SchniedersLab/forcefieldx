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
package ffx.xray.scatter;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * Represents neutron scattering parameters for a given atom, including its name, atomic number,
 * and form factor.
 * <p>
 * This class uses a static cache to store precomputed neutron scattering parameters for atoms
 * and allows for retrieval of these parameters by atom name.
 *
 * @see <a href="https://www.ncnr.nist.gov/resources/n-lengths/list.html" target="_blank">
 * NIST Center for Neutron Research</a>
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000594" target="_blank"> V. F. Sears, Int.
 * Tables Vol. C (2006). Table 4.4.4.1</a>
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000600" target="_blank"> B. T. M. Willis,
 * Int. Tables Vol. C (2006). Chapter 6.1.3</a>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public record NeutronScatteringParameters(String name, int atomicNumber, double[] formFactor) {

    private static final Logger logger = Logger.getLogger(NeutronScatteringParameters.class.getName());

    /**
     * Retrieves the neutron scattering parameters (form factor) for a specified atom.
     * If the atom does not exist in the form factors map, returns null.
     *
     * @param atom the name of the atom for which the form factor is requested
     * @return the neutron scattering parameters for the specified atom,
     * or null if the atom is not found in the form factors map
     */
    public static NeutronScatteringParameters getFormFactor(String atom) {
        if (formfactors.containsKey(atom)) {
            return formfactors.get(atom);
        } else {
            String message = " Form factor for atom: " + atom + " not found!";
            logger.severe(message);
        }
        return null;
    }

    /**
     * Cache of created form factors.
     */
    private static final HashMap<String, NeutronScatteringParameters> formfactors = new HashMap<>();

    private static final String[] atomNames = {
            "H", "D", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
            "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
            "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
            "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
            "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
            "Bi", "Th", "U"
    };
    private static final String[] atomicNumbers = {
            "1_1", "1_2", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
            "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
            "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49",
            "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "62", "63", "64", "65", "66",
            "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82",
            "83", "90", "92"
    };
    private static final double[][] neutronFormFactors = {
            {-3.7390, 0.0},  // H
            {6.671, 0.0},    // D
            {3.26, 0.0},     // He
            {-1.90, 0.0},    // Li
            {7.79, 0.0},     // Be
            {5.30, 0.213},   // B
            {6.6460, 0.0},   // C
            {9.36, 0.0},     // N
            {5.803, 0.0},    // O
            {5.654, 0.0},    // F
            {4.566, 0.0},   // Ne
            {3.63, 0.0},    // Na
            {5.375, 0.0},   // Mg
            {3.449, 0.0},   // Al
            {4.1491, 0.0},  // Si
            {5.13, 0.0},    // P
            {2.847, 0.0},   // S
            {9.5770, 0.0},  // Cl
            {1.909, 0.0},   // Ar
            {3.67, 0.0},    // K
            {4.70, 0.0},    // Ca
            {12.29, 0.0},   // Sc
            {-3.370, 0.0},  // Ti
            {-0.3824, 0.0}, // V
            {3.635, 0.0},   // Cr
            {-3.750, 0.0},  // Mn
            {9.45, 0.0},    // Fe
            {2.49, 0.0},    // Co
            {10.3, 0.0},    // Ni
            {7.718, 0.0},   // Cu
            {5.60, 0.0},    // Zn
            {7.288, 0.0},   // Ga
            {8.185, 0.0},   // Ge
            {6.58, 0.0},    // As
            {7.970, 0.0},   // Se
            {6.795, 0.0},   // Br
            {7.81, 0.0},    // Kr
            {7.09, 0.0},    // Rb
            {7.02, 0.0},    // Sr
            {7.75, 0.0},    // Y
            {7.16, 0.0},    // Zr
            {7.054, 0.0},   // Nb
            {6.715, 0.0},   // Mo
            {6.80, 0.0},    // Tc
            {7.03, 0.0},    // Ru
            {5.88, 0.0},    // Rh
            {5.91, 0.0},    // Pd
            {5.922, 0.0},   // Ag
            {4.87, -0.70},  // Cd
            {4.065, -0.0539}, // In
            {6.225, 0.0},   // Sn
            {5.57, 0.0},    // Sb
            {5.80, 0.0},    // Te
            {5.28, 0.0},    // I
            {4.92, 0.0},    // Xe
            {5.42, 0.0},    // Cs
            {5.07, 0.0},    // Ba
            {8.24, 0.0},    // La
            {4.84, 0.0},    // Ce
            {4.58, 0.0},    // Pr
            {7.69, 0.0},    // Nd
            // No PM
            {0.80, -1.65},  // Sm
            {7.22, -1.26},  // Eu
            {6.5, -13.82},  // Gd
            {7.38, 0.0},    // Tb
            {16.9, -0.276}, // Dy
            {8.01, 0.0},    // Ho
            {7.79, 0.0},    // Er
            {7.07, 0.0},    // Tm
            {12.43, 0.0},   // Yb
            {7.21, 0.0},    // Lu
            {7.77, 0.0},    // Hf
            {6.91, 0.0},    // Ta
            {4.86, 0.0},    // W
            {9.2, 0.0},     // Re
            {10.7, 0.0},    // Os
            {10.6, 0.0},    // Ir
            {9.60, 0.0},    // Pt
            {7.63, 0.0},    // Au
            {12.692, 0.0},  // Hg
            {8.776, 0.0},   // Tl
            {9.405, 0.0},   // Pb
            {8.532, 0.0},   // Bi
            {10.31, 0.0},   // Th
            {8.417, 0.0}    // U
    };

    // Load form factors.
    static {
        for (int i = 0; i < atomNames.length; i++) {
            String number = atomicNumbers[i];
            if (number.contains("_")) {
                number = "1";
            }
            int atomicNumber = Integer.parseInt(number);
            NeutronScatteringParameters parameters = new NeutronScatteringParameters(
                    atomNames[i], atomicNumber, neutronFormFactors[i]);
            formfactors.put(atomicNumbers[i], parameters);
        }
    }

}
