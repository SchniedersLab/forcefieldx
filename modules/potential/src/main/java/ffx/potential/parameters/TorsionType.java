//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.parameters;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;

import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * The TorsionType class defines a torsional angle.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class TorsionType extends BaseType implements Comparator<String> {

    /**
     * Torsion modes include Normal or In-Plane
     */
    public enum TorsionMode {

        NORMAL, IMPROPER
    }

    /**
     * Atom classes that for this Torsion angle.
     */
    public final int[] atomClasses;
    /**
     * Number of terms in the Fourier series.
     */
    public final int terms;
    /**
     * Amplitudes of the Fourier series.
     */
    public final double[] amplitude;
    /**
     * Phases of the Fourier series in degrees.
     */
    public final double[] phase;
    /**
     * Cosine of the phase angle.
     */
    public final double[] cosine;
    /**
     * Sine of the phase angle.
     */
    public final double[] sine;
    /**
     * Periodicity of the Fourier series.
     */
    private final int[] periodicity;
    /**
     * The torsion mode in use.
     */
    private TorsionMode torsionMode;

    /**
     * TorsionType Constructor.
     *
     * @param atomClasses Atom classes.
     * @param amplitude   Amplitudes of the Fourier series.
     * @param phase       Phases of the Fourier series in degrees.
     * @param periodicity Periodicity of the Fourier series.
     */
    public TorsionType(int[] atomClasses, double[] amplitude, double[] phase, int[] periodicity) {
        super(ForceField.ForceFieldType.TORSION, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        int max = 1;
        for (int i1 : periodicity) {
            if (i1 > max) {
                max = i1;
            }
        }
        terms = max;
        if (periodicity.length != max) {
            this.amplitude = new double[max];
            this.phase = new double[max];
            this.periodicity = new int[max];
            for (int i = 0; i < periodicity.length; i++) {
                int j = periodicity[i] - 1;
                this.amplitude[j] = amplitude[i];
                this.phase[j] = phase[i];
                this.periodicity[j] = periodicity[i];
            }
        } else {
            this.amplitude = amplitude;
            this.phase = phase;
            this.periodicity = periodicity;
        }
        cosine = new double[terms];
        sine = new double[terms];
        for (int i = 0; i < terms; i++) {
            double angle = toRadians(this.phase[i]);
            cosine[i] = cos(angle);
            sine[i] = sin(angle);
        }
        torsionMode = TorsionMode.NORMAL;
    }

    /**
     * TorsionType Constructor.
     *
     * @param atomClasses Atom classes.
     * @param amplitude   Amplitudes of the Fourier series.
     * @param phase       Phases of the Fourier series in degrees.
     * @param periodicity Periodicity of the Fourier series.
     */
    public TorsionType(int[] atomClasses, double[] amplitude, double[] phase, int[] periodicity, TorsionMode torsionMode) {
        this(atomClasses, amplitude, phase, periodicity);
        if (torsionMode == TorsionMode.IMPROPER) {
            this.torsionMode = torsionMode;
            forceFieldType = ForceField.ForceFieldType.IMPROPER;
        }
    }

    /**
     * <p>setScaleFactor.</p>
     *
     * @param scale a double.
     */
    void setScaleFactor(double scale) {
        for (int i = 0; i < amplitude.length; i++) {
            amplitude[i] *= scale;
        }
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    public void incrementClasses(int increment) {
        for (int i = 0; i < atomClasses.length; i++) {
            if (atomClasses[i] != 0) {
                atomClasses[i] += increment;
            }
        }
        setKey(sortKey(atomClasses));
    }

    /**
     * Remap new atom classes to known internal ones.
     *
     * @param typeMap a lookup between new atom types and known atom types.
     * @return a {@link ffx.potential.parameters.TorsionType} object.
     */
    public TorsionType patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        int len = atomClasses.length;

        // Look for new TorsionTypes that contain 1 to 3 mapped atom classes.
        for (AtomType newType : typeMap.keySet()) {
            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
                    count++;
                }
            }
        }

        // If found, create a new TorsionType that bridges to known classes.
        if (count == 1 || count == 2 || count == 3) {
            int[] newClasses = copyOf(atomClasses, len);
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < len; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        newClasses[i] = knownType.atomClass;
                    }
                }
            }
            return new TorsionType(newClasses, amplitude, phase, periodicity);
        }
        return null;
    }

    /**
     * <p>average.</p>
     *
     * @param torsionType1 a {@link ffx.potential.parameters.TorsionType} object.
     * @param torsionType2 a {@link ffx.potential.parameters.TorsionType} object.
     * @param atomClasses  an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.TorsionType} object.
     */
    public static TorsionType average(TorsionType torsionType1, TorsionType torsionType2, int[] atomClasses) {
        if (torsionType1 == null || torsionType2 == null || atomClasses == null) {
            return null;
        }
        int len = torsionType1.amplitude.length;
        if (len != torsionType2.amplitude.length) {
            return null;
        }
        double[] amplitude = new double[len];
        double[] phase = new double[len];
        int[] periodicity = new int[len];
        for (int i = 0; i < len; i++) {
            amplitude[i] = (torsionType1.amplitude[i] + torsionType2.amplitude[i]) / 2.0;
            phase[i] = (torsionType1.phase[i] + torsionType2.phase[i]) / 2.0;
            periodicity[i] = (torsionType1.periodicity[i] + torsionType2.periodicity[i]) / 2;
        }

        return new TorsionType(atomClasses, amplitude, phase, periodicity);
    }

    /**
     * This method sorts the atom classes for the torsion.
     *
     * @param c atomClasses
     * @return lookup key
     * @since 1.0
     */
    public static String sortKey(int[] c) {
        if (c == null || c.length != 4) {
            return null;
        }
        if (c[1] < c[2]) {
            // Do nothing.
        } else if (c[2] < c[1]) {
            // Reverse the order.
            int temp = c[0];
            c[0] = c[3];
            c[3] = temp;
            temp = c[1];
            c[1] = c[2];
            c[2] = temp;
        } else if (c[1] == c[2]) {
            if (c[0] > c[3]) {
                // Reverse the order.
                int temp = c[0];
                c[0] = c[3];
                c[3] = temp;
            }
        } else if (c[0] <= c[3]) {
            // Do nothing.
        } else {
            // Reverse the order.
            int temp = c[0];
            c[0] = c[3];
            c[3] = temp;
            temp = c[1];
            c[1] = c[2];
            c[2] = temp;
        }

        return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Torsion angle.
     *
     * @since 1.0
     */
    @Override
    public String toString() {
        StringBuilder torsionBuffer;
        if (torsionMode == TorsionMode.IMPROPER) {
            torsionBuffer = new StringBuilder("improper");
        } else {
            torsionBuffer = new StringBuilder("torsion");
        }
        for (int i : atomClasses) {
            torsionBuffer.append(format(" %5d", i));
        }
        for (int i = 0; i < amplitude.length; i++) {
            if (torsionMode == TorsionMode.NORMAL) {
                torsionBuffer.append(format(" %7.3f %7.3f %1d", amplitude[i], phase[i], periodicity[i]));
            } else {
                torsionBuffer.append(format(" %7.3f %7.3f", amplitude[i], phase[i]));
            }
        }
        return torsionBuffer.toString();
    }

    /**
     * {@inheritDoc}
     *
     * @since 1.0
     */
    @Override
    public int compare(String s1, String s2) {
        String[] keys1 = s1.split(" ");
        String[] keys2 = s2.split(" ");
        int[] c1 = new int[4];
        int[] c2 = new int[4];

        for (int i = 0; i < 4; i++) {
            c1[i] = Integer.parseInt(keys1[i]);
            c2[i] = Integer.parseInt(keys2[i]);
        }

        if (c1[1] < c2[1]) {
            return -1;
        } else if (c1[1] > c2[1]) {
            return 1;
        } else if (c1[2] < c2[2]) {
            return -1;
        } else if (c1[2] > c2[2]) {
            return 1;
        } else if (c1[0] < c2[0]) {
            return -1;
        } else if (c1[0] > c2[0]) {
            return 1;
        } else if (c1[3] < c2[3]) {
            return -1;
        } else if (c1[3] > c2[3]) {
            return 1;
        }

        return 0;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Override the default <code>equals</code> method.
     *
     * @since 1.0
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof TorsionType)) {
            return false;
        }
        TorsionType torsionType = (TorsionType) other;
        for (int i = 0; i < 4; i++) {
            if (torsionType.atomClasses[i] != atomClasses[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Implementation of the <code>hashCode</code> method.
     *
     * @since 1.0
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 89 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }

}
