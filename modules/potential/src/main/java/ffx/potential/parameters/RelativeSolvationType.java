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
package ffx.potential.parameters;

import java.util.Comparator;
import java.util.Objects;

/**
 * A BaseType for relative solvation energies (intended for nonstandard amino
 * acids).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class RelativeSolvationType extends BaseType implements Comparator<String> {
    private final String resName;
    private final double solvEnergy;
    
    public RelativeSolvationType(String resname, double solvEnergy) {
        super(ForceField.ForceFieldType.RELATIVESOLV, resname);
        this.resName = resname;
        this.solvEnergy = solvEnergy;
    }
    
    public double getSolvEnergy() {
        return solvEnergy;
    }
    
    public String getResName() {
        return resName;
    }
    
    @Override
    public boolean equals(Object o) {
        if (o instanceof RelativeSolvationType) {
            if (((RelativeSolvationType) o).getResName().equals(resName)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 23 * hash + Objects.hashCode(this.resName);
        return hash;
    }
    
    @Override
    public String toString() {
        return String.format("relative solvation %10s %8.5f", resName, solvEnergy);
    }

    @Override
    public int compare(String o1, String o2) {
        return o1.compareTo(o2);
    }
}
