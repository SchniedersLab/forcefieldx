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
package ffx.potential;

import java.util.ArrayList;
import java.util.List;

/**
 * <p>PotentialComponent class.</p>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public enum PotentialComponent {

    Topology(null),
    ForceFieldEnergy(Topology),
    VanDerWaals(ForceFieldEnergy),
    Bonded(ForceFieldEnergy),
    Bond(Bonded),
    Angle(Bonded),
    Torsion(Bonded),
    StretchBend(Bonded),
    OutOfPlaneBend(Bonded),
    PiOrbitalTorsion(Bonded),
    TorsionTorsion(Bonded),
    UreyBradley(Bonded),
    RestraintBond(Bonded),
    ImproperTorsion(Bonded),
    NCS(Bonded),
    Restrain(Bonded),
    Electrostatics(ForceFieldEnergy),
    Multipoles(Electrostatics),
    Permanent(Multipoles),
    PermanentRealSpace(Permanent),
    PermanentSelf(Permanent),
    PermanentReciprocal(Permanent),
    Induced(Multipoles),
    InducedRealSpace(Induced),
    InducedSelf(Induced),
    InducedReciprocal(Induced),
    GeneralizedKirkwood(Electrostatics),
    Bias(Topology),
    OSRW(Bias),
    pHMD(Bias),
    Acidostat(pHMD),
    Discretizer(pHMD),
    XRay(Topology),
    ;


    private final PotentialComponent parent;
    private final List<PotentialComponent> children;

    PotentialComponent(PotentialComponent parent) {
        this.parent = parent;
        this.children = new ArrayList<>();
        if (parent != null) {
            parent.addChild(this);
        }
    }

    /**
     * <p>Getter for the field <code>children</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<PotentialComponent> getChildren() {
        return children;
    }

    /**
     * <p>Getter for the field <code>parent</code>.</p>
     *
     * @return a {@link ffx.potential.PotentialComponent} object.
     */
    public PotentialComponent getParent() {
        return parent;
    }

    /**
     * <p>is.</p>
     *
     * @param category a {@link ffx.potential.PotentialComponent} object.
     * @return a boolean.
     */
    public boolean is(PotentialComponent category) {
        if (this == category) {
            return true;
        }
        return (parent != null) && parent.is(category);
    }

    private void addChild(PotentialComponent child) {
        children.add(child);
    }

}
