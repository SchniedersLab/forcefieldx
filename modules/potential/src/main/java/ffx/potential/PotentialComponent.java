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
     * <p>Getter for the field <code>parent</code>.</p>
     *
     * @return a {@link ffx.potential.PotentialComponent} object.
     */
    public PotentialComponent getParent() {
        return parent;
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
     * <p>getAllChildren.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<PotentialComponent> getAllChildren() {
        List<PotentialComponent> allChildren = new ArrayList<>();
        getRecursive(this, allChildren);
        return allChildren;
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
        return (parent != null) ? parent.is(category) : false;
    }

    private void getRecursive(PotentialComponent node, List<PotentialComponent> children) {
        children.addAll(node.getChildren());
        for (PotentialComponent child : children) {
            getRecursive(child, children);
        }
    }

    private void addChild(PotentialComponent child) {
        children.add(child);
    }

}
