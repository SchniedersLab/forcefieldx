package ffx.potential;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 */
public enum PotentialComponent {

	Topology (null),
		ForceFieldEnergy (Topology),
			VanDerWaals (ForceFieldEnergy),
			Bonded (ForceFieldEnergy),
				Bond (Bonded),
				Angle (Bonded),
				Torsion (Bonded),
				StretchBend (Bonded),
				OutOfPlaneBend (Bonded),
				PiOrbitalTorsion (Bonded),
				TorsionTorsion (Bonded),
				UreyBradley (Bonded),
				RestraintBond (Bonded),
				ImproperTorsion (Bonded),
				NCS (Bonded),
				Restrain (Bonded),
			Multipoles (ForceFieldEnergy),
				Permanent (Multipoles),
					PermanentRealSpace (Permanent),
					PermanentSelf (Permanent),
					PermanentReciprocal (Permanent),
				Induced (Multipoles),
					InducedRealSpace (Induced),
					InducedSelf (Induced),
					InducedReciprocal (Induced),
			GeneralizedKirkwood (ForceFieldEnergy),
		Bias (Topology),
			OSRW (Bias),
			pHMD (Bias),
				Acidostat (pHMD),
				Discretizer (pHMD),
		XRay (Topology),
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

	public PotentialComponent getParent() {
		return parent;
	}
	
	public List<PotentialComponent> getChildren() {
		return children;
	}
	
	public List<PotentialComponent> getAllChildren() {
		List<PotentialComponent> allChildren = new ArrayList<>();
		getRecursive(this, allChildren);
		return allChildren;
	}
	
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