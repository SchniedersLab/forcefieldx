package ffx.potential.nonbonded.octree;

/**
 * OctreeParticle:
 * Object class for Octree method presented in the Fast Multipole Method (FMM) tutorial from the Barba Group:
 * https://github.com/barbagroup/FMM_tutorial
 */
public class OctreeParticle extends OctreePoint {

    /**
     * Charge of particle
     */
    private double q = 1.0;
    /**
     * Maximum value for random point coordinate
     * Default = 1.0
     */
    private double domain = 1.0;
    /**
     * Electrostatic potential of particle
     */
    private double phi = 0.0;

    public OctreeParticle(double[] coords, double domain, double q) {
        super(coords, domain);

        setCharge(q);
        addToPhi(phi);
    }

    public void addToPhi(double phi) {
        this.phi += phi;
    }

    public double getCharge() {
        return this.q;
    }

    public void setCharge(double q) {
        this.q = q;
    }

    @Override
    public void setDomain(double domain) {
        this.domain = domain;
    }
}
