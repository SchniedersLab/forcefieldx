package ffx.potential.nonbonded.implicit;

public class OctreeParticle extends OctreePoint{

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

    public OctreeParticle(double[] coords, double domain,double q) {
        super(coords, domain);

        setCharge(q);
        setPhi(phi);
    }

    @Override
    public void setDomain(double domain) {
        this.domain = domain;
    }

    public double getCharge(){
        return this.q;
    }

    public void setCharge(double q){
        this.q = q;
    }

    public void setPhi(double phi){
        this.phi += phi;
    }
}
