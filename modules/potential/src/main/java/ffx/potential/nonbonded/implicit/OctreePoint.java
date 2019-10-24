package ffx.potential.nonbonded.implicit;

import java.util.logging.Logger;

public class OctreePoint {

    private static final Logger logger = Logger.getLogger(OctreePoint.class.getName());

    /**
     * Coordinates of the point
     */
    private double x;
    private double y;
    private double z;

    /**
     * Maximum value for random point coordinate
     * Default = 1.0
     */
    private double domain = 1.0;

    public OctreePoint(double coords[], double domain){

        // Set max for random point coordinate generation
        setDomain(domain);

        // Assign input coordinates, if they are given
        // Otherwise, assign random coordinates between 0 and domain
        if(coords.length > 0){
            if(coords.length == 3){
                this.x = coords[0];
                this.y = coords[1];
                this.z = coords[2];
            }else{
                logger.warning("Coordinate array must have three points");
            }
        } else{
            this.x = this.domain * Math.random();
            this.y = this.domain * Math.random();
            this.z = this.domain * Math.random();
        }
    }

    public void setDomain(double domain){
        this.domain = domain;
    }

    public double distance(OctreePoint other){
        return Math.sqrt(Math.pow((this.x - other.x),2)
                +Math.pow((this.y - other.y),2)
                +Math.pow((this.z-other.z),2));
    }

    public double getX(){ return this.x; }
    public double getY(){ return this.y; }
    public double getZ(){ return this.z; }
}