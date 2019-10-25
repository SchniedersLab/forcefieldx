package ffx.potential.nonbonded.implicit;

import org.apache.commons.lang3.BooleanUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * Octree method based on Fast Multipole Method (FMM) tutorial from the Barba Group:
 * https://github.com/barbagroup/FMM_tutorial
 */
public class Octree {

    private static final Logger logger = Logger.getLogger(Octree.class.getName());

    /**
     * List of cells
     */
    private ArrayList<OctreeCell> cells = new ArrayList<>();
    /**
     * Critical (maximum allowed+1) number of points allowed in any one cell:
     * If a cell contains nCritical points, it needs to be split
     */
    private int nCritical = 10;
    /**
     * List of particles
     */
    private ArrayList<OctreeParticle> particles = new ArrayList<>();
    /**
     * List of all leaf cells
     */
    private ArrayList<OctreeCell> leaves = new ArrayList<>();
    /**
     * Tolerance parameter
     */
    private double theta = 0.5;

    /**
     * Default constructor: only need to pass in a list of particles
     * nCritical and theta set to defaults
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     */
    public Octree(ArrayList<OctreeParticle> particles) {
        this.particles = particles;
        this.nCritical = 10;
        this.theta = 0.5;
    }

    /**
     * Constructor allowing the specification of nCritical, default theta value
     * @param nCritical Critical number of particles; cells must split when they reach nCritical
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     */
    public Octree(int nCritical,ArrayList<OctreeParticle> particles){
        this.nCritical = nCritical;
        this.particles = particles;
        this.theta = 0.5;
    }

    /**
     * Constructor allowing the specification of nCritical and theta
     * @param nCritical Critical number of particles; cells must split when they reach nCritical
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     * @param theta Specifies near field vs far field
     */
    public Octree(int nCritical,ArrayList<OctreeParticle> particles,double theta){
        this.nCritical = nCritical;
        this.particles = particles;
        this.theta = theta;
    }

    public void addChild(int octant, int p){
        OctreeCell tempCell = new OctreeCell(nCritical);

        // Create new cell
        // TODO: should the cells array be passed in or the "this" value from the particular Octree?
        cells.add(tempCell);

        // Last element of cells list is new child, c
        int c = cells.size() - 1;

        // Geometric reference between parent and child
        cells.get(c).setR((cells.get(p).getR())*0.5);
        cells.get(c).setX((cells.get(p).getX())*((octant & 1)*2-1));
        cells.get(c).setY((cells.get(p).getY())*((octant & 2)-1));
        cells.get(c).setZ((cells.get(p).getZ())*((octant & 4)/2-1));

        // Establish mutual reference in cells list
        cells.get(c).setParentIndex(p);
        cells.get(c).setChildren(octant,c);
        cells.get(c).setnChild(cells.get(p).getnChild() | (1 << octant));
    }

    private void splitCell(int p){
        for (int i = 0; i < nCritical; i++){
            int octX = 0;
            int octY = 0;
            int octZ = 0;

            if(particles.get(i).getX() > cells.get(p).getX()){octX = 1;}
            if(particles.get(i).getY() > cells.get(p).getY()){octY = 1;}
            if(particles.get(i).getZ() > cells.get(p).getZ()){octZ = 1;}

            // Find particle's octant - should be an integer from 0 to 7
            int octant = octX + (octY << 1) + (octZ << 2);

            // If there's not a child cell in the particle's octant, create one
            boolean noChildInOctant = BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << octant));
            if(noChildInOctant){
                addChild(octant,p);
            }

            // Reallocate the particle in the child cell
            int c = cells.get(p).getChildAtIndex(octant);
            cells.get(c).setLeaf(cells.get(c).getNumLeaves(),1);
            cells.get(c).setNumLeaves(cells.get(c).getNumLeaves() + 1);

            // Check if child cell reaches nCritical - split recursively if so
            if(cells.get(c).getNumLeaves() >= nCritical){
                splitCell(c);
            }

        }
    }

    public void buildTree(OctreeCell root){
        //ArrayList<OctreeCell> cells = new ArrayList<>();

        // Set root cell
        cells.add(root);

        // Build tree
        int n = particles.size();

        for(int i = 0; i < n; i++){
            int current = 0;

            while(cells.get(current).getNumLeaves() >= nCritical){
                cells.get(current).setNumLeaves(cells.get(current).getNumLeaves() + 1);
                int octX = 0;
                int octY = 0;
                int octZ = 0;

                if(particles.get(i).getX() > cells.get(current).getX()){octX = 1;}
                if(particles.get(i).getY() > cells.get(current).getY()){octY = 1;}
                if(particles.get(i).getZ() > cells.get(current).getZ()){octZ = 1;}

                // Find particle's octant - should be an integer from 0 to 7
                int octant = octX + (octY << 1) + (octZ << 2);

                // If there's not a child cell in the particle's octant, create one
                boolean noChildInOctant = BooleanUtils.toBoolean(cells.get(current).getnChild() & (1 << octant));
                if(noChildInOctant){
                    addChild(octant,current);
                }

                current = cells.get(current).getChildAtIndex(octant);
            }

            // Allocate the particle in the leaf cell
            cells.get(current).setLeaf(cells.get(current).getNumLeaves(),i);
            cells.get(current).setNumLeaves(cells.get(current).getNumLeaves() + 1);

            // Check whether to split cell
            if(cells.get(current).getNumLeaves() >= nCritical){
                splitCell(current);
            }
        }

        //return cells;
    }

    private void getMultipole(int p,ArrayList<Integer> leaves){
        // If the current cell is not a leaf cell, traverse down
        if(cells.get(p).getNumLeaves() >= nCritical){
            for(int c = 0; c < 8; c++){
                if(BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << c))){
                    getMultipole(cells.get(p).getChildAtIndex(c),leaves);
                }
            }
        } else{ // Otherwise, cell p is a leaf cell
            // Loop in leaf particles, do P2M
            for(int i = 0; i < cells.get(p).getNumLeaves();i++){
                int l = cells.get(p).getLeavesValueAtIndex(i);
                double dx = cells.get(p).getX()-particles.get(l).getX();
                double dy = cells.get(p).getY()-particles.get(l).getY();
                double dz = cells.get(p).getZ()-particles.get(l).getZ();

                // Calculate Multipole and fill array
                double charge = particles.get(l).getCharge();
                double[] calculatedMultipole = new double[10];
                // Monopole (charge)
                calculatedMultipole[0] = charge;
                // Dipole
                calculatedMultipole[1] = dx * charge;
                calculatedMultipole[2] = dy * charge;
                calculatedMultipole[3] = dz * charge;
                // Quadropole
                calculatedMultipole[4] = (Math.pow(dx,2)*0.5)*charge;
                calculatedMultipole[5] = (Math.pow(dy,2)*0.5)*charge;
                calculatedMultipole[6] = (Math.pow(dz,2)*0.5)*charge;
                calculatedMultipole[7] = ((dx*dy)*0.5)*charge;
                calculatedMultipole[8] = ((dy*dz)*0.5)*charge;
                calculatedMultipole[9] = ((dz*dx)*0.5)*charge;

                // Set Multipole
                cells.get(p).addToMultipole(calculatedMultipole);

                // TODO: decide if leaves should be an array or ArrayList and adjust accordingly
                leaves.add(p);
            }
        }
    }

    private void M2M(int p, int c){
        double dx = cells.get(p).getX()-cells.get(c).getX();
        double dy = cells.get(p).getY()-cells.get(c).getY();
        double dz = cells.get(p).getZ()-cells.get(c).getZ();

        double[] Dxyz = new double[]{dx, dy, dz};
        double[] Dyzx = new double[]{dy, dz, dx};

        cells.get(p).addToMultipole(cells.get(c).getMultipole());

        double[] currentChildMultipole = cells.get(c).getMultipole();

        // Additional Multipole Terms
        double[] additionalMultipoleTerms = new double[10];
        // Added to charge
        additionalMultipoleTerms[0] = 0;
        // Added to Dipole
        additionalMultipoleTerms[1] = currentChildMultipole[0]*Dxyz[0];
        additionalMultipoleTerms[2] = currentChildMultipole[0]*Dxyz[1];
        additionalMultipoleTerms[3] = currentChildMultipole[0]*Dxyz[2];
        // Added to Quadropole
        additionalMultipoleTerms[4] = currentChildMultipole[1]*Dxyz[0]+0.5*currentChildMultipole[0]*Math.pow(Dxyz[0],2);
        additionalMultipoleTerms[5] = currentChildMultipole[2]*Dxyz[1]+0.5*currentChildMultipole[0]*Math.pow(Dxyz[1],2);
        additionalMultipoleTerms[6] = currentChildMultipole[3]*Dxyz[2]+0.5*currentChildMultipole[0]*Math.pow(Dxyz[2],2);

        additionalMultipoleTerms[7]=0.5*currentChildMultipole[2]*Dyzx[0]+0.5*currentChildMultipole[1]*Dxyz[0]+
                0.5*currentChildMultipole[0]*Dxyz[0]*Dyzx[0];
        additionalMultipoleTerms[8]=0.5*currentChildMultipole[3]*Dyzx[1]+0.5*currentChildMultipole[2]*Dxyz[1]+
                0.5*currentChildMultipole[0]*Dxyz[1]*Dyzx[1];
        additionalMultipoleTerms[9]=0.5*currentChildMultipole[1]*Dyzx[2]+0.5*currentChildMultipole[3]*Dxyz[2]+
                0.5*currentChildMultipole[0]*Dxyz[2]*Dyzx[2];

        cells.get(p).addToMultipole(additionalMultipoleTerms);
    }

    public void upwardSweep(){
        for(int c = cells.size(); c > 0;c++){
            int p = cells.get(c).getParentIndex();
            M2M(p,c);
        }
    }

    public void directSum(){
        for(int i = 0; i < particles.size();i++){
            for(int j = 0; j < particles.size();j++){
                if(j!=i){
                    double r = particles.get(i).distance(particles.get(j));
                    particles.get(j).addToPhi(particles.get(j).getCharge()/r);
                }
            }
        }
        // Reset potential for all particles
        for(int i = 0; i < particles.size(); i++){
            particles.get(i).addToPhi(0);
        }
    }

    public double distance(double[] array, OctreePoint point){
        return Math.sqrt(Math.pow((array[0]-point.getX()),2)
                +Math.pow((array[1]-point.getY()),2)
                +Math.pow((array[2]-point.getZ()),2));
    }

    /**
     * Evaluate potential at one target
     * @param p Index of parent cell
     * @param i Index of target particle
     */
    private void evalAtTarget(int p, int i){

        // Non-leaf cell
        if(cells.get(p).getNumLeaves() >= nCritical){

            // Loop through p's child cells (8 octants)
            for(int oct = 0; oct < 8; oct++){
                if(BooleanUtils.toBoolean(cells.get(p).getnChild() & (1<<oct))){
                    int c = cells.get(p).getChildAtIndex(oct);
                    double r = particles.get(i).distance(cells.get(c));

                    // Near field child cell
                    if(cells.get(c).getR() > theta*r){
                        evalAtTarget(c,i);
                    } else{ // Far field child cell
                        double dx = particles.get(i).getX()-cells.get(c).getX();
                        double dy = particles.get(i).getY()-cells.get(c).getY();
                        double dz = particles.get(i).getZ()-cells.get(c).getZ();
                        double r3 = Math.pow(r,3);
                        double r5 = r3*Math.pow(r,2);

                        // Calculate the weight from each multipole
                        double[] weight = new double[10];
                        weight[0] = 1/r;
                        weight[1] = -dx/r3;
                        weight[2] = -dy/r3;
                        weight[3] = -dz/r3;
                        weight[4] = (3*Math.pow(dx,2))/r5 - (1/r3);
                        weight[5] = (3*Math.pow(dy,2))/r5 - (1/r3);
                        weight[6] = (3*Math.pow(dz,2))/r5 - (1/r3);
                        weight[7] = 3*dx*dy/r5;
                        weight[8] = 3*dy*dz/r5;
                        weight[9] = 3*dz*dx/r5;

                        // Calculate dot product of multipole array and weight array
                        double dotProduct = 0.0;
                        for(int d = 0; d < weight.length; d++){
                            double[] multipoleArray = cells.get(c).getMultipole();
                            dotProduct = dotProduct+multipoleArray[d]*weight[d];
                        }

                        particles.get(i).addToPhi(dotProduct);
                    }
                } else{ // Leaf Cell
                    // Loop in twig cell's particles
                    for(int j = 0; j < cells.get(p).getNumLeaves(); j++){
                        OctreeParticle source = particles.get(cells.get(p).getLeavesValueAtIndex(j));
                        double r = particles.get(i).distance(source);
                        if(r != 0){
                            particles.get(i).addToPhi(source.getCharge()/r);
                        }
                    }
                }
            }
        }
    }

    /**
     * Evaluate potential at all target points
     */
    public void evalPotnetial(){
        for(int i = 0; i < particles.size(); i++){
            evalAtTarget(0,i);
        }
    }

    public void l2Error(double[] phiDirect, double[] phiTree){
        double errorSumNum = 0.0;
        double errorSumDenom = 0.0;
        for(int i = 0; i < phiDirect.length; i++){
            errorSumNum = errorSumNum + Math.pow((phiDirect[i] - phiTree[i]),2);
            errorSumDenom = errorSumDenom + Math.pow(phiDirect[i],2);
        }
        double error = Math.sqrt(errorSumNum/errorSumDenom);
        logger.info("L2 Norm Error: "+error);
    }

    public void readParticle(File file){}

}
