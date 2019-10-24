package ffx.potential.nonbonded.implicit;

import org.apache.commons.lang3.BooleanUtils;

import java.util.ArrayList;
import java.util.logging.Logger;

public class Octree {

    private static final Logger logger = Logger.getLogger(Octree.class.getName());

    /**
     * Reference to one of the eight divisions in 3D
     */
    private int octant = 0;
    /**
     * List of cells
     */
    private ArrayList<OctreeCell> cells = new ArrayList<>();
    /**
     * Index of parent cell in cells list
     */
    private int p = 0;
    /**
     * Index of child cell in cells list
     */
    private int c = 0;
    /**
     * Critical (maximum allowed+1) number of points allowed in any one cell:
     * If a cell contains nCritical points, it needs to be split
     */
    private int nCritical = 0;
    /**
     * List of particles
     */
    private ArrayList<OctreeParticle> particles = new ArrayList<>();
    /**
     * Root cell
     */
    private OctreeCell root;
    /**
     * List of all leaf cells
     */
    private ArrayList<OctreeCell> leaves = new ArrayList<>();
    /**
     * Index of target particle
     */
    private int targetIndex = 0;
    /**
     * Tolerance parameter
     */
    private double theta = 0.5;

    /**
     * Constructor
     */
    public Octree(){
    }

    public void addChild(int octant, int p, ArrayList<OctreeCell> cells, int numCritical){
        OctreeCell tempCell = new OctreeCell(numCritical);

        // Create new cell
        // TODO: should this be passed in or the "this" value from the particular Octree?
        cells.add(tempCell);

        // Last element of cells list is new child, c
        c = cells.size() - 1;

        // Geometric reference between parent and child
        cells.get(c).setR((cells.get(p).getR())*0.5);
        cells.get(c).setX((cells.get(p).getX())*((octant & 1)*2-1));
        cells.get(c).setY((cells.get(p).getY())*((octant & 2)-1));
        cells.get(c).setZ((cells.get(p).getZ())*((octant & 4)/2-1));

        // Establish mutual reference in cells list
        cells.get(c).setParentIndex(p);
        cells.get(c).setChildren(octant,c);
        cells.get(c).setNumChildren(cells.get(p).getNumChildren() | (1 << octant));
    }

    public void splitCell(ArrayList<OctreeParticle> particles, int p, ArrayList<OctreeCell> cells, int nCritical){
        for (int i = 0; i < nCritical; i++){
            int octX = 0;
            int octY = 0;
            int octZ = 0;

            if(particles.get(i).getX() > cells.get(p).getX()){octX = 1;}
            if(particles.get(i).getY() > cells.get(p).getY()){octY = 1;}
            if(particles.get(i).getZ() > cells.get(p).getZ()){octZ = 1;}

            // Find particle's octant - should be an integer from 0 to 7
            octant = octX + (octY << 1) + (octZ << 2);

            // If there's not a child cell in the particle's octant, create one
            boolean noChildInOctant = BooleanUtils.toBoolean(cells.get(p).getNumChildren() & (1 << octant));
            if(noChildInOctant){
                addChild(octant,p,cells,nCritical);
            }

            // Reallocate the particle in the child cell
            c = cells.get(p).getChildAtIndex(octant);
            cells.get(c).setLeaf(cells.get(c).getNumLeaves(),1);
            cells.get(c).setNumLeaves(cells.get(c).getNumLeaves() + 1);

            // Check if child cell reaches nCritical - split recursively if so
            if(cells.get(c).getNumLeaves() >= nCritical){
                splitCell(particles,c,cells,nCritical);
            }

        }
    }

    public ArrayList<OctreeCell> buildTree(ArrayList<OctreeParticle> particles, OctreeCell root, int nCritical){
        ArrayList<OctreeCell> cells = new ArrayList<>();

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

                if(particles.get(i).getX() > cells.get(p).getX()){octX = 1;}
                if(particles.get(i).getY() > cells.get(p).getY()){octY = 1;}
                if(particles.get(i).getZ() > cells.get(p).getZ()){octZ = 1;}

                // Find particle's octant - should be an integer from 0 to 7
                octant = octX + (octY << 1) + (octZ << 2);

                // If there's not a child cell in the particle's octant, create one
                boolean noChildInOctant = BooleanUtils.toBoolean(cells.get(p).getNumChildren() & (1 << octant));
                if(noChildInOctant){
                    addChild(octant,p,cells,nCritical);
                }

                current = cells.get(current).getChildAtIndex(octant);
            }

            // Allocate the particle in the leaf cell
            cells.get(current).setLeaf(cells.get(current).getNumLeaves(),i);
            cells.get(current).setNumLeaves(cells.get(current).getNumLeaves() + 1);

            // Check whether to split cell
            if(cells.get(current).getNumLeaves() >= nCritical){
                splitCell(particles,current,cells,nCritical);
            }
        }

        return cells;
    }

    public void getMultipole(ArrayList<OctreeParticle> particles,int p,ArrayList<OctreeCell> cells,int[] leaves,int nCritical){
        // If the current cell is not a leaf cell, traverse down
        if(cells.get(p).getNumLeaves() >= nCritical){
            for(int c = 0; c < 8; c++){
                if(BooleanUtils.toBoolean(cells.get(p).getNumChildren() & (1 << c))){
                    getMultipole(particles,cells.get(p).getChildAtIndex(c),cells,leaves,nCritical);
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
                // leaves.add(p)
            }
        }

    }

    public void M2M(int p, int c, ArrayList<OctreeCell> cells){
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

    public void upwardSweep(ArrayList<OctreeCell> cells){
        for(int c = cells.size(); c > 0;c++){
            p = cells.get(c).getParentIndex();
            M2M(p,c,cells);
        }
    }

    public void directSum(ArrayList<OctreeParticle> particles){
        for(int i = 0; i < particles.size();i++){
            for(int j = 0; j < particles.size();j++){
                if(j!=i){
                    double r = particles.get(i).distance(particles.get(j));
                    particles.get(j).setPhi(particles.get(j).getCharge()/r);
                }
            }
        }
    }

    public void evalAtTarget(ArrayList<OctreeParticle> particles, int p, int i, ArrayList<OctreeCell> cells, int nCritical, double theta){

        // Non-leaf cell
        if(cells.get(p).getNumLeaves() >= nCritical){

            // Loop through p's child cells (8 octants)
            for(int oct = 0; oct < 8; oct++){
                if(BooleanUtils.toBoolean(cells.get(p).getNumChildren() & (1<<oct))){
                    int c = cells.get(p).getChildAtIndex(oct);
                    double r = 0.0;//particles.get(i).distance(cells.get(c));

                    // Near field child cell
                    if(cells.get(c).getR() > theta*r){
                        evalAtTarget(particles,c,i,cells,nCritical,theta);
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

                        particles.get(i).setPhi(dotProduct);
                    }
                } else{ // Leaf Cell
                    // Loop in twig cell's particles
                    for(int j = 0; j < cells.get(p).getNumLeaves(); j++){
                        //OctreeParticle source = particles.get(cells.get(p).getLea)
                    }
                }
            }
        }
    }

}
