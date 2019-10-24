package ffx.potential.nonbonded.implicit;

public class OctreeCell {

    /**
     * Critical (maximum allowed+1) number of points allowed in any one cell:
     * If a cell contains nCritical points, it needs to be split
     */
    private int nCritical;
    /**
     * Number of leaves
     */
    private int numLeaves = 0;
    /**
     * Array of leaf indices
     */
    private int[] leaves = new int[nCritical];
    /**
     * Number of children
     */
    private int numChildren = 0;
    /**
     * Array of child indices, length 8
     */
    private int[] children = new int[]{0,0,0,0,0,0,0,0};
    /**
     * Parent cell index
     */
    private int parentIndex = 0;

    /**
     * Coordinates for the center of the cell
     */
    private double x = 0.0;
    private double y = 0.0;
    private double z = 0.0;

    /**
     * Radius of cell
     */
    private double r = 0.0;

    /**
     * Multipole array, length 10
     */
    private double[] multipole = new double[]{0,0,0,0,0,0,0,0,0,0};


    public OctreeCell(int nCritical){
        setnCritical(nCritical);
    }

    public void setnCritical(int nCrit){
        this.nCritical = nCrit;
    }

    /**
     * Returns cell radius
     * @return
     */
    public double getR(){ return this.r; }
    /**
     * Sets cell radius
     * @param r
     */
    public void setR(double r){ this.r = r; }
    /**
     * Gets x coordinate of center of cell
     * @return
     */
    public double getX(){ return this.x; }
    /**
     * Sets x coordinate of center of cell
     * @param x
     */
    public void setX(double x){ this.x = x; }
    /**
     * Gets y coordinate of center of cell
     * @return
     */
    public double getY(){ return this.y; }
    /**
     * Sets y coordinate of center of cell
     * @param y
     */
    public void setY(double y){ this.y = y; }
    /**
     * Gets z coordinate of center of cell
     * @return
     */
    public double getZ(){ return this.z; }
    /**
     * Sets z coordinate of center of cell
     * @param z
     */
    public void setZ(double z){ this.z = z; }

    public void setParentIndex(int p){
        this.parentIndex = p;
    }

    public void setChildren(int octant, int c){
        this.children[octant] = c;
    }

    public int getChildAtIndex(int octant){
        return children[octant];
    }

    public int getNumChildren(){
        return this.numChildren;
    }

    public void setNumChildren(int num){
        this.numChildren = num;
    }

    public int getNumLeaves(){ return this.numLeaves; }
    public void setNumLeaves(int num){ this.numLeaves = num; }

    public void setLeaf(int index, int leaf){
        this.leaves[index] = leaf;
    }
    public int getLeavesValueAtIndex(int index){
        return this.leaves[index];
    }

    public double[] getMultipole(){
        return this.multipole;
    }

    public void addToMultipole(double[] calculatedMultipole){
        for(int i = 0; i < 10; i++){
            this.multipole[i] += calculatedMultipole[i];
        }
    }

    public int getParentIndex(){
        return this.parentIndex;
    }
}
