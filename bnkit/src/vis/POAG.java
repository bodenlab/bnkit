/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package poag;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import api.PartialOrderGraph;

/**
 *
 * @author ariane
 */
public class POAG {
    
    PartialOrderGraph poag; 
    PathGen pathGen;
    HashMap<Integer, Node> nodes;
    int x;
    ArrayList<Integer> totalNodeList;
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        //"/home/ariane/Downloads/testPOGraphSmall.dot"; //"
        ///"home/ariane/Documents/bodenlab/bnkit/bnkit/src/test/resources/testPOGraphMed.dot"
        //testPOGraphMSAMed.dot
        //testPOGraphMed.dot
        //testPOGraphLarge.dot
        String dotPath = "/home/ariane/Documents/bodenlab/bnkit/bnkit/src/test/resources/testPOGraphMSAMed.dot"; //"/home/ariane/NetBeansProjects/POAG/src/poag/testPOGraphLarge.dot"; //new PartialOrderGraph("/home/ariane/Documents/stemformatics/biojs_alignment/data/defaultMSA.dot");
        POAG pg = new POAG(dotPath);

        HashMap<Integer, List<Integer>> pathDepths = pg.getPathDepths();

    }
    
    
    private void printNodes() {
        // BE CAREFUL THIS COULD BE BAD LATER ON
        int x = 0;
        Set<Integer> ns = nodes.keySet();
        for (int i: ns) {
            nodes.get(i).setX(x);
            x ++;
            System.out.println(nodes.get(i).toStr());
        }
    }
    public POAG (String dotPath) {
        nodes = new HashMap<>();
        poag = new PartialOrderGraph(dotPath);
        pathGen = new PathGen(poag);
        
    }
    
    /**
     * Returns a map of paths and the depth at which these occur
     * @return paths
     */
    public HashMap<Integer, List<Integer>>  getPathDepths() {
        Integer[] nodeys = poag.getNodeIDs();
        int numNodes = nodeys.length;
        x = numNodes;
        HashMap<Integer, List<Integer>> paths = new HashMap<>();
        int depth = 0;
         // want there to be a distinct x position for each node
        // This gets the main path which will be centered
        pathGen.resetSearchList();
        pathGen.initAStarSearch(pathGen.startID, pathGen.goalID);
        List<Integer> path = pathGen.getMainPath();
        paths.put(depth, path);
        addNodes(path, depth);
        // Now we need to get each of the subsequent paths
        // A depth will be associated for each itteration - this will be
        // visualised as further away from the central line.
        depth ++;
        
        int prevGapStart = 10000000; // Something large (note this needs to be done better)
        int gapEnd = 0;
        int gapStart = 0;
        while (!pathGen.gaps.isEmpty()) {
            // Get the next gap from the FIFO gap queue set up in the POAG
            Integer[] gap = pathGen.gaps.remove();
            gapStart = gap[0];
            gapEnd = gap[1];
            
            // Set up the AStar search environment with the start and end 
            // This also clears the searched nodes in the search map
            pathGen.resetSearchList();
            pathGen.initAStarSearch(gapStart, gapEnd);
            path = pathGen.getMainPath();
            
            // If the path was null we don't want to add it - see code below
            // for alternative method (ugly with nested while loops to step
            // through the gap.
            if (path == null) {
                //System.err.println("depth: " + depth + " Path: was none for gap: " + Arrays.toString(gap));
                continue;
            }
            //System.err.println("Depth: " + depth + ", Path: " + path);  
            paths.put(depth, path);
            addNodes(path, depth);
            // Check if we have reached the end of the first itteration of gaps
            // We can tell this because the gaps will start again from the 
            // End of the path.
            if (gapStart > prevGapStart) {
                depth ++;
            }
            
            prevGapStart = gapStart;
        }
        printNodes();
        // Reset the nodes to have correct x coords
        POAGJson pog = new POAGJson(poag, nodes);
        pog.toJSON();
        return paths;
    }
    
    
    private void addNodes(List<Integer> path, int depth) {
        int id;
       
        for (int i: path) {
            //System.err.println(i);
            // Sanity check that we aren't trying to add duplicate nodes
            if (x < 0) {
                //System.err.println("Error: x < 0");
                return;
            }
            id = i;
            Node n = nodes.get(id);
            if (n == null) {
                // if we don't have the node already then we want to add it
                n = new Node(id, x, depth, poag.getCharacterDistribution(id), poag.getOutEdgeWeights(id), poag.getSeqChars(id));
                x --; // Need to take one away from x
            } else { // Need to reset the x coord of it
                n.setX(x);
            }
            // Otherwise, check if the depth is smaller for this one as we always
            // want the smallest depth, keep the origional x
            if (n.getY() > depth) {
                n.setY(depth);
            }

            nodes.put(id, n);
        }
        
    }
    
   
    
    
    
    
}


