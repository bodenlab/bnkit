package reconstruction;

import api.PartialOrderGraph;
import bn.prob.EnumDistrib;
import dat.POGraph;
import dat.PhyloTree;

import java.io.IOException;
import java.security.Key;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 *
 *  Check to see if information from a reconstruction is isolated to one branch in the phylogenetic tree
 *
 *  @author Gabe
 *
 */
public class BranchIsolation {

    private ASRPOG asr;
    private Map<String, String> ancestralDict;
    private String treePath;
    private String sequencePath;
    private boolean mp;

    public BranchIsolation(ASRPOG asr, Map<String,String> ancestralDict, String treePath, String sequencePath, String node, boolean mp) throws IOException {
        this.asr = asr;
        this.ancestralDict = ancestralDict;
        this.treePath = treePath;
        this.sequencePath = sequencePath;
        this.mp = mp;

        List<Integer> unsupportedLeft = new ArrayList<>();
        List<Integer> unsupportedRight = new ArrayList<>();

        POGraph msa = asr.getMSAGraph();

        System.out.println("\nCurrently at node: " + node);


        //TODO: Can we get children for a node without having to do this conversion to a List
        // Get the children of the node of interest
        Collection<PhyloTree.Node> children = asr.getChildren(node);
        List<PhyloTree.Node> childrenList = new ArrayList(children);

        String leftNode = childrenList.get(0).getLabel().toString();
        String rightNode = childrenList.get(1).getLabel().toString();

        // Get the marginal distributions for each child
        if (ancestralDict.get(leftNode) != null && ancestralDict.get(rightNode) != null) {

            System.out.println("This node is a goer " + node);
        }


//            EnumDistrib[] leftDistrib = getDistrib(childrenList.get(0).getLabel().toString());
//            EnumDistrib[] rightDistrib = getDistrib(childrenList.get(1).getLabel().toString());
//
//
//
//            // Iterate through the MSA distributions
//            for (int k = 0; k < msa.getNumNodes(); k++) {
//                msa.setCurrent(msa.getNodeIDs().get(k));
//                Map<Character, Double> distribution = msa.getCharacterDistribution();
//
//                int currentNode = msa.getNodeIDs().get(k);
//
//                Character maxMSACharacter = Collections.max(distribution.entrySet(), Map.Entry.comparingByValue()).getKey();
//
//                double leftCharacterDistrib = leftDistrib[currentNode].get(getDistibPostion(maxMSACharacter));
//                double rightCharacterDistrib = rightDistrib[currentNode].get(getDistibPostion(maxMSACharacter));
//
//
//
//
//                if (leftCharacterDistrib < 0.1 && rightCharacterDistrib > 0.5) {
//                    unsupportedLeft.add(k);
//                    System.out.println();
//                    System.out.println("We are at node " + node + " in the phylogenetic tree");
//                    System.out.println("Node in the partial order graph is " + k);
//                    System.out.println("Character we're checking is " + maxMSACharacter);
//                    System.out.println("MSA distribution here is " + distribution.get(maxMSACharacter));
//                    System.out.println("Child node " + leftNode + " in the phylogenetic tree has distribution " + leftCharacterDistrib);
//                    System.out.println("Child node " + rightNode + " in the phylogenetic tree has distribution " + rightCharacterDistrib);
//                }
//
//                if (rightCharacterDistrib < 0.1 && leftCharacterDistrib > 0.5) {
//                    unsupportedRight.add(k);
//                    System.out.println();
//                    System.out.println("We are at node " + node + " in the phylogenetic tree");
//                    System.out.println("Node in the partial order graph is " + k);
//                    System.out.println("Character we're checking is " + maxMSACharacter);
//                    System.out.println("MSA distribution here is " + distribution.get(maxMSACharacter));
//                    System.out.println("Child node " + leftNode + " in the phylogenetic tree has distribution " + leftCharacterDistrib);
//                    System.out.println("Child node " + rightNode + " in the phylogenetic tree has distribution " + rightCharacterDistrib);
//                }
//
//
//            }
//
//        }
//
//        else {
//            System.out.println("At least one of the children of this node is an extant");
//        }

        if (!unsupportedLeft.isEmpty() || !unsupportedRight.isEmpty()) {

            System.out.println("\nLeft is unsupported at ");
            for (int pos : unsupportedLeft) {
                if (!unsupportedRight.contains(pos)) {
                    System.out.println(pos);
                }
            }

            System.out.println("\nRight is unsupported at ");
            for (int pos : unsupportedRight) {
                if (!unsupportedLeft.contains(pos)){
                    System.out.println(pos);
                }

            }
        }



    }

    private EnumDistrib[] getDistrib(String node) throws IOException {




            //TODO: Making a marginal reconstruction at an extant node should throw an error

        System.out.println("\nCalculating a marginal reconstruction at " + node);

        final long startTime = System.nanoTime();

        ASRPOG child = new ASRPOG(sequencePath, treePath, sequencePath, node, mp);

        final long duration = System.nanoTime() - startTime;




        System.out.println("Marginal reconstruction at  " + node + " complete");

        System.out.println("Marginal reconstruction took " + TimeUnit.NANOSECONDS.toSeconds(duration) + " seconds");

        EnumDistrib[] distrib = child.getMarginalDistributions();







        return distrib;



    }

    private int getDistibPostion(Character character){

        switch(character) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'D':
                return 2;
            case 'E':
                return 3;
            case 'F':
                return 4;
            case 'G':
                return 5;
            case 'H':
                return 6;
            case 'I':
                return 7;
            case 'K':
                return 8;
            case 'L':
                return 9;
            case 'M':
                return 10;
            case 'N':
                return 11;
            case 'P':
                return 12;
            case 'Q':
                return 13;
            case 'R':
                return 14;
            case 'S':
                return 15;
            case 'T':
                return 16;
            case 'V':
                return 17;
            case 'W':
                return 18;
            case 'Y':
                return 19;
            default:
                System.out.println("Coudn't map this character to a marginal distiburtion: " + character);
        }

        return -1;

    }
}
