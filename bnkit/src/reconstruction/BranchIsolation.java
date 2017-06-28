package reconstruction;

import api.PartialOrderGraph;
import bn.prob.EnumDistrib;
import dat.POGraph;
import dat.PhyloTree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 *
 *  Check to see if information from a reconstruction is isolated to one branch in the phylogenetic tree
 *
 *  @author Gabe
 *
 */
public class BranchIsolation {

    private ASRPOG asr;
    private String treePath;
    private String sequencePath;
    private boolean mp;

    public BranchIsolation(ASRPOG asr, String treePath, String sequencePath, String node, boolean mp) throws IOException {
        this.asr = asr;
        this.treePath = treePath;
        this.sequencePath = sequencePath;
        this.mp = mp;

        POGraph msa = asr.getMSAGraph();


        //TODO: Can we get children for a node without having to do this conversion to a List
        // Get the children of the node of interest
        Collection<PhyloTree.Node> children = asr.getChildren(node);
        List<PhyloTree.Node> childrenList = new ArrayList(children);

        EnumDistrib[] leftDistrib = getDistrib(childrenList.get(0).getLabel().toString());
        EnumDistrib[] rightDistrib = getDistrib(childrenList.get(1).getLabel().toString());


        // Iterate through the MSA distributions
        for (int k = 0; k < msa.getNumNodes(); k++) {
            msa.setCurrent(msa.getNodeIDs().get(k));
            Map<Character, Double> distribution = msa.getCharacterDistribution();

            int currentNode = msa.getNodeIDs().get(k);

            System.out.println("At node " + currentNode);

            for (Character character : distribution.keySet()){


                System.out.println(character + " " +  distribution.get(character));
                System.out.println("Left distrib at this node is " + leftDistrib[currentNode].get(getDistibPostion(character)));
                System.out.println("Right distrib at this node is " + rightDistrib[currentNode].get(getDistibPostion(character)));

            }

        }

    }

    private EnumDistrib[] getDistrib(String node) throws IOException {

        ASRPOG child = new ASRPOG(sequencePath, treePath, sequencePath, node, mp);
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
