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

    public BranchIsolation(ASRPOG asr, String treePath, String sequencePath, String node, boolean mp) throws IOException {

        POGraph msa = asr.getMSAGraph();

        // Iterate through the MSA distributions
        System.out.println("Print all the MSA distributions");

        for (int k = 0; k < msa.getNumNodes(); k++) {
            msa.setCurrent(msa.getNodeIDs().get(k));
            Map<Character, Double> distribution = msa.getCharacterDistribution();
            System.out.println(msa.getNodeIDs().get(k));
            System.out.println(distribution);

        }



        //TODO: Can we get children for a node without having to do this conversion to a List
        // Get the children of the node of interest
        Collection<PhyloTree.Node> children = asr.getChildren(node);
        List<PhyloTree.Node> childrenList = new ArrayList(children);

        String leftNode = childrenList.get(0).getLabel().toString();
        String rightNode = childrenList.get(1).getLabel().toString();



        //TODO: Is it best to recreate ASRPOG objects using the sequencePath / treePath ?
        // Create marginal reconstructions for the children
        ASRPOG leftChild = new ASRPOG(sequencePath, treePath, sequencePath, leftNode, mp);
        ASRPOG rightChild = new ASRPOG(sequencePath, treePath, sequencePath, rightNode.toString(), mp);

        EnumDistrib[] leftDistrib = leftChild.getMarginalDistributions();
        EnumDistrib[] rightDistrib = rightChild.getMarginalDistributions();

        POGraph leftGraph = leftChild.getGraph(leftNode);
        POGraph rightGraph = rightChild.getGraph(rightNode);
        

    }
}
