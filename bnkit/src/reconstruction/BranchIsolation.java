package reconstruction;

import bn.prob.EnumDistrib;
import dat.POGraph;
import dat.PhyloTree;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 *
 *  Check to see if information from a reconstruction is isolated to one branch in the phylogenetic tree
 *
 *  @author Gabe
 *
 */
class BranchIsolation {

    private ASRPOG asr;
    private Map<String, String> ancestralDict;
    private String treePath;
    private String sequencePath;
    private boolean mp;
    private String model;
    private int threads;

    BranchIsolation(ASRPOG asr, Map<String, String> ancestralDict, String treePath, String sequencePath, String node, boolean mp, String model, int threads) throws IOException, InterruptedException {
        this.asr = asr;
        this.ancestralDict = ancestralDict;
        this.treePath = treePath;
        this.sequencePath = sequencePath;
        this.mp = mp;
        this.model = model;
        this.threads = threads;



        List<Integer> unsupportedLeft = new ArrayList<>();
        List<Integer> unsupportedRight = new ArrayList<>();

        POGraph msa = asr.getMSAGraph();



        //TODO: Can we get children for a node without having to do this conversion to a List
        // Get the children of the node of interest
        Collection<PhyloTree.Node> children = asr.getChildren(node);
        List<PhyloTree.Node> childrenList = new ArrayList(children);


        String leftNode = childrenList.get(0).getLabel().toString();
        String rightNode = childrenList.get(1).getLabel().toString();


        // Get the marginal distributions for each child
        if (ancestralDict.get(leftNode) != null && ancestralDict.get(rightNode) != null) {

//        }


            EnumDistrib[] parentDistrib = getDistrib(node);
            EnumDistrib[] leftDistrib = getDistrib(childrenList.get(0).getLabel().toString());
            EnumDistrib[] rightDistrib = getDistrib(childrenList.get(1).getLabel().toString());


            // Iterate through the MSA distributions
            for (int k = 0; k < msa.getNumNodes(); k++) {
                msa.setCurrent(msa.getNodeIDs().get(k));
                Map<Character, Double> distribution = msa.getCharacterDistribution();

                int currentNode = msa.getNodeIDs().get(k);


                Character maxMSACharacter = Collections.max(distribution.entrySet(), Map.Entry.comparingByValue()).getKey();

                Map<String, List<Inference>> ancestralInferences = asr.getAncestralInferences();
                Character maxDistribCharacter = ancestralInferences.get(node).get(k).base;
                Character rootDistribCharacter = ancestralInferences.get("N0").get(k).base;

                if (maxDistribCharacter != '-' && rootDistribCharacter != '-') {


                    //TODO: Add check for Distrib even existing here

                    if (leftDistrib[currentNode] == null) {
                        System.out.println(childrenList.get(0).getLabel().toString() + " doesn't have a graph node at position " + k);
                    }

                    if (rightDistrib[currentNode] == null) {
                        System.out.println(childrenList.get(1).getLabel().toString() + " doesn't have a graph node at position " + k);

                    }

                    if (leftDistrib[currentNode] != null && rightDistrib[currentNode] != null) {

                        Character leftDistribMax = getDistribCharacter(leftDistrib[currentNode].getMaxIndex());
                        Character rightDistribMax = getDistribCharacter(rightDistrib[currentNode].getMaxIndex());


                        double leftCharacterDistrib = leftDistrib[currentNode].get(getDistibPostion(maxDistribCharacter));
                        double rightCharacterDistrib = rightDistrib[currentNode].get(getDistibPostion(maxDistribCharacter));

                        double leftMSADistrib = leftDistrib[currentNode].get(getDistibPostion(maxMSACharacter));
                        double rightMSADistrib = rightDistrib[currentNode].get(getDistibPostion(maxMSACharacter));

                        double leftRootDistrib = leftDistrib[currentNode].get(getDistibPostion(rootDistribCharacter));
                        double rightRootDistrib = rightDistrib[currentNode].get(getDistibPostion(rootDistribCharacter));

                        boolean foundDistrib = false;
                        boolean foundMSA = false;

                        if (leftCharacterDistrib < 0.5 && rightCharacterDistrib > 0.5) {
                            System.out.println("*********YEPPPERS********");
                            foundDistrib = true;
                        } else if (rightCharacterDistrib < 0.5 && leftCharacterDistrib > 0.5) {
                            System.out.println("***********YEPPERS********");
                            foundDistrib = true;
                        }

                        if (leftMSADistrib < 0.5 && rightMSADistrib > 0.5) {
                            System.out.println("*********MSA DISTRIB********");
                            foundMSA = true;
                        } else if (rightMSADistrib < 0.5 && leftMSADistrib > 0.5) {
                            System.out.println("***********MSA DISTRIB********");
                            foundMSA = true;
                        }

                        if (foundDistrib || foundMSA) {

//                if (leftCharacterDistrib < 0.1 && rightCharacterDistrib > 0.5) {
//                    unsupportedLeft.add(k);
                            System.out.println();
                            System.out.println("We are at node " + node + " in the phylogenetic tree");
                            System.out.println("Node in the partial order graph is " + k);
                            System.out.println("Node distrib character we're checking is " + maxDistribCharacter);
                            System.out.println("Root distrib character is " + rootDistribCharacter);
                            System.out.println("MSA highest character is " + maxMSACharacter);

                            System.out.println("Left distrib max character is " + leftDistribMax);
                            System.out.println("Right distrib max character is " + rightDistribMax);


                            System.out.println("MSA distribution of max character here is " + distribution.get(maxMSACharacter));
//                    System.out.println("MSA distribution of character we're checking is " + distribution.get(maxMSACharacter));

                            System.out.println("Child node " + leftNode + " in the phylogenetic tree has distribution " + leftCharacterDistrib + " of the max parent character " + maxDistribCharacter);
                            System.out.println("Child node " + rightNode + " in the phylogenetic tree has distribution " + rightCharacterDistrib + " of the max parent character " + maxDistribCharacter);

                            System.out.println("Child node " + leftNode + " in the phylogenetic tree has distribution " + leftMSADistrib + " of the  max MSA character " + maxMSACharacter);
                            System.out.println("Child node " + rightNode + " in the phylogenetic tree has distribution " + rightMSADistrib + " of the max MSA character " + maxMSACharacter);

                            System.out.println("Child node " + leftNode + " in the phylogenetic tree has distribution " + leftRootDistrib + "of the root character " + rootDistribCharacter);
                            System.out.println("Child node " + rightNode + " in the phylogenetic tree has distribution " + rightRootDistrib + "of the root character " + rootDistribCharacter);
                            foundDistrib = false;
                            foundMSA = false;

                        }
                    }

                }

            }
        }

        else {
            System.out.println("At least one of the children of this node is an extant");
        }

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

    private EnumDistrib[] getDistrib(String node) throws IOException, InterruptedException {




            //TODO: Making a marginal reconstruction at an extant node should throw an error

        System.out.println("\nCalculating a marginal reconstruction at " + node);

        final long startTime = System.nanoTime();


        ASRPOG child = new ASRPOG(sequencePath, treePath, sequencePath, node, mp, model, threads);

        final long duration = System.nanoTime() - startTime;

        System.out.println("Marginal reconstruction at  " + node + " complete");

        System.out.println("Marginal reconstruction took " + TimeUnit.NANOSECONDS.toSeconds(duration) + " seconds");

        return child.getMarginalDistributions();




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

    private Character getDistribCharacter(int pos){

        switch(pos) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'D';
            case 3:
                return 'E';
            case 4:
                return 'F';
            case 5:
                return 'G';
            case 6:
                return 'H';
            case 7:
                return 'I';
            case 8:
                return 'K';
            case 9:
                return 'L';
            case 10:
                return 'M';
            case 11:
                return 'N';
            case 12:
                return 'P';
            case 13:
                return 'Q';
            case 14:
                return 'R';
            case 15:
                return 'S';
            case 16:
                return 'T';
            case 17:
                return 'V';
            case 18:
                return 'W';
            case 19:
                return 'Y';
            default:
                System.out.println("Coudn't map this position to character: " + pos);
        }

        return 'X';

    }
}
