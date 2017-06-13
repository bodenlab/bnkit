//package alignment;
//
//
//import alignment.utilities.MatrixUtils;
////import SubstitutionModels.SubstitutionMatrix;
////import reconstruction.SubstitutionMatrix;
//import alignment.utilities.SubstitutionMatrix;
//import alignment.utilities.HashProfile;
//import java.util.*;
//
//import static java.lang.Math.max;
//
//
//public class Alignment {
//
//    /**
//     * Class for aligning a sequence to a PO Graph
//     */
//
//    private int openGapPenalty;
//    private int extendGapPenalty;
//    private HashProfile profile1;
//    private HashProfile profile2;
//    //    private List<Integer> stringIndexes;
////    private List<Integer> nodeIndexes;
//    private SubstitutionMatrix subMatrix;
//    private boolean MEA;
//    private HashProfile updatedProfile;
//
//    private double[][] scores;
//    private int[][] iBackIndexes;
//    private int[][] jBackIndexes;
//
//
//
//    /**
//     * Constructor that takes sequences as Strings and converts them into HashProfiles
//     * @param seq1 First sequence
//     * @param seq2 Second sequence
//     * @param openGapPenalty Penalty for opening gap
//     * @param extendGapPenalty Penalty for extending gap
//     * @param subMatrix Substitution matrix
//     * @param MEA Boolean representing whether we are performing Maximum Expected Accuracy alignment
//     */
//    public Alignment(String seq1, String seq2, int openGapPenalty, int extendGapPenalty, SubstitutionMatrix subMatrix, boolean MEA) {
//
//        this(new HashProfile(seq1), new HashProfile(seq2), openGapPenalty, extendGapPenalty, subMatrix, MEA);
//    }
//
//    /**
//     * Constructor for aligning two profiles together
//     * @param profile1 First profile to align
//     * @param profile2 Second profile to align
//     * @param openGapPenalty Penalty for opening a gap
//     * @param extendGapPenalty Penalty for extending a gap
//     * @param subMatrix Substituion matrix
//     * @param MEA Boolean representing whether we are performing Maximum Expected Accuracy alignment
//     */
//    public Alignment(HashProfile profile1, HashProfile profile2, int openGapPenalty, int extendGapPenalty, SubstitutionMatrix subMatrix, boolean MEA) {
//        this.profile1 = profile1;
//        this.profile2 = profile2;
//        this.openGapPenalty = openGapPenalty;
//        this.extendGapPenalty = extendGapPenalty;
////        this.stringIndexes = new ArrayList<Integer>();
////        this.nodeIndexes = new ArrayList<Integer>();
//        this.subMatrix = subMatrix;
//        this.MEA = MEA;
//        this.updatedProfile = this.alignProfileToProfile();
//
//    }
//
//
//
//
//
//    /**
//     * @return HashProfile representing the two aligned profiles
//     */
//    public HashProfile getUpdatedProfile() {
//        return this.updatedProfile;
//    }
//
////    /**
////     * Return the scores associated with matching a single character in a PO Graph to an
////     * array of characters.
////     * Used to get the scores for matching each row in the score matrix to each column
////     *
////     * @param base   the character in the PO Graph to be matched against
////     * @param seqvec characters in the sequence to be matched against
////     * @return an int array containing the match scores
////     */
////    public double[] getMatchScore(char base, String seqvec) {
////        double[] matches = new double[seqvec.length()];
////        for (int i = 0; i < seqvec.length(); i++) {
////            double matchScore = subMatrix.getDistance(seqvec.charAt(i), base);
////            matches[i] = matchScore;
////        }
////
////        return matches;
////    }
//
////    public double[] getMEAMatchScore(int i, String seqVec) {
////        double[] matches = new double[seqVec.length()];
////        for (int j = 0; j < seqVec.length(); j++) {
////            double matchScore = matrix[i + 1][j + 1];
////            matches[j] = matchScore;
////
////        }
////
////        return matches;
////    }
//
//    /**
//     * Get the score for making a match in Maximum Expected Accuracy alignment
//     * @param i
//     * @param profile
//     * @return
//     */
//    public double[] getMEAMatchScore(int i, HashProfile profile) {
//        double[] matches = new double[profile.getProfileArray().size()];
//        System.arraycopy(subMatrix.getMatrix()[i + 1], 1, matches, 0, profile.getProfileArray().size());
//
//        return matches;
//    }
//
//    /**
//     *
//     * @param profile1 First profile
//     * @param profile2 Second profile
//     * @param index
//     * @return
//     */
//    public double[] getProfileMatchScore(HashProfile profile1, HashProfile profile2, int index) {
//
//        int profile2Length = profile2.getLength();
//
//        double[] matches = new double[profile2Length];
//
//
//        //TODO: Normalise match score at individual index
//        for (int i = 0; i < profile2Length; i++) {
//            double totalScore = 0;
//            double totalCount;
//            double profile1Count = 0;
//            double profile2Count = 0;
//
//            // Get total count of residues at this position in first profile
//            for (Character residue : profile1.getProfileArray().get(index).keySet()) {
//                profile1Count += profile1.getProfileArray().get(index).get(residue).getValue();
//            }
//
//            // Get total count of residues at this position in second profile
//            for (Character residue : profile2.getProfileArray().get(i).keySet()) {
//                profile2Count += profile2.getProfileArray().get(i).get(residue).getValue();
//            }
//
//            totalCount = profile1Count * profile2Count;
//
//            for (Character name : profile1.getProfileArray().get(index).keySet()) {
//                if (name != '-') {
//
//                    // Get count for specific residue
//                    int profile1Value = profile1.getProfileArray().get(index).get(name).getValue();
//
//                    for (Character name2 : profile2.getProfileArray().get(i).keySet()) {
//                        if (name2 != '-') {
//                            int profile2Value = profile2.getProfileArray().get(i).get(name2).getValue();
//
//
//                            // Get the score of matching these residues
//                            double matchScore = subMatrix.getDistance(name, name2);
//                            totalScore += profile1Value * profile2Value * matchScore;
//                            matches[i] = totalScore / totalCount;
//
//                        }
//                    }
//                }
//
//            }
//        }
//
//        return matches;
//    }
//
//
////    /**
////     * Get list of Node IDs representing current sequence after alignment, null for gaps in sequence
////     * @return list of Node IDS
////     */
////    public List<Integer> getStringIndexes() {
////        return this.stringIndexes;
////    }
////
////
////    /**
////     * Get list of Node IDs that new sequence is being matched to, null for no match
////     * @return list of Node IDS
////     */
////    public List<Integer> getNodeIndexes() {
////        return this.nodeIndexes;
////    }
////
//
////    /**
////     * Returns the indexes of all predecessor nodes
////     * @param node the node to get predecessors for
////     * @param nodeIDtoIndex a mapping of IDs to index positions
////     * @return list of indexes of all predecessor nodes
////     */
////
////    public List<Integer> getPrevIndices(Node node, Map<Integer, Integer> nodeIDtoIndex) {
////        List<Integer> prev = new ArrayList<Integer>();
////
////
////        List<Integer> edgeList = new ArrayList<Integer>();
////
////        for (Integer edge : node.getInEdges().keySet()) {
////            edgeList.add(edge);
////        }
////
////        for (Integer predID : edgeList) {
////            prev.add(nodeIDtoIndex.get(predID));
////        }
////
////        if (prev.size() == 0) {
////            prev.add(-1);
////        }
////
////        return prev;
////    }
//
//    /**
//     * Set up the initial sizes and values for the scores matrix, the matrices to
//     * record the optimal moves through the score matrix, and the two mappings from
//     * node ID to index and index to node ID
//     */
//    public void initialiseDyncamicProgrammingData() {
//        int l1 = this.profile1.getLength();
//        int l2 = this.profile2.getLength();
////        Map<Integer, Integer> nodeIDToIndex = new HashMap<Integer, Integer>();
////        Map<Integer, Integer> nodeIndexToID = new HashMap<Integer, Integer>();
//
////        for (Integer index : this.graph.getNodeDict().keySet()) {
////            nodeIDToIndex.put(this.graph.getNodeIDList().get(index), index);
////            nodeIndexToID.put(index, this.graph.getNodeIDList().get(index));
////        }
//
//        // Create score matrix
//        scores = new double[l1 + 1][l2 + 1];
//
//        // Fill top row with gap penalties
//        for (int i = 0; i < l2; i++) {
//            scores[0][i + 1] = openGapPenalty + (i * extendGapPenalty);
//        }
//
//
//        for (int i = 0; i < l1; ++i) {
//            scores[i + 1][0] = openGapPenalty + (i * extendGapPenalty);
//        }
//
//        // Array to hold the optimal i index positions
//        iBackIndexes = new int[l1 + 1][l2 + 1];
//
//        //Fill first row with optimal previous cell move
//        for (int i = 0; i < l1; i++) {
//            iBackIndexes[i + 1][0] = i;
//        }
//
//        // Array to hold the optimal j index positions
//        jBackIndexes = new int[l1 + 1][l2 + 1];
//
//        // Fill first column with optimal previous cell move
//        for (int i = 0; i < l2; i++) {
//            jBackIndexes[0][i + 1] = i;
//        }
//
//    }
//
//    /**
//     * Trace back through the scores matrix, building the optimal alignment
//     *
//     * @param scores Array representing the filled out dynamic programming matrix
//     * @param jBackIndexes Array representing the j-index to backtrack to
//     * @param iBackIndexes Array representing the i-index to backtrack to
//     * @return HashProfile representing the two aligned profiles
//     */
//
//    public HashProfile backtrack(double[][] scores, int[][] jBackIndexes, int[][] iBackIndexes) {
//
//
//        //TODO: Scores only recording one optimal path - see simpleProfile vs simplerProfile 19/08/2016
//
//        int besti = scores.length;
//        int bestj = scores[0].length;
//
//        besti -= 1;
//        bestj -= 1;
////        List<Integer> terminalIndices = new ArrayList<Integer>();
//
//
////        double bestScore = scores[besti][bestj];
//
////        for (int i = 0; i < terminalIndices.size(); i++) {
////            double score = scores[terminalIndices.get(i)][bestj];
////            if (score > bestScore) {
////                bestScore = score;
////                besti = terminalIndices.get(i);
////            }
////        }
//
//        List<Integer> profile1Indexes = new ArrayList<Integer>();
//        List<Integer> profile2Indexes = new ArrayList<Integer>();
//
//        int profile2Index = 0;
//        int profile1Index = 0;
//
//
//        while (!(besti == 0 && bestj == 0)) {
//            // Get the i and j coordinates of the cell to move to
//            int nexti = iBackIndexes[besti][bestj];
//            int nextj = jBackIndexes[besti][bestj];
//
//            profile2Index = bestj - 1;
//
////            besti = (besti == 0) ? 1 : besti;
//            profile1Index = besti - 1;
//
//            // If the j index to backtrack to is not the same as current j index, there can't be a gap in second profile
//            if (nextj != bestj) {
//                profile2Indexes.add(0, profile2Index);
//
//                // Else, there must be a gap in second profile at this position
//            } else {
//                profile2Indexes.add(0, -1);
//
//            }
//            // If the i index to backtrack to is not the same as current i index, there can't be a gap in first profile
//            if (nexti != besti) {
//                profile1Indexes.add(0, profile1Index);
//
//                // Else, there must be a gap in first profile at this position
//            } else {
//                profile1Indexes.add(0, -1);
//            }
//
//            besti = nexti;
//            bestj = nextj;
//        }
//
//        // Fill out the remaining indexes of each profile
//        while (profile1Indexes.get(0) > 0 && profile1Index > 0) {
//            profile1Indexes.add(0, profile1Index - 1);
//            profile2Indexes.add(0, -1);
//            profile1Index -= 1;
//
//        }
//
//        while (profile2Indexes.get(0) > 0 && profile2Index > 0) {
//            profile2Indexes.add(0, profile2Index - 1);
//            profile1Indexes.add(0, -1);
//            profile2Index -= 1;
//
//        }
//
//        // Update the profiles by adding the gaps to them
//        addGapsToProfiles(profile1Indexes, profile2Indexes);
//
//
//        // Return new profile alignment by joining two updated profiles together
//
//        return new HashProfile(profile1, profile2);
//
//    }
//
////    /**
////     * Align a new sequence to the PO Graph
////     *
////     * @return
////     */
////
////
////    public HashProfile alignSeqToGraph() {
////
////
////        int l1 = this.seq1.length();
////        int l2 = this.seq2.length();
////
////
////        List<Object> initialisedData = this.initialiseDyncamicProgrammingData(l1, l2);
////
////        double[][] scores = (double[][]) initialisedData.get(2);
////        int[][] backStrIdx = (int[][]) initialisedData.get(3);
////        int[][] backGraphIdx = (int[][]) initialisedData.get(4);
////        double[][] insertCost = new double[l1 + 1][l2 + 1];
////        double[][] deleteCost = new double[l1 + 1][l2 + 1];
////
////        for (double[] row : insertCost)
////            Arrays.fill(row, this.openGapPenalty);
////
////        for (double[] row : deleteCost)
////            Arrays.fill(row, this.openGapPenalty);
////
////
////        // Array to keep track of whether we should insert
////        boolean[] inserted = new boolean[l2 + 1];
////        Arrays.fill(inserted, false);
////
////
////        double[] insscores = new double[l2 + 2];
////
////        for (int i = 0; i < this.seq1.length(); i++) {
////
////            // Get character of node
////            char pbase = seq1.charAt(i);
////
////            // Get predecessors of node
//////            List<Integer> predecessors = this.getPrevIndices(this.graph.getNodeDict().get(this.graph.getNodeIDList().get(i)), nodeIDToIndex);
////
////
////            // Get array of scores for matching current node with each position in sequence
////            double[] matchPoints = new double[seq2.length()];
////
////            if (MEA) {
////                matchPoints = this.getMEAMatchScore(i, seq2);
////            } else {
////
////                matchPoints = this.getMatchScore(pbase, seq2);
////
////            }
////
////            // Get array of scores equal to previous row + the cost of a deletion at each position
////            double[] deleteScore = Arrays.copyOfRange(scores[i], 1, scores[0].length);
////            deleteScore = addArrays(deleteScore, Arrays.copyOfRange(deleteCost[i], 1, deleteCost[0].length));
////
////
////            // Array to store the best possible score for deleting at each position
////            int[] bestDelete = new int[l2];
////
////            // Fill array with first predecessor position
////            Arrays.fill(bestDelete, i);
////
////
////            double[] matchScore = addArrays(Arrays.copyOfRange(scores[i], 0, scores[0].length - 1), matchPoints);
////            int[] bestMatch = new int[l2];
////            Arrays.fill(bestMatch, i);
////
////
//////            for (int j = 1; j < predecessors.size(); j++) {
//////                int predecessor = predecessors.get(j);
//////                double[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);
//////
//////
//////                for (int k = 0; k < bestDelete.length; k++) {
//////                    int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
//////                    bestDelete[k] = bestDeleteElement;
//////
//////                    int bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
//////                    deleteScore[k] = bestDeleteScoreElement;
//////                }
//////
//////
//////                double[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);
//////
//////                for (int n = 0; n < bestMatch.length; n++) {
//////                    int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
//////                    bestMatch[n] = bestMatchElement;
//////
//////                    int bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
//////                    matchScore[n] = bestMatchScoreElement;
//////                }
//////            }
////
////
////            boolean[] deleted = new boolean[deleteScore.length];
////
////            for (int k = 0; k < deleteScore.length; k++) {
////                deleted[k] = (deleteScore[k]) >= matchScore[k];
////                scores[i + 1][k + 1] = max(deleteScore[k], matchScore[k]);
////            }
////
////
////            for (int l = 1; l < scores[0].length; l++) {
////                scores[i + 1][l] = max(deleteScore[l - 1], matchScore[l - 1]);
////
////                if (deleted[l - 1]) {
////                    backGraphIdx[i + 1][l] = bestDelete[l - 1];
////                    backStrIdx[i + 1][l] = l;
////                    deleteCost[i + 1][l] = extendGapPenalty;
////
////                } else {
////                    backGraphIdx[i + 1][l] = bestMatch[l - 1];
////                    backStrIdx[i + 1][l] = l - 1;
////                }
////            }
////
////            double[] scoreRow = scores[i + 1];
////            double[] insertRow = insertCost[i + 1];
////            int[] backStrRow = backStrIdx[i + 1];
////            Arrays.fill(inserted, false);
////            insscores = addArrays(scoreRow, insertRow);
////
////            for (int n = 0; n < l2; n++) {
////                if (insscores[n] >= scoreRow[n + 1]) {
////                    scores[i + 1][n + 1] = insscores[n];
////                    scoreRow[n + 1] = insscores[n];
////
////                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
////                    insertRow[n + 1] = this.extendGapPenalty;
////                    backStrIdx[i + 1][n + 1] = n;
////                    backStrRow[n + 1] = n;
////                    inserted[n + 1] = true;
////                    insscores[n + 1] = scores[i + 1][n + 1] + insertCost[i + 1][n + 1];
////                }
////            }
////
////            for (int o = 0; o < inserted.length; o++) {
////                if (inserted[o]) {
////                    insertCost[i + 1][o] = this.extendGapPenalty;
////                    deleteCost[i + 1][o] = this.openGapPenalty;
////                    backGraphIdx[i + 1][o] = i + 1;
////                }
////            }
////        }
//////        System.out.println("Scores matrix is ");
////
////        return backtrack(scores, backStrIdx, backGraphIdx);
////
////    }
//
//    /**
//     * Align a new sequence to the PO Graph
//     * @return Call to backtrack with the the filled out scores, jIndexMatrix, and iIndexMatrix
//     */
//    public HashProfile alignProfileToProfile() {
//
//        int l1 = this.profile1.getProfileArray().size();
//        int l2 = this.profile2.getProfileArray().size();
//
//
//        initialiseDyncamicProgrammingData();
//
//
//
////        List<Object> initialisedData = this.initialiseDyncamicProgrammingData(l1, l2);
////
////        double[][] scores = (double[][]) initialisedData.get(2);
////        int[][] jBackIndexes = (int[][]) initialisedData.get(3);
////        int[][] iBackIndexes = (int[][]) initialisedData.get(4);
//
//        // Arrays to keep track of insert / deletion history and cost
//        double[][] insertCost = new double[l1 + 1][l2 + 1];
//        double[][] deleteCost = new double[l1 + 1][l2 + 1];
//
//        for (double[] row : insertCost)
//            Arrays.fill(row, this.openGapPenalty);
//
//        for (double[] row : deleteCost)
//            Arrays.fill(row, this.openGapPenalty);
//
//
//        // Array to keep track of whether we should insert
//        boolean[] inserted = new boolean[l2 + 1];
//        Arrays.fill(inserted, false);
//
//
//        double[] insscores;
//
//        for (int i = 0; i < l1; i++) {
//
////            profile1.getProfileMatrix()[][i];
//
//
//            double[] matchPoints;
//
//            // If we're performing Maximum Expected Accuracy, get the match score from the posterior probability matrix
//            if (MEA) {
//                matchPoints = this.getMEAMatchScore(i, profile2);
//            }
//
//            // Otherwise get the match score from the chosen substitution matrix
//            else {
//
//                matchPoints = this.getProfileMatchScore(profile1, profile2, i);
//
//            }
//
//            System.out.println("Profile match score ");
//            for (double match : matchPoints){
//                System.out.println(match);
//            }
//
//
//            // Get character of node
////            char pbase = seq1.charAt(i);
//
//            // Get predecessors of node
////            List<Integer> predecessors = this.getPrevIndices(this.graph.getNodeDict().get(this.graph.getNodeIDList().get(i)), nodeIDToIndex);
//
//
//            // Get array of scores for matching current node with each position in sequence
////            double[] matchPoints = new double[seq2.length()];
//
////            if (MEA){
////                matchPoints = this.getMEAMatchScore(i, seq2);
////            }
////
////            else {
////
////                matchPoints = this.getMatchScore(pbase, seq2);
////
////            }
//
//            // Check if there should be a deletion in the move to current row
//
//            // Get array of scores equal to previous row + the cost of moving to a deletion in current row
//            double[] deleteScore = Arrays.copyOfRange(scores[i], 1, scores[0].length);
//            deleteScore = MatrixUtils.addArrays(deleteScore, Arrays.copyOfRange(deleteCost[i], 1, deleteCost[0].length));
//
//
//            // Array to store the best possible score for deleting at each position
//            int[] bestDelete = new int[l2];
//
//            // Fill array with first predecessor position
//            Arrays.fill(bestDelete, i);
//
//
//            // Get array of scores equal to previous row + the cost of moving to a match in current row
//            double[] matchScore = MatrixUtils.addArrays(Arrays.copyOfRange(scores[i], 0, scores[0].length - 1), matchPoints);
//
//            // Array to store the best possible score for matching at each position
//            int[] bestMatch = new int[l2];
//
//            // Fill array with first predecessor position
//            Arrays.fill(bestMatch, i);
//
//
////            for (int j = 1; j < predecessors.size(); j++) {
////                int predecessor = predecessors.get(j);
////                double[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);
////
////
////                for (int k = 0; k < bestDelete.length; k++) {
////                    int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
////                    bestDelete[k] = bestDeleteElement;
////
////                    int bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
////                    deleteScore[k] = bestDeleteScoreElement;
////                }
////
////
////                double[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);
////
////                for (int n = 0; n < bestMatch.length; n++) {
////                    int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
////                    bestMatch[n] = bestMatchElement;
////
////                    int bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
////                    matchScore[n] = bestMatchScoreElement;
////                }
////            }
//
//
//            // Array to indicate whether we deleted at a position
//            boolean[] deleted = new boolean[deleteScore.length];
//
//            // Update the deleted array and put the best score into scores
//            for (int k = 0; k < deleteScore.length; k++) {
//                deleted[k] = (deleteScore[k]) >= matchScore[k];
//                scores[i + 1][k + 1] = max(deleteScore[k], matchScore[k]);
//            }
//
//
//            for (int l = 1; l < scores[0].length; l++) {
//                scores[i + 1][l] = max(deleteScore[l - 1], matchScore[l - 1]);
//
//                /* If there was a deletion, update the backIndexes to reflect predecessor position it came from and
//                   update deleteCost array at this point to be an extension of a gap, not a gap opening
//                */
//                if (deleted[l - 1]) {
//                    iBackIndexes[i + 1][l] = bestDelete[l - 1];
//                    jBackIndexes[i + 1][l] = l;
//                    deleteCost[i + 1][l] = extendGapPenalty;
//
//                    // Otherwise update the backindexes to reflect predecessor position the match came from
//                } else {
//                    iBackIndexes[i + 1][l] = bestMatch[l - 1];
//                    jBackIndexes[i + 1][l] = l - 1;
//                }
//            }
//
//            // Check if there should be an insertion in the current row
//
//            // Get current scores that reflect match / deletion scores
//            double[] scoreRow = scores[i + 1];
//
//            // Get the current cost of insertion at current position
//            double[] insertRow = insertCost[i + 1];
//
//            // Get the optimal predecessor j indexes
//            int[] jBackRow = jBackIndexes[i + 1];
//            Arrays.fill(inserted, false);
//
//            // Get array of scores equal to match scores plus opening an insertion in current row
//            insscores = MatrixUtils.addArrays(scoreRow, insertRow);
//
//            for (int n = 0; n < l2; n++) {
//                /* If opening an insertion has equal or greater score than match or delete, update scores to reflect
//                   opening an insertion
//                  */
//
//                if (insscores[n] >= scoreRow[n + 1]) {
//                    scores[i + 1][n + 1] = insscores[n];
//                    scoreRow[n + 1] = insscores[n];
//
//                    // Also update the insertCost array at this point to be an extension of a gap, not a gap opening
//                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
//                    insertRow[n + 1] = this.extendGapPenalty;
//
//                    //TODO: Think this is where the lack of optional pathways is bleeding from
//                    // Update the backindexes to reflect predecessor position the match came from
//                    jBackIndexes[i + 1][n + 1] = n;
//                    jBackRow[n + 1] = n;
//                    inserted[n + 1] = true;
//
//                    // Update insertScores to be an extension of a gap, not a gap opening
//                    insscores[n + 1] = scores[i + 1][n + 1] + insertCost[i + 1][n + 1];
//                }
//            }
//
//            //TODO: Always favouring insertions over deletions!!
//            // Update the inserted and deleted costs to reflect the insertion and
//            for (int o = 0; o < inserted.length; o++) {
//                if (inserted[o]) {
//                    insertCost[i + 1][o] = this.extendGapPenalty;
//                    deleteCost[i + 1][o] = this.openGapPenalty;
//                    iBackIndexes[i + 1][o] = i + 1;
//                }
//            }
//        }
//
//        return backtrack(scores, jBackIndexes, iBackIndexes);
//
//    }
//
//    /**
//     * Update the two profiles with the gaps given by aligning them
//     * @param profile1Indexes List of integers representing the positions to add gaps in profile 1
//     * @param profile2Indexes List of integers representing the positions to add gaps in profile 2
//     */
//    public void addGapsToProfiles(List<Integer> profile1Indexes, List<Integer> profile2Indexes){
//
//        // Create list of positions with gaps in them
//        List<Integer> gapPos = new ArrayList<Integer>();
//        List<Integer> gapPos2 = new ArrayList<Integer>();
//
//        for (int i = 0; i < profile1Indexes.size(); i++) {
//            if (profile1Indexes.get(i) == -1) {
//                gapPos.add(i);
//            }
//        }
//
//        for (int i = 0; i < profile2Indexes.size(); i++) {
//            if (profile2Indexes.get(i) == -1) {
//                gapPos2.add(i);
//            }
//
//        }
//
//        // Update each profile to have gaps in them
//        if (gapPos.size() > 0) {
//            profile1.addGaps(gapPos);
//
//        }
//
//        if (gapPos2.size() > 0) {
//            profile2.addGaps(gapPos2);
//        }
//
//
//    }
//
//
//}