package bn.reconstruction;

import dat.*;
import dat.substitutionmodels.Blosum62;

import java.io.*;
import java.util.*;

import static java.lang.Math.max;

/**
 * Performs multiple sequence alignment using a partial order graph.
 *
 * Created by gabe, marnie on 20/10/16.
 */
public class MSA {

    private int openGapPenalty = -10;
    private int extendGapPenalty = -1;
    private POGraph graph;

    /**
     * Constructor for the multiple sequence alignment.
     *
     * @param filepath  filepath of the sequence file (.fasta or clustal format)
     */
    public MSA(String filepath) {
        graph = new POGraph();

        try {
            // load the sequence data
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            List<EnumSeq.Gappy<Enumerable>> seqs;
            String line = br.readLine();
            if (line.startsWith("CLUSTAL")) {
                seqs = EnumSeq.Gappy.loadClustal(filepath, Enumerable.aacid);
            } else if (line.startsWith(">")) {
                Character gap = "-".charAt(0);
                seqs = EnumSeq.Gappy.loadFasta(filepath, Enumerable.aacid, gap);
            } else {
                throw new RuntimeException("Alignment should be in Clustal or Fasta format");
            }
            br.close();

            // for each sequence, align with the partial order graph and add to graph
            List<Integer> nodeIds = new ArrayList<>();
            for (int id = 0; id < seqs.get(0).toString().toCharArray().length; id++)
                nodeIds.add(id);
            graph.addSequence(0, seqs.get(0).getName(), seqs.get(0).toString(), nodeIds);
            for (int seqId = 1; seqId < seqs.size(); seqId++) {
                // align sequence with the graph
                List<Integer> alignment = alignSeqToGraph(seqs.get(seqId).toString());
                graph.addSequence(seqId, seqs.get(seqId).getName(), seqs.get(seqId).toString(), alignment);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Save the multiple sequence alignment to a dot file.
     *
     * @param filepath  file path to save graph
     */
    public void saveMSADot(String filepath){
        graph.saveToDot(filepath);
    }

    /**
     * Save the multiple sequence alignment as an .aln file.
     *
     * @param filepath file path to save alignment
     */
    public void saveMSA(String filepath){
        graph.saveSequences(filepath, "clustal");
    }

    /**
     * Get the partial order alignment graph representing the multiple sequence alignment.
     *
     * @return  partial order alignment graph
     */
    public POGraph getMSAGraph() {
        return graph;
    }

    /**
     * Return the scores associated with matching a single character in a PO Graph to an
     * array of characters.
     * Used to get the scores for matching each row in the score matrix to each column
     * @param base the character in the PO Graph to be matched against
     * @param seqvec characters in the sequence to be matched against
     * @return an int array containing the match scores
     */
    private int[] getMatchScore(char base, String seqvec) {
        int[] matches = new int[seqvec.length()];
        for (int i = 0; i < seqvec.length(); i++) {
            int matchScore = Blosum62.getDistance(seqvec.charAt(i), base);
            matches[i] = matchScore;
        }

        return matches;
    }

    /**
     * Set up the initial sizes and values for the scores matrix, the matrices to
     * record the optimal moves through the score matrix, and the two mappings from
     * node ID to index and index to node ID
     * @param   sequence    sequence to set up data for
     * @return              Object containing scores matrix, backSeq matrix, backGraph matrix, nodeID to
     * index, and index to node ID
     */
    private List<Object> initialiseDynamicProgrammingData(String sequence) {
        Integer l1 = this.graph.getNumNodes();
        Integer l2 = sequence.length();
        Map<Integer, Integer> nodeIDToIndex = new HashMap<>();
        Map<Integer, Integer> nodeIndexToID = new HashMap<>();

        int index = 0;
        for (Integer id : this.graph.getNodeIDs()) {
            nodeIDToIndex.put(id, index);
            nodeIndexToID.put(index, id);
            index++;
        }

        // Create score matrix
        int[][] scores = new int[l1+1][l2+1];

        // Fill top row with gap penalties
        for (int i = 0; i < l2; i++)
            scores[0][i+1] = openGapPenalty + (i * extendGapPenalty);

        for (int i = 0; i < this.graph.getNodeIDs().size(); ++i) {
            List<Integer> prevIdxs = this.graph.getPrevIDs();
            scores[0][0] = openGapPenalty - extendGapPenalty;
            int best = scores[(prevIdxs.get(0)) + 1][0];
            for (Integer prevIdx : prevIdxs) {
                best = max(best, scores[prevIdx + 1][0]);
            }

            scores[i+1][0] = best + this.extendGapPenalty;
            scores[0][0] = 0;
        }

        int[][] backStrIdx = new int[l1 + 1][l2 + 1];
        int[][] backGraphIdx = new int[l1 + 1][l2 + 1];

        List<Object> initialisedData = new ArrayList<>();
        initialisedData.add(nodeIDToIndex);
        initialisedData.add(nodeIndexToID);
        initialisedData.add(scores);
        initialisedData.add(backStrIdx);
        initialisedData.add(backGraphIdx);

        return initialisedData;

    }

    /**
     * Trace back through the scores matrix, building the optimal alignment
     *
     * @param scores
     * @param backStrIdx
     * @param backGrphIdx
     * @param nodeIndextoID
     * @return
     */
    private List<Integer> backtrack(int[][] scores, int[][] backStrIdx, int[][] backGrphIdx, Map<Integer, Integer> nodeIndextoID) {
        int besti = scores.length;
        int bestj = scores[0].length;
        besti -= 1;
        bestj -= 1;
        List<Integer> terminalIndices = new ArrayList<>();
        for (Integer nodeId : this.graph.getNodeIDs()) {
            this.graph.setCurrent(nodeId);
            if (this.graph.getNextIDs().isEmpty())
                terminalIndices.add(nodeId);
        }

        int bestScore = scores[besti][bestj];

        for (int i = 0; i < terminalIndices.size(); i++) {
            int score = scores[terminalIndices.get(i)][bestj];
            if (score > bestScore) {
                bestScore = score;
                besti = terminalIndices.get(i);
            }
        }

        List<Integer> matches = new ArrayList<>();
        List<Integer> strIndexes = new ArrayList<>();

        while (!(besti == 0 && bestj == 0)) {
            int nexti = backGrphIdx[besti][bestj];
            int nextj = backStrIdx[besti][bestj];
            int curstrIdx = bestj - 1;

            besti = (besti == 0) ? 1 : besti;
            int curnodeIdx = nodeIndextoID.get(besti - 1);

            if (nextj != bestj)
                strIndexes.add(0, curstrIdx);
            else
                strIndexes.add(0, null);
            if (nexti != besti)
                matches.add(0, curnodeIdx);
            else
                matches.add(0, null);

            besti = nexti;
            bestj = nextj;
        }

        List<Integer> matchesIndex = new ArrayList<>();
        for (int ind = 0; ind < matches.size(); ind++)
            if (matches.get(ind) != null)
                matchesIndex.add(ind);

        return matchesIndex;
    }

    /**
     * Align a new sequence to the PO Graph
     *
     * @param   sequence    sequence to align and add to PO Graph
     * @return
     */
    private List<Integer> alignSeqToGraph(String sequence) {

        int l1 = this.graph.getNumNodes();
        int l2 = sequence.length();

        List<Object> initialisedData = this.initialiseDynamicProgrammingData(sequence);
        HashMap<Integer, Integer> nodeIDToIndex = (HashMap<Integer, Integer>) initialisedData.get(0);
        HashMap<Integer, Integer> nodeIndexToID = (HashMap<Integer, Integer>) initialisedData.get(1);
        int[][] scores = (int[][]) initialisedData.get(2);
        int[][] backStrIdx = (int[][]) initialisedData.get(3);
        int[][] backGraphIdx = (int[][]) initialisedData.get(4);
        int[][] insertCost = new int[l1 + 1][l2 + 1];
        int[][] deleteCost = new int[l1 + 1][l2 + 1];

        for (int[] row : insertCost)
            Arrays.fill(row, this.openGapPenalty);

        for (int[] row : deleteCost)
            Arrays.fill(row, this.openGapPenalty);

        // Array to keep track of whether we should insert
        boolean[] inserted = new boolean[l2 + 1];
        Arrays.fill(inserted, false);

        int[] insscores = new int[l2 + 2];

        for (int i = 0; i < this.graph.getNumNodes(); i++) {
            this.graph.setCurrent(this.graph.getNodeIDs().get(i));

            // Get character of node
            List<Character> bases = new ArrayList<>();
            Character pbase = this.graph.getCurrentBase();
            if (pbase == null) // multiple base characters to consider
                bases.addAll(this.graph.getCurrentBases());
            else
                bases.add(pbase);

            // consider all 'aligned' base characters (set in current node), identify best match of character
            int max = this.openGapPenalty;
            pbase = bases.get(0);
            if (bases.size() > 1)
                for (Character base : bases)
                    for (Integer score : this.getMatchScore(base, sequence))
                        if (score > max) {
                            max = score;
                            pbase = base;
                            break;
                        }

            // Get predecessors of node
            List<Integer> predecessors = this.graph.getPrevIDs();
            if (predecessors.isEmpty())
                predecessors.add(0);

            // Get array of scores for matching current node with each position in sequence
            int[] matchPoints = this.getMatchScore(pbase, sequence);

            // Get array of scores equal to previous row + the cost of a deletion at each position
            int[] deleteScore = Arrays.copyOfRange(scores[predecessors.get(0) + 1], 1, scores[0].length);
            deleteScore = addArrays(deleteScore, Arrays.copyOfRange(deleteCost[predecessors.get(0) + 1], 1, deleteCost[0].length));


            // Array to store the best possible score for deleting at each position
            int[] bestDelete = new int[l2];

            // Fill array with first predecessor position
            Arrays.fill(bestDelete, predecessors.get(0) + 1);

            int[] matchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(0) + 1], 0, scores[0].length - 1), matchPoints);
            int[] bestMatch = new int[l2];
            Arrays.fill(bestMatch, predecessors.get(0) + 1);

            for (int j = 1; j < predecessors.size(); j++) {
                int predecessor = predecessors.get(j);
                int[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);


                for (int k = 0; k < bestDelete.length; k++) {
                    int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
                    bestDelete[k] = bestDeleteElement;

                    int bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
                    deleteScore[k] = bestDeleteScoreElement;
                }

                int[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);

                for (int n = 0; n < bestMatch.length; n++) {
                    int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
                    bestMatch[n] = bestMatchElement;

                    int bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
                    matchScore[n] = bestMatchScoreElement;
                }
            }

            boolean[] deleted = new boolean[deleteScore.length];

            for (int k = 0; k < deleteScore.length; k++) {
                deleted[k] = (deleteScore[k]) >= matchScore[k];
                scores[i + 1][k + 1] = max(deleteScore[k], matchScore[k]);
            }

            for (int l = 1; l < scores[0].length; l++) {
                scores[i + 1][l] = max(deleteScore[l - 1], matchScore[l - 1]);

                if (deleted[l - 1]) {
                    backGraphIdx[i + 1][l] = bestDelete[l - 1];
                    backStrIdx[i + 1][l] = l;
                    deleteCost[i + 1][l] = extendGapPenalty;
                } else {
                    backGraphIdx[i + 1][l] = bestMatch[l - 1];
                    backStrIdx[i + 1][l] = l - 1;
                }
            }

            int[] scoreRow = scores[i + 1];
            int[] insertRow = insertCost[i + 1];
            Arrays.fill(inserted, false);
            insscores = addArrays(scoreRow, insertRow);

            for (int n = 0; n < l2; n++) {
                if (insscores[n] >= scoreRow[n + 1]) {
                    scores[i + 1][n + 1] = insscores[n];
                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
                    backStrIdx[i + 1][n + 1] = n;
                    inserted[n + 1] = true;
                    insscores[n + 1] = scores[i + 1][n + 1] + insertCost[i + 1][n + 1];
                }
            }

            for (int o = 0; o < inserted.length; o++) {
                if (inserted[o]) {
                    insertCost[i + 1][o] = this.extendGapPenalty;
                    deleteCost[i + 1][o] = this.openGapPenalty;
                    backGraphIdx[i + 1][o] = i + 1;
                }
            }
        }

        return backtrack(scores, backStrIdx, backGraphIdx, nodeIndexToID);
    }

    /**
     * Method to add two arrays together by adding the contents of each respective index in each array together.
     * For example
     * [1,4,9,6] + [2,2,3] = [3,7,12,6]
     * @param firstArray first array to add
     * @param secondArray second array to add
     * @return int array with the add
     */
    private int[] addArrays ( int[] firstArray, int[] secondArray){
        int[] shorterArray = (firstArray.length < secondArray.length ? firstArray: secondArray);
        int[] longerArray = (firstArray.length > secondArray.length ? firstArray: secondArray);
        int[] finalArray = new int[longerArray.length];
        for (int i = 0; i < shorterArray.length; i++)
            finalArray[i] = (firstArray[i] + secondArray[i]);
        for (int i = shorterArray.length; i < longerArray.length; i++)
            finalArray[i] = longerArray[i];
        return finalArray;
    }
}
