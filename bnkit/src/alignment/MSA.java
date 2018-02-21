package alignment;

import alignment.utilities.SubstitutionMatrix;
import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static java.lang.Math.max;

/**
 * Performs multiple sequence alignment using a partial order graph.
 *
 * Created by gabe, marnie on 20/10/16.
 */
public class MSA {

    private double openGapPenalty = -4;
    private double extendGapPenalty = -2;
    private POGraph graph;
    private SubstitutionMatrix subMatrix;
    private List<Integer> sortedIDs;


    public MSA(String filepath) throws IOException{

        this(filepath, -4, -2, new SubstitutionMatrix("blosum62"), false, false);
    }

    public MSA(String filepath, String alignmentType) throws IOException{
        this(filepath, 0.01, 0.39, 0.011, 0.05, 0.05, new SubstitutionMatrix("blosum62EstimatedWithX"), alignmentType, true, true);
    }

    /**
     * Constructor for the multiple sequence alignment.
     *
     * @param filepath  filepath of the sequence file (.fasta or clustal format)
     */
    public MSA(String filepath, int openGapPenalty, int extendGapPenalty, SubstitutionMatrix subMatrix, boolean partialOrder, boolean partialOrderTraceback) throws IOException {

        this.openGapPenalty = openGapPenalty;
        this.extendGapPenalty = extendGapPenalty;
        this.subMatrix = subMatrix;

        //try {

            List<EnumSeq.Gappy<Enumerable>> seqs = getSeqs(filepath);
            graph = getGraph(seqs);
            System.out.println("Using Needleman Wunsch alignment");

            for (int seqId = 1; seqId < seqs.size(); seqId++) {
                // align sequence with the graph
                List<Integer> alignment = alignSeqToGraph(seqs.get(seqId).toString(), false, partialOrder, partialOrderTraceback);
                graph.addSequence(seqId, seqs.get(seqId).getName(), seqs.get(seqId).toString(), alignment);
//                graph.saveSequences(filepath + seqId, "fasta");
            }
            saveMSA(filepath);

       // } catch (IOException e) {
        //    e.printStackTrace();
        //}
    }

    public  MSA(POGraph graph,  int openGapPenalty, int extendGapPenalty, SubstitutionMatrix subMatrix, int seqId) throws IOException {

        this.graph = graph;
        this.subMatrix = subMatrix;
        this.openGapPenalty = openGapPenalty;
        this.extendGapPenalty = extendGapPenalty;
    }

    public  void alignMEASequence(int seqId, EnumSeq<Enumerable> seq, boolean partialOrder, boolean partialOrderTraceback){
        List<Integer> alignment = alignSeqToGraph(seq.toString(), true, partialOrder, partialOrderTraceback);
        graph.addSequence(seqId, seq.getName(), seq.toString(), alignment);




    }

    private MSA(String filepath, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type, boolean partialOrder, boolean partialOrderTraceback) throws IOException {
        try {

            List<EnumSeq.Gappy<Enumerable>> seqs = getSeqs(filepath);
            POGraph graph = getGraph(seqs);
            PairHMM pairHMM = new PairHMM(graph, seqs, tau, epsilon, delta, emissionX, emissionY, subMatrix, type, partialOrder, partialOrderTraceback);


            if (type.equals("POViterbi")){
                System.out.println("Using Partial Order Viterbi alignment");
                pairHMM.getViterbiAlignment();

            }

            else if (type.equals("Viterbi")){
                System.out.println("Using Viterbi alignment");
                pairHMM.getViterbiAlignment();
            }

            else if (type.equals("POMEA")){
                System.out.println("Using Partial order maximum expected accuracy alignment");
                pairHMM.getMEAAlignment(1);
            }

            else if (type.equals("MEA")){
                System.out.println("Using Maximum expected accuracy alignment");

//                pairHMM = new PairHMM(graph, seqs, tau, epsilon, delta, emissionX, emissionY, subMatrix, type, false);

                pairHMM.getMEAAlignment(1);
            }





        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private List<EnumSeq.Gappy<Enumerable>> getSeqs(String filepath) throws IOException {
        List<EnumSeq.Gappy<Enumerable>> seqs = new ArrayList<EnumSeq.Gappy<Enumerable>>();

        try {
            // load the sequence data
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            String line = br.readLine();
            if (line.startsWith("CLUSTAL")) {
                seqs = EnumSeq.Gappy.loadClustal(filepath, Enumerable.aacid);
            } else if (line.startsWith(">")) {
                Character gap = "-".charAt(0);
                seqs = EnumSeq.Gappy.loadFasta(filepath, Enumerable.aacid, gap);
            } else {
                throw new RuntimeException("Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
            }
            br.close();

            return seqs;

        } catch (IOException e) {
            e.printStackTrace();
        }
        return seqs;
    }

    private POGraph getGraph(List<EnumSeq.Gappy<Enumerable>> seqs){
        graph = new POGraph();

        // for each sequence, align with the partial order graph and add to graph
        List<Integer> nodeIds = new ArrayList<>();
        int graphLength = seqs.get(0).toString().toCharArray().length;
        for (int id = 0; id < graphLength; id++)
            nodeIds.add(id);
        graph.addSequence(0, seqs.get(0).getName(), seqs.get(0).toString(), nodeIds);
        this.graph = graph;

        return graph;
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
    private void saveMSA(String filepath){
        graph.saveSequences(filepath, "fasta");
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
     * @param seqvec characters in the sequence to be matched against
     * @return an int array containing the match scores
     */
    private double[] getMatchScore(Map<Character, Double> bases, String seqvec, SubstitutionMatrix subMatrix) {
        double[] matches = new double[seqvec.length()];
        for (int i = 0; i < seqvec.length(); i++) {
            double matchScore = 0;
            for (Character base : bases.keySet()) {
                matchScore += subMatrix.getDistance(seqvec.charAt(i), base) * bases.get(base);
            }

            matches[i] = matchScore;
        }

        return matches;
    }

    /**
     * Get the score for making a match in Maximum Expected Accuracy alignment
     * @param i
     * @param profile
     * @return
     */
    private double[] getMEAMatchScore(int i, String profile) {
        double[] matches = new double[profile.length()];
        System.arraycopy(subMatrix.getMatrix()[i + 1], 1, matches, 0, profile.length());

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
        int pos;
        sortedIDs = this.graph.getSortedIDs();

        for (Integer id : sortedIDs) {
            nodeIDToIndex.put(id, index);
            nodeIndexToID.put(index, id);
            index++;
        }

        // Create score matrix
        double[][] scores = new double[l1+1][l2+1];
        int[][] backStrIdx = new int[l1 + 1][l2 + 1];
        int[][] backGraphIdx = new int[l1 + 1][l2 + 1];

        // Fill top row with gap penalties
        for (int i = 0; i < l2; i++)
            scores[0][i+1] = openGapPenalty + (i * extendGapPenalty);


        for (int i = 0; i < this.graph.getNodeIDs().size(); ++i) {
            this.graph.setCurrent(sortedIDs.get(i));
            List<Integer> prevIdxs = this.graph.getPreviousIDs();
            if (prevIdxs.isEmpty() || prevIdxs.get(0) == -1){
                pos = 0;
                scores[i + 1][0] = openGapPenalty;
            } else {
                pos = nodeIDToIndex.get(prevIdxs.get(0)) + 1;
            }
            scores[0][0] = openGapPenalty - extendGapPenalty;



            double best = scores[pos][0] + this.extendGapPenalty;
            for (Integer prevIdx : prevIdxs) {
                if (prevIdx == -1)
                    continue;
                double prevScore = scores[nodeIDToIndex.get(prevIdx) + 1][0] + this.extendGapPenalty;
                //TODO: In event of a tie for score... just defaults to the previous existing best
                if (best <= prevScore) {
                    backGraphIdx[i+1][0] = nodeIDToIndex.get(prevIdx) + 1;
                    scores[i+1][0] = best;


                }

                else {
                    best = prevScore;
                }
            }

            scores[0][0] = 0;
        }


        // Fill first column with optimal previous cell move
        for (int i = 0; i < l2; i++) {
            backStrIdx[0][i + 1] = i;
        }

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
    private List<Integer> backtrack(double[][] scores, int[][] backStrIdx, int[][] backGrphIdx, Map<Integer, Integer> nodeIndextoID) {
        List<Integer> orderedNodes = sortedIDs;


        int besti = scores.length;
        int bestj = scores[0].length;
        besti -= 1;
        bestj -= 1;
        List<Integer> terminalIndices = new ArrayList<>();


        // Get a list of all the indexes that are at the end of the graph
        for (Integer nodeId : this.graph.getNodeIDs()) {
            this.graph.setCurrent(nodeId);
            if (this.graph.getNextIDs().isEmpty())
                terminalIndices.add(nodeId);
        }

        double bestScore = scores[besti][bestj];

        if (terminalIndices.size() > 1) {
            System.out.println("Terminal index size greater than one");


            for (Integer terminalIndice : terminalIndices) {
                double score = scores[terminalIndice][bestj];
                if (score > bestScore) {
                    bestScore = score;
                    besti = terminalIndice;
                }
            }
        }

        List<Integer> profile1Indexes = new ArrayList<Integer>();
        List<Integer> profile2Indexes = new ArrayList<Integer>();

        int profile1Index = 0;
        int profile2Index = 0;

        while (!(besti == 0 && bestj == 0)) {
            int nexti = backGrphIdx[besti][bestj];
            int nextj = backStrIdx[besti][bestj];
            profile1Index = besti - 1;
            profile2Index = bestj - 1;

//            besti = (besti == 0) ? 1 : besti;
//            profile1Index = nodeIndextoID.get(besti - 1);

            if (nextj != bestj)
                profile2Indexes.add(0, profile2Index);
            else
                profile2Indexes.add(0, -2);
            if (nexti != besti)
                profile1Indexes.add(0, profile1Index);
            else
                profile1Indexes.add(0, -2);

            besti = nexti;
            bestj = nextj;
        }

        // Fill out the remaining indexes of each profile
        while (profile1Indexes.get(0) != null && profile1Index > 0) {
            profile1Indexes.add(0, profile1Index - 1);
            profile2Indexes.add(0, -2);
            profile1Index -= 1;

        }

        while (profile2Indexes.get(0) != null && profile2Index > 0) {
            profile2Indexes.add(0, profile2Index - 1);
            profile1Indexes.add(0, -2);
            profile2Index -= 1;

        }

        List<Integer> matchesIndex = new ArrayList<Integer>();

//        if (scores[0].length > scores.length) {
//            matchesIndex = profile1Indexes;
//
//        }
//
//
//        else {
//             matchesIndex = new ArrayList<Integer>();

        //TODO: Does this mess with normal NW alignment?
        for (int i = 0; i < profile2Indexes.size(); i++){

            if (profile2Indexes.get(i) != -2){
                int index = profile1Indexes.get(i) == -2 ? -2: orderedNodes.get(profile1Indexes.get(i));

                matchesIndex.add(index);
            }
        }



        return matchesIndex;
    }

    /**
     * Align a new sequence to the PO Graph
     *
     * @param   sequence    sequence to align and add to PO Graph
     * @return
     */
    private List<Integer> alignSeqToGraph(String sequence, boolean MEA, boolean partialOrder, boolean partialOrderTraceback) {

        int l1 = this.graph.getNumNodes();
        int l2 = sequence.length();



        List<Object> initialisedData = this.initialiseDynamicProgrammingData(sequence);
        HashMap<Integer, Integer> nodeIDToIndex = (HashMap<Integer, Integer>) initialisedData.get(0);
        HashMap<Integer, Integer> nodeIndexToID = (HashMap<Integer, Integer>) initialisedData.get(1);
        double[][] scores = (double[][]) initialisedData.get(2);
        int[][] backStrIdx = (int[][]) initialisedData.get(3);
        int[][] backGraphIdx = (int[][]) initialisedData.get(4);
        double[][] insertCost = new double[l1 + 1][l2 + 1];
        double[][] deleteCost = new double[l1 + 1][l2 + 1];

        for (double[] row : insertCost)
            Arrays.fill(row, this.openGapPenalty);

        for (double[] row : deleteCost)
            Arrays.fill(row, this.openGapPenalty);

        // Array to keep track of whether we should insert
        boolean[] inserted = new boolean[l2 + 1];
        Arrays.fill(inserted, false);

        double[] insscores = new double[l2 + 2];
//        List<Integer> sortedIDs = this.graph.getSortedIDs();



        for (int i = 0; i < l1; i++) {
            this.graph.setCurrent(sortedIDs.get(i));

            // Get character of node
            List<Character> bases = new ArrayList<>();
            Character pbase = this.graph.getCurrentBase();

            if (pbase == null) // multiple base characters to consider
                bases.addAll(this.graph.getSequenceCharacterMapping().values());


//            // consider all 'aligned' base characters (set in current node), identify best match of character
//            Double max = this.openGapPenalty;
////            pbase = bases.get(0);
//            if (bases.size() > 1) {
//                double[] potentialBases;
//                for (Character base : bases) {
//                    if (MEA) {
//                        potentialBases = getMEAMatchScore(base, sequence);
//
//
//                    } else {
//                        potentialBases = this.getMatchScore(base, sequence, subMatrix);
//                    }
//                    for (Double score : potentialBases)
//                        if (score > max) {
//                            max = score;
//                            pbase = base;
//                            break;
//                        }
//                }
//            }

            // Get predecessors of node
            List<Integer> predecessors = this.graph.getPreviousIDs();

            //Get the actual index, not the ID for the predecessors
            for (int j = 0; j < predecessors.size(); j++){
                predecessors.set(j, nodeIDToIndex.get(predecessors.get(j)));
            }

            if (predecessors.isEmpty()) {
                predecessors.add(-1);
            }

            // Get array of scores for matching current node with each position in sequence
            double[] matchPoints;

            if (MEA){
                matchPoints = getMEAMatchScore(i, sequence);
            }

            else {

                matchPoints = this.getMatchScore(this.graph.getCharacterDistribution(), sequence, subMatrix);
            }


            // Get array of scores equal to previous row + the cost of a deletion at each position
            double[] deleteScore = Arrays.copyOfRange(scores[i], 1, scores[0].length);


            deleteScore = addArrays(deleteScore, Arrays.copyOfRange(deleteCost[i ], 1, deleteCost[0].length));



            // Array to store the best possible score for deleting at each position
            int[] bestDelete = new int[l2];


            // Fill array with first predecessor position
            Arrays.fill(bestDelete, i);

            double[] matchScore = addArrays(Arrays.copyOfRange(scores[i], 0, scores[0].length - 1), matchPoints);


            int[] bestMatch = new int[l2];
            Arrays.fill(bestMatch, i);


            if (partialOrderTraceback) {
                //TODO This checks the new score against the existing score, even when they're identical
                for (int j = 0; j < predecessors.size(); j++) {
                    int predecessor = predecessors.get(j);
                    double[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);


                    for (int k = 0; k < bestDelete.length; k++) {
                        //TODO Greater than or equal here would let us select different nodes when values are equal
                        int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
                        bestDelete[k] = bestDeleteElement;

                        double bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
                        deleteScore[k] = bestDeleteScoreElement;
                    }

                    double[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);

                    for (int n = 0; n < bestMatch.length; n++) {
                        int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
                        bestMatch[n] = bestMatchElement;

                        double bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
                        matchScore[n] = bestMatchScoreElement;
                    }
                }
            }

            boolean[] deleted = new boolean[deleteScore.length];

            for (int k = 0; k < deleteScore.length; k++) {
                deleted[k] = (deleteScore[k]) > matchScore[k];
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

            double[] scoreRow = scores[i + 1];
            double[] insertRow = insertCost[i + 1];
            Arrays.fill(inserted, false);
            insscores = addArrays(scoreRow, insertRow);


            for (int n = 0; n < l2; n++) {
                if (insscores[n] >= scoreRow[n + 1]) {
                    scores[i + 1][n + 1] = insscores[n];
                    scoreRow[n + 1] = insscores[n];

                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
                    insertRow[n + 1] = this.extendGapPenalty;

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
    private double[] addArrays ( double[] firstArray, double[] secondArray){
        double[] shorterArray = (firstArray.length < secondArray.length ? firstArray: secondArray);
        double[] longerArray = (firstArray.length > secondArray.length ? firstArray: secondArray);
        double[] finalArray = new double[longerArray.length];
        for (int i = 0; i < shorterArray.length; i++)
            finalArray[i] = (firstArray[i] + secondArray[i]);
        for (int i = shorterArray.length; i < longerArray.length; i++)
            finalArray[i] = longerArray[i];
        return finalArray;
    }
}
