package bn.math;

import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;
import dat.file.AlnWriter;
import dat.substitutionmodels.Blosum62;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Class to perform operations on sets of PO Graphs: union, intersection, merge.
 *
 * Created by marnie on 18/11/16.
 */
public class GraphUtility {
    private Map<String, POGraph> graphs;      // Set of graphs to perform functions on

    /**
     * Constructor for the graph utility
     */
    public GraphUtility() {
        graphs = new HashMap<>();
    }

    /**
     * Constructor for the graph utility. Store the given set of graphs to perform operations on.
     *
     * @param graphs Set of graphs to perform the operations on
     */
    public GraphUtility(Map<String, POGraph> graphs) {
        this();
        this.graphs.putAll(graphs);
    }

    /**
     * Add graph to the set of graphs to consider.
     *
     * @param label label of graph to add
     * @param graph graph to add
     */
    public void addGraph(String label, POGraph graph) {
        if (graphs.keySet().contains(label))
            throw new RuntimeException("Label for partial order graph must be unique: " + label + " already refers to a graph.");
        graphs.put(label, graph);
    }

    /**
     * Perform the union operation on the set of graphs loaded in the class and save as the given filepath.
     *
     * @param filepath  filepath to save resulting structure
     */
    public void savePerformUnion(String filepath) {
        POGraph result = union();
        result.saveSequences(filepath + "_union", "fasta");
        result.saveToDot(filepath + "_union");
    }

    /**
     * Perform the union operation on the set of graphs loaded in the class.
     *
     * @return  Partial order graph representation of the union of all stored graphs
     */
    public POGraph union() {
        if (graphs.isEmpty())
            return null;

        EnumSeq.Gappy<Enumerable>[] seqs;

        String[] labels = new String[graphs.size()];
        graphs.keySet().toArray(labels);
        POGraph[] poGraphs = new POGraph[graphs.size()];
        graphs.values().toArray(poGraphs);

        POGraph g1 = poGraphs[0];
        String labelg1 = labels[0];
        seqs = new EnumSeq.Gappy[g1.getSequences().size()];
        for (int pind = 1; pind < poGraphs.length; pind++) {
            POGraph g2 = poGraphs[pind];
            String labelg2 = labels[pind];

            List<Integer> g1Path = new ArrayList<>();
            List<Integer> g2Path = new ArrayList<>();
            identifyPath(g1, g2, g1Path, g2Path);

            // create 'aligned' sequences to create new graph
            seqs = new EnumSeq.Gappy[g1.getSequences().size() + g2.getSequences().size()];
            Map<Integer, Object[]> mergedSequencesAln1 = constructGappySequences(g1, g1Path, g1Path.size());
            populateEnumSequences(labelg1, g1.getSequences(), mergedSequencesAln1, seqs);
            Map<Integer, Object[]> mergedSequencesAln2 = constructGappySequences(g2, g2Path, g2Path.size());
            populateEnumSequences(labelg2, g2.getSequences(), mergedSequencesAln2, seqs);

            g1 = constructGraph(seqs);
            labelg1 = null;
        }

        return g1;
    }

    /**
     * Perform the intersection operation on the set of graphs loaded in the class and save as the given filepath.
     *
     * @param filepath  filepath to save resulting structure
     */
    public void savePerformIntersection(String filepath) {
        POGraph result = intersection();
        result.saveSequences(filepath + "_intersection", "fasta");
        result.saveToDot(filepath + "_intersection");
    }

    /**
     * Perform the intersection operation on the set of graphs loaded in the class. Identifies the nodes that occur in
     * all stored graphs.
     *
     * @return Partial order graph representation of the intersection of all stored graphs
     */
    public POGraph intersection() {
        if (graphs.isEmpty())
            return null;

        EnumSeq.Gappy<Enumerable>[] seqs;

        String[] labels = new String[graphs.size()];
        graphs.keySet().toArray(labels);
        POGraph[] poGraphs = new POGraph[graphs.size()];
        graphs.values().toArray(poGraphs);

        POGraph g1 = poGraphs[0];
        String labelg1 = labels[0];
        seqs = new EnumSeq.Gappy[g1.getSequences().size()];
        for (int pind = 1; pind < poGraphs.length; pind++) {
            POGraph g2 = poGraphs[pind];
            String labelg2 = labels[pind];

            List<Integer> g1Path = new ArrayList<>();
            List<Integer> g2Path = new ArrayList<>();
            identifyPath(g1, g2, g1Path, g2Path);

            // remove nodes that pair with nulls
            List<Integer> intG1Path = new ArrayList<>();
            List<Integer> intG2Path = new ArrayList<>();
            for (int idx = 0; idx < g1Path.size(); idx++)
                if (g1Path.get(idx) != null && g2Path.get(idx) != null) {
                    intG1Path.add(g1Path.get(idx));
                    intG2Path.add(g2Path.get(idx));
                }

            // create 'aligned' sequences to create new graph
            seqs = new EnumSeq.Gappy[g1.getSequences().size() + g2.getSequences().size()];
            Map<Integer, Object[]> mergedSequencesAln1 = constructGappySequences(g1, intG1Path, intG1Path.size());
            populateEnumSequences(labelg1, g1.getSequences(), mergedSequencesAln1, seqs);
            Map<Integer, Object[]> mergedSequencesAln2 = constructGappySequences(g2, intG2Path, intG2Path.size());
            populateEnumSequences(labelg2, g2.getSequences(), mergedSequencesAln2, seqs);

            g1 = constructGraph(seqs);
            labelg1 = null;
        }

        return g1;
    }

    /**
     * Perform the intersection operation on the set of graphs loaded in the class and save as the given filepath.
     *
     * @param filepath  filepath to save resulting structure
     */
    public void savePerformComplement(String filepath) {
        POGraph result = complement();
        result.saveSequences(filepath + "_complement", "fasta");
        result.saveToDot(filepath + "_complement");
    }

    /**
     * Perform the complement of the given partial order graph and the set of stored graphs. Identifies the nodes that
     * are in the stored graphs, but not in the provided base graph.
     *
     * @return      Partial order graph representing the complement of the set of stored graphs
     */
    public POGraph complement() {
        if (graphs.isEmpty())
            return null;

        // create union and intersection graphs; complement is union - intersect
        POGraph unionGraph = union();
        POGraph intersectGraph = intersection();

        List<Integer> g1Path = new ArrayList<>();
        List<Integer> g2Path = new ArrayList<>();
        identifyPath(unionGraph, intersectGraph, g1Path, g2Path);

        // remove nodes that do not pair with a null (that is, that pair with another node, and is consequently aligned)
        List<Integer> intG1Path = new ArrayList<>();
        List<Integer> intG2Path = new ArrayList<>();
        for (int idx = 0; idx < g1Path.size(); idx++)
            if (g1Path.get(idx) == null || g2Path.get(idx) == null) {
                intG1Path.add(g1Path.get(idx));
                intG2Path.add(g2Path.get(idx));
            }

        // create 'aligned' sequences to create new graph
        EnumSeq.Gappy<Enumerable>[] seqs = new EnumSeq.Gappy[unionGraph.getSequences().size()];
        Map<Integer, Object[]> mergedSequencesAln1 = constructGappySequences(unionGraph, intG1Path, intG1Path.size());
        populateEnumSequences(null, unionGraph.getSequences(), mergedSequencesAln1, seqs);

        return constructGraph(seqs);
    }

    /**
     * Identify path that matches nodes in the graphs by generating s score matrix that represents node similarities
     * between g1 and g2. Finds path that joins graphs, returned as node ID lists for graph 1 and graph 2. Currently
     * greedy search.
     *
     * @param g1       Graph 1
     * @param g2       Graph 2
     * @param g1Pth    Empty list to populate with the path of nodes for graph g1
     * @param g2Pth    Empty list to populate with the path of nodes for graph g2
     */
    private void identifyPath(POGraph g1, POGraph g2, List<Integer> g1Pth, List<Integer> g2Pth) {

        double[][] scores = new double[g1.getNumNodes() + 1][g2.getNumNodes() + 1];

        // based on BLOSUM (collection of base characters)
        List<Integer> g1Ids = g1.getNodeIDs();
        for (int ng1 = 0; ng1 < g1Ids.size(); ng1++) {
            int nodeG1 = g1Ids.get(ng1);
            g1.setCurrent(nodeG1);
            List<Integer> g2Ids = g2.getNodeIDs();
            for (int ng2 = 0; ng2 < g2Ids.size(); ng2++) {
                int nodeG2 = g2Ids.get(ng2);
                g2.setCurrent(nodeG2);
                double baseDist = 0;
                if (g1.getCurrentBase() == null || g2.getCurrentBase() == null) {
                    for (Character bg1 : g1.getSequenceCharacterMapping().values())
                        for (Character bg2 : g2.getSequenceCharacterMapping().values())
                            baseDist += (1.0 * Blosum62.getDistance(bg1, bg2));
                    baseDist /= (g1.getSequenceCharacterMapping().size() * g2.getSequenceCharacterMapping().size());
                } else
                    baseDist = (1.0 * Blosum62.getDistance(g1.getCurrentBase(), g2.getCurrentBase()));
                scores[ng1 + 1][ng2 + 1] = baseDist;
            }
        }

        // Identify path through scores that maximizes scores, only consider adjacent nodes
        // TODO: consider gaps (i.e. nodes can be aligned with non-consecutive nodes)...
        Double[][] adjustedScores = new Double[g1.getNumNodes() + 1][g2.getNumNodes() + 1];
        for (Double[] row : adjustedScores)
            Arrays.fill(row, 0.0);
        g1.reset();
        g2.reset();
        updateAdjustedScores(g1, g2, adjustedScores, scores, 0, 0);

        // find aligned nodes as those with the greatest scores that agree (column/row)
        List<Integer> g2Ids = g2.getNodeIDs();
        Map<Integer, Integer> alignG1 = new HashMap<>();
        Map<Integer, Integer> alignG2 = new HashMap<>();
        for (int ng1 = 0; ng1 < g1Ids.size(); ng1++) {
            int nodeG1 = g1Ids.get(ng1);
            double bestScore = 0.0;
            int bestG2 = 0;
            for (int ng2 = 0; ng2 < g2Ids.size(); ng2++) {
                int nodeG2 = g2Ids.get(ng2);
                if (adjustedScores[ng1 + 1][ng2 + 1] > bestScore) {
                    bestScore = adjustedScores[ng1 + 1][ng2 + 1];
                    bestG2 = nodeG2;
                }
            }
            alignG1.put(nodeG1, bestG2);
        }
        for (int ng2 = 0; ng2 < g2Ids.size(); ng2++) {
            int nodeG2 = g2Ids.get(ng2);
            double bestScore = 0.0;
            int bestG1 = 0;
            for (int ng1 = 0; ng1 < g1Ids.size(); ng1++) {
                int nodeG1 = g1Ids.get(ng1);
                if (adjustedScores[ng1 + 1][ng2 + 1] > bestScore) {
                    bestScore = adjustedScores[ng1 + 1][ng2 + 1];
                    bestG1 = nodeG1;
                }
            }
            alignG2.put(nodeG2, bestG1);
        }

        // pad paths for alignment
        for (Integer g1N : alignG1.keySet()) {
            g1Pth.add(g1N);
            if (alignG2.get(alignG1.get(g1N)) == g1N)
                g2Pth.add(alignG1.get(g1N));
            else
                g2Pth.add(null);
        }
        int idx = 0;
        for (Integer g2N : alignG2.keySet()) {
            if (!g2Pth.contains(g2N)) {
                // find the number of gaps preceding the node
                int gapCount = 0;
                for (int g2ind = 0; g2ind < idx; g2ind++)
                    if (g2Pth.get(g2ind) == null)
                        gapCount++;
                g2Pth.add(idx + gapCount, g2N);
                g1Pth.add(idx + gapCount, null);
            }
            idx++;
        }
    }

    /**
     * Recursively populate adjusted score matrix by traversing edges in the graphs. Checkes graph 1 edges, then graph 2
     * edges, then the combination, for each node.
     *
     * @param g1                graph 1 to align
     * @param g2                graph 2 to align
     * @param adjustedScores    adjusted score matrix to populate
     * @param scores            score matrix representing pairwise base distances
     * @param curN1             current node index position in matrix for graph 1
     * @param curN2             current node index position in matrix for graph 2
     */
    private void updateAdjustedScores(POGraph g1, POGraph g2, Double[][] adjustedScores, double[][] scores, int curN1, int curN2){
        // check only next pointers in g1
        for (Integer nextN : g1.getNextIDs()) {
            int ind = g1.getNodeIDs().indexOf(nextN);
            adjustedScores[ind + 1][curN2 + 1] = Math.max(adjustedScores[ind + 1][curN2 + 1], scores[ind + 1][curN2 + 1]);
        }
        // check only next pointers in g2
        for (Integer nextN : g2.getNextIDs()) {
            int ind = g2.getNodeIDs().indexOf(nextN);
            adjustedScores[curN1 + 1][ind + 1] = Math.max(adjustedScores[curN1 + 1][ind + 1], scores[curN1 + 1][ind + 1]);
        }
        // check combined next pointers, g1 and g2
        for (Integer nextG1 : g1.getNextIDs()) {
            int ind1 = g1.getNodeIDs().indexOf(nextG1);
            for (Integer nextG2 : g2.getNextIDs()) {
                int ind2 = g2.getNodeIDs().indexOf(nextG2);
                adjustedScores[ind1 + 1][ind2 + 1] = Math.max(adjustedScores[ind1 + 1][ind2 + 1], scores[ind1 + 1][ind2 + 1]);
            }
        }
        if (g1.getNumEdgesOut() == 0 && g2.getNumEdgesOut() == 0)
            return;
        if (g1.getNumEdgesOut() == 0)
            for (Integer nextN : g2.getNextIDs()) {
                g2.setCurrent(nextN);
                int ind = g2.getNodeIDs().indexOf(nextN);
                updateAdjustedScores(g1, g2, adjustedScores, scores, curN1, ind);
            }
        else if (g2.getNumEdgesOut() == 0)
            for (Integer nextN : g1.getNextIDs()) {
                g1.setCurrent(nextN);
                int ind = g1.getNodeIDs().indexOf(nextN);
                updateAdjustedScores(g1, g2, adjustedScores, scores, ind, curN2);
            }
        else
            for (Integer nextG1 : g1.getNextIDs())
                for (Integer nextG2 : g2.getNextIDs()) {
                    g2.setCurrent(nextG2);
                    g1.setCurrent(nextG1);
                    int nodeG1 = g1.getNodeIDs().indexOf(nextG1);
                    int nodeG2 = g2.getNodeIDs().indexOf(nextG2);
                    updateAdjustedScores(g1, g2, adjustedScores, scores, nodeG1, nodeG2);
                }
    }

    /**
     * Constructs a gappy sequences of given length for all sequences in graph, based on the node indexes,
     *
     * @param graph         graph to construct sequences from
     * @param nodeIndex     list of ordered node indices
     * @param seqLen        length of alignment
     * @return              map of sequence ID and gappy sequence
     */
    private Map<Integer, Object[]> constructGappySequences(POGraph graph, List<Integer> nodeIndex, int seqLen){
        // construct gappy sequences as strings
        Map<Integer, List<Integer>> graphSeqNodeMap = graph.getSequenceNodeMapping();
        Map<Integer, Object[]> mergedSequencesAln = new HashMap<>();

        for (Integer seqId : graphSeqNodeMap.keySet()) {
            mergedSequencesAln.put(seqId, new Object[seqLen]);
            Arrays.fill(mergedSequencesAln.get(seqId), '-');
            for (Integer nodeId : graphSeqNodeMap.get(seqId))
                if (nodeIndex.contains(nodeId)) {
                    graph.setCurrent(nodeId);
                    mergedSequencesAln.get(seqId)[nodeIndex.indexOf(nodeId)] = graph.getSequenceCharacterMapping().get(seqId);
                }
        }

        return mergedSequencesAln;
    }

    /**
     * Populate EnumSeq with gappy sequence given by sequences
     *
     * @param graphLbl              Label to append to the sequence labels (so unique between graphs)
     * @param sequenceIDLabelMap    Sequence ID to label mapping
     * @param sequences             Map of sequence ID and the gappy sequences
     * @param seqs                  EnumSeq array to populate
     */
    private void populateEnumSequences(String graphLbl, Map<Integer, String> sequenceIDLabelMap, Map<Integer, Object[]> sequences, EnumSeq.Gappy<Enumerable>[] seqs) {
        int offset = 0;
        while (seqs[offset] != null)
            offset++;
        for (Integer seqId : sequenceIDLabelMap.keySet()) {
            seqs[offset + seqId] = new EnumSeq.Gappy<>(Enumerable.aacid_ext);
            seqs[offset + seqId].setInfo(Arrays.deepToString(sequences.get(seqId)));
            seqs[offset + seqId].set(sequences.get(seqId));
            seqs[offset + seqId].setName(graphLbl == null ? sequenceIDLabelMap.get(seqId) : graphLbl + '_' + sequenceIDLabelMap.get(seqId));
        }
    }

    /**
     * Construct partial order graph from aligned sequences.
     *
     * @param seqs  array of aligned sequences
     * @return      partial order graph representation of aligned sequences
     */
    private POGraph constructGraph(EnumSeq.Gappy<Enumerable>[] seqs) {
        POGraph result = null;
        try {
            // write temporary aln file
            File alnfile = new File("tmp.aln");
            AlnWriter aw = new AlnWriter(alnfile);
            aw.save(seqs);
            aw.close();
            result = new POGraph("tmp.aln", "tmp.aln"); // load PO Graph from aln file
            alnfile.delete(); // delete temporary file
        } catch (IOException e) {
            e.printStackTrace();
        }

        return result;
    }
}
