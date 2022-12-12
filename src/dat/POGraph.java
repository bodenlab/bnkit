package dat;

import alignment.utilities.MutableInt;
import bn.Distrib;
import bn.prob.EnumDistrib;
import dat.file.AlnWriter;
import dat.file.DotWriter;
import dat.file.FastaWriter;
import reconstruction.Inference;
import reconstruction.POGEdgeMap;

import java.io.*;
import java.util.*;

/**
 * Implementation of a partial order graph. POGraph stores:
 * 		- list of Nodes (private class implemented below) that contain information:
 * 				base character and a list of Sequence structures that store the ID and base character of sequences that
 * 				pass through the node and pointers to next and previous nodes.
 * 		- current node pointer
 * 		- mapping between sequences and an ordered list of Nodes that the sequence traverses
 * 		- list of starting nodes
 * 		- mapping between sequence ID and sequence label
 *
 * @author Gabe, Marnie, some amendments by Mikael
 *
 */
public class POGraph {
	private Map<Integer, Node> nodes;					// map of node ID and nodes in the structure (for fast node retrieval)
	private Node current;								// pointer to the current node
	private Node initialNode = null;					// initial node (null node), points to the first actual nodes as 'next' nodes
	private Node finalNode = null; 						// final node (null node), 'previous' pointers point to the final actual nodes
	private Map<Integer, String> sequences;				// map of sequence ID and sequence label
	private Map<Integer, Map<Integer, Map<Integer, Integer>>>  edgeCountsMSA; 	// edge counts for MSA POG
	private Map<Integer, Map<Integer, Integer>>  edgeCountsNode; 	// edge counts for node
	private int numSeqsUnderNode; 								// number of sequences under this node

	/**
	 * Constructor to initialise empty graph.
	 */
	public POGraph() {
		sequences = new HashMap<>();
		initialNode = new Node(-1);
		finalNode = new Node();
		nodes = new HashMap<>();
		current = null;
		edgeCountsMSA = new HashMap<>();
		edgeCountsNode = new HashMap<>();
	}

	/**
	 * Constructor to initialise partial order graph using a dot structure and loading the associated sequences from
	 * a fasta file.
	 *
	 * @param structure		Dot representation of partial order graph or filepath to representation, if given aligned
	 *                      sequences (.aln, .fa, .fasta), constructs a graph from the aligned sequences
	 * @param seqPath		File path to sequences
	 */
	public POGraph(String structure, String seqPath) throws IOException {
		this();
		// load sequences
		loadSequencesWithGraphStructure(structure, seqPath);
	}

	public POGraph(List<EnumSeq.Gappy<Enumerable>> seqs) {
		this();

		for (Integer seqId = 0; seqId < seqs.size(); seqId++)
			sequences.put(seqId, seqs.get(seqId).getName());

		current = loadPOGraph(seqs);
	}

	/**
	 * Constructor to initialise a partial order graph using a dot structure.
	 *
	 * @param structure		Dot representation of partial order graph
	 */
	public POGraph(String structure) throws IOException {
		this();
		if (structure.endsWith(".dot"))
			current = loadPOGraph(structure);
		else
			loadSequencesWithGraphStructure(structure, structure);
	}

	/**
	 * Copy constructor for the partial order graph
	 *
	 * @param copy	POGraph to copy from
	 */
	public POGraph(POGraph copy) {
		this();
		this.sequences.putAll(copy.sequences);
		this.initialNode = copy.initialNode.copy();
		for (Node node : this.initialNode.getNextNodes())
			addNode(node);
		this.finalNode = this.nodes.get(this.nodes.size()-1);
		this.nodes.remove(this.nodes.size()-1);
		this.current = initialNode.getNextNodes().get(0);
	}


	/**
	 * Constructor to initialise partial order graph from inferred characters and jumps in POG.
	 * FIXME: make a new sub-class of POGraph for AncestorPOGraph with this as the only constructor?
	 */
	public POGraph(List<Inference> inferred, Map<Integer, String> sequences, Map<Integer, Map<Integer, Integer>> edgeCountsForNode, int numSeqsUnderNode) {
		this.sequences = sequences;
		this.edgeCountsNode = edgeCountsForNode;
		this.numSeqsUnderNode = numSeqsUnderNode;
		initialNode = new Node(-1);
		finalNode = new Node(inferred.size() - 2);
		nodes = new HashMap<>();
		current = null;

		POGEdgeMap pem = new POGEdgeMap();
		for (Inference inf : inferred) {
			int idx1 = inf.pogId;
			Node n1 = null;
			if (idx1 == initialNode.ID)
				n1 = initialNode;
			else if (idx1 == finalNode.ID)
				n1 = finalNode;
			else {
				if (inf.base != '-') {
					n1 = new Node(idx1);
					n1.setBase(inf.base);
					nodes.put(idx1, n1);
				}
			}
			if (inf.base != '-' || n1 == initialNode || n1 == finalNode) {
				for (int idx2 : inf.transitions) {
					if (idx1 < idx2)
						pem.add(idx1, idx2);
					else if (idx1 > idx2)
						pem.add(idx2, idx1);
				}
			}
		}
		Set<POGEdgeMap.POGEdge> deleted = new HashSet<>();
		for (POGEdgeMap.POGEdge pe : pem.getEdges()) {
			int[] idxs = pe.getIndices();
			this.setCurrent(idxs[0]);
			Node n1 = current;
			if (n1 == null) { // index is referring to a deleted position
				deleted.add(pe);
				continue;
			}
			Node n2 = null;
			if (idxs[1] == finalNode.ID)
				n2 = finalNode;
			else {
				n2 = nodes.get(idxs[1]);
				if (n2 == null) { // index is referring to a deleted position
					deleted.add(pe);
					continue;
				}
			}
			// FIXME: look into populating edges with the edge counts
			Edge e1 = new Edge(n2);
			e1.addSequence(0); // temp fix so POGraph.toString can be used
			n1.nextTransitions.add(e1);
			Edge e2 = new Edge(n1);
			e2.addSequence(0); // temp fix so POGraph.toString can be used
			n2.prevTransitions.add(e2);
			boolean recip = pem.isReciprocated(pe);
			e1.setReciprocated(recip);
			e2.setReciprocated(recip);
		}
		for (POGEdgeMap.POGEdge pe : deleted)
			pem.remove(pe);
		// after the removal of POGEdge:s above, the map should be correct
	}

	/**
	 * Allows grasp to access the POG edge counts.
	 * @return
	 */
	public Map<Integer, Map<Integer, Map<Integer, Integer>>> getEdgeCounts() {
		return edgeCountsMSA;
	}

	public Map<Integer, Map<Integer, Integer>> getEdgeCountsNode() {
		return edgeCountsNode;
	}

	public int getNumSeqsUnderNode() {
		return numSeqsUnderNode;
	}

	/**
	 * Store the cost of a particular transition.
	 * @param transitionCost
	 */
	public void setTransitionCost(Map<Integer, Integer> transitionCost) {
		current.transitionCost = transitionCost;
	}

	/**
	 * Ass the cost of a particular transition.
	 * @param edge
	 * @param cost
	 */
	public void addTransitionCost(Integer edge, Integer cost) {
		current.transitionCost.put(edge, cost);
	}

	/**
	 * Remove the cost of a particular transition.
	 * @param edge
	 */
	public void removeTransitionCost(Integer edge) {
		current.transitionCost.remove(edge);
	}

	/**
	 * Gets the cost of a particular transition.
	 * @param transitionId
	 */
	public Integer getTransitionCost(Integer transitionId) {
		return current.transitionCost.get(transitionId);
	}

	/**
	 * Add sequence to the graph with characters in the order specified in nodeIds.
	 *
	 * @param id		Sequence ID
	 * @param label		Sequence label
	 * @param sequence	Sequence
	 * @param nodeIds	List of node ids to add each base character of the sequence
	 */
	public void addSequence(int id, String label, String sequence, List<Integer> nodeIds) {
		ArrayList<Integer> nodeIdsAdded = new ArrayList<>();
		char[] bases = sequence.toCharArray();
		sequences.put(id, label);
		for (int baseInd = 0; baseInd < bases.length; baseInd++) {
			Integer curId = nodeIds.get(baseInd);
			Integer prevId = initialNode.getID();
			if (nodeIds.indexOf(curId) >= 0)
				prevId = (nodeIdsAdded.size() == 0) ? initialNode.getID() : nodeIdsAdded.get(nodeIdsAdded.size()-1);
			if (curId == -2 || !setCurrent(curId)) {
				// base did not align with graph or node doesn't currently exist in the graph, create it
				current = new Node(curId == -2 ? nodes.size() : curId);
				nodes.put(current.getID(), current);
			}
			if (prevId == initialNode.getID()){
				initialNode.addNextNode(current, id);
				current.addPrevNode(initialNode, id);
			} else {
				nodes.get(prevId).addNextNode(current, id);
				current.addPrevNode(nodes.get(prevId), id);
			}
			if (baseInd == bases.length - 1) {
				current.addNextNode(finalNode, id);
				finalNode.addPrevNode(current, id);
			}
			current.addSequence(id, bases[baseInd]);
			nodeIdsAdded.add(current.getID());
		}
	}

	/**
	 * Set the pointer of the current node reference to the node with the specified ID.
	 *
	 * @param nodeID	Node ID to set as current
	 * @return			indication of whether setting the node succeeded or not
	 */
	public boolean setCurrent(Integer nodeID){
		if (nodeID.equals(finalNode.getID()))
			current = finalNode;
		else if (nodeID.equals(initialNode.getID()))
			current = initialNode;
		else
			current = nodes.get(nodeID);
		return (current != null);
	}

	public Node getNode(int index) {
		return nodes.get(index);
	}

	/**
	 * Reset the node pointer to the initial (null) node.
	 *
	 * @return	indication of whether setting the node succeeded or not
	 */
	public boolean reset(){
		current = initialNode;
		return (current != null);
	}

	/**
	 * Get the next out edge weights of the current node. Weights based on total sequence support.
	 *
	 * @return	mapping between next node ID and the edge weight to the next node
	 */
	public Map<Integer, Double> getNextEdgeWeights(){
		Map<Integer, Double> edgeWeights = new HashMap<>();
		for (Edge next : current.getNextTransitions())
			//if (next.getNext() != finalNode)
			edgeWeights.put(next.getNext().getID(), 1.0 * next.getSequences().size() / sequences.size());
		return edgeWeights;
	}

	/**
	 * Get the previous out edge weights of the current node. Weights based on total sequence support.
	 *
	 * @return	mapping between next node ID and the edge weight to the next node
	 */
	public Map<Integer, Double> getPreviousEdgeWeights(){
		Map<Integer, Double> edgeWeights = new HashMap<>();
		for (Edge next : current.getPreviousTransitions())
			//if (next.getNext() != finalNode)
			edgeWeights.put(next.getNext().getID(), 1.0 * next.getSequences().size() / sequences.size());
		return edgeWeights;
	}

	/**
	 * Get the mapping of sequence ID and sequence label.
	 *
	 * @return	map of sequence ID and label
	 */
	public Map<Integer, String> getSequences(){ return sequences; }

	/**
	 * Get the mapping of sequence ID and the list of ordered Node IDs that it traverses.
	 *
	 * @return	map of sequence ID and ordered node IDs
	 */
	public Map<Integer, List<Integer>> getSequenceNodeMapping() {
		Map<Integer, List<Integer>> seqNodeMapping = new HashMap<>();

		for (Integer seqId : sequences.keySet())
			for (Edge out : initialNode.getNextTransitions())
				if (out.getSequences().contains(seqId)) {
					ArrayList<Integer> ids = new ArrayList<>();
					ids.add(out.getNext().getID());
					seqNodeMapping.put(seqId, getNodeIdSequence(seqId, ids, out.getNext()));
					break;
				}
		return seqNodeMapping;
	}

	/**
	 * Recursively identify the sequence of node IDs for sequence with the given ID.
	 *
	 * @param seqId		ID of the sequence
	 * @param nodeIds	List of node IDs to populate recursively
	 * @param node		node to check
	 * @return			Ordered list of node IDs
	 */
	private List<Integer> getNodeIdSequence(int seqId, List<Integer> nodeIds, Node node) {
		if (node.getNextTransitions().isEmpty())
			return nodeIds;
		Node next = null;
		for (Edge out : node.getNextTransitions())
			if (out.getSequences().contains(seqId)) {
				next = out.getNext();
				break;
			}
		if (next == null)
			return nodeIds;
		nodeIds.add(next.getID());
		return getNodeIdSequence(seqId, nodeIds, next);
	}

	/**
	 * Get the mapping between sequence ID and base character for the current node.
	 *
	 * @return	mapping <seqID, base character>
	 */
	public Map<Integer, Character> getSequenceCharacterMapping(){ return (current == null) ? null : current.getSeqCharMapping(); }

	/**
	 * Set the inferred base character of the current node.
	 *
	 * @param base	inferred base character
	 */
	public void setBase(char base){
		if (current != null)
			current.setBase(base);
	}

	/**
	 * Get the unique possible base characters in the current node.
	 *
	 * @return	list of base characters
	 */
	public List<Character> getCurrentBases() {
		if (current == null)
			return null;
		List<Character> bases = new ArrayList<>();
		for (Character base : current.getSeqCharMapping().values())
			if (!bases.contains(base))
				bases.add(base);
		return bases;
	}

	/**
	 * Get a list of base characters in the current node (note: can have repeats).
	 *
	 * @return	List of base characters in current node.
	 */
	public List<Character> getCurrentBasesNotUnique() {
		return new ArrayList<>(current.getSeqCharMapping().values());
	}

	/**
	 * Get the number of nodes in the graph.
	 *
	 * @return	number of nodes in the structure.
	 */
	public int getNumNodes(){ return nodes.size(); }

	/**
	 * Get a sorted list of IDs in the graph
	 *
	 * @return topologically sorted list of IDs
	 */
	public List<Integer> getSortedIDs(){
		List<Integer> sortedIDs = topologicalSort();
		return sortedIDs;

	}

	/**
	 * Get the ID of the current node.
	 *
	 * @return	id of the current node, -1 if no node is set
	 */
	public Integer getCurrentId(){ return current.getID(); }

	/**
	 * Get the base character of the current node. Returns null if it is not set.
	 *
	 * @return	base character or null
	 */
	public Character getCurrentBase() {
		if (current == null || current.getBase() == null)
			return null;
		return current.getBase();
	}

	/**
	 * Get the node IDs of previous nodes from the current node.
	 *
	 * @return	previous node IDs
	 */
	public List<Integer> getPreviousIDs(){
		ArrayList<Integer> prevIDs = new ArrayList<>();
		if (current == null)
			return prevIDs;
		for (Node node : current.getPreviousNodes())
			prevIDs.add(node.getID());
		return prevIDs;
	}

	/**
	 * Get the node IDs of next nodes from the current node.
	 *
	 * @return	next node IDs
	 */
	public List<Integer> getNextIDs(){
		if (current == null)
			return null;
		ArrayList<Integer> nextIDs = new ArrayList<>();
		for (Node node : current.getNextNodes())
			nextIDs.add(node.getID());
		return nextIDs;
	}

	/**
	 * Sets the edge between the current node and the node with the provided nextId as reciprocated.
	 *
	 * @param nextId	ID of the next node that the reciprocated edge points to
	 */
	public void setReciprocated(Integer nextId) {
		for (Edge e : current.getNextTransitions())
			if (e.getNext().getID().equals(nextId)) {
				e.setReciprocated(true);
				for (Edge p : e.getNext().getPreviousTransitions())
					if (p.getNext().getID().equals(current.getID()))
						p.setReciprocated(true);
			}
	}

	/**
	 * Get the set of next node IDs that have reciprocated edges from the current node.
	 *
	 * @return	List of next node IDs with reciprocated edges
	 */
	public List<Integer> getReciprocatedNextIDs() {
		ArrayList<Integer> next = new ArrayList<>();
		for (Edge e : current.getNextTransitions())
			if (e.getReciprocated())
				next.add(e.getNext().getID());
		return next;
	}

	/**
	 * Remove transition to next node with the provided ID
	 *
	 * @param removeId	ID of the 'next' node transition to remove
	 */
	public void removeNextTransition(Integer removeId) {
		for (Edge edge : current.getNextTransitions())
			if (edge.getNext().getID().equals(removeId)) {
				current.removeNextNode(edge.getNext());
				if (edge.getNext().getPreviousNodes().isEmpty()) {
					Node tmp = current;
					setCurrent(edge.getNext().getID());
					removeNode();
					current = tmp;
				}
				return;
			}
	}

	public Integer getFinalNodeID(){
		return finalNode.getID();
	}

	public Integer getInitialNodeID() { return initialNode.getID(); }

	/**
	 * Remove transition to the previous node with the provided ID
	 *
	 * @param removeId	ID of the 'previous' node transition to remove
	 */
	public void removePreviousTransition(Integer removeId) {
		for (Edge edge : current.getPreviousTransitions())
			if (edge.getNext().getID().equals(removeId)) {
				current.removePrevNode(edge.getNext());
				if (current.getPreviousNodes().isEmpty()) {
					removeNode();
					current = initialNode;
				}
				return;
			}
	}

	/**
	 * Removes the current node from the graph.
	 */
	public void removeNode() {
		if (current == null)
			return;
		// set all previous nodes to point to the next node, for the sequences in the edges
		List<Edge> previousEdges = new ArrayList<>(current.getPreviousTransitions());
		for (Edge prev : previousEdges) {
			// find the next node for each sequence in the previous edge and update pointers
			List<Integer> seqs = new ArrayList<>(prev.getSequences());
			for (Integer seqId : seqs) {
				// find the next edge/node
				Edge nextEdge = null;
				for (Edge next : current.getNextTransitions())
					if (next.getSequences().contains(seqId)) {
						nextEdge = next;
						break;
					}
				// if there is an edge, join nodes otherwise just delete this sequence
				prev.removeSequence(seqId);
				// re-order transitions so that the transitions stay in order of max number of sequence traversal
				if (nextEdge != null) {
					// update seqId to go from the previous node to the next
					prev.getNext().addNextNode(nextEdge.getNext(), null);
					nextEdge.getNext().addPrevNode(prev.getNext(), null);
				}
			}
			prev.getNext().removeNextNode(current);
		}
		// remove current node from the 'previous' edges of next nodes
		List<Edge> nextEdges = new ArrayList<>(current.getNextTransitions());
		for (Edge next : nextEdges) {
			next.getNext().removePrevNode(current);
			if (next.getNext().getPreviousNodes().isEmpty()) {
				Node tmp = current;
				setCurrent(next.getNext().getID());
				removeNode();
				current = tmp;
			}
		}

		nodes.remove(current.getID());
		current = initialNode.getNextNodes().get(0);

	}

	/**
	 * Get a list of IDs of the nodes in the graph.
	 *
	 * @return	List of structure node IDs.
	 */
	public List<Integer> getNodeIDs(){
		return new ArrayList<>(nodes.keySet());
	}

	/**
	 * Save sequences represented in the partial order graph in fasta or clustal format.
	 *
	 * @param filepath	file path to save sequences
	 * @param format	file format: "fasta" or "clustal"
	 */
	public void saveSequences(String filepath, String format) {
		try {
			// determine the order of the nodes in the graph
			List<Integer> orderedNodeIds = topologicalSort();
			EnumSeq[] seqs = new EnumSeq[sequences.size()];

			for (Integer seqId : sequences.keySet()) {
				List<Character> characters = new ArrayList<>();
				String sequence  = getGappySequence(seqId, initialNode, "", orderedNodeIds);
				for (Character c : sequence.toCharArray())
					characters.add(c);
				EnumSeq seq = new EnumSeq(Enumerable.aacid_ext);
				seq.setName(sequences.get(seqId));
				seq.setInfo(characters.toString());
				Object[] chars = new Object[characters.size()];
				characters.toArray(chars);
				seq.set(chars);
				seqs[seqId] = seq;
			}
			if (format.equalsIgnoreCase("fasta")) {
				// fasta writer
				FastaWriter writer = new FastaWriter(filepath + ".fa");
				writer.save(seqs);
				writer.close();
			} else if (format.equalsIgnoreCase("clustal")) {
				// aln writer
				AlnWriter writer = new AlnWriter(filepath + ".aln");
				writer.save(seqs);
				writer.close();
			} else
				System.err.print("Incorrect file type. Must be FASTA or Clustal.");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Recursively get the character sequence for a sequence with the specified ID
	 *
	 * @param seqId		ID of sequence
	 * @param node		Initial node
	 * @param sequence	Sequence string to append to
	 * @return			Sequence characters
	 */
	private String getSequence(Integer seqId, Node node, String sequence) {
		if (node.getNextTransitions().isEmpty())
			return sequence;

		Node next = null;
		for (Edge edge : node.getNextTransitions())
			if (edge.getSequences().contains(seqId)) {
				next = edge.getNext();
				break;
			}

		if (next == null)
			return sequence;

		return getSequence(seqId, next, sequence + next.getBase());
	}

	/**
	 * Recursively get the gappy character sequence for a sequence with the specified ID.
	 *
	 * @param seqId				ID of the sequence
	 * @param node				Initial node
	 * @param sequence			Sequence string to append to
	 * @param orderedNodeIds	List of ordered Node IDs
	 * @return					Sequence characters (gappy)
	 */
	private String getGappySequence(Integer seqId, Node node, String sequence, List<Integer> orderedNodeIds) {

		Node next = null;
		String characters = (node.getSeqCharMapping().get(seqId) == null ? "" : node.getSeqCharMapping().get(seqId)) + "";

		// find next node
		for (Edge edge : node.getNextTransitions())
			if (edge.getSequences().contains(seqId)) {
				next = edge.getNext();
				break;
			}

		if (next.getNextNodes().isEmpty()) { // end of sequence
			// check if there are gaps at the end..
			int numGaps = orderedNodeIds.size() - orderedNodeIds.indexOf(node.getID()) - 1;
			characters = node.getSeqCharMapping().get(seqId) + "";
			for (int i = 0; i < numGaps; i++)
				characters += '-';
			return sequence + characters;
		}

		// sequence did not traverse in a path, find the next ordered node with the sequence in it
		if (next == null)
			for (int i = orderedNodeIds.indexOf(node.getID()); i < orderedNodeIds.size(); i++)
				if (nodes.get(orderedNodeIds.get(i)).getSeqCharMapping().keySet().contains(seqId)) {
					next = nodes.get(orderedNodeIds.get(i));
					break;
				}

		if (next == null) { // sequence finishes in this node
			// check if there are gaps at the end..
			int numGaps = orderedNodeIds.size() - orderedNodeIds.indexOf(node.getID()) - 1;
			characters = node.getSeqCharMapping().get(seqId) + "";
			for (int i = 0; i < numGaps; i++)
				characters += '-';
			return sequence + characters;
		}

		// identify the number of gaps between the nodes (based on the sorted node IDs)
		int numGaps = orderedNodeIds.indexOf(next.getID()) - orderedNodeIds.indexOf(node.getID()) - 1;
		for (int i = 0; i < numGaps; i++)
			characters += '-';

		return getGappySequence(seqId, next, sequence + characters, orderedNodeIds);
	}

	/**
	 * Get a histogram of the characters in the current node.
	 *
	 * @return	Map of base character to histogram count
	 */
	public Map<Character, MutableInt> getCurrentBaseCounts() {
		Map<Character, MutableInt> bases = new HashMap<>();
		for (Character b : current.getSeqCharMapping().values())
			if (!bases.containsKey(b))
				bases.put(b, new MutableInt());
			else
				bases.get(b).increment();
		return bases;
	}

	/**
	 * Get the out path of sequences from the current node.
	 *
	 * @return	Map of (next) nodeId and a list of sequences that traverse to that node
	 */
	public Map<Integer, List<Integer>> getSequencesOutEdges() {
		return getSequencesOutEdges(current);
	}

	/**
	 * Get the number of out edges of the current node.
	 *
	 * @return number of out edges
	 */
	public int getNumEdgesOut() { return current.getNextNodes().size(); }

	/**
	 * Get the number of in edges of the current node.
	 *
	 * @return number of in edges
	 */
	public int getNumEdgesIn() { return current.getPreviousNodes().size(); }

	/**
	 * Set the probability distribution of the current node.
	 *
	 * @param dist	Character probability distribution {(Character, Double)}
	 * @deprecated use setDistrib instead
	 */
	public void setCharacterDistribution(Map<Character, Double> dist) {
		current.setCharacterDistribution(dist);
	}

	public void setDistrib(EnumDistrib d) {
		current.setDistrib(d);
	}

	/**
	 * Get the character probability distribution of the current node.
	 *
	 * @return	probability distribution of the characters in the current node {(Character, Double)}
	 */
	public Map<Character, Double> getCharacterDistribution() {
		return current.getDistribution();
	}

	public Distrib getDistrib() {
		Node n = current;
		if (n == null)
			return null;
		return n.getDistrib();
	}
	public Distrib getDistrib(int index) {
		Node n = getNode(index);
		if (n == null)
			return null;
		return n.getDistrib();
	}
	/**
	 * Get the out path of sequences from the node.
	 *
	 * @param 	node	Node to get sequence out-path from
	 * @return	Map of (next) nodeId and a list of sequence IDs that traverse to that node
	 */
	private Map<Integer, List<Integer>> getSequencesOutEdges(Node node){
		Map<Integer, List<Integer>> nextNodeSeqs = new HashMap<>();
		for (Edge next : node.getNextTransitions()) {
			if (next.getNext() == finalNode)
				continue;
			if (!nextNodeSeqs.containsKey(next.getNext().getID()))
				nextNodeSeqs.put(next.getNext().getID(), new ArrayList<>());
			for (Integer seqId : next.getSequences())
				nextNodeSeqs.get(next.getNext().getID()).add(seqId);
		}
		return nextNodeSeqs;
	}

	/**
	 * Load provided partial order graph structure and sequences.
	 *
	 * @param structure		Dot representation of partial order graph or filepath to representation, if given aligned
	 *                      sequences (.aln, .fa, .fasta), constructs a graph from the aligned sequences
	 * @param seqPath		File path to sequences
	 */
	private void loadSequencesWithGraphStructure(String structure, String seqPath) throws FileNotFoundException, IOException{
		List<EnumSeq.Gappy<Enumerable>> seqs = new ArrayList<>();
		BufferedReader seqfile = new BufferedReader(new FileReader(seqPath));
		String line = null;
		line = seqfile.readLine();
		if (line.startsWith("CLUSTAL")) {
			seqs = EnumSeq.Gappy.loadClustal(seqPath, Enumerable.aacid_ext);
		} else if (line.startsWith(">")) {
			seqs = EnumSeq.Gappy.loadFasta(seqPath, Enumerable.aacid_ext, '-');
		} else {
			throw new RuntimeException("Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
		}
		seqfile.close();

		for (Integer seqId = 0; seqId < seqs.size(); seqId++)
			sequences.put(seqId, seqs.get(seqId).getName());

		if (structure.endsWith(".aln") || structure.endsWith(".fa") || structure.endsWith(".fasta")) {
			int seqLen = seqs.get(0).toString().length();
			for (Integer seqId = 0; seqId < seqs.size(); seqId++)
				// check that all sequences have the same length (i.e. check that they are aligned)
				if (seqs.get(seqId).toString().length() != seqLen)
					throw new RuntimeException("Aligned sequences must have the same length.");
			// load graph from aligned sequences
			current = loadPOGraph(seqs);
		} else
			// load graph from a dot file
			current = loadPOGraph(structure);
	}

	/**
	 * Uses Dijkstra's algorithm to calculate the minimum path distance (number of edges traversed) between the node
	 * with the given ID and every other node in the graph.
	 *
	 * @param nodeId	ID of starting node
	 * @return	minimum number of edges to traverse from node (ID) to every other node <node ID, distance>
	 */
	public Map<Integer, Integer> minDistance(int nodeId) {
		Map<Integer, Integer> edgeCount = new HashMap<>();
		final int MAXDIST = Integer.MAX_VALUE;

		// initialise distance
		for (Integer node : nodes.keySet())
			if (node == nodeId)
				edgeCount.put(node, 0);	// 0 distance from nodeId to nodeId
			else
				edgeCount.put(node, MAXDIST);

		findDist(nodes.get(nodeId), edgeCount);

		return edgeCount;
	}

	public Map<String, Object> getNextMapping() {
		Map<String, Object> extantNextIDs = new HashMap<>();
		for (Edge next : current.getNextTransitions())
			for (Integer seq : next.getSequences())
				extantNextIDs.put(sequences.get(seq), next.getNext().getID());
		return extantNextIDs;
	}

	public Map<String, Object> getPrevMapping() {
		Map<String, Object> extantPrevIDs = new HashMap<>();
		for (Edge next : current.getPreviousTransitions())
			for (Integer seq : next.getSequences())
				extantPrevIDs.put(sequences.get(seq), next.getNext().getID());
		return extantPrevIDs;
	}

	/**
	 * Returns ordered list of next transitions, order by extant support.
	 *
	 * @return	ordered list of next transitions
	 */
	public ArrayList<Integer> getOrderedNext() {
		Integer[] ids = new Integer[current.getNextTransitions().size()];
		Map<Integer, Integer> weights = new HashMap<>();
		int i = 0;
		for (Edge edge : current.getNextTransitions()) {
			ids[i++] = edge.getNext().getID();
			weights.put(edge.getNext().getID(), edge.getSequences().size());
		}
		quickSortIDsExtantSupport(ids, weights, true, 0, ids.length-1);

		return new ArrayList<>(Arrays.asList(ids));
	}

	/**
	 * Returns ordered list of previous transitions, order by extant support.
	 *
	 * @return	ordered list of previous transitions
	 */
	public ArrayList<Integer> getOrderedPrev() {
		Integer[] ids = new Integer[current.getPreviousTransitions().size()];
		Map<Integer, Integer> weights = new HashMap<>();
		int i = 0;
		for (Edge edge : current.getPreviousTransitions()) {
			ids[i++] = edge.getNext().getID();
			weights.put(edge.getNext().getID(), edge.getSequences().size());
		}
		quickSortIDsExtantSupport(ids, weights, true, 0, ids.length-1);

		return new ArrayList<>(Arrays.asList(ids));
	}

	/**
	 * Quick sort ids based on weights in the specified direction
	 *
	 * @param ids
	 * @param weights
	 * @param highToLow
	 * @param left
	 * @param right
	 */
	private void quickSortIDsExtantSupport(Integer[] ids, Map<Integer, Integer> weights, boolean highToLow, int left, int right) {
		// pivot is middle number between left and right

		if (ids.length <= 1)
			return;

		int pivotId = ids[left + (right - left) / 2];
		int pivot = weights.get(pivotId);

		int i = left;
		int j = right;

		while (i < j) {
			// find next array index (i) who's associated extant weight is less than the pivot (sorting high -> low)
			// if highToLow is true, otherwise greater than the pivot (sorting low -> high)
			// (moving from left to wards the pivot)
			while ((highToLow && weights.get(ids[i]) > pivot) || (!highToLow && weights.get(ids[i]) < pivot))
				i++;
			// find next array index (j) who's associated extant weight is greater/less than than the pivot (moving from
			// right towards the pivot)
			while ((highToLow && weights.get(ids[j]) < pivot) || (!highToLow && weights.get(ids[j]) > pivot))
				j--;
			// switch i and j over pivot
			Integer tmp = ids[i];
			ids[i] = ids[j];
			ids[j] = tmp;
			// increment indices and continue sorting
			i++;
			j--;
		}

		if (left < j)
			quickSortIDsExtantSupport(ids, weights, highToLow, left, j);
		if (i < right)
			quickSortIDsExtantSupport(ids, weights, highToLow, i, right);
	}

	/**
	 * Recursive function to calculate minimum distance to node.
	 *
	 * @param currNode	Node to cound distance from
	 * @param distance	Map to keep track of minimum distance for the node IDs
	 */
	private void findDist(Node currNode, Map<Integer,Integer> distance) {
		for (Node nextNode : currNode.getNextNodes())
			if (distance.get(nextNode.getID()) > distance.get(currNode.getID()) + 1)
				distance.put(nextNode.getID(), distance.get(currNode.getID()) + 1);
		for (Node nextNode : currNode.getNextNodes())
			findDist(nextNode, distance);
	}


	/**
	 * Gets the consensus sequences using an A star search algorithm.
	 *
	 * @param gappy
	 * @return
	 */
	public String getSupportedSequence(boolean gappy) {
		// FIXME: used currently for N0, needs to be updated to new consensus as in GraspCmd
		return null;
	}


	/**
	 * Comparator class for ensuring the nodes with the largest cost are
	 * kept at the front of the queue.
	 *
	 * Here, we want to maximise the cost as we are trying to get either:
	 * 		a. the longest sequence, or
	 * 		b. the path which the greatest number of sequences agree on.
	 */
	public class NodeComparator implements Comparator<Node>
	{
		@Override
		public int compare(Node x, Node y)
		{
			if (x.getCost() < y.getCost())
			{
				return -1;
			}
			if (x.getCost() > y.getCost())
			{
				return 1;
			}
			return 0;
		}
	}

	/**
	 * Helper class to store the path. Keeps track of the edge and the node
	 * that lead to a particular path. We need to keep track of the edge to
	 * be able to set it as 'consensus' and the node to be able to get the
	 * character.
	 */
	public class Path {
		private Node node;
		private Edge edge;

		public Path(Node node, Edge edge) {
			this.edge = edge;
			this.node = node;
		}

		public Node getNode() { return this.node; }

		public Edge getEdge() { return this.edge; }
	}

	/**
	 * Traverses the graph structure to construct the most supported sequence of characters.
	 *
	 * @return	most supported sequence of base characters
	 */
	public String getSupportedSequenceOld(boolean gappy) {

		String sequence = "";
		Node current = initialNode;
		while (current != finalNode) {
			current.setConsensus(true);
			// find next edge as first reciprocated edge in the ordered extant list, if no reciprocated, then
			// default to the first edge (extant support)
			Edge next = null;
			for (int n = 0; n < current.getNextTransitions().size(); n++) {
				if (current.getNextTransitions().get(n).reciprocated) {
					next = current.getNextTransitions().get(n);
					break;
				}
			}
			next = (next == null) ? current.getNextTransitions().get(0) : next;
			next.setConsensus(true);
			if (gappy) {
				Integer n = next.getNext().getID() - current.getID() - 1;
				if (n > 0)
					for (int g = 0; g < n; g++)
						sequence += '-';
			}
			if (next.getNext() != finalNode)
				sequence += next.getNext().getBase();
			current = next.getNext();
		}

		return sequence;
	}


	/**
	 * Get indication of if the current node is part of the consensus path.
	 *
	 * @return	indication of consensus membership
	 */
	public boolean getCurrentConsensusFlag() {
		return current.getConsensus();
	}

	/**
	 * Get the ID of the next node in the consensus path.
	 *
	 * @return	ID of the next consensus node
	 */
	public Integer getNextConsensusID() {
		for (Edge edge : current.getNextTransitions())
			if (edge.getConsensus())
				return edge.getNext().getID();
		return null;
	}

	/**
	 * Save partial order alignment graph in a dot format in the given directory.
	 *
	 * @param filepath	pathname to save dot product as
	 */
	public void saveToDot(String filepath) {
		try {
			DotWriter dw = new DotWriter(filepath + ".dot", "directed");
			// build rank lists
			List<Node> rankedNodes = new ArrayList<>();
			for (Node node : nodes.values())
				if (!rankedNodes.contains(node) && node.getAlignedNodes() != null) {
					String alignedNodes = "[";
					rankedNodes.add(node);
					for (Node alignedNode : node.getAlignedNodes()) {
						alignedNodes += alignedNode.getID() + " ";
						rankedNodes.add(alignedNode);
					}
					alignedNodes += "]";
					dw.writeGraphOption("rank", alignedNodes);
				}
			dw.writeGraphOption("rankdir", "\"LR\"");

			// if the base value of the node has not been instantiated, replace the base value with all possible values in the node
			Map<Node, String> nodeToLabel = getNodeLabels();
			for (Node node : nodes.values()) {
				if (node == finalNode)
					continue;
				// format node character distribution, if available
				EnumDistrib dist = node.getDistrib();
				String distStr = "\"";
				if (dist != null)
					for (Object base : dist.getDomain().getValues())
						distStr += base + ":" + String.format("%.0e", dist.get(base)) + " ";
				distStr = distStr.trim() + "\"";

				StringBuilder sb = new StringBuilder();

				// save node information
				for (Integer seqId : node.getSeqCharMapping().keySet())
					sb.append(sequences.get(seqId) + ":" + node.getSeqCharMapping().get(seqId) + ";");

				sb.replace(sb.length()-1, sb.length(),"");
				dw.writeNode(Integer.toString(node.getID()), "label", "\"" + nodeToLabel.get(node) + "\"", "fontsize", 15, "style", "\"filled\"", "fillcolor",
						"\"" + (node.getBase()==null?"#FFFFFF":dat.colourschemes.Clustal.getColour(node.getBase())) + "\"", "distribution", distStr, "sequences", "\"" + sb.toString() + "\"");

				// save edge information
				for (Edge out : node.getNextTransitions()) {
					if (out.getNext() == finalNode)
						continue;
					sb = new StringBuilder();
					for (Integer seqId : out.getSequences()){
						sb.append(sequences.get(seqId));
						sb.append(",");
					}
					sb.replace(sb.length()-1, sb.length(), "");
					float percent = (100.1f * out.getSequences().size() / sequences.size());
					dw.writeEdge(Integer.toString(node.getID()), Integer.toString(out.getNext().getID()), "fontsize", 12,
							"fontcolor", "darkgray", "penwidth", (out.getSequences().size() > 20 ? 8 : out.getSequences().size()/3 + 1), "dir", "forward", "label",
							String.format("\"%.0f", percent) + "%\"", "sequences", "\""+sb.toString()+"\"");
				}
			}
			dw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Load the partial order graph using the specified structure and sequences.
	 *
	 * @param structure		Dot representation of partial order graph or filepath to representation, if given aligned
	 *                      sequences (.aln, .fa, .fasta), constructs a graph from the aligned sequences
	 * @param seqPath		File path to sequences
	 */
	private void loadSequencesWithStructure(String structure, String seqPath) {
		List<EnumSeq.Gappy<Enumerable>> seqs = new ArrayList<>();
		try {
			BufferedReader seqfile = new BufferedReader(new FileReader(seqPath));
			String line = null;
			line = seqfile.readLine();
			if (line.startsWith("CLUSTAL")) {
				seqs = EnumSeq.Gappy.loadClustal(seqPath, Enumerable.aacid_ext);
			} else if (line.startsWith(">")) {
				seqs = EnumSeq.Gappy.loadFasta(seqPath, Enumerable.aacid_ext, '-');
			} else {
				throw new RuntimeException("Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
			}
			seqfile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (Integer seqId = 0; seqId < seqs.size(); seqId++)
			sequences.put(seqId, seqs.get(seqId).getName());

		if (structure.endsWith(".aln") || structure.endsWith(".fa") || structure.endsWith(".fasta")) {
			int seqLen = seqs.get(0).toString().length();
			for (Integer seqId = 0; seqId < seqs.size(); seqId++)
				// check that all sequences have the same length (i.e. check that they are aligned)
				if (seqs.get(seqId).toString().length() != seqLen)
					throw new RuntimeException("Aligned sequences must have the same length.");
			// load graph from aligned sequences
			current = loadPOGraph(seqs);
		} else
			// load graph from a dot file
			current = loadPOGraph(structure);
	}

	/**
	 * Generate a string representation of each node as the base character (if not empty), the sequence of unique characters
	 * stored in the node or "X" if all are empty.
	 *
	 * @return	Mapping of nodes and their labels
	 */
	private Map<Node, String> getNodeLabels() {
		Map<Node, String> nodeToLabel = new HashMap<>();
		for (Node node : nodes.values())
			if (node.getBase() == null) {
				String label = "";
				for (Character base : node.getSeqCharMapping().values())
					if (!label.contains(base.toString()))
						label += base;
				nodeToLabel.put(node, (node.getSeqCharMapping().values().isEmpty() ? "X" : label));
			} else {
				nodeToLabel.put(node, node.getBase()+"");
			}
		return nodeToLabel;
	}

	/**
	 * Recursively checks if the given sequence can be traversed in the partial order graph.
	 *
	 * @param sequence	sequence to find in the graph
	 * @return	if the sequence can be constructed by traversing the graph
	 */
	public boolean checkSequenceMembership(String sequence) {
		char[] bases = sequence.toCharArray();
		// find starting node
		Node next = null;
		for (Node node : nodes.values())
			if (node.getPreviousNodes().isEmpty() && node.getBase() == bases[0])
				next = node;
		// if next is still null, haven't found the starting character in the set of starting nodes
		if (next == null)
			return false;
		return findSequencePath(bases, next, 0);
	}

	/**
	 * Get the label of the current node. If a single character label is not populated, gets a string of the unique
	 * characters in the node.
	 *
	 * @return	node label
	 */
	public String getCurrentLabel() {
		Character base = this.getCurrentBase();
		String output = (base == null) ? "" : base + "";
		if (base == null)
			for (Character c : this.getCurrentBases())
				output += c;
		return output;
	}

	/**
	 * String representation of the partial order alignment graph.
	 *
	 * @return	String representation as dot format
	 */
	public String toString() {
		String sb = "digraph {\n\trankdir=\"LR\";\n";
		// build rank lists
		List<Node> rankedNodes = new ArrayList<>();
		for (Node node : nodes.values())
			if (!rankedNodes.contains(node) && node.getAlignedNodes() != null) {
				String alignedNodes = "[";
				rankedNodes.add(node);
				for (Node alignedNode : node.getAlignedNodes()) {
					alignedNodes += alignedNode.getID() + " ";
					rankedNodes.add(alignedNode);
				}
				alignedNodes += "]";
				sb += "rank="+alignedNodes+";\n";
			}
		// if the base value of the node has not been instantiated, replace the base value with all possible values in the node
		Map<Node, String> nodeToLabel = getNodeLabels();
		for (Node node : nodes.values()) {
			if (node == finalNode)
				continue;
			sb += "\t\"" + Integer.toString(node.getID()) + "\"[label=\"" + nodeToLabel.get(node) + "\"];\n";
			for (Edge edge : node.getNextTransitions()) {
				if (edge.getNext() == finalNode)
					continue;
				// find the number of sequences that traverse to the next node and calculate the weighting and percentage
				StringBuilder seqb = new StringBuilder();
				for (Integer seqId : edge.getSequences()) {
					seqb.append(sequences.get(seqId));
					seqb.append(",");
				}
				seqb.replace(seqb.length()-1, seqb.length(), "");
				float percent = (100.1f * edge.getSequences().size() / sequences.size());
				sb += "\t\"" + Integer.toString(node.getID()) + "\"->\"" + Integer.toString(edge.getNext().getID()) + "\"[dir=forward," + "label=" +
						String.format("\"%.1f", percent) + "%\"," + "sequences=\"" + seqb.toString() + "\"];\n";
			}
		}
		sb += "}";
		return sb;
	}


	/**
	 * Load partial order alignment structure from filename or string format
	 *
	 * @param 	structure string representation of structure or filepath to structure representation
	 * @return	first node in the graph
	 */
	private Node loadPOGraph(String structure) {
		String[] lines = null;			// array of lines if not reading from a file
		int lineCount = 0;				// counter for iterating through lines if not reading from a file
		BufferedReader reader = null;
		String line = "";
		if (structure.endsWith(".dot")) {
			try {
				reader = new BufferedReader(new FileReader(new File(structure)));
				line = reader.readLine();
			} catch (FileNotFoundException e) {
				throw new RuntimeException("Cannot find partial order graph file.");
//				System.exit(1);
			} catch (IOException e) {
				throw new RuntimeException("Incorrect partial order graph file format.");
//				System.exit(1);
			}
		} else {
			// structure is a string representation
			lines = structure.split("\n");
			line = lines[lineCount];
		}
		// parse graph structure
		try {
			// load all nodes
			while (line != null) {
				line = line.replace("\t", "");
				if (!line.contains("->")) {
					// node information, store distribution and sequence characters
					String[] elements = line.split("[\\[]+");
					if (elements.length > 1) {
						String nodeId = elements[0].replace("\"","");
						nodeId = nodeId.replaceAll("[^\\d]", "");
						int pogId = Integer.parseInt(nodeId);

						Node node = new Node(pogId);
						if (nodes.containsKey(pogId))
							node = nodes.get(pogId);
						else
							nodes.put(pogId, node);

						elements = elements[1].split("[,]+");
						for (String el : elements)
							if (el.contains("label")) {
								// load node label
								elements = el.split("[\"]+");
								String label = elements[1].replaceAll("\"", "");
								if (label.length() == 1)
									node.setBase(label.toCharArray()[0]);
							} else if (el.contains("distribution")) {
								// load node character distribution, expects "char:prob char:prob ... "
								elements = el.split("[\" ]+");
								Map<Object, Double> dist = new HashMap<>();
								for (String cp : elements)
									if (cp.contains(":")) {
										String[] els = cp.split("[:]+");
										dist.put(els[0].toCharArray()[0], Double.parseDouble(els[1]));
									}
								EnumDistrib d = new EnumDistrib(dist);
								//node.setCharacterDistribution(dist);
								node.setDistrib(d);
							} else if (el.contains("sequences")) {
								el = el.replace("\"", "");
								String seqs = el.split("sequences=")[1];
								for (String seq : seqs.split("[;]+")) {
									Integer seqId = getSequenceID(seq.split("[:]+")[0]);
									if (seqId == null) {
										seqId = sequences.size()+1;
										sequences.put(seqId, seq.split("[:]+")[0]);
									}
									node.addSequence(seqId, seq.split("[:]+")[1].toCharArray()[0]);
								}
							}
					}
				} else {
					// edge information, store pointers
					String[] elements = line.split("->");
					String id = elements[0].replace("\"","");
					id = id.replaceAll("[^\\d]", "");
					int fromId = Integer.parseInt(id);

					elements = elements[1].split("[\\[]+");
					id = elements[0].replace("\"","");
					id = id.replaceAll("[^\\d]", "");
					int toId = Integer.parseInt(id);

					Node nodeFrom = new Node(fromId);
					Node nodeTo = new Node(toId);

					// check if the node has already been defined, if so use that instead
					if (nodes.containsKey(fromId))
						nodeFrom = nodes.get(fromId);
					else
						nodes.put(fromId, nodeFrom);
					if (nodes.containsKey(toId))
						nodeTo = nodes.get(toId);
					else
						nodes.put(toId, nodeTo);

					// put sequences in the edge
					elements = elements[1].split("[\"]+");
					for (int el = 0; el < elements.length ; el++)
						if (elements[el].contains("sequences")) {
							elements = elements[el+1].split("[,]+");
							for (String seq : elements) {
								Integer seqId = getSequenceID(seq);
								if (seqId == null) {
									seqId = sequences.size()+1;
									sequences.put(seqId, seq);
								}
								nodeFrom.addNextNode(nodeTo, seqId);
								nodeTo.addPrevNode(nodeFrom, seqId);
							}
							break;
						}
				}
				if (reader != null)
					line = reader.readLine();
				else if (lineCount + 1 == lines.length)
					line = null;
				else
					line = lines[++lineCount];
			}

			if (reader != null)
				reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// find starting nodes
		for (Node node : nodes.values()) {
			if (node.getPreviousNodes().isEmpty())
				for (Integer seqId : node.getSeqIds()) {
					initialNode.addNextNode(node, seqId);
					node.addPrevNode(initialNode, seqId);
				}
			if (node.getNextTransitions().isEmpty())
				for (Integer seqId : node.getSeqIds()) {
					finalNode.ID = nodes.size();
					finalNode.addPrevNode(node, seqId);
					node.addNextNode(finalNode, seqId);
				}
		}

		return initialNode.getNextNodes().get(0);
	}

	/**
	 * Create PO graph based on a multiple sequence alignment .aln file
	 *
	 * @param seqs			list of sequences to load structure
	 */
	private Node loadPOGraph(List<EnumSeq.Gappy<Enumerable>> seqs) {
		// input:	seqLabel	sequence(with gap character)
		// each Node is a column of the aln file
		initialNode = new Node(-1);
		finalNode = new Node(seqs.get(0).toString().length());
		nodes = new HashMap<>();
		Map<Integer, List<Node>> seqNodeMap = new HashMap<>();

		int numNodes = seqs.get(0).toString().toCharArray().length;
		for (int nodeId = 0; nodeId < numNodes; nodeId++)
			nodes.put(nodeId, new Node(nodeId));

		// for each sequence, iterate through characters, if gap character, don't add to node
		for (int seqId = 0; seqId < seqs.size(); seqId++) {
			seqNodeMap.put(seqId, new ArrayList<>());

			char[] bases = seqs.get(seqId).toString().toCharArray();
			//List<Character> filteredBases = new ArrayList<>();
			for (int baseInd = 0; baseInd < bases.length; baseInd++)
				if (bases[baseInd] != '-') {
					//filteredBases.add(bases[baseInd]);
					nodes.get(baseInd).addSequence(seqId, bases[baseInd]);
					seqNodeMap.get(seqId).add(nodes.get(baseInd));
				}
			//Character[] chars = new Character[filteredBases.size()];
			//filteredBases.toArray(chars);
			//seqs.get(seqId).set(chars);
		}

		// iterate through the lists of nodes in the sequence mapping and assign edges to other nodes
		for (Integer seqId : seqNodeMap.keySet())
			for (int nodeInd = 0; nodeInd + 1 < seqNodeMap.get(seqId).size(); nodeInd++) {
				seqNodeMap.get(seqId).get(nodeInd).addNextNode(seqNodeMap.get(seqId).get(nodeInd + 1), seqId);
				addEdgeToSeq(seqNodeMap.get(seqId).get(nodeInd).ID, seqNodeMap.get(seqId).get(nodeInd + 1).ID, seqId );

				seqNodeMap.get(seqId).get(nodeInd + 1).addPrevNode(seqNodeMap.get(seqId).get(nodeInd), seqId);
			}

		List<Node> nodelist = new ArrayList<>(nodes.values());
		for (Node node : nodelist) {
			// if no sequences in node, remove (i.e. was a full gap column in aln)
			if (node.getSeqCharMapping().size() == 0) {
				setCurrent(node.getID());
				throw new RuntimeException("A column in your multiple sequence alignment contained only gaps. Remove this column or realign.");
			}

			// for all nodes, if only one unique character, set as base
			Character baseChar = null;
			for (Character b : node.getSeqCharMapping().values())
				if (baseChar == null)
					baseChar = b;
				else if (b != baseChar) {
					baseChar = null;
					break;
				}
			if (node.getBase() == null)
				node.setBase(baseChar);
		}

		// find starting and ending nodes
		for (Integer seqId : seqNodeMap.keySet()) {
			initialNode.addNextNode(seqNodeMap.get(seqId).get(0), seqId);
			seqNodeMap.get(seqId).get(0).addPrevNode(initialNode, seqId);
			finalNode.addPrevNode(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).size()-1), seqId);
			seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).size()-1).addNextNode(finalNode, seqId);

			addEdgeToSeq(initialNode.ID, seqNodeMap.get(seqId).get(0).ID, seqId );
			addEdgeToSeq(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).size()-1).ID, finalNode.ID, seqId );


		}

		return initialNode.getNextNodes().get(0);
	}


	/**
	 * Recursively populates node list.
	 *
	 * @param node	node to add to list
	 */
	private void addNode(Node node){
		if (node == finalNode)
			return;
		if (this.nodes.containsKey(node.getID()))
			return;
		this.nodes.put(node.getID(), node);
		for (Node next : node.getNextNodes())
			addNode(next);
	}

	/**
	 * Gets the ID of the extant sequence, i.e. position in the extant sequence list
	 *
	 * @param 	seqLabel	label of extant sequence
	 * @return	position of sequence in extant sequence list
	 */
	private Integer getSequenceID(String seqLabel) {
		for (Integer seqId : sequences.keySet())
			if (sequences.get(seqId).equalsIgnoreCase(seqLabel))
				return seqId;
		return null;
	}

	/**
	 * Recursively search through graph to identify if the sequence can be traversed fully
	 *
	 * @param sequence	sequence of characters to check
	 * @param node		current node to check
	 * @param baseInd	current character to check
	 * @return			if the current character is the current node and we can keep traversing the graph from this node
	 */
	private boolean findSequencePath(char[] sequence, Node node, int baseInd) {
		boolean found = node.getBase() == sequence[baseInd];
		boolean exploredNext = true;
		if (baseInd + 1 < sequence.length) {
			exploredNext = false;
			for (Node next : node.getNextNodes())
				if (next.getBase() == sequence[baseInd + 1]) {
					found &= findSequencePath(sequence, next, baseInd + 1);
					exploredNext = true;
				}
		}
		return found & exploredNext;
	}

	/**
	 * Performs a topological sort to ensure the Nodes are in the correct order
	 *
	 * @return 	list of IDs of sorted nodes
	 */
	public List<Integer> topologicalSort() {
		List<Integer> sortedList = new ArrayList<>();
		Set<Integer> completed = new HashSet<>();

		while (sortedList.size() < this.getNumNodes()) {
			Integer found = null;
			for (Integer node : this.nodes.keySet())
				if (!completed.contains(node)) {
					found = node;
					break;
				}
			if (found != null)
				depthFirstSearch(found, completed, sortedList);
		}

		return sortedList;
	}

	/**
	 * Performs a depth first search to get all subsequent Nodes from a Node in a PO Graph
	 *
	 * @param start id of Node to start search at
	 * @param completed set of Nodes that have already been sorted
	 * @param sortedList sorted order of Nodes
	 */
	private void depthFirstSearch(Integer start, Set<Integer> completed, List<Integer> sortedList) {
		List<Integer> successors = new ArrayList<>();
		List<Integer> stack = new ArrayList<>();
		stack.add(start);
		Set<Integer> started = new HashSet<>();

		while (!stack.isEmpty()) {
			Integer nodeID = stack.get(stack.size() - 1);
			stack.remove(stack.size() - 1);

			if (completed.contains(nodeID))
				continue;

			if (started.contains(nodeID)) {
				completed.add(nodeID);
				sortedList.add(0, nodeID);
				started.remove(nodeID);
				continue;
			}
			successors.clear();

			for (Node next : nodes.get(nodeID).getNextNodes())
				if (!next.getNextNodes().isEmpty() && !completed.contains(next.getID()))
					successors.add(0, next.getID());

			started.add(nodeID);
			stack.add(nodeID);
			stack.addAll(successors);
		}
	}

	
	private void addEdgeToSeq(Integer start, Integer end, Integer seqID){
		Map<Integer, Integer> edgeCount = new HashMap<>();
		edgeCount.put(end, 1);

		if (this.edgeCountsMSA.get(seqID) != null){

			this.edgeCountsMSA.get(seqID).put(start, edgeCount);

		}
		else {
			this.edgeCountsMSA.put(seqID, new HashMap(){{put(start, edgeCount);}});

		}


	}

	/**
	 * Find the shortest path from the start node to the end node using breadth first search.
	 *
	 * @param start		Node to start the search
	 * @param end		Node to finish
	 * @param ignore	Node to ignore (if applicable)
	 * @return			List of nodes to traverse
	 */
	private List<Node> findShortestPathToNode(Node start, Node end, Node ignore) {
		Queue<ArrayList<Node>> paths = new ArrayDeque<>();

		// create adjacency map of nodes, where adjacency is defined as nodes that are a 'next' node of the previous
		// considered node (or parent node)
		Map<Node, List<Node>> adjacency = new HashMap<>();
		ArrayList<Node> allNodes = new ArrayList<>();
		allNodes.add(initialNode);
		allNodes.addAll(nodes.values());
		allNodes.add(finalNode);

		for (Node node : allNodes) {
			if (!node.getNextNodes().isEmpty() && !adjacency.keySet().contains(node))
				adjacency.put(node, new ArrayList<>());
			for (Node next : node.getNextNodes()) {
				if (!adjacency.keySet().contains(next))
					adjacency.put(next, new ArrayList<>());
				for (Node adjacent : node.getNextNodes())
					if (!adjacency.get(next).contains(adjacent))
						adjacency.get(next).add(adjacent);
			}
		}

		List<Integer> sortedIds = getSortedIDs();

		ArrayList<Node> initialPath = new ArrayList<>();
		initialPath.add(start);
		paths.add(initialPath);

		ArrayList<Node> path = null;
		while (!paths.isEmpty()) {
			path = paths.remove();
			if (path.contains(end))
				break;
			// if there are no adjacent nodes, find a 'next' node
			if (path.get(path.size()-1).getNextNodes().isEmpty())
				continue;
			Node next = path.get(path.size()-1).getNextNodes().get(0);
			if (!adjacency.containsKey(next) || adjacency.get(next).isEmpty())
				adjacency.get(next);
			for (Node adjacent : adjacency.get(next)) {
				if (path.contains(adjacent) || adjacent.getNextNodes().isEmpty() || (!end.getNextNodes().isEmpty() && sortedIds.indexOf(adjacent.getID()) > sortedIds.indexOf(end.getID()))
						|| adjacent == ignore || (sortedIds.indexOf(adjacent.getID()) < sortedIds.indexOf(start.getID())))
					continue;
				ArrayList<Node> newpath = new ArrayList<>();
				newpath.addAll(path);
				newpath.add(adjacent);
				if (!newpath.contains(end))
					paths.add(newpath);
				else if (newpath.size() > 2)
					return newpath;
			}
		}
		if (!path.contains(end)) // couldn't find an alternative path to the node
			return null;
		return path;
	}

	private Node copyInitialNode(Node initial) {
		return initial.copy();
	}

	/**
	 * Node for encapsulating a partial order graph alignment of nodes. 'Aligned' POG nodes are combined to represent
	 * a single node.
	 */
	public class Node {
		private Integer ID = null;								// alignment ID
		private Character base;									// base character
		private List<Edge> nextTransitions;						// transitions to next nodes
		private List<Edge> prevTransitions;						// transitions to previous nodes
		private List<Node> alignedTo = null;					// list of nodes that are aligned with this node
		private Map<Integer, Character> seqChars;				// map of sequence Ids and their base character
		private Map<Character, Double> distribution = null;		// probability distribution of inferred character (deprecated, use distrib instead)
		private EnumDistrib distrib = null;						// probability distribution of character state
		private boolean consensus = false; 						// flag to indicate if belongs to the consensus path
		private long cost = 10000;								// cost to  reach this node from the start
		private Map<Integer, Integer> transitionCost; // Keeps track of the cost of a particular transition

		/**
		 * Constructor
		 */
		public Node() {
			this.prevTransitions = new ArrayList<>();
			this.nextTransitions = new ArrayList<>();
			this.seqChars = new HashMap<>();
			this.transitionCost = new HashMap<>();
		}

		/**
		 * Constructor for when we're recreating the objects from JSON objs.
		 *
		 * @param id
		 * @param base
		 * @param consensus
		 */
		public Node(int id, char base, boolean consensus) {
			this.consensus = consensus;
			this.ID = id;
			this.base = base;
			this.prevTransitions = new ArrayList<>();
			this.nextTransitions = new ArrayList<>();
			this.seqChars = new HashMap<>();
			this.transitionCost = new HashMap<>();
		}


		/**
		 * Gets the maximum parsimony score of transitioning from one node to another.
		 * @param nodeId
		 * @return
		 */
		public Integer getTransitionCost(Integer nodeId) {
			return transitionCost.get(nodeId);
		}

		/**
		 * Constructor
		 *
		 * @param ID Node id
		 */
		public Node(Integer ID) {
			this();
			this.ID = ID;
		}

		/**
		 * Sets the cost based on a heuristic.
		 * @param cost
		 */
		public void setCost(long cost) {
			this.cost = cost;
		}

		/**
		 * Gets the cost to this node.
		 * @return
		 */
		public long getCost() {
			return this.cost;
		}

		/**
		 * Returns a deep copy of the node structure.
		 *
		 * @param map	Map structure to keep track of visited nodes
		 * @return		deep copy of node
		 */
		private Node copy(Map<Node, Node> map) {
			Node copy = map.get(this);
			if (copy == null) {
				copy = new Node(this.getID());
				copy.setBase(this.base);
				if (this.getAlignedNodes() != null)
					for (Node node : this.alignedTo)
						copy.addAlignedNode(node.copy(map));
				for (Integer sc : this.seqChars.keySet())
					copy.addSequence(sc, seqChars.get(sc));
				for (Node next : this.getNextNodes())
					next.copy(map);
				map.put(this, copy);
			}
			return copy;
		}

		public Node copy() {
			Map<Node,Node> map = new IdentityHashMap<>();
			initialNode.copy(map);
			// link prev/next pointers
			for (Node fromOriginal : map.keySet())
				for (Edge toEdge : fromOriginal.getNextTransitions()) {
					Node fromCopy = map.get(fromOriginal);
					Node toCopy = map.get(toEdge.getNext());
					for (Integer seqId : toEdge.getSequences()) {
						fromCopy.addNextNode(toCopy, seqId);
						toCopy.addPrevNode(fromCopy, seqId);
					}
				}
			// find node with ID
			for (Node node : map.values())
				if (node.getID().equals(this.getID()))
					return node;
			return null;
		}

		/**
		 * Add transition from this to the given next node.
		 *
		 * @param next	Next node where transition points
		 * @param seq	Sequences that traverse to the next node
		 */
		private void addNextNode(Node next, Integer seq) {
			Edge nextT = null;
			for (Edge edge : nextTransitions)
				if (edge.getNext().getID().equals(next.getID())) {
					// edge already exists, just add sequence
					nextT = edge;
					this.nextTransitions.remove(nextT);
					break;
				}
			// edge doesn't already exist, create new transition
			if (nextT == null)
				nextT = new Edge(next);
			if (seq != null)
				nextT.addSequence(seq);
			// find index to put next edge
			// order transitions based on #seqs
			int index = 0;
			int seqCount = nextT.getSequences().size();
			for (int i = 0; i < this.nextTransitions.size(); i++)
				if (seqCount < this.nextTransitions.get(i).getSequences().size())
					index = i+1;
			this.nextTransitions.add(index, nextT);
		}

		/**
		 * Add transition from this to the given previous node.
		 *
		 * @param prev	Previous node where transition points
		 * @param seq	Sequence that traverses to the next node
		 */
		private void addPrevNode(Node prev, Integer seq) {
			for (Edge edge : prevTransitions)
				if (edge.getNext().getID().equals(prev.getID())) {
					// edge already exists, just add sequence
					if (seq != null)
						edge.addSequence(seq);
					return;
				}
			// edge doesn't already exist, create new transition
			Edge prevT = new Edge(prev);
			if (seq != null)
				prevT.addSequence(seq);
			this.prevTransitions.add(prevT);
		}

		/**
		 * Get the ID of the alignment
		 *
		 * @return	Alignment ID
		 */
		private Integer getID(){
			return this.ID;
		}

		/**
		 * Get a list of nodes that are aligned with this node.
		 *
		 * @return	List of aligned nodes. Null if not applicable.
		 */
		public List<Node> getAlignedNodes() { return this.alignedTo; }

		/**
		 * Add an aligned node.
		 *
		 * @param node	aligned node
		 */
		public void addAlignedNode(Node node) {
			if (alignedTo == null)
				alignedTo = new ArrayList<>();
			alignedTo.add(node);
		}

		/**
		 * Set the probability distribution of characters in this node.
		 *
		 * @param dist	distribution {(Character, probability)}
		 * @deprecated use setDistrib instead
		 */
		public void setCharacterDistribution(Map<Character, Double> dist) {
			this.distribution = new HashMap<>(dist);
			if (this.distrib == null)
                this.distrib = new EnumDistrib(Enumerable.aacid);
            for (Character ch : dist.keySet())
                if (this.distrib.isValid(ch))
                    this.distrib.set(ch, dist.get(ch));
		}

		public void setDistrib(EnumDistrib d) {
			this.distrib = d;
			// below is a FIX to work with GRASP GUI
            if (this.distribution == null) {
                this.distribution = new HashMap<>();
                for (Object ch : this.distrib.getDomain().getValues()) {
                    this.distribution.put((Character)ch, this.distrib.get(ch));
                }
            }
		}

		/**
		 * Get the distribution of characters in this node.
         * Actually, it SETs the distribution from seqChars if it has not be set before...
		 *
		 * @return distribution {(Character, probability)}
		 * @deprecated use getDistrib instead
		 */
		public Map<Character, Double> getDistribution() {
		    // MB: I think the first if-clause is when the distribution has been set from a marginal recon, or if the code below has been executed previously
			if (distribution != null)
				return distribution;
			// MB: I think the below code is ONLY executed if this is an extant node
			distribution = new HashMap<>();
			for (Character b : seqChars.values())
				if (distribution.containsKey(b))
					distribution.put(b, distribution.get(b) + 1.0);
				else
					distribution.put(b, 1.0);
			for (Character b : distribution.keySet())
				distribution.put(b, distribution.get(b) / seqChars.size());
			if (distrib == null) {
			    distrib = new EnumDistrib(Enumerable.aacid);
			    for (Object ch : distrib.getDomain().getValues())
			        if (distribution.containsValue(ch))
			            distrib.set(ch, distribution.get(ch));
            }
			return this.distribution;
		}

		public EnumDistrib getDistrib() {
			return distrib;
		}

		/**
		 * Remove the previous node from the list.
		 *
		 * @param prev	previous node
		 */
		private void removePrevNode(Node prev){
			for (Edge edge : prevTransitions)
				if (edge.getNext().getID().equals(prev.getID())) {
					prevTransitions.remove(edge);
					prev.removeNextNode(this);
					return;
				}
		}

		/**
		 * Remove the next node from the list.
		 *
		 * @param next	next node
		 */
		private void removeNextNode(Node next){
			for (Edge edge : nextTransitions)
				if (edge.getNext().getID().equals(next.getID())) {
					nextTransitions.remove(edge);
					next.removePrevNode(this);
					return;
				}
		}

		/**
		 * Store sequence that traverses through this node.
		 *
		 * @param seqId	ID of extant sequence
		 * @param base	Base character of sequence
		 */
		private void addSequence(int seqId, char base){
			for (Integer sqId : seqChars.keySet())
				if (sqId == seqId)
					return;
			this.seqChars.put(seqId, base);
		}

		/**
		 * Set the consensus flag.
		 *
		 * @param flag	Indication of whether node belongs to the consensus path.
		 */
		private void setConsensus(boolean flag) { this.consensus = flag; }

		/**
		 * Get the consensus flag.
		 *
		 * @return	flag indicating consensus membership
		 */
		private boolean getConsensus() { return this.consensus; }

		/**
		 * Get list of next nodes.
		 *
		 * @return	list of next nodes
		 */
		private List<Node> getNextNodes(){
			List<Node> nodes = new ArrayList<>();
			for (Edge edge : nextTransitions)
				nodes.add(edge.getNext());
			return nodes;
		}

		/**
		 * Get list of previous nodes.
		 *
		 * @return	list of previous nodes
		 */
		private List<Node> getPreviousNodes() {
			List<Node> nodes = new ArrayList<>();
			for (Edge edge : prevTransitions)
				nodes.add(edge.getNext());
			return nodes;
		}

		/**
		 * Get the list of in edges.
		 *
		 * @return	List of in transition edges.
		 */
		private List<Edge> getPreviousTransitions() {
			return this.prevTransitions;
		}

		/**
		 * Get the list of out edges.
		 *
		 * @return	List of out transition edges.
		 */
		private List<Edge> getNextTransitions() {
			return this.nextTransitions;
		}

		/**
		 * Set the inferred base character of this Node.
		 *
		 * @param base	inferred base character
		 */
		private void setBase(Character base) {
			this.base = base;
		}


		/**
		 * Get the inferred base character.
		 *
		 * @return	inferred base character
		 */
		private Character getBase(){ return base; }


		/**
		 * Generates string representation of the Node.
		 *
		 * @return		String representation in reduced dot format
		 */
		public String toString() {
			String sb = "\"" + Integer.toString(ID) + "\"" + "[label=\"" + base + "\"];\n";
			for (Node nextNode : getNextNodes())
				if (!nextNode.getNextNodes().isEmpty())
					sb += "\"" + Integer.toString(ID) + "\"->\"" + Integer.toString(nextNode.getID())+ "\"" + "[]\n;";
			return sb;
		}

		/**
		 * Gets a list of the sequence IDs that traverse this node.
		 *
		 * @return	List of sequence IDs
		 */
		private List<Integer> getSeqIds(){
			List<Integer> ids = new ArrayList<>();
			for (Integer sc : seqChars.keySet())
				ids.add(sc);
			return ids;
		}

		/**
		 * Get the mapping between sequence ID and base character.
		 *
		 * @return	Mapping of sequence ID and base character for the current node
		 */
		private Map<Integer,Character> getSeqCharMapping(){
			Map<Integer, Character> mapping = new HashMap<>();
			for (Integer sc : seqChars.keySet())
				mapping.put(sc, seqChars.get(sc));
			return mapping;
		}
	}

	/**
	 * Edge for storing transitions between nodes. This is required to easily track sequence paths when edges and nodes
	 * are removed from the graph.
	 */
	public class Edge {
		private Node next = null;
		private boolean consensus = false;
		private boolean reciprocated = false;
		private List<Integer> sequences;

		// Used once we already have the sequences and we no longer retain this info, instead converting
		// it to a weight.
		private double weight = 0.0;

		public Edge() {
			this.sequences = new ArrayList<>();
		}

		public Edge(Node nextNode) {
			this();
			this.next = nextNode;
		}

		/**
		 * Initialiser for when we already have the edge saved.
		 * @param nextNode
		 * @param weight
		 * @param reciprocated
		 */
		public Edge(Node nextNode, double weight, boolean reciprocated) {
			this.next = nextNode;
			this.weight = weight;
			this.reciprocated = reciprocated;
		}

		public void setConsensus(boolean flag) { this.consensus = flag; }

		public boolean getConsensus() { return this.consensus; }

		public void setReciprocated(boolean flag) { this.reciprocated = flag; }

		public boolean getReciprocated() { return this.reciprocated; }

		private void addSequence(int seqId) {
			this.sequences.add(seqId);
		}

		private void removeSequence(int seqId) {
			this.sequences.remove((Integer)seqId);
		}

		private List<Integer> getSequences() {
			return this.sequences;
		}

		public Node getNext() {
			return this.next;
		}

		public void setNext(Node nextNode) {
			this.next = nextNode;
		}
	}

}