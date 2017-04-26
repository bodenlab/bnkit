package dat;

import dat.file.AlnWriter;
import dat.file.DotWriter;
import dat.file.FastaWriter;

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
 * @author Gabe, Marnie
 *
 */
public class POGraph {
	private Map<Integer, Node> nodes;				// map of node ID and nodes in the structure (for fast node retrieval)
	private Node current;							// pointer to the current node
	private Map<Integer, List<Node>> seqNodeMap;	// map of sequence ID and the order of structure nodes it traverses
	private Node initialNode = null;				// initial node (null node), points to the first actual nodes as 'next' nodes
	private Map<Integer, String> sequences;			// map of sequence ID and sequence label

	/**
	 * Constructor to initialise empty graph.
	 */
	public POGraph() {
		sequences = new HashMap<>();
		initialNode = new Node();
		nodes = new HashMap<>();
		seqNodeMap = new HashMap<>();
		current = null;
	}

	/**
	 * Constructor to initialise partial order graph using a dot structure and loading the associated sequences from
	 * a fasta file.
	 *
	 * @param structure		Dot representation of partial order graph or filepath to representation, if given aligned
	 *                      sequences (.aln, .fa, .fasta), constructs a graph from the aligned sequences
	 * @param seqPath		File path to sequences
	 */
	public POGraph(String structure, String seqPath) {
		this();
		// load sequences
		loadSequencesWithGraphStructure(structure, seqPath);
	}

	/**
	 * Constructor to initialise a partial order graph using a dot structure.
	 *
	 * @param structure		Dot representation of partial order graph
	 */
	public POGraph(String structure) {
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
		IdentityHashMap<Node,Node> map = new IdentityHashMap<>();
		this.initialNode = copy.initialNode.copy(map);
		for (Node node : this.initialNode.getNextNodes())
			addNode(node);
		for (Integer seqId : copy.seqNodeMap.keySet()) {
			this.seqNodeMap.put(seqId, new ArrayList<>());
			for (Node node : copy.seqNodeMap.get(seqId))
				this.seqNodeMap.get(seqId).add(this.nodes.get(node.getID()));
		}
		this.current = initialNode.getNextNodes().get(0);
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
		char[] bases = sequence.toCharArray();
		sequences.put(id, label);
		List<Node> seqnodes = new ArrayList<>();
		for (int baseInd = 0; baseInd < bases.length; baseInd++) {
			if (baseInd >= nodeIds.size() || !setCurrent(nodeIds.get(baseInd))) { // node doesn't currently exist in the graph, create and add to graph
				// find next ID number
				current = new Node(nodes.size());
				nodes.put(current.getID(), current);
				if (baseInd == 0) {
					initialNode.addNextNode(current);
					current.addPrevNode(initialNode);
				}
			}
			current.addSequence(id, bases[baseInd]);
			seqnodes.add(current);
		}
		// update node pointers
		seqnodes.get(0).addNextNode(seqnodes.get(1));
		seqnodes.get(seqnodes.size()-1).addPrevNode(seqnodes.get(seqnodes.size()-2));
		for (int nodeind = 1; nodeind < seqnodes.size() - 1; nodeind++) {
			seqnodes.get(nodeind).addPrevNode(seqnodes.get(nodeind-1));
			seqnodes.get(nodeind).addNextNode(seqnodes.get(nodeind+1));
		}
		seqNodeMap.put(id, seqnodes);
	}

	/**
	 * Set the pointer of the current node reference to the node with the specified ID.
	 * 
	 * @param nodeID	Node ID to set as current
	 * @return			indication of whether setting the node succeeded or not
	 */
	public boolean setCurrent(int nodeID){
		current = nodes.get(nodeID);
		return (current != null);
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
	 * Get the out edge weights of the current node. Weights based on total sequence support.
	 *
	 * @return	mapping between next node ID and the edge weight to that node
	 */
	public Map<Integer, Double> getEdgeWeights(){
		HashMap<Integer, Double> edgeWeights = new HashMap<>();
		for (Integer nextId : getNextIDs())
			edgeWeights.put(nextId, 1.0 * getSequencesOutEdges().get(nextId).size() / sequences.size());
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
		for (Integer seqId : seqNodeMap.keySet()) {
			seqNodeMapping.put(seqId, new ArrayList<>());
			for (Node node : seqNodeMap.get(seqId))
				seqNodeMapping.get(seqId).add(node.getID());
		}
		return seqNodeMapping;
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
	 * Get the number of nodes in the graph.
	 *
	 * @return	number of nodes in the structure.
	 */
	public int getNumNodes(){ return nodes.size(); }

	/**
	 * Get the ID of the current node.
	 * 
	 * @return	id of the current node, -1 if no node is set
	 */
	public int getCurrentId(){ return (current == null) ? -1 : current.getID(); }

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
	public List<Integer> getPrevIDs(){
		if (current == null)
			return null;
		ArrayList<Integer> prevIDs = new ArrayList<>();
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
	 * Removes the current node from the graph.
	 */
	public void removeNode() {
		if (current == null)
			return;
		for (Integer seqId : current.getSeqIds()) {
			// set the previous node and next node pointers
			if (seqNodeMap.get(seqId).indexOf(current) == 0) {
				// first node, only affects next node pointers
				Node nextNode = seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1);
				nextNode.removePrevNode(current);
			} else if (seqNodeMap.get(seqId).indexOf(current) == seqNodeMap.get(seqId).size() - 1) {
				// final node, only affect previous node pointers
				Node prevNode = seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) - 1);
				prevNode.removeNextNode(current);
			} else {
				// somewhere in between, affects both previous and next node pointers
				Node prevNode = seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) - 1);
				Node nextNode = seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1);
				prevNode.addNextNode(nextNode);
				prevNode.removeNextNode(current);
				nextNode.addPrevNode(prevNode);
				nextNode.removePrevNode(current);
			}
			seqNodeMap.get(seqId).remove(current);
		}

		nodes.remove(current.getID(), current);
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
	 * Get a list of base characters in the current node (note: can have repeats).
	 *
	 * @return	List of base characters in current node.
	 */
	public List<Character> getCurrentBasesNotUnique() {
		return new ArrayList<>(current.getSeqCharMapping().values());
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
			for (Integer seqId : seqNodeMap.keySet()) {
				List<Node> seqNodes = seqNodeMap.get(seqId);
				// step through sorted nodes and add gap or character, depending on if the node appears in the
				// node list for the sequence
				List<Character> characters = new ArrayList<>();
				for (Integer nodeId : orderedNodeIds) {
					boolean found = false;
					for (Node node : seqNodes)
						if (node.getID() == nodeId) {
							characters.add(node.getSeqCharMapping().get(seqId));
							found = true;
							break;
						}
					if (!found)
						characters.add('-');
				}
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
				System.err.print("Incorrect file type. Must be 'fasta' or 'clustal'.");
		} catch (IOException e) {
			e.printStackTrace();
		}
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
	 */
	public void setCharacterDistribution(HashMap<Character, Double> dist) {
		current.setCharacterDistribution(dist);
	}

	/**
	 * Get the character probability distribution of the current node.
	 *
	 * @return	probability distribution of the characters in the current node {(Character, Double)}
	 */
	public Map<Character, Double> getCharacterDistribution() {
		return current.getDistribution();
	}

	/**
	 * Get the out path of sequences from the node.
	 *
	 * @param 	node	Node to get sequence out-path from
	 * @return	Map of (next) nodeId and a list of sequence IDs that traverse to that node
	 */
	private Map<Integer, List<Integer>> getSequencesOutEdges(Node node){
		Map<Integer, List<Integer>> nextNodeSeqs = new HashMap<>();
		Map<Integer, List<Integer>> seqNodeMap = getSequenceNodeMapping();
		for (Integer seqId : node.getSeqIds())
			for (Node nextNode : node.getNextNodes())
				// check if sequence is in the next node
				if (nextNode.getSeqIds().contains(seqId))
					// check that there are no nodes between node and nextNode for seqId
					for (int i = 0; i < seqNodeMap.get(seqId).size() - 1; i++)
						if (seqNodeMap.get(seqId).get(i) == node.getID() && seqNodeMap.get(seqId).get(i+1) ==
								nextNode.getID()) {
							if (!nextNodeSeqs.containsKey(nextNode.getID()))
								nextNodeSeqs.put(nextNode.getID(), new ArrayList<>());
							nextNodeSeqs.get(nextNode.getID()).add(seqId);
							break;
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
	private void loadSequencesWithGraphStructure(String structure, String seqPath) {
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
				throw new RuntimeException("Alignment should be in Clustal or Fasta format");
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
	 * Traverses the graph structure to construct the most supported sequence of characters.
	 *
	 * @return	most supported sequence of base characters
	 */
	public String getSupportedSequence() {
		int maxId = 0;
		for (Integer nodeId : nodes.keySet())
			if (nodeId > maxId)
				maxId = nodeId;
		String seq = "";
		current = null;

		while (current == null || !current.getNextNodes().isEmpty()) {
			// find the next node based on how many sequences traverse to that node
			HashMap<Node, Integer> nodeCount = new HashMap<>();
			for (Integer seqId : seqNodeMap.keySet())
				if (current == null) {
					if (!nodeCount.containsKey(seqNodeMap.get(seqId).get(0)))
						nodeCount.put(seqNodeMap.get(seqId).get(0), 0);
					nodeCount.put(seqNodeMap.get(seqId).get(0), nodeCount.get(seqNodeMap.get(seqId).get(0))+1);
				} else if (seqNodeMap.get(seqId).contains(current) && seqNodeMap.get(seqId).indexOf(current) + 1 < seqNodeMap.get(seqId).size()) {
					if (!nodeCount.containsKey(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1)))
						nodeCount.put(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1), 0);
					nodeCount.put(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1), nodeCount.get(seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(current) + 1))+1);
				}

			Node next = null;
			for (Node node : nodeCount.keySet())
				if (next == null || nodeCount.get(node) > nodeCount.get(next))
					next = node;

			current = next;
			seq += current.getBase();
		}
		return seq;
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
				// format node character distribution, if available
				HashMap<Character, Double> dist = node.getDistribution();
				String distStr = "\"";
				if (dist != null)
					for (Character base : dist.keySet())
						distStr += base + ":" + String.format("%.0e", dist.get(base)) + " ";
				distStr = distStr.trim() + "\"";
				StringBuilder sb = new StringBuilder();
				for (Integer seqId : seqNodeMap.keySet())
					if (seqNodeMap.get(seqId).contains(node))
						sb.append(sequences.get(seqId) + ":" + node.seqChars.get(seqId) + ";");
				sb.replace(sb.length()-1, sb.length(),"");
				dw.writeNode(Integer.toString(node.getID()), "label", "\"" + nodeToLabel.get(node) + "\"", "fontsize", 15, "style", "\"filled\"", "fillcolor",
							"\"" + (node.getBase()==0?"#FFFFFF":dat.colourschemes.Clustal.getColour(node.getBase())) + "\"", "distribution", distStr, "sequences", "\"" + sb.toString() + "\"");
				for (Node next : node.getNextNodes()) {
					// find the number of sequences that traverse to the next node and calculate the weighting and percentage
					sb = new StringBuilder();
					int numSeqs = 0;
					for (Integer seqId : seqNodeMap.keySet())
						if (seqNodeMap.get(seqId).contains(node) && seqNodeMap.get(seqId).indexOf(node)+1 < seqNodeMap.get(seqId).size() &&
								seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(node)+1) == next) {
							numSeqs++;
							sb.append(sequences.get(seqId));
							sb.append(",");
						}
					if (numSeqs > 0)
						sb.replace(sb.length()-1, sb.length(), "");
					float percent = (100.1f * numSeqs / seqNodeMap.keySet().size());
					dw.writeEdge(Integer.toString(node.getID()), Integer.toString(next.getID()), "fontsize", 12,
							"fontcolor", "darkgray", "penwidth", (numSeqs > 20 ? 8 : numSeqs/3 + 1), "dir", "forward", "label",
							String.format("\"%.0f", percent) + "%\"", "sequences", "\""+sb.toString()+"\"");
				}
			}
			dw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
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
			sb += "\t\"" + Integer.toString(node.getID()) + "\"[label=\"" + nodeToLabel.get(node) + "\"];\n";
			for (Node next : node.getNextNodes()) {
				// find the number of sequences that traverse to the next node and calculate the weighting and percentage
				StringBuilder seqb = new StringBuilder();
				int numSeqs = 0;
				for (Integer seqId : seqNodeMap.keySet())
					if (seqNodeMap.get(seqId).contains(node) && seqNodeMap.get(seqId).indexOf(node)+1 < seqNodeMap.get(seqId).size() &&
							seqNodeMap.get(seqId).get(seqNodeMap.get(seqId).indexOf(node)+1) == next) {
						numSeqs++;
						seqb.append(sequences.get(seqId));
						seqb.append(",");
					}
				if (numSeqs > 0)
					seqb.replace(seqb.length()-1, seqb.length(), "");
				float percent = (100.1f * numSeqs / seqNodeMap.keySet().size());
				sb += "\t\"" + Integer.toString(node.getID()) + "\"->\"" + Integer.toString(next.getID()) + "\"[dir=forward," + "label=" +
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
				System.err.println("Cannot find partial order graph file.");
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Incorrect partial order graph file format.");
				System.exit(1);
			}
		} else {
			// structure is a string representation
			lines = structure.split("\n");
			line = lines[lineCount];
		}
		// parse graph structure
		Map<String, Integer> inputNodeToPONode = new HashMap<>();				// Mapping of input nodes and which PO graph nodes they are stored in
		Map<Integer, Map<String, Character>> nodeSeqCharMap = new HashMap<>();	// Mapping of input nodes, the sequences in that node and their base characters from the aln
		try {
			// load all nodes
			while (line != null) {
				line = line.replace("\t", "");
				if (!line.contains("->")) {
					String[] elements = line.split("[\\[]+");
					if (elements.length > 1) {
						String nodeId = elements[0].replace("\"","");
						int pogId = Integer.parseInt(nodeId);
						HashMap<Character, Double> dist = null;
						Character base = null;
						elements = elements[1].split("[,]+");
						for (String el : elements)
							if (el.contains("label")) {
								// load node label
								elements = el.split("[\"]+");
								String label = elements[1].replaceAll("\"", "");
								if (label.length() == 1)
									base = label.toCharArray()[0];
							} else if (el.contains("distribution")) {
								// load node character distribution, expects "char:prob char:prob ... "
								elements = el.split("[\" ]+");
								dist = new HashMap<>();
								for (String cp : elements)
									if (cp.contains(":")) {
										String[] els = cp.split("[:]+");
										dist.put(els[0].toCharArray()[0], Double.parseDouble(els[1]));
									}
							} else if (el.contains("sequences")) {
								nodeSeqCharMap.put(pogId, new HashMap<>());
								el = el.replace("\"", "");
								String seqs = el.split("[=]")[1];
								for (String seq : seqs.split("[;]+"))
									nodeSeqCharMap.get(pogId).put(seq.split("[:]+")[0], seq.split("[:]+")[1].toCharArray()[0]);
							}
						nodes.put(pogId, new Node(pogId));
						nodes.get(pogId).setBase(base);
						if (dist != null)
							nodes.get(pogId).setCharacterDistribution(dist);
						inputNodeToPONode.put(nodeId, pogId);
					}
				}
				if (reader != null)
					line = reader.readLine();
				else if (lineCount + 1 == lines.length)
					line = null;
				else
					line = lines[++lineCount];
			}

			// update pointers and sequence character mappings using edge information
			if (reader != null) {
				reader.close();
				reader = new BufferedReader(new FileReader(new File(structure)));
				line = reader.readLine();
			} else {
				lineCount = 0;
				line = lines[++lineCount];
			}
			while (line != null) {
				line = line.replace("\t", "");
				if (line.contains("->")) {
					String[] elements = line.split("[->]+");
					String fromId = elements[0].replace("\"","");
					int fromNodeId = Integer.parseInt(fromId);
					elements = elements[1].split("[\\[]+");
					String toId = elements[0].replace("\"","");
					int toNodeId = Integer.parseInt(toId);
					inputNodeToPONode.put(toId, toNodeId);
					int toPOGID = inputNodeToPONode.get(toId);
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
								// TODO: Base character for MSA (if base is null)
								nodes.get(inputNodeToPONode.get(fromId)).addSequence(seqId, nodeSeqCharMap.get(fromNodeId).get(seq));
								nodes.get(inputNodeToPONode.get(toId)).addSequence(seqId, nodeSeqCharMap.get(toNodeId).get(seq));
								if (!seqNodeMap.containsKey(seqId))
									seqNodeMap.put(seqId, new ArrayList<>());
								if (!seqNodeMap.get(seqId).contains(nodes.get(fromNodeId)))
									seqNodeMap.get(seqId).add(nodes.get(fromNodeId));
								if (!seqNodeMap.get(seqId).contains(nodes.get(toPOGID)))
									seqNodeMap.get(seqId).add(nodes.get(toPOGID));
							}
							break;
						}
					nodes.get(fromNodeId).addNextNode(nodes.get(toNodeId));
					nodes.get(toNodeId).addPrevNode(nodes.get(fromNodeId));
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

		// sort node order for each sequence and add the seqChar mapping for the final nodes
		for (Integer seqId : seqNodeMap.keySet()) {
			List<Node> sorted = new ArrayList<>();
			List<Node> unsorted = seqNodeMap.get(seqId);
			while (!unsorted.isEmpty()) {
				Node current = unsorted.get(0);
				for (Node node : unsorted) {
					if (current.getPreviousNodes().contains(node))
						current = node;
					if (current.getPreviousNodes().isEmpty())
						// found first node, add to sorted list
						break;
				}
				sorted.add(current);
				unsorted.remove(current);
			}
			seqNodeMap.get(seqId).addAll(sorted);
		}

		// find starting nodes
		for (Node node : nodes.values())
			if (node.getPreviousNodes().isEmpty()) {
				initialNode.addNextNode(node);
				node.addPrevNode(initialNode);
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
		initialNode = new Node();
		nodes = new HashMap<>();
		seqNodeMap = new HashMap<>();

		int numNodes = seqs.get(0).toString().toCharArray().length;
		for (int nodeId = 0; nodeId < numNodes; nodeId++)
			nodes.put(nodeId, new Node(nodeId));

		// for each sequence, iterate through characters, if gap character, don't add to node
		for (int seqId = 0; seqId < seqs.size(); seqId++) {
			seqNodeMap.put(seqId, new ArrayList<>());

			char[] bases = seqs.get(seqId).toString().toCharArray();
			List<Character> filteredBases = new ArrayList<>();
			for (int baseInd = 0; baseInd < bases.length; baseInd++)
				if (bases[baseInd] != '-') {
					filteredBases.add(bases[baseInd]);
					nodes.get(baseInd).addSequence(seqId, bases[baseInd]);
					seqNodeMap.get(seqId).add(nodes.get(baseInd));
				}
			Character[] chars = new Character[filteredBases.size()];
			filteredBases.toArray(chars);
			seqs.get(seqId).set(chars);
		}

		// iterate through the lists of nodes in the sequence mapping and assign pointers in nodes
		for (List<Node> nodeSeqs : seqNodeMap.values()) {
			nodeSeqs.size();
			for (int nodeInd = 0; nodeInd < nodeSeqs.size(); nodeInd++)
				if (nodeInd + 1 < nodeSeqs.size()) {
					nodeSeqs.get(nodeInd).addNextNode(nodeSeqs.get(nodeInd + 1));
					nodeSeqs.get(nodeInd + 1).addPrevNode(nodeSeqs.get(nodeInd));
				}
		}

		for (Node node : nodes.values()) {
			// find starting nodes
			if (node.getPreviousNodes().isEmpty()) {
				initialNode.addNextNode(node);
				node.addPrevNode(initialNode);
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

		return initialNode.getNextNodes().get(0);
	}


	/**
	 * Recursively populates node list.
	 *
	 * @param node	node to add to list
	 */
	private void addNode(Node node){
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
				if (!completed.contains(next.getID()))
					successors.add(0, next.getID());

			started.add(nodeID);
			stack.add(nodeID);
			stack.addAll(successors);
		}
	}

	/**
     * Node for encapsulating a partial order graph alignment of nodes. 'Aligned' POG nodes are combined to represent
	 * a single node.
     */
    private class Node {
    	private Integer ID = null;								// alignment ID
    	private Character base;									// base character
    	private List<Node> prevNodes;							// list of neighbouring previous nodes
    	private List<Node> nextNodes;							// list of neighbouring next nodes
		private List<Node> alignedTo = null;					// list of nodes that are aligned with this node
		private HashMap<Integer, Character> seqChars;			// map of sequence Ids and their base character
		private HashMap<Character, Double> distribution = null;	// probability distribution of inferred character
    	
    	/**
    	 * Constructor
    	 */
    	public Node() {
			this.prevNodes = new ArrayList<>();
			this.nextNodes = new ArrayList<>();
			this.seqChars = new HashMap<>();
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
    	 * Returns a deep copy of the node structure.
    	 * 
    	 * @param map	Map structure to keep track of visited nodes
    	 * @return		deep copy of node
    	 */
		public Node copy(Map<Node, Node> map) {
			Node copy = map.get(this);
			if (copy == null) {
				copy = new Node(this.getID());
				copy.setBase(this.base);
				if (this.getAlignedNodes() != null)
					for (Node node : this.alignedTo)
						copy.addAlignedNode(node.copy(map));
				for (Integer sc : this.seqChars.keySet())
					copy.addSequence(sc, seqChars.get(sc));
				map.put(this, copy);
				for (Node node: this.prevNodes)
					copy.addPrevNode(node.copy(map));
				for (Node node: this.nextNodes)
					copy.addNextNode(node.copy(map));
			}
			return copy;
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
		 */
		public void setCharacterDistribution(HashMap<Character, Double> dist) {
			this.distribution = new HashMap<>(dist);
		}

		/**
		 * Get the distribution of characters in this node.
		 *
		 * @return distribution {(Character, probability)}
		 */
		public HashMap<Character, Double> getDistribution() {
			return this.distribution;
		}

    	/**
    	 * Checks if the previous node is already stored, if not, add to list of previous nodes.
    	 * 
    	 * @param prev		Pointer to the prev node
    	 */
    	private void addPrevNode(Node prev) {
    		for (Node previous : prevNodes)
    			if (previous.getID() == prev.getID())
    				return;
			prevNodes.add(prev);
    	}
    	
    	/**
    	 * Checks if the next node is already stored, if not, add to list of next nodes.
    	 * 
    	 * @param next		Pointer to the next node
    	 */
    	private void addNextNode(Node next) {
			for (Node nextN : nextNodes)
				if (nextN.getID() == next.getID())
					return;
			nextNodes.add(next);
    	}
    	
    	/**
    	 * Remove the previous node from the list.
    	 * 
    	 * @param prev	previous node
    	 */
		private void removePrevNode(Node prev){
    		for (Node node : prevNodes)
    			if (node.getID() == prev.getID()) {
    				prevNodes.remove(node);
    				return;
    			}
    	}
    	
    	/**
    	 * Remove the next node from the list.
    	 * 
    	 * @param next	next node
    	 */
		private void removeNextNode(Node next){
    		for (Node node : nextNodes)
    			if (node.getID() == next.getID()) {
    				nextNodes.remove(node);
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
    	 * Get list of next nodes.
    	 * 
    	 * @return	list of next nodes
    	 */
		private List<Node> getNextNodes(){ return this.nextNodes; }
    	
    	/**
    	 * Get list of previous nodes.
    	 * 
    	 * @return	list of previous nodes
    	 */
		private List<Node> getPreviousNodes(){
    		return this.prevNodes;
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
    	 *  Generates string representation of the Node.
    	 *  
    	 * @return		String representation in reduced dot format
    	 */
    	public String toString() {
    		String sb = "\"" + ((ID == null) ? "null" : Integer.toString(ID)) + "\"" + "[label=\"" + base + "\"];\n";
    		for (Node nextNode : getNextNodes())
    			sb += "\"" + ((ID == null) ? "null" : Integer.toString(ID)) + "\"->\"" + Integer.toString(nextNode.getID())+ "\"" + "[]\n;";
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
}
