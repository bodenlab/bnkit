package dat.pog;

import asr.ASRException;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class POAGraph extends IdxEdgeGraph<SeqEdge> {

    private Enumerable domain;

    public POAGraph(Enumerable domain, int nNodes) {
        super(nNodes, false, true);
        this.domain = domain;
    }

    public POAGraph(dat.EnumSeq.Alignment<Enumerable> aln) {
        this(aln.getDomain(), aln.getWidth());
        for (int i = 0; i < aln.getWidth(); i ++) {
            // before adding, make sure there is at least one sequence with content in this column
            if (aln.getOccupancy(i) > 0) {
                EnumNode node = addNode(i); // Create the node
                for (int j = 0; j < aln.getHeight(); j ++) {
                    EnumSeq.Gappy<Enumerable> gseq = aln.getEnumSeq(j);
                    Object sym = gseq.get(i);
                    if (sym != null) { // when a character is present, count it
                        node.add(sym);
                        node.setXLabel(i + 1);
                    }
                }
            }
        }
        for (int j = 0; j < aln.getHeight(); j ++) {
            EnumSeq.Gappy<Enumerable> gseq = aln.getEnumSeq(j);
            int start = 0, from = -1, to = -1;
            if (isTerminated()) {
                for (; start < aln.getWidth(); start ++)
                    if (gseq.get(start) != null)
                        break;
                SeqEdge e = this.getEdge(-1, start);
                if (e == null) {
                    e = new SeqEdge(domain);
                    this.addEdge(-1, start, e);
                    e.setTotal(aln.getHeight());
                }
                e.add(gseq);
            }
            for (int i = start; i < aln.getWidth(); i ++) {
                if (gseq.get(i) == null)
                    continue;
                else if (to == -1)
                    to = i;
                else {
                    from = to;
                    to = i;
                    SeqEdge e = this.getEdge(from, to);
                    if (e == null) {
                        e = new SeqEdge(domain);
                        this.addEdge(from, to, e);
                        e.setTotal(aln.getHeight());
                    }
                    e.add(gseq, gseq.get(from), gseq.get(to));
                }
            }
            if (isTerminated()) {
                from = to;
                to = this.maxsize();
                SeqEdge e = this.getEdge(from, to);
                if (e == null) {
                    e = new SeqEdge(domain);
                    this.addEdge(from, to, e);
                    e.setTotal(aln.getHeight());
                }
                e.add(gseq);
            }
        }
    }

    public Enumerable getDomain() {
        return domain;
    }

    /**
     * Modify the graph by adding a node at a specified index.
     * Note: does not connect the node.
     * @param nid node index, which must be valid, i.e. between 0 and N - 1, where N is the size of the possible/valid node-set
     */
    public synchronized EnumNode addNode(int nid) {
        EnumNode node = new EnumNode(this.domain);
        addNode(nid, node);
        return node;
    }

    /**
     * Modify the graph by adding a node.
     * Note: does not connect the node.
     * @param node the node
     * @return the index assigned to the node, if indices are exhausted, a runtime exception is thrown
     */
    public synchronized int addNode(EnumNode node) {
        return super.addNode(node);
    }

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        out.println("Usage: asr.GRASP \n" +
                "\t[-a | --aln <filename>]\n" +
                "\t[-o | --out <filename>]\n" +
                "\t{-n | --nwk <filename>}\n" +
                "\t{--ancestor <ancestor-label>}"
        );
        if (msg != null)
            out.println(msg);
        System.exit(error);
    }

    public static void main(String[] args) {
        String ALN_FILE = null;
        String NWK_FILE = null;
        String OUT_FILE = null;
        Integer MARG_NODE = null;
        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if ((arg.equalsIgnoreCase("-aln") || arg.equalsIgnoreCase("a")) && args.length > a + 1) {
                    ALN_FILE = args[++a];
                } else if ((arg.equalsIgnoreCase("-nwk") || arg.equalsIgnoreCase("n")) && args.length > a + 1) {
                    NWK_FILE = args[++a];
                } else if ((arg.equalsIgnoreCase("-out") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    OUT_FILE = args[++a];
                } else if (arg.equalsIgnoreCase("-ancestor") && args.length > a + 1) {
                    String ancid = args[++a];
                    if (ancid.startsWith("N"))
                        ancid = ancid.substring(1);
                    try {
                        MARG_NODE = Integer.parseInt(ancid);
                    } catch (NumberFormatException e) {
                        usage(2, args[a] + " is not a valid ancestor name (use <number>, or \"N<number>\", where <number> starts with 0 at root, depth-first). Tip: perform joint reconstruction first to check branch point numbering in tree.");
                    }
                }
            }
        }
        if (ALN_FILE == null)
            usage(1, "No alignment file given");
        if (OUT_FILE == null)
            usage(3, "No output file given");
        if (NWK_FILE == null && MARG_NODE != null)
            usage(2, "No tree file given, so cannot determine where ancestor \"" + MARG_NODE + "\" is");
        try {
            EnumSeq.Alignment aln = Utils.loadAlignment(ALN_FILE, Enumerable.aacid);
            if (NWK_FILE != null) {
                List<EnumSeq.Gappy> select = new ArrayList<>();
                Tree tree = Utils.loadTree(NWK_FILE);
                Utils.checkData(aln, tree);
                int bpidx = 0; // default root
                if (MARG_NODE != null)
                bpidx = tree.getIndex(MARG_NODE);
                if (bpidx >= 0) {
                    String[] names = aln.getNames();
                    for (int idx : tree.getLeaves(bpidx)) {
                        Object label = tree.getLabel(idx);
                        for (int i = 0; i < names.length; i ++) {
                            if (names[i].equals(label.toString())) {
                                EnumSeq.Gappy seq = aln.getEnumSeq(i);
                                select.add(seq);
                            }
                        }
                    }
                    aln = new EnumSeq.Alignment(select);
                }
            }
            POAGraph poag = new POAGraph(aln);
            poag.saveToDOT(OUT_FILE);
        } catch (IOException e) {
            usage(5, e.getMessage());
        } catch (ASRException e) {
            usage(4, e.getMessage());
        }
    }
}
