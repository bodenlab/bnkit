package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;

import java.io.IOException;

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

    public static void main(String[] args) {
        try {
            EnumSeq.Alignment aln = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/dp16_poag.aln", Enumerable.aacid));
            POAGraph poag = new POAGraph(aln);
            poag.saveToDOT("/Users/mikael/simhome/ASR/dp16_poag.dot");
            poag.saveToMatrix("/Users/mikael/simhome/ASR/dp16_poag.m");
        } catch (IOException e) {
            System.err.println(e);
        }
    }
}
