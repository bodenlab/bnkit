package dat.file;

import asr.ASRException;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.Tree;
import dat.pog.Edge;
import dat.pog.POGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * File utilities intended for GRASP in particular.
 */
public class Utils {

    public static <E extends Edge> E toEdge(String edgetype, Object[] input) {
        return null;
    }

    /**
     * Check out stuff about a file that can be read
     */
    public static class SneakPeek {
        public Format format = null;
        public Enumerable alphabet = null;
        /**
         * Check out this file
         * @param filename file
         * @throws IOException if the file cannot be opened or read
         */
        public SneakPeek(String filename) throws IOException {
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line = reader.readLine();
            // first scan until a line with content is found
            while (line != null) {
                if (line.trim().length() > 0)
                    break;
                line = reader.readLine();
            }
            // we make a call what general format it is by parsing header info
            if (line.startsWith("CLUSTAL"))
                format = Format.CLUSTAL;
            else if (line.startsWith(">"))
                format = Format.FASTA;
            else if (line.startsWith("#")) {
                String follows = line.substring(1).trim();
                if (follows.startsWith("STOCKHOLM"))
                    format = Format.STOCKHOLM;
                else if (follows.startsWith("NEXUS"))
                    format = Format.NEXUS;
            } else if (line.startsWith("("))
                format = Format.NEWICK;
            // if sequences, work out what alphabet
            if (format == Format.CLUSTAL || format == Format.FASTA) {
                // todo
            }
            reader.close();
        }

        public Format getFormat() {
            return format;
        }
    }

    /**
     * File formats
     */
    public enum Format {
        FASTA,
        CLUSTAL,
        STOCKHOLM,
        NEWICK,
        NEXUS
    }

    public static String format2string(Format format) {
        if (format == null) return "Invalid";
        switch (format) {
            case CLUSTAL:   return "CLUSTAL";
            case FASTA:     return "FASTA";
            case STOCKHOLM: return "STOCKHOLM";
            case NEWICK:    return "NEWICK";
            case NEXUS:     return "NEXUS";
            default:        return "Unknown";
        }
    }

    /**
     * Read the first line of a file to figure out its format.
     * @param filename the filename
     * @return file format
     * @throws IOException if the file can't be read
     */
    public static Format getFormat(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        return new SneakPeek(filename).getFormat();
    }

    public static dat.EnumSeq.Alignment<Enumerable> loadAlignment(String filename, Enumerable alphabet) throws IOException, ASRException {
        SneakPeek sp = new SneakPeek(filename);
        Format format = sp.getFormat();
        if (format == null) throw new ASRException("Format of alignment is unknown");
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        List<EnumSeq.Gappy<Enumerable>> seqs = null;         // List of sequences (characters)
        if (format == Format.CLUSTAL) {
            seqs = EnumSeq.Gappy.loadClustal(reader, alphabet);
        } else if (format == Format.FASTA) {
            seqs = EnumSeq.Gappy.loadFasta(reader, alphabet, '-');
        } else {
            throw new ASRException("Format is not an alignment. Detected format " + format2string(format));
        }
        try {
            EnumSeq.Alignment aln = new EnumSeq.Alignment(seqs);
            return aln;
        } catch (RuntimeException e) {
            throw new ASRException(e.getMessage());
        }
    }

    public static Tree loadTree(String filename) throws IOException, ASRException {
        SneakPeek sp = new SneakPeek(filename);
        Format format = sp.getFormat();
        if (format == null) throw new ASRException("Format of tree is unknown");
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        Tree tree = null;
        if (format == Format.NEWICK) {
            tree = Newick.load(filename);
        } else if (format == Format.NEXUS) {
            throw new ASRException("NEXUS file reader not implemented.");
        } else {
            throw new ASRException("Format is not a tree format. Detected format " + format2string(format));
        }
        return tree;
    }

    public static void checkData(EnumSeq.Alignment aln, Tree tree) throws ASRException {
        ASRException.FileIssues exception = new ASRException.FileIssues();
        // Check if there are duplicate extant node names in the phylogenetic tree
        // Duplicate extant node names not allowed - will influence reconstruction outcomes
        Set<String> names_in_tree = new HashSet<>();
        Set<String> duplicates_in_tree = new HashSet<>();
        for (int idx : tree.getLeaves()) {
            Object name = tree.getBranchPoint(idx).getID();
            if (names_in_tree.contains(name.toString()))
                duplicates_in_tree.add(name.toString());
            names_in_tree.add(name.toString());
        }
        if (duplicates_in_tree.size() > 0)
            exception.add("Duplicate identifiers in tree: ", duplicates_in_tree);
        // Check if there are duplicate sequence names in the extant sequences
        // Duplicate sequence names not allowed - will influence reconstruction outcomes
        Set<String> names_in_aln = new HashSet<>();
        Set<String> duplicates_in_aln = new HashSet<>();
        for (String name : aln.getNames()) {
            if (names_in_aln.contains(name))
                duplicates_in_aln.add(name);
            names_in_aln.add(name);
        }
        if (duplicates_in_aln.size() > 0)
            exception.add("Duplicate identifiers in alignment: ", duplicates_in_aln);
        // Check if the provided extant sequences match up to the provided tree
        // save sequence information in internal nodes of the phylogenetic tree
        Set<String> missing_in_tree = new HashSet<>();
        for (String seqname : names_in_aln)
            if (!names_in_tree.contains(seqname))
                missing_in_tree.add(seqname);
        if (missing_in_tree.size() > 0)
            exception.add("Identifiers found in alignment but missing in tree", missing_in_tree);
        Set<String> missing_in_aln = new HashSet<>();
        for (String nodename : names_in_tree)
            if (!names_in_aln.contains(nodename))
                missing_in_aln.add(nodename);
        if (missing_in_aln.size() > 0)
            exception.add("Identifiers found in tree but missing in alignment", missing_in_aln);
        // if any problems were flagged above, we'll "throw" the exception
        if (exception.isLoaded())
            throw new ASRException(exception);
    }


}
