/**
 * 
 */
package dat.file;

import dat.EnumSeq;
import dat.Enumerable;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * Reader for ClustalW multiple alignment files (.aln)
 * @author m.boden
 */
public class AlnReader {

    final File file;
    final Enumerable alpha;
    BufferedReader reader;

    /**
     * Construct a reader for Clustal files.
     *
     * @param filename the name of the file
     * @param alpha the alphabet that the alignment uses
     */
    public AlnReader(String filename, Enumerable alpha) throws IOException {
        this.file = new File(filename);
        reader = new BufferedReader(new FileReader(file));
        this.alpha = alpha;
    }

    /**
     * Loads a file on the Clustal format, requiring a header starting with
     * "CLUSTAL", then expects at least one empty line before a block of
     * (partial) alignments. The longest space- or tab-delimited string on a
     * line (excluding the first) is the alignment. The tokens preceding the
     * alignment is the label of the sequence. Later blocks are used to extend
     * the alignments. Here, the labels should be in agreement with the ones
     * given in the first block OR if not, appear in the same order.
     *
     * @return the sequences annotated with alignment data
     * @throws java.io.IOException
     * @see seq.file.SeqReader#load()
     */
    public EnumSeq.Gappy[] load() throws IOException {
        Map<String, StringBuffer> seq = new HashMap<>();
        Map<Integer, String> order = new HashMap<>();
        boolean collectNames = true; // true if still finding new line names
        int rowCount = 0;
        int precEmptyLineCount = 0;
        boolean uniqueNames = true;
        /* Example lines in an .aln file:
         CLUSTAL W (1.83) multiple sequence alignment


         P25101          ---METLCLRASFWLALVG---CVISDN---PERYSTNLSNHVDDFTTFRGTELS-----
         P24530          MQPPPSLCGRALVALVLACGLSRIWGEERGFPPDRATPLLQTAEIMTPPTKTLWPKGSNA
         P30556          ---------------------------------------MILNSSTEDG-----------
         :             

         P25101          FLVTTHQPTNLVLP--SNGS-----MHNYCPQQTKITSAFKYINTVISCTIFIVGMVGNA
         P24530          SLARSLAPAEVPKGDRTAGSPPRTISPPPCQGPIEIKETFKYINTVVSCLVFVLGIIGNS
         P30556          --IKRIQDDCPK---------------------AGRHNYIFVMIPTLYSIIFVVGIFGNS
         :::: :: :*: *: **:
         */
        String line = reader.readLine();
        String titleLine = null;
        int count = 0;
        int rowInBlock = 0;
        while (line != null) {
            rowCount++;
            line = line.trim();
            if (line.length() < 1) { // empty line
                precEmptyLineCount++;
            } else {
                if (precEmptyLineCount > 0) {
                    rowInBlock = 0;
                } else {
                    rowInBlock++;
                }
                if (titleLine == null) {
                    titleLine = line.toUpperCase();
                    if (!titleLine.startsWith("CLUSTAL")) {
                        throw new RuntimeException("Not a CLUSTAL file: \"" + file.getAbsolutePath() + "\". First row should contain \"CLUSTAL\".");
                    }
                } else {
                    if (Character.isLetter(line.charAt(0)) || Character.isDigit(line.charAt(0))) { 			// this is a proper line starting with an identifier
                        StringTokenizer tok = new StringTokenizer(line, " \t");
                        int nTokens = tok.countTokens();
                        String[] allTokens = new String[nTokens];
                        int longestIndex = 1;
                        for (int i = 0; i < nTokens; i++) {
                            allTokens[i] = tok.nextToken();
                            if (i > 0) { // the sequence can't be the first token so don't even consider it even if it's the longest string on the row
                                if (allTokens[i].length() >= allTokens[longestIndex].length()) {
                                    longestIndex = i;
                                }
                            }
                        }
                        if (nTokens > 1) { 								// this line has contents... 
                            // the first token(s) constitute the sequence identifier (not necessarily unique)
                            String strid = allTokens[0];
                            for (int i = 1; i < longestIndex; i++) { // we use all tokens up until the one that has most letter (this is the alignment)
                                strid = strid.concat(" " + allTokens[i]);
                            }
                            strid = strid + "\\" + rowInBlock; 			// check out the first token on the line
                            String seqstr = allTokens[longestIndex];
                            StringBuffer oldbuf = seq.get(strid);
                            if (oldbuf == null) { // Did not find identifier
                                // Could still be there. Look more carefully...
                                String orig = order.get(new Integer(rowInBlock + 1)); // this is the name that appeared on this row in the previous block (if there's one)
                                if (orig == null) {
                                    if (collectNames) { // we're still in the first block so must be a new identifier
                                        count++;
                                        seq.put(strid, new StringBuffer(seqstr));
                                        order.put(count, strid);
                                    } else {
                                        throw new RuntimeException("Invalid identifier \"" + strid + "\" on row " + rowCount);
                                    }
                                } else { // did find it
                                    collectNames = false;
                                    oldbuf = seq.get(orig);
                                    oldbuf.append(seqstr);
                                }
                            } else { // found the identifier first attempt 
                                collectNames = false;
                                oldbuf.append(seqstr);
                            }
                        }
                    }
                }
                precEmptyLineCount = 0;
            }
            line = reader.readLine();
        }
        EnumSeq.Gappy[] aseqs = new EnumSeq.Gappy[seq.size()];
        String[] names = new String[seq.size()];
        int i = 0;
        List<Integer> sorted = new ArrayList<>(order.keySet());
        Collections.sort(sorted);
        for (Integer index : sorted) {
            names[i++] = order.get(index);
        }
        for (int a = 0; a < aseqs.length; a++) {
            StringBuffer strbuf = seq.get(names[a]);
            if (strbuf != null) {
                int idx = names[a].indexOf("\\");
                String id = names[a];
                if (idx > 0) {
                    id = names[a].substring(0, idx);
                }
                idx = id.indexOf("/");
                if (idx > 0) {
                    id = id.substring(0, idx);
                }
                String syms = strbuf.toString();
                List<Character> symlist = new ArrayList<>();
                for (int index = 0; index < strbuf.length(); index++) {
                    Character sym = strbuf.charAt(index);
                    if (sym == '-')
                        symlist.add(null);
                    else if (alpha.isValid(sym))
                        symlist.add(sym);
                    else
                        throw new RuntimeException("Invalid symbol " + sym + " in sequence " + id + " at index " + index);
                }
                Character[] symarr = new Character[symlist.size()];
                symlist.toArray(symarr);
                // remove the numbering from the name
                aseqs[a] = new EnumSeq.Gappy<>(alpha);
                aseqs[a].set(symarr);
                aseqs[a].setName(id);
            }
        }
        return aseqs;
    }

    /**
     * Closes a file. 
     * @throws IOException if the file can not be closed
     */
    public void close() throws IOException {
        reader.close();
    }
    

    /**
     * @param args
     */
    public static void main(String[] args) {
        String filename = "/Users/mikael/Desktop/cyp3.aln";
        if (args.length > 0)
            filename = args[0];
        try {
            AlnReader aln = new AlnReader(filename, Enumerable.aacid);
            EnumSeq.Gappy[] e = aln.load();
            int nameLen = 0;
            for (int i = 0; i < e.length; i++) {
                if (nameLen < e[i].getName().length()) {
                    nameLen = e[i].getName().length();
                }
            }

            for (int i = 0; i < e.length; i++) {
                StringBuffer empty = new StringBuffer("                                                       ");
                empty.setLength(nameLen - e[i].getName().length());
                
                System.out.println(e[i].getName() + empty.toString() + "\t[" + e[i].length() + "]\t" + e[i].toString());
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

}
