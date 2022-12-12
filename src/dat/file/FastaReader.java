/*
 * bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014  M. Boden et al.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package dat.file;

/**
 * @author mikael
 */

import dat.EnumSeq;
import dat.Enumerable;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FastaReader {
    /* FASTA defline
     * fields are normally separated with '|' (PIPE), however, to avoid a complicated process of identifying which is the label-field
     * the separaor we use is SPACE, TAB and COMMA. See below for examples for various databanks
        GenBank                           gi|gi-number|gb|accession|locus
        EMBL Data Library                 gi|gi-number|emb|accession|locus
        DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
        NBRF PIR                          pir||entry
        Protein Research Foundation       prf||name
        SWISS-PROT                        sp|accession|name
        Brookhaven Protein Data Bank (1)  pdb|entry|chain
        Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
        Patents                           pat|country|number 
        GenInfo Backbone Id               bbs|number 
        General database identifier       gnl|database|identifier
        NCBI Reference Sequence           ref|accession|locus
        Local Sequence identifier         lcl|identifier
    */
    
    /** Separator between sequences, default the FASTA separator */
    public static String SEQ_FILE_SEPARATOR="\\>";

    private List<EnumSeq> sequences=null; // holds all sequences in store

    public Enumerable getAlphabet() {
        return alphabet;
    }

    private Enumerable alphabet=null; // the alphabet type by which all sequences conform to

    /** if an error occurs, the current line number will be stored here. */
    public int CURRENT_LINE=0;
    
    /** the reference to the file which is currently read */
    private BufferedReader reader=null;
    private Character gapSymbol = null;
    public boolean CONVERT_TO_UPPERCASE = true;

    /**
     * The alphabets that will be tested for compliance if sequences are read without supplied alphabet.
     * Presented in order of preference.
     */
    public static Enumerable[] checkAlphabetsInOrder = new Enumerable [] {
            Enumerable.nacid,
            Enumerable.nacidwn,
            Enumerable.nacidRNA,
            Enumerable.nacidwnRNA,
            Enumerable.aacid,
            Enumerable.aacidwx};

    /**
     * Opens a file on the FASTA format.
     * This constructor makes an educated guess on what the alphabet is; it checks in order
     * DNA, DNA with N, RNA, RNA with N, Amino Acid and the resorts to create its own.
     * It will recognise a (separate) gap character ('-' by default), and not include it as part of the alphabet.
     * @param filename the filename
     * @throws IOException if the file can not be found
     */
    public FastaReader(String filename) throws IOException {
        this(filename, '-');
    }

    /**
     * Takes an already open file i.e. a buffered reader and does the same as above.
     *
     * This constructor makes an educated guess on what the alphabet is; it checks in order
     * DNA, DNA with N, RNA, RNA with N, Amino Acid and the resorts to create its own.
     * It will recognise a (separate) gap character ('-' by default), and not include it as part of the alphabet.
     * @param br the buffered reader
     * @throws IOException if the file can not be found
     */
    public FastaReader(BufferedReader br) throws IOException {
        this(br, '-');
    }

    /**
     * Opens a file on the FASTA format.
     * This constructor makes an educated guess on what the alphabet is; it checks in order
     * DNA, DNA with N, RNA, RNA with N, Amino Acid and the resorts to create its own.
     * It will recognise a (separate) gap character ('-' by default), and not include it as part of the alphabet.
     * @param reader the buffered reader
     * @param gapSymbol the gap character
     * @throws IOException if the file can not be found
     */
    public FastaReader(BufferedReader reader, Character gapSymbol) throws IOException {
        this.gapSymbol = gapSymbol;
        try {
            String line = reader.readLine();
            Set<Object> charset = new HashSet<>();
            while (line != null) {
                String myline = line.trim();
                if (!myline.startsWith(">"))
                    for (Character ch : myline.toCharArray())
                        charset.add(ch);
                line = reader.readLine();
            }
            reader.close();
            for (Enumerable alpha : checkAlphabetsInOrder) {
                boolean valid = true;
                for (Object ch : charset) {
                    if (ch != gapSymbol && !Enumerable.aacid.isValid(ch)) {
                        valid = false;
                        break;
                    }
                }
                if (valid) { // we are happy with this alphabet
                    this.alphabet = alpha;
                    break;
                }
            }
            if (this.alphabet == null)
                this.alphabet = new Enumerable(charset.toArray());
            sequences=new ArrayList<>();
        } catch (FileNotFoundException e) {
            reader=null;
            throw new IOException(e.getMessage());
        }
    }

    /**
     * Opens a file on the FASTA format.
     * This constructor makes an educated guess on what the alphabet is; it checks in order
     * DNA, DNA with N, RNA, RNA with N, Amino Acid and the resorts to create its own.
     * It will recognise a (separate) gap character ('-' by default), and not include it as part of the alphabet.
     * @param filename the filename
     * @param gapSymbol the gap character
     * @throws IOException if the file can not be found
     */
    public FastaReader(String filename, Character gapSymbol) throws IOException {
        this.gapSymbol = gapSymbol;
        try {
            FileReader freader=new FileReader(filename);
            reader=new BufferedReader(freader);
            String line = reader.readLine();
            Set<Object> charset = new HashSet<>();
            while (line != null) {
                String myline = line.trim();
                if (!myline.startsWith(">"))
                    for (Character ch : myline.toCharArray())
                        charset.add(ch);
                line = reader.readLine();
            }
            reader.close();
            freader.close();
            for (Enumerable alpha : checkAlphabetsInOrder) {
                boolean valid = true;
                for (Object ch : charset) {
                    if (ch != gapSymbol && !Enumerable.aacid.isValid(ch)) {
                        valid = false;
                        break;
                    }
                }
                if (valid) { // we are happy with this alphabet
                    this.alphabet = alpha;
                    break;
                }
            }
            if (this.alphabet == null)
                this.alphabet = new Enumerable(charset.toArray());
            freader=new FileReader(filename);
            reader=new BufferedReader(freader);
            sequences=new ArrayList<>();
        } catch (FileNotFoundException e) {
            reader=null;
            throw new IOException(e.getMessage());
        }
    }

    /**
     * Assigns the buffered reader for use in GRASP.
     *
     * @param br the open file already buffered.
     * @param alpha the symbol alphabet to use
     * @throws IOException if the file can not be found
     */
    public FastaReader(BufferedReader br, Enumerable alpha, Character gapSymbol) throws IOException {
        sequences=new ArrayList<>();
        alphabet=alpha;
        this.gapSymbol = gapSymbol;
        reader = br;
    }


    /**
     * Opens a file on the FASTA format.
     * @param filename the filename
     * @param alpha the symbol alphabet to use
     * @throws IOException if the file can not be found
     */
    public FastaReader(String filename, Enumerable alpha) throws IOException {
        this(new File(filename), alpha, null);
    }

    /**
     * Opens a file on the FASTA format. 
     * @param filename the filename
     * @param alpha the symbol alphabet to use
     * @param gapSymbol the symbol used for gaps
     * @throws IOException if the file can not be found
     */
    public FastaReader(String filename, Enumerable alpha, Character gapSymbol) throws IOException {
        this(new File(filename), alpha, gapSymbol);
    }
    
    /**
     * Opens a file on the FASTA format. 
     * @param file the file to be read
     * @param alpha the symbol alphabet to use
     * @param gapSymbol the symbol used for gaps
     * @throws IOException if the file can not be found
     */
    public FastaReader(File file, Enumerable alpha, Character gapSymbol) throws IOException {
        sequences=new ArrayList<>();
        alphabet=alpha;
        this.gapSymbol = gapSymbol;
        try {
            FileReader freader=new FileReader(file);
            reader=new BufferedReader(freader);
        } catch (FileNotFoundException e) {
            reader=null;
            throw new IOException(e.getMessage());
        }
    }

    protected EnumSeq extract(String[] seqDef, boolean gappy) {
//        boolean gappy = false;
        if (seqDef.length > 0) {
            StringTokenizer tok = new StringTokenizer(seqDef[0], " \t,;");
            if (tok.hasMoreTokens()) {
                // find name and collect annotations from first line
                String first = tok.nextToken();
                String name = first.substring(1, first.length()); // sequence name
                String info = seqDef[0];
                // there's no universal way annotating a FASTA entry, but all info is found on the first line
                for (int a = 1; tok.hasMoreTokens(); a++) {
                    tok.nextToken(); // currently not used
                }
                // now add symbols
                List<Character> syms = new ArrayList<>();
                for (int i = 1; i < seqDef.length; i++) { // start on second row and onwards...
                    String parse = seqDef[i].trim();
                    if (parse.length() > 0) {
                        for (int s = 0; s < parse.length(); s++) {
                            Character current = parse.charAt(s);
                            if (CONVERT_TO_UPPERCASE) // default setting is true
                                current = Character.toUpperCase(current);  // FASTA allows symbols to be specified in lower case
                            if (current.equals(gapSymbol)) {
                                syms.add(null);
                                gappy = true;
                            } else if (alphabet.isValid(current)) { // sequence component
                                syms.add(current);
                            } else {
                                throw new RuntimeException("Unrecognised symbol: " + current + " at index " + s);

                            }
                        }
                    }
                }
                Character[] arr = new Character[syms.size()];
                syms.toArray(arr);
                EnumSeq seq;
                if (gappy) 
                    seq = new EnumSeq.Gappy(alphabet);
                else
                    seq = new EnumSeq(alphabet);
                seq.setName(name);
                seq.setInfo(info);
                seq.set(arr);
                return seq;
            } else {
                throw new RuntimeException("Invalid Sequence format (not proper FASTA)");
            }
        } else {
            throw new RuntimeException("Invalid Sequence format (not proper FASTA)");
        }
    }

    /**
     * Reads a file on the FASTA format. 
     * @return all sequences in the file
     * @throws IOException if the file can not be read
     */
    public EnumSeq[] load() throws IOException, RuntimeException {
        Pattern pattern = Pattern.compile(SEQ_FILE_SEPARATOR);
        String line;
        line = reader.readLine();
        CURRENT_LINE++;
        while (line != null) {
            Matcher matcher = pattern.matcher(line);
            if (matcher.find()) {
                List<String> def = new ArrayList<>();
                def.add(line);
                line = reader.readLine();
                CURRENT_LINE++;
                while (line != null) {
                    matcher = pattern.matcher(line);
                    if (matcher.find()) {
                        break;
                    }
                    String trimmed = line.trim();
                    if (trimmed.length() > 0) {
                        def.add(trimmed);
                    }
                    line = reader.readLine();
                    CURRENT_LINE++;
                }
                String[] seqdef = new String[def.size()];
                for (int i = 0; i < def.size(); i++) {
                    seqdef[i] = (String) def.get(i);
                }
                try {
                    EnumSeq sbuf = extract(seqdef, false);
                    sequences.add(sbuf);
                } catch (RuntimeException e) {
                    throw new RuntimeException("Sequence \"" + seqdef[0].substring(1) + "\" is using an invalid symbol. " + e.getMessage());
                }
            } else {
                line = reader.readLine();
                CURRENT_LINE++;
            }
        }
        EnumSeq[] seqarr = new EnumSeq[sequences.size()];
        sequences.toArray(seqarr);
        return seqarr;
    }

    /**
     * Reads a file on the FASTA format.
     * @return all sequences in the file
     * @throws IOException if the file can not be read
     */
    public EnumSeq.Gappy[] loadGappy() throws IOException {
        Pattern pattern = Pattern.compile(SEQ_FILE_SEPARATOR);
        String line;
        line = reader.readLine();
        CURRENT_LINE++;
        while (line != null) {
            Matcher matcher = pattern.matcher(line);
            if (matcher.find()) {
                List<String> def = new ArrayList<>();
                def.add(line);
                line = reader.readLine();
                CURRENT_LINE++;
                while (line != null) {
                    matcher = pattern.matcher(line);
                    if (matcher.find()) {
                        break;
                    }
                    String trimmed = line.trim();
                    if (trimmed.length() > 0) {
                        def.add(trimmed);
                    }
                    line = reader.readLine();
                    CURRENT_LINE++;
                }
                String[] seqdef = new String[def.size()];
                for (int i = 0; i < def.size(); i++) {
                    seqdef[i] = (String) def.get(i);
                }
                try {
                    EnumSeq sbuf = extract(seqdef, true);
                    //FIXME - untested
                    sequences.add(sbuf);
                } catch (RuntimeException e) {
                    throw new RuntimeException("Sequence \"" + seqdef[0].substring(1) + "\" is using an invalid symbol. " + e.getMessage());
                }
            } else {
                line = reader.readLine();
                CURRENT_LINE++;
            }
        }
        EnumSeq.Gappy[] seqarr = new EnumSeq.Gappy[sequences.size()];
        sequences.toArray(seqarr);
        return seqarr;
    }

    /**
     * Closes a file. 
     * @throws IOException if the file can not be closed
     */
    public void close() throws IOException {
        reader.close();
    }
    
    public static void main(String[] args) {
        try {
            FastaReader freader=new FastaReader("/Users/mikael/simhome/ChIA-PET/duplx_500.fa", Enumerable.nacid);
            EnumSeq[] seqs=freader.load();
            System.out.println("Loaded " + seqs.length + " sequences.");
            for (int i = 0; i < seqs.length; i ++) {
                System.out.println(seqs[i].getName());
                Object[] content = seqs[i].get();
                for (int j = 0; j < content.length; j ++) {
                    if (content[j] == null)
                        System.out.print('-');
                    System.out.print(content[j]);
                }
                System.out.println();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
