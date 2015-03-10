/*
 * Created on 8/05/2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package dat.file;

import dat.EnumSeq;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;

/**
 * @author mikael
 */

public class FastaWriter {

    public static int LINE_WIDTH = 60;
    private BufferedWriter writer=null;
    private boolean fileExists=false;
    
    /**
     * Create a file for storing sequences.
     * @param filename the name of the file 
     * @throws IOException if the operation fails 
     */
    public FastaWriter(String filename) throws IOException  {
        this(new File(filename));
    }

    /**
     * Create a file for storing sequences.
     * @param file the file to be opened for writing
     * @throws IOException if the operation fails 
     */
    public FastaWriter(File file) throws IOException {
        fileExists=file.exists();
        try {
            FileWriter fwriter=new FileWriter(file);
            writer=new BufferedWriter(fwriter);
        } catch (IOException e) {
            writer=null;
            throw new IOException(e.getMessage());
        }
    }
    
    /**
     * Retrieves the (Buffered)Writer that has been opened for writing.
     * @return the Writer to which sequences can be written
     */
    protected BufferedWriter getWriter() {
        return writer;
    }
    
    /**
     * Checks if this file existed previously. Useful if one wants to caution the user of overwriting.
     * @return true if the file existed before the current save.
     */
    public boolean exists() {
        return fileExists;
    }
    
    /**
     * Writes the sequence name on FASTA format to a string.
     * @return the defline string (FASTA)
     */
    public static String defline(EnumSeq seq) {
        StringBuffer sbuf=new StringBuffer();
        sbuf.append(">"+seq.getName()+" ");
        return sbuf.toString().trim();
    }

    /**
     * Saves the sequences to the file.
     * It uses the Sequence.write method to generate the String[] to write to the file
     * So if you want to make your own version of the FASTA entry, extend Sequence and overwrite "write".
     * Don't forget to {@see close} the file after all sequences have been stored.
     * @param collection the collection of sequences
     * @throws IOException if the write operation fails
     */
    public void save(EnumSeq[] collection) throws IOException {
        for (int s=0; s<collection.length; s++) {
            Object[] str=collection[s].get();
            try { 
                writer.write(defline(collection[s]));
                for (int i=0; i<str.length; i++) {
                    if (i % LINE_WIDTH == 0)
                        writer.newLine();
                    writer.write(str[i].toString());
                }
                writer.newLine();
            } catch (IOException e) {
                throw new IOException("Error in writing sequence "+collection[s]+" (index "+s+")");
            }
        }
    }
    
    /**
     * Closes the file so that it can be read by others.
     * @throws BioSeqIOException if the close operation fails
     */
    public void close() throws IOException {
        try {
            writer.flush();
            writer.close();
        } catch (IOException e) {
            throw new IOException("Error while closing");
        }
    }
    
}
