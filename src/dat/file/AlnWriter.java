package dat.file;

import dat.EnumSeq;

import java.io.*;

/**
 * Created by marnie on 24/10/16.
 */
public class AlnWriter {
    public static int LINE_WIDTH = 60;
    private BufferedWriter writer=null;
    private boolean fileExists = false;

    /**
     * Create a file for storing sequences.
     *
     * @param filename the name of the file
     * @throws IOException if the operation fails
     */
    public AlnWriter(String filename) throws IOException  {
        this(new File(filename));
    }

    /**
     * Create a file for storing sequences.
     *
     * @param file the file to be opened for writing
     * @throws IOException if the operation fails
     */
    public AlnWriter(File file) throws IOException {
        try {
            fileExists = file.exists();
            writer = new BufferedWriter(new FileWriter(file));
        } catch (IOException e) {
            writer=null;
            throw new IOException(e.getMessage());
        }
    }

    /**
     * Retrieves the (Buffered)Writer that has been opened for writing.
     *
     * @return the Writer to which sequences can be written
     */
    protected BufferedWriter getWriter() {
        return writer;
    }

    /**
     * Checks if this file existed previously.
     *
     * @return true if the file existed before the current save.
     */
    public boolean exists() {
        return fileExists;
    }

    /**
     * Saves the sequences to the file.
     *
     * @param collection the collection of sequences
     * @throws IOException if the write operation fails
     */
    public void save(EnumSeq[] collection) throws IOException {
        try {
            writer.write("CLUSTAL");
            writer.newLine();
            writer.newLine();
            int line = 0;
            while (line < collection[0].length()) {
                for (int s = 0; s < collection.length; s++) {
                    EnumSeq seq = collection[s];
                    writer.write(seq.getName() + "\t");
                    char[] sequence = seq.toString().toCharArray();
                    for (int i = line; i < line + LINE_WIDTH && i < sequence.length; i++)
                        writer.write(sequence[i]);
                    writer.newLine();
                }
                line += LINE_WIDTH;
                writer.newLine();
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException("Error writing sequences.");
        }
    }

    /**
     * Closes the file so that it can be read by others.
     *
     * @throws IOException if the close operation fails
     */
    public void close() throws IOException {
        try {
            writer.flush();
            writer.close();
        } catch (IOException e) {
            throw new IOException("Error while closing.");
        }
    }
}
