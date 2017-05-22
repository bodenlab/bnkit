package dat.file;


import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;

/**
 * Writes graph structure to a dot file.
 *
 * Created by marnie on 26/09/2016.
 */
public class DotWriter {

    private enum GRAPHTYPE { DIRECTED, UNDIRECTED };
    private BufferedWriter writer = null;
    private GRAPHTYPE graphType = null;

    /**
     * Create a file for storing the graph structure.
     *
     * @param filename      filepath and name
     * @param type          type of graph to write, "digraph" for directed graph, "graph" for undirected
     * @throws IOException
     */
    public DotWriter(String filename, String type) throws IOException {
        this(new File(filename), type);
    }

    /**
     * Write graph option.
     *
     * @param parameter    graph parameter to set
     * @param option       setting
     * @throws IOException
     */
    public void writeGraphOption(Object parameter, Object option) throws IOException {
        writer.write("\t" + parameter + "=");
        writer.write(option + ";");
        writer.newLine();
    }

    /**
     * Open a writer for the given file.
     *
     * @param file          File to write to
     * @param type          type of graph to write, "directed" for directed graph, "undirected" for undirected
     * @param options       list of options for the graph
     * @throws IOException
     */
    public DotWriter(File file, String type, Object... options) throws IOException {
        if (type.equalsIgnoreCase("directed"))
            graphType = GRAPHTYPE.DIRECTED;
        else
            graphType = GRAPHTYPE.UNDIRECTED;
        try {
            writer = new BufferedWriter(new FileWriter(file));
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
        if (graphType == GRAPHTYPE.DIRECTED)
            writer.write("digraph ");
        else
            writer.write("graph ");
        writer.write("{");
        writer.newLine();
    }

    /**
     * Write node element.
     *
     * @param id    graph node ID
     */
    public void writeNode(Object id) throws IOException {
        writeNode(id, null);
    }

    /**
     * Write node element information. Options is a list of information to include, e.g.:
     *      label, fontname, style, fillcolor, fontsize, etc.
     * See dot file formatting for a comprehensive list: http://www.graphviz.org/doc/info/attrs.html
     *
     * @param id        graph node ID
     * @param options   List of options for the node
     */
    public void writeNode(Object id, Object... options) throws IOException {
        writer.write("\t\"" + id + "\"");
        writeOptions(options);
        writer.write(";");
        writer.newLine();
    }

    /**
     * Write edge between nodes. Options is a list of information to include, e.g.:
     *      label, fontname, style, fillcolor, fontsize, etc.
     * See dot file formatting for a comprehensive list: http://www.graphviz.org/doc/info/attrs.html
     *
     * @param idFrom    ID of initial node
     * @param idTo      ID of final node
     * @param options   List of options for the graph edge
     */
    public void writeEdge(Object idFrom, Object idTo, Object... options) throws IOException {
        String from = "";
        try {
            from = Integer.toString((int)idFrom);
        } catch (Exception e ){
            from = idFrom.toString();
        }
        String to = "";
        try {
            to = Integer.toString((int)idTo);
        } catch (Exception e ){
            to = idTo.toString();
        }
        writer.write("\t\"" + from + "\"");
        if (graphType == GRAPHTYPE.DIRECTED)
            writer.write("->");
        else
            writer.write("--");
        writer.write("\"" + to + "\"");
        writeOptions(options);
        writer.write(";");
        writer.newLine();
    }


    /**
     * Close the writer stream.
     *
     * @throws IOException
     */
    public void close() throws IOException {
        writer.write("}");
        try {
            writer.flush();
            writer.close();
        } catch (IOException e) {
            throw new IOException("Error while closing");
        }
    }

    /**
     * Iterates though the list of options and writes to the file.
     *
     * @param options       list of options
     * @throws IOException
     */
    private void writeOptions(Object... options) throws IOException {
        if (options == null)
            return;
        writer.write("[");
        StringBuilder sb = new StringBuilder();
        for (int opt = 0; opt < options.length; opt += 2) {
            sb.append(options[opt]);
            sb.append("=");
            sb.append(options[opt+1]);
            sb.append(", ");
        }
        sb.replace(sb.length()-2,sb.length(), "");
        writer.write(sb.toString());
        writer.write("]");
    }

}
