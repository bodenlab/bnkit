package vis;

public class Defines {

    /**
     * Defines for the Nodes in the MSA.
     *
     * Each position in an array corresponds to a particular value. An array as such is used
     * to save on space in the database and front end.
     */
    public static final int G_LABEL = 0;
    public static final int G_ID = 1;
    public static final int G_X = 2;
    public static final int G_MUTANTS = 3;
    public static final int G_CONSENSUS = 4;
    public static final int G_SEQ = 5;
    public static final int G_GRAPH = 6;

    // For sequence values
    public static final int G_VALUE = 1;

    // Where we place the characters in the sequence object / array
    public static final int G_CHAR = 0;

    /**
     * Defines for the Edges.
     */
    public static final int E_CONSENSUS = 0;
    public static final int E_RECIPROCATED = 1;
    public static final int E_FROM = 2;
    public static final int E_TO = 3;
    public static final int E_WEIGHT = 4;
    public static final int E_SINGLE = 5;

    /**
     * General Defines.
     */
    public static final int TRUE = 1;
    public static final int FALSE = 0;
}
