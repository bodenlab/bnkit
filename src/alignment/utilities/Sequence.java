package alignment.utilities;

/**
 * Class for modelling a sequence with id and character information
 */
public class Sequence {

    private String id;
    private String seq;

    public Sequence(String id, String seq){
        this.id = id;
        this.seq = seq;
    }

    String getID(){
        return this.id;
    }

    public String getSeq(){
        return this.seq;
    }
}
