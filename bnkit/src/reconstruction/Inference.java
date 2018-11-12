package reconstruction;

import java.util.List;
import json.JSONArray;

public class Inference {
    /**
     * Helper class to store changes to an ancestral graph node
     *
     * Information:
     * 		- POG structure index
     * 		- Inferred base character: base character that is inferred or '-' to represent a gap (i.e. that the node needs to be deleted when updating the structure)
     */
    public Integer pogId;
    public char base;
    public List<Integer> transitions;

    public Inference(Integer id, char ch, List<Integer> tr){
        pogId = id;
        base = ch;
        transitions = tr;
    }
    public String toString(){
        return pogId + "->" + base;
    }

    /**
     * Return a JSON representation so we can save it to the database
     * @return
     */
    public JSONArray getAsJSON() {
        JSONArray infDetails = new JSONArray();
        infDetails.put(pogId);
        infDetails.put(base);
        JSONArray transitionsArr = new JSONArray();
        for (Integer transition : transitions)
            transitionsArr.put(transition);
        infDetails.put(transitionsArr);
        return infDetails;
    }
}
