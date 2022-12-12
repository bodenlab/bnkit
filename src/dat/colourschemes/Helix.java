package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Helix{


    private static final Map<String, String> helix;
    static {
        helix = new HashMap<String, String>();
        helix.put("A", "#E718E7");
        helix.put("R", "#6F906F");
        helix.put("N", "#1BE41B");
        helix.put("D", "#778877");
        helix.put("C", "#23DC23");
        helix.put("Q", "#926D92");
        helix.put("E", "#FF00FF");
        helix.put("G", "#00FF00");
        helix.put("H", "#758A75");
        helix.put("I", "#8A758A");
        helix.put("L", "#AE51AE");
        helix.put("K", "#A05FA0");
        helix.put("M", "#EF10EF");
        helix.put("F", "#986798");
        helix.put("P", "#00FF00");
        helix.put("S", "#36C936");
        helix.put("T", "#47B847");
        helix.put("W", "#8A758A");
        helix.put("Y", "#21DE21");
        helix.put("V", "#857A85");
        helix.put("B", "#797986");
        helix.put("X", "#758A75");
        helix.put("Z", "#C936C9");
    }



    public static String getColour(char base) {
        return helix.get(Character.toString(base));
    }
}