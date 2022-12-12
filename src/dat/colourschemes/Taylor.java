package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Taylor{


    private static final Map<String, String> taylor;
    static {
        taylor = new HashMap<String, String>();
        taylor.put("A", "#CCFF00");
        taylor.put("R", "#0000FF");
        taylor.put("N", "#CC00FF");
        taylor.put("D", "#FF0000");
        taylor.put("C", "#FFFF00");
        taylor.put("Q", "#FF00CC");
        taylor.put("E", "#FF0066");
        taylor.put("G", "#FF9900");
        taylor.put("H", "#0066FF");
        taylor.put("I", "#66FF00");
        taylor.put("L", "#33FF00");
        taylor.put("K", "#6600FF");
        taylor.put("M", "#00FF00");
        taylor.put("F", "#00FF66");
        taylor.put("P", "#FFCC00");
        taylor.put("S", "#FF3300");
        taylor.put("T", "#FF6600");
        taylor.put("W", "#00CCFF");
        taylor.put("Y", "#00FFCC");
        taylor.put("V", "#99FF00");
        taylor.put("B", "#FFFFFF");
        taylor.put("X", "#FFFFFF");
        taylor.put("Z", "#FFFFFF");
    }



    public static String getColour(char base) {
        return taylor.get(Character.toString(base));
    }
}