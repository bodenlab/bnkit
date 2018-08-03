package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Zappo{


    private static final Map<String, String> zappo;
    static {
        zappo = new HashMap<String, String>();
        zappo.put("A", "#FFAFAF");
        zappo.put("R", "#6464FF");
        zappo.put("N", "#00FF00");
        zappo.put("D", "#FF0000");
        zappo.put("C", "#FFFF00");
        zappo.put("Q", "#00FF00");
        zappo.put("E", "#FF0000");
        zappo.put("G", "#FF00FF");
        zappo.put("H", "#6464FF");
        zappo.put("I", "#FFAFAF");
        zappo.put("L", "#FFAFAF");
        zappo.put("K", "#6464FF");
        zappo.put("M", "#FFAFAF");
        zappo.put("F", "#FFC800");
        zappo.put("P", "#FF00FF");
        zappo.put("S", "#00FF00");
        zappo.put("T", "#00FF00");
        zappo.put("W", "#FFC800");
        zappo.put("Y", "#FFC800");
        zappo.put("V", "#FFAFAF");
        zappo.put("B", "#FFFFFF");
        zappo.put("X", "#FFFFFF");
        zappo.put("Z", "#FFFFFF");
    }



    public static String getColour(char base) {
        return zappo.get(Character.toString(base));
    }
}